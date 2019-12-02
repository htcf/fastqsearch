#!/usr/bin/env python3

import argparse
import gzip
from multiprocessing import Pool, Queue, Process
import re
import sys
import time

def get_seqs(fastq_file, chunk_size=40):
    seqs = []
    seq_re = re.compile(b'@.*\n(.*)\n', re.MULTILINE)
    with gzip.open(fastq_file, 'rb') as f:
        data = b''
        size = 0 
        while True:
            data += f.read(int(1024*1024*chunk_size))
            if len(data) == 0: break
            for m in seq_re.finditer(data):
                seqs.append(m.groups(0)[0])
                size += 1
                if size % 1000 == 0:
                    yield seqs
                    size = 0
                    seqs = []
            data = data[m.end():]
        if seqs:
            yield seqs

def worker(task_queue, done_queue, keys):
    lookup = dict([(k,0) for k in keys])
    for seqs, seq_length, key_length in iter(task_queue.get, 'STOP'):
        for seq in seqs:
            for i in range(0, seq_length-key_length):
                try:
                    lookup[seq[i:i+key_length]] += 1
                except:
                    pass
    done_queue.put(lookup)
    
def go(fastq_file, guide_file, num_workers):
    task_queue = Queue(50)
    done_queue = Queue()
    lookup = dict([(x.strip(),0) for x in open(guide_file, 'rb').readlines()])
    workers = []

    for i in range(num_workers):
        p = Process(target=worker, args=(task_queue,done_queue,lookup.keys()))
        p.start()
        workers.append(p)

    key_length = len(next(iter(lookup.keys())))
    seq_length = None
    for i, seqs in enumerate(get_seqs(fastq_file)):
        if seq_length is None:
            seq_length = len(seqs[0])
        task_queue.put((seqs, seq_length, key_length))
    task_queue.put((seqs, seq_length, key_length))
    

    for i in range(num_workers):
        task_queue.put('STOP')

    for i in range(num_workers):
        for k,v in done_queue.get().items():
            lookup[k] += v
    return lookup


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('guide_file')
    p.add_argument('fastq_file')
    p.add_argument('-w', '--workers', default=4, type=int)
    args = p.parse_args()
    
    lookup = go(args.fastq_file, args.guide_file, args.workers)
    for k in lookup.keys():
        sys.stdout.write("{} {}\n".format(k.decode(), lookup[k]))
