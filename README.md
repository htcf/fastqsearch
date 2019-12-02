# fastqsearch

Slurm sbatch example:

```
#!/bin/bash

#SBATCH --cpus-per-task=5

./fastqsearch.py -w $(( SLURM_CPUS_PER_TASK - 1 )) guide.txt big_R1.fastq.gz
```
