#!/bin/bash
#BSUB -J variant_calls                # Job name
#BSUB -P acc_oscarlr                    # Project allocation
#BSUB -q express                       # Queue name
#BSUB -n 8                              # 8 compute cores
#BSUB -R "rusage[mem=10000]"             # 10 GB per core → 80 GB total
#BSUB -R "span[hosts=1]"                # All cores on the same node
#BSUB -W 12:00                          # 85 hour wall-time limit
#BSUB -o variant_calls.%J.out       # STDOUT log
#BSUB -eo variant_calls.%J.err      # STDERR log
#BSUB -L /bin/bash


## Cause any big bugs to fail out of the script
set -eou pipefail 
set -x
  

## Conda 
export JAVA_HOME="/hpc/packages/minerva-rocky9/java/21.0.4/jdk-21.0.4"
export JAVA_LD_LIBRARY_PATH="${JAVA_HOME}/lib"
eval "$(micromamba shell hook --shell bash)"

module load micromamba/1.5.3-0

# load module and environment; export this environment into a yml file and add to the repo for reproducibility; conda env export --name variant_calls > variant_calls_env.yml

module load deepvariant/1.9.0
module

## Directories 
data=/sc/arion/work/willij115/projects/HAI/data/2026-04-28_variant_calls
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-04-28_variant_calls
results=/sc/arion/work/willij115/projects/HAI/results/2026-04-28_variant_calls
results1=/sc/arion/work/willij115/projects/HAI/results/2026-04-23_align_pacbio_reads_to_reference


