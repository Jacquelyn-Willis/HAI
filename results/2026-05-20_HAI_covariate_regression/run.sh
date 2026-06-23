#!/bin/bash
#BSUB -J variant_calls                # Job name
#BSUB -P acc_oscarlr                    # Project allocation
#BSUB -q express                  # Queue name
#BSUB -n 8                              # 8 compute cores
#BSUB -R "rusage[mem=50000]"             # 10 GB per core → 80 GB total
#BSUB -R "span[hosts=1]"                # All cores on the same node
#BSUB -W 12:00                          # 85 hour wall-time limit
#BSUB -o "variant_calls.%J.out.txt"       # STDOUT log
#BSUB -eo "variant_calls.%J.err.txt"      # STDERR log
#BSUB -L /bin/bash


## Cause any big bugs to fail out of the script
set -eou pipefail 
set -x
  

# --- Conda ---
source "/hpc/packages/minerva-rocky9/miniforge3/26.1.1-3/miniforge/etc/profile.d/conda.sh"



# load module and environment; export this environment into a yml file and add to the repo for reproducibility; conda env export --name variant_calls > variant_calls_env.yml
source "/hpc/packages/minerva-rocky9/miniforge3/26.1.1-3/miniforge/etc/profile.d/conda.sh"
module load anaconda3/latest

## Directories 
data=/sc/arion/work/willij115/projects/HAI/data/2026-05-20_HAI_covariate_regression
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-05-20_HAI_covariate_regression
results=/sc/arion/work/willij115/projects/HAI/results/2026-05-20_HAI_covariate_regression



#1. merge immunespace HAI data tables into one table 


#2. group by study and arm, and summarize the number of participants, events, and assays per study and arm


#3. merge the summarized table with the immunespace HAI data tables to create a final table with all relevant information for each study and arm and save the final table to a CSV file in the results directory


#4. run the covariate regression analysis using the final table and save the results to a CSV file in the results directory











### Function calls