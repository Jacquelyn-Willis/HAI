#!/bin/bash
#BSUB -P acc_oscarlr         
#BSUB -q interactive              
#BSUB -n 1                      # one core is plenty for downloading
#BSUB -W 12:00                    
#BSUB -J HAI              # descriptive job name
#BSUB -o "variant.%J.out.txt"       # stdout 
#BSUB -e "variant.%J.err.txt"       # stderr 
#BSUB -Is /bin/bash



#to cause any big bugs to fail out of the script
set -eou pipefail 
set -x

## Directories 
data=/sc/arion/work/willij115/projects/HAI/data/2026-05-20_HAI_covariate_regression
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-05-20_HAI_covariate_regression


get_immunespace_HAI_data_tables () {


    wget -c -O "${scratch}/immunespaceHAI_studies_tables.csv" 'https://immunespace.org/api_kb/get_full_studies_data/?assay_curie_filter=OBI:0000875+OBI:0000872'
    wget -c -O "${scratch}/immunespaceHAI_arms_tables.csv" 'https://immunespace.org/api_kb/get_full_arms_data/?assay_curie_filter=OBI:0000875+OBI:0000872'
    wget -c -O "${scratch}/immunespaceHAI_participants_tables.csv" 'https://immunespace.org/api_kb/get_full_participants_data/?assay_curie_filter=OBI:0000875+OBI:0000872'
    wget -c -O "${scratch}/immunespaceHAI_events_tables.csv" 'https://immunespace.org/api_kb/get_full_events_data/?assay_curie_filter=OBI:0000875+OBI:0000872'
    wget -c -O "${scratch}/immunespaceHAI_assays_tables.csv" 'https://immunespace.org/api_kb/get_full_assays_data/?assay_curie_filter=OBI:0000875+OBI:0000872'

}



##function calls

get_immunespace_HAI_data_tables