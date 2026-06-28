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


    wget -c -O "${scratch}/immunespaceHAI_studies_tables.csv" 'https://immunespace.org/api_kb/get_full_studies_data/?assay_curie_filter=OBI%3A0000875+OBI%3A0000872&ordering_studies=investigation_id&ordering_studies_dir=desc&ordering_tab=study_tab'
    wget -c -O "${scratch}/immunespaceHAI_arms_tables.csv" 'https://immunespace.org/api_kb/get_full_arms_data/?assay_curie_filter=OBI%3A0000875+OBI%3A0000872'
    wget -c -O "${scratch}/immunespaceHAI_participants_tables.csv" 'https://immunespace.org/api_kb/get_full_participants_data/?assay_curie_filter=OBI%3A0000875+OBI%3A0000872'
    wget -c -O "${scratch}/immunespaceHAI_events_tables.csv" 'https://immunespace.org/api_kb/get_full_events_data/?assay_curie_filter=OBI%3A0000875+OBI%3A0000872'
    wget -c -O "${scratch}/immunespaceHAI_assays_tables.csv" 'https://immunespace.org/api_kb/get_full_assays_data/?assay_curie_filter=OBI%3A0000875+OBI%3A0000872'

    
}


export_immunespace_datatools_tables () {
    #export the data tables to the data directory
    cp "/Users/jwillis/Downloads/demographics_2026-06-28_13-34-03.csv" "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/datatools_demographic_Table.csv"
    cp "/Users/jwillis/Downloads/hai_2026-06-28_13-38-11.csv" "/Users/jwillis/minerva_scratch/projects/HAI/2026-05-20_HAI_covariate_regression/datatools_HAI_Table.csv"
   
}

##function calls

#get_immunespace_HAI_data_tables
export_immunespace_datatools_tables