#!/bin/bash
#BSUB -P acc_oscarlr         
#BSUB -q interactive              
#BSUB -n 1                      # one core is plenty for downloading
#BSUB -W 12:00                    
#BSUB -J variant              # descriptive job name
#BSUB -o variant.%J.out       # stdout 
#BSUB -e variant.%J.err       # stderr 
#BSUB -Is /bin/bash



#to cause any big bugs to fail out of the script
set -eou pipefail 
set -x

## Directories 
data=/sc/arion/work/willij115/projects/HAI/data/2026-04-28_variant_calls
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-04-28_variant_calls



get_bed_file () {
    wget -c https://raw.githubusercontent.com/Watson-IG/immune_receptor_genomics/refs/heads/main/251106/gene.bed\
        -O "${scratch}/gene.bed"
}

get_bam_files () {


	cp /sc/arion/projects/oscarlr/rodrio10/to_share/for_jacquelyn/2026-05-04_EE_8_serum_donors/samples.zip  /sc/arion/work/willij115/projects/HAI/data/2026-04-28_variant_calls/

} 





##function calls

#get_bed_file
get_bam_files
