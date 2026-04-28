#!/bin/bash
#BSUB -P acc_oscarlr         
#BSUB -q interactive              
#BSUB -n 1                      # one core is plenty for downloading
#BSUB -W 12:00                    
#BSUB -J pacbio_reads              # descriptive job name
#BSUB -o pacbio_reads.%J.out       # stdout 
#BSUB -e pacbio_reads.%J.err       # stderr 
#BSUB -Is /bin/bash



#to cause any big bugs to fail out of the script
set -eou pipefail 
set -x



work=/sc/arion/work/willij115/projects/HAI/data/2026-04-23_align_pacbio_reads_to_reference
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-04-23_align_pacbio_reads_to_reference



# --- Get the 1 kg VCF genome files ---
get_index_and_manifest_file () {

    cd $work || exit 1. #do this command or exit if it fails

    wget -c https://raw.githubusercontent.com/franklin-hiciano/RodrioData/refs/heads/main/datasets/2026-Light_EE_NatComm/2026-Light_EE_NatComm.std.index
    wget -c "https://raw.githubusercontent.com/franklin-hiciano/RodrioData/refs/heads/main/datasets/2026-Light_EE_NatComm/metadata-15346978-processed-ok%20(2).tsv" samtools faidx reference.fasta
}

parse_files () {
    
    #awk: A tool for processing text files row by row and column by column.
    #Sets the field separator to a tab, So each line is split into columns:
    #$6 ~ /^(CE0005604|CE0007821|CE0005779|CE0006233|CE0006435|CE0006481|CE0007070|CE0006634)/' means: If the 6th column matches regex any of these sample IDs (regular expression), keep the line.    


    awk -F'\t' ' $6 ~ /^(CE0005604_genomic_1|CE0007821_genomic_1|CE0005779_genomic_1|CE0006233_genomic_1|CE0006435_genomic_1|CE0006481_genomic_1|CE0007070_genomic_1|CE0006634_genomic_1)$/' "metadata-15346978-processed-ok (2).tsv" > metadata_filtered.tsv 

    cut -f1 metadata_filtered.tsv > sample_ids.txt

    awk -F'\t' 'NR==FNR {ids[$1]; next} $1 in ids' sample_ids.txt 2026-Light_EE_NatComm.std.index > index_filtered.std.index

}

# ---------- HTTP -----------
function HTTP_download_one_file {
    #-O saves to current directory; -o allows you to specify the output file name and path; -L follows redirects; -c continues a previous download if it was interrupted.

    curl -L -c "$1" -o "$2"
}


# This function processes and feeds the links from the .std.index files into the download functions one by one, which would be inconvenient to have to code up every time.
function download_files {
        index_file="$(realpath "$1")"
        out_dir="$(realpath $2)"
        method="$3"
        region="${4:-}"
        cut -f9 "${index_file}" | tr ';' '\n' | while read -r url; do
            echo ${method}
            ${method}_download_one_file $url "${out_dir}/$(basename ${url})" $region
        done
}

get_ref_fasta () {

    cd $work || exit 1. #do this command or exit if it fails
    wget -c http://immunogenomics.louisville.edu/wasp/251106/reference.fasta
    wget -c https://raw.githubusercontent.com/Watson-IG/immune_receptor_genomics/refs/heads/main/251106/reference.fasta.fai

}
#function calls

#get_index_and_manifest_file 
#parse_files
#download_files index_filtered.std.index ${scratch} HTTP
get_ref_fasta
