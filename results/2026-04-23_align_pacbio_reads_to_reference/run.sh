#!/bin/bash
#BSUB -J pacbio_reads                # Job name
#BSUB -P acc_oscarlr                    # Project allocation
#BSUB -q express                       # Queue name
#BSUB -n 8                              # 8 compute cores
#BSUB -R "rusage[mem=10000]"             # 10 GB per core → 80 GB total
#BSUB -R "span[hosts=1]"                # All cores on the same node
#BSUB -W 12:00                          # 85 hour wall-time limit
#BSUB -o pacbio_reads   .%J.out       # STDOUT log
#BSUB -eo pacbio_reads   .%J.err      # STDERR log
#BSUB -L /bin/bash


#to cause any big bugs to fail out of the script
set -eou pipefail 
set -x


# --- Conda ---
source "/hpc/packages/minerva-rocky9/miniforge3/26.1.1-3/miniforge/etc/profile.d/conda.sh"
module load anaconda3/latest

module load sratoolkit/3.0.1
module load samtools/1.21
module load minimap2/2.24

 #export this environment into a yml file and add to the repo for reproducibility; conda env export --name pacbio_reads > pacbio_reads_env.yml



work=/sc/arion/work/willij115/projects/HAI/data/2026-04-23_align_pacbio_reads_to_reference
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-04-23_align_pacbio_reads_to_reference
results=/sc/arion/work/willij115/projects/HAI/results/2026-04-23_align_pacbio_reads_to_reference


convert_sra_to_fastq () {

        for file in SRR33916016 SRR33916038 SRR33916044  SRR33916105   SRR33916134  SRR33916143  SRR33916144  SRR33916146 ; do
            fasterq-dump --outdir $scratch $file
        done

}

create_ref_minimizer_index () {
    minimap2 -x map-hifi -d ${scratch}/ref.mmi -t 8 ${work}/reference.fasta


}


align_donor_reads_to_reference () {
    for file in $scratch/*.fastq ; do
        
        minimap2 -t 8 -ax map-hifi ${scratch}/ref.mmi $file \
        | samtools sort -@ 8 -O BAM -o ${results}/$(basename "$file" .fastq).sorted.bam

        samtools index ${results}/$(basename "$file" .fastq).sorted.bam

    done
}


#function calls

#convert_sra_to_fastq
create_ref_minimizer_index
align_donor_reads_to_reference
