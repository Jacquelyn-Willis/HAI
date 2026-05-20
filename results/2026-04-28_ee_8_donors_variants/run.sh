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
  

# --- Conda ---
source "/hpc/packages/minerva-rocky9/miniforge3/26.1.1-3/miniforge/etc/profile.d/conda.sh"



# load module and environment; export this environment into a yml file and add to the repo for reproducibility; conda env export --name variant_calls > variant_calls_env.yml



## Directories 
data=/sc/arion/work/willij115/projects/HAI/data/2026-04-28_ee_8_donors_variants
data1=/sc/arion/work/willij115/projects/HAI/data/2026-04-23_align_pacbio_reads_to_reference
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-04-28_ee_8_donors_variants
results=/sc/arion/work/willij115/projects/HAI/results/2026-04-28_ee_8_donors_variants
results1=/sc/arion/work/willij115/projects/HAI/results/2026-04-23_align_pacbio_reads_to_reference


##functions

##1 call variants with deepvariant
make_bed_file () {
    cat <<EOF > "${data}/ig_loci.bed"
chr2	88807162	90310090	igk
chr22	22414484	23392002	igl
igh	0	1193129	igh
ighc	0	401750	
EOF
}

call_variants_w_deepvariant () {    

    module load singularity/3.11.0
    module load deepvariant/1.9.0

    mkdir -p "${results}/deep_variants"

    for file in "${results1}"/*.sorted.bam ; do
        sample_name=$(basename "$file" .sorted.bam)

        singularity exec \
            --cleanenv \
            --no-home \
            -B "${results1}:${results1},${results}:${results},${data1}:${data1},${data}:${data}" \
            "$DEEPVARIANT_SIF" \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref="${data1}/reference.fasta" \
            --reads="${file}" \
            --regions="${data}/ig_loci.bed" \
            --output_vcf="${results}/deep_variants/${sample_name}.vcf.gz" \
            --output_gvcf="${results}/deep_variants/${sample_name}.g.vcf.gz" \
             --sample_name="$sample_name" \
            --num_shards=${LSB_DJOB_NUMPROC:-8}

    done

    module purge

}




merge_vcf_files () {

    module load singularity/3.11.0
    module load glnexus

    mkdir -p "${results}/deep_variants_merged"

    gvcfs=$(ls "${results}/deep_variants/"*.g.vcf.gz)

    # clean previous GLnexus run
    rm -rf GLnexus.DB

    singularity exec \
        -B "${results}/deep_variants:/data,${data}:/bed" \
        docker://quay.io/mlin/glnexus:v1.2.2 \
        glnexus_cli \
        --config DeepVariantWGS \
        --bed /bed/ig_loci.bed \
        $gvcfs \
    | bcftools view - | bgzip -c > "${results}/deep_variants_merged/deepvariant.cohort.vcf.gz"

    module purge
}




##2 call variants with IGentotyper

index_bam_files () { 
    #sample download only came with .pbi files, so we need to index the bam files before we can use them for downstream analyses
    module load samtools

    for file in "${scratch}/aligned_bams"/*.bam ; do
        samtools index "$file"
    done

    module purge samtools
}

create_and_copy_reference_sa_index () {

    #IGentotyper needs the reference genome to be indexed with sawriter, so we will create the index and copy the reference fasta and index files to the IGenotyper data directory
    export SJOB_DEFALLOC=NONE

    set +u
    conda activate /sc/arion/work/willij115/test_env/envs/IGv2
    set -u
    
    #sawriter ${data1}/reference.fasta 
    cp ${data1}/reference.fasta* /sc/arion/work/willij115/test_env/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/
   
    
}


phase_bam_files_w_igenotyper () {

    #IGentotyper can phase variants in the bam files, so we will use it to phase the bam files and output phased bam files for downstream analyses
    export SJOB_DEFALLOC=NONE

    set +u
    conda activate /sc/arion/work/willij115/test_env/envs/IGv2
    set -u    

    module load minimap2
    
    sed -i -e '/chr7_.*/d' /sc/arion/work/willij115/test_env/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/target_regions.bed

    rm -r "${scratch}/igenotyper_run"
    mkdir -p "${scratch}/igenotyper_run"


    for file in "${data}/samples"/*.bam; do
        sample=$(basename "$file" .bam)
        outdir="${scratch}/igenotyper_run/${sample}"
        mkdir -p "$outdir"
        bsub -J "variant_calls" -P acc_oscarlr -q express -n 8 -W 12:00 -R "rusage[mem=100000] span[hosts=1]" -o "${outdir}/variant_calls.${sample}.%J.out.txt" -e "${outdir}/variant_calls.${sample}.%J.err.txt"  \
        IG phase "$file" "$outdir" --threads 8
       
    done

    
}

phase_bam_files_w_igenotyper_repeat () {

    #IGentotyper can phase variants in the bam files, so we will use it to phase the bam files and output phased bam files for downstream analyses
    export SJOB_DEFALLOC=NONE

    set +u
    conda activate /sc/arion/work/willij115/test_env/envs/IGv2
    set -u    

    module load minimap2
    
    sed -i -e '/chr7_.*/d' /sc/arion/work/willij115/test_env/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/target_regions.bed


    for file in "${data}/samples/m64152e_230602_034553.bc1046--bc1046.bam" "${data}/samples/m64407e_230722_040625.bc1057--bc1057.bam"; do
        sample=$(basename "$file" .bam)
        outdir="${scratch}/igenotyper_run/${sample}"
        mkdir -p "$outdir"
        bsub -J "variant_calls" -P acc_oscarlr -q express -n 8 -W 12:00 -R "rusage[mem=100000] span[hosts=1]" -o "${outdir}/variant_calls_phase.${sample}.%J.out.txt" -e "${outdir}/variant_calls_phase.${sample}.%J.err.txt"  \
        IG phase "$file" "$outdir" --threads 8 --tmp "${scratch}/igenotyper_run/${sample}/tmp_final"
       
    done

    
}
    

assemble_bam_files_w_igenotyper () {

    #IGentotyper can phase variants in the bam files, so we will use it to phase the bam files and output phased bam files for downstream analyses
    export SJOB_DEFALLOC=NONE

    set +u
    conda activate /sc/arion/work/willij115/test_env/envs/IGv2
    set -u    

    for file in "${data}/samples"/*.bam; do
            sample=$(basename "$file" .bam)
            outdir="${scratch}/igenotyper_run/${sample}"
            mkdir -p "$outdir"
            bsub -J "variant_calls" -P acc_oscarlr -q express -n 8 -W 12:00 -R "rusage[mem=100000] span[hosts=1]" -o "${outdir}/variant_calls_assembly.${sample}.%J.out.txt" -e "${outdir}/variant_calls_assembly.${sample}.%J.err.txt"  \
            IG assembly "$outdir" --threads 8
    done

}



IG_detect_w_igenotyper () {

    #IGentotyper can phase variants in the bam files, so we will use it to phase the bam files and output phased bam files for downstream analyses
    export SJOB_DEFALLOC=NONE

    set +u
    conda activate /sc/arion/work/willij115/test_env/envs/IGv2
    set -u    

    for file in "${data}/samples"/*.bam; do
            sample=$(basename "$file" .bam)
            outdir="${scratch}/igenotyper_run/${sample}"
            mkdir -p "$outdir"
            bsub -J "variant_calls" -P acc_oscarlr -q express -n 8 -W 12:00 -R "rusage[mem=100000] span[hosts=1]" -o "${outdir}/variant_calls_detect.${sample}.%J.out.txt" -e "${outdir}/variant_calls_detect.${sample}.%J.err.txt"  \
            IG detect "$outdir" 
    done


}

IG_detect_w_igenotyper_wo_2_samples () {

    #IGentotyper can phase variants in the bam files, so we will use it to phase the bam files and output phased bam files for downstream analyses
    export SJOB_DEFALLOC=NONE

    set +u
    conda activate /sc/arion/work/willij115/test_env/envs/IGv2
    set -u   

    for file in "${data}/samples"/*.bam; do
        sample=$(basename "$file" .bam)

        # skip these samples
        if [[ "$sample" == "m64152e_230602_034553.bc1046--bc1046" || "$sample" == "m64407e_230722_040625.bc1057--bc1057" ]]; then
            continue
        fi

        outdir="${scratch}/igenotyper_run/${sample}"
        mkdir -p "$outdir"

        bsub -J "variant_calls" \
            -P acc_oscarlr \
            -q express \
            -n 8 \
            -W 12:00 \
            -R "rusage[mem=100000] span[hosts=1]" \
            -o "${outdir}/variant_calls_detect.${sample}.%J.out.txt" \
            -e "${outdir}/variant_calls_detect.${sample}.%J.err.txt" \
            IG detect "$outdir" 
    done
}



## 3. Analyze deep variant output 

get_vcf_stats_of_deepvariant_output () {
    # here we can add code to analyze the output from deepvariant, such as calculating the number of variants called, the distribution of variant types, etc.
    module load bcftools

    bcftools stats "${results}/deep_variants_merged/deepvariant.cohort.vcf.gz" > "${scratch}/deepvariant.cohort.vcf.stats.txt"


    module purge
}



summarize_merged_vcf () {

  module load bcftools
  local vcf="$1"
  local out="${2:-summary.tsv}"
  local tmpdir
  tmpdir="$(mktemp -d)"

  echo "Working in $tmpdir"

  bcftools norm -m -any -Oz -o "$tmpdir/split.vcf.gz" "$vcf"
  bcftools index -t "$tmpdir/split.vcf.gz"

  bcftools +fill-tags "$tmpdir/split.vcf.gz" -Oz -o "$tmpdir/split.maf.vcf.gz" -- -t MAF
  bcftools index -t "$tmpdir/split.maf.vcf.gz"

  bcftools query -l "$tmpdir/split.maf.vcf.gz" > "$tmpdir/samples.txt"

  echo -e "sample\tSNP\tINDEL" > "$out"

  bcftools query -f '%REF\t%ALT[\t%GT]\n' "$tmpdir/split.maf.vcf.gz" | \
  awk -v samples_file="$tmpdir/samples.txt" '
    function nonref(gt,   a,n,i) {
      if (gt == "./." || gt == ".|." || gt == "." || gt == "") return 0
      n = split(gt, a, /[\/|]/)
      for (i = 1; i <= n; i++) {
        if (a[i] != "0" && a[i] != ".") return 1
      }
      return 0
    }
    function vartype(ref, alt) {
      if (length(ref) == 1 && length(alt) == 1) return "SNP"
      return "INDEL"
    }
    BEGIN {
      while ((getline s < samples_file) > 0) samples[++n] = s
      close(samples_file)
    }
    {
      type = vartype($1, $2)
      for (i = 3; i <= NF; i++) {
        if (nonref($i)) {
          sample = samples[i-2]
          cnt[sample, type]++
        }
      }
    }
    END {
      for (i = 1; i <= n; i++) {
        s = samples[i]
        print s "\t" cnt[s, "SNP"]+0 "\t" cnt[s, "INDEL"]+0
      }
    }' >> "$out"

  maf_count=$(bcftools view -i 'MAF>0.4' -H "$tmpdir/split.maf.vcf.gz" | wc -l)
  printf "# SNP_sites_MAF_gt_0.4\t%s\n" "$maf_count" >> "$out"

  rm -rf "$tmpdir"
  echo "$out"
}

plot_merged_vcf_summary () {
    module load python

    local summary_tsv="$1"
    local out_png="${2:-summary_counts.png}"

    python - "$summary_tsv" "$out_png" <<'PY'
import sys
import pandas as pd
import matplotlib.pyplot as plt

summary_tsv = sys.argv[1]
out_png = sys.argv[2]

df = pd.read_csv(summary_tsv, sep="\t", comment="#")
df = df[df["sample"].notna()]

ax = df.set_index("sample")[["SNP", "INDEL"]].plot(kind="bar", figsize=(12, 5))
ax.set_ylabel("Count")
ax.set_xlabel("Donor")

plt.tight_layout()
plt.savefig(out_png, dpi=300)
PY
}


##4. Analyze IGentotyper output


count_snps_per_sample_for_IG_loci () {
    local out_tsv="${scratch}/igenotyper_run/snv_counts.tsv"
    echo -e "sample\tIGH_all\tIGH_hetero_cnt\tIGK_all\tIGK_hetero_cnt\tIGL_all\tIGL_hetero_cnt" > "$out_tsv"

    # Strip commas from coordinates
    local IGH_CHR="igh"
    local IGH_START="1"
    local IGH_END="1193129"

    local IGK_CHR="chr2"
    local IGK_START="88807162"
    local IGK_END="90310090"

    local IGL_CHR="chr22"
    local IGL_START="22414484"
    local IGL_END="23392002"

    for SAMPLE_DIR in "${scratch}/igenotyper_run"/*; do
        [ -d "$SAMPLE_DIR" ] || continue
        SAMPLE=$(basename "$SAMPLE_DIR")

        VCF=$(find "$SAMPLE_DIR/variants/" -name "snvs_assembly.vcf" -o -name "snvs_assembly.vcf.gz" 2>/dev/null | head -1)
        [ -z "$VCF" ] && continue

        if [[ "$VCF" == *.vcf.gz ]]; then
            READ_CMD="zcat"
        else
            READ_CMD="cat"
        fi

        IGH_ALL=$($READ_CMD "$VCF" | awk -v chr="$IGH_CHR" -v start="$IGH_START" -v end="$IGH_END" '
            !/^#/ && $1==chr && $2>=start && $2<=end && length($4)==1 && length($5)==1 { c++ }
            END { print c+0 }
        ')

        IGH_MAF=$($READ_CMD "$VCF" | awk -v chr="$IGH_CHR" -v start="$IGH_START" -v end="$IGH_END" '
            !/^#/ && $1==chr && $2>=start && $2<=end && length($4)==1 && length($5)==1 {
                split($9, fmt, ":")
                split($10, sample, ":")
                gt_i = 0
                for (i=1; i<=length(fmt); i++) {
                    if (fmt[i] == "GT") { gt_i = i; break }
                }
                if (gt_i == 0) next
                gt = sample[gt_i]
                gsub(/\|/, "/", gt)
                if (gt == "0/1" || gt == "1/0") c++
            }
            END { print c+0 }
        ')

        IGK_ALL=$($READ_CMD "$VCF" | awk -v chr="$IGK_CHR" -v start="$IGK_START" -v end="$IGK_END" '
            !/^#/ && $1==chr && $2>=start && $2<=end && length($4)==1 && length($5)==1 { c++ }
            END { print c+0 }
        ')

        IGK_MAF=$($READ_CMD "$VCF" | awk -v chr="$IGK_CHR" -v start="$IGK_START" -v end="$IGK_END" '
            !/^#/ && $1==chr && $2>=start && $2<=end && length($4)==1 && length($5)==1 {
                split($9, fmt, ":")
                split($10, sample, ":")
                gt_i = 0
                for (i=1; i<=length(fmt); i++) {
                    if (fmt[i] == "GT") { gt_i = i; break }
                }
                if (gt_i == 0) next
                gt = sample[gt_i]
                gsub(/\|/, "/", gt)
                if (gt == "0/1" || gt == "1/0") c++
            }
            END { print c+0 }
        ')

        IGL_ALL=$($READ_CMD "$VCF" | awk -v chr="$IGL_CHR" -v start="$IGL_START" -v end="$IGL_END" '
            !/^#/ && $1==chr && $2>=start && $2<=end && length($4)==1 && length($5)==1 { c++ }
            END { print c+0 }
        ')

        IGL_MAF=$($READ_CMD "$VCF" | awk -v chr="$IGL_CHR" -v start="$IGL_START" -v end="$IGL_END" '
            !/^#/ && $1==chr && $2>=start && $2<=end && length($4)==1 && length($5)==1 {
                split($9, fmt, ":")
                split($10, sample, ":")
                gt_i = 0
                for (i=1; i<=length(fmt); i++) {
                    if (fmt[i] == "GT") { gt_i = i; break }
                }
                if (gt_i == 0) next
                gt = sample[gt_i]
                gsub(/\|/, "/", gt)
                if (gt == "0/1" || gt == "1/0") c++
            }
            END { print c+0 }
        ')

        echo -e "${SAMPLE}\t${IGH_ALL}\t${IGH_MAF}\t${IGK_ALL}\t${IGK_MAF}\t${IGL_ALL}\t${IGL_MAF}" >> "$out_tsv"
    done
}


plot_ig_locus_summary () {
    module load python

    local summary_tsv="${scratch}/igenotyper_run/snv_counts.tsv"
    local out_png="${2:-ig_loci_snv_counts.png}"

    python - "$summary_tsv" "${scratch}/igenotyper_run/${out_png}" <<'PY'
import sys
import pandas as pd
import matplotlib.pyplot as plt

summary_tsv = sys.argv[1]
out_png = sys.argv[2]

df = pd.read_csv(summary_tsv, sep="\t", comment="#")
df = df[df["sample"].notna()].copy()

total_cols = ["IGH_all", "IGK_all", "IGL_all"]
maf_cols = ["IGH_hetero_cnt", "IGK_hetero_cnt", "IGL_hetero_cnt"]

for c in total_cols + maf_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

# Top: total SNP counts
df.set_index("sample")[total_cols].plot(
    kind="bar",
    ax=axes[0],
    width=0.85
)
axes[0].set_ylabel("Total SNPs")
axes[0].set_title("SNP counts in IG loci")

# Bottom: MAF > 0.4 proxy counts
df.set_index("sample")[maf_cols].plot(
    kind="bar",
    ax=axes[1],
    width=0.85
)
axes[1].set_ylabel("Het SNPs")
axes[1].set_xlabel("Donor")
axes[1].set_title("SNPs with MAF > 0.4 proxy")

for ax in axes:
    ax.legend(title="", frameon=False)

axes[1].tick_params(axis="x", rotation=45)
for label in axes[1].get_xticklabels():
    label.set_horizontalalignment("right")

plt.tight_layout()
plt.savefig(out_png, dpi=300)
PY
}


##function calls 

#make_bed_file
#call_variants_w_deepvariant
#merge_vcf_files

#index_bam_files
#create_and_copy_reference_sa_index 
#phase_bam_files_w_igenotyper
#phase_bam_files_w_igenotyper_repeat
#assemble_bam_files_w_igenotyper 
#IG_detect_w_igenotyper
#IG_detect_w_igenotyper_wo_2_samples


#get_vcf_stats_of_deepvariant_output
#summarize_merged_vcf "${results}/deep_variants_merged/deepvariant.cohort.vcf.gz"  "${scratch}/merged_vcf_summary.tsv"
#plot_merged_vcf_summary "${scratch}/merged_vcf_summary.tsv" "${scratch}/merged_vcf_summary.png"

#count_snps_per_sample_for_IG_loci
plot_ig_locus_summary