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
data=/sc/arion/work/willij115/projects/HAI/data/2026-04-28_ee_8_donors_variants
data1=/sc/arion/work/willij115/projects/HAI/data/2026-04-23_align_pacbio_reads_to_reference
scratch=/sc/arion/scratch/willij115/projects/HAI/2026-04-28_ee_8_donors_variants
results=/sc/arion/work/willij115/projects/HAI/results/2026-04-28_ee_8_donors_variants
results1=/sc/arion/work/willij115/projects/HAI/results/2026-04-23_align_pacbio_reads_to_reference


##functions

##1 call variants with deepvariant


align_bam_files_to_reference_for_deepvariant_input123 () {

    conda activate /sc/arion/work/willij115/test_env/envs/pbmm2_env

    #pbmm2 index -j 8 ${data1}/reference.fasta ${scratch}/ref.mmi 

    rm -rf "${scratch}/aligned_bams"
    mkdir -p "${scratch}/aligned_bams"

    for file in "${data}"/samples/*.bam ; do

        # Extract the base name without the .bam extension
        local sample_name=$(basename "$file" .bam)
        
        echo "Aligning ${sample_name}..."
            
        bsub -J "pbmm2_${sample_name}" \
                -P acc_oscarlr \
                -q express \
                -n 8 \
                -W 12:00 \
                -R "rusage[mem=50000] span[hosts=1]" \
                -o "${scratch}/aligned_bams/pbmm2_align.${sample_name}.%J.out.txt" \
                -e "${scratch}/aligned_bams/pbmm2_align.${sample_name}.%J.err.txt"  \
                pbmm2 align  "${scratch}/ref.mmi"  "$file" "${scratch}/aligned_bams/${sample_name}_aligned.sorted.bam"  --sort -j 8

    done

}

align_bam_files_to_reference_for_deepvariant_input () {

    conda activate /sc/arion/work/willij115/test_env/envs/pbmm2_env

    rm -rf "${scratch}/aligned_bams"
    mkdir -p "${scratch}/aligned_bams"

    for file in "${data}"/samples/*.bam ; do
        local sample_name
        sample_name=$(basename "$file" .bam)

        echo "Aligning ${sample_name}..."

        bsub -J "pbmm2_${sample_name}" \
             -P acc_oscarlr \
             -q express \
             -n 8 \
             -W 12:00 \
             -R "rusage[mem=50000] span[hosts=1]" \
             -o "${scratch}/aligned_bams/pbmm2_align.${sample_name}.%J.out.txt" \
             -e "${scratch}/aligned_bams/pbmm2_align.${sample_name}.%J.err.txt" \
             bash -lc "
                 set -euo pipefail
                 export TMPDIR='${scratch}/tmp/${sample_name}'
                 mkdir -p \"\$TMPDIR\"
                 pbmm2 align '${scratch}/ref.mmi' '$file' \
                     '${scratch}/aligned_bams/${sample_name}_aligned.sorted.bam' \
                     --sort -j 8
             "
    done
}

make_bed_file () {
    cat <<EOF > "${data}/ig_loci.bed"
chr2	88807162	90310090	igk
chr22	22414484	23392002	igl
igh	0	1193129	igh
ighc	0	401750	
EOF
}

call_variants_w_deepvariant () {    

    # call in SLEEP job

    module load singularity/3.11.0
    module load deepvariant/1.9.0

    rm -rf "${scratch}/deep_variant_calls"
    mkdir -p "${scratch}/deep_variant_calls"

    for file in "${scratch}/aligned_bams"/*_aligned.sorted.bam ; do
        sample_name=$(basename "$file" _aligned.sorted.bam)


        singularity exec \
            --cleanenv \
            --no-home \
            -B "${results1}:${results1},${results}:${results},${data1}:${data1},${data}:${data},${scratch}:${scratch}" \
            "$DEEPVARIANT_SIF" \
            /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref="${data1}/reference.fasta" \
            --reads="${file}" \
            --regions="${data}/ig_loci.bed" \
            --output_vcf="${scratch}/deep_variant_calls/${sample_name}.vcf.gz" \
            --output_gvcf="${scratch}/deep_variant_calls/${sample_name}.g.vcf.gz" \
             --sample_name="$sample_name" \
            --num_shards=${LSB_DJOB_NUMPROC:-8}

    done

    module purge

}


merge_deepvariant_vcf_files () {

    module load singularity/3.11.0
    module load bcftools
    module load htslib

    mkdir -p "${results}/deepvariants_merged_new_bam"
    rm -rf GLnexus.DB

    # 1. Get just the filenames and add the /data/ prefix to them
    gvcfs=$(cd "${results}/deep_variant_calls" && ls *.g.vcf.gz | sed 's|^|/data/|')

    # 2. Run Singularity normally
    singularity exec \
        -B "${results}/deep_variant_calls:/data,${data}:/bed" \
        "/hpc/users/willij115/glnexus_v1.2.2.sif" \
        glnexus_cli \
        --config DeepVariantWGS \
        --bed /bed/ig_loci.bed \
        $gvcfs \
    | bcftools view - | bgzip -c > "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"

    bcftools index -t "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"

    module purge
}











## 2. Analyze deep variant output 

get_vcf_stats_of_deepvariant_output () {
    # here we can add code to analyze the output from deepvariant, such as calculating the number of variants called, the distribution of variant types, etc.
    module load bcftools

    bcftools stats "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz" > "${scratch}/deepvariant.cohort.vcf.stats.txt"


    module purge
}




summarize_merged_deepvariant_vcf_split () {
  module load bcftools

  local vcf="$1"
  local out="${2:-summary.tsv}"
  local tmpdir
  tmpdir="$(mktemp -d)"
  trap 'rm -rf "$tmpdir"' RETURN

  local IGH_CHR="igh"
  local IGH_START="1"
  local IGH_END="1193129"

  local IGK_CHR="chr2"
  local IGK_START="88807162"
  local IGK_END="90310090"

  local IGL_CHR="chr22"
  local IGL_START="22414484"
  local IGL_END="23392002"

  rm -f "${vcf}.tbi" "${vcf}.csi"
  bcftools index -f -c "$vcf"

  echo "Working in $tmpdir" >&2
  echo -e "locus\tsample\tSNP\tINDEL\tSNP_sites_MAF_gt_0.4" > "$out"

  for locus in IGH IGK IGL; do
    local chr start end region
    local region_vcf split_vcf maf_vcf samples_file maf_count

    case "$locus" in
      IGH)
        chr="$IGH_CHR"; start="$IGH_START"; end="$IGH_END"
        ;;
      IGK)
        chr="$IGK_CHR"; start="$IGK_START"; end="$IGK_END"
        ;;
      IGL)
        chr="$IGL_CHR"; start="$IGL_START"; end="$IGL_END"
        ;;
    esac

    region="${chr}:${start}-${end}"
    region_vcf="$tmpdir/${locus}.region.vcf.gz"
    split_vcf="$tmpdir/${locus}.split.vcf.gz"
    maf_vcf="$tmpdir/${locus}.split.maf.vcf.gz"
    samples_file="$tmpdir/${locus}.samples.txt"

    bcftools view -r "$region" -Oz -o "$region_vcf" "$vcf"
    bcftools index -t "$region_vcf"

    bcftools norm -m -any -Oz -o "$split_vcf" "$region_vcf"
    bcftools index -t "$split_vcf"

    bcftools +fill-tags "$split_vcf" -Oz -o "$maf_vcf" -- -t MAF
    bcftools index -t "$maf_vcf"

    bcftools query -l "$maf_vcf" > "$samples_file"

    maf_count=$(bcftools view -i 'MAF>0.4' -H "$maf_vcf" | wc -l | awk '{print $1}')

    bcftools query -f '%REF\t%ALT[\t%GT]\n' "$maf_vcf" | \
    awk -v locus="$locus" -v maf_count="$maf_count" -v samples_file="$samples_file" '
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
          print locus "\t" s "\t" cnt[s, "SNP"]+0 "\t" cnt[s, "INDEL"]+0 "\t" maf_count+0
        }
      }' >> "$out"
  done

  echo "$out"
}

plot_merged_deepvariant_vcf_summary () {
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
df = df[df["sample"].notna()].copy()

loci = [x for x in ["IGH", "IGK", "IGL"] if x in set(df["locus"])]
if not loci:
    loci = sorted(df["locus"].dropna().unique())

fig, axes = plt.subplots(len(loci), 1, figsize=(14, max(4, 4 * len(loci))), sharex=False)
if len(loci) == 1:
    axes = [axes]

# define one consistent sample order
sample_order = sorted(df["sample"].unique())

for ax, locus in zip(axes, loci):
    sub = (
        df[df["locus"] == locus]
        .set_index("sample")
        .reindex(sample_order)[["SNP", "INDEL"]]
    )

    sub.plot(kind="bar", ax=ax)
    ax.set_title(locus)
    ax.set_ylabel("Count")
    ax.set_xlabel("Sample")
    ax.legend(loc="upper right")

plt.tight_layout()
plt.savefig(out_png, dpi=300)
PY
}




#summarize_merged_deepvariant_vcf_split "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"  "${scratch}/merged_vcf_summary_split_new_bam.tsv"
#plot_merged_deepvariant_vcf_summary "${scratch}/merged_vcf_summary_split_new_bam.tsv" "${scratch}/merged_vcf_summary_split_new_bam.png"








##3 call variants with IGentotyper

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
    
    #sed -i -e '/chr7_.*/d' /sc/arion/work/willij115/test_env/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/target_regions.bed

    rm -r "${scratch}/igenotyper_run"
    mkdir -p "${scratch}/igenotyper_run"


    for file in "${data}/samples"/*.bam; do
        sample=$(basename "$file" .bam)
        outdir="${scratch}/igenotyper_run/${sample}"
        mkdir -p "$outdir"
        bsub -J "variant_calls_phase" -P acc_oscarlr -q express -n 8 -W 12:00 -R "rusage[mem=50000] span[hosts=1]" -o "${outdir}/variant_calls_phase.${sample}.%J.out.txt" -e "${outdir}/variant_calls_phase.${sample}.%J.err.txt"  \
        IG phase "$file" "$outdir" --threads 8 --tmp "${scratch}/igenotyper_run/${sample}/tmp_final"
       
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
            bsub -J "variant_calls_assembly" -P acc_oscarlr -q express -n 8 -W 12:00 -R "rusage[mem=50000] span[hosts=1]" -o "${outdir}/variant_calls_assembly.${sample}.%J.out.txt" -e "${outdir}/variant_calls_assembly.${sample}.%J.err.txt"  \
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
            -R "rusage[mem=50000] span[hosts=1]" \
            -o "${outdir}/variant_calls_detect.${sample}.%J.out.txt" \
            -e "${outdir}/variant_calls_detect.${sample}.%J.err.txt" \
            IG detect "$outdir" 
    done
}


#phase_bam_files_w_igenotyper
#assemble_bam_files_w_igenotyper
IG_detect_w_igenotyper_wo_2_samples





##4. Analyze IGentotyper output

merge_igenotyper_vcf_files () {

    module load singularity/3.11.0
    module load bcftools
    module load htslib

    mkdir -p "${results}/deepvariants_merged_new_bam"
    rm -rf GLnexus.DB

    # 1. Get just the filenames and add the /data/ prefix to them
    gvcfs=$(cd "${results}/deep_variant_calls" && ls *.g.vcf.gz | sed 's|^|/data/|')

    # 2. Run Singularity normally
    singularity exec \
        -B "${results}/deep_variant_calls:/data,${data}:/bed" \
        "/hpc/users/willij115/glnexus_v1.2.2.sif" \
        glnexus_cli \
        --config DeepVariantWGS \
        --bed /bed/ig_loci.bed \
        $gvcfs \
    | bcftools view - | bgzip -c > "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"

    bcftools index -t "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"

    module purge
}







count_igenotyper_snps_per_sample_for_IG_loci () {
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

        IGH_hetero=$($READ_CMD "$VCF" | awk -v chr="$IGH_CHR" -v start="$IGH_START" -v end="$IGH_END" '
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

        IGK_hetero=$($READ_CMD "$VCF" | awk -v chr="$IGK_CHR" -v start="$IGK_START" -v end="$IGK_END" '
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

        IGL_hetero=$($READ_CMD "$VCF" | awk -v chr="$IGL_CHR" -v start="$IGL_START" -v end="$IGL_END" '
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

        echo -e "${SAMPLE}\t${IGH_ALL}\t${IGH_hetero}\t${IGK_ALL}\t${IGK_hetero}\t${IGL_ALL}\t${IGL_hetero}" >> "$out_tsv"
    done
}


plot_igenotyper_ig_locus_summary () {
    module load python

    local summary_tsv="${scratch}/igenotyper_run/snv_counts.tsv"
    local out_png="${2:-ig_loci_snv_counts_igenotyper.png}"

    python - "$summary_tsv" "${scratch}/igenotyper_run/${out_png}" <<'PY'
import sys
import pandas as pd
import matplotlib.pyplot as plt

summary_tsv = sys.argv[1]
out_png = sys.argv[2]

df = pd.read_csv(summary_tsv, sep="\t", comment="#")
df = df[df["sample"].notna()].copy()

total_cols = ["IGH_all", "IGK_all", "IGL_all"]
hetero_cols = ["IGH_hetero_cnt", "IGK_hetero_cnt", "IGL_hetero_cnt"]

for c in total_cols + hetero_cols:
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

# Bottom: hetero counts
df.set_index("sample")[hetero_cols].plot(
    kind="bar",
    ax=axes[1],
    width=0.85
)
axes[1].set_ylabel("Het SNPs")
axes[1].set_xlabel("Donor")
axes[1].set_title("count of Heterozygous SNPs")

for ax in axes:
    ax.legend(title="", frameon=False)

axes[1].tick_params(axis="x", rotation=45)
for label in axes[1].get_xticklabels():
    label.set_horizontalalignment("right")

plt.tight_layout()
plt.savefig(out_png, dpi=300)
PY
}





##5. Histogram of read coverage across IG loci for each sample to identify SV from phased and assembled contigs from IGentotyper


make_SV_bed_file () {
    cat <<EOF > "${data}/SV_ig_loci.bed"
igh	157688	167669	IGHV7-4-1.CH17
igh	210158	257471	IGHV5-10-1.CH17
igh	394658	426627	IGHV3-23.region.ABC9
igh	484922	559858	IGHV3-30.region.ABC11
igh	609413	682906	IGHV4-38-2.region.mixFosmids
igh	955770	1033787	IGHV1-69.region.CH17
igh	1149129	1194129	IGHV1-8.GRCh37
EOF
}


phased_histogram_of_read_coverage_across_IG_loci_for_each_sample () {
    module load samtools

    for SAMPLE_DIR in "${scratch}/igenotyper_run"/*; do
        [ -d "$SAMPLE_DIR" ] || continue
        SAMPLE=$(basename "$SAMPLE_DIR")

        BAM=$(find "$SAMPLE_DIR/alignments/" -name "contigs_to_ref_phased.sorted.bam" 2>/dev/null | head -1)
        [ -z "$BAM" ] && continue

        # Extract Read Groups directly into an array to avoid subshell complexities
        mapfile -t READ_GROUPS < <(samtools view -H "$BAM" | grep '^@RG' | sed 's/.*ID:\([^\t]*\).*/\1/')

        for RG in "${READ_GROUPS[@]}"; do
            [ -z "$RG" ] && continue

            while read -r chr start end name; do
                [ -z "$chr" ] && continue
                igh_region="${chr}:${start}-${end}"
                
                # Stream the specific RG and region directly into samtools coverage without writing to disk
                samtools view -b -r "$RG" "$BAM" "$igh_region" | \
                samtools coverage -A -D --histogram - > "${SAMPLE_DIR}/${name}.${RG}_coverage.txt"
                # Instead of samtools coverage, extract raw per-base depth
                samtools depth -r "$igh_region" "$RG_BAM" > "${SAMPLE_DIR}/${name}.${RG}_depth.txt"
            done < "${data}/SV_ig_loci.bed"
        done
    done

    module purge
}


plot_depth_histogram_for_each_sample () {
    module load python

    python - "$1" "${2}" "${3}" <<'PY'


    # This function will read the depth files generated in the previous step and create histograms of read depth across the IG loci for each sample and RG, which can help identify potential SVs based on coverage patterns.
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Usage: python plot_depth.py input_depth.txt output_image.png "Region Name"
depth_file = sys.argv[1]
output_img = sys.argv[2]
region_title = sys.argv[3]

# Read samtools depth output (Chromosome, Position, Depth)
df = pd.read_csv(depth_file, sep="\t", names=["Chr", "Pos", "Depth"])

plt.figure(figsize=(12, 4))
sns.lineplot(data=df, x="Pos", y="Depth", color="crimson", linewidth=1.5)

# Fill under the curve to easily spot homozygous deletions (drops to 0) or duplications (doubling of depth)
plt.fill_between(df["Pos"], df["Depth"], color="crimson", alpha=0.2)

plt.title(f"Read Depth Across {region_title}")
plt.xlabel("Genomic Position")
plt.ylabel("Depth of Coverage")
plt.tight_layout()

plt.savefig(output_img, dpi=300)
plt.close()
PY
}

#phased_histogram_of_read_coverage_across_IG_loci_for_each_sample
#plot_depth_histogram_for_each_sample 




phased_histogram_of_read_coverage_across_IG_loci_for_each_sample2 () {
    module load samtools
    module load python

    for SAMPLE_DIR in "${scratch}/igenotyper_run"/*; do
        [ -d "$SAMPLE_DIR" ] || continue
        SAMPLE=$(basename "$SAMPLE_DIR")

        BAM=$(find "$SAMPLE_DIR/alignments/" -name "contigs_to_ref_phased.sorted.bam" 2>/dev/null | head -1)
        [ -z "$BAM" ] && continue

        # Safely extract Read Groups into an array
        mapfile -t READ_GROUPS < <(samtools view -H "$BAM" | grep '^@RG' | sed 's/.*ID:\([^\t]*\).*/\1/')

        for RG in "${READ_GROUPS[@]}"; do
            [ -z "$RG" ] && continue

            # Reset standard input descriptor so the internal python heredoc doesn't swallow the BED file lines
            while read -r chr start end name <&3; do
                [ -z "$chr" ] && continue
                igh_region="${chr}:${start}-${end}"
                
                DEPTH_TXT="${SAMPLE_DIR}/${name}.${RG}_depth.txt"
                PLOT_IMG="${SAMPLE_DIR}/${name}.${RG}_SV_coverage.png"

                # Stream specific RG + region from the master BAM straight into depth
                samtools view -b -r "$RG" "$BAM" "$igh_region" | samtools depth - > "$DEPTH_TXT"

                # Skip plotting if no reads map to this haplotype/region (prevents python pandas errors)
                if [ ! -s "$DEPTH_TXT" ]; then
                    rm -f "$DEPTH_TXT"
                    continue
                fi

                # Inline Python execution (matching your structural style)
                python - "$DEPTH_TXT" "$PLOT_IMG" "${SAMPLE} - ${name} (Hap/RG: ${RG})" <<'PY'
import sys
import pandas as pd
import matplotlib.pyplot as plt

depth_file = sys.argv[1]
output_img = sys.argv[2]
title_str = sys.argv[3]

# Read the streamed per-base coverage data (samtools depth outputs: Chr, Pos, Depth)
# Since we used a pipe stream, chromosome names might match or be standard input placeholders
df = pd.read_csv(depth_file, sep="\t", names=["Chr", "Pos", "Depth"])

# Drop rows with missing values just in case
df.dropna(subset=["Pos", "Depth"], inplace=True)

if df.empty:
    sys.exit(0)

# Build the structural variant visualization plot
plt.figure(figsize=(14, 4.5))
plt.plot(df["Pos"], df["Depth"], color="#1f77b4", linewidth=1.2, label="Read Depth")
plt.fill_between(df["Pos"], df["Depth"], color="#1f77b4", alpha=0.15)

plt.title(title_str, fontsize=12, fontweight="bold")
plt.xlabel("Genomic Position (bp)", fontsize=10)
plt.ylabel("Depth of Coverage", fontsize=10)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.xlim(df["Pos"].min(), df["Pos"].max())
plt.ylim(bottom=0)

plt.tight_layout()
plt.savefig(output_img, dpi=200)
plt.close()
PY

                # Clean up per-base files to protect storage space
                rm -f "$DEPTH_TXT"

            done 3< "${data}/SV_ig_loci.bed" # Using File Descriptor 3 keeps loops stable
        done
    done

    module purge
}








##function calls 

#align_bam_files_to_reference_for_deepvariant_input
#make_bed_file
#call_variants_w_deepvariant
#merge_deepvariant_vcf_files

#index_bam_files
#create_and_copy_reference_sa_index 
#phase_bam_files_w_igenotyper
#phase_bam_files_w_igenotyper_repeat
#assemble_bam_files_w_igenotyper 
#IG_detect_w_igenotyper
#IG_detect_w_igenotyper_wo_2_samples


#get_vcf_stats_of_deepvariant_output
#summarize_merged_deepvariant_vcf_split "${results}/deep_variants_merged/deepvariant.cohort.vcf.gz"  "${scratch}/merged_vcf_summary_split.tsv"
#plot_merged_deepvariant_vcf_summary "${scratch}/merged_vcf_summary_split.tsv" "${scratch}/merged_vcf_summary_split.png"


#count_igenotyper_snps_per_sample_for_IG_loci
#plot_igenotyper_ig_locus_summary

#make_SV_bed_file
#phased_histogram_of_read_coverage_across_IG_loci_for_each_sample

#phased_histogram_of_read_coverage_across_IG_loci_for_each_sample2



