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

index_bam_files () { 
    #sample download only came with .pbi files, so we need to index the bam files before we can use them for downstream analyses
    module load samtools

    for file in "${scratch}/aligned_bams"/*.bam ; do
        samtools index "$file"
    done

    module purge samtools
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
    gvcfs=$(cd "${scratch}/deep_variant_calls" && ls *.g.vcf.gz | sed 's|^|/data/|')

    # 2. Run Singularity normally
    singularity exec \
        -B "${scratch}/deep_variant_calls:/data,${data}:/bed" \
        "/hpc/users/willij115/glnexus_v1.2.2.sif" \
        glnexus_cli \
        --config DeepVariantWGS \
        --bed /bed/ig_loci.bed \
        $gvcfs \
    | bcftools view - | bgzip -c > "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"

    bcftools index -t "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"

    module purge
}



## vcf stats func
get_vcf_stats_of_deepvariant_output () {
    # here we can add code to analyze the output from deepvariant, such as calculating the number of variants called, the distribution of variant types, etc.
    module load bcftools

    bcftools stats "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz" > "${scratch}/deepvariant.cohort.vcf.stats.txt"


    module purge
}







##2 call variants with IGentotyper

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



merge_igenotyper_snv_vcf_files () {
    module load bcftools
    module load htslib
    set -euo pipefail

    local outdir="${scratch}/igenotyper_merged_sample_vcf"
    local workdir="${scratch}/igenotyper_merged_sample_vcf/tmp"
    mkdir -p "$outdir" "$workdir"

    local vcfs=()
    for file in "${scratch}/igenotyper_run"/*/variants/snvs_assembly.vcf; do
        sample_dir=$(basename "$(dirname "$(dirname "$file")")")
        sample_name=$(grep -oE 'bc[0-9]+' <<< "$sample_dir" | head -n1) || {
            echo "ERROR: could not parse sample name from $sample_dir" >&2
            return 1
        }


        local base="${workdir}/${sample_name}"
        printf "sample\t%s\n" "$sample_name" > "${base}.samples"

        # Step 1: rename sample, add contig lengths from .fai, in one reheader pass
        bcftools reheader \
            -s "${base}.samples" \
            -f "${data1}/reference.fasta.fai" \
            "$file" \
            -o "${base}.reheadered.vcf.gz"

        # Step 2: sort
        bcftools sort -Oz -o "${base}.sorted.vcf.gz" "${base}.reheadered.vcf.gz"
        bcftools index -c "${base}.sorted.vcf.gz"

        vcfs+=("${base}.sorted.vcf.gz")
    done

    bcftools merge \
        -Oz \
        -o "${outdir}/merged_igenotyper.vcf.gz" \
        "${vcfs[@]}"

    bcftools index -c "${outdir}/merged_igenotyper.vcf.gz"
}






### 4. analyze and compare igentoyper vs deepvariant outputs 

## ---------------------------------------------------------------------------
## Utility: (re)index a VCF with a CSI index, forcing a fresh index
## ---------------------------------------------------------------------------
index_vcf () {
    local vcf="$1"
    rm -f "${vcf}.tbi" "${vcf}.csi"
    bcftools index -f -c "$vcf"
}

## ---------------------------------------------------------------------------
## Extract a single locus region from a VCF (no multiallelic splitting).
## Writes: ${outdir}/${locus}.region.vcf.gz (+ .csi index)
## Prints the path to the region VCF on stdout.
## ---------------------------------------------------------------------------
extract_locus_region () {
    local vcf="$1"
    local locus="$2"
    local region="$3"
    local outdir="$4"

    local region_vcf="${outdir}/${locus}.region.vcf.gz"

    bcftools view -r "$region" -Oz -o "$region_vcf" "$vcf"
    bcftools index -t "$region_vcf"

    local n_sites
    n_sites=$(bcftools view -H "$region_vcf" | wc -l)
    if [[ "$n_sites" -eq 0 ]]; then
        echo "WARNING: 0 sites for $locus at $region — check contig naming in $vcf" >&2
    fi

    echo "$region_vcf"
}

## ---------------------------------------------------------------------------
## Split multiallelic sites (SNPs and indels) into separate biallelic records.
## Writes: ${outdir}/${locus}.split.vcf.gz (+ .csi index)
## Prints the path to the split VCF on stdout.
## ---------------------------------------------------------------------------
split_multiallelic_sites () {
    local region_vcf="$1"
    local locus="$2"
    local outdir="$3"

    local split_vcf="${outdir}/${locus}.split.vcf.gz"

    bcftools norm -m -any -Oz -o "$split_vcf" "$region_vcf"
    bcftools index -t "$split_vcf"

    echo "$split_vcf"
}

## ---------------------------------------------------------------------------
## Annotate a (split, biallelic) VCF with MAF and count sites with MAF > threshold.
## Writes: ${outdir}/${locus}.split.maf.vcf.gz (+ .csi index)
## Prints: "<maf_vcf_path>\t<maf_count>" on stdout
## ---------------------------------------------------------------------------
compute_locus_maf_count () {
    local split_vcf="$1"
    local locus="$2"
    local outdir="$3"
    local maf_threshold="${4:-0.4}"

    local maf_vcf="${outdir}/${locus}.split.maf.vcf.gz"

    bcftools +fill-tags "$split_vcf" -Oz -o "$maf_vcf" -- -t MAF
    bcftools index -t "$maf_vcf"

    local maf_count
    maf_count=$(bcftools view -i "MAF>${maf_threshold}" -H "$maf_vcf" | wc -l | awk '{print $1}')

    printf "%s\t%s\n" "$maf_vcf" "$maf_count"
}

## ---------------------------------------------------------------------------
## Count per-sample total heterozygous genotype calls from the RAW
## (pre-multiallelic-split) region VCF for one locus. This must run on the
## original genotype calls, since splitting a multiallelic site (e.g. 1/2)
## into biallelic records changes how each resulting genotype looks and would
## distort a het count taken afterwards.
## Writes ${outdir}/${locus}.het_counts.tsv with columns: sample, het_count.
## Prints that file's path on stdout.
## ---------------------------------------------------------------------------
count_heterozygous_per_sample () {
    local region_vcf="$1"
    local locus="$2"
    local outdir="$3"

    local samples_file="${outdir}/${locus}.samples.txt"
    local het_file="${outdir}/${locus}.het_counts.tsv"
    bcftools query -l "$region_vcf" > "$samples_file"

    bcftools query -f '[%GT\t]\n' "$region_vcf" | \
    awk -v samples_file="$samples_file" '
      # true if genotype is heterozygous: at least 2 distinct non-missing alleles
      # (handles multiallelic genotypes like 1/2, 0/2, etc. correctly, since
      # this runs BEFORE any multiallelic splitting)
      function ishet(gt,   a,n,i,seen,distinct,allele) {
        if (gt == "./." || gt == ".|." || gt == "." || gt == "") return 0
        n = split(gt, a, /[\/|]/)
        distinct = 0
        delete seen
        for (i = 1; i <= n; i++) {
          allele = a[i]
          if (allele == ".") continue
          if (!(allele in seen)) { seen[allele] = 1; distinct++ }
        }
        return (distinct >= 2) ? 1 : 0
      }
      BEGIN {
        while ((getline s < samples_file) > 0) samples[++n] = s
        close(samples_file)
      }
      {
        for (i = 1; i <= NF; i++) {
          if ($i == "") continue
          sample = samples[i]
          if (ishet($i)) hetcnt[sample]++
        }
      }
      END {
        for (i = 1; i <= n; i++) {
          s = samples[i]
          print s "\t" hetcnt[s]+0
        }
      }' > "$het_file"

    echo "$het_file"
}

## ---------------------------------------------------------------------------
## Count per-sample SNP / INDEL from a (split, biallelic) VCF for one locus.
## Heterozygosity is NOT computed here — it is computed separately on the
## pre-split VCF by count_heterozygous_per_sample and merged in afterwards.
## Appends rows to $out (locus, sample, SNP, INDEL, het_count,
## locus_SNP_sites_MAF_gt_<threshold>).
## ---------------------------------------------------------------------------
count_sample_variants_per_locus () {
    local maf_vcf="$1"
    local locus="$2"
    local maf_count="$3"
    local het_file="$4"
    local outdir="$5"
    local out="$6"

    local samples_file="${outdir}/${locus}.samples.txt"
    bcftools query -l "$maf_vcf" > "$samples_file"

    bcftools query -f '%REF\t%ALT[\t%GT]\n' "$maf_vcf" | \
    awk -v locus="$locus" -v maf_count="$maf_count" \
        -v samples_file="$samples_file" -v het_file="$het_file" '
      # true if genotype carries at least one non-ref, non-missing allele
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
        while ((getline line < het_file) > 0) {
          split(line, parts, "\t")
          het[parts[1]] = parts[2]
        }
        close(het_file)
      }
      {
        type = vartype($1, $2)
        for (i = 3; i <= NF; i++) {
          sample = samples[i-2]
          if (nonref($i)) cnt[sample, type]++
        }
      }
      END {
        for (i = 1; i <= n; i++) {
          s = samples[i]
          print locus "\t" s "\t" \
            cnt[s, "SNP"]+0 "\t" cnt[s, "INDEL"]+0 "\t" \
            het[s]+0 "\t" \
            maf_count+0
        }
      }' >> "$out"
}

## ---------------------------------------------------------------------------
## Orchestrator: indexes the input VCF, then for each Ig locus:
##   1. extracts the region (raw, pre-split)
##   2. counts per-sample heterozygous genotypes on the RAW region VCF
##   3. splits multiallelic sites into biallelic records
##   4. computes per-locus MAF count on the split VCF
##   5. counts per-sample SNP/INDEL on the split VCF
## and writes one combined TSV to $out.
## ---------------------------------------------------------------------------
summarize_merged_deepvariant_vcf_split () {
    module load bcftools
    set -euo pipefail
    local vcf="$1"
    local out="${2:-summary.tsv}"
    local maf_threshold="${3:-0.4}"
    local tmpdir
    tmpdir="$(mktemp -d)"
    trap 'rm -rf "$tmpdir"' RETURN

    local IGH_CHR="igh";   local IGH_START="1";        local IGH_END="1193129"
    local IGK_CHR="chr2";  local IGK_START="88807162"; local IGK_END="90310090"
    local IGL_CHR="chr22"; local IGL_START="22414484"; local IGL_END="23392002"

    index_vcf "$vcf"

    echo "Working in $tmpdir" >&2
    echo -e "locus\tsample\tSNP\tINDEL\theterozygous_count\tlocus_SNP_sites_MAF_gt_${maf_threshold}" > "$out"

    for locus in IGH IGK IGL; do
        local chr start end region region_vcf het_file split_vcf maf_result maf_vcf maf_count
        case "$locus" in
            IGH) chr="$IGH_CHR"; start="$IGH_START"; end="$IGH_END" ;;
            IGK) chr="$IGK_CHR"; start="$IGK_START"; end="$IGK_END" ;;
            IGL) chr="$IGL_CHR"; start="$IGL_START"; end="$IGL_END" ;;
        esac
        region="${chr}:${start}-${end}"

        # 1. extract region (raw, pre-split)
        region_vcf=$(extract_locus_region "$vcf" "$locus" "$region" "$tmpdir")

        # 2. heterozygous count on RAW genotypes, before multiallelic split
        het_file=$(count_heterozygous_per_sample "$region_vcf" "$locus" "$tmpdir")

        # 3. split multiallelic sites
        split_vcf=$(split_multiallelic_sites "$region_vcf" "$locus" "$tmpdir")

        # 4. per-locus MAF count on the split VCF
        maf_result=$(compute_locus_maf_count "$split_vcf" "$locus" "$tmpdir" "$maf_threshold")
        maf_vcf=$(cut -f1 <<< "$maf_result")
        maf_count=$(cut -f2 <<< "$maf_result")

        # 5. per-sample SNP/INDEL count on the split VCF, merged with the het count
        count_sample_variants_per_locus "$maf_vcf" "$locus" "$maf_count" "$het_file" "$tmpdir" "$out"
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
import matplotlib.gridspec as gridspec
import numpy as np

summary_tsv = sys.argv[1]
out_png = sys.argv[2]

df = pd.read_csv(summary_tsv, sep="\t", comment="#")
df = df[df["sample"].notna()].copy()

metrics = ["SNP", "INDEL", "heterozygous_count"]
maf_col = next((c for c in df.columns if c.startswith("locus_SNP_sites_MAF")), None)

loci = [x for x in ["IGH", "IGK", "IGL"] if x in set(df["locus"])]
if not loci:
    loci = sorted(df["locus"].dropna().unique())

sample_order = sorted(df["sample"].unique())
colors = {"SNP": "#4C72B0", "INDEL": "#DD8452", "heterozygous_count": "#55A868"}
width = 0.25
x = np.arange(len(sample_order))

fig = plt.figure(figsize=(5.5 * len(loci), 9))
gs = gridspec.GridSpec(2, len(loci), height_ratios=[2, 1], figure=fig)

for col, locus in enumerate(loci):
    ax = fig.add_subplot(gs[0, col])
    sub = df[df["locus"] == locus].set_index("sample").reindex(sample_order)

    for i, metric in enumerate(metrics):
        ax.bar(x + (i - 1) * width, sub[metric].fillna(0), width,
               label=metric, color=colors[metric])

    ax.set_title(locus, fontsize=13, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=45, ha="right")
    ax.set_ylabel("Count" if col == 0 else "")
    if col == 0:
        ax.legend(loc="upper right", fontsize=9)

    if maf_col is not None and sub[maf_col].notna().any():
        maf_val = sub[maf_col].dropna().iloc[0]
        ax.text(0.98, 0.02, f"MAF>0.4 sites: {int(maf_val)}",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=9, color="dimgray",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="lightgray", alpha=0.8))

ax_bottom = fig.add_subplot(gs[1, :])
totals = df.groupby("locus")[metrics].sum().reindex(loci)
xt = np.arange(len(loci))
for i, metric in enumerate(metrics):
    ax_bottom.bar(xt + (i - 1) * width, totals[metric], width,
                  label=metric, color=colors[metric])
ax_bottom.set_xticks(xt)
ax_bottom.set_xticklabels(loci)
ax_bottom.set_ylabel("Total count\n(summed across samples)")
ax_bottom.set_title("Per-locus totals", fontsize=12, fontweight="bold")
ax_bottom.legend(loc="upper right", fontsize=9)

plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved plot to {out_png}")
PY
}




















































##function calls 

###deep variant calls
#align_bam_files_to_reference_for_deepvariant_input
#index_bam_files
#make_bed_file
#call_variants_w_deepvariant
#merge_deepvariant_vcf_files
#get_vcf_stats_of_deepvariant_output 

###igenotyper calls
#create_and_copy_reference_sa_index 
#phase_bam_files_w_igenotyper
#phase_bam_files_w_igenotyper_repeat
#assemble_bam_files_w_igenotyper 
#IG_detect_w_igenotyper
#IG_detect_w_igenotyper_wo_2_samples
#merge_igenotyper_snv_vcf_files


### analysis and comparison of deepvariant vs igentotyper plots
#summarize_merged_deepvariant_vcf_split "${results}/deepvariants_merged_new_bam/deepvariant.cohort.vcf.gz"  "${scratch}/merged_deepvariant_summary_split_new_bam.tsv"
#summarize_merged_deepvariant_vcf_split "${scratch}/igenotyper_merged_sample_vcf/merged_igenotyper.vcf.gz"  "${scratch}/merged_igenotyper_summary_split.tsv"
#plot_merged_deepvariant_vcf_summary "${scratch}/merged_deepvariant_summary_split_new_bam.tsv" "${scratch}/deepvariant_summary_counts.png"
#plot_merged_deepvariant_vcf_summary "${scratch}/merged_igenotyper_summary_split.tsv" "${scratch}/igenotyper_summary_counts.png"    



