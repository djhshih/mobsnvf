#!/usr/bin/env bash

set -eo pipefail

## User-provided inputs | Change as necessary
sample_id=""
damage_type="ffpe"
fp_cut="1e-08"
use_phi="true"

## Prefix and suffix for VCF and BAM files
pvcf=""
svcf=""
pbam=""
sbam=""
ref_name="Homo_sapiens_assembly38.fasta"

# Process command-line arguments for phi and sample id
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--i|--id|--sample-id)
            sample_id="$2"
            shift 2
            ;;
        --use-phi)
            use_phi="$2"
            shift 2
            ;;
        --prefix-vcf|--pvcf|-pv)
            pvcf="$2"
            shift 2
            ;;
        --suffix-vcf|--svcf|-sv)
            svcf="$2"
            shift 2
            ;;
        --prefix-bam|--pbam|-pb)
            pbam="$2"
            shift 2
            ;;
        --suffix-bam|--sbam|-sb)
            sbam="$2"
            shift 2
            ;;
        --damage-type|--dt|-d)
            damage_type="$2"
            shift 2
            ;;
        --reference-genome|--ref|-r)
            ref_name="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1
            ;;
    esac
done

# Verify that required parameters are provided
if [[ -z "$sample_id" ]]; then
    echo "Error: sample_id is required. Use --i <sample_id>"
    exit 1
fi

echo -e "Using Parameters:"
echo -e "Sample ID: $sample_id"
echo -e "Use phi: $use_phi\n"

## Set up directories
#### Don't use absolute paths
base_dir="."

if [[ "$use_phi" == "true" ]]; then
  out_dir="${base_dir}/results/${sample_id}/known_phi/"
else
  out_dir="${base_dir}/results/${sample_id}/unknown_phi/"
fi
mkdir -p "${out_dir}"

bam="${base_dir}/bam/${pbam}${sample_id}${sbam}.bam"
bai="${base_dir}/bam/${pbam}${sample_id}${sbam}.bai"
ref="${base_dir}/ref/${ref_name}"
vcf="${base_dir}/vcf/${pvcf}${sample_id}${svcf}.vcf"
rscript="${base_dir}/R/fdr-failed.R"

(
echo -e "Inputs:"
echo -e "1) Sample_ID     = ${sample_id}"
echo -e "2) BAM           = ${bam}"
echo -e "3) BAM_Index     = ${bai}"
echo -e "4) REF           = ${ref}"
echo -e "5) VCF           = ${vcf}"
echo -e "6) VCF_index     = ${vcf}.idx\n\n"


## 1) Task: bam_phi_estimation
if [[ "$use_phi" == "true" ]]; then
  echo -e "==== 1) Estimating phi from BAM using hts-mobsnvf quantify ===="
# Use hts-mobsnvf to estimates the phi parameter based on the BAM file.
  hts-mobsnvf quantify \
    -M freq \
    -t "${damage_type}" \
    -f "${ref}" \
    -b "${bam}" \
    -J "${out_dir}${damage_type}_obquant.json"

  phi_json="${out_dir}${damage_type}_obquant.json"
  echo -e "phi_json created at: ${phi_json}\n"
else
  echo -e "==== 1) Continuing without phi estimation ===="
fi


## 2) Task: snv_mobsnvf_filter
echo -e "==== 2) Filtering SNVs using hts-mobsnvf identify ===="
# The phi value is parsed from the JSON.
if [[ "$use_phi" == "true" ]]; then
  phi_n=$(grep -v 'bam_file' "${phi_json}" | grep "\"phi\"" | sed -E 's/.*"phi":.([0-9.e+\-]+),?/\1/')
  echo -e "Parsed phi value: $phi_n"
else
  phi_n="0.0"
fi

# Phi opt flag for hts-mobsnvf identify is generated. If phi is zero, it is not used for downstream damage identificaton.
if [[ "$phi_n" =~ 0\.0+([eE]00)? ]]; then
  echo -e "Estimated phi is $phi_n$; using unknown phi"
  phi_opts=""
else
  echo -e "Estimated phi is $phi_n; using fixed phi"
  phi_opts="--phi ${phi_n} --fixed-phi"
fi

# Perform variant identification:
## -g is the p-value cutoff. By setting it to zero we are disabling the p-value filter.
hts-mobsnvf identify \
  -M freq \
  -t "${damage_type}" \
  -b "${bam}" \
  -V "${vcf}" \
  -g 0 \
  ${phi_opts} \
  -o "${out_dir}${sample_id}_${damage_type}.snv" \
  > "${out_dir}${sample_id}_snv_${damage_type}_mobsnvf.tsv"

annotated_snv="${out_dir}${sample_id}_${damage_type}.snv"
echo -e "Annotated SNV file: ${annotated_snv}\n"


## 3) Task: snv_fdr_filter
echo -e "==== 3) FDR filter on SNVs using R script ===="
# Applies an FDR-based filter to the SNV file.
Rscript "${rscript}" \
  "${annotated_snv}" \
  -o "${out_dir}${sample_id}_failed.snv" \
  --fp-cut "${fp_cut}"

failed_snv="${out_dir}${sample_id}_failed.snv"
echo -e "Failed SNV file: ${failed_snv}\n"


## 4) Task: vcf_to_header
echo -e "==== 4) Extracting the header from the original VCF ===="
grep '^#' "${vcf}" | sed 's/INFO\t.*/INFO/' > "${out_dir}${sample_id}.vcf.header"
vcf_header="${out_dir}${sample_id}.vcf.header"
echo -e "Header file: ${vcf_header}\n"


## 5) Task: snv_to_vcf
echo -e "==== 5) Converting SNVs to VCF ===="
out_name="${out_dir}${sample_id}_removed"
cat "${vcf_header}" > "${out_name}.vcf"

# Create a new VCF by combining info from the SNV data.
paste -d'\t' \
  <(cut -f 1-2 "${failed_snv}" | sed 's/$/\t./') \
  <(cut -f 3-4 "${failed_snv}" | sed 's/$/\t.\tartifact\t./') >> "${out_name}.vcf"

# Index the new VCF
gatk IndexFeatureFile \
  -I "${out_name}.vcf" \
  -O "${out_name}.vcf.idx"

removed_vcf="${out_name}.vcf"
removed_vcf_index="${out_name}.vcf.idx"
echo -e "Mask VCF (failed SNVs): ${removed_vcf}"
echo -e "Mask VCF index: ${removed_vcf_index}\n"


## 6) Task: vcf_mask_variants
echo -e "==== 6) Masking variants in the original VCF ===="
# Apply the mask 'removed_vcf' to the original VCF to identify differences in varaints
gatk VariantFiltration \
  -V "${vcf}" \
  --invalidate-previous-filters false \
  --mask "${removed_vcf}" \
  --mask-name "${damage_type}" \
  -O "${out_dir}${sample_id}_masked.vcf"

masked_vcf="${out_dir}${sample_id}_masked.vcf"
masked_vcf_index="${out_dir}${sample_id}_masked.vcf.idx"
echo -e "Masked VCF: ${masked_vcf}"
echo -e "Masked VCF index: ${masked_vcf_index}\n"


## 7) Task: vcf_select_variants
echo -e "==== 7) Selecting unfiltered variants ===="
gatk SelectVariants \
  -V "${masked_vcf}" \
  --exclude-filtered \
  -O "${out_dir}${sample_id}_selected.vcf"

selected_vcf="${out_dir}${sample_id}_selected.vcf"
selected_vcf_index="${out_dir}${sample_id}_selected.vcf.idx"

echo -e "Selected VCF (final): ${selected_vcf}"
echo -e "Selected VCF index (final): ${selected_vcf_index}\n"


## Remove Intermediate files
# Users can comment the lines if they want them to be kept
rm "${vcf_header}"
rm "${masked_vcf}"
rm "${masked_vcf_index}"
rm "${out_dir}${sample_id}_snv_${damage_type}_mobsnvf.tsv"


## Final output summary
echo -e "-------------------"
echo -e "Filtering completed!"
echo -e "Outputs:"
echo -e "1) phi_json            = ${phi_json}"
echo -e "2) annotated_snv       = ${annotated_snv}"
echo -e "3) removed_vcf         = ${removed_vcf}"
echo -e "4) removed_vcf_index   = ${removed_vcf_index}"
echo -e "5) selected_vcf        = ${selected_vcf}"
echo -e "6) selected_vcf_index  = ${selected_vcf_index}"

## Log the output
) 2>&1 | tee "${out_dir}${sample_id}.log"
