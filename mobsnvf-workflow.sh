#!/usr/bin/env bash

set -euo pipefail

## User-provided inputs | Change as necessary
sample_id=""
damage_type="ffpe"
fp_cut="0.5"
use_phi="true"

## File paths

base_dir="."
out_dir=""
ref_path=""
bam_path=""
vcf_path=""



# Process command-line arguments for phi and sample id
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--i|--id|--sample-id)
            sample_id="$2"
            shift 2
            ;;
        -b|--b|--bam|--bam-path)
            bam_path="$2"
            shift 2
            ;;
        -v|--v|--vcf|--vcf-path)
            vcf_path="$2"
            shift 2
            ;;
        -r|--r|--ref|--ref-path|--reference-genome-path)
            ref_path="$2"
            shift 2
            ;;
        -o|--o|-out|--out-dir)
            use_phi="$2"
            shift 2
            ;;
        --use-phi)
            use_phi="$2"
            shift 2
            ;;
        --damage-type|--dt|-d)
            damage_type="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            exit 1
            ;;
    esac
done




## Set up directories
#### Don't use absolute paths
base_dir="."

## Assign outdir
if [[ -z "$out_dir" ]]; then
  out_dir="${base_dir}/result"
fi

## Assign bam
if [[ -z "$bam_path" ]]; then
  echo "BAM file is required for mobsnvf"
  exit 1
elif [[ "$bam_path" == *"/"* ]]; then
  bam="$bam_path"
else
  bam="${base_dir}/bam/${sample_id}.bam"
fi

## Check if BAM file exists
if [[ ! -f "$bam" ]]; then
  echo "$bam does not exist"
  exit 1
fi

## Check and assign bam index
if [[ -f "${bam%.*}.bai" ]]; then
  bai="${bam%.*}.bai"
elif [[ -f "${bam}.bai" ]]; then
  bai="${bam}.bai"
else
  echo "$bam is not indexed. Please index the bam file"
  echo -e "Example: \n samtools index $bam"
  exit 1
fi

## Assign vcf
if [[ -z "$vcf_path" ]]; then
  echo "VCF file is required for mobsnvf"
  exit 1
elif [[ "$vcf_path" == *"/"* ]]; then
  vcf="$vcf_path"
else
  vcf="${base_dir}/vcf/${vcf_path}"
fi

## Check if vcf exists
if [[ ! -f "$vcf" ]]; then
  echo "$vcf does not exist"
  exit 1
fi

## Check and assign vcf index
if [[ "$vcf" =~ \.vcf\.gz$ || "$vcf" =~ \.bcf$ ]]; then
  idx_tbi=${vcf}.tbi
  idx_csi=${vcf}.csi
fi

if [[ -f $idx_tbi ]]; then
  vcf_idx="$idx_tbi"
elif [[ -f "$idx_csi" ]]; then
  vcf_idx="$idx_csi"
else
  echo "$vcf is is compressed but not indexed. Please index your vcf to proceed"
  echo -e "Example: \n bcftools index -t $vcf"
  exit 1
fi

# Assign reference genome
if [[ -z "$ref_path" ]]; then
  echo "Reference genome path is required for mobsnvf"
  exit 1
elif [[ "$ref_path" == *"/"* ]]; then
  ref="$ref_path"
else
  ref="${base_dir}/ref/${ref_path}"
fi

## Assign Sample ID
if [[ -z "$sample_id" ]]; then
  sample_id="${vcf}"
  if [[ $vcf =~ \.vcf\.gz$ ]]; then
    sample_id="${sample_id%.*}"
  fi
  sample_id=$(basename "${sample_id%.*}")

fi

## Assign final out path based on the sample ID and if Phi is used or not
if [[ "$use_phi" == "true" ]]; then
  out_dir="${out_dir}/${sample_id}/known_phi"
else
  out_dir="${out_dir}/${sample_id}/unknown_phi"
fi

mkdir -p "${out_dir}"

echo -e "Using Parameters:"
echo -e "Sample ID: $sample_id"
echo -e "Use phi: $use_phi\n"

rscript="${base_dir}/R/fdr-failed.R"

(
echo -e "Inputs:"
echo -e "1) Sample_ID     = ${sample_id}"
echo -e "2) BAM           = ${bam}"
echo -e "3) BAM_Index     = ${bai}"
echo -e "4) REF           = ${ref}"
echo -e "5) VCF           = ${vcf}"
echo -e "6) VCF_index     = ${vcf_idx}\n\n"


## 1) Task: bam_phi_estimation
if [[ "$use_phi" == "true" ]]; then
  echo -e "==== 1) Estimating phi from BAM using hts-mobsnvf quantify ===="
# Use hts-mobsnvf to estimates the phi parameter based on the BAM file.
  hts-mobsnvf quantify \
    -M freq \
    -t "${damage_type}" \
    -f "${ref}" \
    -b "${bam}" \
    -J "${out_dir}/${damage_type}_obquant.json"

  phi_json="${out_dir}/${damage_type}_obquant.json"
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
  -o "${out_dir}/${sample_id}_${damage_type}.snv" \
  > "${out_dir}/${sample_id}_snv_${damage_type}_mobsnvf.tsv"

annotated_snv="${out_dir}/${sample_id}_${damage_type}.snv"
echo -e "Annotated SNV file: ${annotated_snv}\n"


## 3) Task: snv_fdr_filter
echo -e "==== 3) FDR filter on SNVs using R script ===="
# Applies an FDR-based filter to the SNV file.
Rscript "${rscript}" \
  "${annotated_snv}" \
  -o "${out_dir}/${sample_id}_failed.snv" \
  --fp-cut "${fp_cut}"

failed_snv="${out_dir}/${sample_id}_failed.snv"
echo -e "Failed SNV file: ${failed_snv}\n"


## 4) Task: vcf_to_header
echo -e "==== 4) Extracting the header from the original VCF ===="
if [[ $vcf =~ \.vcf\.gz$ ]]; then
  zcat "${vcf}" | grep '^#' | sed 's/INFO\t.*/INFO/' > "${out_dir}/${sample_id}.vcf.header"
elif [[ $vcf =~ \.bcf$ ]]; then
  bcftools view -h "${vcf}" | sed 's/INFO\t.*/INFO/' > "${out_dir}/${sample_id}.vcf.header"
else
  grep '^#' "${vcf}" | sed 's/INFO\t.*/INFO/' > "${out_dir}/${sample_id}.vcf.header"
fi

vcf_header="${out_dir}/${sample_id}.vcf.header"
echo -e "Header file: ${vcf_header}\n"

## 5) Task: snv_to_vcf
echo -e "==== 5) Converting SNVs to VCF ===="
out_name="${out_dir}/${sample_id}_artifacts"
cat "${vcf_header}" > "${out_name}.vcf"

# Create a new VCF by combining info from the SNV data.
paste -d'\t' \
  <(cut -f 1-2 "${failed_snv}" | sed 's/$/\t./') \
  <(cut -f 3-4 "${failed_snv}" | sed 's/$/\t.\tartifact\t./') >> "${out_name}.vcf"

# Index the new VCF
gatk IndexFeatureFile \
  -I "${out_name}.vcf" \
  -O "${out_name}.vcf.idx"

artifacts_vcf="${out_name}.vcf"
artifacts_vcf_index="${out_name}.vcf.idx"
echo -e "Mask VCF (failed SNVs): ${artifacts_vcf}"
echo -e "Mask VCF index: ${artifacts_vcf_index}\n"


## 6) Task: vcf_mask_variants
echo -e "==== 6) Masking variants in the original VCF ===="
# Apply the mask 'artifacts_vcf' to the original VCF to identify differences in varaints
gatk VariantFiltration \
  -V "${vcf}" \
  --invalidate-previous-filters false \
  --mask "${artifacts_vcf}" \
  --mask-name "${damage_type}" \
  -O "${out_dir}/${sample_id}_masked.vcf"

masked_vcf="${out_dir}/${sample_id}_masked.vcf"
masked_vcf_index="${out_dir}/${sample_id}_masked.vcf.idx"
echo -e "Masked VCF: ${masked_vcf}"
echo -e "Masked VCF index: ${masked_vcf_index}\n"


## 7) Task: vcf_select_variants
echo -e "==== 7) Selecting unfiltered variants ===="

filtered_vcf="${out_dir}/${sample_id}_filtered.vcf"
filtered_vcf_index="${out_dir}/${sample_id}_filtered.vcf.idx"

gatk SelectVariants \
  -V "${masked_vcf}" \
  --exclude-filtered \
  -O "${filtered_vcf}"


gatk IndexFeatureFile \
  -I "${filtered_vcf}" \
  -O "${filtered_vcf_index}"

echo -e "Selected VCF (final): ${filtered_vcf}"
echo -e "Selected VCF index (final): ${filtered_vcf_index}\n"


## Remove Intermediate files
# Users can comment the lines if they want them to be kept
rm "${vcf_header}"
rm "${masked_vcf}"
rm "${masked_vcf_index}"
rm "${out_dir}/${sample_id}_snv_${damage_type}_mobsnvf.tsv"


## Final output summary
echo -e "-------------------"
echo -e "Filtering completed!"
echo -e "Outputs:"

n=1
if [[ "$use_phi" == "true" ]]; then
  echo -e "${n}) phi_json            = ${phi_json}"
  ((n++))
fi
echo -e "${n}) annotated_snv       = ${annotated_snv}"
((n++))
echo -e "${n}) artifacts_vcf         = ${artifacts_vcf}"
((n++))
echo -e "${n}) artifacts_vcf_index   = ${artifacts_vcf_index}"
((n++))
echo -e "${n}) selected_vcf        = ${filtered_vcf}"
((n++))
echo -e "${n}) selected_vcf_index  = ${filtered_vcf_index}"

## Log the output
) 2>&1 | tee "${out_dir}/${sample_id}.log"