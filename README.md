# HTSPAN MOBSNVF Filter

The mobsnvf filter is a module of htspan for removing artifacts from high throughput sequencing data. Supported damage types are:

*   FFPE
*   OXOG

This repository provides a pipeline for using the mobsnvf module of htspan for a production workflow.

## Dependencies

*   htspan (see **[djhshih/htspan](https://github.com/djhshih/setup-linux)** for instructions on setting up htspan)
*   R 4.4.3
*   Python 3.12.10
*   GATK 4.6.2.0
*   libcurl-dev 8.14.1
*   liblzma-dev 5.4.5
*   bzip2 1.22
*   cmake 4.0.3
*   g++ 13.3.0
*   gcc 13.3.0
*   zlib1g-dev
*   libbz2-dev

### R libraries:
*   `stringr`
*   `argparser`

### Python libraries:
*   polars
*   pysam
*   numpy
*   seaborn
*   matplotlib

### Optional dependencies:
*   gsutil

## Setup Procedure

If your system does git installed, please setup and configure git first. You may take a look here for a **[demonstration](https://github.com/djhshih/setup-linux)**.

```bash
sudo apt install git
```

### Install htspan

Please refer to **[djhshih/htspan](https://github.com/djhshih/htspan)** for instructions on setting up **htspan** and a minimalistic use case.

### Clone mobsnvf repo

Then clone this `mobsnvf` repo:

```bash
git clone https://github.com/djhshih/mobsnvf.git
cd mobsnvf
```

## Installing dependencies

The dependencies for **htspan** and ways to install them are mentioned in it's respective repostory:

This pipeline also relies on R, R packages and GATK.

GATK can be installed according to the [instructions provided here](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4).

The R packages may be installed as follows:

```bash
sudo apt install r-base r-base-dev
R -e 'install.packages(c("argparser","stringr"), repos="https://cloud.r-project.org")'
```

### Alternate method

Alternatively, you may use the `Dockerfile` included with this repo to build and run a docker container. This solves all dependencies for compiling htspan and running the pipeline.

If docker is not installed on your system, follow instructions **[here]** to install docker or contact your system administrator if you don't have the privileges.

Then navigate to the cloned `mobsnvf/` repo and use the following commands to make a Docker image and run an interactive Docker container:

```bash
docker build -t mobsnvf .
docker run --rm -it -v .:/work -w /work mobsnvf
```

This will create and run an interactive Docker container, mounting the current directory into the `/work` directory inside the container and setting it as the working directory.

## How to use this pipeline

The workflow requires explicit file paths for inputs.

### Command-Line Arguments

| Argument | Alias(es) | Required? | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--bam-path` | `-b`, `--bam` | **Yes** | | Path to the input BAM file. Must be indexed (`.bai`). |
| `--vcf-path` | `-v`, `--vcf` | **Yes** | | Path to the input VCF file. Must be indexed (`.tbi` or `.csi` if gzipped). |
| `--ref-path` | `-r`, `--ref` | **Yes** | | Path to the reference genome FASTA file. Must be indexed (`.fai`). |
| `--sample-id` | `-i`, `--id` | No | Inferred from VCF | A unique identifier for the sample. If not provided, it will be automatically derived from the VCF filename. |
| `--use-phi` | | No | `true` | Set to `true` to perform phi estimation or `false` to skip it. |
| `--damage-type`| `-d`, `--damage` | No | `ffpe` | The type of DNA damage to model. Supported values are `ffpe` and `oxog`. |
| `--out-dir` | `-o`, `--out` | No | `./result` | The base directory where results will be saved. |
| `--fp-cut` | | No | `0.5` | False positive cutoff value for the FDR filter step. |

### Example Usage

Below is a typical command to run the pipeline for a single sample. You must provide the paths to your BAM, VCF, and reference genome.

**Note:** It is important that the reference fasta provided is the same as the one used for sequence alignment. 

```bash
bash mobsnvf-workflow.sh \
    --sample-id "sample_name" \
    --bam-path "/path/to/data/bams/sample01.bam" \
    --vcf-path "/path/to/data/vcfs/sample01.vcf.gz" \
    --ref-path "/path/to/genomes/bam.specific.fasta" \
    --damage-type "ffpe" \
    --use-phi "true" \
    --out-dir "/path/to/output"
```

A more barebones command without optional parameters would look like this:

```bash
bash mobsnvf-workflow.sh \
    --bam-path "/path/to/data/bams/sample01.bam" \
    --vcf-path "/path/to/data/vcfs/sample01.vcf.gz" \
    --ref-path "/path/to/genomes/bam.specific.fasta"
```

In this case the sample id will be inferred from the VCF filename, and the output will be saved to `./result/<sample_id>/known_phi/`.

### Reference Genome

If you are sure that your BAM uses a standard hg38 reference genome, you can use the `get.sh` script included in the ref directory to download it:

**Note:** This requires gsutil to be installed on your system. If you don't have gsutil, you can install it by following the instructions [here](https://cloud.google.com/storage/docs/gsutil_install).

```bash
cd ref
bash get.sh
```

### Testing the Pipeline with Example Data

This repository includes an `example-data` directory containing a simple ffpe BAM and VCF file along with a matched normal, allowing you to perform a test run to ensure the pipeline is configured correctly.

You can run the test using the following command from the root of the repository. **Remember to replace `/path/to/your/GRCh38.d1.vd1.fa`** with the actual path to your reference genome file.

```bash
bash mobsnvf-workflow.sh \
    --bam-path "example-data/bam/TCGA-A6-6650-01B-02D-A270-10_Illumina_gdc_realn_compressed.bam" \
    --vcf-path "example-data/vcf/a82846b3-c3df-443e-b9e6-836380fa60e3.vcf.gz" \
    --ref-path "example-data/vcf/tp53_hg38.fasta"
```

Upon successful completion, you will find the results in the `example-output/TCGA-A6-6650-test/known_phi/` directory.

### Phi Parameter

The parameter **phi** is an estimate of the extent of artificial DNA damage.

Setting the `--use-phi` flag to `true` (the default) will estimate phi with `hts-mobsnvf quantify` and use this single value when identifying artifacts across all sites with `hts-mobsnvf identify`. This is the recommended approach as it generally improves accuracy.

When `--use-phi` is set to `false`, the initial phi estimation step is skipped, and `hts-mobsnvf identify` will instead estimate the extent of DNA damage (phi) individually for each SNV site.

### Results

The results for the workflow will be saved to the directory specified by `--out-dir`, organized by sample ID and whether phi was used. For example: `/path/to/output/sample01/known_phi/`.

The primary output files are:

*   **`{sample_id}_{damage_type}.snv`**: A tab-separated file containing all SNVs annotated by `hts-mobsnvf identify`.
*   **`{sample_id}_artifacts.vcf`**: A VCF file containing the variants that were identified as artifacts.
*   **`{sample_id}_artifacts.vcf.idx`**: The index for the artifacts VCF.
*   **`{sample_id}_filtered.vcf`**: The final, filtered VCF file with artifacts removed. **This is the main result file.**
*   **`{sample_id}_filtered.vcf.idx`**: The index for the final filtered VCF.
*   **`{sample_id}.log`**: A log file containing the complete standard output of the workflow run.
*   **`{damage_type}_obquant.json`**: (Only if `--use-phi` is `true`) A JSON file containing the estimated phi value.

## Running multiple samples

To process multiple samples, we recommend creating a manifest file and looping through it.

1.  **Create a manifest file**

    Create a tab-separated values (TSV) file (e.g., `samples.tsv`) with columns for the sample ID and the absolute paths to the corresponding BAM and VCF files.

    **Example `samples.tsv`:**

    ```tsv
    /path/to/bams/sample01.bam	/path/to/vcfs/sample01.vcf.gz
    /path/to/bams/sample02.bam	/path/to/vcfs/sample02.vcf.gz
    /path/to/bams/sample03.bam	/path/to/vcfs/sample03.vcf.gz
    ```

2.  **Run the workflow using a loop**

    You can now use a simple `bash` loop to read the manifest file and create a script for each sample. Then you may run these scripts in batch mode.

    ```bash
    #!/bin/bash
    
    ref="/path/to/your/hg38.fasta"
    out="/path/to/your/results"
    
    mkdir -p job

    # Read the manifest and generate a command for each sample
    while IFS=$'\t' read -r  bam_path vcf_path; do
      # Skip header or empty lines
      [[ "$bam_path" == "bam_path" || -z "$bam_path" ]] && continue
    
      echo "bash mobsnvf-workflow.sh \
        --bam-path \"${bam_path}\" \
        --vcf-path \"${vcf_path}\" \
        --ref-path \"${ref}\" \
        --out-dir \"${out}\"" > "job/$(basename "${vcf_path%%.*}").sh"
    done < samples.tsv
    
    ```

    These scripts can then be executed in batch, for example using [dlazy](https://github.com/djhshih/dlazy/tree/main).

    ```bash
    dlazy job
    ```

## Evaluation

An evaluation python script is included with this repository to analyse the mutation profile, mutation counts, and mutation signature plots for the VCF files.

### Example Usage

To run the evaluation script, you can use the following command:

```bash
python3 evaluate.py \
    --var-path "/path/to/your/variants.vcf" \
    --ref-path "/path/to/your/reference.fasta" \
    --outdir "/path/to/your/output" \
    --sample-name "your_sample_name"
```

Usage using example data:

```bash
    python evaluate.py \
        -r example-data/ref/tp53_hg38.fasta \
        -v result/ffpe.calls/known_phi/ffpe.calls_filtered.vcf
```

By default, the output sample id will be inferred from the VCF filename, and the results will be saved to the same directory as the input VCF file.

### Outputs

Running the evaluation script will generate the following files in the specified output directory:

-  **`{your_sample_name}_96c_mutations_profiles.tsv`**: A tab-separated file containing the mutation profile.
-  **`{your_sample_name}_96c_mutations_count.tsv`**: A tab-separated file containing the mutation counts.
-  **`{your_sample_name}_96c_mutations_signature.pdf`**: A PDF image of the 96 channel mutation signature plot.
