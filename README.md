# HTSPAN MOBSNVF Filter

The mobsnvf filter is a module of htspan for removing artifacts from high throughput sequencing data. Supported damage types are:

*   FFPE
*   OXOG

This repository provides a pipeline for using the mobsnvf module of htspan for a production workflow.

## Dependencies

*   htspan (see **[djhshih/htspan](https://github.com/djhshih/htspan)** for instructions on setting up htspan)
*   R 4.4.3
*   python 3.12.10
*   gatk 4.6.2.0

### R libraries:
*   stringr 1.5.1
*   argparser 0.7.2

### Python libraries:
*   polars 1.31.0
*   pysam 0.23.3
*   numpy 2.3.0
*   seaborn 0.13.2
*   matplotlib 3.10.3

### Optional dependencies:
*   gsutil
*   dlazy (see **[djhshih/dlazy](https://github.com/djhshih/dlazy/tree/main)** for instructions on setting up dlazy)

## Setup Procedure

If your system does **git** installed, please setup and configure **git** first. You may take a look here for a **[demonstration](https://github.com/djhshih/setup-linux)**.

```bash
sudo apt install git
```

### Install htspan

Please refer to **[djhshih/htspan](https://github.com/djhshih/htspan)** for instructions on setting up **htspan** and a minimalistic use case.

### Clone mobsnvf repo

Then clone this `mobsnvf` repo:

```bash
git clone https://github.com/djhshih/mobsnvf.git
```

## Installing dependencies

The dependencies for **htspan** and its installation guide is available in its own repository (**[djhshih/htspan](https://github.com/djhshih/htspan)**), as mentioned previously. The remaining dependencies can be installed as mentioned below.

This main artifact filtering pipeline relies on R, two R packages, and GATK.

GATK can be installed according to the [instructions provided here](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4).

The R packages may be installed as follows in a debian-based system (e.g., Ubuntu):

```bash
sudo apt install r-base r-base-dev
R -e 'install.packages(c("argparser","stringr"), repos="https://cloud.r-project.org")'
```

The evaluation of the data requires python and some additional python libraries.

Python can be installed as follows in a debian-based system (e.g., Ubuntu):

```bash
sudo apt update && sudo apt install python3 python3-pip
```

Then the python libraries can be installed using pip:

```bash
pip install polars pysam numpy seaborn matplotlib
```

### Alternate methods

#### Docker

Alternatively, you may use the `Dockerfile` included with this repo to build and run a docker container. This solves all dependencies for compiling htspan and running this pipeline.

If docker is not installed on your system, follow instructions **[here](https://docs.docker.com/engine/install/ubuntu/)** to install docker or contact your system administrator if you don't have the privileges.

Then navigate to the cloned `mobsnvf/` repo and use the following commands to make a Docker image and run an interactive Docker container:

```bash
docker build -t mobsnvf .
docker run --rm -it -v .:/home/work -w /home/work mobsnvf
```

This will create and run an interactive Docker container, mounting the current directory into the `/home/work` directory inside the container.

_**Re-emphasis**_: The docker image includes the dependencies for compiling **htspan** but not **htspan** itself. Therefore, it still needs to be compiled from source as described in **[djhshih/htspan](https://github.com/djhshih/setup-linux)**. **Dlazy** is also not included in the Docker image, so you will need to install it separately if you want to use it.

#### Conda

Except for htspan, the rest of the dependencies can also be easily installed using Conda. This is useful if you don't have administrative privileges to install system-wide packages or if you prefer an isolated environment.

If you don't have Conda installed, you can follow the instructions [here](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-terminal-installer).

Once Conda is installed, create a new environment and install the required packages with the following commands:

```bash
conda create -n mobsnvf_test -c conda-forge -c bioconda -c defaults \
r-base==4.4.2 r-argparser=0.7.2 r-stringr=1.5.1 \
python=3.12.11 pysam=0.23.3 numpy=2.3.0 polars=1.31.0 matplotlib=3.10.3 seaborn=0.13.2 \
gatk4=4.6.2.0
```

Then activate this environment with:

```bash
conda activate mobsnvf
```
And perform the tasks within this environment.

_**Re-emphasis**_: The **htspan** and its dependency is not available via Conda, so you will still need to install its dependencies separately and compile it from source as described in the **Install htspan** section above.

## How to use this pipeline

### Command-Line Arguments

| Argument | Alias(es) | Required? | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--bam-path` | `-b`, `--bam` | **Yes** | | Path to the input BAM file. Must be indexed (`.bai`). |
| `--vcf-path` | `-v`, `--vcf` | **Yes** | | Path to the input VCF file. Must be indexed (`.tbi` or `.csi` if gzipped). |
| `--ref-path` | `-r`, `--ref` | **Yes** | | Path to the reference genome FASTA file. Must be indexed (`.fai`). |
| `--sample-id` | `-i`, `--id` | No | Inferred from VCF | A unique identifier for the sample. If not provided, it will be automatically derived from the VCF filename. |
| `--use-phi` | `-p` | No | `true` | Set to `true` to perform phi estimation or `false` to skip it. |
| `--damage-type`| `-d`, `--damage` | No | `ffpe` | The type of DNA damage to model. Supported values are `ffpe` and `oxog`. |
| `--out-dir` | `-o`, `--out` | No | `./result` | The base directory where results will be saved. |
| `--fp-cut` | `-c`, `--cut` | No | `0.5` | False positive cutoff value for the FDR filter step. |

### Example Usage

Below is a typical command to run the pipeline for a single sample. You must provide the paths to your BAM, VCF, and reference genome.

The workflow requires explicit file paths for inputs unless they are placed in the `ref`, `bam` and `vcf` directories within the repository. In this case they can be just referred to by their filenames. The input files i.e. the BAM, VCF and Reference Fenome must be indexed.

**Note:** It is important that the Reference Fasta provided is the same as the one used for sequence alignment.

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

**Note:** This requires gsutil to be installed on your system. If you don't have gsutil, you can install it by following the instructions [here](https://cloud.google.com/storage/docs/gsutil_install). This already comes prepackaged if you are following the **conda** or **docker** methods mentioned above.

```bash
cd ref
bash get.sh
```

### Phi Parameter

The parameter **phi** is an estimate of the extent of artificial DNA damage.

Setting the `--use-phi` flag to `true` (the default) will estimate phi with `hts-mobsnvf quantify` and use this single value when identifying artifacts across all sites with `hts-mobsnvf identify`. This is the recommended approach as it generally improves accuracy.

When `--use-phi` is set to `false`, the initial phi estimation step is skipped, and `hts-mobsnvf identify` will instead estimate the extent of DNA damage (phi) individually for each SNV site.

### Testing the Pipeline with Example Data

This repository includes an `example-data` directory containing a simple ffpe BAM and VCF file along with a matched normal VCF. This allows you to perform a test run to ensure that your environment is configured correctly and the pipeline works as expected.

You can run the test using the following command from the root of the repository.

```bash
bash mobsnvf-workflow.sh \
    --ref example-data/ref/tp53_hg38.fasta \
    --bam example-data/bam/ffpe.bam \
    --vcf example-data/vcf/ffpe.calls.vcf.gz
```

Upon successful completion, you will find the results in the `result/ffpe.calls/known_phi/` directory.


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

### Running multiple samples

To process multiple samples, we recommend creating a manifest file and looping through it. However, you can feel free to adapt to your own needs.

1.  **Create a manifest file**

    Create a tab-separated values (TSV) file (e.g., `samples.tsv`) with absolute paths to the corresponding BAM and VCF files.

    **Example `samples.tsv`:**

    ```tsv
    /path/to/bams/sample01.bam	/path/to/vcfs/sample01.vcf.gz
    /path/to/bams/sample02.bam	/path/to/vcfs/sample02.vcf.gz
    /path/to/bams/sample03.bam	/path/to/vcfs/sample03.vcf.gz
    ```

2.  **Run the workflow in batch**

    You can now use a simple `bash` loop to read the manifest file and create a script for each sample. Then you may run these scripts in batch mode.

    ```bash
    #!/bin/bash
    
    ref="/path/to/your/reference.fasta"
    out="/path/to/your/output/results"
    
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

    Dlazy prepares logs for each job and makes it easy to resume or rerun jobs in case of failures.

## Evaluation

An evaluation python script (`evaluate.py`) is included with this repository to analyse the mutation profile, mutation counts, and mutation signature plots for the vartiant (VCF/SNV) files.

### Example Usage

To run the evaluation script, you can use the following command:

```bash
python evaluate.py \
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

This script can also be used to evaluate tabular variant data such as .snv files. This will work as long as the file is tab delimited and includes the **[`chrom`, `pos`, `ref`, `alt`]** columns.

By default, the output sample id will be inferred from the VCF filename, and the results will be saved to the same directory as the input VCF file.

### Outputs

Running the evaluation script will generate the following files in the specified output directory:

-  **`{your_sample_name}_96c_mutations_profiles.tsv`**: A tab-separated file containing the mutation profile.
-  **`{your_sample_name}_96c_mutations_count.tsv`**: A tab-separated file containing the mutation counts.
-  **`{your_sample_name}_96c_mutations_signature.pdf`**: A PDF showing the 96 channel mutation signature plot along with total variant count and total C>T mutation count.


## Issues to be resolved

Currently, htspan mobsnvf does not handle multiallelic sites very well. As of right now multiallelic sites are only handled as intended if the **ref** and first **alt** allele makes up for a **`C>T`** mutation. In this case the site is split into multiple biallelic sites. In all other multiallelic cases only the first alt allele is retained in the output and the successive alleles besides are discarded. It also discards all INDELs from its output.

### How it affects the workflow

The output from the mobsnvf module, i.e. `{sample_id}_{damage_type}.snv`, is only used to identify and mask artifacts in the VCF, which is used to generate the `{sample_id}_artifacts.vcf` and `{sample_id}_filtered.vcf` files. Therefore, the issue is mostly mitigated and only affects cases where the ref and the successive alt alleles in a multiallelic site result in a `C>T` mutation. In such cases, these mutations are not analyzed for artifact identification. For example:

| CHR | POS | REF | ALT      |
|-----|-----|-----|----------|
| ... | ... | C   | A, T     |
| ... | ... | C   | G, T     |
| ... | ... | C   | GAA, T   |
| ... | ... | G   | C, A     |
| ... | ... | G   | TC, A    |

In cases like this the C>T mutations are not considered for artifact identification.

### Workaround
If your VCF contains multiallelic sites, you may need to preprocess it to split these sites into biallelic ones before running the mobsnvf workflow.

An easy way to do this is by using **bcftools**:

```bash
bcftools norm -m - -o output.vcf input.vcf
```
