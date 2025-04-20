# HTSPAN MOBSNVF Filter

The mobsnvf filter is a modue of htspan for removing artifacts from high thoughput sequencing data. Supported damage types are:

* FFPE
* OXOG

This repository provides a pipeline for using the mobsnvf module of htspan for a production workflow.

## Setup Procedure

Please refer to [djhshih/htspan](https://github.com/djhshih/htspan) for setup procedure and example of a minimalistic use case.

## Dependencies
* R >= 4.0
* R Libraries: srtringr, argparser
* GATK >= 4.2.2.0
* gsutil >= 5.00
* gcc >= 4.8
* bzip2 >= 1.0 (for htslib)
* liblzma-dev >= 5.2.4 (for htslib)
* libcurl-dev >= 7.64 (for htslib, any implementation)

## Installing dependencies

Follow the guideline in [GitHub: djhshih/htspan](https://github.com/djhshih/htspan) to install the dependencies for htspan.

Follow the instructions below to install the R, R package and GATK dependencies for this pipeline:

```bash
sudo apt install r-base r-base-dev
R -e 'install.packages(c("argparser","stringr"), repos="https://cloud.r-project.org")'
conda install gatk4
```

### Alternate method:

Alternatively, you may use the included Dockerfile to build and run a docker container which includes all the dependencies. This is covered at [GitHub: djhshih/htspan](https://github.com/djhshih/htspan).

## How to use this pipeline

### Sample Preparation

This pipeline requires sequence alignment (__BAM__) files and called variants (__VCF__) files. Hence, place the BAM and VCF of your samples in their respective _`BAM/`_ and _`VCF/`_ directory.

__Important:__ Make sure the __BAM__ and __VCF__ files have the same name for each sample. Example: __`sample01.bam`__ and __`sample01.vcf`__.

In case your bam and vcf have prefixes or suffixes, you may use the `--prefix-bam` and `--suffix-bam` for declaring prefix and suffix to your bam file. Similarly, you may add the prefix and suffix to the vcf files using the `--prefix-vcf` and `--suffix-vcf` parameters.

__Example:__
If your bam has the name `sample_01_aligned.bam` and your vcf have the name `sample_01_called.vcf`. You may declare this using: `--sample-id sample_01 --suffix-bam _aligned --suffix-vcf _called`

### Get reference genome

Run the `get_reference_genome.sh` script to download the reference genomes in the _`ref/`_ dirctory

```bash
bash get_reference_genome.sh
```

### Running samples for filtering

After placing the samples in the _`bcf/`_ and _`vcf/`_ directory, the pipeline may be run for FFPE artifact filtering.

```bash
bash workflow-mobsnvf --use-phi 'true' --sample_id {your-sample-name} 
```

The `--sample-id` field takes the sample name of the __BAM__ and __VCF__ files. For example, the `--sample-id` input for __`sample01.bam`__ and __`sample01.vcf`__ would be `--sample-id sample01`.

The parameter __phi__ is an estimate of the extent of artificial DNA damage. The `--use-phi` flag being set to `true` will calculate phi and use it in when identifying artifacts. When this flag is set to `false` it will skip the phi calculation and proceed with artifact identification without fixing phi.

#### Results

The result for this workflow will be saved to the `results/` directory, followed by the __sample name__ and __phi__. Example `./results/sample01/known_phi/`.

The outputs files are:

- `{sample_name}_ffpe.snv` : This contains the list of SNV files with their P-value.
- `{sample_name}_selected.vcf` : This contains the genetic variants with FFPE Artifacts filtered out.
- `{sample_name}_selected.vcf.idx` : Index of the `{sample_name}_selected.vcf`.
- `{sample_name}_removed.vcf` : This contains the genetic variants which are presumably FFPE artifacts.
- `{sample_name}_removed.vcf.idx` : Index of the `{sample_name}_removed.vcf`.

## Running multiple samples

An easy way to run multiple samples using this workflow is to use Dlazy. Check __[GitHub: djhshih/dlazy](https://github.com/djhshih/dlazy)__ for instructions on installation and basic use.

You may use dlazy to run your samples in as follows:

- Run the `list_samples.py` script to obtain a list of your samples:

```bash
python list_samples.py
```
This will generate a samples.txt file with listing the sample name of all your samples by referring to the `bam/` direcrtory.

Then you may run your samples using:
```
djobs workflow-mobsnvf --use-phi 'true' --sample_id
dlazy job
```















