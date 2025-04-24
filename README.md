# HTSPAN MOBSNVF Filter

The mobsnvf filter is a module of htspan for removing artifacts from high throughput sequencing data. Supported damage types are:

* FFPE
* OXOG

This repository provides a pipeline for using the mobsnvf module of htspan for a production workflow.


## Dependencies
* R >= 4.0
* R Libraries: srtringr, argparser
* GATK >= 4.2.2.0
* gsutil >= 5.00
* bzip2 >= 1.0 (for htslib)
* cmake
* gcc >= 4.8
* g++
* libbz2-dev
* libcurl-dev >= 7.64 (for htslib, any implementation)
* liblzma-dev >= 5.2.4 (for htslib)
* zlib1g-dev


## Setup Procedure

If your system does not have git, install git:

```bash
sudo apt install git
```

### Compile htspan
Please refer to __[djhshih/htspan](https://github.com/djhshih/htspan)__ for setting up __htspan__ and for an example of a minimalistic workflow.

### Clone mobsnvf repo
Then clone this `mobsnvf` repo:
```bash
git clone https://github.com/djhshih/mobsnvf.git
```


## Installing dependencies

The dependencies for __htspan__ may be installed in Debian/Ubuntu based distributions using this command:
```bash
sudo apt install gcc bzip2 liblzma-dev libcurl4-openssl-dev
```

This pipeline also relies on R, R packages and GATK. These may be installed as follows:

```bash
sudo apt install r-base r-base-dev
R -e 'install.packages(c("argparser","stringr"), repos="https://cloud.r-project.org")'
conda install -c bioconda gatk4
```

If your system does not have __Conda__, you may follow these __[installation instructions](https://www.anaconda.com/docs/getting-started/miniconda/install#linux)__ or install GATK4 through other means.

### Alternate method

Alternatively, you may use the `Dockerfile` included with this repo to build and run a docker container. This solves all dependencies for compiling htspan and running the pipeline.

If docker is not installed on your system, follow instructions __[here](https://docs.docker.com/engine/install/)__ to install docker or contact your system administrator if you don't have the privileges.

Then navigate to the cloned `mobsnvf/` repo and use the following commands to make a Docker image and run an interactive Docker container:

```bash
docker build -t mobsnvf .
docker run --rm -it mobsnvf
```

Alternatively, you may use the `-v` flag to mount your current working directory to your docker container:

```bash
docker run --rm -it -v .:/home/ubuntu mobsnvf
```

This will create and run an interactive Docker container, mounting the current working directory and including all the dependencies.

## How to use this pipeline

### Sample Preparation

This pipeline requires sequence alignment (__BAM__) files and called variants (__VCF__) files. Hence, place the BAM and VCF of your samples in their respective _`bam/`_ and _`vcf/`_ directory along with their associated index files.

__Easiest method:__ Make sure the __BAM__ and __VCF__ files have the same file name for each sample. Example: __`sample01.bam`__ and __`sample01.vcf`__.

In case your bam and vcf have prefixes or suffixes, you may use the `--prefix-bam` and `--suffix-bam` for declaring prefix and suffix for your bam file. Similarly, you may declare the prefix and suffix for the vcf files using the `--prefix-vcf` and `--suffix-vcf` parameters.

__Example:__
If your bam has the name `sample_01_aligned.bam` and your vcf has the name `sample_01_called.vcf`. You may declare this using: `--sample-id sample_01 --suffix-bam _aligned --suffix-vcf _called`

### Reference genome

#### Using your own:
If you are going to use your own reference genome, place the __fasta__ files in the `ref/` directory and make sure they are __indexed__.

Custom reference genome needs to be declared to the pipeline using the `--ref {file-name}`. Replace _{file-name}_ with actual the file name.

#### Downloading a reference genome:
By default the pipeline will look for `Homo_sapiens_assembly38.fasta` which may be obtained by running the `get_reference_genome.sh` script to download this reference genome into the _`ref/`_ directory. 

Downloading the reference genome requires __gsutil__. If this is not installed, you may obtain this by running:

```bash
conda install -c conda-forge gsutil
```
Then download the reference genome using:

```bash
bash get_reference_genome.sh
```

### Running samples for filtering

After placing the samples in the _`bcf/`_ and _`vcf/`_ directory, the pipeline may be run for FFPE artifact filtering.

```bash
bash mobsnvf-workflow.sh --use-phi 'true' --sample-id {your-sample-name} 
```

The `--sample-id` field takes the sample name of the __BAM__ and __VCF__ files. For example, the `--sample-id` input for __`sample01.bam`__ and __`sample01.vcf`__ would be `--sample-id sample01`.

#### Phi
The parameter __phi__ is an estimate of the extent of artificial DNA damage. 

Setting the `--use-phi` flag to `true` will estimate phi with `hts-mobsnvf quantify` and use it in when identifying artifacts with `hts-mobsnvf identify`. 

When `--use-phi` is set to `false`, phi estimation will be skipped and `hts-mobsnvf identify` will estimate the extent of DNA damage (phi) individually for each SNV site.

Typically, accuracy is bit better if phi is estimated using `quantify` first, and then provided to `identify`. Hence, it is generally recommended to set `--use-phi` to **true**.

#### Results

The result for this workflow will be saved to the `results/` directory, followed by the __sample name__ and __phi__. Example `./results/sample01/known_phi/`.

The outputs files are:

- `{sample_name}_ffpe.snv` : This contains the list of SNVs with their P-value.
- `{sample_name}_selected.vcf` : This contains the genetic variants with FFPE Artifacts filtered out.
- `{sample_name}_selected.vcf.idx` : Index of the `{sample_name}_selected.vcf`.
- `{sample_name}_removed.vcf` : This contains the genetic variants which are presumably FFPE artifacts.
- `{sample_name}_removed.vcf.idx` : Index of the `{sample_name}_removed.vcf`.

## Running multiple samples

An easy way to run multiple samples using this workflow is to use __Dlazy__. Check __[djhshih/dlazy](https://github.com/djhshih/dlazy)__ for instructions on installation and basic use.

You may use dlazy to run your samples in as follows:

- Run the `list_samples.py` script to obtain a list of your samples:

```bash
python list_samples.py
```
This will generate a `samples.txt` file with listing the sample name of all your samples by referring to the _`bam/`_ directory. 

Then the pipeline can be run using:
```
djobs workflow-mobsnvf --use-phi 'true' --sample-id
dlazy job
```

### Dealing with prefixes and suffixes

If your bam file name includes suffixes or prefixes along with the sample name, they may be declared using the `--suffix` and `--prefix` arguments to omit them when generating the sample list.

#### Example
Imagine you have a bam file named `prefix_sample_01_suffix.bam` and a vcf file named `prefix_sample_01_suffix.vcf`. 

As the sample list is generated based on the bam files present in the `bam/` directory, the prefix and suffix can be omitted during sample list generation using __`--suffix _suffix --prefix prefix_`__ to just keep the sample name i.e. __sample_01__:

 ```bash
python list_samples.py --suffix _suffix --prefix prefix_
 ```

Then the pipeline with multiple samples can be run with **Dlazy** using:

```
djobs ./mobsnvf-workflow.sh --use-phi 'true' --prefix-bam prefix_ --prefix-vcf prefix_ --suffix-bam _suffix --suffix-vcf _suffix --sample-id

dlazy job
```



















