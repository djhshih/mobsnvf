# htspan mobsnvf filter

The htspan mobsnvf filter is a toolkit for removing artifacts from high thoughput sequencing data. Supported damage types are:

* FFPE
* OXOG

## Setup Procedure

### Dependencies
* R >= 4.0
* gcc >= 4.8
* bzip2 >= 1.0 (for htslib)
* liblzma SDK >= 5.2.4 (for htslib)
* libcurl SDK >= 7.64 (for htslib, any implementation)

### Optional dependencies for testing and development

* samtools >= 1.8
* python3 >= 3.6
* boost >= 1.69
* bcftools >= 1.3.1

## Install Conda

If your linux system doesn't have conda installed, download and install the latest version of miniconda. Skip this step is you already have conda installed on your system.

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the instructions and until conda is successfully installed in your system. This typically goes as follows:

- When prompted Enter `yes` to agree to the TOS.

- Press **Return** to accept the default install location (`PREFIX=/Users/<USER>/miniconda3`), or enter another file path to specify an alternate installation directory. The installation might take a few minutes to complete.

- When choosing initialization option, type `Yes`. This will conda modify your shell configuration to initialize conda whenever you open a new shell and to recognize conda commands automatically.


## Create environment with necessary dependencies

Now create a conda environment for running the hts-mobsnvf filter

```bash
conda create --name hts-mobsnvf
```
Aftewards, we activate this newly created environment and install the necessary dependencies for running hts-mobsnvf:

```bash
conda activate hts-mobsnvf

conda install -c conda-forge -c bioconda r-base r-essentials r-stringr r-argparse bzip2 gatk4 git htslib
```


## Install hts-mobsnvf filter

Now we can install the hts-mobsnvf filter in our system. Clone the github repo (make sure you have access since this is a private repo)

```bash
git clone https://github.com/djhshih/htspan
cd htspan
git submodule update --init --recursive
make
```

If the compilation is successful, the binaries will appear in `bin`.

The binaries may be added to `$PATH` if you want them to be accessible from any directory:

- Open `~/.bashrc` file

```bash
nano ~/.bashrc
```

Now navigate to the end and add the following line to the file. 
```txt
export PATH="{path-to-bin-folder}:$PATH"
```

Replace `{path-to-bin-folder}` with the actual path to where htspan binaries are compiled. You may obtain this by navigating to the _`bin/`_ directory and entering `pwd` in the terminal.

Press `Ctrl + O` and then `Enter` to save the file. And then press `Ctrl + X` to exit.

Now re-initialize the shell by entering `exec bash`.

Now __hts-mobsnvf__ should be accessible from any directory. You may test this by typing `hts-mobsnvf` in the terminal.

## Running the FFPE Artifact Filtering Workflow

### Sample Preparation

This pipeline requires sequence alignment __BAM__ files and called variants __VCF__. Therefore, place the BAM and VCF for your samples in their respective _`BAM/`_ and _`VCF/`_ directory. 

__Important:__ Make sure the __BAM__ and __VCF__ files have the same name for each sample. Example: sample01.bam and sample01.vcf.

If they don't you can add a prefix or suffix to the bam files using the  `--add-prefix-bam` and `--add-suffix-bam` respectively. Similarly, you may add the prefix and suffix to the vcf files using the `--add-prefix-vcf` and `-add-suffix-vcf` parameters.

It also requires a reference genome to align against. The reference genomes may be obtained by navigating to the _`ref/`_ directory and running the `get.sh` script:

```bash
bash get.sh
```

### Running samples for filtering

After placing the samples in the _`bcf/`_ and _`vcf/`_ directory, the pipeline may finally be run for applying FFPE artifact filtering to our samples.

```bash
bash workflow-mobsnvf --sample_id {your-sample-name} --use-phi 'true'
```

The `--sample-id` field takes the sample name of the __BAM__ and __VCF__ files. For example, the `--sample-id` input for sample01.bam and sample01.vcf would be `--sample-id sample01`.

The `--use-phi` flag being set to `true` will calculate phi and identify artifacts with a fixed phi. When this flag is set to `false` it will skip the phi calculation and identify aritfacts without fixing phi.









