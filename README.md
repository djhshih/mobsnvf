# htspan mobsnvf filter

The mobsnvf filter is a modue of htspan for removing artifacts from high thoughput sequencing data. Supported damage types are:

* FFPE
* OXOG

## Setup Procedure

### Dependencies
* R >= 4.0
* R Libraries: stringr, argparse
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

If the dependencies for this pipeline are not installed on your system, a simple way to get them is through conda.

If your linux system doesn't have conda installed, download and install the latest version of miniconda. Skip this step if you already have conda installed on your system.

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

## Minimalisitc Use Case

The most straightforward way to use mobsnvf would be to call the `identify` program and pass in the __*statistics model*__, __*damage type*__, __*bam file*__, __*vcf file*__, __*p-value cutoff*__, and __*output name*__. 

```
hts-mobsnvf identify \
  -M freq \
  -t {damage_type} \
  -b {path-to-bam} \
  -V {path-to-vcf} \
  -g 0 \
  -o "{output-name}.snv" > {output-name}.tsv
```

Make sure your vcf and bam files are paired with their index files within their respective directories.

**Parameters:**

- `-M`: Model to use for identification or quantification. Choices are 'freq' or 'bayes'.

- `-t`:  Type of damage to identify or quantify (e.g. what variant type to analyze). Choices are **'ffpe'** or **'oxog'**.

- `-b`: Path to the BAM data being examined for damage. __Must be paired-end__.

- `-V`: Path to the list of SNVs to be examined for damage identification. Valid formats are plain **TSV** or **VCF/BCF** (may be compressed with gzip or bgzip).

- `-g`: P-value threshold for a variant to be considered damaged. Setting this to `0` will disable automatic filtering but statistics will still be recorded in the output file.

- `-o`: Path for the output list of filtered SNVs.


## Proper FFPE Artifact Filtering Workflow

### Additional Dependencies

Make sure GATK4 is installed and added to `$PATH`.

* GATK4

We can start by making a clone of the mobsnvf GitHub repository and navigating to the directory:

```bash
git clone git@github.com:djhshih/mobsnvf.git
cd mobsnvf
```

### Sample Preparation

This pipeline requires sequence alignment __BAM__ and called variants __VCF__ files. Therefore, place the BAM and VCF of your samples in their respective _`BAM/`_ and _`VCF/`_ directory.

__Important:__ Make sure the __BAM__ and __VCF__ files have the same name for each sample. Example: __`sample01.bam`__ and __`sample01.vcf`__.

In case they don't you may add a prefix or suffix to the bam files using the `--prefix-bam` and `--suffix-bam` respectively. Similarly, you may add the prefix and suffix to the vcf files using the `--prefix-vcf` and `--suffix-vcf` parameters.

It also requires a reference genome to align against. The reference genomes may be obtained by navigating to the _`ref/`_ directory and running the `get.sh` script:

```bash
bash get.sh
```

### Running samples for filtering

After placing the samples in the _`bcf/`_ and _`vcf/`_ directory, the pipeline may finally be run for applying FFPE artifact filtering to our samples.

```bash
bash workflow-mobsnvf --use-phi 'true' --sample_id {your-sample-name} 
```

The `--sample-id` field takes the sample name of the __BAM__ and __VCF__ files. For example, the `--sample-id` input for sample01.bam and sample01.vcf would be `--sample-id sample01`.

The parameter __phi__ is an estimate of the extent of artificial DNA damage. The `--use-phi` flag being set to `true` will calculate phi and use it in when identifying artifacts. When this flag is set to `false` it will skip the phi calculation and proceed with artifact identification without fixing phi.

#### Results

The result for this workflow will be saved to the `results/` directory, followed by the __sample name__ and __phi__. Example `./results/sample01/known_phi/`.

The outputs files are:

- `{sample_name}_ffpe.snv` : This contains the list of SNV files with their P-value.
- `{sample_name}_selected.vcf` : This contains the genetic variants with FFPE Artifacts filtered out.
- `{sample_name}_selected.vcf.idx` : Index of the `{sample_name}_selected.vcf`.
- `{sample_name}_removed.vcf` : This contains the genetic variants which are presumably FFPE artifacts.
- `{sample_name}_removed.vcf.idx` : Index of the `{sample_name}_removed.vcf`.













