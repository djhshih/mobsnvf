# htspan mobsnvf filter

The htspan mobsnvf filter is a toolkit for removing artifacts from high thoughput sequencing data. Supported damage types are

* FFPE
* OXOG

## Dependencies

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
Here, we created a new conda environment named hts-mobsnvf. Now we can install the necessary dependencies for running hts-mobsnvf:

```bash
conda install -c conda-forge -c bioconda r-base r-essentials r-stringr r-argparse bzip2 gatk4 git
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

Add the binaries to `$PATH` if you want them to be accessible from any directory:

- Open `~/.bashrc` file

```
nano ~/.bashrc
```

Navigate to the end and add the following line to the file:

```
export PATH="{path-to-bin-folder}:$PATH"
```

Replace `{path-to-bin-folder}` with the actual path to where htspan binaries are compiled. You may obtain this by navigating to the `bin` folder and entering `pwd` in the terminal.

Press `Ctrl + O` and then `Enter` to save the file. And then press `Ctrl + X` to exit.

Now re-initialize the shell by entering `exec bash`.

Now `hts-mobsnvf` should be accessible from any directory.


