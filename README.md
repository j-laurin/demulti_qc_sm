# QC and demultiplexing of the libraries for Vilardo lab

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow for QC and demultiplexing of Vilardo libraries


## Installation

### Download this Workflow

Download this Snakemake workflow, e.g., with `git clone`:

  ```bash
  git clone git@github.com:j-laurin/demulti_qc_sm.git
  ```


### Fill in the sample list in config folder and transfer raw fastq.gz files into the data folder. Fill in the table as usual if using MS Excel. 

Sample file you can download here: https://github.com/j-laurin/demulti_qc_sm/blob/main/config/samples.csv

or use the file from the repository. If downloading as an extra file, replace the samples.csv with edited one your cloned repository.

The file name and the indexes used are seprated by a single comma in texteditors.


| file_name   | indexes | 
| ----------- | --------|
| EV230201    | 1234    | 
| EV230202    | 5678    | 

### Following dependencies are required for the pipeline:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/#getting-started)
- cutadapt
- fastqc
- multiqc
- seaborn


To install dependecies, create a conda environment for this pipeline: 

```bash
conda create -n lib_qc -c conda-forge -c bioconda snakemake cutadapt fastqc multiqc seaborn
```

Then activate the environment:

```bash
conda activate lib_qc
```

And run the pipeline:

```bash
snakemake -c 2 "get_all_results.txt"
```

 Replace '2' with the number of cores you wish to use, at least 2.

To see how many cores you have available on linux machines, run:
  
  ```bash
  nproc
  ```

### All your results will be then located in the results folder. 

### Settings:

You can modify the minimum length of sequences kept by modifying the 'length_cutoff' 
in the config.yaml file in the config folder.
The default is set to 19 nt (lowest nt count after adaptor removal can be 10nt). 

  ```bash
  length_cutoff: 19
  ```
