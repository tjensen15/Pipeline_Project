# Pipeline_Project

This is a pipeline that includes various functions and options to compare transcriptomes based on a given condition. For the test data included, the pipeline will compare the Human cytomegalovirus (HCMV) transcriptomes 2 and 6 days post infection (dpi).

## Installation
You can clone the repository for all scripts and example data with the command below:
```
git clone https://github.com/tjensen15/Pipeline_Project
```
### Dependencies
To run the wrapper.py script, the following tools and libraries are necessary:

1. R: R with sleuth library downloaded. Here is guidance on how to download it: [Sleuth R Package Download Information] (https://pachterlab.github.io/sleuth/download)
2. Python: Must have the following libraries downloaded: argparse, logging, subprocess, os, Biopython, collections, and re. 
3. Command line tools: Must have the following tools downloaded (documentation hyperlinked): [kallisto] (https://pachterlab.github.io/kallisto/download), [bowtie2] (https://github.com/BenLangmead/bowtie2), [SPAdes] (https://github.com/ablab/spades?tab=readme-ov-file), and [BLAST+] (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)

## Parameters description
Main parameters:
```
Usage: python wrapper.py --acc <ACCESSION #> --email <YOUR EMAIL> --srr <LIST OF SRR, SPACE DELIM: SRRX SRRX SRRX ...> --samples_file samples.txt --step [OPTIONS]

    - fetch
    - index
    - quant
    - analyze
    - bowtie2-build
    - bowtie2
    - spades
    - blast
    - all       will run all steps above

```
Need to have steps/outputs from any steps above the one you want to run individually.

### General usage:
```
### To be run inside the github repo folder
python wrapper.py --acc NC_006273.2 --email johnsmith@gmail.com --srr SRR5660030_sample SRR5660033_sample SRR5660044_sample SRR5660045_sample --samples_file samples.txt --step all
```
Samples text file is needed for parsing ID names, conditions and for an input for sleuth R script. If a samples file is not included, empty dictionaries will be returned and sleuth will not run.

General format:
```
sample  condition       path
SRR5660030_sample       2dpi    results/SRR5660030_sample
SRR5660033_sample       6dpi    results/SRR5660033_sample
SRR5660044_sample       2dpi    results/SRR5660044_sample
SRR5660045_sample       6dpi    results/SRR5660045_sample
```
This specific pipeline uses HCMV genome and nucleotide database of the Betaherpesvirinae subfamily. If this is not the same genomes and databases you are analyzing, they will need to be downloaded using NCBI's command line tools and need to replace the genomic.fna file and betaherpesvirinae folder.