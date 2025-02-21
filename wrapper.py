import argparse
import logging
import subprocess
import os
import time
from Bio import Entrez, SeqIO
import pandas as pd
from collections import defaultdict
import re

def fetch_genome(acc, email, output_fasta):
    '''Fetch genome from NCBI and extract CDS sequences'''
    Entrez.email = email # email from argparse entry
    
    try:
        # fetching from nucleotide database with id entered in argparse
        handle = Entrez.efetch(db='nucleotide', id=acc, rettype='gb', retmode='text')
        record = SeqIO.read(handle, "genbank") # getting sequence
        handle.close()
    except Exception as e:
        logging.error(f"Error fetching genome: {e}") # prints to log if there was an error getting genome
        return None

    # sums the CDS by counting how many features contain CDS
    cds_count = sum(1 for feature in record.features if feature.type == "CDS") 
    logging.info(f"The genome {acc} has {cds_count} CDS.") # outputs to log file

    # writes CDS to an output fasta file
    with open(output_fasta, "w") as out_file:
        for feature in record.features:
            if feature.type == "CDS": # if the feature type is CDS
                cds_sequence = feature.extract(record.seq) # get the sequence
                protein_id = feature.qualifiers.get("protein_id", ["unknown"])[0] # get protein ID
                out_file.write(f">{protein_id}\n{cds_sequence}\n") # write the ID and CDS seq to fasta
    
    return output_fasta

def build_kallisto_index(fasta_file, index_file):
    '''Build kallisto index from CDS sequences'''
    command = ["kallisto", "index", "-i", index_file, fasta_file] # command for kallisto using CDS seqs we just got
    try:
        subprocess.run(command, check=True) # run the command
    except subprocess.CalledProcessError as e:
        logging.error(f"Error building Kallisto index: {e}") # log the error if we can't build index

def run_kallisto_quant(index_file, srr_ids):
    '''Run Kallisto quantification for each SRR ID'''
    subprocess.run(["mkdir", "-p", "results"]) # make a results directory for following commands
    
    for id in srr_ids:
        output_dir = f"results/{id}" # output directory as results/id 
        # command to run kallisto quant using index, 10 bootstraps, 2 processors, and the paired end fastq files
        command = f"kallisto quant -i {index_file} -o {output_dir} -b 10 -t 2 SRA/{id}_1.fastq SRA/{id}_2.fastq"
        try:
            subprocess.run(command, shell=True, check=True) # run via shell
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running Kallisto quant for {id}: {e}") # if there is an error, print to log file

def process_abundance(samples, log_file):
    '''Process abundance.tsv files and log TPM statistics'''
    with open(log_file, "a") as log:
        log.write("Sample\tCondition\tMin TPM\tMedian TPM\tMean TPM\tMax TPM\n") # create header in log file
        for sample_id, condition in samples.items(): # for each sample id and condition
            file_path = f"results/{sample_id}/abundance.tsv" # get the file path to abundance file
            try:
                # adding delay to ensure the file is fully written
                time.sleep(2)
                if not os.path.exists(file_path): 
                    # if the path does not exist, output to log file that it does not exist
                    logging.error(f"File not found: {file_path}")
                    continue
                df = pd.read_csv(file_path, sep="\t") # read in tsv file as dataframe
                df = df.drop(['length', 'eff_length', 'est_counts'], axis=1) # drop columns we are not looking at
                # get the sample id, condition, min tpm, median tpm, mean tpm, and max tpm
                result_line = f"{sample_id}\t{condition}\t{df['tpm'].min()}\t{df['tpm'].median()}\t{df['tpm'].mean()}\t{df['tpm'].max()}\n"
                log.write(result_line) # write the result to the log file
            except Exception as e:
                logging.error(f"Error processing {sample_id}: {e}") # if there is an error with any sample, output to log file

def run_r_analysis(log_file):
    '''Run R script for differential expression analysis.
    Need to make own samples.txt table for sleuth to read, details in documentation'''
    try:
        subprocess.run(["Rscript", "sleuth_script.R"], check=True) # run r script
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running R script: {e}") # if there is an error, write to log file
        return
     # check if the output file exists before trying to read it
    if not os.path.exists("sleuth_results.txt"):
        logging.error("sleuth_results.txt not found. Skipping analysis.")
        return

    try:
        df = pd.read_csv('sleuth_results.txt', sep = ' ') # space separated file
        df = df[['target_id','test_stat','pval','qval']] # extracting only columns of interest

        df_filtered = df[df['qval'] < 0.05] # filter significant results, want anything below 0.05
        print(df_filtered)
        # append results as tab separated with header
        df_filtered.to_csv(log_file, sep="\t", index=False, mode="a", header=True)
    except Exception as e:
        logging.error(f"Error writing to file: {e}") # if there is an error, write to log file

def run_bowtie2_build():
    '''Build Bowtie2 index'''
    command = ["bowtie2-build", "genomic.fna", "bowtie_index"] # command to build index with sample database
    try:
        subprocess.run(command, check=True) # run command
    except subprocess.CalledProcessError as e:
        logging.error(f"Error building Bowtie2 index: {e}") # if there is an error, write it to log file

def run_bowtie2_mapping(ids, samples, patient_ids, log_file):
    '''Run Bowtie2 mapping for given SRR IDs.'''
    subprocess.run(["mkdir", "-p", "bowtie2_results"]) # make a bowtie2_results directory
    for id in ids:
        # for each id, map paired end fastq files to the bowtie index, output aligned sam file, but also an aligned fastq file for counting # of reads later
        command = f"bowtie2 -x bowtie_index -1 SRA/{id}_1.fastq -2 SRA/{id}_2.fastq -S bowtie2_results/{id}_aligned.sam --al-conc bowtie2_results/{id}_aligned.fastq -p 2 > {id}_bowtie2.log 2>&1 "
        try:
            subprocess.run(command, shell=True, check=True) # run via shell
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running Bowtie2 for {id}: {e}") # if there is an error, write to log file
    def count_reads(filename):
        '''Counts number of reads in a FASTQ file'''
        with open(filename, 'r') as file:
            # sum the number of lines that contain @SRR, since this indicates a read
            return sum(1 for line in file if re.match(r'^@SRR', line)) # re.match uses regular expression

    with open(log_file, 'a') as log:
        for sample, condition in samples.items(): # for each sample and condition
            total_reads_before = count_reads(f"SRA/{sample}_1.fastq") # get total reads before mapping
            total_reads_after = count_reads(f"bowtie2_results/{sample}_aligned.1.fastq") # get total reads after mapping from file with only aligned reads

            donor = "Unknown" # initialize donor
            # if the condition is 2dpi, donor is donor 1
            if condition == "2dpi":
                donor = "Donor 1"
            elif condition == "6dpi": # if condition is 6dpi, donor is donor 3
                donor = "Donor 3"
            else: # if none of these conditions are true, output to log file that this needs to be changed
                log.write('Please alter run_bowtie2_mapping code to match your conditions.')

            # append to log
            log.write(f"{donor} ({condition}) had {total_reads_before} read pairs before Bowtie2 filtering and {total_reads_after} read pairs after.\n")

def concatenate_and_run_spades(patient_ids, log_file):
    '''Concatenate FASTQ files and run SPAdes assembly dynamically for multiple patients'''
    bowtie2_dir = os.path.abspath("bowtie2_results") # get absolute path of bowtie2_results
    results_dir = os.path.abspath("results") # get absolute path of results 

    with open(log_file, "a") as log:
        # loop through each patient group
        for donor, ids in patient_ids.items():
            # get file paths for one set of paired end reads 1 as a list
            aligned_1_files = [os.path.join(bowtie2_dir, f"{id}_aligned.1.fastq") for id in ids]
            # get file paths for other set of paired end reads 2 as a list
            aligned_2_files = [os.path.join(bowtie2_dir, f"{id}_aligned.2.fastq") for id in ids]
            # output file names for concatenated files, one per donor
            concatenated_1 = f"{donor}_aligned.1.fastq"
            concatenated_2 = f"{donor}_aligned.2.fastq"

            # concatenate all the aligned 1 files into a single file
            with open(concatenated_1, 'w') as outfile:
                for fname in aligned_1_files: # iterate over each input file
                    with open(fname) as infile: # open each file in read mode
                        outfile.write(infile.read()) # read the file content and write it to the output file

            # concatenate all the aligned 2 files into a single file, same process as above
            with open(concatenated_2, 'w') as outfile:
                for fname in aligned_2_files:
                    with open(fname) as infile:
                        outfile.write(infile.read())

            # adding a small delay to ensure files are fully written
            time.sleep(2)

            # run spades using kmer size of 77, adjust if needed, 2 processors and only assembler mode using the concatenated files
            spades_command = f"spades.py -k 77 -t 2 --only-assembler -1 {concatenated_1} -2 {concatenated_2} -o {donor}_assembly/"
            log.write(f"{spades_command}\n") # write the spades command used to the log file
            try:
                subprocess.run(spades_command, shell=True, check=True) # run command via shell
            except subprocess.CalledProcessError as e:
                logging.error(f"Error running SPAdes for {donor}: {e}") # log the error
                log.write(f"Error running SPAdes for {donor}: {e}\n") # write the specific error to log file
                
def parse_samples(file_path):
    '''Parses a text file to generate patient ID lists and a samples dictionary'''
    
    # read the text file into a DataFrame
    df = pd.read_csv(file_path, sep="\t", header=0, names=["sample", "condition", "path"])
    
    # create the samples dictionary 
    samples = {row["path"].split("/")[-1]: row["condition"] for _, row in df.iterrows()}

    # create patient ID groups based on condition using a defaultdict list to make sure no repeats
    condition_groups = defaultdict(list)
    
    for _, row in df.iterrows():
        sample_id = row["path"].split("/")[-1]  # extract sample ID from the path
        condition_groups[row["condition"]].append(sample_id) # append sample id to condition groups list based on condition
    
    # create patient ID lists for each sample ID and their condition
    patient_ids = {f"patient{i+1}_ids": ids for i, ids in enumerate(condition_groups.values())}

    return patient_ids, samples


def main():
    '''Main function to parse arguments and execute selected pipeline steps'''
    parser = argparse.ArgumentParser(description="RNA-seq pipeline wrapper.")
    parser.add_argument("--acc", help="NCBI accession number for genome.")
    parser.add_argument("--email", help="Your email (required by NCBI).")
    parser.add_argument("--srr", nargs="+", help="List of SRR IDs.")
    parser.add_argument("--samples_file", help="Path to the samples.txt file.")
    parser.add_argument("--step", choices=["fetch", "index", "quant", "analyze", "bowtie2-build", "bowtie2", "spades", "all"], required=True, help="Specify the pipeline step to run.")
    
    args = parser.parse_args()
    logging.basicConfig(filename="Pipeline_Project.log", level=logging.INFO, format="%(message)s")
    fasta_file = "CDS_sequences.fasta"
    index_file = "index.idx"
    # parse samples.txt file
    if args.samples_file:
        patient_ids, samples = parse_samples(args.samples_file)
    else:
        patient_ids, samples = {}, {} # creates empty dictionaries if sample file not included

    if args.step == "fetch":
        if not args.acc or not args.email:
            logging.error("Fetching genome requires --acc and --email.")
        else:
            fetch_genome(args.acc, args.email, fasta_file)
    elif args.step == "index":
        build_kallisto_index(fasta_file, index_file)
    elif args.step == "quant":
        if not args.srr:
            logging.error("Quantification requires --srr IDs.")
        else:
            run_kallisto_quant(index_file, args.srr)
    elif args.step == "process_abundance":
        process_abundance(samples, "Pipeline_Project.log")
    elif args.step == "analyze":
        run_r_analysis()
    elif args.step == "bowtie2-build":
        run_bowtie2_build()
    elif args.step == "bowtie2":
        if not args.srr:
            logging.error("Bowtie2 mapping requires --srr IDs.")
        else:
            run_bowtie2_mapping(args.srr, samples, patient_ids, 'Pipeline_Project.log')
    elif args.step == "spades":
        # Now dynamically use the patient_id dictionary
        concatenate_and_run_spades(patient_ids, "Pipeline_Project.log")
    elif args.step == "all": # if all commands are wanting to be run:
        if not args.acc or not args.email or not args.srr:
            logging.error("Running all steps requires --acc, --email, and --srr IDs.")
        else:
            fetch_genome(args.acc, args.email, fasta_file)
            build_kallisto_index(fasta_file, index_file)
            run_kallisto_quant(index_file, args.srr)
            process_abundance(samples, "Pipeline_Project.log")
            run_r_analysis("Pipeline_Project.log")
            run_bowtie2_build()
            run_bowtie2_mapping(args.srr, samples, patient_ids, 'Pipeline_Project.log')
            concatenate_and_run_spades(patient_ids, "Pipeline_Project.log")

if __name__ == "__main__":
    main()


# python wrapper2.py --acc NC_006273.2 --email hjensen2@luc.edu --srr SRR5660030_sample SRR5660033_sample SRR5660044_sample SRR5660045_sample --samples_file samples.txt --step all