import os
import datetime
import pandas as pd
from snakemake.utils import validate,available_cpu_count
import time

report: "../report/workflow.rst"

configfile: 'config/config.yaml'
# validate(config, schema="schemas/config_schema.yaml")

# Setup vars
config_threads = int(config["threads"])
run_name = os.path.basename(os.getcwd())

# Load in metadata
metadata_file = "config/metadata.tsv"
metadata_df = pd.read_csv(metadata_file, sep = "\t").set_index('sample_name', drop=False)
# validate(metadata_df, schema="schemas/metadata_schema.yaml")

# Setup samplesheet
samples = metadata_df.index.tolist()

# Convert df to dict for trt/ctrl pairing
metadata_dict = metadata_df.to_dict('index')
chip_samples = [sample for sample in metadata_dict if metadata_dict[sample]['type'] == "chip"]
se_samples = [sample for sample in metadata_dict if metadata_dict[sample]['end_type'] == "se"]
pe_samples = [sample for sample in metadata_dict if metadata_dict[sample]['end_type'] == "pe"]

# Get the raw fastq files to start with
def get_fastq_se(wildcards):
	return metadata_df.loc[(wildcards.sample), "r1"]

def get_fastq_pe(wildcards):
	return metadata_df.loc[(wildcards.sample), ["r1", "r2"]]

def get_sortbam_input(wildcards):
	if metadata_df.loc[(wildcards.sample), "end_type"] == "se":
		return "results/bam/{}.bowtie.bam".format(wildcards.sample)
	elif metadata_df.loc[(wildcards.sample), "end_type"] == "pe":
		return "results/bam/{}.bowtie2.bam".format(wildcards.sample)

def get_bamfilter_params(wildcards):
	if metadata_df.loc[(wildcards.sample), "end_type"] == "se":
		return "-F {remove}".format(remove = config['filter_bam']['se_remove_reads_flags'])
	elif metadata_df.loc[(wildcards.sample), "end_type"] == "pe":
		return "-F {remove} -f {keep}".format(remove = config['filter_bam']['pe_remove_reads_flags'], keep = config['filter_bam']['pe_keep_reads_flags'])
def get_macs2_params(wildcards):
	if metadata_df.loc[(wildcards.chip_sample), "end_type"] == "se":
		return "--format BAM --gsize {macs2_genome} --keep-dup all --nomodel --extsize 200".format(macs2_genome = config['macs2']['gsize'])
	elif metadata_df.loc[(wildcards.chip_sample), "end_type"] == "pe":
		if metadata_df.loc[(metadata_df.loc[(wildcards.chip_sample), "input"]), "end_type"] == "pe":
			return "--format BAMPE --gsize {macs2_genome} --keep-dup all --nomodel".format(macs2_genome = config['macs2']['gsize'])
		elif metadata_df.loc[(metadata_df.loc[(wildcards.chip_sample), "input"]), "end_type"] == "se":
			return "--format BAM --gsize {macs2_genome} --keep-dup all --nomodel --extsize 200".format(macs2_genome = config['macs2']['gsize'])

def get_input_bam(wildcards):
	return "results/bam/{}.sorted.remdup.nonblklst.filt.resort.bam".format(metadata_df.loc[(wildcards.chip_sample), "input"])

def get_read_extension_length(wildcards):
	if metadata_df.loc[(wildcards.sample), "end_type"] == "se":
		return config['bamcoverage_bw']['read_extension_length']
	elif metadata_df.loc[(wildcards.sample), "end_type"] == "pe":
		return ''


