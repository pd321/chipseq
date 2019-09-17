import os
import datetime
import pandas as pd
from snakemake.utils import validate,available_cpu_count
import time

report: "report/workflow.rst"

configfile: 'config.yaml'
# validate(config, schema="schemas/config_schema.yaml")

# Setup vars
threads_high = available_cpu_count()
threads_mid = int(threads_high/2)
threads_low = int(threads_high/4)
run_name = os.path.basename(os.getcwd())
bdg_types = ['treat_pileup', 'control_lambda', 'input_minus']

# Load in metadata
metadata_file = "metadata.tsv"
metadata_df = pd.read_csv(metadata_file, sep = "\t").set_index('sample_name', drop=False)
# validate(metadata_df, schema="schemas/metadata_schema.yaml")

# Setup samplesheet
samples = metadata_df.index.tolist()

# Convert df to dict for tumor/normal pairing
metadata_dict = metadata_df.to_dict('index')
chip_samples = [sample for sample in metadata_dict if metadata_dict[sample]['type'] == "chip"]

# Get the raw fastq files to start with
def get_fastq(wildcards):
    return metadata_df.loc[(wildcards.sample), "r1"]

def get_input_bam(wildcards):
    return "results/bam/{}.bam".format(metadata_df.loc[(wildcards.chip_sample), "input"])
