#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import pandas as pd


def main(args):
	logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)
	logging.info("Will include samples from the following dirs: {dirs}".format(dirs=",".join(args.sample_dirs)))

	metadata = {}

	for sample_dir in args.sample_dirs:
		files = glob.glob(sample_dir + '/**/*' + args.read_extension, recursive=True)
		for file in files:
			sample_name = os.path.basename(os.path.dirname(file))
			metadata.setdefault(sample_name, {})
			logging.info("Using the {file} as input for {sample}".format(file=file, sample=sample_name))
			metadata[sample_name]['r1'] = file

	metadata_df = pd.DataFrame(metadata).transpose()
	metadata_df['type'] = ""
	metadata_df['input'] = ""
	metadata_df.to_csv("metadata.tsv", sep="\t", index_label="sample_name")


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Script to generate metadata.tsv")
	parser.add_argument('-d', '--dirs', dest='sample_dirs', action='append')
	parser.add_argument('-f', '--read1', dest='read_extension', default='_1.fq.gz',
						choices=['_1.fq.gz', '_R1_001.fastq.gz', 'fastq.gz', 'fq.gz'])
	cmd_args = parser.parse_args()
	main(cmd_args)
