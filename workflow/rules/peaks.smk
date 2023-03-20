rule macs2:
	input:
		trt_bam = "results/bam/{chip_sample}.sorted.remdup.nonblklst.filt.resort.bam",
		cnt_bam = get_input_bam
	output:
		multiext("results/peaks/{chip_sample}/{chip_sample}",
				 "_peaks.xls",
				 "_peaks.narrowPeak",
				 "_summits.bed"
				 )
	log:
		"logs/macs2/{chip_sample}.log"
	threads: config_threads
	params:	
		extra = "--format BAMPE --gsize {macs2_genome} --keep-dup all --nomodel".format(macs2_genome = config['macs2']['gsize'])
	wrapper:
		"v1.24.0/bio/macs2/callpeak"
