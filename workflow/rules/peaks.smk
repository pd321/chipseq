rule macs2:
	input:
		treatment = "results/bam/{chip_sample}.sorted.remdup.nonblklst.filt.resort.bam",
		control = get_input_bam
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
		extra = get_macs2_params
	wrapper:
		"v1.24.0/bio/macs2/callpeak"
