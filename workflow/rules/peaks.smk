rule macs2:
	input:
		trt_bam = "results/bam/{chip_sample}_remdup.bam",
		cnt_bam = get_input_bam
	output:
		peaks = temp("results/peaks/{chip_sample}/{chip_sample}_peaks.narrowPeak")
	conda:
		"envs/macs2.yaml"
	log:
		"logs/macs2/{chip_sample}.log"
	threads: threads_mid
	params:	
		qvalue = config['macs2']['qvalue'],
		gsize = config['macs2']['gsize']
	shell:
		'macs2 callpeak '
		'--treatment {input.trt_bam} '
		'--control {input.cnt_bam} '
		'--format AUTO '
		'--gsize {params.gsize} '
		'--keep-dup all '
		'--outdir results/peaks/{wildcards.chip_sample} '
		'--nomodel '
		'--extsize 200 '
		'--qvalue {params.qvalue} '
		'--name {wildcards.chip_sample} 2> {log}'

rule blklist_filt:
	input:
		peaks = rules.macs2.output.peaks,
		blklist_regions = config['blklist_regions']
	output:
		"results/peaks/{chip_sample}/{chip_sample}_peaks.filt.narrowPeak"
	conda:
		"envs/bedtools.yaml"
	threads: threads_low
	shell:
		'bedtools intersect '
		'-v -a {input.peaks} -b {input.blklist_regions} > {output}'
