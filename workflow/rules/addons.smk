rule homer_annot:
	"""
	Annotate the peaks with homer
	annotatePeaks.pl 02macs/A78/A78_peaks.broadPeak hg19 > 05homerAnnot/raw/A78_homer.xls
	"""
	input:
		rules.blklist_filt.output
	output:
		temp("results/peaks/{chip_sample}/{chip_sample}_peaks.filt.annot.bed")
	threads: threads_low
	params:
		genome = config['genome']
	shell:
		'annotatePeaks.pl {input} {params.genome} > {output}'

rule clean_merge_homer:
	"""
	This will clean up the homer output and add in details of the narrowpeak
	"""
	input:
		rules.homer_annot.output
	output:
		"results/peaks/{chip_sample}/{chip_sample}_peaks_annot.xls"
	threads: threads_low
	log: 'logs/addons/clean_merge_homer/{chip_sample}.log'
	script:
		"R/homer_clean.R"

rule findmotifs:
    input:
        rules.blklist_filt.output
    output:
        known = "results/motif/{chip_sample}/knownResults.html",
        denovo = "results/motif/{chip_sample}/homerResults.html"
    params:
        genome = config['genome']
    threads: threads_mid
    shell:
        'findMotifsGenome.pl {input} {params.genome} results/motif/{wildcards.chip_sample}'
