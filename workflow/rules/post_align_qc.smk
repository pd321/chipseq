rule remdup:
	input:
		bams = rules.sortbam.output
	output:
		bam = temp("results/bam/{sample}.sorted.remdup.bam"),
		bai = temp("results/bam/{sample}.sorted.remdup.bai"),
		metrics = "results/qc/remdup/{sample}.metrics.txt"
	threads: 
		config_threads
	log:
		"logs/qc/remdup/{sample}.log"
	params:
		extra="--REMOVE_DUPLICATES true --CREATE_INDEX true",
	wrapper:
		"v1.23.5-48-gf27313f0/bio/picard/markduplicates"

rule remove_blacklist_reads:
	input:
		left = rules.remdup.output,
		right = config['blklist_regions']
	output:
		temp("results/bam/{sample}.sorted.remdup.nonblklst.bam")
	log:
		"logs/qc/remove_blacklist_reads/{sample}.log"
	params:
		extra="-v"
	wrapper:
		"v1.24.0/bio/bedtools/intersect"

rule filter_bam:
	input:
		rules.remove_blacklist_reads.output
	output:
		bam = temp("results/bam/{sample}.sorted.remdup.nonblklst.filt.bam"),
		bai = temp("results/bam/{sample}.sorted.remdup.nonblklst.filt.bai")
	log:
		"logs/qc/filter_bam/{sample}.log"
	threads: 
		config_threads
	params:
		extra="-F 1804 -f 2"
	wrapper:
		"v1.23.5-48-gf27313f0/bio/samtools/view"

rule re_sort_bam:
	input: rules.filter_bam.output.bam
	output: 
		bam = "results/bam/{sample}.sorted.remdup.nonblklst.filt.resort.bam",
		bai = "results/bam/{sample}.sorted.remdup.nonblklst.filt.resort.bai"
	log:
		"logs/samtools/resortbam/{sample}.log"
	threads: config_threads
	wrapper:
		"v1.23.5-48-gf27313f0/bio/samtools/sort"

rule flagstat:
	input:
		rules.re_sort_bam.output
	output:
		"results/qc/flagstat/{sample}.txt"
	log:
		"logs/qc/flagstat/{sample}.log",
	threads: config_threads
	wrapper:
		"v1.24.0/bio/samtools/flagstat"

rule multiqc:
	input:
		expand(["results/qc/fastqc/{se_sample}_fastqc.html", 
			"results/qc/fastqc/{pe_sample}_{group}_fastqc.html", 
		 "results/qc/flagstat/{sample}.txt"], sample = samples, 
		 se_sample = se_samples, pe_sample = pe_samples, group = ["r1", "r2"],)
	output:
		report("results/qc/multiqc/multiqc_report.html", caption ="report/multiqc.rst", category = "Quality control")
	conda:
		"../envs/multiqc.yaml"
	threads: 1
	log:
		"logs/qc/multiqc/multiqc.log"
	shell:
		'multiqc '
		'--force '
		'--outdir results/qc/multiqc '
		'--zip-data-dir '
		'logs results 2>&1 | tee {log}'
