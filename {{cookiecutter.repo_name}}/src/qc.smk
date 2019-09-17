rule fastqc:
	input:
		get_fastq
	output:
		html="results/qc/fastqc/{sample}_fastqc.html",
		zip="results/qc/fastqc/{sample}_fastqc.zip"
	threads: threads_mid
	wrapper:
		"0.38.0/bio/fastqc"

rule flagstat:
	input:
		rules.remdup.output.bam
	output:
		"results/qc/flagstat/{sample}.txt"
	conda:
		"envs/samtools.yaml"
	threads: threads_mid
	shell:
		'samtools flagstat '
		'--threads {threads} '
		'{input} > {output}'

rule multiqc:
	input:
		expand(["results/qc/fastqc/{sample}_fastqc.html", "results/qc/flagstat/{sample}.txt", 
			"results/qc/remdup/{sample}.metrics.txt"], sample = samples)
	output:
		report("results/qc/multiqc/multiqc_report.html", caption="report/multiqc.rst", category="Quality control")
	conda:
		"envs/multiqc.yaml"
	threads: threads_high
	log:
		"logs/multiqc/multiqc.log"
	shell:
		'multiqc '
		'--force '
		'--outdir results/qc/multiqc '
		'--zip-data-dir '
		'. 2>&1 | tee {log}'
