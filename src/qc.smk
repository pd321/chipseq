rule fastqc_se:
	input:
		lambda wildcards: metadata_dict[wildcards.se_sample]["r1"]
	output:
		html="results/qc/fastqc/{se_sample}_fastqc.html",
		zip="results/qc/fastqc/{se_sample}_fastqc.zip"
	threads: threads_mid
	wrapper:
		"0.38.0/bio/fastqc"

rule fastqc_pe:
	input:
		lambda wildcards: metadata_dict[wildcards.pe_sample][wildcards.group]
	output:
		html="results/qc/fastqc/{pe_sample}_{group}_fastqc.html",
		zip="results/qc/fastqc/{pe_sample}_{group}_fastqc.zip"
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
		expand(["results/qc/fastqc/{se_sample}_fastqc.html", 
			"results/qc/fastqc/{pe_sample}_{group}_fastqc.html", 
		 "results/qc/flagstat/{sample}.txt"], sample = samples, 
		 se_sample = se_samples, pe_sample = pe_samples, group = ["r1", "r2"],)
	output:
		report("results/qc/multiqc/multiqc_report.html", caption ="report/multiqc.rst", category = "Quality control")
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
		'logs results 2>&1 | tee {log}'
