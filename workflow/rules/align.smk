rule trimgalore_se:
	input: get_fastq_se
	output: temp("results/bam/{sample}_trimmed.fq.gz")
	conda: "envs/trimgalore.yaml"
	log: "logs/trimgalore/{sample}.log"
	params:
		quality = config['trimgalore']['quality'],
		stringency = config['trimgalore']['stringency'],
		e = config['trimgalore']['e'],
		#trimgalore does not recommend using more than 4 threads
		threads_actual = threads_mid if threads_mid < 4 else 4
	# setting it to mid ensures not no many trimgalore jobs launched
	threads: threads_mid
	shell:
		'trim_galore '
		'--quality {params.quality} '
		'--stringency {params.stringency} '
		'-e {params.e} '
		'--gzip '
		'--output_dir results/bam/ '
		'--cores {params.threads_actual} '
		'--basename {wildcards.sample} '
		'--no_report_file '
		'{input} 2>&1 | tee {log}'

rule trimgalore_pe:
	input:
		get_fastq_pe
	output:
		r1 = temp("results/bam/{sample}_val_1.fq.gz"),
		r2 = temp("results/bam/{sample}_val_2.fq.gz")
	conda:
		"envs/trimgalore.yaml"
	log:
		"logs/trimgalore/{sample}.log"
	params:
		quality = config['trimgalore']['quality'],
		stringency = config['trimgalore']['stringency'],
		e = config['trimgalore']['e'],
		#trimgalore does not recommend using more than 4 threads
		threads_actual = threads_mid if threads_mid < 4 else 4
	# setting it to mid ensures not no many trimgalore jobs launched
	threads: threads_mid
	shell:
		'trim_galore '
		'--quality {params.quality} '
		'--stringency {params.stringency} '
		'-e {params.e} '
		'--gzip '
		'--output_dir results/bam/ '
		'--cores {params.threads_actual} '
		'--basename {wildcards.sample} '
		'--paired --no_report_file '
		'{input[0]} {input[1]} 2>&1 | tee {log}'

rule bowtie:
	input: 
		rules.trimgalore_se.output
	output: 
		temp("results/bam/{sample}_bowtie.bam")
	conda: 
		'envs/bowtie.yaml'
	threads: 
		threads_mid
	log: 
		'logs/bowtie/{sample}.log'
	params:
		idx = config['bowtie']['idx'],
		max_multiple_aln = config['bowtie']['max_multiple_aln'],
		max_aln_per_read = config['bowtie']['max_aln_per_read'],
		chunkmbs = config['bowtie']['chunkmbs'],
		platform = config['platform'],
		labname = '--sam-RG CN:{}'.format(config['labname']) if config['labname'] else '',
		date = time.strftime("%Y-%m-%d")
	shell:
		'zcat {input} | '
		'bowtie '
		'-m {params.max_multiple_aln} '
		'-k {params.max_aln_per_read} '
		'--best --strata '
		'--sam '
		'--chunkmbs {params.chunkmbs} '
		'--time '
		'--threads {threads} '
		'--sam-RG SM:{wildcards.sample} '
		'--sam-RG LB:{wildcards.sample} '
		'--sam-RG ID:{wildcards.sample} '
		'--sam-RG PL:{params.platform} '
		'{params.labname} '
		'--sam-RG DT:{params.date} '
		'{params.idx} - 2> {log} | '
		'samtools view -@ {threads} -bS - | '
		'samtools sort -@ {threads} -o {output} -'

rule bowtie2:
	input:
		r1 = rules.trimgalore_pe.output.r1,
		r2 = rules.trimgalore_pe.output.r2
	output:
		temp("results/bam/{sample}_bowtie2.bam")
	conda:
		"envs/bowtie2.yaml"
	log:
		"logs/bowtie2/{sample}.log"
	params:	
		idx = config['bowtie2']['idx'],
		maxins = config['bowtie2']['maxins']
	threads: threads_high
	shell:
		'bowtie2 '
		'--very-sensitive '
		'--maxins {params.maxins} '
		'--no-mixed '
		'--no-discordant '
		'--time '
		'--threads {threads} '
		'-x {params.idx} '
		'-1 {input[0]} '
		'-2 {input[1]} '
		'2> {log} | '
		'samtools view '
		'--threads {threads} '
		'-b -| '
		'samtools sort '
		'--threads {threads} '
		'-o {output}'

rule remdup:
	input:
		get_remdup_input
	output:
		bam = "results/bam/{sample}_remdup.bam",
		bai = "results/bam/{sample}_remdup.bai",
		metrics = "results/qc/remdup/{sample}.metrics.txt"
	conda:
		"envs/picard.yaml"
	threads: 
		threads_mid
	log:
		"logs/remdup/{sample}.log"
	shell:
		'picard MarkDuplicates '
		'CREATE_INDEX=true '
		'INPUT={input} '
		'OUTPUT={output.bam} '
		'REMOVE_DUPLICATES=true '
		'METRICS_FILE={output.metrics} 2>&1 | tee {log}'
