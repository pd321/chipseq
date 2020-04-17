rule bowtie:
	input: 
		rules.trimgalore.output if config['addons']['trimgalore'] else get_fastq
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

rule remdup:
	input:
		rules.bowtie.output
	output:
		bam = "results/bam/{sample}.bam",
		bai = "results/bam/{sample}.bai",
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
		'REMOVE_DUPLICATES= true '
		'METRICS_FILE={output.metrics} 2>&1 | tee {log}'
