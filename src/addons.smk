rule trimgalore:
	input: get_fastq
	output: temp("results/bam/{sample}_1_trimmed.fq.gz")
	conda: "envs/trimgalore.yaml"
	log: "logs/trimgalore/{sample}.log"
	params:
		quality = config['trimgalore']['quality'],
		stringency = config['trimgalore']['stringency'],
		e = config['trimgalore']['e']
	# trimgalore does not recommend using more than 4 threads
	threads: threads_mid if threads_mid < 4 else 4
	shell:
		'trim_galore '
		'--quality {params.quality} '
		'--stringency {params.stringency} '
		'-e {params.e} '
		'--gzip '
		'--output_dir results/bam/ '
		'--cores {threads} '
		'--basename {wildcards.sample} '
		'--no_report_file '
		'{input} 2>&1 | tee {log}'
