rule rose_main:
	input:
		trt_bam = "results/bam/{chip_sample}.bam",
		cnt_bam = get_input_bam,
		peaks = rules.blklist_filt.output
	output:
		enhancers_table = "results/se/{chip_sample}/{chip_sample}_peaks_AllEnhancers.table.txt",
		enhancers_plot = "results/se/{chip_sample}/{chip_sample}_Plot_points.png"
	threads: threads_mid
	log: "logs/se/{chip_sample}_rose_main.log"
	params:
		rose_dir = config['se']['rose_dir'],
		genome = config['genome'],
		tss_dist = config['se']['tss_dist'],
		rundir = os.getcwd()
	shell:
		'cd {params.rose_dir} && '
		'python2 ROSE_main.py '
		'--i \'{params.rundir}/{input.peaks}\' '
		'--rankby \'{params.rundir}/{input.trt_bam}\' '
		'--control \'{params.rundir}/{input.cnt_bam}\' '
		'--genome {params.genome} '
		'--out \'{params.rundir}/results/se/{wildcards.chip_sample}\' '
		'--tss {params.tss_dist} > {params.rundir}/{log}'

rule rose_annot:
	input:
		rules.rose_main.output.enhancers_table
	output:
		"results/se/{chip_sample}/{chip_sample}_peaks_AllEnhancers_ENHANCER_TO_GENE.txt"
	log: "logs/se/{chip_sample}_rose_annot.log"
	threads: threads_low
	params:
		rose_dir = config['se']['rose_dir'],
		genome = config['genome'],
		rundir = os.getcwd()
	shell:
		'cd {params.rose_dir} && '
		'python2 ROSE_geneMapper.py '
		'-i \'{params.rundir}/{input}\' '
		'--genome {params.genome} > {params.rundir}/{log}'
