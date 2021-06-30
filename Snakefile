include: "src/common.smk"
include: "src/align.smk"
include: "src/peaks.smk"
include: "src/bw.smk"
include: "src/se.smk"
include: "src/qc.smk"
include: "src/addons.smk"


peaks = expand("results/peaks/{chip_sample}/{chip_sample}_peaks.filt.bed", chip_sample = chip_samples)
bw = expand("results/bw/{sample}.bw", sample = samples)
qc = ["results/qc/multiqc/multiqc_report.html"]
out_files = peaks + bw + qc

if config['addons']['se']:
	out_files += expand("results/se/{chip_sample}/{chip_sample}_peaks_AllEnhancers_ENHANCER_TO_GENE.txt", chip_sample = chip_samples)

if config['addons']['homer_annot']:
	out_files += expand("results/peaks/{chip_sample}/{chip_sample}_peaks_annot.xls", chip_sample = chip_samples)

if config['addons']['homer_motifs']:
	out_files += expand("results/motif/{chip_sample}/homerResults.html", chip_sample = chip_samples)

rule all:
	input:
		out_files
		
		
onsuccess:
	shell("if command -v telegram-notify; then telegram-notify --success --text \'snakemake:chipseq:{} completed\'; fi".format(run_name))

onerror:
	shell("if command -v telegram-notify; then telegram-notify --error --text \'snakemake:chipseq:{} failed\'; fi".format(run_name))
