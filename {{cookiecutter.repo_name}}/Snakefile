include: "src/common.smk"
include: "src/addons.smk"
include: "src/align.smk"
include: "src/peaks.smk"
include: "src/bw.smk"
include: "src/se.smk"
include: "src/qc.smk"


peaks = expand("results/peaks/{chip_sample}/{chip_sample}_peaks.filt.bed", chip_sample = chip_samples)
bw = expand("results/bw/{chip_sample}-{bdg_type}.bw", chip_sample = chip_samples, bdg_type = bdg_types)
qc = ["results/qc/multiqc/multiqc_report.html"]
out_files = peaks + bw + qc

if config['addons']['se']:
	out_files += expand("results/se/{chip_sample}/{chip_sample}_peaks_AllEnhancers_ENHANCER_TO_GENE.txt", chip_sample = chip_samples)

rule all:
	input:
		out_files
		
		
onsuccess:
	shell("if command -v telegram-notify; then telegram-notify --success --text \'snakemake:chipseq:{} completed\'; fi".format(run_name))

onerror:
	shell("if command -v telegram-notify; then telegram-notify --error --text \'snakemake:chipseq:{} failed\'; fi".format(run_name))
