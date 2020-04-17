rule subtract_input:
	input:
		trt_bdg = rules.macs2.output.trt_bdg,
		cnt_bdg = rules.macs2.output.cnt_bdg
	output:
		temp("results/peaks/{chip_sample}/{chip_sample}_subtract.bdg")
	conda:
		"envs/macs2.yaml"
	threads: threads_low
	shell:
		'macs2 bdgcmp '
		'--tfile {input.trt_bdg} '
		'--cfile {input.cnt_bdg} '
		'--ofile {output} '
		'--method subtract '

rule fix_neg_vals:
	input:
		rules.subtract_input.output
	output:
		temp("results/peaks/{chip_sample}/{chip_sample}_input_minus.bdg")
	conda:
		"envs/gawk.yaml"
	threads: threads_low
	shell:
		'awk \'{{OFS="\\t" ; if($4 < 0) {{print $1,$2,$3,0}} else {{print $0}}}}\' '
		'{input} > {output}'


rule bdg2bw:
	input:
		bdg = "results/peaks/{chip_sample}/{chip_sample}_{bdg_type}.bdg",
		genome_sizes = config['genome_sizes']
	output:
		bw = "results/bw/{chip_sample}-{bdg_type}.bw",
		clip = temp("results/peaks/{chip_sample}/{chip_sample}_{bdg_type}.bdg.clip"),
		sort_clip = temp("results/peaks/{chip_sample}/{chip_sample}_{bdg_type}.bdg.sort.clip"),
	conda:
		"envs/bdg2bw.yaml"
	threads: threads_mid
	shell:
		"""
		bedtools slop -i {input.bdg} -g {input.genome_sizes} -b 0 | bedClip stdin {input.genome_sizes} {output.clip}
		LC_COLLATE=C sort -k1,1 -k2,2n {output.clip} > {output.sort_clip}
		bedGraphToBigWig {output.sort_clip} {input.genome_sizes} {output.bw}
		"""
