import random
import math
import os 

#################
## PEAKCALLING ##
#################

# BE CAREFUL when using other protocols then you might need to define a different bitflag!!
rule mate_reads_fitlering:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	output:
		bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}.bam",
		bam_pos=MATEFILTER_OUTDIR + "/{sample}_{replicate}_pos.bam",
		bam_neg=MATEFILTER_OUTDIR + "/{sample}_{replicate}_neg.bam",
		bam_bai_pos=MATEFILTER_OUTDIR + "/{sample}_{replicate}_pos.bam.bai",
		bam_bai_neg=MATEFILTER_OUTDIR + "/{sample}_{replicate}_neg.bam.bai",
		sorted_bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam",
		sorted_bam_bai=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam.bai"
	threads: 2
	conda:
		config["conda_envs"] + "/samtools.yml"
	shell:
		"if [ ! -d {MATEFILTER_OUTDIR} ]; then mkdir {MATEFILTER_OUTDIR}; fi"
		"&& samtools view -b -f 0x0080 {input} > {output.bam}"
		"&& samtools view -b -F 0x0010 {output.bam} > {output.bam_pos}"
		"&& samtools view -b -f 0x0010 {output.bam} > {output.bam_neg}"
		"&& samtools sort {output.bam} > {output.sorted_bam}"
		"&& samtools index {output.sorted_bam}"
		"&& samtools index {output.bam_pos}"
		"&& samtools index {output.bam_neg}"

if ( control == "yes" ):

	if ( peakcaller == "PanPeaker" ):
	rule panpeaker:
			input:
		    	clip_bam=expand(MATEFILTER_OUTDIR + "/{sample}.bam", sample=REPLICATES_CLIP),
		    	clip_bai=expand(MATEFILTER_OUTDIR + "/{sample}.bai", sample=REPLICATES_CLIP),
		    	ctl_bam=expand(MATEFILTER_OUTDIR + "/{sample}.bam", sample=REPLICATES_CONTROL),
		    	ctl_bai=expand(MATEFILTER_OUTDIR + "/{sample}.bai", sample=REPLICATES_CONTROL),
		    	genome_fasta=GENOME_FASTA,
		    	chr_sizes=GENOME_SIZES
			output:
				PEAKCALLING_OUTDIR + "/robust_peaks.bed"
			threads: 10 
			shell:
				"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
				"&& echo {config[panpeaker]} >> {file_tool_params}"
				"&& source activate panpeaker"
				"&& python3 {config[panpeaker_script]}/Panpeaker.py {config[panpeaker]} -nt {threads} "
				"-i {input.clip_bam} -b {input.clip_bai} -c {input.ctl_bam} -k {input.ctl_bai} -o {PEAKCALLING_OUTDIR} -g {input.genome_fasta} "
				"--chr_sizes {input.chr_sizes} "
				"&& source deactivate"

	if ( peakcaller == "PureCLIP" ):
		rule pureclip:
			input:
				experiment=MATEFILTER_OUTDIR + "/{sample_exp}_{replicate_exp}_sorted.bam",
				experiment_bai=MATEFILTER_OUTDIR + "/{sample_exp}_{replicate_exp}_sorted.bam.bai",
				control=MATEFILTER_OUTDIR + "/{sample_ctl}_{replicate_ctl}_sorted.bam",
				control_bai=MATEFILTER_OUTDIR + "/{sample_ctl}_{replicate_ctl}_sorted.bam.bai",
				genome_fasta=GENOME_FASTA
			output:
				crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_crosslinkind_sites.bed",
				binding_regions=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed",
			threads: 10
			params:
				tmp=PEAKCALLING_OUTDIR + "/tmp/",
				parameters=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_parameters.txt"
			conda:
				config["conda_envs"] + "/pureclip.yml"
			shell:
				"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
				"&& echo {config[pureclip]} >> {file_tool_params}"
				"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
				"-ibam {input.control} -ibai {input.control_bai} -or {output.binding_regions} -p {params.parameters} -nt {threads} {config[pureclip]} "

		rule peaks_extend_frontiers:
			input:
				bed=PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_binding_regions.bed",
				genome=GENOME_SIZES
			output:
				PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"
			shell:
				"echo {config[peaks_extend_frontiers]} >> {file_tool_params}"
				"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

	rule find_robust_peaks:
		input:
			expand(PEAKCALLING_OUTDIR + "/{sample_exp}_{replicate_exp}_{sample_ctl}_{replicate_ctl}_peaks_extended.bed", 
				sample_exp=SAMPLES[0], replicate_exp=REP_NAME_CLIP, sample_ctl=SAMPLES[1], replicate_ctl=REP_NAME_CONTROL)
		output:
			ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed"
		threads: 2
		conda:
			config["conda_envs"] + "/bedtools.yml"
		params:
			input_folder=PEAKCALLING_OUTDIR,
			output_folder=ROBUSTPEAKS_OUTDIR
		shell:
			"if [ ! -d {ROBUSTPEAKS_OUTDIR} ]; then mkdir {ROBUSTPEAKS_OUTDIR}; fi"
			"&& {config[find_robust_intersections]}/robust_intersections.sh {params.input_folder} {params.output_folder} bed"

else:

	if ( peakcaller == "Piranha" ):
		rule piranha:
			input:
				MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam"
			output:
				PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks.bed"
			conda:
				config["conda_envs"] + "/piranha.yml"
			threads: 2 
			shell:
				"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
				"&& echo {config[piranha]} >> {file_tool_params}"
				"&& Piranha {config[piranha]} {input} -o {output}"

		rule peaks_extend_frontiers:
			input:
				bed=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks.bed",
				genome=GENOME_SIZES
			output:
				PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"
			shell:
				"echo {config[peaks_extend_frontiers]} >> {file_tool_params}"
				"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

	elif ( peakcaller == "PureCLIP" ):
		rule pureclip:
			input:
				experiment=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam",
				experiment_bai=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam.bai",
				genome_fasta=GENOME_FASTA
			output:
				crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_crosslinkind_sites.bed",
				binding_regions=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed"
			threads: 10
			params:
				tmp=PEAKCALLING_OUTDIR + "/tmp/",
				parameters=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_parameters.txt"
			conda:
				config["conda_envs"] + "/pureclip.yml"
			shell:
				"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
				"&& echo {config[pureclip]} >> {file_tool_params}"
				"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
				"-or {output.binding_regions} -p {params.parameters} -nt {threads} {config[pureclip]} "

		rule peaks_extend_frontiers:
			input:
				bed=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed",
				genome=GENOME_SIZES
			output:
				PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
			threads: 2
			conda:
				config["conda_envs"] + "/bedtools.yml"
			shell:
				"echo {config[peaks_extend_frontiers]} >> {file_tool_params}"
				"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

	rule find_robust_peaks:
		input:
			expand(PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
		output:
			ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed"
		threads: 2
		conda:
			config["conda_envs"] + "/bedtools.yml"
		params:
			input_folder=PEAKCALLING_OUTDIR,
			output_folder=ROBUSTPEAKS_OUTDIR
		shell:
			"if [ ! -d {ROBUSTPEAKS_OUTDIR} ]; then mkdir {ROBUSTPEAKS_OUTDIR}; fi"
			"&& {config[find_robust_intersections]}/robust_intersections.sh {params.input_folder} {params.output_folder} bed"
