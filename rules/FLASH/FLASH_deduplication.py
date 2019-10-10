import random
import math
import os

###################
## DEDUPLICATION ##
###################

rule fastqc_before_dedup:
    input:
    	PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
    output:
    	FASTQC_BEFORE_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html"
    threads: 2
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_BEFORE_DEDUP_OUTDIR} ]; then mkdir {FASTQC_BEFORE_DEDUP_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_BEFORE_DEDUP_OUTDIR}"

rule deduplication:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
	output:
		bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}.bam",
		log=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_log.txt",
		sorted_bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	threads: 2
	conda:
		config["conda_envs"] + "/umi.yml"  # samtools alreads included
	shell:
		"if [ ! -d {DEDUPLICAITON_OUTDIR} ]; then mkdir {DEDUPLICAITON_OUTDIR}; fi "
		"&& echo {config[umi_dedup]} >> {file_tool_params}"
		"&& umi_tools dedup {config[umi_dedup]} -I {input} -S {output.bam} -L {output.log}"
		"&& samtools sort {output.bam} > {output.sorted_bam}"
		"&& samtools index {output.sorted_bam}"

rule sam_name_sort:
	input:
		sorted_bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
	output:
		name_sorted_bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_name_sorted.bam"
	threads: 2
	conda:
		config["conda_envs"] + "/samtools.yml"
	params:
		tmp=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}"
	shell:
		"samtools sort -n -T {params.tmp} -o {output.name_sorted_bam} {input.sorted_bam}"
