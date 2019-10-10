import random
import math
import os

rule fastqc_beginning:
    input:
    	RENAMING + "/{sample}_{replicate}_{pair}.fastq"
    output:
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.html",
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.zip"
    threads: 2
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_BEG_OUTDIR} ]; then mkdir {FASTQC_BEG_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_BEG_OUTDIR}"

# rule cutadapt_first_read_clip:
# 	input:
# 		first=RENAMING + "/{sample}_{replicate}_r1.fastq",
# 		second=RENAMING + "/{sample}_{replicate}_r2.fastq"
# 	output:
# 		seq_first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastq",
# 		seq_second=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r2.fastq",
# 		log=CUTADAPT_OUTDIR + "/{sample}_{replicate}.txt"
# 	threads: 2
# 	conda:
# 		config["conda_envs"] + "/cutadapt.yml"
# 	shell:
# 		"if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; fi"
# 		"&& echo {config[cutadapt]} >> {file_tool_params}"
# 		"&& cutadapt -j {threads} {config[cutadapt]} "
# 		"--paired-output={output.seq_second} --output={output.seq_first} {input.first} {input.second} > {output.log}"

if ( config["remove_tail"] != "0" ):
	rule remove_tail:
		input:
			first=RENAMING + "/{sample}_{replicate}_r1.fastq",
			second=RENAMING + "/{sample}_{replicate}_r2.fastq"
		output:
			first=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
			second=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastq"
		threads: 2
		conda:
			config["conda_envs"] + "/bctools.yml"
		shell:
			"if [ ! -d {REMOVE_TAIL_OUTDIR} ]; then mkdir {REMOVE_TAIL_OUTDIR}; fi"
			"&& echo {config[remove_tail]} >> {file_tool_params}"
			"&& python {config[bctools]}/remove_tail.py {input.first} {config[remove_tail]} > {output.first}"
			"&& cp {input.second} {output.second}"
else:
	rule remove_tail:
		input:
			first=RENAMING + "/{sample}_{replicate}_r1.fastq",
			second=RENAMING + "/{sample}_{replicate}_r2.fastq"
		output:
			first=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
			second=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastq"
		threads: 2
		conda:
			config["conda_envs"] + "/bctools.yml"
		shell:
			"if [ ! -d {REMOVE_TAIL_OUTDIR} ]; then mkdir {REMOVE_TAIL_OUTDIR}; fi"
			"&& echo {config[remove_tail]} >> {file_tool_params}"
			"&& cp {input.first} {output.first}"
			"&& cp {input.second} {output.second}"

rule fastqc_after_adapter_removal:
    input:
    	REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastq"
    output:
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html",
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.zip"
    threads:
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_ADAPT_OUTDIR} ]; then mkdir {FASTQC_ADAPT_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_ADAPT_OUTDIR}"
