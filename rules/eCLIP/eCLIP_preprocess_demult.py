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

rule cutadapt_first_read_clip:
	input:
		first=RENAMING + "/{sample}_{replicate}_r1.fastq",
		second=RENAMING + "/{sample}_{replicate}_r2.fastq"
	output:
		seq_first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastq",
		seq_second=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r2.fastq",
		log=CUTADAPT_OUTDIR + "/{sample}_{replicate}.txt"
	threads: 2
	conda:
		config["conda_envs"] + "/cutadapt.yml"
	shell:
		"if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; fi"
		"&& echo {config[cutadapt]} >> {file_tool_params}"
		"&& cutadapt -j {threads} {config[cutadapt]} "
		"--paired-output={output.seq_second} --output={output.seq_first} {input.first} {input.second} > {output.log}"

rule remove_tail:
	input:
		first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastq",
		second=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r2.fastq"
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

rule fastqc_after_adapter_removal:
    input:
    	REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastq"
    output:
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html",
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.zip"
    threads: 2
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_ADAPT_OUTDIR} ]; then mkdir {FASTQC_ADAPT_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_ADAPT_OUTDIR}"

# Activate if you need to remove UMIs and barcodes
rule got_umis:
	input:
		first=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
		second=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastq"
	output:
		first=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
		second=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastq",
		log=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_log.txt"
	threads: 2
	conda:
		config["conda_envs"] + "/umi.yml"
	shell:
		"if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; fi "
		"&& echo {config[got_umis_1]} >> {file_tool_params}"
		"&& umi_tools extract {config[got_umis_1]} -I {input.second} -S {output.second} --read2-in {input.first} --read2-out {output.first} -L {output.log}"

rule demultiplex:
	input:
		first=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
		second=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_r2_trimmed.fastq",
		barcodes=config["barcodes"]
	output:
		diag=DEMULTI_OUTDIR + "/{sample}_{replicate}_diag.log",
		files=[DEMULTI_OUTDIR + "/{sample}_{replicate}" + x for x in BARCODE_FILES]
	threads: 2
	conda:
		config["conda_envs"] + "/je_demultiplex.yml"
	params:
		outdir=DEMULTI_OUTDIR
	shell:
		"if [ ! -d {DEMULTI_OUTDIR} ]; then mkdir {DEMULTI_OUTDIR}; fi "
		"&& echo {config[demultiplex]} >> {file_tool_params}"
		"&& je demultiplex F1={input.first} F2={input.second} BARCODE_FILE={input.barcodes} {config[demultiplex]} "
		"OUTPUT_DIR={params.outdir} BARCODE_DIAG_FILE={output.diag} METRICS_FILE_NAME={params.outdir}/metric.txt"

rule renaming_demultiplex:
	input:
		#BARCODE_FILES
		[DEMULTI_OUTDIR + "/" + MULTIPLEX_SAMPLE_NAME + "_rep1_" + x for x in BARCODE_FILES]
	output:
		BARCODE_NEWFILES
	run:
		for i in range(0, len(BARCODE_FILES)):
			shell("cp " + input[i] + " " + output[i])
