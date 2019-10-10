import random
import math
import os 

#############
## MAPPING ##
#############

rule star_generate_index_for_genome:
	input:
		fasta=GENOME_FASTA,
		annotation=GENOME_GTF
	output:
		REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab"
	threads: 4
	conda:
		config["conda_envs"] + "/star.yml"
	shell:
		"if [ -d {config[sample_data_dir]}/STAR_tmp_Index ]; then rm -r {config[sample_data_dir]}/STAR_tmp_Index; fi "
		"&& STAR --outTmpDir {config[sample_data_dir]}/STAR_tmp_Index --runThreadN {threads} --runMode genomeGenerate --genomeDir {REF_GENOME_DIR} "
		"--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.annotation}"

rule bowtie2_generate_index_for_genome:
	input:
		fasta=GENOME_FASTA
	output:
		REF_GENOME_DIR + "/bowtie2index.1.bt2"
	threads: 4
	conda:
		config["conda_envs"] + "/bowtie2.yml"
	params:
		basename=REF_GENOME_DIR + "/bowtie2index"
	shell:
		"bowtie2-build --seed 123 --threads {threads} {input.fasta} {params.basename}"

MAPPING_INDIR = ""
if ( demultiplexed == "yes" ): 
	MAPPING_INDIR = REMOVE_TAIL_OUTDIR
else:
	MAPPING_INDIR = DEMULTI_OUTDIR

if ( mapper == "STAR" ):

	if ( PROTOCOL == "eCLIP" ):

		rule star:
			input:
				REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab",
				first_read=MAPPING_INDIR + "/{sample}_{replicate}_r2_trimmed.fastq",
				second_read=MAPPING_INDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
			output:
				log=MAPPING_OUTDIR + "/{sample}_{replicate}.txt",
				bam=MAPPING_OUTDIR + "/{sample}_{replicate}.bam",
				final=MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out"
			threads: 4
			conda:
				config["conda_envs"] + "/star.yml"
			params:
				output_folder=MAPPING_OUTDIR + "/{sample}_{replicate}"	
			shell:
				"if [ ! -d {MAPPING_OUTDIR} ]; then mkdir {MAPPING_OUTDIR}; fi "
				"&& TIME=$(date +%N) "
				"&& echo {config[star_all]} >> {file_tool_params}"
				"&& echo {config[star_indi]} >> {file_tool_params}"
				"&& echo {config[star_frag]} >> {file_tool_params}"
				"&& if [ -d {config[sample_data_dir]}/STAR_tmp_$TIME ]; then rm -r {config[sample_data_dir]}/STAR_tmp_$TIME; fi "
				"&& STAR --runThreadN {threads} --genomeLoad NoSharedMemory --genomeDir {REF_GENOME_DIR} "   
				"--readFilesIn {input.first_read} {input.second_read} --outTmpDir {config[sample_data_dir]}/STAR_tmp_$TIME  --outFileNamePrefix {params.output_folder}_ "  
				"{config[star_all]} {config[star_indi]} {config[star_frag]} > {output.log} "
				"&& mv {params.output_folder}_Aligned.sortedByCoord.out.bam {output.bam} "
				"&& rm -r {config[sample_data_dir]}/STAR_tmp_$TIME"
	else:

		rule star:
			input:
				REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab",
				first_read=MAPPING_INDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
				second_read=MAPPING_INDIR + "/{sample}_{replicate}_r2_trimmed.fastq",
			output:
				log=MAPPING_OUTDIR + "/{sample}_{replicate}.txt",
				bam=MAPPING_OUTDIR + "/{sample}_{replicate}.bam",
				final=MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out"
			threads: 4
			conda:
				config["conda_envs"] + "/star.yml"
			params:
				output_folder=MAPPING_OUTDIR + "/{sample}_{replicate}"	
			shell:
				"if [ ! -d {MAPPING_OUTDIR} ]; then mkdir {MAPPING_OUTDIR}; fi "
				"&& TIME=$(date +%N) "
				"&& echo {config[star_all]} >> {file_tool_params}"
				"&& echo {config[star_indi]} >> {file_tool_params}"
				"&& echo {config[star_frag]} >> {file_tool_params}"
				"&& if [ -d {config[sample_data_dir]}/STAR_tmp_$TIME ]; then rm -r {config[sample_data_dir]}/STAR_tmp_$TIME; fi "
				"&& STAR --runThreadN {threads} --genomeLoad NoSharedMemory --genomeDir {REF_GENOME_DIR} "   
				"--readFilesIn {input.first_read} {input.second_read} --outTmpDir {config[sample_data_dir]}/STAR_tmp_$TIME  --outFileNamePrefix {params.output_folder}_ "  
				"{config[star_all]} {config[star_indi]} {config[star_frag]} > {output.log} "
				"&& mv {params.output_folder}_Aligned.sortedByCoord.out.bam {output.bam} "
				"&& rm -r {config[sample_data_dir]}/STAR_tmp_$TIME"

if ( mapper == "Bowtie2"):

	if ( PROTOCOL == "eCLIP" ):

		rule bowtie2:
			input:
				REF_GENOME_DIR + "/bowtie2index.1.bt2",
				first_read=MAPPING_INDIR + "/{sample}_{replicate}_r2_trimmed.fastq",
				second_read=MAPPING_INDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
			output:
				sam=MAPPING_OUTDIR + "/{sample}_{replicate}.sam",
				log=MAPPING_OUTDIR + "/{sample}_{replicate}_log.txt"
			threads: 4
			conda:
				config["conda_envs"] + "/bowtie2.yml" # samtools already included
			params:
				basename=REF_GENOME_DIR + "/bowtie2index"
			shell:
				"if [ ! -d {MAPPING_OUTDIR} ]; then mkdir {MAPPING_OUTDIR}; fi "
				"&& echo {config[bowtie2]} >> {file_tool_params}"
				"&& bowtie2 -p {threads} -x {params.basename} -1 {input.first_read} -2 {input.second_read} {config[bowtie2]} -S {output.sam} 2> {output.log}"

	else:

		rule bowtie2:
			input:
				REF_GENOME_DIR + "/bowtie2index.1.bt2",
				first_read=MAPPING_INDIR + "/{sample}_{replicate}_r1_trimmed.fastq",
				second_read=MAPPING_INDIR + "/{sample}_{replicate}_r2_trimmed.fastq",
			output:
				sam=MAPPING_OUTDIR + "/{sample}_{replicate}.sam",
				log=MAPPING_OUTDIR + "/{sample}_{replicate}_log.txt"
			threads: 4
			conda:
				config["conda_envs"] + "/bowtie2.yml" # samtools already included
			params:
				basename=REF_GENOME_DIR + "/bowtie2index"
			shell:
				"if [ ! -d {MAPPING_OUTDIR} ]; then mkdir {MAPPING_OUTDIR}; fi "
				"&& echo {config[bowtie2]} >> {file_tool_params}"
				"&& bowtie2 -p {threads} -x {params.basename} -1 {input.first_read} -2 {input.second_read} {config[bowtie2]} -S {output.sam} 2> {output.log}"

	rule bowtie2_convert:
		input:
			MAPPING_OUTDIR + "/{sample}_{replicate}.sam"
		output:
			MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
		threads: 2
		conda:
			config["conda_envs"] + "/samtools.yml"
		shell:
			"samtools sort -O bam -o {output} {input}"

rule indexing:
	input:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai"
	threads: 2
	conda:
		config["conda_envs"] + "/samtools.yml"
	shell:
		"samtools index {input}"
