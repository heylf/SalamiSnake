
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

MAPPING_INDIR = ""
if ( demultiplexed == "yes" ): 
	MAPPING_INDIR = REMOVE_TAIL_OUTDIR
else:
	MAPPING_INDIR = DEMULTI_OUTDIR

rule star:
	input:
		REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab",
		first_read=MAPPING_INDIR + "/{sample}_{replicate}_r1_trimmed.fastq"
	output:
		log=MAPPING_OUTDIR + "/{sample}_{replicate}.txt",
		bam=MAPPING_OUTDIR + "/{sample}_{replicate}.bam",
		qual=MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out"
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
		"--readFilesIn {input.first_read} --outTmpDir {config[sample_data_dir]}/STAR_tmp_$TIME  --outFileNamePrefix {params.output_folder}_ "  
		"{config[star_all]} {config[star_indi]} > {output.log} "
		"&& mv {params.output_folder}_Aligned.sortedByCoord.out.bam {output.bam} "
		"&& rm -r {config[sample_data_dir]}/STAR_tmp_$TIME"

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
