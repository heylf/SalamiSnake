import random
import math
import os 

SRC_PATH = os.getcwd()

configfile: "config.yml"

REF_GENOME_DIR= config['ref_genome_dir']
GENOME_FASTA = REF_GENOME_DIR + config['genome_fasta']
GENOME_2BIT = REF_GENOME_DIR + config['genome_2bit']
GENOME_GTF = REF_GENOME_DIR + config['genome_gtf']
GENOME_SIZES = REF_GENOME_DIR + config['genome_sizes']

###############
## FUNCTIONS ##
###############

def list_all_values_of_dict(dictionary):
	for values in dictionary.values():
	    if isinstance(values, dict):
	      yield from list_all_values_of_dict(values)
	    else:
	      yield values

def find_all_values_on_key(key, dictionary):
	for k, v in dictionary.items():
		if k == key:
			yield v
		if isinstance(v, dict):
			yield from find_all_values_on_key(key, v)

##############################
## WRITE TOOL PARAMS TO LOG ##
##############################

file_tool_params = open(config['sample_data_dir'] + "/tool_params.txt", "w")
for i in config['all_tool_params']:
	file_tool_params.write(i + "\n")
file_tool_params.close()

######################
## SAMPLE VARIABLES ##
######################

#check if control samples exist
FIRST_READS = list(find_all_values_on_key("r1", config['clip_samples']))
ALL_SAMPLES = FIRST_READS

REPLICATES_CLIP = list(config['clip_samples'].keys())
ALL_REPLICATES = REPLICATES_CLIP

print("[NOTE] Loaded Samples.")
print("All samples: " + str(ALL_SAMPLES))
print("First read samples: " + str(FIRST_READS))

PAIR = ["r1"]
REP_NAME_CLIP = ["rep" + str(x+1) for x in range(0,len(REPLICATES_CLIP))]
SAMPLES = [REPLICATES_CLIP[0].split("_")[0]]

print(SAMPLES)

######################
## PATH DEFINITIONS ##
######################

RENAMING = config['sample_data_dir'] + "/" + config['renaming_outdir']
FASTQC_BEG_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_beg_outdir']
FASTQC_ADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_adapt_outdir']
CUTADAPT_OUTDIR = config['sample_data_dir'] + "/" + config['cutadapt_outdir']
REMOVE_TAIL_OUTDIR = config['sample_data_dir'] + "/" + config['remove_tail_outdir']
TRIMGALORE_OUTDIR = config['sample_data_dir'] + "/" + config['trim_galore_outdir']
MAPPING_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_outdir']
MAPPING_QUALITY_OUTDIR = config['sample_data_dir'] + "/" + config['mapping_quality_outdir']
PRE_FOR_UMI_OUTDIR = config['sample_data_dir'] + "/" + config['preprocessing_for_umi_outdir']
DEDUPLICAITON_OUTDIR = config['sample_data_dir'] + "/" + config['deduplication_outdir']
FASTQC_DEDUP_OUTDIR = config['sample_data_dir'] + "/" + config['fastqc_dedup_outdir']
COVERAGE_OUTDIR = config['sample_data_dir'] + "/" + config['coverage_outdir']
PEAKCALLING_OUTDIR = config['sample_data_dir'] + "/" + config['peakcalling_outdir']
POSTPROCESSING_OUTDIR = config['sample_data_dir'] + "/" + config['post_processing_outdir']
ROBUSTPEAKS_OUTDIR = config['sample_data_dir'] + "/" + config['robust_peak_search_outdir']
MOTIF_DETECTION_OUTDIR = config['sample_data_dir'] + "/" + config['motif_detection_outdir']
MOTIF_SEARCH_OUTDIR = config['sample_data_dir'] + "/" + config['motif_search_outdir']
MATEFILTER_OUTDIR = config['sample_data_dir'] + "/" + config['matefilter_outdir']
MULTIQC_OUTDIR = config['sample_data_dir'] + "/" + config['multiqc_outdir']

###################
## PREPROCESSING ##
###################

rule all:
	input: 
		expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
		expand(CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
		REF_GENOME_DIR + "/sjdbList.fromGTF.out.tab",
		expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png",
		MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_plot.png",
		expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(PEAKCALLING_OUTDIR + "/{sample}_{replicate}_crosslinkind_sites.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed",
		expand(COVERAGE_OUTDIR + "/{sample}_{replicate}_check_sorted.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_pos_strand.bedgraph", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
		expand(COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_coverage_{type}_strand.bigwig", sample=SAMPLES[0], replicate=REP_NAME_CLIP, type=["pos", "neg", "both"]),
		MOTIF_DETECTION_OUTDIR + "/meme_chip/meme-chip.html",
		MULTIQC_OUTDIR + "/multiqc_report.html"

ALL_NEW_FILE_NAMES = ["none"] * len(REPLICATES_CLIP)
i = 0 
for j in REP_NAME_CLIP:
	for k in PAIR:
		ALL_NEW_FILE_NAMES[i] = RENAMING + "/" + SAMPLES[0] + "_" + j + "_" + k + ".fastqsanger"
		i += 1

rule renaming:
	input:
		expand(config['sample_data_dir'] + "/{all}.fastqsanger", all=ALL_SAMPLES)
	output:
		file=ALL_NEW_FILE_NAMES,
		log=RENAMING + "/renaming_log.txt"
	run:
		shell("if [ ! -d {RENAMING} ]; then mkdir {RENAMING}; fi "
		"&& echo 'original_name \t new_name' > {output.log}")
		for i in range(0, len(input)): 
			shell("cp " + input[i] + " " + output.file[i] + 
			"&& echo " + input[i] + " \t " + output.file[i] + " >> {output.log}")

rule fastqc_beginning:
    input:
    	RENAMING + "/{sample}_{replicate}_{pair}.fastqsanger"
    output:
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.html",
    	FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}.fastqsanger_fastqc.zip"
    threads: 2
    conda:
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_BEG_OUTDIR} ]; then mkdir {FASTQC_BEG_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_BEG_OUTDIR}"

# remove adapter from first and second read
rule cutadapt_adapter_trimming:
	input:
		first=RENAMING + "/{sample}_{replicate}_r1.fastqsanger"
	output:
		seq_first=CUTADAPT_OUTDIR + "/{sample}_{replicate}_r1.fastqsanger",
		log=CUTADAPT_OUTDIR + "/{sample}_{replicate}.txt"
	threads: 8 
	conda:
		"envs/cutadapt.yml"
	shell:
		"if [ ! -d {CUTADAPT_OUTDIR} ]; then mkdir {CUTADAPT_OUTDIR}; fi"
		"&& cutadapt -j {threads} {config[cutadapt]} " 
		"--output={output.seq_first} {input.first} > {output.log}"

# trim of low quality bases at the ned of the parclip data 
rule remove_tail:
	input:
		first=CUTADAPT_OUTDIR + "/{samples}_{replicate}_r1.fastqsanger"
	output:
		first=REMOVE_TAIL_OUTDIR + "/{samples}_{replicate}_r1_trimmed.fastqsanger"
	threads: 2
	conda:
		"envs/bctools.yml"
	shell:
		"if [ ! -d {REMOVE_TAIL_OUTDIR} ]; then mkdir {REMOVE_TAIL_OUTDIR}; fi"
		"&& python {config[bctools]}/remove_tail.py {input.first} {config[remove_tail]} > {output.first}"

rule fastqc_after_adapter_removal:
    input:
    	REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger"
    output:
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.html",
    	FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed.fastqsanger_fastqc.zip"
    threads: 2
    conda:
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_ADAPT_OUTDIR} ]; then mkdir {FASTQC_ADAPT_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_ADAPT_OUTDIR}"

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
		"envs/star.yml"
	shell:
		"if [ -d {config[sample_data_dir]}/STAR_tmp_Index ]; then rm -r {config[sample_data_dir]}/STAR_tmp_Index; fi "
		"&& STAR --outTmpDir {config[sample_data_dir]}/STAR_tmp_Index --runThreadN {threads} --runMode genomeGenerate --genomeDir {REF_GENOME_DIR} "
		"--genomeFastaFiles {input.fasta} --sjdbGTFfile {input.annotation}"

rule star:
	input:
		first_read=REMOVE_TAIL_OUTDIR + "/{sample}_{replicate}_r1_trimmed.fastqsanger"
	output:
		log=MAPPING_OUTDIR + "/{sample}_{replicate}.txt",
		bam=MAPPING_OUTDIR + "/{sample}_{replicate}.bam",
		qual=MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out"
	threads: 4
	conda:
		"envs/star.yml"
	params:
		output_folder=MAPPING_OUTDIR + "/{sample}_{replicate}"	
	shell:
		"if [ ! -d {MAPPING_OUTDIR} ]; then mkdir {MAPPING_OUTDIR}; fi "
		"&& TIME=$(date +%N) "
		"&& if [ -d {config[sample_data_dir]}/STAR_tmp_$TIME ]; then rm -r {config[sample_data_dir]}/STAR_tmp_$TIME; fi "
		"&& STAR --runThreadN {threads} --genomeLoad NoSharedMemory --genomeDir {REF_GENOME_DIR} "   
		"--readFilesIn {input.first_read} --outTmpDir {config[sample_data_dir]}/STAR_tmp_$TIME  --outFileNamePrefix {params.output_folder}_ "  
		"{config[star]} > {output.log} "
		"&& mv {params.output_folder}_Aligned.sortedByCoord.out.bam {output.bam} "
		"&& rm -r {config[sample_data_dir]}/STAR_tmp_$TIME"

rule indexing:
	input:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam.bai"
	threads: 2
	conda:
		"envs/samtools.yml"
	shell:
		"samtools index {input}"

########################
## POST-MAP FILTERING ##
########################

# -F 0x100 do not include secondary alignments 
rule unique_reads_fitlering:
	input:
		MAPPING_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_unique_reads_fitlering.bam"
	threads: 2
	conda:
		"envs/samtools.yml"
	shell:
		"if [ ! -d {PRE_FOR_UMI_OUTDIR} ]; then mkdir {PRE_FOR_UMI_OUTDIR}; fi"
		"&& samtools view -h -F 0x100 {input} "
		"| awk '$0 ~ /^@/{{ print }} ($2 == 0 || $2 == 16)&&$0!~/XS:i/{{ print }} ' "  
		"| samtools view -bSh > {output}"

rule filter_out_unlocalized_regions_for_later_genome_versions:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_unique_reads_fitlering.bam"
	output:
		bam=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam",
		bai=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam.bai"
	threads: 2
	conda:
		"envs/samtools.yml"
	shell:
		"samtools view -h {input} "
		"""| awk -F "\t" 'BEGIN {{ OFS = FS }} {{ if ($0 ~ /^@/) {{print $0;}} else {{ if ($3 ~/^chr/) {{print $0;}} }} }}' """
		"| samtools view -bSh > {output.bam}"
		"&& samtools index {output.bam}"
 
###################
## DEDUPLICATION ##
###################

# rule deduplication:
# 	input:
# 		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
# 	output:
# 		bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}.bam",
# 		log=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_log.txt",
# 		sorted_bam=DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
# 	threads: 2
# 	conda:
# 		"envs/umi.yml"
# 	params:
# 		paired="--paired",
# 		single=""
# 	shell:
# 		"if [ ! -d {DEDUPLICAITON_OUTDIR} ]; then mkdir {DEDUPLICAITON_OUTDIR}; fi "
# 		"&& umi_tools dedup {params.paired} {config[umi_dedup]} -I {input} -S {output.bam} -L {output.log}"
# 		"&& samtools sort {output.bam} > {output.sorted_bam}"
# 		"&& samtools index {output.sorted_bam}"

#####################
## MAPPING QUALITY ##
#####################

rule finger_print_plot:
	input:
		expand(PRE_FOR_UMI_OUTDIR + "/{all}_got_umis_unlocalized_check.bam", all=ALL_REPLICATES)
	output:
		MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png"
	threads: 2
	conda:
		"envs/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& plotFingerprint --numberOfProcessors {threads} --bamfiles {input} {config[fingerprint]} "
		"--labels {ALL_REPLICATES} --plotFile {output}  --plotFileFormat 'png'"

rule correlating_bam_files_plot:
	input:
		expand(PRE_FOR_UMI_OUTDIR + "/{all}_got_umis_unlocalized_check.bam", all=ALL_REPLICATES)
	output:
		bamsummary=MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_summary.txt",
		plot=MAPPING_QUALITY_OUTDIR + "/correlating_bam_files_plot.png"
	threads: 4
	conda:
		"envs/deeptools.yml"
	shell:	
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"	
		"&& multiBamSummary bins --numberOfProcessors {threads} --outFileName {output.bamsummary} --bamfiles {input} "
		"--labels {ALL_REPLICATES} {config[multiBamSummary]} "
		"&& plotCorrelation {config[plotCorrelation]} --corData {output.bamsummary} --plotFile {output.plot}"

############################
## SECOND QUALITY CONTROL ##
############################

rule fastqc_after_dedup:
    input:
    	PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
    output:
    	FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html"
    threads: 2
    conda:
    	"envs/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_DEDUP_OUTDIR} ]; then mkdir {FASTQC_DEDUP_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_DEDUP_OUTDIR}"

#################
## PEAKCALLING ##
#################

# BE CAREFUL when using other protocols then you might need to define a different bitflag!!
# rule mate_reads_fitlering:
# 	input:
# 		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
# 	output:
# 		bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}.bam",
# 		sorted_bam=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam",
# 		sorted_bam_bai=MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam.bai"
# 	threads: 2
# 	conda:
# 		"envs/samtools.yml"
# 	shell:
# 		"if [ ! -d {MATEFILTER_OUTDIR} ]; then mkdir {MATEFILTER_OUTDIR}; fi"
# 		"&& samtools view -b -F 0x100 {input} > {output.bam}"
# 		"&& samtools sort {output.bam} > {output.sorted_bam}"
# 		"&& samtools index {output.sorted_bam}"

# rule piranha:
# 	input:
# 		MATEFILTER_OUTDIR + "/{sample}_{replicate}_sorted.bam"
# 	output:
# 		PEAKCALLING_OUTDIR + "/{sample}_{replicate}_peaks.bed"
# 	conda:
# 		"envs/piranha.yml"
# 	threads: 2 
# 	shell:
# 		"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi"
# 		"&& Piranha {config[prianha]} {input} -o {output}"

rule pureclip:
	input:
		experiment=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam",
		experiment_bai=PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam.bai",
		genome_fasta=GENOME_FASTA
	output:
		crosslinking_sites=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_crosslinkind_sites.bed",
		binding_regions=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed"
	threads: 4
	params:
		tmp=PEAKCALLING_OUTDIR + "/tmp/",
		parameters=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_parameters.txt"
	conda:
		"envs/pureclip.yml"
	shell:
		"if [ ! -d {PEAKCALLING_OUTDIR} ]; then mkdir {PEAKCALLING_OUTDIR}; fi "
		"&& pureclip -i {input.experiment} -bai {input.experiment_bai} -g {input.genome_fasta} -o {output.crosslinking_sites} -tmp {params.tmp} "
		"-or {output.binding_regions} -p {params.parameters} -nt {threads} -nta {threads} {config[pureclip]} "

###############################
## COUNTING / COVERAGE FILES ##
###############################

rule bam_to_bed:
	input:
		PRE_FOR_UMI_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check.bam"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}.bed"
	threads: 2
	conda:
		"envs/bedtools.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR} ]; then mkdir {COVERAGE_OUTDIR}; fi"
		"&& bedtools bamtobed -i {input} > {output}"

rule correct_bed_format_check:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}.bed"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check.bed"
	threads: 2
	shell:
		"""awk -v OFS="\t" -v FS="\t" '$3 != 0 {{ print $0 }}' {input} > {output} """

rule sort_beds:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check.bed"
		#COVERAGE_OUTDIR + "/{sample}_{replicate}.bed"
	output:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check_sorted.bed"
	threads: 2
	conda:
		#"envs/bedtools.yml"
		"envs/bedops.yml"
	shell:
		"&& sort-bed {input} > {output}"
		#"&& bedtools sort -i {input} > {output}"

rule calculate_coverage:
	input:
		COVERAGE_OUTDIR + "/{sample}_{replicate}_check_sorted.bed"
	output:
		pos=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_pos_strand.bedgraph",
		neg=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_neg_strand.bedgraph",
		bot=COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_both_strand.bedgraph"
	threads: 4
	conda:
		"envs/bedtools.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bedgraph ]; then mkdir {COVERAGE_OUTDIR}/bedgraph; fi"
		"&& genomeCoverageBed -i {input} -g {GENOME_SIZES} -bg -strand + > {output.pos}" 
		"&& genomeCoverageBed -i {input} -g {GENOME_SIZES} -bg -strand - > {output.neg}" 
		"&& genomeCoverageBed -i {input} -g {GENOME_SIZES} -bg > {output.bot}" 

rule calculate_coverage_bigwig:
	input:
		COVERAGE_OUTDIR + "/bedgraph/{sample}_{replicate}_coverage_{type}_strand.bedgraph"
	output:
		COVERAGE_OUTDIR + "/bigwig/{sample}_{replicate}_coverage_{type}_strand.bigwig"
	threads: 4
	conda:
		"envs/bedGraphToBigWig.yml"
	shell:
		"if [ ! -d {COVERAGE_OUTDIR}/bigwig ]; then mkdir {COVERAGE_OUTDIR}/bigwig; fi "
		"&& bedGraphToBigWig {input} {GENOME_SIZES} {output}"

##############################
## POST-PROCESSING OF PEAKS ##
##############################

rule peaks_extend_frontiers:
	input:
		bed=PEAKCALLING_OUTDIR + "/{sample}_{replicate}_binding_regions.bed",
		genome=GENOME_SIZES
	output:
		POSTPROCESSING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed"
	threads: 2
	conda:
		"envs/bedtools.yml"
	shell:
		"if [ ! -d {POSTPROCESSING_OUTDIR} ]; then mkdir {POSTPROCESSING_OUTDIR}; fi"
		"&& bedtools slop {config[peaks_extend_frontiers]} -i {input.bed} -g {input.genome} > {output}"

rule find_robust_peaks:
	input:
		expand(POSTPROCESSING_OUTDIR + "/{sample}_{replicate}_peaks_extended.bed", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
	output:
		ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed"
	threads: 2
	conda:
		"envs/bedtools.yml"
	params:
		input_folder=POSTPROCESSING_OUTDIR,
		output_folder=ROBUSTPEAKS_OUTDIR
	shell:
		"if [ ! -d {ROBUSTPEAKS_OUTDIR} ]; then mkdir {ROBUSTPEAKS_OUTDIR}; fi"
		"&& {config[find_robust_intersections]}/robust_intersections.sh {params.input_folder} {params.output_folder} bed"

#####################
## MOTIF DETECTION ##
#####################

rule extract_genomic_DNA_dreme:
	input:
		ROBUSTPEAKS_OUTDIR + "/robust_between_all.bed"
	output:
		MOTIF_DETECTION_OUTDIR + "/peaks.fa"
	threads: 2
	shell:
		"if [ ! -d {MOTIF_DETECTION_OUTDIR} ]; then mkdir {MOTIF_DETECTION_OUTDIR}; fi"
		"""&& awk '{{ print $1"\t"$2"\t"$3 }}' {input} > {MOTIF_DETECTION_OUTDIR}/peaks.bed3 """
		"&& python " + config["extract_genomic_dna"] + "/fetch_DNA_sequence.py -o {output} {MOTIF_DETECTION_OUTDIR}/peaks.bed3 {GENOME_FASTA}"

rule meme_chip:
	input: 
		MOTIF_DETECTION_OUTDIR + "/peaks.fa"
	output:
		MOTIF_DETECTION_OUTDIR + "/meme_chip/meme-chip.html"
	threads: 4
	conda:
		"envs/meme_suite.yml"
	shell:
		"meme-chip {input} -oc {MOTIF_DETECTION_OUTDIR}/meme_chip {config[meme_chip]}"

#################
## END QUALITY ##
#################

rule multiqc:
    input:
    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
    output:
    	MULTIQC_OUTDIR + "/multiqc_report.html"
    threads: 2
    conda:
    	"envs/multiqc.yml"
    shell:
    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"