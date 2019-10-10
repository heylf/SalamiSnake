import random
import math
import os

#####################
## MAPPING QUALITY ##
#####################

rule genome_bed_conversion:
	input:
		GENOME_GTF
	output:
		MAPPING_QUALITY_OUTDIR + "/reference_genome.bed"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"""&& awk -v FS="\t" -v OFS="\t" '$0 !~/^#/ {{print $1,$4,$5,$3,".",$7}}' {input} > {output}"""

rule orientation:
	input:
		bam=MAPPING_OUTDIR + "/{sample}_{replicate}.bam",
		genome_bed=MAPPING_QUALITY_OUTDIR + "/reference_genome.bed"
	output:
		MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_orientation.txt"
	threads: 2
	conda:
		config["conda_envs"] + "/orientation.yml"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& echo {config[orientation]} >> {file_tool_params}"
		"&& python2 {config[infer_experiment]}/infer_experiment.py -i {input.bam} -r {input.genome_bed} --sample-size 200000 --mapq 30 > {output}"

rule compute_gc_bias_plots:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		plot=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_gc_bias_plot.png",
		file=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_gc_bias_data.txt"
	threads: 2
	conda:
		config["conda_envs"] + "/deeptools.yml"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& echo {config[gc_bias]} >> {file_tool_params}"
		"&& computeGCBias --numberOfProcessors {threads} --bamfile {input} --GCbiasFrequenciesFile {output.file} {config[gc_bias]} "
		"--genome {GENOME_2BIT} --biasPlot {output.plot} --plotFileFormat png"

# for picard please conda install R
# Swap in /scratch/bi01/heylf/miniconda3/envs/picard/etc/conda/activate.d/java_home.sh the first line
# with the second line.
rule estimate_insert_size:
	input:
		DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}.bam"
	output:
		plot=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_insert_size_plot.pdf",
		report=MAPPING_QUALITY_OUTDIR + "/{sample}_{replicate}_insert_size_report.txt"
	threads: 2
	conda:
		config["conda_envs"] + "/picard.yml"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& echo {config[estimate_insert_size]} >> {file_tool_params}"
		"&& source activate picard"
		"&& picard CollectInsertSizeMetrics INPUT={input} OUTPUT={output.report} "
		"HISTOGRAM_FILE={output.plot} {config[estimate_insert_size]} "
		"&& source deactivate"

rule finger_print_plot:
	input:
		expand(DEDUPLICAITON_OUTDIR + "/{sample}_sorted.bam", sample=ALL_REPLICATES)
	output:
		MAPPING_QUALITY_OUTDIR + "/fingerprint_plot.png"
	threads: 2
	conda:
		config["conda_envs"] + "/deeptools.yml"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& echo {config[fingerprint]} >> {file_tool_params}"
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
		config["conda_envs"] + "/deeptools.yml"
	shell:
		"if [ ! -d {MAPPING_QUALITY_OUTDIR} ]; then mkdir {MAPPING_QUALITY_OUTDIR}; fi"
		"&& echo {config[multiBamSummary]} >> {file_tool_params}"
		"&& multiBamSummary bins --numberOfProcessors {threads} --outFileName {output.bamsummary} --bamfiles {input} "
		"--labels {ALL_REPLICATES} {config[multiBamSummary]} "
		"&& plotCorrelation {config[plotCorrelation]} --corData {output.bamsummary} --plotFile {output.plot}"

rule fastqc_after_dedup:
    input:
    	DEDUPLICAITON_OUTDIR + "/{sample}_{replicate}_sorted.bam"
    output:
    	FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_sorted_fastqc.html"
    threads: 2
    conda:
    	config["conda_envs"] + "/fastqc.yml"
    shell:
    	"if [ ! -d {FASTQC_DEDUP_OUTDIR} ]; then mkdir {FASTQC_DEDUP_OUTDIR}; fi"
    	"&& fastqc {input} --outdir {FASTQC_DEDUP_OUTDIR}"

#################
## END QUALITY ##
#################


if ( PROTOCOL == "ChIP" or PROTOCOL == "PARCLIP" ):

	mapping_suffix = ""
	if ( mapper == "STAR" ):
		mapping_suffix = "_Log.final.out"
	if ( mapper == "Bowtie2" ):
		mapping_suffix = "_log.txt"

	if ( control == "yes" ):
		if ( demultiplexed == "yes" ):
			rule multiqc:
			    input:
			    	begin=[expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    		   expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, pair=PAIR)],
			    	trim=[expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    		  expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, pair=PAIR)],
			    	dedup=[expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    		   expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
			    	mapping=[expand(MAPPING_OUTDIR + "/{sample}_{replicate}" + mapping_suffix, sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    			 expand(MAPPING_OUTDIR + "/{sample}_{replicate}" + mapping_suffix, sample=SAMPLES[1], replicate=REP_NAME_CONTROL)]
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
		else:
			rule multiqc:
			    input:
			    	begin=expand(FASTQC_BEG_OUTDIR + "/{s}_{r}_{pair}_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{s}_{r}_{pair}_trimmed_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	dedup=[expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    		   expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
			    	mapping=[expand(MAPPING_OUTDIR + "/{sample}_{replicate}" + mapping_suffix, sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    			 expand(MAPPING_OUTDIR + "/{sample}_{replicate}" + mapping_suffix, sample=SAMPLES[1], replicate=REP_NAME_CONTROL)]
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
	else:
		if ( demultiplexed == "yes" ):
			rule multiqc:
			    input:
			    	begin=expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}" + mapping_suffix, sample=SAMPLES[0], replicate=REP_NAME_CLIP)
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
		else:
			rule multiqc:
			    input:
			    	begin=expand(FASTQC_BEG_OUTDIR + "/{s}_{r}_{pair}_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{s}_{r}_{pair}_trimmed_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}" + mapping_suffix, sample=SAMPLES[0], replicate=REP_NAME_CLIP)
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"

else:

	if ( control == "yes" ):
		if ( demultiplexed == "yes" ):
			rule multiqc:
			    input:
			    	begin=[expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    		   expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, pair=PAIR)],
			    	trim=[expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    		  expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL, pair=PAIR)],
			    	be_dedup=[expand(FASTQC_BEFORE_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    		      expand(FASTQC_BEFORE_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
			    	dedup=[expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_sorted_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    		   expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_sorted_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
			    	mapping=[expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    			 expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)]
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {FASTQC_BEFORE_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
		else:
			rule multiqc:
			    input:
			    	begin=expand(FASTQC_BEG_OUTDIR + "/{s}_{r}_{pair}_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{s}_{r}_{pair}_trimmed_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	be_dedup=[expand(FASTQC_BEFORE_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    		      expand(FASTQC_BEFORE_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
			    	dedup=[expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_sorted_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    		   expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_sorted_fastqc.html", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)],
			    	mapping=[expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    			 expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[1], replicate=REP_NAME_CONTROL)]
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {FASTQC_BEFORE_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
	else:
		if ( demultiplexed == "yes" ):
			rule multiqc:
			    input:
			    	begin=expand(FASTQC_BEG_OUTDIR + "/{sample}_{replicate}_{pair}_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{sample}_{replicate}_{pair}_trimmed_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP, pair=PAIR),
			    	be_dedup=expand(FASTQC_BEFORE_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {FASTQC_BEFORE_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
		else:
			rule multiqc:
			    input:
			    	begin=expand(FASTQC_BEG_OUTDIR + "/{s}_{r}_{pair}_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	trim=expand(FASTQC_ADAPT_OUTDIR + "/{s}_{r}_{pair}_trimmed_fastqc.html", s=MULTIPLEX_SAMPLE_NAME, r="rep1", pair=PAIR),
			    	be_dedup=expand(FASTQC_BEFORE_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    	dedup=expand(FASTQC_DEDUP_OUTDIR + "/{sample}_{replicate}_got_umis_unlocalized_check_fastqc.html", sample=SAMPLES[0], replicate=REP_NAME_CLIP),
			    	mapping=expand(MAPPING_OUTDIR + "/{sample}_{replicate}_Log.final.out", sample=SAMPLES[0], replicate=REP_NAME_CLIP)
			    output:
			    	MULTIQC_OUTDIR + "/multiqc_report.html"
			    threads: 2
			    conda:
			    	config["conda_envs"] + "/multiqc.yml"
			    shell:
			    	"if [ ! -d {MULTIQC_OUTDIR} ]; then mkdir {MULTIQC_OUTDIR}; fi"
			    	"&& multiqc -s {FASTQC_BEG_OUTDIR} {FASTQC_ADAPT_OUTDIR} {FASTQC_DEDUP_OUTDIR} {FASTQC_BEFORE_DEDUP_OUTDIR} {input.mapping} --outdir {MULTIQC_OUTDIR}"
