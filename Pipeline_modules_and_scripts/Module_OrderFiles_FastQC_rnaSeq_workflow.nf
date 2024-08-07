#!/usr/bin/env nextflow

/*
##################################################################################################
############################### Module order files in directories ################################
##################################################################################################
*/

nextflow.enable.dsl=2

/*

############################### Module order files in directories ################################

    This module sorts the R1 and R2 documents of the same sample in a subfolder and performs the fastq quality control with the FastQC tool, finally it returns the directory of the created subfolder.

*/

if (params.help) {
	log.info ''
	exit 0
}

process Create_dir {

conda '/home/bio/soft/miniconda/iper-condria/envs/ipt-fastqc'

	errorStrategy 'ignore'
	
	input:
	val R1_file
	
	output:
	stdout
	
	script:
	"""
	path=\$(basename \$(echo $R1_file) | awk -F "_R1" '{print \$1}')
	mkdir -p "${params.directory}/\${path}"
	mkdir -p "${params.directory}/\${path}/FastQC_results"
	
	fastqc --outdir "${params.directory}/\${path}/FastQC_results" -f fastq $R1_file > /dev/null 2>&1
	R2_file=\$(echo $R1_file | sed 's/R1/R2/')
	fastqc --outdir "${params.directory}/\${path}/FastQC_results" -f fastq \${R2_file} > /dev/null 2>&1
	
	mv "$R1_file" "${params.directory}/\${path}"
	mv "\${R2_file}" "${params.directory}/\${path}"
	
	echo "${params.directory}/\${path}"
	"""
}
