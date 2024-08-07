#!/usr/bin/env nextflow

/*
###################################################################################################
##################################### Module counts_matrix #######################################
###################################################################################################
*/

nextflow.enable.dsl = 2

/*

##################################### Module counts_matrix #######################################

    This module makes the count matrix of all bam files created by STAR with two different tools, htseq and featurecounts.

*/

if (params.help) {
	log.info ''
	exit 0
}    

process Htseq_module {

conda '/home/nbarco/.conda/envs/ipx-htseq'
	
	errorStrategy 'ignore'
	
	input:
	val bam
	
	script:
	"""
	mkdir -p ${params.directory}/Htseq_result
	htseq-count --mode intersection-strict ${params.directory}/*/*_split.bam ${params.Ensembl}/Homo_sapiens.GRCh38.104.chr.gtf > ${params.directory}/Htseq_result/Htseq_counts.txt
	"""
}

process FeatureCounts_module {

conda '/home/nbarco/.conda/envs/ipx-subread'

	errorStrategy 'ignore'
	
	input:
	val bam	
	
	script:
	"""
	mkdir -p ${params.directory}/FeatureCounts_result
	featureCounts -p -a ${params.Ensembl}/Homo_sapiens.GRCh38.104.chr.gtf -o ${params.directory}/FeatureCounts_result/featureCounts_counts.txt ${params.directory}/*/*_split.bam
	"""
}
