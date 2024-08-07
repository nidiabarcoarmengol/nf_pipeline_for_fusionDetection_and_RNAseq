#!/usr/bin/env nextflow

/*
###################################################################################################
##################################### Module multiQC_report #######################################
###################################################################################################
*/

nextflow.enable.dsl = 2

/*

##################################### Module multiQC_report #######################################

    This module makes quality control report with multiQC tool.

*/

if (params.help) {
	log.info ''
	exit 0
}    

process MultiQC_module {

conda '/home/nbarco/.conda/envs/ipa-multiqc-v1.20'
	
	errorStrategy 'ignore'
	
	input:
	val path_bam_stat
	val out_fastqc
	val out_inf_exp
	val out_inn_dist
	val out_junc_ann
	val out_junc_sat
	val out_qualimap
	val out_read_dis
	val out_read_dup
	val out_tin
	
	script:
	"""
	multiqc -f -n \$(echo $path_bam_stat)_multiqc_report \
	${params.directory}/\$(echo $path_bam_stat)/bam_stat_results/ \
	${params.directory}/\$(echo $path_bam_stat)/FastQC_results/ \
	${params.directory}/\$(echo $path_bam_stat)/infer_experiment_results/ \
	${params.directory}/\$(echo $path_bam_stat)/inner_distance_results/ \
	${params.directory}/\$(echo $path_bam_stat)/junction_annotation_results/ \
	${params.directory}/\$(echo $path_bam_stat)/junction_saturation_results/ \
	${params.directory}/\$(echo $path_bam_stat)/qualimap_result/ \
	${params.directory}/\$(echo $path_bam_stat)/read_distribution_results/ \
	${params.directory}/\$(echo $path_bam_stat)/read_duplication_results/ \
	${params.directory}/\$(echo $path_bam_stat)/tin_results/ \
	-o ${params.directory}/\$(echo $path_bam_stat)
	"""
}
