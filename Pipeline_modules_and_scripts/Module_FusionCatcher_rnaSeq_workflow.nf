#!/usr/bin/env nextflow

/*
###################################################################################################
#################################### Module Fusion_catcher #######################################
###################################################################################################
*/

nextflow.enable.dsl=2

/*

#################################### Module Fusion_catcher #######################################

    This module performs the fusion detection with the FusionCatcher tool starting as input from the subdirectory returned by the CreateDir process (Module_OrderFiles_FastQC_rnaSeq_workflow.nf) and returns the subdirectory containing the results of the process.		

*/

if (params.help) {
	log.info ''
	exit 0
}    

process Fusion_catcher_module {

conda '/home/nbarco/.conda/envs/ipx-fusioncatcher2'
	
	errorStrategy 'ignore'
	
	input:
	val directory
	
	output:
	stdout
	
	script:
	"""
	R1=\$(ls \$(echo $directory)/*_R1_*.fastq.gz)
	R2=\$(ls \$(echo $directory)/*_R2_*.fastq.gz)
	
	mkdir -p \$(echo $directory)/FusionCatcher_results
	
	fusioncatcher	-d ${params.FusionCatcher_data} \
			-i \$(echo "\${R1},\${R2}") \
			-o \$(echo $directory)/FusionCatcher_results \
			--limitSjdbInsertNsj ${params.limitSjdbInsertNsj}
			
	echo \$(echo $directory)/FusionCatcher_results
	"""
}

