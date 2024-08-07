#!/usr/bin/env nextflow

/*
###################################################################################################
##################################### Module Fusions_result ######################################
###################################################################################################
*/

nextflow.enable.dsl = 2

/*

##################################### Module Fusions_result ######################################

    This module makes the resume of fusions.

*/

if (params.help) {
	log.info ''
	exit 0
}    

process Fusions {
	
	errorStrategy 'ignore'
	
	input:
	val directory
	val Arriba
	val CICERO
	val STAR
	val FusionCatcher
	
	output:
	stdout
	
	script:
	"""
	R1=\$(ls \$(echo $directory)/*_R1_*.fastq.gz)
	path=\$(basename \$(echo \${R1}) | awk -F "_R1" '{print \$1}')
	
	mkdir -p \$(echo $directory)/Fusions_Results
	
	${params.mix_fusions}/mix_fusions_3.sh -a "\$(echo $directory)/Arriba_results/\${path}_fusions.tsv" \
	-f "\$(echo $directory)/FusionCatcher_results/final-list_candidate-fusion-genes.txt" \
	-u "\$(echo $directory)/STAR_Fusion_results/star-fusion.fusion_predictions.tsv" \
	-w "\$(echo $directory)/CICERO_results/CICERO_DATADIR/\${path}_aligned_merged_sorted_chr/final_fusions.txt" \
	-b "\${path}" -d ${params.Sample_type} -x ${params.Disease_name} \
	-o "\$(echo $directory)/Fusions_Results"
	"""
}
