#!/usr/bin/env nextflow

/*
###################################################################################################
####################################### Module CICERO ############################################
###################################################################################################
*/

nextflow.enable.dsl = 2

/*

####################################### Module CICERO ############################################

    This module performs the fusion detection with the CICERO tool starting as input from the BAM file created by the STAR aligner and returns the subdirectory containing the results of the process.

*/

if (params.help) {
	log.info ''
	exit 0
}    

process CICERO_Fusion_module {

conda '/home/nbarco/.conda/envs/ipt-samtools-1.9'
	
	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout
	
	script:
	"""
	files_dir=\$(dirname \$(echo $bam))
	path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
	
	mkdir -p \${files_dir}/CICERO_results
	
	samtools view -H \${files_dir}/\${path}_aligned_merged_sorted.bam | sed -e 's/SN:\\([0-9XY]\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - \${files_dir}/\${path}_aligned_merged_sorted.bam > \${files_dir}/CICERO_results/\${path}_aligned_merged_sorted_chr.bam
	
	samtools index \${files_dir}/CICERO_results/\${path}_aligned_merged_sorted_chr.bam
	
	docker run -u root -v ${params.CICERO_data}:/reference -v \${files_dir}/CICERO_results:/data magearna:5000/iper/cicero Cicero.sh -b /data/\${path}_aligned_merged_sorted_chr.bam -g GRCh38_no_alt -r /reference -o /data -n ${params.ThreadN}
	
	echo \${path}
	
	"""
}

