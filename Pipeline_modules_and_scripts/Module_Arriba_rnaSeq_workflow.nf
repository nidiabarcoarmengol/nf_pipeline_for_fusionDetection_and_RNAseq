#!/usr/bin/env nextflow

/*
##################################################################################################
###################################### Module Arriba #############################################
##################################################################################################
*/

nextflow.enable.dsl=2

/*

#################################### Module STAR_Arriba ##########################################

    This module performs the merger detection with the Arriba tool starting as input from the subdirectory returned by the CreateDir process (Module_OrderFiles_FastQC_rnaSeq_workflow.nf) and returns the subdirectory containing the results of the process.		

*/

if (params.help) {
	log.info ''
	exit 0
}    

process Star_Arriba_module {
	
conda '/home/nbarco/.conda/envs/ipx-arriba'

	errorStrategy 'ignore'
	
	input:
	val directory
	
	output:
	stdout
	
	label 'Heavy_Tools'
		
	script:
	"""
	R1=\$(ls \$(echo $directory)/*_R1_*.fastq.gz)
	R2=\$(ls \$(echo $directory)/*_R2_*.fastq.gz)
	path=\$(basename \$(echo \${R1}) | awk -F "_R1" '{print \$1}')
	
	mkdir -p \$(echo $directory)/Arriba_results
	
	STAR	--runThreadN 8\
		--genomeDir ${params.Ensembl}\
		--genomeLoad NoSharedMemory\
		--readFilesIn \${R1} \${R2}\
		--readFilesCommand zcat\
		--outStd BAM_Unsorted\
		--outSAMtype BAM Unsorted\
		--outSAMunmapped Within\
		--outBAMcompression 0\
		--outFilterMultimapNmax 50\
		--peOverlapNbasesMin 10\
		--alignSplicedMateMapLminOverLmate 0.5\
		--alignSJstitchMismatchNmax 5 -1 5 5\
		--chimSegmentMin 10\
		--chimOutType WithinBAM HardClip\
		--chimJunctionOverhangMin 10\
		--chimScoreDropMax 30\
		--chimScoreJunctionNonGTAG 0\
		--chimScoreSeparation 1\
		--chimSegmentReadGapMax 3\
		-chimMultimapNmax 50|
	arriba	-x /dev/stdin\
		-o "\$(echo $directory)/Arriba_results/\${path}_fusions.tsv"\
		-O "\$(echo $directory)/Arriba_results/\${path}_fusions_discarded.tsv"\
		-a "${params.Ensembl}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"\
		-g "${params.Ensembl}/Homo_sapiens.GRCh38.104.chr.gtf"\
		-b "${params.Arriba_data}/blacklist_hg38_GRCh38_v2.1.0.tsv.gz"\
		-k "${params.Arriba_data}/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz"\
		-t "${params.Arriba_data}/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz"\
		-p "${params.Arriba_data}/protein_domains_hg38_GRCh38_v2.1.0.gff3" > "\$(echo $directory)/Arriba_results/Arriba_stdout.txt"
		
	echo \${path}
	"""
}

