#!/usr/bin/env nextflow

/*
###################################################################################################
################################### Module STAR_Alignment #######################################
###################################################################################################
*/

nextflow.enable.dsl = 2

/*

################################### Module STAR_Alignment #######################################

    This module contains three processes that perform an alignment with STAR from the subdirectory created with the Create_dir process (Module_OrderFiles_FastQC_rnaSeq_workflow.nf) and returns the path to the bam document generated from the alignment.

*/

if (params.help) {
	log.info ''
	exit 0
}    

process cutAdapt {

conda '/home/bio/soft/miniconda/iper-condria/envs/ipt-cutadapt'

	errorStrategy 'ignore'
	
	input:
	val directory
	
	output:
	val directory
	
	script:
	"""
	R1=\$(ls \$(echo $directory)/*_R1_*.fastq.gz)
	R2=\$(ls \$(echo $directory)/*_R2_*.fastq.gz)
	dir_to_files=\$(dirname \${R1})
	path=\$(basename \$(echo \${R1}) | awk -F "_R1" '{print \$1}')
	
	cutadapt	-a TATAAGAGACAG\
		  	-A TATAAGAGACAG\
		  	-u "${params.OutHead}"\
		  	-U "${params.OutHead}"\
		  	-m "${params.Length}"\
		  	-q "${params.Qual}"\
		  	--nextseq-trim "${params.Gtrim}"\
		  	-o "\$(echo $directory)/\${path}_R1trim_cut.fastq.gz"\
		  	-p "\$(echo $directory)/\${path}_R2trim_cut.fastq.gz"\
		  	\${R1} \${R2}\
		  	--too-short-output "\$(echo $directory)/\${path}_shortR1_cuti.fastq.gz"\
		  	--too-short-paired-output "\$(echo $directory)/\${path}_shortR2_cuti.fastq.gz" > "\$(echo $directory)/\${path}_summary.cutadapt"
	"""
}

process STAR_alignment {

conda '/home/nbarco/.conda/envs/ipt-star'

	errorStrategy 'ignore'
	
	input:
	val directory
	
	output:
	val directory
	
	script:
	"""
	R1=\$(ls \$(echo $directory)/*_R1_*.fastq.gz)
	path=\$(basename \$(echo \${R1}) | awk -F "_R1" '{print \$1}')
			  	
	STAR --genomeDir ${params.Ensembl} --runThreadN ${params.ThreadN} --readFilesIn "\$(echo $directory)/\${path}_R1trim_cut.fastq.gz" "\$(echo $directory)/\${path}_R2trim_cut.fastq.gz" --readFilesCommand "gunzip -c" --sjdbOverhang 75 --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix "\$(echo $directory)/\${path}"
							
	java -jar -Djava.io.tmpdir="\$(echo $directory)" /home/bio/soft/picard/build/libs/picard.jar AddOrReplaceReadGroups I="\$(echo $directory)/\${path}Aligned.sortedByCoord.out.bam" O="\$(echo $directory)/\${path}_aligned_merged_sorted.bam" SO=coordinate RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
	"""
}

process STAR_samtools {

conda '/home/nbarco/.conda/envs/ipt-samtools-1.20'
	
	errorStrategy 'ignore'
	
	input:
	val directory
	
	output:
	stdout
	
	script:
	"""
	R1=\$(ls \$(echo $directory)/*_R1_*.fastq.gz)
	path=\$(basename \$(echo \${R1}) | awk -F "_R1" '{print \$1}')
	
	samtools index "\$(echo $directory)/\${path}_aligned_merged_sorted.bam" > /dev/null 2>&1
	
	java -Xmx15g -jar -Djava.io.tmpdir="\$(echo $directory)" /home/bio/soft/GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R "${params.Ensembl}/Homo_sapiens.GRCh38.dna.primary_assembly.fa" -I "\$(echo $directory)/\${path}_aligned_merged_sorted.bam" -o "\$(echo $directory)/\${path}_split.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60  -U ALLOW_N_CIGAR_READS > /dev/null 2>&1
	
	samtools index "\$(echo $directory)/\${path}_split.bam" > /dev/null 2>&1
	
	echo "\$(echo $directory)/\${path}_split.bam"
	"""
}

