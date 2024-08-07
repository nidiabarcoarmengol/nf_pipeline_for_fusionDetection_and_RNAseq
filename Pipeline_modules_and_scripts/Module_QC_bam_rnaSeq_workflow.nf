#!/usr/bin/env nextflow

/*
###################################################################################################
###################################### Module QC_bam_file #######################################
###################################################################################################
*/

nextflow.enable.dsl = 2

/*

###################################### Module QC_bam_file #######################################

    This module contains different processes to perform quality control of STAR and Kallisto BAM files with the RSeQC and Preseq tools.

*/

if (params.help) {
	log.info ''
	exit 0
}    

process QC_bam_stat {

conda ("/home/nbarco/.conda/envs/ipa-rseqc")

	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout

	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/bam_stat_results
	
	bam_stat.py -i \${bam_al} > \${file_dir}/bam_stat_results/\${path}.bam_stat.txt
	
	echo \${path}
	"""
}

process QC_infer_experiment {

conda ("/home/nbarco/.conda/envs/ipa-rseqc")
	
	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout
	
	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
		echo "Entra"	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		echo \${path}
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		echo "NO Entra"
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		echo \${path}
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/infer_experiment_results
	
	infer_experiment.py -i \${bam_al} -r ${params.bed_for_bamQC} > \${file_dir}/infer_experiment_results/\${path}.infer_experiment.txt
	"""
}

process QC_inner_distance {
	
conda ("/home/nbarco/.conda/envs/ipa-rseqc")

	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout

	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
		echo "Entra"	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		echo \${path}
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		echo "NO Entra"
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		echo \${path}
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/inner_distance_results
	
	inner_distance.py -i \${bam_al} -o \${file_dir}/inner_distance_results/\${path}.inner_distance -r ${params.bed_for_bamQC} > \${file_dir}/inner_distance_results/\${path}.inner_distance_mean.txt
	"""
}

process QC_junction_annotation {

conda ("/home/nbarco/.conda/envs/ipa-rseqc")

	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout

	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
		echo "Entra"	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		echo \${path}
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		echo "NO Entra"
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		echo \${path}
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/junction_annotation_results
	
	junction_annotation.py -i \${bam_al} -r ${params.bed_for_bamQC} -o \${file_dir}/junction_annotation_results/\${path}.junction_annotation > \${file_dir}/junction_annotation_results/\${path}.junction_annotation.log
	"""
}

process QC_junction_saturation {

conda ("/home/nbarco/.conda/envs/ipa-rseqc")

	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout

	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
		echo "Entra"	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		echo \${path}
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		echo "NO Entra"
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		echo \${path}
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/junction_saturation_results
	
	junction_saturation.py -i \${bam_al} -r ${params.bed_for_bamQC} -o \${file_dir}/junction_saturation_results/\${path}.junction_saturation
	"""
}

process QC_read_distribution {

conda ("/home/nbarco/.conda/envs/ipa-rseqc")

	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout

	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
		echo "Entra"	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		echo \${path}
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		echo "NO Entra"
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		echo \${path}
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/read_distribution_results
	
	read_distribution.py -i \${bam_al} -r ${params.bed_for_bamQC} > \${file_dir}/read_distribution_results/\${path}.read_distribution.txt
	"""
}

process QC_read_duplication {

conda ("/home/nbarco/.conda/envs/ipa-rseqc")

	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout

	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
		echo "Entra"	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		echo \${path}
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		echo "NO Entra"
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		echo \${path}
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/read_duplication_results
	
	read_duplication.py -i \${bam_al} -o \${file_dir}/read_duplication_results/\${path}.read_duplication
	"""
}

process QC_tin {

conda ("/home/nbarco/.conda/envs/ipa-rseqc")

	errorStrategy 'ignore'
	
	input:
	val bam
	
	output:
	stdout

	script:
	"""
	file_dir=\$(dirname \$(echo $bam))
	
	if [[ \$(basename \$(echo $bam)) != "pseudoalignments.bam" ]]; then
		echo "Entra"	
		path=\$(basename \$(echo $bam) | awk -F "_split" '{print \$1}')
		echo \${path}
		bam_al=\${file_dir}/\${path}_aligned_merged_sorted.bam
	else
		echo "NO Entra"
		dir=\$(dirname \$(dirname \$(echo $bam)))
		R1=\$(ls \${dir}/*_R1_*.fastq.gz)
		path=\$(basename \${R1} | awk -F "_R1" '{print \$1}')
		echo \${path}
		bam_al=\$(echo $bam)
	fi
	
	mkdir -p \${file_dir}/tin_results
	
	tin.py -i \${bam_al} -r ${params.bed_for_bamQC} > \${file_dir}/tin_results/\${path}.tin.txt
	"""
}

