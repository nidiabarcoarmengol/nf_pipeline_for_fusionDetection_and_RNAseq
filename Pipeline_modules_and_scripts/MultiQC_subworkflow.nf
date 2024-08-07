#!/usr/bin/env nextflow

/*
###################################################################################################
################################## RNA_Seq global workflow #######################################
###################################################################################################
*/

nextflow.enable.dsl=2

params.OutHead = 0
params.Length = 25
params.Qual = 10
params.Gtrim = 25
params.ThreadN = 15
params.limitSjdbInsertNsj = 2500000
params.Ensembl = "/datos3/DB/HG38/Ensembl"
params.STAR_fusion_data = "/datos3/DB/HG38/CTAT/v33"
params.Arribe_data = "/home/bio/soft/arriba_v2.1.0/database"
params.FusionCatcher_data = "/home/nbarco/test_fusion_detection_tools/data_FusionCatcher/data_for_FusionCatcher/human_v102"
params.idx_kallisto = "/home/nbarco/test_fusion_detection_tools/gencode.v38.transcripts.idx"
params.CICERO_data = "/home/nbarco/test_fusion_detection_tools/reference/"
params.bed_for_bamQC = "/home/nbarco/test_fusion_detection_tools/Homo_sapiens.GRCh38.104.chr.bed"

params.help = false

log.info """

################################## RNA_Seq global workflow #######################################

    Usage: ./RNASeq_Global_workflow.nf --directory <directory> 
    
    Options:
    --directory <directory>    Path to the directory where the R1 and R2 files are stored.
    """

include { Create_dir } from './test_Module_OrderFiles_FastQC_rnaSeq_workflow.nf'

include {cutAdapt; STAR_alignment; STAR_samtools} from './test_Module_STAR_alignment_rnaSeq_workflow.nf'

include {qualimap_module} from './test_Module_qualimap_rnaSeq_workflow.nf'

include {QC_bam_stat; QC_infer_experiment; QC_inner_distance; QC_junction_annotation; QC_junction_saturation; QC_read_distribution; QC_read_duplication; QC_tin} from './test_Module_QC_bam_rnaSeq_workflow.nf'

include {MultiQC_module} from './Module_MultiQC_report_rnaSeq_workflow.nf'


workflow STAR_subworkflow {
	
	take:
	directory
	
	main:
	cutAdapt(directory)
	STAR_alignment(cutAdapt.out)
	STAR_samtools(STAR_alignment.out)
	
	emit:
	STAR_samtools.out
}

workflow QC_bam_subworkflow_STAR {
	
	take:
	bam
	
	main:
	QC_bam_stat(bam)
	QC_infer_experiment(bam)
	QC_inner_distance(bam)
	QC_junction_annotation(bam)
	QC_junction_saturation(bam)
	QC_read_distribution(bam)
	QC_read_duplication(bam)
	QC_tin(bam)
	
	emit:
	bam_stat_out = QC_bam_stat.out
	infer_experiment_out = QC_infer_experiment.out
	inner_distance_out = QC_inner_distance.out
	junction_annotation_out = QC_junction_annotation.out
	junction_saturation_out = QC_junction_saturation.out
	read_distribution_out = QC_read_distribution.out
	read_duplication_out = QC_read_duplication.out
	tin_out = QC_tin.out
}

workflow {

	main:
	if (params.directory) {
		path_r1 = "${params.directory}/*_R1_*.fastq.gz" 
		R1_files = channel.fromPath(path_r1)
		Create_dir(R1_files)
		STAR_subworkflow(Create_dir.out)
		qualimap_module(STAR_subworkflow.out)
		qc_bam_out = QC_bam_subworkflow_STAR(STAR_subworkflow.out)
		MultiQC_module(qc_bam_out.bam_stat_out, Create_dir.out, qc_bam_out.infer_experiment_out, qc_bam_out.inner_distance_out, qc_bam_out.junction_annotation_out, qc_bam_out.junction_saturation_out, qualimap_module.out, qc_bam_out.read_distribution_out, qc_bam_out.read_duplication_out, qc_bam_out.tin_out)
	}
}
