#!/usr/bin/env nextflow

/*
#################################################################################################
##################################### Global workflow ###########################################
#################################################################################################
*/

nextflow.enable.dsl=2

params.mix_fusions = "/path/to/directory/when/we/found/mix.fusions_scpript"
params.Sample_type = "NA"
params.Disease_name = "NA"
params.OutHead = 0
params.Length = 25
params.Qual = 10
params.Gtrim = 25
params.ThreadN = 15
params.limitSjdbInsertNsj = 2500000
params.Ensembl = "/path/to/directory/with/ensembl/reference"
params.STAR_fusion_data = "/path/to/directory/with/reference/data/for/STAR_fusion"
params.Arriba_data = "/path/to/directory/with/reference/data/for/Arriba"
params.FusionCatcher_data = "/path/to/directory/with/reference/data/for/FusionCatcher"
params.idx_kallisto = "/path/to/kallisto.idx"
params.CICERO_data = "/path/to/directory/with/reference/data/for/CICERO"
params.bed_for_bamQC = "/path/to/reference/file/for/QC/reference.bed"

params.help = false

log.info """

##################################### Global workflow ###########################################

    Usage: ./RNASeq_Global_workflow.nf --directory <directory> 
    
    Options:
    --directory <directory>    Path to the directory where the R1 and R2 files are stored.
    """

include { Create_dir } from './Module_OrderFiles_FastQC_rnaSeq_workflow.nf'

include {cutAdapt; STAR_alignment; STAR_samtools} from './Module_STAR_alignment_rnaSeq_workflow.nf'

include {STAR_Fusion_module} from './Module_STAR_Fusion_rnaSeq_workflow.nf'

include {Star_Arriba_module} from './Module_Arriba_rnaSeq_workflow.nf'

include {Fusion_catcher_module} from './Module_FusionCatcher_rnaSeq_workflow.nf'

include {Kallisto_alignment_module} from './Module_Kallisto_rnaSeq_workflow.nf'

include {CICERO_Fusion_module} from './Module_CICERO_rnaSeq_workflow.nf'

include {qualimap_module} from './Module_qualimap_rnaSeq_workflow.nf'

include {QC_bam_stat; QC_infer_experiment; QC_inner_distance; QC_junction_annotation; QC_junction_saturation; QC_read_distribution; QC_read_duplication; QC_tin} from './Module_QC_bam_rnaSeq_workflow.nf'

include {Htseq_module; FeatureCounts_module} from './Module_Counts_matrix_rnaSeq_workflow.nf'

include {MultiQC_module} from './Module_MultiQC_report_rnaSeq_workflow.nf'

include {Fusions} from './Module_Fusions_resume.nf'


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
		STAR_Fusion_module(Create_dir.out)
		Star_Arriba_module(Create_dir.out)
		Fusion_catcher_module(Create_dir.out)
		Kallisto_alignment_module(Create_dir.out)
		STAR_subworkflow(Create_dir.out)
		Htseq_module(STAR_subworkflow.out.collect())
		FeatureCounts_module(STAR_subworkflow.out.collect())
		CICERO_Fusion_module(STAR_subworkflow.out)
		qualimap_module(STAR_subworkflow.out)
		qc_bam_out = QC_bam_subworkflow_STAR(STAR_subworkflow.out)
		MultiQC_module(qc_bam_out.bam_stat_out, Create_dir.out, qc_bam_out.infer_experiment_out, qc_bam_out.inner_distance_out, qc_bam_out.junction_annotation_out, qc_bam_out.junction_saturation_out, qualimap_module.out, qc_bam_out.read_distribution_out, qc_bam_out.read_duplication_out, qc_bam_out.tin_out)
		Fusions(Create_dir.out, Star_Arriba_module.out, CICERO_Fusion_module.out, STAR_Fusion_module.out, Fusion_catcher_module.out)
	}
}

