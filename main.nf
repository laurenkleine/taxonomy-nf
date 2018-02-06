#!/usr/bin/env nextflow

params.output = "./test" 
params.reads = "/home/laurenrk/projects/nexus/data/*_{1,2}.fastq"

Channel
	.fromFilePairs( params.reads, flat: true) /*When true the matching files are produced as sole elements in the emitted tuples*/
	.ifEmpty { return "Error" }
	.into { read_pairs; qc_pairs }

process RunPreFastQC {
	publishDir "${params.output}/FastQCResults_Pre", mode: 'move'

	tag { dataset_id }

	input:
		set dataset_id, file(forward), file(reverse) from (read_pairs)
	
	output:
		set dataset_id, file('*_fastqc.{html,zip}') into (fastqc_results_pre)

	"""
	mkdir temp
	fastqc -f fastq $forward $reverse -o temp
	mv temp/*.{html,zip} .
	"""
}

process RunCutAdapt {
	publishDir "${params.output}/RunCutAdapt", mode: 'copy'

	tag { dataset_id }

	input:
		set dataset_id, file(forward), file(reverse) from (qc_pairs)

	output:
		set dataset_id, file("${dataset_id}.R1.fastq"), file("${dataset_id}.R2.fastq") into (cut_adapt_results1, cut_adapt_results2)

	"""
	cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -A AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT -G CTGTCTCTTATACACATCT -q 30,30 --minimum-length 80 -u 1 -o ${dataset_id}.R1.fastq -p ${dataset_id}.R2.fastq $forward $reverse
	"""	
}

process cd_hit_est {
	publishDir "${params.output}/cd_hit_est_Results", mode: 'copy'

	tag { dataset_id }

	input:
		set dataset_id, file("${dataset_id}.R1.fastq"), file("${dataset_id}.R2.fastq") from (cut_adapt_results1)
	
	output:
		set dataset_id, file("${dataset_id}_R12_30_f.fa.cdhit") into (cd_hit_est_results)

	"""
	$baseDir/bin/fastq_to_fasta < ${dataset_id}.R1.fastq > ${dataset_id}_R1_f.fa
	$baseDir/bin/fastq_to_fasta < ${dataset_id}.R2.fastq > ${dataset_id}_R2_f.fa
	$baseDir/bin/concat_fasta_records -5 30 ${dataset_id}_R1_f.fa ${dataset_id}_R2_f.fa > ${dataset_id}_R12_30_f.fa
	cd-hit-est -c 0.96 -i ${dataset_id}_R12_30_f.fa  -o ${dataset_id}_R12_30_f.fa.cdhit  -T 12 -M 0 
	"""
}


/* fromPath */


process reconcile_reads {
	publishDir "${params.output}/reconcile_reads_Results", mode: 'copy'

	tag { dataset_id }

	input:
		set dataset_id, file(forward), file(reverse) from (cut_adapt_results2)
		set dataset_id, file("${dataset_id}_R12_30_f.fa.cdhit") from (cd_hit_est_results)
		
	
	output:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") into (paired_filtered_fastq)

	"""
	$baseDir/bin/reconcile_fastq_to_fasta ${dataset_id}.R1.fastq $baseDir/test/cd_hit_est_Results/${dataset_id}_R12_30_f.fa.cdhit  > ${dataset_id}_R1_fu.fastq
	$baseDir/bin/reconcile_fastq_to_fasta ${dataset_id}.R2.fastq $baseDir/test/cd_hit_est_Results/${dataset_id}_R12_30_f.fa.cdhit  > ${dataset_id}_R2_fu.fastq
	"""
}
