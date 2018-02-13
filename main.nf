#!/usr/bin/env nextflow

params.output = "./test" // Output subdirectory
params.reads = "/home/laurenrk/projects/nexus/data/*_{1,2}.fastq" // Location of forward and reverse read pairs
params.btindex = "/home/databases/cow/cow_genome" // Location of bowtie2 index of cow genome
params.output_suffix = "cow_genome" // output suffix of files created
params.num_cpu = "12" // num_cpu for running bowtie2
params.min_score = "60" // minimum alignment score -a score of 60 corresponds to approximately a perfect match over 30 nt

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


process reconcile_reads {
	
	publishDir "${params.output}/reconcile_reads_Results", mode: 'copy'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}.R1.fastq"), file("${dataset_id}.R2.fastq") from (cut_adapt_results2)
		set dataset_id, file("${dataset_id}_R12_30_f.fa.cdhit") from (cd_hit_est_results)
		
	
	output:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") into (reconcile_reads_results1, reconcile_reads_results2)	
	"""
	$baseDir/bin/reconcile_fastq_to_fasta ${dataset_id}.R1.fastq ${dataset_id}_R12_30_f.fa.cdhit  > ${dataset_id}_R1_fu.fastq \
	$baseDir/bin/reconcile_fastq_to_fasta ${dataset_id}.R2.fastq ${dataset_id}_R12_30_f.fa.cdhit  > ${dataset_id}_R2_fu.fastq
	"""
}


process RunPostFastQC {
	
	publishDir "${params.output}/FastQCResults_Post", mode: 'move'

	tag { dataset_id }

	input:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") from reconcile_reads_results1
	
	output:
		set dataset_id, file('*_fastqc.{html,zip}') into fastqc_results_post

	"""
	mkdir temp
	fastqc -f fastq $baseDir/test/reconcile_reads_Results/${dataset_id}_R1_fu.fastq $baseDir/test/reconcile_reads_Results/${dataset_id}_R2_fu.fastq -o temp
	mv temp/*.{html,zip} .
	"""
}

process align_reads {
	
	publishDir "${params.output}/aligned_reads_Results", mode: 'copy'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") from reconcile_reads_results2
		
	output:
		set dataset_id, file('${dataset_id}_R2_fuh1.fastq') into align_results
		
	"""
	bowtie2 -x $btindex --local -q -U ${dataset_id}_R1_fu.fastq --sensitive --score-min C,${params.min_score},0 --time -p $params.num_cpu -S ${dataset_id}_R1_fu.fastq.${params.output_suffix}_bt.sam 2> ${dataset_id}_R1_fu.fastq.${parmas.output_suffix}_bt.log
	$baseDir/bin/fasta_from_sam -f ${dataset_id}_R1_fu.fastq -r ${dataset_id}_R1_fu.fastq.${params.output_suffix}_bt.sam > ${dataset_id}_R1_fuh1.fastq
	$baseDir/bin/reconcile_read2_file ${dataset_id}_R1_fuh1.fastq ${dataset_id}_R2_fu.fastq f2 > ${dataset_id}_R2_fuh1.fastq
	"""
