#!/usr/bin/env nextflow

params.output = "./test" // Output subdirectory
params.reads = "/home/laurenrk/projects/nexus/data/*_R{1,2}.fastq" // Location of forward and reverse read pairs /*  sed 's/\/1*$/ 1/g' SRR532663_1.fastq */
params.btindex = "/home/databases/croc/croc_genome" // Location of bowtie2 index of cow genome
params.output_suffix = "croc_genome" // output suffix of files created
params.num_cpu = "12" // num_cpu for running bowtie2
params.min_score = "60" // minimum alignment score -a score of 60 corresponds to approximately a perfect match over 30 nt

Channel
	.fromFilePairs( params.reads, flat: true) /*When true the matching files are produced as sole elements in the emitted tuples*/
	.ifEmpty { return "Error" }
	.into { read_pairs; qc_pairs }


/*dumby process*/

process RunPreFastQC {
	
	publishDir "${params.output}/FastQCResults_Pre", mode: 'move'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file(forward), file(reverse) from (qc_pairs)
	
	output:
		set dataset_id, file('*_fastqc.{html,zip}') into (fastqc_results_pre)
	"""
	mkdir temp
	fastqc -f fastq $forward $reverse -o temp
	mv temp/*.{html,zip} .
	"""
}


process RunCutAdapt {
	
	publishDir "${params.output}/RunCutAdapt", mode: 'link'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file(forward), file(reverse) from (read_pairs)
	output:
		set dataset_id, file("${dataset_id}.R1.fastq"), file("${dataset_id}.R2.fastq") into (cut_adapt_results1, cut_adapt_results2)
	"""
	cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -A AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT -G CTGTCTCTTATACACATCT -q 30,30 --minimum-length 80 -u 1 -o ${dataset_id}.R1.fastq -p ${dataset_id}.R2.fastq $forward $reverse
	"""	
}

process cd_hit_est {
	
	publishDir "${params.output}/cd_hit_est_Results", mode: 'link'
	
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
	publishDir "${params.output}/reconcile_reads_Results", mode: 'link'

	tag { dataset_id }

	input:
		set dataset_id, file("${dataset_id}.R1.fastq"), file("${dataset_id}.R2.fastq") from cut_adapt_results2
		set dataset_id, file("${dataset_id}_R12_30_f.fa.cdhit") from (cd_hit_est_results)
	
	output:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") into (reconcile_reads_results1,reconcile_reads_results2)

	"""
	$baseDir/bin/reconcile_fastq_to_fasta ${dataset_id}.R1.fastq ${dataset_id}_R12_30_f.fa.cdhit  > ${dataset_id}_R1_fu.fastq
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

/*
TODO:
*write align_reads process
*write check input process
*maybe write check_single_or_paired process
*/ 



process align_reads1 {
	
	publishDir "${params.output}/aligned_reads_Results", mode: 'link'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") from reconcile_reads_results2
		
	output:
		set dataset_id, file("${dataset_id}_R1_fu_bt.sam"), file("${dataset_id}_R2_fu_bt.sam") into alignment_results
		set dataset_id, file("${dataset_id}_R1_bt.log"), file("${dataset_id}_R2_bt.log") into bt_log
		set dataset_id, file("${dataset_id}_R2_fuh.fastq"), file("${dataset_id}_R1_fuh.fastq") into host_filted_results
		
	"""
	bowtie2 -x $params.btindex --local -q -U ${dataset_id}_R1_fu.fastq --sensitive --score-min C,${params.min_score},0 --time -p $params.num_cpu -S ${dataset_id}_R1_fu_bt.sam \
	2> ${dataset_id}_R1_bt.log
	
	$baseDir/bin/fasta_from_sam -f ${dataset_id}_R1_fu.fastq -r ${dataset_id}_R1_fu_bt.sam > ${dataset_id}_R1_fuh1.fastq

	$baseDir/bin/reconcile_read2_file ${dataset_id}_R1_fuh1.fastq ${dataset_id}_R2_fu.fastq f2 > ${dataset_id}_R2_fuh1.fastq
	
	bowtie2 -x $params.btindex --local -q -U ${dataset_id}_R2_fu.fastq --sensitive --score-min C,${params.min_score},0 --time -p $params.num_cpu -S ${dataset_id}_R2_fu_bt.sam \
	2> ${dataset_id}_R2_bt.log

	$baseDir/bin/fasta_from_sam -f ${dataset_id}_R2_fuh1.fastq -r ${dataset_id}_R2_fu_bt.sam > ${dataset_id}_R2_fuh.fastq

	$baseDir/bin/reconcile_read2_file ${dataset_id}_R2_fuh.fastq ${dataset_id}_R1_fuh1.fastq > ${dataset_id}_R1_fuh.fastq
	"""
}	


process runSPAdes {
	
	publishDir "${params.output}/SPAdes_Results", mode: 'link'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_R2_fuh.fastq"), file("${dataset_id}_R1_fuh.fastq") from host_filted_results
	
	output:
		set dataset_id, file("${dataset_id}.spades") into SPAdes_results
		
	"""
	spades.py -o ${dataset_id}.spades --pe1-1 ${dataset_id}_R1_fuh.fastq --pe1-2 ${dataset_id}_R2_fuh.fastq -t 24 -m 150
	"""
}

process filterContigs {
	
	publishDir "${params.output}/SPAdes_Results", mode: 'link'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}.spades") from SPAdes_results
		
	output:
		set dataset_id, file("${dataset_id}_spade_contigs_gt_150.fa") into filtered_contig_fastas
		
	"""
	$baseDir/bin/fasta_to_fasta ${dataset_id}.spades/contigs.fasta > ./${dataset_id}_spade_contigs.fa
	
	$baseDir/bin/filter_fasta_by_size -b 150 ${dataset_id}_spade_contigs.fa | $baseDir/bin/fasta_to_fasta >  ${dataset_id}_spade_contigs_gt_150.fa
	"""
}	
