#!/usr/bin/env nextflow

min_len = params.min_len

btindex = params.btindex

num_cpu = params.num_cpu
min_score = params.min_score

db = params.db

/* fix DIAMOND index to work correctly */



Channel
	.fromFilePairs( params.reads, flat: true) /*When true the matching files are produced as sole elements in the emitted tuples*/
	.ifEmpty { return "Error" }
	.into { read_pairs; qc_pairs }


/*process that establishes correct fastq version*/

process RunPreFastQC {
	
	publishDir "${params.output}/FastQCResults/Pre", mode: 'link'
	
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
	
	publishDir "${params.output}/CutAdaptResults", mode: 'link'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file(forward), file(reverse) from (read_pairs)
	output:
		set dataset_id, file("${dataset_id}.R1.fastq"), file("${dataset_id}.R2.fastq") into (cut_adapt_results1, cut_adapt_results2)
	"""
	cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -A AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT -G CTGTCTCTTATACACATCT -q 30,30 --minimum-length ${min_len} -u 1 -o ${dataset_id}.R1.fastq -p ${dataset_id}.R2.fastq $forward $reverse
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
	
	publishDir "${params.output}/FastQCResults/Post", mode: 'link'

	tag { dataset_id }

	input:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") from reconcile_reads_results1
	
	output:
		set dataset_id, file('*_fastqc.{html,zip}') into fastqc_results_post

	"""
	mkdir temp
	fastqc -f fastq ${dataset_id}_R1_fu.fastq ${dataset_id}_R2_fu.fastq -o temp
	mv temp/*.{html,zip} .
	"""
}


process runBowtie2 {
	
	publishDir "${params.output}/Bowtie2_Results", mode: 'link'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_R1_fu.fastq"), file("${dataset_id}_R2_fu.fastq") from reconcile_reads_results2
		
	output:
		set dataset_id, file("${dataset_id}_R1_fu_bt.sam"), file("${dataset_id}_R2_fu_bt.sam") into alignment_results
		set dataset_id, file("${dataset_id}_R1_bt.log"), file("${dataset_id}_R2_bt.log") into bt_log
		set dataset_id, file("${dataset_id}_R2_fuh.fastq"), file("${dataset_id}_R1_fuh.fastq") into (host_filted_results1, host_filted_results2, host_filted_results3)
		
	"""
	bowtie2 -x ${btindex} --local -q -U ${dataset_id}_R1_fu.fastq --sensitive --score-min C,${min_score},0 --time -p ${num_cpu} -S ${dataset_id}_R1_fu_bt.sam \
	2> ${dataset_id}_R1_bt.log
	
	$baseDir/bin/fasta_from_sam -f ${dataset_id}_R1_fu.fastq -r ${dataset_id}_R1_fu_bt.sam > ${dataset_id}_R1_fuh1.fastq

	$baseDir/bin/reconcile_read2_file ${dataset_id}_R1_fuh1.fastq ${dataset_id}_R2_fu.fastq f2 > ${dataset_id}_R2_fuh1.fastq
	
	bowtie2 -x ${btindex} --local -q -U ${dataset_id}_R2_fu.fastq --sensitive --score-min C,${min_score},0 --time -p ${num_cpu} -S ${dataset_id}_R2_fu_bt.sam \
	2> ${dataset_id}_R2_bt.log

	$baseDir/bin/fasta_from_sam -f ${dataset_id}_R2_fuh1.fastq -r ${dataset_id}_R2_fu_bt.sam > ${dataset_id}_R2_fuh.fastq

	$baseDir/bin/reconcile_read2_file ${dataset_id}_R2_fuh.fastq ${dataset_id}_R1_fuh1.fastq > ${dataset_id}_R1_fuh.fastq
	"""
}	


process runSPAdes {
	
	publishDir "${params.output}/SPAdes_Results", mode: 'link'
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_R2_fuh.fastq"), file("${dataset_id}_R1_fuh.fastq") from host_filted_results1
	
	output:
		set dataset_id, file("${dataset_id}.spades") into SPAdesDir
		set dataset_id, file("${dataset_id}_spade_contigs.fa") into (SPAdes_results1, SPAdes_results2, SPAdes_results3, SPAdes_results4)
		
	"""
	spades.py -o ${dataset_id}.spades --pe1-1 ${dataset_id}_R1_fuh.fastq --pe1-2 ${dataset_id}_R2_fuh.fastq -t 24 -m 150
	
	$baseDir/bin/fasta_to_fasta ${dataset_id}.spades/contigs.fasta > ./${dataset_id}_spade_contigs.fa
	"""
}


process contigBowtie {

	publishDir "${params.output}/contigBowtie_Results", mode: "link"
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs.fa") from (SPAdes_results1)
		set dataset_id, file("${dataset_id}_R2_fuh.fastq"), file("${dataset_id}_R1_fuh.fastq") from host_filted_results2
		
	
	output:
		set dataset_id, file("${dataset_id}_spade_contigs_R1.sam"), file("${dataset_id}_spade_contigs_R2.sam") into (fuh_sam_results1, fuh_sam_results2)
		set dataset_id, file("${dataset_id}_spade_contigs_R1.bt_log"), file("${dataset_id}_spade_contigs_R2.bt_log") into fuh_btlog_results
	
	"""
	bowtie2-build ${dataset_id}_spade_contigs.fa ${dataset_id}_spade_contigs
	bowtie2 -x ${dataset_id}_spade_contigs --local --score-min C,120,1 -q -U ${dataset_id}_R1_fuh.fastq -p 12 -S ${dataset_id}_spade_contigs_R1.sam 2> ${dataset_id}_spade_contigs_R1.bt_log
	bowtie2 -x ${dataset_id}_spade_contigs --local --score-min C,120,1 -q -U ${dataset_id}_R2_fuh.fastq -p 12 -S ${dataset_id}_spade_contigs_R2.sam 2> ${dataset_id}_spade_contigs_R2.bt_log
	"""
}


process tallySAMsubjects {

	publishDir "${params.output}/tallySAMsubjects_Results", mode: "link"
	
	tag { dataset_id }

	input:
		set dataset_id, file("${dataset_id}_spade_contigs_R1.sam"), file("${dataset_id}_spade_contigs_R2.sam") from (fuh_sam_results1)
		
	
	output:
		set dataset_id, file("${dataset_id}_spade_contigs_R1.sam.subject_hits") into (sam_subject_hits_resultsR1_1, sam_subject_hits_resultsR1_2)
		set dataset_id, file("${dataset_id}_spade_contigs_R2.sam.subject_hits") into (sam_subject_hits_resultsR2)
	
	"""
	$baseDir/bin/tally_sam_subjects ${dataset_id}_spade_contigs_R1.sam > ${dataset_id}_spade_contigs_R1.sam.subject_hits
	$baseDir/bin/tally_sam_subjects ${dataset_id}_spade_contigs_R2.sam > ${dataset_id}_spade_contigs_R2.sam.subject_hits
	"""
}	
	
process compositeSAMfile {

	publishDir "${params.output}/compositeSAMfile_Results", mode: "link"
	
	tag {dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs_R1.sam"), file("${dataset_id}_spade_contigs_R2.sam") from (fuh_sam_results2)
		set dataset_id, file("${dataset_id}_R2_fuh.fastq"), file("${dataset_id}_R1_fuh.fastq") from host_filted_results3
		
		
	output:
		set dataset_id, file("${dataset_id}_R1_fuhs.fastq"), file("${dataset_id}_R2_fuhs.fastq") into compositeSAMfile_results
		
	"""
	cat ${dataset_id}_spade_contigs_R1.sam ${dataset_id}_spade_contigs_R2.sam > ${dataset_id}_R12_fuh.fastq.sam
	$baseDir/bin/fasta_from_sam -r -f ${dataset_id}_R1_fuh.fastq ${dataset_id}_R12_fuh.fastq.sam > ${dataset_id}_R1_fuhs.fastq
	$baseDir/bin/fasta_from_sam -r -f ${dataset_id}_R2_fuh.fastq ${dataset_id}_R12_fuh.fastq.sam > ${dataset_id}_R2_fuhs.fastq
	"""
	
}

process BLASTcontigs_vs_nt {
	
	publishDir "${params.output}/BLAST_Results", mode: "link"
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs.fa") from (SPAdes_results2)
	
	output:
		set dataset_id, file("${dataset_id}.fa.bn_nt") into (BLASTresults1, BLASTresults2)
		set dataset_id, file("${dataset_id}_spade_contigs_n.fa") into (BLASTresults)
		
	"""
	blastn -query ${dataset_id}_spade_contigs.fa -db ${db} -num_threads 24 -evalue 1e-8 -task megablast -outfmt 6 | $baseDir/bin/consolidate_blast_output > ${dataset_id}.fa.bn_nt
	
	$baseDir/bin/fasta_from_blast -f ${dataset_id}_spade_contigs.fa -r ${dataset_id}.fa.bn_nt > ${dataset_id}_spade_contigs_n.fa
	"""
}

process taxonomic_Assessment_nt {

	publishDir ".", mode: "link"
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs_R1.sam.subject_hits") from (sam_subject_hits_resultsR1_1)
		set dataset_id, file("${dataset_id}.fa.bn_nt") from BLASTresults1
		
	output:
		set dataset_id, file("${dataset_id}.fa.bn_nt.tally") into bn_nt_tally_results
		set dataset_id, file("${dataset_id}.bn_nt.desc_tally") into desc_tally_results
		set dataset_id, file("${dataset_id}.bn_nt.tab_tree_tally") into tab_tree_tally_results
	
	"""
	$baseDir/bin/tally_hits_universal.pl -lca -w ${dataset_id}_spade_contigs_R1.sam.subject_hits ${dataset_id}.fa.bn_nt > ${dataset_id}.fa.bn_nt.tally
	$baseDir/bin/tally_hits_universal.pl -lca -w ${dataset_id}_spade_contigs_R1.sam.subject_hits -d -o desc_tally ${dataset_id}.fa.bn_nt > ${dataset_id}.bn_nt.desc_tally
	$baseDir/bin/tally_hits_universal.pl -lca -w ${dataset_id}_spade_contigs_R1.sam.subject_hits -t -ti -o tab_tree_tally ${dataset_id}.fa.bn_nt > ${dataset_id}.bn_nt.tab_tree_tally
	"""
}	

/* my $tax_dir = "/home/databases/NCBI_Taxonomy"; */

process virusDerivedReads_nt {

	publishDir ".", mode: "link"
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs.fa") from (SPAdes_results3)
		set dataset_id, file("${dataset_id}.fa.bn_nt") from BLASTresults2
	
	output:
		set dataset_id, file('*.fa') into reads_results
		
	"""
	$baseDir/bin/distribute_fasta_by_blast_taxid.pl -v ${dataset_id}_spade_contigs.fa ${dataset_id}.fa.bn_nt
	"""
}
	
process runDIAMOND {
	
	publishDir "${params.output}/DIAMOND_Results", mode: "link"
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs_n.fa") from (BLASTresults)
	
	output:
		set dataset_id, file("${dataset_id}_spade_contigs_n.fa.dmd_nr") into runDIAMOND_results1, runDIAMOND_results2
		set dataset_id, file("${dataset_id}_spade_contigs_nn.fa") into fastaDIAMOND_results1
		
	"""
	diamond blastx --db /home/databases/nr_nt/nr.dmnd --threads 16 --out ${dataset_id}_spade_contigs_n.fa.dmd_nr --outfmt 6 --query ${dataset_id}_spade_contigs_n.fa --unal 0 --evalue 1e-3
	$baseDir/bin/fasta_from_blast -r -f ${dataset_id}_spade_contigs_n.fa ${dataset_id}_spade_contigs_n.fa.dmd_nr > ${dataset_id}_spade_contigs_nn.fa
	"""
}

process taxonomic_Assessment_nr {

	publishDir ".", mode: "link"
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs_R1.sam.subject_hits") from (sam_subject_hits_resultsR1_2)
		set dataset_id, file("${dataset_id}_spade_contigs_n.fa.dmd_nr") from runDIAMOND_results1
		
	output:
		set dataset_id, file("${dataset_id}.dmd_nr.tally") into nr_tally_results
		set dataset_id, file("${dataset_id}.dmd_nr.desc_tally") into nr_desc_tally_results
		set dataset_id, file("${dataset_id}.dmd_nr.tab_tree_tally") into nr_tab_tree_tally_results
		
	"""
	$baseDir/bin/tally_hits_universal.pl -lca -w ${dataset_id}_spade_contigs_R1.sam.subject_hits ${dataset_id}_spade_contigs_n.fa.dmd_nr > ${dataset_id}.dmd_nr.tally
	$baseDir/bin/tally_hits_universal.pl -lca -w ${dataset_id}_spade_contigs_R1.sam.subject_hits -d -o desc_tally ${dataset_id}_spade_contigs_n.fa.dmd_nr > ${dataset_id}.dmd_nr.desc_tally
	$baseDir/bin/tally_hits_universal.pl -lca -w ${dataset_id}_spade_contigs_R1.sam.subject_hits -t -ti -o tab_tree_tally ${dataset_id}_spade_contigs_n.fa.dmd_nr > ${dataset_id}.dmd_nr.tab_tree_tally
	"""
}

process virusDerivedReads_nr {

	publishDir ".", mode: "link"
	
	tag { dataset_id }
	
	input:
		set dataset_id, file("${dataset_id}_spade_contigs.fa") from (SPAdes_results4)
		set dataset_id, file("${dataset_id}_spade_contigs_n.fa.dmd_nr") from runDIAMOND_results2
	
	output:
		set dataset_id, file("*.fa") into final_channel
	
	"""
	$baseDir/bin/distribute_fasta_by_blast_taxid.pl -v ${dataset_id}_spade_contigs.fa ${dataset_id}_spade_contigs_n.fa.dmd_nr
	"""
}
