#!/usr/bin/env nextflow

/*  
	
	To run this practice nf script with Docker: 
	
	nextflow run main_practice_docker.nf -profile docker

*/

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


/*dumby process*/

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
        fastqc -f fastq ${dataset_id}_R1_fu.fastq {dataset_id}_R2_fu.fastq -o temp
        mv temp/*.{html,zip} .
        """
}
