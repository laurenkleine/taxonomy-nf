#!/usr/bin/env nextflow

params.reads = "$baseDir/data/read_pairs/*_{1,2}.fastq"

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_pairs_ch } 

process fastqc {
    tag "FASTQC on $sample_id"
    
    input:
    set sample_id, file(reads) from read_pairs_ch
   
    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch
    
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
} 