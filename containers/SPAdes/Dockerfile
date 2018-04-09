FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>


RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    python-pip \
    libpython2.7-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    python2.7-minimal \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*


RUN wget http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz -O /opt/SPAdes-3.5.0-Linux.tar.gz && \
    tar -xzf /opt/SPAdes-3.5.0-Linux.tar.gz -C /opt/ && \
    cd /opt/SPAdes-3.5.0-Linux/bin/ && \
    ln -s /opt/SPAdes-3.5.0-Linux/bin/* /usr/local/bin && \
    rm -rf /opt/SPAdes-3.5.0-Linux.tar.gz
    
########################################################################################

N E X T F L O W  ~  version 0.28.1
Launching `main_bt.nf` [boring_franklin] - revision: d9f57c93d3
[warm up] executor > local
[70/dcbd8f] Cached process > RunCutAdapt (Croc)
[ee/4fbdac] Cached process > RunPreFastQC (Croc)
[65/13bba8] Cached process > cd_hit_est (Croc)
[32/5749c7] Cached process > reconcile_reads (Croc)
[b6/9481aa] Cached process > runBowtie2 (Croc)
[79/74a66b] Cached process > RunPostFastQC (Croc)
[26/8d14f4] Submitted process > runSPAdes (Croc)
ERROR ~ Error executing process > 'runSPAdes (Croc)'

Caused by:
  Process `runSPAdes (Croc)` terminated with an error exit status (1)

Command executed:

  spades.py -o Croc.spades --pe1-1 Croc_R1_fuh.fastq --pe1-2 Croc_R2_fuh.fastq -t 24 -m 150
  
  /home/lauren/Desktop/taxonomy-nf/bin/fasta_to_fasta Croc.spades/contigs.fasta > ./Croc_spade_contigs.fa

Command exit status:
  1

Command output:
     0:00:01.066   16M /    3G   INFO  General                 (kmer_index.hpp            : 353)   Merging temporary buckets.
     0:00:01.067   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 476)   Building perfect hash indices
     0:00:01.068   16M /    3G   INFO  General                 (kmer_index.hpp            : 371)   Merging final buckets.
     0:00:01.069   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 515)   Index built. Total 1228 bytes occupied (2.86916 bits per kmer).
     0:00:01.069   16M /    3G   INFO DeBruijnExtensionIndexBu (kmer_extension_index.hpp  : 382)   Building k-mer extensions from k+1-mers
     0:00:01.070   16M /    3G   INFO DeBruijnExtensionIndexBu (kmer_extension_index.hpp  : 386)   Building k-mer extensions from k+1-mers finished.
     0:00:01.070   16M /    3G   INFO Early tip clipping       (early_simplification.hpp  : 236)   Early tip clipping
     0:00:01.071   16M /    3G   INFO Early tip clipping       (early_simplification.hpp  : 239)   0 22-mers were removed by early tip clipper
     0:00:01.071   16M /    3G   INFO  General                 (graph_construction.hpp    : 145)   Condensing graph
     0:00:01.071   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 337)   Extracting unbranching paths
     0:00:01.072   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 351)   Extracting unbranching paths finished. 51 sequences extracted
     0:00:01.072   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 311)   Collecting perfect loops
     0:00:01.073   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 326)   Collecting perfect loops finished. 0 loops collected
     0:00:01.073   16M /    3G   INFO  General                 (graph_construction.hpp    : 150)   Building index with from graph
     0:00:01.073   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 467)   Building kmer index
     0:00:01.073   16M /    3G   INFO  General                 (kmer_splitters.hpp        : 289)   Splitting kmer instances into 16 buckets. This might take a while.
     0:00:01.073   16M /    3G   INFO  General                 (file_limit.hpp            :  29)   Open file limit set to 1048576
     0:00:01.074   16M /    3G   INFO  General                 (kmer_splitters.hpp        : 307)   Memory available for splitting buffers: 49.9948 Gb
     0:00:01.074   16M /    3G   INFO  General                 (kmer_splitters.hpp        : 312)   Using cell size of 4194304
     0:00:01.075  592M /    3G   INFO  General                 (kmer_splitters.hpp        : 330)   Used 51 sequences.
     0:00:01.075   16M /    3G   INFO  General                 (kmer_index.hpp            : 345)   Starting k-mer counting.
     0:00:01.076   16M /    3G   INFO  General                 (kmer_index.hpp            : 351)   K-mer counting done. There are 3403 kmers in total.
     0:00:01.076   16M /    3G   INFO  General                 (kmer_index.hpp            : 353)   Merging temporary buckets.
     0:00:01.077   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 476)   Building perfect hash indices
     0:00:01.077   16M /    3G   INFO  General                 (kmer_index.hpp            : 371)   Merging final buckets.
     0:00:01.078   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 515)   Index built. Total 1217 bytes occupied (2.861 bits per kmer).
     0:00:01.078   16M /    3G   INFO  General                 (edge_index_builders.hpp   :  21)   Collecting k-mer coverage information from graph, this takes a while.
     0:00:01.079   16M /    3G   INFO  General                 (graph_construction.hpp    : 182)   Filling coverage index
     0:00:01.079   16M /    3G   INFO  General                 (edge_index_builders.hpp   :  99)   Collecting k-mer coverage information from reads, this takes a while.
     0:00:01.080   16M /    3G   INFO  General                 (graph_construction.hpp    : 184)   Filling coverage and flanking coverage from index
     0:00:01.080   16M /    3G   INFO  General                 (construction.cpp          :  33)   Figured out: read length = 150
     0:00:01.080   16M /    3G   INFO StageManager             (stage.cpp                 : 125)   STAGE == EC Threshold Finding
  === Stack Trace ===
  [0x406f4a]
  [0x4012c2]
  [0x5ee96b]
  [0x431c0c]
  [0x614fc2]
  [0x626258]
  [0x7b3ea0]
  [0x4001b9]
  Verification of expression 'cov_.size() > 10' failed in function 'void cov_model::KMerCoverageModel::Fit()'. In file '/home/yasha/gitrep2/algorithmic-biology/assembler/src/debruijn/kmer_coverage_model.cpp' on line 184. Message 'Invalid kmer coverage histogram'.
  Verification of expression 'cov_.size() > 10' failed in function 'void cov_model::KMerCoverageModel::Fit()'. In file '/home/yasha/gitrep2/algorithmic-biology/assembler/src/debruijn/kmer_coverage_model.cpp' on line 184. Message 'Invalid kmer coverage histogram'.
  spades: /home/yasha/gitrep2/algorithmic-biology/assembler/src/debruijn/kmer_coverage_model.cpp:184: void cov_model::KMerCoverageModel::Fit(): Assertion `cov_.size() > 10' failed.
  
  
  == Error ==  system call for: "['/opt/SPAdes-3.5.0-Linux/bin/spades', 'Croc.spades/K21/configs/config.info']" finished abnormally, err code: -6
  
  In case you have troubles running SPAdes, you can write to spades.support@bioinf.spbau.ru
  Please provide us with params.txt and spades.log files from the output directory.

Command wrapper:
     0:00:01.066   16M /    3G   INFO  General                 (kmer_index.hpp            : 353)   Merging temporary buckets.
     0:00:01.067   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 476)   Building perfect hash indices
     0:00:01.068   16M /    3G   INFO  General                 (kmer_index.hpp            : 371)   Merging final buckets.
     0:00:01.069   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 515)   Index built. Total 1228 bytes occupied (2.86916 bits per kmer).
     0:00:01.069   16M /    3G   INFO DeBruijnExtensionIndexBu (kmer_extension_index.hpp  : 382)   Building k-mer extensions from k+1-mers
     0:00:01.070   16M /    3G   INFO DeBruijnExtensionIndexBu (kmer_extension_index.hpp  : 386)   Building k-mer extensions from k+1-mers finished.
     0:00:01.070   16M /    3G   INFO Early tip clipping       (early_simplification.hpp  : 236)   Early tip clipping
     0:00:01.071   16M /    3G   INFO Early tip clipping       (early_simplification.hpp  : 239)   0 22-mers were removed by early tip clipper
     0:00:01.071   16M /    3G   INFO  General                 (graph_construction.hpp    : 145)   Condensing graph
     0:00:01.071   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 337)   Extracting unbranching paths
     0:00:01.072   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 351)   Extracting unbranching paths finished. 51 sequences extracted
     0:00:01.072   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 311)   Collecting perfect loops
     0:00:01.073   16M /    3G   INFO UnbranchingPathExtractor (debruijn_graph_constructor: 326)   Collecting perfect loops finished. 0 loops collected
     0:00:01.073   16M /    3G   INFO  General                 (graph_construction.hpp    : 150)   Building index with from graph
     0:00:01.073   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 467)   Building kmer index
     0:00:01.073   16M /    3G   INFO  General                 (kmer_splitters.hpp        : 289)   Splitting kmer instances into 16 buckets. This might take a while.
     0:00:01.073   16M /    3G   INFO  General                 (file_limit.hpp            :  29)   Open file limit set to 1048576
     0:00:01.074   16M /    3G   INFO  General                 (kmer_splitters.hpp        : 307)   Memory available for splitting buffers: 49.9948 Gb
     0:00:01.074   16M /    3G   INFO  General                 (kmer_splitters.hpp        : 312)   Using cell size of 4194304
     0:00:01.075  592M /    3G   INFO  General                 (kmer_splitters.hpp        : 330)   Used 51 sequences.
     0:00:01.075   16M /    3G   INFO  General                 (kmer_index.hpp            : 345)   Starting k-mer counting.
     0:00:01.076   16M /    3G   INFO  General                 (kmer_index.hpp            : 351)   K-mer counting done. There are 3403 kmers in total.
     0:00:01.076   16M /    3G   INFO  General                 (kmer_index.hpp            : 353)   Merging temporary buckets.
     0:00:01.077   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 476)   Building perfect hash indices
     0:00:01.077   16M /    3G   INFO  General                 (kmer_index.hpp            : 371)   Merging final buckets.
     0:00:01.078   16M /    3G   INFO K-mer Index Building     (kmer_index.hpp            : 515)   Index built. Total 1217 bytes occupied (2.861 bits per kmer).
     0:00:01.078   16M /    3G   INFO  General                 (edge_index_builders.hpp   :  21)   Collecting k-mer coverage information from graph, this takes a while.
     0:00:01.079   16M /    3G   INFO  General                 (graph_construction.hpp    : 182)   Filling coverage index
     0:00:01.079   16M /    3G   INFO  General                 (edge_index_builders.hpp   :  99)   Collecting k-mer coverage information from reads, this takes a while.
     0:00:01.080   16M /    3G   INFO  General                 (graph_construction.hpp    : 184)   Filling coverage and flanking coverage from index
     0:00:01.080   16M /    3G   INFO  General                 (construction.cpp          :  33)   Figured out: read length = 150
     0:00:01.080   16M /    3G   INFO StageManager             (stage.cpp                 : 125)   STAGE == EC Threshold Finding
  === Stack Trace ===
  [0x406f4a]
  [0x4012c2]
  [0x5ee96b]
  [0x431c0c]
  [0x614fc2]
  [0x626258]
  [0x7b3ea0]
  [0x4001b9]
  Verification of expression 'cov_.size() > 10' failed in function 'void cov_model::KMerCoverageModel::Fit()'. In file '/home/yasha/gitrep2/algorithmic-biology/assembler/src/debruijn/kmer_coverage_model.cpp' on line 184. Message 'Invalid kmer coverage histogram'.
  Verification of expression 'cov_.size() > 10' failed in function 'void cov_model::KMerCoverageModel::Fit()'. In file '/home/yasha/gitrep2/algorithmic-biology/assembler/src/debruijn/kmer_coverage_model.cpp' on line 184. Message 'Invalid kmer coverage histogram'.
  spades: /home/yasha/gitrep2/algorithmic-biology/assembler/src/debruijn/kmer_coverage_model.cpp:184: void cov_model::KMerCoverageModel::Fit(): Assertion `cov_.size() > 10' failed.
  
  
  == Error ==  system call for: "['/opt/SPAdes-3.5.0-Linux/bin/spades', 'Croc.spades/K21/configs/config.info']" finished abnormally, err code: -6
  
  In case you have troubles running SPAdes, you can write to spades.support@bioinf.spbau.ru
  Please provide us with params.txt and spades.log files from the output directory.

Work dir:
  /home/lauren/Desktop/taxonomy-nf/work/26/8d14f444c98fb0377a2b52f76d216f