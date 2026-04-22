======================================
ATAC-seq Data Processing and Analysis
======================================


Reference genome and annotations
=================================

The following were used:


* ``danRer11`` with ``alt`` contigs removed;

* ENSEMBL gene models;

* TF binding motifs were ``JASPAR2022_CORE_vertebrates_non-redundant``


Raw data processing
=======================

* ``Fastq`` files were processed using ``nf-core ataseq`` pipeline (1.2.1) (https://nf-co.re/atacseq/2.1.1/);

* ``bam`` files were filtered to remove fragments mapped to MT chromosome;


Data analysis
=======================

* peaks were called using ``MACS2`` (2.2.6);

* peaks detected in at least two replicates per tissue were retained in the analysis;

* TF footprinting was performed on bam files merged per tissue using TOBIAS (0.16.0) (https://github.com/loosolab/TOBIAS)

* 


