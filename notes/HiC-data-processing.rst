======================================
HiC Data Processing and Analysis
======================================


Reference genome and annotations
=================================

The following were used:


* ``danRer11`` with ``alt`` contigs removed;

* ENSEMBL gene models;



Raw data processing
=======================

* ``Fastq`` files were processed and summarised using ``Juicer`` 1.6 (https://github.com/aidenlab/juicer). ``BWA 0.7.17-r1188`` was used for read mapping.

* Files ``juicer.sh`` and ``split_rmdups.awk`` were modified to fit job scheduling configuration and module system on HPC Rackham (`Uppmax <https://www.uu.se/centrum/uppmax/>`_); the modified code is at https://github.com/agata-sm/juicer-rackham.

* To generate merged matrices (megamap in Juicer nomenclature) modified ``mega.sh`` was used; the modified code is at https://github.com/agata-sm/juicer-rackham.


Data analysis
=================



TADs
-------


* TADs were detected using ``OnTAD`` v. 1.4 (gcc 8.4.0), on merged replicate matrices, 10 kb bin size.

	Parameters used:

	* ``-lsize`` Lsize 7 (The local region size that used to determine local minimum)

	* ``-penalty`` 0.075 (The penalty applied in scoring function to select positive TADs. Higher penalty score will result in fewer TADs)

	* ``-maxsz`` 300 (The maximum size of TADs can be called. The size is determined by number of bins covered in the contact matrix)


* Consensus and differential TAD boundaries between matrices (10 kb bin size) were identified using ``TADCompare``.
	
	* Consensus: between individual replicate matrices (within tissue): tissue specific OnTAD boundaries were queried;

	* Differential: between merged replicate matrices (between tissues): union of tissue specific OnTAD boundaries was queried;


Loops
-------

* Loops were detected using ``hiccups`` GPU::

	java -jar /proj/snic2019-30-13/nobackup/nbis5931/scripts/juicer_2023/juicer_tools.2.20.00.jar hiccups --threads 20 -r 5000,10000,25000 -k KR -f .1,.1,.1 -p 4,2,1 -i 7,5,3 -t 0.02,1.5,1.75,2 -d 20000,20000,50000 $infile_hic $outdir1

* Loops consistent between replicates were identified::

	java -jar /proj/snic2019-30-13/nobackup/juicer/scripts/common/juicer_tools.jar compare -m 25000 0 $chromsizes $loops102 $loops104 102_vs_104

* Loops common between LEC and BEC: at least 0.4 of the loop region length from both lists must overlap::

	bedtools pairtopair -f 0.4 -a $loops_LEC -b $loops_BEC > LEC_vs_BEC.common.0_4.overlap.bedpe

	bedtools pairtopair -f 0.4 -a $loops_BEC -b $loops_LEC > BEC_vs_LEC.common.0_4.overlap.bedpe

	#output files were concatenated and further filtered to remove redundancies

* Loops different between LEC and BEC (tissue-specific): loops within at least 40kb distance from each other; at least one end of A mustn't overlap B::

	bedtools pairtopair -type notboth -slop 40000 -a $loops_LEC -b $loops_BEC > LEC_vs_BEC.common.nooverlap_40kb.bedpe

	bedtools pairtopair -type notboth -slop 40000 -a $loops_BEC -b $loops_LEC > BEC_vs_LEC.common.nooverlap_40kb.bedpe


