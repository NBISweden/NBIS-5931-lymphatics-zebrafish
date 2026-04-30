#!/bin/bash


# for BINDdetect2nodes

annot_dir=""

peaks_file="/path/to/ATAC_peaks_annot_closest.ensembl.consensus_replicate_merged_MACS2_broad_TMM_paired.14xi2023.tsv"
	
binddetect="/path/to//BINDetect_jaspar_mafb/BINDetect_mafb_all_jaspar/LEC_BEC_diff/BINDdetect"

motifs_dir="/path/to/motifs_for_network"

# for select_nodes

dir_rnaseq="/path/to/TF_network/rnaseq"
	
file_abundance="average_expression_5days.csv"
	
file_lFC="240202_fold_change_5days.csv"
	
file_GO="/path/to//path/to/genes_endothelium.txt"


# outdir
outdir="/path/to/mafb_altmotifs_subnetworks"


# RNAseq cutoffs

rnaseq_lfc="0.5"; #cutoff for log2FC of RNA-seq abundance

rnaseq_frac1="0.2"; #cutoff for pct cells expressing a gene
rnaseq_frac2="0.3"; #cutoff for pct cells expressing a gene


# TF for network generation

TF_tst="mafba"

mkdir -p $outdir

######

motif_files="${motifs_dir}/*txt"

for motifs_file in $motif_files
do

	motifs_fname=${motifs_file##*/}

	echo "processing file $motifs_fname"

	motifs_pth="${motifs_dir}/${motifs_fname}"
	motif_base="${motifs_fname%%.txt}"
	resdir="${outdir}/${motif_base}" 
	outfile_Bd2nodes="${resdir}/LEC_BEC_TF_0."$motif_base".footprints_to_nodes.tsv"


	perl BINDdetect2nodes.v0.7.pl --infile_motifs ${motifs_pth} \
	--infile_ortho ${annot_dir}/Drer_hsa.orthology-alliance.tsv \
	--annotated_peaks $peaks_file \
	--indir_BINDdetect $binddetect \
	--outdir $resdir \
	--pref ${motif_base}

	echo "creating subnetworks"

	for rnaseq_frac in $rnaseq_frac1 $rnaseq_frac2
	do
		echo "RNA frac expr co $rnaseq_frac"
	
		outdir_net="${resdir}/subnetworks/${TF_tst}/RNAseq_lfc${rnaseq_lfc}_expr${rnaseq_frac}"

		perl select_nodes_v0.5.pl --infile_edges $outfile_Bd2nodes\
		 --infile_abundance ${dir_rnaseq}/${file_abundance} \
		 --infile_log2FC ${dir_rnaseq}/${file_lFC} \
		 --TF $TF_tst \
		 --outdir $outdir_net \
		 --pref ${motif_base} \
		 --rnaseq_lfc ${rnaseq_lfc} \
		 --rnaseq_frac ${rnaseq_frac} \
		 --smpl LEC \
		 --infile_genelist $file_GO

		perl select_nodes_v0.5.pl --infile_edges $outfile_Bd2nodes\
		 --infile_abundance ${dir_rnaseq}/${file_abundance} \
		 --infile_log2FC ${dir_rnaseq}/${file_lFC} \
		 --TF $TF_tst \
		 --outdir $outdir_net \
		 --pref ${motif_base} \
		 --rnaseq_lfc ${rnaseq_lfc} \
		 --rnaseq_frac ${rnaseq_frac} \
		 --smpl BEC \
		 --infile_genelist $file_GO

		perl select_nodes_v0.5.pl --infile_edges $outfile_Bd2nodes\
		 --infile_abundance ${dir_rnaseq}/${file_abundance} \
		 --infile_log2FC ${dir_rnaseq}/${file_lFC} \
		 --TF $TF_tst \
		 --outdir $outdir_net \
		 --pref ${motif_base} \
		 --rnaseq_lfc ${rnaseq_lfc} \
		 --rnaseq_frac ${rnaseq_frac} \
		 --smpl LEC 

		 perl select_nodes_v0.5.pl --infile_edges $outfile_Bd2nodes\
		 --infile_abundance ${dir_rnaseq}/${file_abundance} \
		 --infile_log2FC ${dir_rnaseq}/${file_lFC} \
		 --TF $TF_tst \
		 --outdir $outdir_net \
		 --pref ${motif_base} \
		 --rnaseq_lfc ${rnaseq_lfc} \
		 --rnaseq_frac ${rnaseq_frac} \
		 --smpl BEC
	done

done

echo "*** all done ***"

