#!c:/perl/bin/perl.exe


# script to parse output of BINDetect (TOBIAS) to get data required for TF network reconstruction

# output sites bound in only one of the sample and sites which pass the log2FC cutoff if bound in both samples

# input files

#######
# from Jaspar, processed by TOBIAS into one file

# all_motifs.txt
# >MA0804.1	MA0804.1.TBX19

# also this version exists in directories of TOBIAS output

# all_motifs.txt
# >MA0737.1	GLIS3



#######
# from https://www.alliancegenome.org/downloads and further parsed to contain hsa and drer only

# Drer_hsa.orthology-alliance.tsv

# # Species: Homo sapiens, Rattus norvegicus, Mus musculus, Danio rerio, Xenopus tropicalis, Xenopus laevis, Drosophila melanogaster, Caenorhabditis elegans, Saccharomyces cerevisiae
# HGNC:4232	GDNF	NCBITaxon:9606	Homo sapiens	ZFIN:ZDB-GENE-010226-1	gdnfa	NCBITaxon:7955	Danio rerio	OrthoInspector|OMA|OrthoFinder|Hieranoid|SonicParanoid|ZFIN|PhylomeDB|PANTHER|InParanoid	9	10	Yes	Yes
# HGNC:5105	HOXA4	NCBITaxon:9606	Homo sapiens	ZFIN:ZDB-GENE-000823-4	hoxa4a	NCBITaxon:7955	Danio rerio	Ensembl Compara|ZFIN|PANTHER|OrthoFinder	4	10	Yes	Yes
# HGNC:2384	CRY1	NCBITaxon:9606	Homo sapiens	ZFIN:ZDB-GENE-010426-3	cry1b	NCBITaxon:7955	Danio rerio	OrthoInspector|OrthoFinder|SonicParanoid|ZFIN|PhylomeDB|PANTHER|Ensembl Compara|InParanoid	8	10	No	Yes
# HGNC:400	ALCAM	NCBITaxon:9606	Homo sapiens	ZFIN:ZDB-GENE-990415-30	alcama	NCBITaxon:7955	Danio rerio	OrthoInspector|OMA|OrthoFinder|Hieranoid|SonicParanoid|ZFIN|PhylomeDB|PANTHER|InParanoid	9	10	No	Yes
# HGNC:2681	DAXX	NCBITaxon:9606	Homo sapiens	ZFIN:ZDB-GENE-010110-3	daxx	NCBITaxon:7955	Danio rerio	OrthoInspector|OMA|OrthoFinder|Hieranoid|SonicParanoid|ZFIN|PhylomeDB|PANTHER|Ensembl Compara|InParanoid	10	10	Yes	Yes



#######
# output of TOBIAS BINDdetect

# directories named <TF>_<motif.version>
# file <TF>_<motif.version>_overview.txt

#TFBS_chr	TFBS_start	TFBS_end	TFBS_name	TFBS_score	TFBS_strand	peak_chr	peak_start	peak_end	additional_1	additional_2	additional_3	additional_4	additional_5	LEC_merged_replicates_footprints_score	BEC_merged_replicates_footprints_score	LEC_merged_replicates_footprints_bound	BEC_merged_replicates_footprints_bound	LEC_merged_replicates_footprints_BEC_merged_replicates_footprints_log2fc
#chr1	18810	18820	Ascl2_MA0816.1	7.75194	-	chr1	18610	19386	macs2_broad_atac_peak_5	.	.	ENSDARG00000102097	nfkbiz	0.85775	0.71893	0	0	0.21452
#chr1	18956	18966	Ascl2_MA0816.1	8.71182	+	chr1	18610	19386	macs2_broad_atac_peak_5	.	.	ENSDARG00000102097	nfkbiz	0.35128	0.55749	0	0	-0.49965

# all fields
#TFBS_chr	TFBS_start	TFBS_end	TFBS_name	TFBS_score	TFBS_strand	peak_chr	peak_start	peak_end	additional_1	additional_2	additional_3	additional_4	additional_5	LEC_merged_replicates_footprints_score	BEC_merged_replicates_footprints_score	LEC_merged_replicates_footprints_bound	BEC_merged_replicates_footprints_bound	LEC_merged_replicates_footprints_BEC_merged_replicates_footprints_log2fc

# fields to consider
#additional_5	LEC_merged_replicates_footprints_score	BEC_merged_replicates_footprints_score	LEC_merged_replicates_footprints_bound	BEC_merged_replicates_footprints_bound	LEC_merged_replicates_footprints_BEC_merged_replicates_footprints_log2fc


########
# annoated peaks - to get peaks in promoters

# ATAC_peaks_annot_closest.ensembl.consensus_replicate_merged_MACS2_broad_TMM_paired.14xi2023.tsv

# logFC	logCPM	LR	PValue	FDR	Chr	Start	End	PeakID	seqnames	start	end	width	strand	gc	annotation.closest	geneChr.closest	geneStart.closest	geneEnd.closest	geneLength.closest	geneStrand.closest	geneId.closest	transcriptId.closest	distanceToTSS.closest	ENTREZID.closest	SYMBOL.closest	GENENAME.closest
# -2.36637693349842	4.57286806187301	176.602720586521	2.67437659610388e-40	2.31844381492841e-35	chr7	33732998	33733703	macs2_broad_atac_peak_78922	chr7	33732998	33733703	706	*	0.424929178470255	Distal Intergenic	937	33645652	33684632	38981	2	ENSDARG00000069006	ENSDART00000130553	-48366	100007463	tle3b	TLE family member 3, transcriptional corepressor b
# -1.91985769668591	5.72878314105868	172.08120950317	2.59791748446956e-39	1.12608032323075e-34	chr16	27306094	27307036	macs2_broad_atac_peak_26007	chr16	27306094	27307036	943	*	0.422057264050901	Distal Intergenic	946	27345383	27369338	23956ENSDARG00000055854	ENSDART00000078250	-38347	548604	nr4a3	nuclear receptor subfamily 4, group A, member 3
# -1.9607082070662	5.2375031639393	159.460811404984	1.484074679771e-36	3.57921263684291e-32	chr9	10747618	10749089	macs2_broad_atac_peak_86094	chr9	10747618	10749089	1472	*	0.385190217391304	Exon (ENSDART00000102472/ENSDARG00000070043, exon 11 of 18)	939	10758339	10769097	10759	2	ENSDARG00000109373	ENSDART00000189967	20008	566927	LOC566927	protein lifeguard 3-like
# 1.85315526852375	5.85811231147094	159.248366816488	1.65148060898728e-36	3.57921263684291e-32	chr17	32936210	32937306	macs2_broad_atac_peak_30449	chr17	32936210	32937306	1097	*	0.469462169553327	Distal Intergenic	947	32822588	32865788	43201	2	ENSDARG00000055158	ENSDART00000077476	-70422	30679	prox1a	prospero homeobox 1a
# -1.93862876658958	5.25491515574702	157.961756256377	3.1549856287221e-36	5.47017718279095e-32	chr5	3169503	3170944	macs2_broad_atac_peak_68111	chr5	3169503	3170944	1442	*	0.395284327323162	Promoter (<=1kb)	935	3170736	3170826	91	1	ENSDARG00000082167	ENSDART00000118394	0	100033664	mir142a	microRNA 142a

# fields to consider
# 9 PeakID
# 16 annotation.closest

# uses match to "Promoter" which gives the following annotations
# Promoter (<=1kb)
# Promoter (1-2kb)
# Promoter (2-3kb)



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Path qw(make_path);
use List::Util qw( min max );
use List::MoreUtils qw(any); 

my $script_name="BINDdetect2nodes";

## set this
my $score_FC_cutoff=0;



if ($ARGV[0]eq qw '-h'){
	print "this script lists TF - target interactions based on TF footprinting results from TOBIAS \n";
	print "as well as various annotations\n";
	print "if in doubt, please check this script header for file format examples\n";
	print "please note this is a PROJECT SPECIFIC TOOL\n\n";

	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--BINDdetect_dir: /path/to/dir/BINDdetect_results\n";
	print "--infile_ortho: /path/to/file/Drer_hsa.orthology-alliance.tsv from Genome Alliance orthology resource\n";
	print "--infile_motifs: /path/to/file/all_motifs.txt used for footprinting (from Jaspar)\n";
	print "--annot_peaks: /path/to/file/annotated_peaks\n";

	print "--outdir: /path/to/outdir (will be created if does not exist)\n";
	print "--pref: file name prefix\n";
	print "-h prints this message\n";
}


else{
	my $parameters=join(' ', @ARGV);

	GetOptions(
		'infile_ortho=s'		=>	\(my $infile_orthol),
		'infile_motifs=s'		=>	\(my $infile_tf_motifs),
		'annotated_peaks=s'		=>	\(my $infile_annotpeaks),
		'indir_BINDdetect=s'		=>	\(my $indir_BINDdetect),
		'pref=s'		=>	\(my $prefix),
		'outdir=s'		=>	\(my $outdir)	
	) or die "Error in command line arguments";

	print "Input files: orthology info $infile_orthol\nmotifs $infile_tf_motifs\nannotated peaks $infile_annotpeaks\nBINDdetect $indir_BINDdetect\n";


	# other variables
	my $score_FC_cutoff=0; #cutoff on footprint score fild change in LEC vs BEC

	############
	#### outfiles
	eval { make_path($outdir); 1};
	if ($@) {
  		print "Couldn't create $outdir: $@";
	}

	my $outfile=$outdir."/LEC_BEC_TF_".$score_FC_cutoff.".".$prefix.".footprints_to_nodes.tsv";

	#############
	### get annotations

	#### gene orthology
	my %gene_names_hsa_drer_array;
	open (INFILE_GENES, "<","$infile_orthol") or die "Cannot open input file $infile_orthol: $!"; 
	while(<INFILE_GENES>){

		chomp $_;
		my @line=split /\t/;

		unless ($_=~m/^#/){
			if ( ($line[2] eq qw "NCBITaxon:9606") && ($line[6] eq qw "NCBITaxon:7955") ){
				my $hsa_gene=uc($line[1]);
				push @{$gene_names_hsa_drer_array{$hsa_gene}}, $line[5]; #because >1 drer orthologues to some hsa genes
			}
		}
	}

	close(INFILE_GENES);


	#### peaks in promoters
	my %peaks_promoters;
	open (INFILE_PEAKS, "<","$infile_annotpeaks") or die "Cannot open input file $infile_annotpeaks: $!"; 
	while(<INFILE_PEAKS>){

		chomp $_;
		my @line=split /\t/;

		unless ($_=~m/^logFC/){
			if ($line[15] =~m/Promoter/){
				$peaks_promoters{$line[8]}=$line[15];
				#print "$line[8]\t$line[15]\n";
			}
		}


	}
	close(INFILE_PEAKS);


	################
	### motifs used in the footprinting
	# while parsing this information, count TFs in Jaspar which have and do not have annotated orthologue, lost those that do not have


	open (INFILE_MOT, "<","$infile_tf_motifs") or die "Cannot open input file $infile_tf_motifs: $!"; 

	my $counter_all=0;
	my $counter_orthol=0;
	my $counter_no_orthol=0;
	my $counter_fam=0;

	my %TF_drer_hsa;
	my %TF_hsa_drer;
	my %TF_hsa_motifid;


	while(<INFILE_MOT>){
		
		chomp $_;

		# motif input in jaspar format
		if ($_ =~m/^>/){

			++$counter_all;

			my @line=split /\t/;

			# for motif line >MA0737.1	GLIS3
			my $Tf=$line[1]; # Original capitalisation
			my $TF=uc($Tf);
			my $motif=$line[0];
			$motif=~s/>//;

			if( exists $gene_names_hsa_drer_array{$TF}){
				++$counter_orthol;

				foreach my $drer_TF_paralog (@{$gene_names_hsa_drer_array{$TF}}){

					# create data structures with only TFs (also for drer arrary)
					#because >1 hsa orthologues to some drer genes is possible
					push @{$TF_drer_hsa{$drer_TF_paralog}},$TF;
					push @{$TF_hsa_drer{$TF}},$drer_TF_paralog;

					$TF_hsa_motifid{$Tf}=$motif; #keep the original case for traversing dir structure

				}
			}else{
				if($TF =~m/::/){ #motif family / dimer rather than individual TF e.g. ERF::FOXI1
					++$counter_fam;
				}
				else{
					++$counter_no_orthol;
				}
			}

		}


		# motif input in meme format
		if ($_=~m/^MOTIF/){

			++$counter_all;

			my @line=split /\s{1,}/;

			#motif name line MOTIF MA0117.1 Mafb_117.1
			my $Tf=$line[2]; # we use the original TF notation (some are hsa some are mmu, hence differences in capitalisation)
			my $motif=$line[1];
			my $TF=uc($Tf);

			if( exists $gene_names_hsa_drer_array{$TF}){
				++$counter_orthol;

				foreach my $drer_TF_paralog (@{$gene_names_hsa_drer_array{$TF}}){

					# create data structures with only TFs (also for drer arrary)
					#because >1 hsa orthologues to some drer genes is possible
					push @{$TF_drer_hsa{$drer_TF_paralog}},$TF;
					push @{$TF_hsa_drer{$TF}},$drer_TF_paralog;

					#was
					$TF_hsa_motifid{$Tf}=$motif; #keep the original case for traversing dir structure

				}
			}else{
				if($TF =~m/::/){ #motif family / dimer rather than individual TF e.g. ERF::FOXI1
					++$counter_fam;
				}
				else{
					++$counter_no_orthol;
				}
			}

		}

	}
	close(INFILE_MOT);

	print "Total TFs: $counter_all\nTFs with orthologues: $counter_orthol\nTFs without orthologues: $counter_no_orthol\nTF / Motif families: $counter_fam\n";


	########################################
	### save the footprints
	##### in promoters

	my $inter_type="TF_binding";

	open (OUTFILE, ">","$outfile") or die "Cannot open input file $outfile: $!"; 

	my $header="Source\tTarget\tType\tScore_LEC\tScore_BEC\tlog2FC_score";
	print OUTFILE "$header\n";



	foreach my $Tf_hsa (keys(%TF_hsa_motifid)) {

		my %target_highestlFC; # storing TF - target highest scores
		#i.e. for given "source" $TF_drer_source, the hash contains keys "targets" $target_TF with values "largest abs(lFC)"


		my $TF_motif=$Tf_hsa . "_" . $TF_hsa_motifid{$Tf_hsa};

		my $overview_pth=$indir_BINDdetect . "/" . $TF_motif . "/" . $TF_motif . "_overview.txt";
		

		if (-e $overview_pth) {

			# parsing round 1
			# parse the txt file and save info only targets which are TFs
			# for multiple instances of the same target select one with largest log2FC(score)

			my %target_score_lFC; #HoAs; keys are targets and array of scores - to select one largest per pair

			open (FH, "<", "$overview_pth") or die "Cannot open overview file $overview_pth: $!";
			while(<FH>){
				chomp $_;
				my @infile_line=split/\t/; #TFBS_chr	TFBS_start	TFBS_end	TFBS_name	TFBS_score	TFBS_strand	peak_chr	peak_start	peak_end	additional_1	additional_2	additional_3	additional_4	additional_5	LEC_merged_replicates_footprints_score	BEC_merged_replicates_footprints_score	LEC_merged_replicates_footprints_bound	BEC_merged_replicates_footprints_bound	LEC_merged_replicates_footprints_BEC_merged_replicates_footprints_log2fc

				my $annot_gene_drer=$infile_line[-6];
				my $annot_gene_score_log2FC=$infile_line[-1];

				my $bound_LEC=$infile_line[-3];
				my $bound_BEC=$infile_line[-2];

				my $peak_id=$infile_line[9];

				if (exists $TF_drer_hsa{$annot_gene_drer}){ #if target is a TF

					if (exists $peaks_promoters{$peak_id}){
						## proceed if bound in any of the conditions
						if( ($bound_LEC==1) | ($bound_BEC==1) ){
							push @{$target_score_lFC{$annot_gene_drer}}, $annot_gene_score_log2FC;
						}
					}
				}
			}
			close(FH);


			# my %target_highestlFC; # storing TF - target highest scores 
			# i.e. for given "source" $Tf (hsa /mmu in original capitalisation) (from the outer loop of hash with motif ids);
			# this hash contains keys "targets" $target_TF (drer) with values "largest abs(lFC)"
			# the $Tf is converetd to $TF_drer_source downstream

			for my $target_TF (keys %target_score_lFC){

				my @pairlog2FC=@{$target_score_lFC{$target_TF}};
				#print "$Tf_hsa\t$target_TF\t@pairlog2FC\n";

				my $min_lfc=min(@pairlog2FC);
				my $max_lfc=max(@pairlog2FC);

				if (abs($min_lfc) >$max_lfc ){
					#print "$min_lfc\n";
					$target_highestlFC{$target_TF}=$min_lfc;
					
				}else{
					#print "$max_lfc\n";
					$target_highestlFC{$target_TF}=$max_lfc;

				}

			}


			my $TF_hsa=uc($Tf_hsa);
			my $TF_drer_source; #drer orthologue(s) of the TF used for footprinting


			# parsing round 2
			# to save selected interaction to file

			open (FH, "<", "$overview_pth") or die "Cannot open overview file $overview_pth: $!";
			while(<FH>){
				chomp $_;
				my @infile_line=split/\t/; #TFBS_chr	TFBS_start	TFBS_end	TFBS_name	TFBS_score	TFBS_strand	peak_chr	peak_start	peak_end	additional_1	additional_2	additional_3	additional_4	additional_5	LEC_merged_replicates_footprints_score	BEC_merged_replicates_footprints_score	LEC_merged_replicates_footprints_bound	BEC_merged_replicates_footprints_bound	LEC_merged_replicates_footprints_BEC_merged_replicates_footprints_log2fc
				my $TF_hsa_motif_id=$infile_line[3];

				my $annot_gene_drer=$infile_line[-6];
				my $annot_gene_score_log2FC=$infile_line[-1];

				my $bound_LEC=$infile_line[-3];
				my $bound_BEC=$infile_line[-2];

				my $score_LEC=$infile_line[-5];
				my $score_BEC=$infile_line[-4];

				my $peak_id=$infile_line[9];


				if (exists $TF_drer_hsa{$annot_gene_drer}){ #if target is a TF

					if (exists $peaks_promoters{$peak_id}){# if the peak is in the promoter

						## proceed if bound in any of the conditions
						if( ($bound_LEC==1) | ($bound_BEC==1) ){
							
							## check if score differencs is the highest (selected TFBS!)
							if ($annot_gene_score_log2FC == $target_highestlFC{$annot_gene_drer}){

								## assign ns to score if not bound in either condition
								if($bound_LEC==0){
									$score_LEC="ns";
								}					

								if($bound_BEC==0){
									$score_BEC="ns";
								}


								### resolve hsa-drer ortholgy
								my @TFs_drer=@{$TF_hsa_drer{$TF_hsa}}; #drer orthologue(s) of TF footprinted

								foreach my $TF_drer (@TFs_drer){

									my $line="$TF_drer\t$annot_gene_drer\t$inter_type\t$score_LEC\t$score_BEC\t$annot_gene_score_log2FC";

									### print line to file

									## if bound in one condition only just print
									if( ($bound_LEC==1 && $bound_BEC==0) | ($bound_LEC==0 && $bound_BEC==1) ){
										print OUTFILE "$line\n";

									}
								
									## if bound in both conditions check if the score ratio > cutoff
									elsif($bound_LEC==1 && $bound_BEC==1){

										if($annot_gene_score_log2FC>=$score_FC_cutoff){
											print OUTFILE "$line\n";
										}
									}
								}
							}
						}
					}
				}
			}
			close(FH);
		}
	}
	close(OUTFILE);
}
exit;
