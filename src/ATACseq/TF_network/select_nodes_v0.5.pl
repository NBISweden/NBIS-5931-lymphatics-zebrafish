#!c:/perl/bin/perl.exe

# script to parse output of BINDdetect2nodes.v0.x.pl (processed output of TOBIAS BINDdetect)
# select subnetwork of a selected TF
# incorporate expression data from scRNA-seq: average expression in LEC and VEC and log2FC VEC vs LEC


# in v 0.5
# prefix in command line
# RNA-seq related cutoffs exposed at the command line


# input files

#######
# head LEC_BEC_TF_lfc0_footprints_to_nodes.tsv                   
# Source	Target	Type	Score_LEC	Score_BEC	log2FC_score
# gmeb1	bhlhe40	TF_binding	2.66119	1.75145	0.56483
# gmeb1	hlfb	TF_binding	ns	1.66878	-1.15752
# gmeb1	nkx6.2	TF_binding	ns	1.03238	-0.45036
# gmeb1	lbx2	TF_binding	1.15742	ns	0.66231
# gmeb1	bach2a	TF_binding	ns	0.96179	-0.69509
# gmeb1	prox1a	TF_binding	3.4133	ns	2.24959




######
#### RNA-seq
## LEC is cluster 0
## VEC is cluster 1


#fold_change_5days.csv
#LEC vs VEC
# clust 0 vs clust 1

# "","avg_log2FC","pct.1","pct.2"
# "slc35a5",1.26146088841453,0.102,0.058
# "ccdc80",0.747049574946611,0.867,0.868
# "nrf1",-0.0534979700418461,0.389,0.416

###
# average_expression_5days.csv
# linear scale
# LEC - cluster 0 VEC - cluster 1

#"","X0","X1","X2","X3","X4","X5","X6","X7","X8"
#"prox1a",132.343180616607,5.12183114182662,4.05317773409395,2.25705428981188,1.59983068169138,3.8915492477986,1.29293139782549,32.7912068732523,649.035753444582



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Path qw(make_path);
use File::Basename qw(basename);



my $script_name="select_nodes";

## set this for results suffix - related to TFscoreLFC used to pull source - target edges file
my $suffix="TFscoreLFC0";

my $level="2"; #for file name only; atm levels hardcoded

my $interaction_type="TF_binding"; #for output files

my $abund_pseudocount=0.1; #avoid log(0)



if ($ARGV[0]eq qw '-h'){
	print "this script selects nodes in a TF interaction network using parsed footprinting and scRNAseq results\n";
	print "and other annotations\n";

	print "if in doubt, please check this script header for file format examples\n";
	print "please note this is a PROJECT SPECIFIC TOOL\n\n";

	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	
	print "--infile_edges: /path/to/file/LEC_BEC_TF_0.footprints_to_nodes.tsv\n";
	print "--infile_log2FC: /path/to/file/240202_fold_change_5days.csv\n";
	print "--infile_abundance: /path/to/file/average_expression_5days.csv\n";
	print "--infile_genelist: /path/to/file/genelist for genes of interest e.g. GO terms associated\n";

	print "--TF: name of TF starting pint e.g. prox1a nr2f2 znf423\n";
	print "--smpl: LEC or BEC\n";

	print "--outdir: /path/to/outdir (will be created if does not exist)\n";
	print "--pref: file name prefix\n";

	print "--rnaseq_lfc: log2FC cutoff for RNAseq data\n";
	print "--rnaseq_frac: fraction cells expressed cutoff for RNAseq data\n";

	print "-h prints this message\n";
}


else{
	my $parameters=join(' ', @ARGV);

	GetOptions(
		'infile_edges=s'		=>	\(my $infile_edges),
		'infile_abundance=s'		=>	\(my $infile_abundance),
		'infile_log2FC=s'		=>	\(my $infile_log2FC),
		'infile_genelist=s'		=>	\(my $infile_genelist),
		'TF=s'		=>	\(my $TF),
		'smpl=s'		=>	\(my $smpl),
		'pref=s'		=>	\(my $prefix),
		'rnaseq_lfc=s'		=>	\(my $rna_log2FC_co),
		'rnaseq_frac=s'		=>	\(my $perc_expr),
		'outdir=s'		=>	\(my $outdir)	
	) or die "Error in command line arguments";


	########################
	### formatting subroutines

	sub format_line # pass array of strings, return tab-delimited line
	{
		my $line_TF=join "\t", @_;
		return($line_TF);
	}


	#log in perl is ln
	sub log2
	{ 
	    my $n = shift; 
	    # using pre-defined log function 
	    return log($n) / log(2); 
	} 




	############
	#### outfiles
	eval { make_path($outdir); 1};
	if ($@) {
  		print "Couldn't create $outdir: $@";
	}

	my $outfile_nt;
	my $outfile_abund;
	my $outfile_pref=$TF . "_" . $smpl.".".$prefix.".RNAlfc_".$rna_log2FC_co."_RNAfracexpr_".$perc_expr."_lvl".$level."_".$suffix;

	if(defined $infile_genelist){
		my $genelist_fname=basename($infile_genelist);
		my $genelist_basename=basename($genelist_fname,".txt");

		$outfile_nt=$outdir."/".$outfile_pref.".".$genelist_basename."_network.tsv";
		$outfile_abund=$outdir."/".$outfile_pref.".".$genelist_basename."_abundance.tsv";

	}else{
		$outfile_nt=$outdir."/".$outfile_pref."_network.tsv";
		$outfile_abund=$outdir."/".$outfile_pref."_abundance.tsv";

	}


	# file with abundance of all TFs in the footprinting analysis - for cytoscape when working with more subnetworks
	# not really that useful as it needs to be imported for each network anyway
	my $outfile_abund_all="all_gene_abundance.tsv";


	##################
	##################
	#### code



	##################
	#### read in genes from gene list

	my %genelist;
	if(defined $infile_genelist){
		open(GENELIST, "<", $infile_genelist) or die "Cannot open file $infile_genelist: $!";
		while(<GENELIST>){
			chomp $_;
			my @line=split/\t/;
			$genelist{$line[0]}="gene";
		}
		close(GENELIST);

	}

	##############
	### read in edges

	my %pair_TFBS_score_smpl; # source-target > score
	my %pair_TFBS_lfc; # source-target > lfc LEC vs BEC

	open (INFILE_EDGES, "<","$infile_edges") or die "Cannot open input file $infile_edges: $!"; 

	while(<INFILE_EDGES>){
		chomp $_;

		unless($_=~m/^Source/){
			my @line=split/\t/;

			my $source_target=$line[0] . "::" . $line[1];

			if ($smpl eq "LEC"){
				$pair_TFBS_score_smpl{$source_target}=$line[3];
			}


			if ($smpl eq "BEC"){
				$pair_TFBS_score_smpl{$source_target}=$line[4];
			}

		$pair_TFBS_lfc{$source_target}=$line[5];

		}
	}
	close(INFILE_EDGES);



	##############
	### read in expression data (linear scale avg expression per cluster)
	### this is not used & will be removed in the final version

	my %gene_abundance_clust0; #abundance on log2 scale  #VEC
	my %gene_abundance_clust1; #abundance on log2 scale  #LEC


	open(INFILE_ABUNDANCE, "<", "$infile_abundance") or die "Cannot open input file $infile_abundance: $!"; 

	while(<INFILE_ABUNDANCE>){
		chomp $_;

		unless($_=~m/""/){
			my @line=split/,/;
			my $gene_name=$line[0];
			$gene_name=~s/"//g;
			$line[1]=$line[1]+$abund_pseudocount;
			$line[2]=$line[2]+$abund_pseudocount;

			my $log2_clust0=log2($line[1]);
			my $log2_clust1=log2($line[2]);

			$gene_abundance_clust0{$gene_name}=$log2_clust0;
			$gene_abundance_clust1{$gene_name}=$log2_clust1;



		}
	}

	close(INFILE_ABUNDANCE);




	##############
	### read in data with log2FC and fraction of cells where the gene is expressed
	### select genes to consider for network based on fraction cells expressed & other criteria

	open(OUTFILE_ABN_ALL, ">", "$outfile_abund_all") or die "Cannot open output file $outfile_abund_all: $!";
	my $header_abn_all="Gene\tlog2FC\tlogcnts_clust0\tlogcnts_clust1";
	print OUTFILE_ABN_ALL "$header_abn_all\n";


	my %gene_lFC;
	my %gene_fraction_expressed;
	my %gene_abundance_line;

	open(INFILE_LFC, "<", "$infile_log2FC") or die "Cannot open input file $infile_log2FC: $!"; 

	while(<INFILE_LFC>){
		chomp $_;

		unless($_=~m/^""/){
			my @line=split/,/;

			my $gene_name=$line[0];
			$gene_name=~s/"//g;
		
			my $log2FC=$line[1];
			my $frac_expr_LEC_BEC=$line[2]."::".$line[3]; ## frac_LEC::frac_BEC

			$gene_lFC{$gene_name}=$log2FC;
			$gene_fraction_expressed{$gene_name}=$frac_expr_LEC_BEC;

			#for printing abundance to a separate file
			my $goi_abund_clust1;
			if(exists $gene_abundance_clust1{$gene_name}){
				$goi_abund_clust1=$gene_abundance_clust1{$gene_name};
			}else{$goi_abund_clust1="NA"}
			my $goi_abund_clust0;
			if(exists $gene_abundance_clust0{$gene_name}){
				$goi_abund_clust0=$gene_abundance_clust0{$gene_name};
			}else{$goi_abund_clust0="NA"}


			#line to print in all TF abundance data in one file
			my $line_abundance=format_line($gene_name,$log2FC,$goi_abund_clust0,$goi_abund_clust1);
			print OUTFILE_ABN_ALL "$line_abundance\n";

			$gene_abundance_line{$gene_name}=$line_abundance;


		}

	}
	close(INFILE_LFC);

	close(OUTFILE_ABN_ALL);




	##############
	### select genes fulfilling criteria: fraction expressed, presence in a gene list, absolute log2(fold change) of $rna_log2FC_co

	my %genes_nt_selected;

	while ( my ($gene_id, $fraction_expr) = each (%gene_fraction_expressed)) {

		my($frac_expr_LEC,$frac_expr_BEC)=split(/::/,$fraction_expr);

		if( ($frac_expr_LEC>=$perc_expr) | ($frac_expr_BEC>=$perc_expr) ){

			my $TF_lFC;
			if(exists $gene_lFC{$gene_id}){
				$TF_lFC=$gene_lFC{$gene_id};
			}else{$TF_lFC="NA"}

			if($TF_lFC ne qw "NA") {
				if (abs($TF_lFC) >= $rna_log2FC_co){#only consider TFs with absolute log2(fold change) of $rna_log2FC_co
					if(defined $infile_genelist){
						if (exists $genelist{$gene_id}){
							$genes_nt_selected{$gene_id}="selected";
						}
					}else{# no gene list defined so all passing genes are included
						$genes_nt_selected{$gene_id}="selected";

					}
				}
			}
		}
	}



	##############
	##############
	### construct subnetwork

	my %genes_in_nt; # for propagating down the network
	my %pairs_in_nt; # for output file

	##############
	### define functions

	### subroutine to get abundance and log2FC from RNAseq data for given TF

	sub get_exprs
	{
		my $goi = shift; 

		my $TF_abund_clust0;
		if(exists $gene_abundance_clust0{$goi}){
			$TF_abund_clust0=$gene_abundance_clust0{$goi};
		}else{$TF_abund_clust0="NA"}
		
		my $TF_abund_clust1;
		if(exists $gene_abundance_clust1{$goi}){
			$TF_abund_clust1=$gene_abundance_clust1{$goi};
		}else{$TF_abund_clust1="NA"}
		
		my $TF_lFC;
		if(exists $gene_lFC{$goi}){
			$TF_lFC=$gene_lFC{$goi};
		}else{$TF_lFC="NA"}

		my @results=($TF_lFC,$TF_abund_clust0,$TF_abund_clust1);

		return(@results);
	}


	#### subroutine to get a list of network subnodes (targets) to be passed further for final file printing
	#### the list is transformed to a hash because the keys are unique, to avoid gene duplicates
	#### the subroutine can be repeated as many times as disired; new subnodes will be added to the hash upon each iteration

	sub get_subnodes
	{
		my $gene_node = shift; 


		my @genes_in_nt;
		my @pairs_in_nt;

		if (exists $genes_nt_selected{$gene_node}){ #gene fulfilled criteria to be included in the network
		
			foreach my $interacting_pair (keys(%pair_TFBS_score_smpl)){
				my ($source_TF, $target_TF)=split/::/,$interacting_pair;

				if($source_TF eq $gene_node){

					my $score_TFBS=$pair_TFBS_score_smpl{$interacting_pair};

					if($score_TFBS ne "ns"){

						if (exists $genes_nt_selected{$target_TF}){ #gene fulfilled criteria to be included in the network

							push @genes_in_nt, $target_TF;
							push @pairs_in_nt, $interacting_pair;
						}

					}

				}
			}
		}
		my @res=(\@genes_in_nt,\@pairs_in_nt);

		return(@res);
	}



	##############
	### targets of the initial TF (1st level sources)
	
	$genes_in_nt{$TF}="source";

	my @node_results=get_subnodes($TF);

	my $ref1=$node_results[0];
	foreach my $nt_node (@$ref1){
		#print "$nt_node\n";
		$genes_in_nt{$nt_node}="target";
	}

	my $ref2=$node_results[1];
	foreach my $nt_pair (@$ref2){
		#print "$nt_pair\n";
		$pairs_in_nt{$nt_pair}="pair";
	}


	##############
	### targets of the primary targets (2nd level sources)

	foreach my $TF_subnode (keys %genes_in_nt){
		print "\n$TF_subnode\n";
		my @node_results_lvl=get_subnodes($TF_subnode);
		
		my $ref1_lvl=$node_results_lvl[0];
		foreach my $nt_node (@$ref1_lvl){
		print "$nt_node\n";
			$genes_in_nt{$nt_node}="target";
		}
		
		my $ref2_lvl=$node_results_lvl[1];
		foreach my $nt_pair (@$ref2_lvl){
			$pairs_in_nt{$nt_pair}="pair";
		}
	}



	foreach my $TF_pair (keys %pairs_in_nt){
		print "$TF_pair\t$pairs_in_nt{$TF_pair}\n";
	}

	##############
	### output network files
	### outfiles with subnetwork and its abundance information


	open(OUTFILE_NT, ">", "$outfile_nt") or die "Cannot open output file $outfile_nt: $!"; 
	my $header="Source\tTarget\tType\tScore"."_".$smpl;
	print OUTFILE_NT "$header\n";


	open(OUTFILE_ABN, ">", "$outfile_abund") or die "Cannot open output file $outfile_abund: $!"; 
	my $header_abn="Gene\tlog2FC\tlogcnts_clust0\tlogcnts_clust1";
	print OUTFILE_ABN "$header_abn\n";

	#line for main TF at the root of the nt
	my $TF1_abund_clust0=$gene_abundance_clust0{$TF};
	my $TF1_abund_clust1=$gene_abundance_clust1{$TF};
	my $TF1_lFC=$gene_lFC{$TF};

	my $line_TF1=$TF . "\t" . $TF1_lFC . "\t" . $TF1_abund_clust0 . "\t" . $TF1_abund_clust1;
	print OUTFILE_ABN "$line_TF1\n";


	foreach my $interacting_pair (keys(%pairs_in_nt)){
		my ($source_TF, $target_TF)=split/::/,$interacting_pair;
		my $score_TFBS=$pair_TFBS_score_smpl{$interacting_pair};
		my $lfc_score_TFBS=$pair_TFBS_lfc{$interacting_pair};
		my ($TF_lFC,$TF_abund_clust0,$TF_abund_clust1)=get_exprs($target_TF);
		my $line_TF=format_line($target_TF,$TF_lFC,$TF_abund_clust0,$TF_abund_clust1);
		print OUTFILE_ABN "$line_TF\n";
		my $line_nt=format_line($source_TF,$target_TF,$interaction_type,$score_TFBS);
		print OUTFILE_NT "$line_nt\n";
	}

	close(OUTFILE_ABN);
	close(OUTFILE_NT);


}
exit;