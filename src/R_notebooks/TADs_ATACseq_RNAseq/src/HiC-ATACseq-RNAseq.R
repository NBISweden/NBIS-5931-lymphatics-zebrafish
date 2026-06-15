



## ---- prep-env



library(renv)

renv::activate()
##renv::restore()

library(tidyverse)
library(dplyr)
library(reshape2)
library(kableExtra)

library(scales)
library(GenomicRanges)
library(regioneR)
library(regioneReloaded)
library(ChIPseeker)
#library(orthogene)

library(AnnotationHub)
library(org.Dr.eg.db)
library(BSgenome.Drerio.UCSC.danRer11)

library(viridis)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(eulerr)

#library(rtracklayer)
#library(HiContacts)


##renv::snapshot()


tab_n_sv=0
tab_n=0 # table labels
fig_n=0
fig_n_cap=0

# for file names
file_pref="9vi2026"

# colour scheme

LEC_light="#10963085"
LEC_dark="#109630"
BEC_light="#762A8385"
BEC_dark="#762A83"



options(scipen=999)



#matrix bin size
binsize=10000


## ---- feature-TAD-overlap-cutoffs

# how many bins are tolerated to be called the same boundary and for feature overlap
# bins to consider TADs to be identical in LEC and BEC
ontad_wiggle_bins=3


## cutoffs for TAD-feature positional overlaps

## in kb
#loop_wiggle_distance=10
atacpeak_wiggle_distance=10

atac_FDR_CO=0.05
#atac_FDR_CO=0.01
#atac_lFC_CO=0
atac_lFC_CO=0.5
atac_lCPM_CO=2

#n permutations for overlap test
n_perms=500

#n permutations for producion run
#n_perms=5000



## ---- functions


#save tables in tab-delimted format
save.table=function(df=df, file=file, dir=dir) {
	write.table(df, file.path(dir, file), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, fileEncoding = "")
}


#save tables in tab-delimted format
save.bed=function(df=df, file=file, dir=dir) {
	write.table(df, file.path(dir, file), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, fileEncoding = "")
}




# save plots in pdf and eps
save.plots=function(plot_ts=plot_ts,savedir=savedir,pref=pref) {
        fname_pdf=paste(pref,"pdf",sep=".")
        fname_eps=paste(pref,"eps",sep=".")

        plot_outdir=file.path(savedir,pref)
        dir.create(plot_outdir, recursive=TRUE)

        graphics.off()
        pdf(file.path(plot_outdir,fname_pdf))
                print(plot_ts)
        dev.off()

        graphics.off()
        setEPS()                                       
        postscript(file.path(plot_outdir,fname_eps))             
                print(plot_ts)
        dev.off()   

}



## ---- dirs
projrootdir="/Users/agatasmialowska/NBIS/NBIS_proj/5931_fishlymphatics/2026/manuscript_figures"

datadir=file.path(projrootdir,"data")

resdir=file.path(projrootdir,"results",paste0("HiC_ATACseq_RNAseq_",file_pref))
dir.create(resdir, recursive=TRUE)



## ---- gene-models

annotdir=file.path(projrootdir,"reference","annotation")

ensembl_annot_trxdb=file.path(annotdir,"Ensembl.txdb.Drer11.Rdata")
txdb_ens=loadDb(ensembl_annot_trxdb)


# change IDs of the flanking genes
# this is a bit of an roundabout manner, but allows for adding more annotations to flanking genes, if required and keeps everything in df format
# to be used by apply
ah=AnnotationHub(localHub=FALSE)
drer11Db=query(ah, pattern = c("Danio Rerio", "EnsDb"))
drer11.ensDb=drer11Db[["AH116824"]]
drer11.ensDb


# gene names from biomart
fname="drer11_gene_names.tab"
gene_annot_drer11=read.delim(file.path(annotdir,fname), sep="\t", header=TRUE, quote = "")

# get coordinates for all genes; if external_gene_name not available, sub ensembl_gene_id
genes_all=gene_annot_drer11|>
        dplyr::select(ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position,strand)|>
        dplyr::distinct(ensembl_gene_id, external_gene_name, .keep_all = TRUE) |>
        dplyr::filter(!str_detect(chromosome_name, "CHR_ALT"))|>        #rm genes on alt contigs
        mutate_if(is.character, list(~na_if(.,"")))|>
        mutate(external_gene_name = coalesce(external_gene_name, ensembl_gene_id))

genes_all$strand=gsub("^1","+1",genes_all$strand,perl=TRUE)


## ---- gene-lists

###################################
## genes - BEC / LEC markers

## based on RNAseq, same cutoffs as for ATAC peak overlaps
rnaseq_log2FC_CO=0.3
rnaseq_padj_CO=0.1
rnaseq_pct_expr_CO=0.3

#Files sent by Marleen (29 Oct 2024):
rnaseq_datadir=file.path(datadir,"RNAseq","RNA_seq_Marleen_29x2024")

DE_d5=file.path(rnaseq_datadir,"fullmarkers_5dpf.csv")

DE_rnaseq_d5=read.csv(DE_d5)

DE_rnaseq_d5.sig=DE_rnaseq_d5|>
        dplyr::filter(p_val_adj<rnaseq_padj_CO & abs(avg_log2FC)>=rnaseq_log2FC_CO)|>
        dplyr::filter(pct.1>=rnaseq_pct_expr_CO | pct.2>=rnaseq_pct_expr_CO)
nrow(DE_rnaseq_d5.sig)


rnaseq_fc.d5.lec=DE_rnaseq_d5.sig|>
        dplyr::filter(avg_log2FC>0 & pct.1>=rnaseq_pct_expr_CO)
nrow(rnaseq_fc.d5.lec)

rnaseq_fc.d5.bec=DE_rnaseq_d5.sig|>
        dplyr::filter(avg_log2FC<0 & pct.2>=rnaseq_pct_expr_CO)
nrow(rnaseq_fc.d5.bec)


genes_lec_dpf5.ls=rnaseq_fc.d5.lec$X
genes_bec_dpf5.ls=rnaseq_fc.d5.bec$X

genes_lec_dpf5.ls.sorted=genes_lec_dpf5.ls[order(genes_lec_dpf5.ls)]
genes_bec_dpf5.ls.sorted=genes_bec_dpf5.ls[order(genes_bec_dpf5.ls)]

# get coordinates for the genes of interest

genes_lec=gene_annot_drer11|>dplyr::filter(external_gene_name%in%genes_lec_dpf5.ls)|>
        dplyr::select(ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position)|>
        dplyr::distinct(ensembl_gene_id, external_gene_name, .keep_all = TRUE)|>
        dplyr::filter(!str_detect(chromosome_name, "CHR_ALT"))#rm genes on alt contigs

genes_bec=gene_annot_drer11|>dplyr::filter(external_gene_name%in%genes_bec_dpf5.ls) |>
        dplyr::select(ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position)|>
        dplyr::distinct(ensembl_gene_id, external_gene_name, .keep_all = TRUE)|>
        dplyr::filter(!str_detect(chromosome_name, "CHR_ALT"))#rm genes on alt contigs


genes_bec_incl=genes_bec$external_gene_name[order(genes_bec$external_gene_name)]
genes_lec_incl=genes_lec$external_gene_name[order(genes_lec$external_gene_name)]


###################################
## vascular markers
## 

## Gene IDs converted to D.rerio previously using combination of orthogene::convert_orthologs and manual conversion
vascular_markers=c("ace","afdna","akt3a","anxa3a","anxa3b","aplnrb","batf3","CABZ01044099.1","CABZ01084131.1","calcrla","calcrlb","cav1","cav2","cavin1a","cavin1b","cavin2a","cavin2b","cd151","cd248a","cd248b","cd36","cd44a","cd44b","cd82a","cd82b","cdh5","cdh6","cldn5a","cldn5b","erg","esama","esamb","ets1","f3a","f3b","flt1","flt4","ghra","ghrb","il11ra","il13ra1","il1rapl1b","il6r","itga4","itga5","ITGB1","itgb1a","itgb1b","itgb1b.1","itgb1b.2","itgb3a","itgb3b","kdr","kdrl","kita","kitb","lama4","lama5","lamb2","ldb2a","ldb2b","lyve1b","mcama","mcamb","mef2ca","mef2cb","mrc1a","myl12.2","notch1a","notch1b","podxl","prox1a","ptk2aa","ptk2ab","pxna","ramp2","RAMP3","robo4","s1pr1","sash1a","sash1b","sele","selp","si:ch1073-391i24.1","si:ch211-156j16.1","si:ch211-286o17.1","si:ch211-8c17.3","si:ch73-22o12.1","si:dkey-237i9.8","sox18","st6galnac3","stab1","stab2","stard13a","stard13b","tek","tie1","tjp1b","vcam1b","vwf","zgc:153911","zmp:0000000936")

# get coordinates for the genes of interest
# gene names from biomart

genes_vasc=gene_annot_drer11|>dplyr::filter(external_gene_name%in%vascular_markers)|>
        dplyr::select(ensembl_gene_id,external_gene_name,chromosome_name,start_position,end_position)|>
        dplyr::distinct(ensembl_gene_id, external_gene_name, .keep_all = TRUE)|>
        dplyr::filter(!str_detect(chromosome_name, "CHR_ALT"))#rm genes on alt contigs


genes_vasc_incl=genes_vasc$external_gene_name[order(genes_vasc$external_gene_name)]





## ---- onTAD-functions


# add info on TAD boundary hierarchy
# for situations of same boundary, multiple levels (when inner TADs share the boundary with outer TADs) assign the outermost level

read_ontad_raw=function(dirpth=dirpth, sample_id=sample_id){
        require(dplyr)

        tads_smpl=list()

        #25 chr
        for (i in c(1:25)){
                fpth_i=file.path(dirpth,paste0("chr_",i,".tad"))
                tad_i=read.csv(fpth_i,sep="\t",header=FALSE)
                tad_i$chr=i
                tads_smpl[[i]]=tad_i
        }
        tads_smpl_df=bind_rows(tads_smpl, .id = "df_label")
        tads_smpl_df=tads_smpl_df[,!names(tads_smpl_df)%in%c("df_label")]
        colnames(tads_smpl_df)=c("startpos","endpos","TAD_level","TADmean","TADscore","chr")

        tads_smpl_df=tads_smpl_df|>
                        dplyr::filter(TAD_level>0)|>
                dplyr::mutate(start=startpos*binsize, end=endpos*binsize)|>
                        dplyr::mutate(tad_id=paste(chr,startpos*binsize,endpos*binsize,sep="_"))|>
                dplyr::mutate(bnd_id_5=paste(chr,start,sep="_"))|>
                dplyr::mutate(bnd_id_3=paste(chr,end,sep="_"))|>
                dplyr::mutate(smpl=sample_id)

        #assign boundary level
        tads_smpl_df=tads_smpl_df|>
                group_by(bnd_id_5)|>
                dplyr::mutate(bnd_5_level=min(TAD_level))|>
                group_by(bnd_id_3)|>
                dplyr::mutate(bnd_3_level=min(TAD_level))

        return(as.data.frame(tads_smpl_df))
}




## ---- TADcompare-functions


## basic formatting of TADcompare results

format_TADcompare=function(TADcompare_res=TADcompare_res){
        require(dplyr)

        colnames(TADcompare_res)[6]=c("Status")
        TADcompare_res=TADcompare_res|>dplyr::mutate(bnd_id=paste(chr,Boundary,sep="_"))

}





## ---- data-ontad

## megamap TADs
# using unprocessed OnTAD output
mega_TAD_root=file.path(datadir,"OnTAD/raw/megamap")
mega_TAD_pref="Lsize_7_penalty_0.075_maxsize_300"

ontad_LEC=file.path(mega_TAD_root,"LEC",mega_TAD_pref)
ontad_BEC=file.path(mega_TAD_root,"BEC",mega_TAD_pref)

tads_LEC=read_ontad_raw(dirpth=ontad_LEC, sample_id="LEC")
tads_BEC=read_ontad_raw(dirpth=ontad_BEC, sample_id="BEC")





## ---- data-TADcompare-megaTADs


# TADcompare using sample matrices and mega TADs

tadcompare_datadir=file.path(datadir,"OnTAD/TADcompare_28iii2025/ontad_samples_consensus_28iii2025")
ontad_pref.m="Lsize_7_penalty_0.075_maxsize_300"

#BEC
test_name="BEC_mega_TADs.101_vs_103"

wins="50"
mega_bnds.101_103.ontad.50 = readRDS(file.path(tadcompare_datadir,paste0("TADcompare.ontad.diff_zCO2_5.",test_name,".",ontad_pref.m,".windows_",wins,".18iii2025.rds")))
mega_bnds_allsites.101_103.ontad.50 = readRDS(file.path(tadcompare_datadir,paste0("TADcompare.ontad.BoundaryScore_allpos_zCO2_5.",test_name,".",ontad_pref.m,".windows_",wins,".18iii2025.rds")))


#LEC
test_name="LEC_mega_TADs.102_vs_104"

wins="50"
mega_bnds.102_104.ontad.50 = readRDS(file.path(tadcompare_datadir,paste0("TADcompare.ontad.diff_zCO2_5.",test_name,".",ontad_pref.m,".windows_",wins,".18iii2025.rds")))
mega_bnds_allsites.102_104.ontad.50 = readRDS(file.path(tadcompare_datadir,paste0("TADcompare.ontad.BoundaryScore_allpos_zCO2_5.",test_name,".",ontad_pref.m,".windows_",wins,".18iii2025.rds")))


mega_bnds.102_104.ontad.50=format_TADcompare(TADcompare_res=mega_bnds.102_104.ontad.50)
mega_bnds.101_103.ontad.50=format_TADcompare(TADcompare_res=mega_bnds.101_103.ontad.50)

mega_bnds_allsites.102_104.ontad.50=format_TADcompare(TADcompare_res=mega_bnds_allsites.102_104.ontad.50)
mega_bnds_allsites.101_103.ontad.50=format_TADcompare(TADcompare_res=mega_bnds_allsites.101_103.ontad.50)



# add info for missing boundaries

get_all_bnds=function(df_tads=df_tads,df_bnds=df_bnds,df_allsites=df_allsites){
        require(dplyr)

        all_bnds_tads=unique(c(df_tads$bnd_id_5,df_tads$bnd_id_3))
        missing_bnds=setdiff(all_bnds_tads,df_bnds$bnd_id)
        missing_bnds_allpos=df_allsites|>dplyr::filter(bnd_id%in%missing_bnds)
        all_bnds=rbind(df_bnds,missing_bnds_allpos)
        return(all_bnds)
}

mega_bnds.102_104.ontad.50.all=get_all_bnds(df_tads=tads_LEC,df_bnds=mega_bnds.102_104.ontad.50,df_allsites=mega_bnds_allsites.102_104.ontad.50)
mega_bnds.101_103.ontad.50.all=get_all_bnds(df_tads=tads_BEC,df_bnds=mega_bnds.101_103.ontad.50,df_allsites=mega_bnds_allsites.101_103.ontad.50)



## ---- data-tadcompare-wins

tadcompare_dir=file.path(datadir,"TADcompare/2i2026/ontad_samples_consensus_v2_2i2026")

# these in form to be used with paste0
pref_smpl="LEC_BEC_mega_TADs.LEC_vs_BEC.windows_"
suf_run=".Lsize_7_penalty_0.075_maxsize_300.2i2026.rds"

pref_allpos="TADcompare.ontad.BoundaryScore_allpos."
pref_bnds="TADcompare.ontad.diff."



append_score=function(bnds=bnds,win=win){
  require(dplyr)

  colname=paste0("Gap_Score.",win)
  bnds=bnds|>dplyr::mutate(colname1=Gap_Score)
  names(bnds)[names(bnds) == "colname1"]=colname

  return(bnds)
}



read_TADcompare=function(tadcompare_dir=tadcompare_dir,win=win){

  fname_allpos=paste0(pref_allpos,pref_smpl,win,suf_run)
  fpath_allpos=file.path(tadcompare_dir,win,fname_allpos)

  allpos.df=readRDS(fpath_allpos)
  allpos.df=format_TADcompare(TADcompare_res=allpos.df)
  allpos.df=append_score(bnds=allpos.df,win=win)

  fname_bnds=paste0(pref_bnds,pref_smpl,win,suf_run)
  fpath_bnds=file.path(tadcompare_dir,win,fname_bnds)

  bnds.df=readRDS(fpath_bnds)
  bnds.df=format_TADcompare(TADcompare_res=bnds.df)
  bnds.df=append_score(bnds=bnds.df,win=win)


  res_out=list()
  res_out[["bnds"]]=bnds.df
  res_out[["allpos"]]=allpos.df

  return(res_out)
}



window_sizes=c(20,35,50)


all_bnd_calls=list()

for (win in window_sizes ){

  all_bnd_calls[[win]]=read_TADcompare(tadcompare_dir=tadcompare_dir,win=win)

}




## ---- published_TADs

#Yang 2020

datadir_yang=file.path(datadir,"published/Yang2020/TADs_Yang2020")

# brain
brain_tads_pth=file.path(datadir_yang,"brain_TADs-Table1.tsv")
brain_tads=read.csv(brain_tads_pth, sep="\t", header=FALSE)

brain_tads=brain_tads[,c(1:3)]
colnames(brain_tads)=c("chr","start","end")
brain_tads$length=brain_tads$end-brain_tads$start+1


muscle_tads_pth=file.path(datadir_yang,"Muscle_TADs-Table1.tsv")
muscle_tads=read.csv(muscle_tads_pth, sep="\t", header=FALSE)

muscle_tads=muscle_tads[,c(1:3)]
colnames(muscle_tads)=c("chr","start","end")
muscle_tads$length=muscle_tads$end-muscle_tads$start+1

#Franke 2021

datadir_franke=file.path(datadir,"published/Franke2021/TADs/domains")

file1=file.path(datadir_franke,"boundaries_24h_wt_clean.domains.bed")
#file2=file.path(datadir_franke,"boundaries_48h_wt_clean.domains.bed")


Franke_1_tads=read.csv(file1, sep="\t", header=FALSE)

tads=Franke_1_tads

tads=tads[,c(1:3)]
colnames(tads)=c("chr","start","end")
tads$length=tads$end-tads$start

Franke_1_tads=tads


## ---- CTCF-data

# data from Franke et al 2021 (GSE156096)



ctcf_chip_datadir=file.path(datadir,"published/Franke2021/CTCF_chip/GSE156096_RAW")

ctcf_files=list.files(path=ctcf_chip_datadir, pattern = "\\.narrowPeak$", full.names = TRUE)

ctcf_chip_data=list()
ctcf_chip_smpl=list()

for (i in 1:length(ctcf_files)){

        path=ctcf_files[[i]]
        fname=basename(path)
        smpl=gsub("_peaks.narrowPeak","",fname,perl=TRUE)
        geo_acc_smpl_id=unlist(strsplit(smpl,"_ChIP_",fixed=TRUE))
        geo_acc=geo_acc_smpl_id[1]
        smpl_id=geo_acc_smpl_id[2]
    df=read.csv(path,header=FALSE,sep="\t")
    df=df|>dplyr::filter(grepl('chr', V1))|>dplyr::filter(!grepl('Un|alt', V1))
    df$V1=gsub("chr","",df$V1)

    ctcf_chip_data[[geo_acc]]=df
    ctcf_chip_smpl[[geo_acc]]=smpl_id
}



ctcf_wt_24_r1=ctcf_chip_data[["GSM4724551"]]
ctcf_wt_24_r2=ctcf_chip_data[["GSM5344491"]]

ctcf_wt_48_r1=ctcf_chip_data[["GSM5344494"]]
ctcf_wt_48_r2=ctcf_chip_data[["GSM5344495"]]


ctcf_wt_24_r1.gr=regioneR::toGRanges(ctcf_wt_24_r1[,c(1:4)])
ctcf_wt_24_r2.gr=regioneR::toGRanges(ctcf_wt_24_r2[,c(1:4)])
ctcf_wt_48_r1.gr=regioneR::toGRanges(ctcf_wt_48_r1[,c(1:4)])
ctcf_wt_48_r2.gr=regioneR::toGRanges(ctcf_wt_48_r2[,c(1:4)])

ctcf_wt_24.gr=c(ctcf_wt_24_r1.gr,ctcf_wt_24_r2.gr)
ctcf_wt_24.gr=GenomicRanges::reduce(ctcf_wt_24.gr)
names=paste0("ctcf_24_peak_",c(1:length(ctcf_wt_24.gr)))
mcols(ctcf_wt_24.gr)$id=names

ctcf_wt_48.gr=c(ctcf_wt_48_r1.gr,ctcf_wt_48_r2.gr)
ctcf_wt_48.gr=GenomicRanges::reduce(ctcf_wt_48.gr)
names=paste0("ctcf_48_peak_",c(1:length(ctcf_wt_48.gr)))
mcols(ctcf_wt_48.gr)$id=names



## ---- data-loops

BEC_loops=file.path(datadir,"loops/juicer_compare_lists/LEC_vs_BEC/BEC/BEC_vs_LEC.common_AB25kb.nooverlap_40kb.bedpe")
LEC_loops=file.path(datadir,"loops/juicer_compare_lists/LEC_vs_BEC/LEC/LEC_vs_BEC.common_AB25kb.nooverlap_40kb.bedpe")

common_loops1=file.path(datadir,"loops/juicer_compare_lists/LEC_vs_BEC/intersection/BEC_vs_LEC.hiccups.common_AB25kb.0_4.overlap.bedpe")
common_loops2=file.path(datadir,"loops/juicer_compare_lists/LEC_vs_BEC/intersection/LEC_vs_BEC.hiccups.common_AB25kb.0_4.overlap.bedpe")


loops_bec=read.csv(BEC_loops, sep="\t", skip = 2,header=FALSE)
loops_lec=read.csv(LEC_loops, sep="\t", skip = 2,header=FALSE)

loops_common1=read.csv(common_loops1, sep="\t", skip = 2,header=FALSE)
loops_common2=read.csv(common_loops2, sep="\t", skip = 2,header=FALSE)

loops_common=rbind(loops_common1,loops_common2)

# this table is LEC loops and BEC loops after rbind corresponding loops
loops_common_uniq_coords=loops_common|>distinct(V1,V2,V3,V4,V5,V6,V13,V14,V15,V16,V17,V18 , .keep_all = TRUE)

# > nrow(loops_common_uniq_coords)
# [1] 7061


loops_common_fmt=as.data.frame(loops_common_uniq_coords|>rowwise()|>
        dplyr::mutate(start1=min(V2,V14),end1=max(V3,V15),start2=min(V5,V17),end2=max(V6,V18))|>
        dplyr::select(V1,start1,end1,V4,start2,end2,V7,V8,V9,V10,V22)|>
        dplyr::mutate(V24="common")|>
        dplyr::rename(V2=start1,V3=end1,V5=start2,V6=end2))


## ---- data-atac

atac_resdir=file.path(projrootdir,"results",paste0("ATACseq_",file_pref),
        "consensus_intrareplicate_merged","DA_all_peaks")
atac_res_DA_pth=file.path(atac_resdir,paste("ATAC_peaks_annot_closest.ensembl.DA_groups",file_pref,"tsv",sep="."))

atacseq_res=read.delim(atac_res_DA_pth, sep="\t", header=TRUE, quote = "")

#atacseq_res$SYMBOL.closest=dplyr::na_if(atacseq_res$SYMBOL.closest, "<NA>")

atacseq_res.lec=atacseq_res|>
        dplyr::filter(FDR<atac_FDR_CO & abs(logFC)>=atac_lFC_CO)|>
        dplyr::filter(logFC>0)|>
        dplyr::filter(logCPM>atac_lCPM_CO)

atacseq_res.bec=atacseq_res|>
        dplyr::filter(FDR<atac_FDR_CO & abs(logFC)>=atac_lFC_CO)|>
        dplyr::filter(logFC<0)|>
        dplyr::filter(logCPM>atac_lCPM_CO)


# format change to be used with dplyranges
# objects for overlaps with TADs

atacseq_res$Chr=gsub("chr","",atacseq_res$Chr)

atac_da_df=atacseq_res|>
        dplyr::select(Chr,Start,End,PeakID,logFC,FDR)|>
        dplyr::rename(id=PeakID)

atac_peaks.gr=regioneR::toGRanges(as.data.frame(atac_da_df))

atac_peaks.up.gr=atac_peaks.gr|>
        plyranges::filter(id%in%atacseq_res.lec$PeakID)

atac_peaks.dn.gr=atac_peaks.gr|>
        plyranges::filter(id%in%atacseq_res.bec$PeakID)



## ---- hic-loops-summary

loopplotdir=file.path(resdir,"loop_summary")
dir.create(loopplotdir,recursive=TRUE)


set.seed(314)
overlaps.loops= euler(c("LEC"=nrow(loops_lec),"BEC"=nrow(loops_bec),
                                 "LEC&BEC"=nrow(loops_common_fmt)), shape = "ellipse")

venn.overlaps.loops=plot(overlaps.loops,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

fpref="loops_venn_LEC_BEC"
save.plots(plot_ts=venn.overlaps.loops,savedir=loopplotdir, pref=fpref)


loops_lec.len=loops_lec|>
	rowwise()|>
	dplyr::mutate(length=(V6-V2)/1000)|>
	dplyr::mutate(smpl="LEC")|>
	dplyr::select(length,smpl)

loops_bec.len=loops_bec|>
	rowwise()|>
	dplyr::mutate(length=(V6-V2)/1000)|>
	dplyr::mutate(smpl="BEC")|>
	dplyr::select(length,smpl)

loops_common.len=loops_common_fmt|>
	rowwise()|>
	dplyr::mutate(length=(V6-V2)/1000)|>
	dplyr::mutate(smpl="common")|>
	dplyr::select(length,smpl)

loops.all=rbind(loops_lec.len,loops_bec.len,loops_common.len)

scale_col3=c(BEC_dark,"#404040",LEC_dark)
scale_fil3=c(BEC_light,"#BABABA",LEC_light)


plot_title="Loop length"
hist_loops=ggplot(loops.all, aes(x=length, color=smpl, fill=smpl)) +
                geom_histogram(alpha=0.3, bins=50, position="identity") +
                theme_bw() +
                theme(aspect.ratio = 1/1.617)+  
                xlab("Length in kbps") +
                ggtitle(plot_title) +
                scale_color_manual(values=scale_col3)+
                scale_fill_manual(values=scale_fil3)+
                scale_x_continuous(labels = comma)+
                theme(text = element_text(size = 20))

hist_loops.1=ggplot(loops.all, aes(x=length, color=smpl, fill=smpl)) +
                geom_histogram(alpha=0.3, bins=50, position="identity") +
                theme_bw() +
                theme(aspect.ratio = 1/1.617)+  
                xlab("Length in kbps") +
                ggtitle(plot_title) +
                scale_color_manual(values=scale_col3)+
                scale_fill_manual(values=scale_fil3)+
                scale_x_continuous(labels = comma)+
                theme(text = element_text(size = 20))+
                xlim(0,2000)


fpref="loops_histogram_LEC_BEC_common"
save.plots(plot_ts=hist_loops.1,savedir=loopplotdir, pref=fpref)



## ---- section-reproducible-tads

### before
### section-consensus-tads

########################
########################
###########    tissue consensus TADs i.e. reproducible TADs within each tissue
########################


## ---- TADs-bnds-functions



get_TADs_from_bounds=function(tads_df=tads_df, bnds=bnds, type=c("both","any")){
        require(dplyr)

        if(type=="both"){
                tads_df_res=tads_df|>
                        dplyr::filter(bnd_id_5%in%bnds & bnd_id_3%in%bnds)
        }
        if(type=="any"){
                tads_df_res=tads_df|>
                        dplyr::filter(bnd_id_5%in%bnds | bnd_id_3%in%bnds)
        }

        return(tads_df_res)
}

#for "dif" the tad_df must contain column `dif_bnd` of values: bnd_5_3 bnd_5 bnd_3
#this is not checked at the moment

# returns bounds_df which correspond to TADs in tads_df$bnd_id_5 and tads_df$bnd_id_3

filter_bounds_from_TADs=function(tads_df=tads_df, bnds_df=bnds_df, type=c("both","dif")){
        require(dplyr)

        if(type=="both"){
                bnds_lst=unique(c(tads_df$bnd_id_5,tads_df$bnd_id_3))
                bnds_df_res=bnds_df|>
                        dplyr::filter(bnd_id%in%bnds_lst)
        }
        if(type=="dif"){
                tads_df.bnds.3=tads_df|>
                        dplyr::filter(dif_bnd=="bnd_3" | dif_bnd=="bnd_5_3")|>
                        pull(bnd_id_3)

                tads_df.bnds.5=tads_df|>
                        dplyr::filter(dif_bnd=="bnd_5" | dif_bnd=="bnd_5_3")|>
                        pull(bnd_id_5)

                bnds_lst=unique(c(tads_df.bnds.5,tads_df.bnds.3))

                bnds_df_res=bnds_df|>
                        dplyr::filter(bnd_id%in%bnds_lst)
        }

        return(bnds_df_res)
}


#returns bounds list present in tads_df
get_bounds_from_TADs=function(tads_df=tads_df, type=c("both","dif")){
        require(dplyr)

        if(type=="both"){
                bnds_lst=unique(c(tads_df$bnd_id_5,tads_df$bnd_id_3))

        }
        if(type=="dif"){
                tads_df.bnds.3=tads_df|>
                        dplyr::filter(dif_bnd=="bnd_3" | dif_bnd=="bnd_5_3")|>
                        pull(bnd_id_3)

                tads_df.bnds.5=tads_df|>
                        dplyr::filter(dif_bnd=="bnd_5" | dif_bnd=="bnd_5_3")|>
                        pull(bnd_id_5)

                bnds_lst=unique(c(tads_df.bnds.5,tads_df.bnds.3))
        }

        return(bnds_lst)
}





## ---- TADcompare-cons-megaTADs

#for overlaps
file_string1="OnTAD_megamap"
txt_string1="OnTAD megamap"

file_string1c="OnTAD_megamap_reproducible"
txt_string1c="OnTAD megamap-TADs reproducible"


LEC_bnds.cons=mega_bnds.102_104.ontad.50.all|>dplyr::filter(Status=="Non-Differential")|>pull(bnd_id)
tads_LEC.cons=get_TADs_from_bounds(tads_df=tads_LEC,bnds=LEC_bnds.cons,type="both")

BEC_bnds.cons=mega_bnds.101_103.ontad.50.all|>dplyr::filter(Status=="Non-Differential")|>pull(bnd_id)
tads_BEC.cons=get_TADs_from_bounds(tads_df=tads_BEC,bnds=BEC_bnds.cons,type="both")

format_consTADs=function(tads_df=tads_df,ident=ident){
        tads_df=tads_df|>
                dplyr::select(!c(startpos,endpos))|>
                dplyr::mutate(tad_len_win=(end-start)/binsize)|>
                dplyr::mutate(smpl=ident)

        return(tads_df)
}


tads_LEC.cons=format_consTADs(tads_df=tads_LEC.cons,ident="LEC")
tads_BEC.cons=format_consTADs(tads_df=tads_BEC.cons,ident="BEC")


## ---- TADcompare-cons-megaTADs-combine

##### these three following functions perform TAD overlaps resulting in df of "equivalent" and unique TADs
##### wrapped into one proj specific function below

#get overlapping TADs (identical and with boundaries within wiggle room)
overlap_TADs=function(tads_df.1=tads_df.1,tads_df.2=tad_df.2, wiggle_bins=wiggle_bins){
        require(regioneR)

        id.1=tads_df.1$smpl[1]
        id.2=tads_df.2$smpl[1]
        id.1_2=paste(id.1,id.2,sep="_")

        tads_1.gr=regioneR::toGRanges(tads_df.1|>
                                dplyr::select(chr,start,end,tad_id) )

        tads_2.gr=regioneR::toGRanges(tads_df.2|>
                                dplyr::select(chr,start,end,tad_id) )

        overlap=regioneR::overlapRegions(tads_1.gr,tads_2.gr, colA=c("tad_id"),colB=c("tad_id"))
        colnames(overlap)[c(6,7)]=c(paste0("tad_id_",id.1),paste0("tad_id_",id.2))

        #assign a combo tad_id encompassing longest region - for list comparisons only
        overlap=overlap|>
                rowwise()|>
                dplyr::mutate(start_min=min(startA,startB))|>
                dplyr::mutate(end_max=max(endA,endB))|>
                dplyr::mutate(tad_id_wg=paste(chr,start_min,end_max,sep="_"))|>
                dplyr::select(!c(start_min,end_max))#|>
                #dplyr::mutate(smpl=id.1_2)

        #identical tads
        overlap.filt.ident=overlap|>
                dplyr::filter(type=="equal")

        #check if boundary in wiggle boundaries
        overlap.filt.wg=overlap|>
                dplyr::filter(abs(startA-startB)<=ontad_wiggle_bins*binsize)|>
                dplyr::filter(abs(endA-endB)<=ontad_wiggle_bins*binsize)

        df_lst=list()
        df_lst[["ident"]]=overlap.filt.ident
        df_lst[["wg"]]=overlap.filt.wg

        return(df_lst)
}

wiggle_bins=3
tads_overlap.unfilt=overlap_TADs(tads_df.1=tads_LEC,tads_df.2=tads_BEC,wiggle_bins=wiggle_bins)

tads_overlap.unfilt.ident=tads_overlap.unfilt[["ident"]]
tads_overlap.unfilt.wg=tads_overlap.unfilt[["wg"]]

tads_overlap.cons=overlap_TADs(tads_df.1=tads_LEC.cons,tads_df.2=tads_BEC.cons,wiggle_bins=wiggle_bins)

tads_overlap.cons.ident=tads_overlap.cons[["ident"]]
tads_overlap.cons.wg=tads_overlap.cons[["wg"]]



## ---- TADcompare-cons-megaTADs-fmt

# get tads and boundary_ids (mega coords) demarcating conserved TADs (for all TADs)


# for table with overlaps
n_TADs_raw_mega_LEC=nrow(tads_LEC)
n_TADs_raw_mega_BEC=nrow(tads_BEC)

n_TADs_101_103_cons_BEC=nrow(tads_LEC.cons)
n_TADs_102_104_cons_LEC=nrow(tads_BEC.cons)

n_TADs_raw_mega_LEC.1=nrow(tads_LEC|>dplyr::filter(TAD_level==1))
n_TADs_raw_mega_BEC.1=nrow(tads_BEC|>dplyr::filter(TAD_level==1))

n_TADs_101_103_cons_BEC.1=nrow(tads_BEC.cons|>dplyr::filter(TAD_level==1))
n_TADs_102_104_cons_LEC.1=nrow(tads_LEC.cons|>dplyr::filter(TAD_level==1))


TAD_cons_summary_df=data.frame(rbind(
        c("unfiltered TADs",n_TADs_raw_mega_LEC,n_TADs_raw_mega_BEC,n_TADs_raw_mega_LEC.1,n_TADs_raw_mega_BEC.1),
        c("reproducible TADs",n_TADs_102_104_cons_LEC,n_TADs_101_103_cons_BEC,n_TADs_102_104_cons_LEC.1,n_TADs_101_103_cons_BEC.1)
        ))
colnames(TAD_cons_summary_df)=c("Type","LEC","BEC", "LEC lv1","BEC lv1")

## save consensus TADs as tables and bed files

save_TADs_bed=function(tads_df=tads_df,fname_str=fname_str,outdir=outdir){
        
        fname=paste(fname_str,"tsv",sep=".")
        save.table(df=tads_df,file=fname,dir=outdir)



        tads_bed=tads_df|>
                dplyr::select(chr,start,end,tad_id,TADscore)|>
                dplyr::mutate(strand=".")

        fname_bed=paste(fname_str,"bed",sep=".")
        save.bed(df=tads_bed,dir=out_resdir,file=fname_bed)

        #save bed with each TAD level, from 1 to n
        for (i in c(1:max(tads_df$TAD_level))){

                tads_bed.i=tads_df|>
                        dplyr::filter(TAD_level==i)|>
                        dplyr::select(chr,start,end,tad_id,TADscore)|>
                        dplyr::mutate(strand=".")

                tad_lvl=paste0("lvl_",i)
                fname_bed=paste(fname_str,tad_lvl,"bed",sep=".")
                save.bed(df=tads_bed.i,dir=out_resdir,file=fname_bed)
        }

}

fname_str1="LEC_reproducible_TADs"
out_resdir=file.path(resdir,fname_str1)
dir.create(out_resdir,recursive=TRUE)
fname_str=paste(fname_str1,file_string1,file_pref,sep=".")
save_TADs_bed(tads_df=tads_LEC.cons,fname_str=fname_str,outdir=out_resdir)

fname_str1="BEC_reproducible_TADs"
out_resdir=file.path(resdir,fname_str1)
dir.create(out_resdir,recursive=TRUE)
fname_str=paste(fname_str1,file_string1,file_pref,sep=".")
save_TADs_bed(tads_df=tads_BEC.cons,fname_str=fname_str,outdir=out_resdir)





## ---- smpl-tads-venn

#for overlaps
file_string1="OnTAD_megamap"
txt_string1="OnTAD megamap"

file_string2="reproducible" #or unfilt
txt_string2="reproducible"

file_string3="identical" #or wiggle_n_bins
txt_string3="identical"


get_overlaps_euler_2=function(df_tads.lec=df_tads.lec,df_tads.bec=df_tads.bec, df_overlap=df_overlap, level=c("all","1")){

        if (level!="all"){
        df_tads.lec=df_tads.lec|>dplyr::filter(TAD_level==level)
        df_tads.bec=df_tads.bec|>dplyr::filter(TAD_level==level)
        }

  df_overlap=df_overlap|>dplyr::filter(tad_id_LEC%in%df_tads.lec$tad_id | tad_id_BEC%in%df_tads.bec$tad_id)

  LEC_TADs=df_tads.lec|>dplyr::filter(!tad_id%in%df_overlap$tad_id_LEC)
  BEC_TADs=df_tads.bec|>dplyr::filter(!tad_id%in%df_overlap$tad_id_BEC)

  n_LEC=length(unique(LEC_TADs$tad_id))
  n_BEC=length(unique(BEC_TADs$tad_id))
  n_LEC_BEC=length(unique(df_overlap$tad_id_wg))

  set.seed(314)
  tad_intersections = euler(c("LEC" = n_LEC, "BEC" = n_BEC,
                                 "LEC&BEC" = n_LEC_BEC), shape = "ellipse")

  return(tad_intersections)

}


# all levels 


venn_overlaps.cons_tads_LEC_BEC.wg=get_overlaps_euler_2(df_tads.lec=tads_LEC.cons,df_tads.bec=tads_BEC.cons,df_overlap=tads_overlap.cons.wg,level="all")
venn_overlaps.cons_tads_LEC_BEC.ident=get_overlaps_euler_2(df_tads.lec=tads_LEC.cons,df_tads.bec=tads_BEC.cons,df_overlap=tads_overlap.cons.ident,level="all")

venn_overlaps.unfilt_tads_LEC_BEC.wg=get_overlaps_euler_2(df_tads.lec=tads_LEC,df_tads.bec=tads_BEC,df_overlap=tads_overlap.unfilt.wg,level="all")
venn_overlaps.unfilt_tads_LEC_BEC.ident=get_overlaps_euler_2(df_tads.lec=tads_LEC,df_tads.bec=tads_BEC,df_overlap=tads_overlap.unfilt.ident,level="all")




venn.cons.wg=plot(venn_overlaps.cons_tads_LEC_BEC.wg,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

venn.cons.ident=plot(venn_overlaps.cons_tads_LEC_BEC.ident,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

venn.unfilt.wg=plot(venn_overlaps.unfilt_tads_LEC_BEC.wg,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

venn.unfilt.ident=plot(venn_overlaps.unfilt_tads_LEC_BEC.ident,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

# level 1 TADs



venn_overlaps.cons_tads_LEC_BEC.wg.1=get_overlaps_euler_2(df_tads.lec=tads_LEC.cons,df_tads.bec=tads_BEC.cons,df_overlap=tads_overlap.cons.wg,level="1")
venn_overlaps.cons_tads_LEC_BEC.ident.1=get_overlaps_euler_2(df_tads.lec=tads_LEC.cons,df_tads.bec=tads_BEC.cons,df_overlap=tads_overlap.cons.ident,level="1")

venn_overlaps.unfilt_tads_LEC_BEC.wg.1=get_overlaps_euler_2(df_tads.lec=tads_LEC,df_tads.bec=tads_BEC,df_overlap=tads_overlap.unfilt.wg,level="1")
venn_overlaps.unfilt_tads_LEC_BEC.ident.1=get_overlaps_euler_2(df_tads.lec=tads_LEC,df_tads.bec=tads_BEC,df_overlap=tads_overlap.unfilt.ident,level="1")


venn.cons.wg.1=plot(venn_overlaps.cons_tads_LEC_BEC.wg.1,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

venn.cons.ident.1=plot(venn_overlaps.cons_tads_LEC_BEC.ident.1,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

venn.unfilt.wg.1=plot(venn_overlaps.unfilt_tads_LEC_BEC.wg.1,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))

venn.unfilt.ident.1=plot(venn_overlaps.unfilt_tads_LEC_BEC.ident.1,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))


#save them as files

wiggle_n_bins=paste0("wiggle_",ontad_wiggle_bins)

file_string2="reproducible" #or unfilt
file_string3="identical" #or wiggle_n_bins=paste0("wiggle_",ontad_wiggle_bins)


vennplotdir=file.path(resdir,"venn_LEC_BEC_overlaps")
dir.create(vennplotdir,recursive=TRUE)

file_string2="reproducible" #or unfilt
file_string3=wiggle_n_bins
fpref=paste("venn_LEC_BEC",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.cons.wg,savedir=vennplotdir, pref=fpref)

file_string2="unfilt" #or unfilt
file_string3=wiggle_n_bins
fpref=paste("venn_LEC_BEC",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.unfilt.wg,savedir=vennplotdir, pref=fpref)

file_string2="reproducible" #or unfilt
file_string3="identical"
fpref=paste("venn_LEC_BEC",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.cons.ident,savedir=vennplotdir, pref=fpref)

file_string2="unfilt" #or unfilt
file_string3="identical"
fpref=paste("venn_LEC_BEC",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.unfilt.ident,savedir=vennplotdir, pref=fpref)

## lv1
file_string2="reproducible" #or unfilt
file_string3=wiggle_n_bins
fpref=paste("venn_LEC_BEC","lvl_1",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.cons.wg.1,savedir=vennplotdir, pref=fpref)

file_string2="unfilt" #or unfilt
file_string3=wiggle_n_bins
fpref=paste("venn_LEC_BEC","lvl_1",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.unfilt.wg.1,savedir=vennplotdir, pref=fpref)

file_string2="reproducible" #or unfilt
file_string3="identical"
fpref=paste("venn_LEC_BEC","lvl_1",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.cons.ident.1,savedir=vennplotdir, pref=fpref)

file_string2="unfilt" #or unfilt
file_string3="identical"
fpref=paste("venn_LEC_BEC","lvl_1",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=venn.unfilt.ident.1,savedir=vennplotdir, pref=fpref)


## table

get_cnts_table_tads=function(tads.lec,tads.bec,overlap_df) {

                #all levels
                tad_cnts.all=get_overlaps_euler_2(
                        df_tads.lec=tads.lec,df_tads.bec=tads.bec,
                        df_overlap=overlap_df, level="all")

                tad_cnts.all.v=tad_cnts.all$original

                #lvl1
                tad_cnts.all.1=get_overlaps_euler_2(
                        df_tads.lec=tads.lec,df_tads.bec=tads.bec,
                        df_overlap=overlap_df, level="1")

                tad_cnts.all.1.v=tad_cnts.all.1$original


                tad_cnts=c(tad_cnts.all.v,tad_cnts.all.1.v)

        return(tad_cnts)
}


n.unfilt.ident=get_cnts_table_tads(tads.lec=tads_LEC,tads.bec=tads_BEC,overlap_df=tads_overlap.unfilt.ident)
n.cons.ident=get_cnts_table_tads(tads.lec=tads_LEC.cons,tads.bec=tads_BEC.cons,overlap_df=tads_overlap.cons.ident)
n.unfilt.wg=get_cnts_table_tads(tads.lec=tads_LEC,tads.bec=tads_BEC,overlap_df=tads_overlap.unfilt.wg)
n.cons.wg=get_cnts_table_tads(tads.lec=tads_LEC.cons,tads.bec=tads_BEC.cons,overlap_df=tads_overlap.cons.wg)

TAD_overlap_summary_df=data.frame(cbind(
                        c("Tissue / TAD","LEC, all levels","BEC, all levels","common, all levels","LEC, level 1","BEC, level 1","common, level 1"),
                        c("unfiltered TADs, identical coordinates", n.unfilt.ident),
                        c("unfiltered TADs, wiggle coordinates", n.unfilt.wg),
                        c("reproducible TADs, identical coordinates", n.cons.ident),
                        c("reproducible TADs, wiggle coordinates", n.cons.wg) ))

colnames(TAD_overlap_summary_df)=TAD_overlap_summary_df[1,]




## ---- TAD-len-hist


histplotdir=file.path(resdir,"histogram_TAD_length")
dir.create(histplotdir,recursive=TRUE)


## LEC and BEC superimposed

# saves the plots
#  for crop xlim(0,3.5e6)


plot_TADlen_hist2=function(tads_lec=tads_lec,tads_bec=tads_bec,plot_title=plot_title,
        plotdir=plotdir,file_str_pl=file_str_pl){
        
        tads_lec=tads_lec|>
                dplyr::select(chr,start,end,TAD_level)|>
                dplyr::mutate(length=end-start)|>
                dplyr::mutate(smpl="LEC")

        tads_lec.1=tads_lec|>dplyr::filter(TAD_level==1)

        tads_bec=tads_bec|>
                dplyr::select(chr,start,end,TAD_level)|>
                dplyr::mutate(length=end-start)|>
                dplyr::mutate(smpl="BEC")

        tads_bec.1=tads_bec|>dplyr::filter(TAD_level==1)


        allTADs=rbind(tads_lec, tads_bec)
        allTADs.1=rbind(tads_lec.1, tads_bec.1)
        
        plot_title_all=paste0(plot_title,"\n","all TAD levels")
        plot_title_lv1=paste0(plot_title,"\n","TAD level 1")

        hist_all_tads=ggplot(allTADs, aes(x=length, color=smpl, fill=smpl)) +
                geom_histogram(alpha=0.3, bins=50, position="identity") +
                theme_bw() +
                theme(aspect.ratio = 1/1.617)+  
                xlab("Length in bps") +
                ggtitle(plot_title_all) +
                scale_color_manual(values=scale_col)+
                scale_fill_manual(values=scale_fil)+
                scale_x_continuous(labels = comma)+
                theme(text = element_text(size = 20))

        hist_all_tads.1=ggplot(allTADs.1, aes(x=length, color=smpl, fill=smpl)) +
                geom_histogram(alpha=0.3, bins=50, position="identity") +
                theme_bw() +
                theme(aspect.ratio = 1/1.617)+  
                xlab("Length in bps") +
                ggtitle(plot_title_lv1) +
                scale_color_manual(values=scale_col)+
                scale_fill_manual(values=scale_fil)+
                scale_x_continuous(labels = comma)+
                theme(text = element_text(size = 20))

        plots_list=list()
        plots_list[["all"]]=hist_all_tads
        plots_list[["lvl1"]]=hist_all_tads.1
        plots_list[["all-crop"]]=hist_all_tads + xlim(0,3.5e6)
        plots_list[["lvl1-crop"]]=hist_all_tads.1 + xlim(0,3.5e6)

        # save to files
        for (i in names(plots_list)){

                fname=paste(file_str_pl,i,sep=".")

                for (devplot in c("pdf","eps")){
                        fname_full=paste(fname,devplot,sep=".")
                        ggsave(plot=plots_list[[i]], filename=fname_full, path=histplotdir, device=devplot)
                }
        }

        plots_list=list()
        plots_list[["all"]]=hist_all_tads+theme(text = element_text(size = 15)) + theme(legend.position="bottom")+scale_x_continuous(labels = comma)
        plots_list[["lvl1"]]=hist_all_tads.1+theme(text = element_text(size = 15)) + theme(legend.position="bottom")+scale_x_continuous(labels = comma)
        plots_list[["all-crop"]]=hist_all_tads+theme(text = element_text(size = 15))+theme(legend.position="bottom")+scale_x_continuous(labels = comma)+xlim(0,3.5e6)
        plots_list[["lvl1-crop"]]=hist_all_tads.1+theme(text = element_text(size = 15))+theme(legend.position="bottom")+scale_x_continuous(labels = comma)+xlim(0,3.5e6)

        return(plots_list)
}



scale_col=c(BEC_dark,LEC_dark)
scale_fil=c(BEC_light,LEC_light)


plot_title1="Reproducible TADs"
file_str_pl=paste("hist_TADlength_superimposed",file_string1,"reproducible",sep=".")
hist_LEC_BEC.cons=plot_TADlen_hist2(tads_lec=tads_LEC.cons,tads_bec=tads_BEC.cons,
        plot_title=plot_title1,plotdir=histplotdir,file_str_pl=file_str_pl)


plot_title1="Unfiltered TADs"
file_str_pl=paste("hist_TADlength_superimposed",file_string1,"unfiltered",sep=".")
hist_LEC_BEC.unfilt=plot_TADlen_hist2(tads_lec=tads_LEC,tads_bec=tads_BEC,
        plot_title=plot_title1,plotdir=histplotdir,file_str_pl=file_str_pl)



## per sample separately

# saves the plots
#  for crop xlim(0,3.5e6)


plot_TADlen_hist=function(tads=tads,plot_title=plot_title,
        plotdir=plotdir,file_str_pl=file_str_pl){
        
        tads=tads|>
                dplyr::select(chr,start,end,TAD_level)|>
                dplyr::mutate(length=end-start)

 
        hist_tads.fullx=ggplot(tads, aes(x=length)) +
                geom_histogram(color=hcol,fill=hfill, alpha=0.5, bins=30) +
                theme_bw() +
                theme(aspect.ratio = 1/1.617)+  
                xlab("Length in bps") + ylab(paste0("Count, n=",nrow(tads)))+
                ggtitle(plot_title)+
                theme(text = element_text(size = 20))+
                                scale_x_continuous(labels = comma)


        hist_tads.cropx=hist_tads.fullx+ 
                        xlim(0,3.5e6)

        tads.1=tads|>dplyr::filter(TAD_level==1)
        
        plot_title=paste0(plot_title,"\nlevel 1")
        hist_tads.1.fullx=ggplot(tads.1, aes(x=length)) +
                geom_histogram(color=hcol,fill=hfill, alpha=0.5, bins=30) +
                theme_bw() +
                theme(aspect.ratio = 1/1.617)+  
                xlab("Length in bps") + ylab(paste0("Count, n=",nrow(tads)))+
                ggtitle(plot_title)+
                theme(text = element_text(size = 20))+
                                scale_x_continuous(labels = comma)


        hist_tads.1.cropx=hist_tads.1.fullx+ 
                        xlim(0,3.5e6)


        plots_list=list()
        plots_list[["all"]]=hist_tads.fullx
        plots_list[["lvl1"]]=hist_tads.1.fullx
        plots_list[["all-crop"]]= hist_tads.cropx
        plots_list[["lvl1-crop"]]=hist_tads.1.cropx


        # save to files
        for (i in names(plots_list)){

                fname=paste(file_str_pl,i,sep=".")

                for (devplot in c("pdf","eps")){
                        fname_full=paste(fname,devplot,sep=".")
                        ggsave(plot=plots_list[[i]], filename=fname_full, path=histplotdir, device=devplot)
                }
        }

        plots_list=list()
        plots_list[["all"]]=hist_tads.fullx+theme(text = element_text(size = 15)) + theme(legend.position="bottom")+scale_x_continuous(labels = comma)
        plots_list[["lvl1"]]=hist_tads.1.fullx+theme(text = element_text(size = 15)) + theme(legend.position="bottom")+scale_x_continuous(labels = comma)
        plots_list[["all-crop"]]=hist_tads.cropx+xlim(0,3.5e6)+theme(text = element_text(size = 15))+theme(legend.position="bottom")+scale_x_continuous(labels = comma)
        plots_list[["lvl1-crop"]]=hist_tads.1.cropx+xlim(0,3.5e6)+theme(text = element_text(size = 15))+theme(legend.position="bottom")+scale_x_continuous(labels = comma)

        return(plots_list)
}


scale_col=c(BEC_dark,LEC_dark)
scale_fil=c(BEC_light,LEC_light)

hcol=LEC_dark
hfill=LEC_light

plot_title1="LEC Reproducible TADs"
file_str_pl=paste("hist_TADlength_LEC",file_string1,"reproducible",sep=".")
hist_LEC.cons=plot_TADlen_hist(tads=tads_LEC.cons,
        plot_title=plot_title1,plotdir=histplotdir,file_str_pl=file_str_pl)

plot_title1="LEC Unfiltered TADs"
file_str_pl=paste("hist_TADlength_LEC",file_string1,"unfiltered",sep=".")
hist_LEC.unfilt=plot_TADlen_hist(tads=tads_LEC,
        plot_title=plot_title1,plotdir=histplotdir,file_str_pl=file_str_pl)

hcol=BEC_dark
hfill=BEC_light

plot_title1="BEC Reproducible TADs"
file_str_pl=paste("hist_TADlength_BEC",file_string1,"reproducible",sep=".")
hist_BEC.cons=plot_TADlen_hist(tads=tads_BEC.cons,
        plot_title=plot_title1,plotdir=histplotdir,file_str_pl=file_str_pl)

plot_title1="BEC Unfiltered TADs"
file_str_pl=paste("hist_TADlength_BEC",file_string1,"unfiltered",sep=".")
hist_BEC.unfilt=plot_TADlen_hist(tads=tads_BEC,
        plot_title=plot_title1,plotdir=histplotdir,file_str_pl=file_str_pl)



plot_figure1=list()
plot_figure1[["LEC"]]=hist_LEC.unfilt[["lvl1"]]
plot_figure1[["BEC"]]=hist_BEC.unfilt[["lvl1"]]
plot_figure1[["LEC-xlim"]]=hist_LEC.unfilt[["lvl1-crop"]]+xlim(0,3.5e6)
plot_figure1[["BEC-xlim"]]=hist_BEC.unfilt[["lvl1-crop"]]+xlim(0,3.5e6)


plot_figure2=list()
plot_figure2[["LEC"]]=hist_LEC.unfilt[["all"]]
plot_figure2[["BEC"]]=hist_BEC.unfilt[["all"]]
plot_figure2[["LEC-xlim"]]=hist_LEC.unfilt[["all-crop"]]+xlim(0,3.5e6)
plot_figure2[["BEC-xlim"]]=hist_BEC.unfilt[["all-crop"]]+xlim(0,3.5e6)

plot_figure3=list()
plot_figure3[["LEC"]]=hist_LEC.cons[["lvl1"]]
plot_figure3[["BEC"]]=hist_BEC.cons[["lvl1"]]
plot_figure3[["LEC-xlim"]]=hist_LEC.cons[["lvl1-crop"]]+xlim(0,3.5e6)
plot_figure3[["BEC-xlim"]]=hist_BEC.cons[["lvl1-crop"]]+xlim(0,3.5e6)


plot_figure4=list()
plot_figure4[["LEC"]]=hist_LEC.cons[["all"]]
plot_figure4[["BEC"]]=hist_BEC.cons[["all"]]
plot_figure4[["LEC-xlim"]]=hist_LEC.cons[["all-crop"]]+xlim(0,3.5e6)
plot_figure4[["BEC-xlim"]]=hist_BEC.cons[["all-crop"]]+xlim(0,3.5e6)



## ---- TAD-len-hist-pub


Yang_m_l="#33999975"
Yang_m_d="#339999"

Yang_b_l="#FF990075"
Yang_b_d="#FF9900"

Franke_l="#FF99CC75"
Franke_d="#FF99CC"



plot_title1="TAD length in different data sets\n(unfiltered TADs)"
plot_title2="TAD length in different data sets\n(reproducible TADs)"



Franke_1_tads$smpl="embryo_24h"
brain_tads$smpl="brain"
lec_tads=tads_LEC|>dplyr::filter(TAD_level==1)|>dplyr::select(chr,start,end)|>dplyr::mutate(length=end-start)|>dplyr::mutate(smpl="LEC")
bec_tads=tads_BEC|>dplyr::filter(TAD_level==1)|>dplyr::select(chr,start,end)|>dplyr::mutate(length=end-start)|>dplyr::mutate(smpl="BEC")
lec_tads_cons=tads_LEC.cons|>dplyr::filter(TAD_level==1)|>dplyr::select(chr,start,end)|>dplyr::mutate(length=end-start)|>dplyr::mutate(smpl="LEC")
bec_tads_cons=tads_BEC.cons|>dplyr::filter(TAD_level==1)|>dplyr::select(chr,start,end)|>dplyr::mutate(length=end-start)|>dplyr::mutate(smpl="BEC")

allTADs=rbind(lec_tads,bec_tads,Franke_1_tads,brain_tads)
allTADs_cons=rbind(lec_tads_cons,bec_tads_cons,Franke_1_tads,brain_tads)


lab_LEC=paste0("LEC n=",nrow(lec_tads))
lab_BEC=paste0("BEC n=",nrow(bec_tads))
lab_Fr=paste0("embryo n=",nrow(Franke_1_tads))
lab_Yb=paste0("brain n=",nrow(brain_tads))

lab_LEC_c=paste0("LEC n=",nrow(lec_tads_cons))
lab_BEC_c=paste0("BEC n=",nrow(bec_tads_cons))



hist_all_tads=ggplot(allTADs, aes(x=length, color=smpl, fill=smpl)) +
    geom_histogram(alpha=0.3, bins=50, position="identity") +
    theme_bw() +
    theme(aspect.ratio = 1/1.617)+  
    xlab("Length in bps") +
    ggtitle(plot_title1) +
    scale_color_manual(values=c(BEC_dark, Yang_b_d,Franke_d,LEC_dark))+
    scale_fill_manual(values=c(BEC_light,Yang_b_l,Franke_l,LEC_light))+
    scale_x_continuous(labels = comma)

hist_all_tads1=hist_all_tads+ 
        annotate("text", x=5500000, y=1000,label=lab_LEC, color=LEC_dark)+
        annotate("text", x=5500000, y=800,label=lab_BEC, color=BEC_dark)+
        annotate("text", x=5500000, y=600,label=lab_Fr, color=Franke_d)+
        annotate("text", x=5500000, y=400,label=lab_Yb, color=Yang_b_d)+
        scale_x_continuous(labels = comma)

hist_all_tads2=hist_all_tads+
        coord_cartesian(xlim=c(0,3e6))

hist_all_tads3=hist_all_tads+
        coord_cartesian(xlim=c(0,3e6))+ 
        annotate("text", x=2000000, y=1000,label=lab_LEC, color=LEC_dark)+
        annotate("text", x=2000000, y=800,label=lab_BEC, color=BEC_dark)+
        annotate("text", x=2000000, y=600,label=lab_Fr, color=Franke_d)+
        annotate("text", x=2000000, y=400,label=lab_Yb, color=Yang_b_d)





plot_title="TAD length in different data sets\n(reproducible TADs)"

hist_all_tads_c=ggplot(allTADs_cons, aes(x=length, color=smpl, fill=smpl)) +
    geom_histogram(alpha=0.3, bins=50, position="identity") +
    theme_bw() +
    theme(aspect.ratio = 1/1.617)+  
    xlab("Length in bps") +
    ggtitle(plot_title2) +
    scale_color_manual(values=c(BEC_dark, Yang_b_d,Franke_d,LEC_dark))+
    scale_fill_manual(values=c(BEC_light,Yang_b_l,Franke_l,LEC_light))+
    scale_x_continuous(labels = comma)


hist_all_tads_c1=hist_all_tads_c+ 
        annotate("text", x=5500000, y=1000,label=lab_LEC_c, color=LEC_dark)+
        annotate("text", x=5500000, y=800,label=lab_BEC_c, color=BEC_dark)+
        annotate("text", x=5500000, y=600,label=lab_Fr, color=Franke_d)+
        annotate("text", x=5500000, y=400,label=lab_Yb, color=Yang_b_d)+
        scale_x_continuous(labels = comma)

hist_all_tads_c2=hist_all_tads_c+ 
        coord_cartesian(xlim=c(0,3e6))

hist_all_tads_c3=hist_all_tads_c+ 
        coord_cartesian(xlim=c(0,3e6)) +
        annotate("text", x=2000000, y=1000,label=lab_LEC, color=LEC_dark)+
        annotate("text", x=2000000, y=800,label=lab_BEC, color=BEC_dark)+
        annotate("text", x=2000000, y=600,label=lab_Fr, color=Franke_d)+
        annotate("text", x=2000000, y=400,label=lab_Yb, color=Yang_b_d)


histplotdir2=file.path(resdir,"histogram_TAD_length_pubdata")
dir.create(histplotdir2,recursive=TRUE)

save_hist=function(gplot=gplot,namepref=namepref){
        
        fname=paste0(file_string1,namepref,".Franke_Yang.TAD_len_histogram")
        for (devplot in c("pdf","eps")){
                fname_full=paste(fname,devplot,sep=".")
                ggsave(plot=gplot, filename=fname_full, path=histplotdir2, device=devplot)
        }

}


save_hist(gplot=hist_all_tads,namepref="unfilt_lvl1_TADs.fl_scale")
save_hist(gplot=hist_all_tads1,namepref="unfilt_lvl1_TADs.fl_scale_n")
save_hist(gplot=hist_all_tads2,namepref="unfilt_lvl1_TADs.crop_scale")
save_hist(gplot=hist_all_tads3,namepref="unfilt_lvl1_TADs.crop_scale_n")

save_hist(gplot=hist_all_tads_c,namepref="reproducible_lvl1_TADs.fl_scale")
save_hist(gplot=hist_all_tads_c1,namepref="reproducible_lvl1_TADs.fl_scale_n")
save_hist(gplot=hist_all_tads_c2,namepref="reproducible_lvl1_TADs.crop_scale")
save_hist(gplot=hist_all_tads_c3,namepref="reproducible_lvl1_TADs.crop_scale_n")



plot_figure3=list()
plot_figure3[["full"]]=hist_all_tads
plot_figure3[["crop"]]=hist_all_tads3
plot_figure3[["full-cons"]]=hist_all_tads_c
plot_figure3[["crop-cons"]]=hist_all_tads_c3



## ---- section-differential-tads

########################
########################
###########    differential reproducible TADs
########################


## ---- diff-TADs-TADcompare



##### prep TAD object


# consensus TADs
# may remove TADs < 20 windows
tad_len_CO_bins=1

tads_LEC.cons.flt=tads_LEC.cons|>dplyr::filter(tad_len_win>=tad_len_CO_bins)
tads_BEC.cons.flt=tads_BEC.cons|>dplyr::filter(tad_len_win>=tad_len_CO_bins)

## select boundaries of consensus TADs to use for subsetting dif boundaries
cons_bnds=c(tads_LEC.cons.flt$bnd_id_5,tads_LEC.cons.flt$bnd_id_3,tads_BEC.cons.flt$bnd_id_5,tads_BEC.cons.flt$bnd_id_3)
cons_bnds=unique(cons_bnds)

# df of  boundary scores - each boundary interrogated in each of the three windows
#provisional df for joins
bnds_gap_scores=data.frame(bnd_id=unique(c(all_bnd_calls[[20]][["bnds"]]$bnd_id,all_bnd_calls[[35]][["bnds"]]$bnd_id,all_bnd_calls[[50]][["bnds"]]$bnd_id)))

for (win in window_sizes ){

  bnds_gap_score=all_bnd_calls[[win]][["bnds"]]|>
    dplyr::select(bnd_id, starts_with("Gap_Score."))

  bnds_gap_scores=full_join(bnds_gap_scores,bnds_gap_score, by="bnd_id")

}

bnds_gap_scores.dif=bnds_gap_scores|>dplyr::filter(if_any(starts_with("Gap_Score."), ~abs(.x)>1.8))
bnds_gap_scores.dif.cons=bnds_gap_scores.dif|>dplyr::filter(bnd_id%in%cons_bnds)

bnds_gap_scores.dif.cons.20=bnds_gap_scores.dif.cons|>dplyr::filter(abs(Gap_Score.20)>1.8)
bnds_gap_scores.dif.cons.35=bnds_gap_scores.dif.cons|>dplyr::filter(abs(Gap_Score.35)>1.8)
bnds_gap_scores.dif.cons.50=bnds_gap_scores.dif.cons|>dplyr::filter(abs(Gap_Score.50)>1.8)



## select appropriate window for boundary gap score based on TAD length
## should be min in LEC and BEC
## if given boundary delimits TADs of different levels (having different length), select max TAD level i.e. min length

## do this in each tissue separately


# get dif TADs from dif bounds in three TAD length strata (match the TAD length to appropriate boundary score window)
# starting from consensus TADs filtered by length


get_difTADs=function(bnds_gapscores=bnds_gapscores, cons_tads=cons_tads, win=win, all_bnd_calls=all_bnd_calls){

        #get boundary gap scores for defined window
        colname=paste0("Gap_Score.",win)
        
        names(bnds_gapscores)[names(bnds_gapscores) == colname]="Gap_Score"
        bnds_gapscores=bnds_gapscores|>dplyr::filter(abs(Gap_Score)>1.8)

        cons_tads.dif=cons_tads|>
                dplyr::filter(bnd_id_5%in%bnds_gapscores$bnd_id | bnd_id_3%in%bnds_gapscores$bnd_id)|>
                dplyr::mutate(dif_bnd=case_when(bnd_id_5%in%bnds_gapscores$bnd_id & !bnd_id_3%in%bnds_gapscores$bnd_id ~ "bnd_5",
                        !bnd_id_5%in%bnds_gapscores$bnd_id & bnd_id_3%in%bnds_gapscores$bnd_id ~ "bnd_3",
                        bnd_id_5%in%bnds_gapscores$bnd_id & bnd_id_3%in%bnds_gapscores$bnd_id ~ "bnd_5_3",
                 TRUE ~ "NA"))

        tadcompare_res_win_bnds=all_bnd_calls[[win]][["bnds"]]
        tadcompare_res_win_allpos=all_bnd_calls[[win]][["allpos"]]
        
        tadcompare_res_win=get_all_bnds(df_tads=cons_tads.dif,df_bnds=tadcompare_res_win_bnds,df_allsites=tadcompare_res_win_allpos)

        tadcompare_res_win=tadcompare_res_win|>dplyr::select(bnd_id,Gap_Score,Status,Type)|>
        mutate(Type=case_when(abs(Gap_Score)>1.8 & Status=="Non-Differential" ~ "other",
                                        abs(Gap_Score)>1.8 & Status=="Differential" & is.na(Type) ~ "other",
            TRUE ~ Type))|>
        mutate(Status=case_when(abs(Gap_Score)>1.8 ~ "Differential",
            TRUE ~ Status))

    # tad oriented df

   cons_tads.dif3=cons_tads.dif|>
    dplyr::filter(dif_bnd=="bnd_3" | dif_bnd=="bnd_5_3")|>
    dplyr::inner_join(tadcompare_res_win,join_by("bnd_id_3"=="bnd_id"))|>
        dplyr::mutate(win_gapscore=win)

   cons_tads.dif5=cons_tads.dif|>
    dplyr::filter(dif_bnd=="bnd_5" | dif_bnd=="bnd_5_3")|>
    dplyr::inner_join(tadcompare_res_win,join_by("bnd_id_5"=="bnd_id"))|>
        dplyr::mutate(win_gapscore=win)

  cons_tads.dif.35=rbind(cons_tads.dif3,cons_tads.dif5)

    # bnd oriented df
   difbnds_df=tadcompare_res_win|>
        dplyr::filter(bnd_id%in%bnds_gapscores$bnd_id)|>
        dplyr::mutate(win_gapscore=win)


  res_lst=list()
  res_lst[["tads"]]=as.data.frame(cons_tads.dif.35)
  res_lst[["bnds"]]=as.data.frame(difbnds_df)

  return(res_lst)
}




# stratify by TAD length and obtain gap score appropriate for the TAD length
tads_LEC.cons.flt.20=tads_LEC.cons.flt|>dplyr::filter(tad_len_win<35)
tads_BEC.cons.flt.20=tads_BEC.cons.flt|>dplyr::filter(tad_len_win<35)
tads_LEC.cons.flt.35=tads_LEC.cons.flt|>dplyr::filter(tad_len_win>=35, tad_len_win<100)
tads_BEC.cons.flt.35=tads_BEC.cons.flt|>dplyr::filter(tad_len_win>=35, tad_len_win<100)
tads_LEC.cons.flt.50=tads_LEC.cons.flt|>dplyr::filter(tad_len_win>=100)
tads_BEC.cons.flt.50=tads_BEC.cons.flt|>dplyr::filter(tad_len_win>=100)

lec_tads_win=list()
lec_tads_win[[20]]=get_difTADs(cons_tads=tads_LEC.cons.flt.20, win=20, all_bnd_calls=all_bnd_calls, bnds_gapscores=bnds_gap_scores)
lec_tads_win[[35]]=get_difTADs(cons_tads=tads_LEC.cons.flt.35, win=35, all_bnd_calls=all_bnd_calls, bnds_gapscores=bnds_gap_scores)
lec_tads_win[[50]]=get_difTADs(cons_tads=tads_LEC.cons.flt.50, win=50, all_bnd_calls=all_bnd_calls, bnds_gapscores=bnds_gap_scores)

bec_tads_win=list()
bec_tads_win[[20]]=get_difTADs(cons_tads=tads_BEC.cons.flt.20, win=20, all_bnd_calls=all_bnd_calls, bnds_gapscores=bnds_gap_scores)
bec_tads_win[[35]]=get_difTADs(cons_tads=tads_BEC.cons.flt.35, win=35, all_bnd_calls=all_bnd_calls, bnds_gapscores=bnds_gap_scores)
bec_tads_win[[50]]=get_difTADs(cons_tads=tads_BEC.cons.flt.50, win=50, all_bnd_calls=all_bnd_calls, bnds_gapscores=bnds_gap_scores)


lec_diftads=list()
lec_difbnds=list()

for (i in window_sizes){

        diftads.i=lec_tads_win[[i]][["tads"]]
        difbnds.i=lec_tads_win[[i]][["bnds"]]

        lec_diftads[[i]]=diftads.i
        lec_difbnds[[i]]=difbnds.i
}

bec_diftads=list()
bec_difbnds=list()

for (i in window_sizes){

        diftads.i=bec_tads_win[[i]][["tads"]]
        difbnds.i=bec_tads_win[[i]][["bnds"]]

        bec_diftads[[i]]=diftads.i
        bec_difbnds[[i]]=difbnds.i
}

cons_tads.dif.lec=dplyr::bind_rows(lec_diftads)
cons_tads.dif.bec=dplyr::bind_rows(bec_diftads)

## rbind to obtain full list of all dif TADs in their appropriate gap score window
cons_tads.dif.combined=rbind(cons_tads.dif.bec,cons_tads.dif.lec)


#only one df for bounds
lec_bec_difbnds.df=dplyr::bind_rows(bec_difbnds)

max_gapscores=lec_bec_difbnds.df|>
        group_by(bnd_id)|>
        dplyr::mutate(max_abs_gapscore=max(abs(Gap_Score)))|>
        rowwise()|>
        dplyr::filter(abs(Gap_Score)==max_abs_gapscore)

lec_bec_difbnds.df.final=max_gapscores|>
        dplyr::select(bnd_id,win_gapscore,max_abs_gapscore)|>
        dplyr::rename(win_max_gapscore=win_gapscore)|>
        left_join(lec_bec_difbnds.df,by="bnd_id")|>
        dplyr::filter(win_gapscore==win_max_gapscore)|>
        dplyr::mutate(Type=case_when(is.na(Type) ~ "other",
                TRUE ~ Type))




## add info on smpl identity for each tad_id
## wiggle & identical TADs

cons_tads.dif.combined.wg=cons_tads.dif.combined|>
                dplyr::mutate(smpl=case_when(tad_id%in%tads_overlap.cons.wg$tad_id_LEC | tad_id%in%tads_overlap.cons.wg$tad_id_BEC ~ "LEC_BEC",
                                tad_id%in%tads_LEC$tad_id & !tad_id%in%tads_overlap.cons.wg$tad_id_LEC ~ "LEC",
                                tad_id%in%tads_BEC$tad_id & !tad_id%in%tads_overlap.cons.wg$tad_id_BEC ~ "BEC",
                                TRUE ~ "NA"
                        ))

cons_tads.dif.combined.ident=cons_tads.dif.combined|>
                dplyr::mutate(smpl=case_when(tad_id%in%tads_overlap.cons.ident$tad_id_LEC | tad_id%in%tads_overlap.cons.ident$tad_id_BEC ~ "LEC_BEC",
                                tad_id%in%tads_LEC$tad_id & !tad_id%in%tads_overlap.cons.ident$tad_id_LEC ~ "LEC",
                                tad_id%in%tads_BEC$tad_id & !tad_id%in%tads_overlap.cons.ident$tad_id_BEC ~ "BEC",
                                TRUE ~ "NA"
                        ))
                

## save as bed

diftads.bed=cons_tads.dif.combined.ident|>
        dplyr::mutate(strand=".",score=abs(Gap_Score))|>
        dplyr::select(chr, start, end, tad_id, score, strand)

fname_str1="diffTADs"
fnamebed=paste(fname_str1,file_string1,file_pref,"bed",sep=".")

outdir_dif_TADs=file.path(resdir,"bed_dif_TADs")
dir.create(outdir_dif_TADs,recursive=TRUE)
save.bed(df=diftads.bed,dir=outdir_dif_TADs,file=fnamebed)


# save separately for LEC and BEC

diftads.lec.bed=cons_tads.dif.combined.ident|>
        dplyr::filter(smpl=="LEC" | smpl=="LEC_BEC")|>
        dplyr::mutate(strand=".",score=abs(Gap_Score))|>
    dplyr::select(chr, start, end, tad_id, score, strand)

diftads.bec.bed=cons_tads.dif.combined.ident|>
        dplyr::filter(smpl=="BEC" | smpl=="LEC_BEC")|>
        dplyr::mutate(strand=".",score=abs(Gap_Score))|>
    dplyr::select(chr, start, end, tad_id, score, strand)

fname_str1="diffTADs_LEC"
fnamebed=paste(fname_str1,file_string1,file_pref,"bed",sep=".")
save.bed(df=diftads.lec.bed,dir=outdir_dif_TADs,file=fnamebed)

fname_str1="diffTADs_BEC"
fnamebed=paste(fname_str1,file_string1,file_pref,"bed",sep=".")
save.bed(df=diftads.bec.bed,dir=outdir_dif_TADs,file=fnamebed)




## ---- diff-bnds-TADcompare


#####
### dif bnds formatting & summarising


# differential boundaries per summarising window for report table

dif_bnds_by_win=as.data.frame(table(lec_bec_difbnds.df.final$win_gapscore,lec_bec_difbnds.df.final$Type))

dif_bnds_by_win=dif_bnds_by_win|>
        dplyr::rename(sum_window=Var1,Type=Var2)|>
        tidyr::pivot_wider(names_from=c(Type),values_from=Freq)

dif_bnds_by_win$sum_window=as.character(dif_bnds_by_win$sum_window)
dif_bnds_by_win=rbind(dif_bnds_by_win,c("all",as.vector(table(lec_bec_difbnds.df.final$Type))))

all_difbounds_per_win=as.data.frame(table(lec_bec_difbnds.df.final$win_gapscore))
all_difbounds_per_win=all_difbounds_per_win|>
        dplyr::rename(sum_window=Var1)|>
        tidyr::pivot_wider(names_from=c(sum_window),values_from=Freq)
all_difbounds_per_win=as.data.frame(all_difbounds_per_win)       
all_difbounds_per_win=c(as.character(unname(all_difbounds_per_win[1,])),length(unique(lec_bec_difbnds.df.final$bnd_id)))

dif_bnds_by_win=cbind(dif_bnds_by_win,all_difbounds_per_win)
colnames(dif_bnds_by_win)[length(colnames(dif_bnds_by_win))]="all"

## bed file for visualisation

dif_bounds_all.bed=lec_bec_difbnds.df|>
        tidyr::separate(bnd_id,c("chr","start"),sep="_", remove=FALSE)|>
        dplyr::mutate(strand=".",end=as.numeric(start)+binsize,score=abs(Gap_Score))|>
        dplyr::select(chr, start, end, bnd_id, score, strand)

fname_str1="diffTAD_boundaries_LEC_BEC"
fnamebed=paste(fname_str1,file_string1,file_pref,"bed",sep=".")

outdir_dif_bnds=file.path(resdir,"bed_dif_bounds")
dir.create(outdir_dif_bnds,recursive=TRUE)
save.bed(df=dif_bounds_all.bed,dir=outdir_dif_bnds,file=fnamebed)





## ---- difTADs-euler

### OBS! something is not correct here, check later

get_overlaps_euler=function(overlap_tads=overlap_tads,dif_tads=dif_tads,tads_lec=tads_lec,tads_bec=tads_bec){

        dif_tads=dif_tads|>dplyr::distinct(tad_id,smpl)

        tads_lec=tads_lec|>dplyr::filter(!tad_id%in%overlap_tads$tad_id_LEC)|>pull(tad_id)
        tads_bec=tads_bec|>dplyr::filter(!tad_id%in%overlap_tads$tad_id_BEC)|>pull(tad_id)
        tads_lec_bec=overlap_tads$tad_id_wg

        n_all_lec=length(unique(tads_lec))
        n_all_bec=length(unique(tads_bec))
        n_all_lec_bec=length(unique(tads_lec_bec))

        diftads_lec=dif_tads|>dplyr::filter(smpl=="LEC")|>pull(tad_id)
        diftads_bec=dif_tads|>dplyr::filter(smpl=="BEC")|>pull(tad_id)
        diftads_lec_bec=dif_tads|>dplyr::filter(smpl=="LEC_BEC")|>pull(tad_id)

        n_dif_lec=length(unique(diftads_lec))
        n_dif_bec=length(unique(diftads_bec))
        n_dif_lec_bec=length(unique(diftads_lec_bec))

        n_lec_nodif=n_all_lec-n_dif_lec
        n_bec_nodif=n_all_bec-n_dif_bec
        n_lec_bec_nodif=n_all_lec_bec-n_dif_lec_bec


        set.seed(314)
        tad_all_diff_intersections = euler(c("LEC" = n_lec_nodif, "BEC" = n_bec_nodif, "diff"=0,
                                 "LEC&BEC" = n_lec_bec_nodif, "LEC&diff" = n_dif_lec, "BEC&diff" = n_dif_bec,
                                 "LEC&BEC&diff" = n_dif_lec_bec),shape = "ellipse")

  return(tad_all_diff_intersections)

}



euler_overlaps.cons.wg=get_overlaps_euler(dif_tads=cons_tads.dif.combined.wg,overlap_tads=tads_overlap.cons.wg,tads_lec=tads_LEC.cons,tads_bec=tads_BEC.cons)

euler.cons.wg=plot(euler_overlaps.cons.wg,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#FF69b450","#80808060"),
  edges=c(LEC_dark,BEC_dark,"#FF69b4","#808080"))

euler_overlaps.cons.ident=get_overlaps_euler(dif_tads=cons_tads.dif.combined.ident,overlap_tads=tads_overlap.cons.ident,tads_lec=tads_LEC.cons,tads_bec=tads_BEC.cons)

euler.cons.id=plot(euler_overlaps.cons.ident,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#FF69b450","#80808060"),
  edges=c(LEC_dark,BEC_dark,"#FF69b4","#808080"))



# save them as files

wiggle_n_bins=paste0("wiggle_",ontad_wiggle_bins)

eulerplotdir=file.path(resdir,"euler_dif_LEC_BEC_overlaps")
dir.create(eulerplotdir,recursive=TRUE)

file_string2="consensus" #or unfilt
file_string3=wiggle_n_bins
fpref=paste("euler_diff_LEC_BEC",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=euler.cons.wg,savedir=eulerplotdir, pref=fpref)

file_string2="consensus" #or unfilt
file_string3="identical"
fpref=paste("euler_diff_LEC_BEC",file_string1,file_string2,file_string3,sep=".")
save.plots(plot_ts=euler.cons.id,savedir=eulerplotdir, pref=fpref)


## ---- diff-bnds-pie

pieplotdir=file.path(resdir,"pie_dif_bnds")
dir.create(pieplotdir,recursive=TRUE)


## consensus boundaries


lec_bec_difbnds.df.final.cons=lec_bec_difbnds.df.final|>
        dplyr::filter()
lec_bec_difbnds.df.final.cons=lec_bec_difbnds.df.final|>
        dplyr::filter(bnd_id%in%cons_bnds)

tab=table(lec_bec_difbnds.df.final.cons$Type)


pie_df=as.data.frame(cbind(as.character(c(names(tab[1]),names(tab[2]),names(tab[4]),names(tab[5]),names(tab[6]),names(tab[3]))),
      as.numeric(c(unname(tab[1]),unname(tab[2]),unname(tab[4]),unname(tab[5]),unname(tab[6]),unname(tab[3])))))
colnames(pie_df)=c("type","count")
pie_df$count=as.numeric(pie_df$count)

# Compute the position of labels
pie_df=pie_df |>
  arrange(desc(type)) |>
  mutate(prop = count / sum(pie_df$count) *100) |>
  mutate(ypos = cumsum(prop)- 0.5*prop )

pie_cons=ggplot(pie_df|>arrange(desc(type)), aes(x="", y=prop, fill=type))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  theme_void()+
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#99999940","#8A2BE2","#008080","#BC8F8F"))+
  geom_text(aes(y = ypos, label = count), color = "white", size=4)

fname=paste0(file_string1,".difTAD_pie_reproducible")
for (devplot in c("pdf","eps")){
        fname_full=paste(fname,devplot,sep=".")
        ggsave(plot=pie_cons, filename=fname_full, path=pieplotdir, device=devplot)
}

pie_df



## ---- section-dif-bnds-annot

########################
########################
###########     Dif Boundaries TSS annotation
########################


## ---- annot-difbnds-TSS


# change IDs of the flanking genes
swap_gene_ids<-function(x){ 

    ensids=unlist(strsplit(x, split=";"))

    #rm NAs > substitute with empty strings
    ensids[is.na(ensids)]<-""
    ens2symbol=mapIds(drer11.ensDb,keys=ensids,column="SYMBOL", keytype="GENEID", multiVals="first")

    flank_geneIDs=data.frame(SYMBOL=ens2symbol,ensembl_gene_id=ensids)
    flank_geneIDs=flank_geneIDs|>mutate_if(is.character, list(~na_if(.,"")))%>%
        mutate(SYMBOL = coalesce(SYMBOL, ensembl_gene_id))

    flank_geneIDs.str=paste(flank_geneIDs$SYMBOL,collapse="; ")

    return(flank_geneIDs.str)
}

# change the delimiter by adding space for distances - to make tables nicer on the eye
change_str_delim<-function(x){
    str_sep=unlist(strsplit(x, split=";"))
    str_modified=paste(str_sep,collapse="; ")
    return(str_modified)
}


lec_bec_difbnds.df.2=lec_bec_difbnds.df.final|>
        tidyr::separate(bnd_id,c("chr","start"),sep="_", remove=FALSE)|>
        dplyr::mutate(start = as.numeric(start))|>
        dplyr::mutate(end=start+binsize, .after="start")|>
        dplyr::select(chr,start,end,bnd_id,win_gapscore,Gap_Score,Type)


dif_boundary.gr=regioneR::toGRanges(as.data.frame(lec_bec_difbnds.df.2))

dif_boundary.annot=as.data.frame(ChIPseeker::annotatePeak(dif_boundary.gr, tssRegion=c(-3000, 3000),TxDb=txdb_ens, annoDb="org.Dr.eg.db", level="gene", addFlankGeneInfo = TRUE, flankDistance=50000))


#for table
dif_boundary.annot$gene_symbol=unlist(lapply(dif_boundary.annot$flank_geneIds,FUN=swap_gene_ids))
dif_boundary.annot$flank_gene_distances=unlist(lapply(dif_boundary.annot$flank_gene_distances,FUN=change_str_delim))

dif_boundary_all.annot=dif_boundary.annot|>
        dplyr::rename(start_bnd=start,end_bnd=end,chr_bnd=seqnames)|>
        dplyr::select(!c("width","strand","geneChr"))|>
    dplyr::distinct(bnd_id, .keep_all=TRUE)

fname_str1="Differential_TAD_boundaries_ontad_all_lvls.LEC_vs_BEC"
fname=paste(fname_str1,file_string1,file_pref,"tsv",sep=".")
save.table(df=dif_boundary_all.annot, file=fname,dir=outdir_dif_bnds)

dif_boundary_all.annot.tbl=dif_boundary_all.annot|>
        dplyr::select(bnd_id,Gap_Score,Type,gene_symbol,flank_gene_distances)|>
        dplyr::arrange(bnd_id)|>
        dplyr::rename(Boundary_ID=bnd_id)|>
        dplyr::distinct(Boundary_ID, .keep_all=TRUE)



## ---- difboundaries-overlaps-fcn

#functions

#overlap of features (body) with dif TADs

feature_body_overlap=function(test_tads.gr=test_tads.gr,feature.gr=feature.gr){

  feature.test_tads=regioneR::overlapRegions(feature.gr,test_tads.gr, colA=c("id"),colB=colnames(mcols(test_tads.gr)))
  colnames(feature.test_tads)[c(6:ncol(feature.test_tads))]=c("id",colnames(mcols(test_tads.gr)),"overlap_type")

  feature.test_tads.red_gene=feature.test_tads|> 
    group_by(id)|>
    dplyr::mutate(tad_ids = paste0(tad_id,collapse=","))|> 
        distinct(id, tad_ids, .keep_all = TRUE)


  feature.test_tads.overlaps=list()
  feature.test_tads.overlaps[["complete"]]=feature.test_tads
  feature.test_tads.overlaps[["genes"]]=as.data.frame(feature.test_tads.red_gene)


  return(feature.test_tads.overlaps)

}

## test the significance of overlap of dif TADs with a selected gene list vs all genes

# gene_sel.ls list of selected features, subset of genes.all
# genes.all
# overlaps_dif and overlaps.all are output of feature_body_overlap[["genes"]]

# > colnames(overlaps.dif)
#  [1] "chr"          "startA"       "endA"         "startB"       "endB"         "id"           "tad_id"       "bnd_id_5"    
#  [9] "bnd_id_3"     "dif_bnd"      "smpl"         "overlap_type" "tad_ids"     

# genes.all=genes_all$external_gene_name
# genes.sel=genes_lec$external_gene_name
# overlaps.dif=gene_annot_difTAD[["genes"]]
# overlaps.all=gene_annot_consTAD[["genes"]]


get_overlaps_sig=function(genes.sel=genes.sel,genes.all=genes.all,overlaps.dif=overlaps.dif,overlaps.all=overlaps.all){


  genes_all.difTAD.ls=unique(overlaps.dif$id)
  genes_sel.difTAD.ls=genes_all.difTAD.ls[genes_all.difTAD.ls%in%genes.sel]

  genes_all.allTAD.ls=unique(overlaps.all$id)
  genes_sel.allTAD.ls=genes_all.allTAD.ls[genes_all.allTAD.ls%in%genes.sel]


  n_genes.TAD=length(genes_all.allTAD.ls)
  n_genes.TAD_sel=length(genes_all.difTAD.ls)
  n_genes.TAD_nosel=n_genes.TAD-n_genes.TAD_sel

  n_genes_sel.TAD=length(genes_sel.allTAD.ls)  
  n_genes_sel.TAD_sel=length(genes_sel.difTAD.ls)
  n_genes_sel.TAD_nosel=n_genes_sel.TAD-n_genes_sel.TAD_sel

  n_genes_nosel.TAD=n_genes.TAD-n_genes_sel.TAD
  n_genes_nosel.TAD_sel=n_genes.TAD_sel-n_genes_sel.TAD_sel
  n_genes_nosel.TAD_nosel=n_genes.TAD-n_genes.TAD_sel-n_genes_sel.TAD

  n_gens_nosel.TAD=n_genes.TAD-n_genes_sel.TAD

  contingency_table_feat_tad=rbind(c("TADs / Features",  "diff features",  "no-diff features",         "Sum"),
                              c("diff TADs",         n_genes_sel.TAD_sel,   n_genes_nosel.TAD_sel,    n_genes.TAD_sel ),
                              c("no-diff TADs",      n_genes_sel.TAD_nosel, n_genes_nosel.TAD_nosel,  n_genes.TAD_nosel) ,
                              c("Sum",               n_genes_sel.TAD,       n_genes_nosel.TAD,        n_genes.TAD) )



  mat.fisher=matrix(c(n_genes.TAD - (n_genes_sel.TAD+n_genes.TAD_sel),n_genes_sel.TAD_nosel,
    n_genes_nosel.TAD_sel, n_genes_sel.TAD_sel), nrow=2) 

  p.fisher=fisher.test(mat.fisher, alternative="greater") 

  list_res=list()
  list_res[["p.fisher"]]=p.fisher
  list_res[["cont_tab"]]=contingency_table_feat_tad

  return(list_res)

}

#############
#region association by permutation test

## chr lengths for permutation test
lengths = seqlengths(BSgenome.Drerio.UCSC.danRer11)
lengths_df=as.data.frame(lengths)
lengths_df$chr=rownames(lengths_df)
lengths_df=lengths_df|>dplyr::filter(grepl('chr', chr))|>dplyr::filter(!grepl('Un|alt', chr))
lengths_df$chr=gsub("chr","",lengths_df$chr)
lengths_df=lengths_df|>dplyr::rename(end=lengths)|>dplyr::mutate(start=1)|>
  dplyr::select(chr,start,end)



## ---- difbnds-overlaps-ctcf

## permutation test
# consensus bounds
bnds_LEC.cons.5=tads_LEC.cons|>
      dplyr::select(chr,start,bnd_id_5,bnd_5_level)|>
      dplyr::mutate(end=start+binsize, .after="start")|>
      dplyr::rename(bnd_id=bnd_id_5,bnd_level=bnd_5_level)
bnds_LEC.cons.3=tads_LEC.cons|>
      dplyr::select(chr,start,bnd_id_3,bnd_3_level)|>
      dplyr::mutate(end=start+binsize, .after="start")|>
      dplyr::rename(bnd_id=bnd_id_3,bnd_level=bnd_3_level)

bnds_BEC.cons.5=tads_BEC.cons|>
      dplyr::select(chr,start,bnd_id_5,bnd_5_level)|>
      dplyr::mutate(end=start+binsize, .after="start")|>
      dplyr::rename(bnd_id=bnd_id_5,bnd_level=bnd_5_level)
bnds_BEC.cons.3=tads_BEC.cons|>
      dplyr::select(chr,start,bnd_id_3,bnd_3_level)|>
      dplyr::mutate(end=start+binsize, .after="start")|>
      dplyr::rename(bnd_id=bnd_id_3,bnd_level=bnd_3_level)

all_cons_bounds.gr=rbind(bnds_LEC.cons.5,bnds_LEC.cons.3,bnds_BEC.cons.5,bnds_BEC.cons.3)
all_cons_bounds.gr=all_cons_bounds.gr|>
  dplyr::distinct(bnd_id, .keep_all=TRUE)

all_cons_bounds.gr=regioneR::toGRanges(as.data.frame(all_cons_bounds.gr))

# dif bounds

difbounds.gr=lec_bec_difbnds.df.final|>
          tidyr::separate(bnd_id,c("chr","start"),sep="_", remove=FALSE)|>
          dplyr::mutate_at(vars(chr,start), as.numeric)|>
          dplyr::mutate(end=start+binsize,.after="start")|>
          dplyr::select(chr,start,end,bnd_id)

difbounds.gr=regioneR::toGRanges(as.data.frame(difbounds.gr))

nodif_bounds.gr=setdiff(all_cons_bounds.gr, difbounds.gr, ignore.strand=TRUE) 
nodif_bounds_df=as.data.frame(nodif_bounds.gr)
names=paste(nodif_bounds_df$seqnames,nodif_bounds_df$start,sep="_")
mcols(nodif_bounds.gr)$bnd_id=names


## use regioneReloded
all_cons_bounds.gr.1=all_cons_bounds.gr|>plyranges::filter(bnd_level==1)
all_cons_bounds.gr.2_3=all_cons_bounds.gr|>plyranges::filter(bnd_level==2 | bnd_level==3)
all_cons_bounds.gr.4=all_cons_bounds.gr|>plyranges::filter(bnd_level>=4)

difbounds.gr.1=difbounds.gr|>plyranges::filter(bnd_id%in%all_cons_bounds.gr.1$bnd_id)
difbounds.gr.2_3=difbounds.gr|>plyranges::filter(bnd_id%in%all_cons_bounds.gr.2_3$bnd_id)
difbounds.gr.4=difbounds.gr|>plyranges::filter(bnd_id%in%all_cons_bounds.gr.4$bnd_id)


ctcf_bnds.lst=list()
ctcf_bnds.lst[["ctctf-48h"]]=ctcf_wt_48.gr
ctcf_bnds.lst[["cons-1"]]=all_cons_bounds.gr.1
ctcf_bnds.lst[["cons-2_3"]]=all_cons_bounds.gr.2_3
ctcf_bnds.lst[["cons-4"]]=all_cons_bounds.gr.4
ctcf_bnds.lst[["diff-1"]]=difbounds.gr.1
ctcf_bnds.lst[["diff-2_3"]]=difbounds.gr.2_3
ctcf_bnds.lst[["diff-4"]]=difbounds.gr.4
ctcf_bnds.lst[["cons_bnds"]]=all_cons_bounds.gr
ctcf_bnds.lst[["dif_bnds"]]=difbounds.gr
ctcf_bnds.lst[["nodif_bnds"]]=nodif_bounds.gr



set.seed(2357)
bnds_ctcf.perm=crosswisePermTest(Alist=ctcf_bnds.lst, Blist=ctcf_bnds.lst, sampling=FALSE,
    gen=lengths_df,ntimes=n_perms,alternative="greater",force.parallel=TRUE, mc.set.seed=FALSE)

bnds.ctcf48_perm_df=bnds_ctcf.perm@multiOverlaps
bnds.ctcf48_perm_res_full=bind_rows(bnds.ctcf48_perm_df, .id="setA")|>
  dplyr::select(!order.id)|>
  dplyr::rename(setB=name)|>
  dplyr::filter(!setA==setB)|>
  dplyr::select(setA,setB,p_value,adj.p_value,n_regionA,n_regionB,n_overlaps,z_score,norm_zscore,mean_perm_test,sd_perm_test)

pertest_ctcf_dir=file.path(resdir,"permutation_tests_difBnds_CTCF")
dir.create(pertest_ctcf_dir,recursive=TRUE)

fname_pref_bnds_ctcf=paste0("boundaries_CTCF_48.permtest_regioneReloaded_nperm_",n_perms)

fname_perm_ctcf=paste0(fname_pref_bnds_ctcf,".tsv")
save.table(df=bnds.ctcf48_perm_res_full,file=fname_perm_ctcf,dir=pertest_ctcf_dir)

fname_perm_ctcf_rds=paste0(fname_pref_bnds_ctcf,".rds")
fpath_rds=file.path(pertest_ctcf_dir,fname_perm_ctcf_rds)
saveRDS(bnds_ctcf.perm, file=fpath_rds)



## report table

bnds.ctcf48_perm_res_tbl=bnds.ctcf48_perm_res_full|>
  dplyr::select(setA,setB,n_regionA,n_regionB,n_overlaps,p_value,adj.p_value,z_score,norm_zscore)|>
  dplyr::mutate_at(vars(p_value,adj.p_value,norm_zscore,z_score), as.numeric)|>
  dplyr::filter(!setA==setB)|>
  dplyr::filter(setB=="ctctf-48h")|>
  dplyr::mutate_if(is.numeric, format, digits=4,nsmall = 4)|>
  dplyr::rename(Boundaries=setA,CTCF_48=setB,n_boundaries=n_regionA,n_CTCF_48_peaks=n_regionB)


## plot

bnds.ctcf48_perm_res_tbl.sel=bnds.ctcf48_perm_res_tbl|>
        dplyr::mutate_at(vars(p_value,adj.p_value,norm_zscore), as.numeric)|>
        dplyr::mutate(logPval=-log10(adj.p_value))|>
        dplyr::filter(!Boundaries=="nodif_bnds")|>
        dplyr::mutate(Boundary_level=case_when(
                str_detect(Boundaries, regex("-1", ignore_case=TRUE)) ~ "1",
                str_detect(Boundaries, regex("2_3", ignore_case=TRUE)) ~ "2_3",
                str_detect(Boundaries, regex("-4", ignore_case=TRUE)) ~ "4",
                TRUE ~ "all"))|>
        dplyr::mutate(Status=case_when(
                str_detect(Boundaries, regex("cons", ignore_case=TRUE)) ~ "reproducible",
                str_detect(Boundaries, regex("dif", ignore_case=TRUE)) ~ "differential",
                TRUE ~ "na"))

#https://waldyrious.net/viridis-palette-generator/
#scale_col4=c("#440154","#365c8d","#4ac16d","#fde725")
scale_col4=c("#440154","#2a788e","#4ac16d","#fde725")

bnds_ctcf48.pl=ggplot(bnds.ctcf48_perm_res_tbl.sel,aes(x=norm_zscore,y=logPval, color=Boundary_level,shape=Status)) +
        geom_point(size=3)+
        scale_color_manual(values=scale_col4)+
        scale_fill_manual(values=scale_col4)+
        scale_shape_manual(values = c(15,16))+ # c(17,2)
        geom_hline(yintercept=-log10(0.05),color="black",linetype=2)+
        theme_bw()+
        theme(aspect.ratio = 1/1.617)+
        xlab("Normalised Z score\n(permutation test)") +
        ylab("-log10 (p adjusted)")

fname=paste0("summary_plot.",fname_pref_bnds_ctcf)
for (devplot in c("pdf","eps")){
        fname_full=paste(fname,devplot,sep=".")
        ggsave(plot=bnds_ctcf48.pl, filename=fname_full, path=pertest_ctcf_dir, device=devplot)
}






## ---- diffTADs-overlaps-fcn


get_seltad_overlaps <- function(sel_tads.gr=sel_tads.gr,goi_genes_df=goi_genes_df){

  goi_genes_df_sel=goi_genes_df|>dplyr::select(chromosome_name,start_position,end_position,external_gene_name)|>
        dplyr::rename(chr=chromosome_name,start=start_position,end=end_position,gene=external_gene_name)
  goi_genes.gr=regioneR::toGRanges(goi_genes_df_sel)

  goi_genes.sel_tads=regioneR::overlapRegions(goi_genes.gr,sel_tads.gr, colA=c("gene"),colB=colnames(mcols(sel_tads.gr)))
  colnames(goi_genes.sel_tads)[c(6:ncol(goi_genes.sel_tads))]=c("gene_name",colnames(mcols(sel_tads.gr)),"overlap_type")

  goi_genes.sel_tads.red=goi_genes.sel_tads|>
    dplyr::group_by(gene_name,chr,startA,endA )|>
    dplyr::summarise(tad_ids=paste0(tad_id,collapse=","), diff_status=paste0(dif_bnd,collapse=",") )|>
    dplyr::rename(start_gene=startA,end_gene=endA)


  df_list=list()
  df_list[[1]]=as.data.frame(goi_genes.sel_tads.red)
  df_list[[2]]=as.data.frame(goi_genes.sel_tads)

  return(df_list)
}

#save bed for region plotting & gene ids for reference
# df=df_list[[2]]

save_overlap_genes<-function(df=df,dir=dir,pref1=pref1,pref2=pref2,cons_tads=cons_tads){
    
    cons_tads=cons_tads|>dplyr::select(tad_id,Type,Gap_Score)|>
        dplyr::distinct(tad_id,Type,Gap_Score)

    df_file=df|>
    dplyr::left_join(cons_tads,by="tad_id")|>
    group_by(gene_name)|>
    dplyr::distinct(tad_id)
    
    fname=paste("gene_TAD_overlap",pref1,pref2,"tsv",sep=".")
    fname_bed=paste("gene_TAD_overlap",pref1,pref2,"bed",sep=".")

    save.table(df=df_file,file=fname,dir=dir)

    # bed file
    bed.genes=df|>dplyr::mutate(strand=".")|>
        dplyr::mutate(score=".")|>
        dplyr::select(chr,startA,endA,gene_name,score,strand)

    
    save.bed(df=bed.genes,dir=dir,file=fname_bed)

}



format_goi_tads_fortable<-function(overlaps=overlaps){
  overlaps=overlaps|>
    dplyr::rename(
            start_gene=startA,end_gene=endA) |>
    dplyr::arrange(gene_name)

  return(as.data.frame(overlaps))

}




## ---- diffTADs-foroverlaps

txt_string_seltads="differential TADs"

pref_seltads="difTADs"

subresdir_name=paste0("overlaps-",pref_seltads)

subresdir=file.path(resdir,subresdir_name)
dir.create(subresdir,recursive=TRUE)

cons_tads.dif.combined.wg.uniq=cons_tads.dif.combined.wg|>dplyr::distinct(tad_id,smpl, .keep_all=TRUE)

# > colnames(cons_tads.dif.combined.wg)
#  [1] "TAD_level"    "TADmean"      "TADscore"     "chr"          "start"        "end"          "tad_id"       "bnd_id_5"    
#  [9] "bnd_id_3"     "smpl"         "bnd_5_level"  "bnd_3_level"  "tad_len_win"  "dif_bnd"      "Gap_Score"    "Status"      
# [17] "Type"         "win_gapscore"

dif_tads.gr=cons_tads.dif.combined.wg.uniq|>dplyr::select(chr,start,end,tad_id,bnd_id_5,bnd_id_3,dif_bnd,smpl)
dif_tads.gr=regioneR::toGRanges(as.data.frame(dif_tads.gr))


### for overlap significance

#all consensus TADs
#the original tad_id is subbed to tad_id_wg for LEC_BEC tads to avoid feature number inflation
tads_LEC.cons.df=tads_LEC.cons|>
        dplyr::select(chr,start,end,tad_id,smpl,TAD_level)|>
        left_join(tads_overlap.cons.wg|>dplyr::select(!chr),join_by("tad_id"=="tad_id_LEC"))|>
        dplyr::mutate(smpl=case_when(tad_id%in%tads_overlap.cons.wg$tad_id_LEC ~ "LEC_BEC",
                                !tad_id%in%tads_overlap.cons.wg$tad_id_LEC ~ "LEC",
                                TRUE ~ "NA"))|>
        dplyr::mutate(tad_id_gr=case_when(!is.na(tad_id_wg) ~ tad_id_wg,
                                TRUE ~ tad_id))|>
        dplyr::select(!c(startA,startB,endA,endB,tad_id_BEC,tad_id_wg,type))


tads_BEC.cons.df=tads_BEC.cons|>
        dplyr::select(chr,start,end,tad_id,smpl,TAD_level)|>
        left_join(tads_overlap.cons.wg|>dplyr::select(!chr),join_by("tad_id"=="tad_id_BEC"))|>
        dplyr::mutate(smpl=case_when(tad_id%in%tads_overlap.cons.wg$tad_id_BEC ~ "LEC_BEC",
                                !tad_id%in%tads_overlap.cons.wg$tad_id_BEC ~ "BEC",
                                TRUE ~ "NA"))|>
        dplyr::mutate(tad_id_gr=case_when(!is.na(tad_id_wg) ~ tad_id_wg,
                                TRUE ~ tad_id))|>
        dplyr::select(!c(startA,startB,endA,endB,tad_id_LEC,tad_id_wg,type))

tads_all.cons.df=rbind(tads_LEC.cons.df,tads_BEC.cons.df)
tads_all.cons.df=tads_all.cons.df|>
        dplyr::select(!tad_id)|>
        dplyr::distinct(tad_id_gr,smpl,.keep_all=TRUE)|>
        dplyr::rename(tad_id=tad_id_gr)

all_consTADs.gr=regioneR::toGRanges(as.data.frame(tads_all.cons.df))


# for child
dif_bounds.gr=dif_boundary.gr
sel_tads.gr=dif_tads.gr



## ---- foo





