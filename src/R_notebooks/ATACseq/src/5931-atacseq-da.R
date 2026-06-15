
## ---- prep-env

library(renv)

## renv::restore()
renv::activate()


library(edgeR)
library(GenomicAlignments)
library(GenomicFeatures)
library(ChIPseeker)
library(Hmisc)
library(sva)

library(org.Dr.eg.db)

library(tidyverse)
library(tibble)
library(dplyr)
library(DT)
library(htmltools)
library(magrittr)

library(plotly)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggupset)
library(ggimage)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(wesanderson)
library(circlize)
library(dendextend)
library(eulerr)

library(knitr)
library(kableExtra)

#renv::snapshot()


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


colour_scale_hm_div=circlize::colorRamp2(c(-2, 0, 2), c("#00868B", "#FFFAF0",  "#CD0000"))

colour_scale_hm=viridis(100)

ann_colors.1=list(
    tissue=c(BEC="#762A83",LEC="#109630")
    )



## ---- cutoffs

FDR_CO=0.05
lFC_CO=0.5
#lFC_CO=0



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

resdir=file.path(projrootdir,"results",paste0("ATACseq_",file_pref))
dir.create(resdir, recursive=TRUE)


## ---- annots

## annotations
annotdir=file.path(projrootdir,"reference","annotation")

ff=FaFile(file.path(projrootdir,"reference","fasta","danRer11noALT.fa"))

ensembl_annot_trxdb=file.path(annotdir,"Ensembl.txdb.Drer11.Rdata")
txdb_ens=loadDb(ensembl_annot_trxdb)
#append chr to sequence names
newSeqNames=paste('chr', seqlevels(txdb_ens), sep = '')
names(newSeqNames)=seqlevels(txdb_ens)
txdb_ens=renameSeqlevels(txdb_ens, newSeqNames )

# gene names from biomart
fname="drer11_gene_names.tab"
gene_annot_drer11=read.delim(file.path(annotdir,fname), sep="\t", header=TRUE, quote = "")



## ---- data

data.pref.txt="Consensus peaks within replicates"
dat_prefix="consensus_intrareplicate_merged"
DAstring="DA_groups" # file naming

res_subdir=file.path(resdir,dat_prefix)
dir.create(res_subdir, recursive=TRUE)


atac_datadir=file.path(datadir,"ATACseq")

cnt_table.pth=file.path(atac_datadir,"consensus_merged_macs2broad.tsv")
cnt_table=read.table(cnt_table.pth, sep="\t", header=TRUE, blank.lines.skip=TRUE)
rownames(cnt_table)=cnt_table$Geneid
colnames(cnt_table)=c(colnames(cnt_table)[1:6],gsub(".filtered.sorted.bam","",colnames(cnt_table)[7:12]))
colnames(cnt_table)=c(colnames(cnt_table)[1:6],gsub("P17204_","",colnames(cnt_table)[7:12]))

cnt_table=cnt_table[!cnt_table$Chr=="chrM",]

groups=factor(c(rep("LEC",3),rep("BEC",3)))


# get GC content
gr=GRanges(seqnames=cnt_table$Chr, ranges=IRanges(cnt_table$Start, cnt_table$End), strand="*", mcols=data.frame(peakID=cnt_table$Geneid))
peakSeqs=getSeq(x=ff, gr)
gcContentPeaks=letterFrequency(peakSeqs, "GC",as.prob=TRUE)[,1]
gcGroups=Hmisc::cut2(gcContentPeaks, g=20)
mcols(gr)$gc=gcContentPeaks



## ---- peaks_profile

#profile around TSS
promoter=getPromoters(TxDb=txdb_ens, upstream=3000, downstream=3000)
tagMatrix=getTagMatrix(gr, windows=promoter)

TSS_profile=plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")



## ---- peak_annotation

fname_annot_closest=paste("ATAC_peaks_annot_closest.ensembl",DAstring,file_pref,"tsv", sep=".")

peakAnno=annotatePeak(gr, tssRegion=c(-3000, 3000),TxDb=txdb_ens, annoDb="org.Dr.eg.db")
peakAnnoplot=upsetplot(peakAnno, vennpie=TRUE)


peakAnno.closest=as.data.frame(peakAnno)
colnames(peakAnno.closest)=c(colnames(peakAnno.closest)[1:7],paste(colnames(peakAnno.closest)[8:19],"closest",sep="."))

peakAnno.closest$PeakID=peakAnno.closest$mcols.peakID
rm_cols=!names(peakAnno.closest)%in%c("mcols.peakID")
peakAnno.closest=peakAnno.closest[,rm_cols ,drop = FALSE]


fname_annot_flanking=paste("ATAC_peaks_annot_flanking.ensembl",DAstring,file_pref,"tsv", sep=".")

#only 3´genes
peakAnno.3=as.data.frame(annotatePeak(gr, tssRegion=c(-3000, 3000),TxDb=txdb_ens, annoDb="org.Dr.eg.db", ignoreUpstream=TRUE,ignoreOverlap=TRUE))
colnames(peakAnno.3)=c(colnames(peakAnno.3)[1:7],paste(colnames(peakAnno.3)[8:19],"genes3",sep="."))

#only 5´genes
#gr object has 1st 2 rows removed (the 5' most peak on chr1 does not have any 5' annotation) to avoid:  Unknown ID type, gene annotation will not be added...
peakAnno.5=as.data.frame(annotatePeak(gr[-c(1,2)], tssRegion=c(-3000, 3000),TxDb=txdb_ens, annoDb="org.Dr.eg.db", ignoreDownstream=TRUE,ignoreOverlap=TRUE))
colnames(peakAnno.5)=c(colnames(peakAnno.5)[1:7],paste(colnames(peakAnno.5)[8:19],"genes5",sep="."))
peakAnno.5=peakAnno.5[,c(6,8:19)]

# closest incl overlaps
peakAnno.overlap=as.data.frame(annotatePeak(gr, tssRegion=c(-3000, 3000),TxDb=txdb_ens, annoDb="org.Dr.eg.db", ignoreOverlap=FALSE))
colnames(peakAnno.overlap)=c(colnames(peakAnno.overlap)[1:7],paste(colnames(peakAnno.overlap)[8:19],"genesOverlap",sep="."))

# replace values with NA when no overlap found
peakAnno.overlap=peakAnno.overlap[,c(6,8:19)]
peakAnno.overlap_mod=peakAnno.overlap|>
	mutate(distanceToTSS.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", distanceToTSS.genesOverlap))|>
 	mutate(transcriptId.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", transcriptId.genesOverlap))|>
 	mutate(annotation.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", annotation.genesOverlap))|>
 	mutate(geneChr.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", geneChr.genesOverlap))|>
 	mutate(geneStart.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", geneStart.genesOverlap))|>
 	mutate(geneEnd.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", geneEnd.genesOverlap))|>
 	mutate(geneLength.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", geneLength.genesOverlap))|>
 	mutate(geneStrand.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", geneStrand.genesOverlap))|>
 	mutate(geneId.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", geneId.genesOverlap))|>
 	mutate(transcriptId.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", transcriptId.genesOverlap))|>
 	mutate(distanceToTSS.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", distanceToTSS.genesOverlap))|>
 	mutate(ENTREZID.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", ENTREZID.genesOverlap))|>
 	mutate(SYMBOL.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", SYMBOL.genesOverlap))|>
 	mutate(GENENAME.genesOverlap=ifelse(distanceToTSS.genesOverlap !=0, "NA", GENENAME.genesOverlap))


#merge and clean up the data frame
all_annots=list(peakAnno.3,peakAnno.5,peakAnno.overlap_mod)
all_annots_df=all_annots|>reduce(full_join, by='mcols.peakID')

all_annots_df$PeakID=all_annots_df$mcols.peakID
rm_cols=!names(all_annots_df)%in%c("mcols.peakID")
all_annots_df=all_annots_df[,rm_cols ,drop = FALSE]




## ---- TMM-norm


reads.peak=cnt_table[,c(7:12)]

groups=factor(c(rep("LEC",3),rep("BEC",3)))
sorting=c(1,2,3,1,2,3)
sorting=factor(sorting)
design=model.matrix(~sorting+groups)
rownames(design)=colnames(reads.peak)

design


reads.dge=DGEList(counts=reads.peak, group=groups)
reads.dge=calcNormFactors(reads.dge)
keep=filterByExpr(reads.dge)

summary(keep)

reads.dge=reads.dge[keep,,keep.lib.sizes=FALSE]

reads.TMM=cpm(reads.dge, normalized.lib.sizes = TRUE,log = FALSE) # for plots
logreads.TMM=cpm(reads.dge, normalized.lib.sizes = TRUE,log = TRUE) # for plots

#subset counts
count_table.sub=cnt_table[keep,]

# remove batch effect from sorting for plotting
metadata=data.frame(cbind(as.character(groups),sorting))
colnames(metadata)=c("groups","sorting")

design.combat=model.matrix(~groups)
reads.TMM.combat=ComBat(reads.TMM, metadata$sorting, mod = design.combat)
logreads.TMM.combat=ComBat(logreads.TMM, metadata$sorting, mod = design.combat)


#subset GRanges object for logFC binning
gr.sub=gr[keep,]
gcGroups.sub=gcGroups[keep]

# for child documents
reads.lognorm_plots=reads.TMM.combat
report_res_phrase="paired"


## ---- pca

reads_pca=logreads.TMM.combat

pca=prcomp(t(reads_pca),center = TRUE)

df_pca=as.data.frame(pca$x)
df_pca$library=rownames(df_pca)

pca=ggplot(df_pca, aes(x=PC1,y=PC2, colour=library)) +
	geom_point(size=6) +
	theme(text = element_text(size = 6)) +
  	theme_bw() + theme(aspect.ratio = 1)+
  	scale_color_manual(values=wesanderson::wes_palette("Zissou1", 6, "continuous"))+
  	geom_text_repel(size=4, label=df_pca$library) 

## ---- da

###### DA using TMM

reads.dge=estimateDisp(reads.dge, design)
fit=glmFit(reads.dge, design)
lrt.tmm=glmLRT(fit)

DA_res=as.data.frame(topTags(lrt.tmm, nrow(lrt.tmm$table)))

DA.res.sig.dn=DA_res|>
		dplyr::filter(FDR<FDR_CO & abs(logFC)>=lFC_CO)|>
        dplyr::filter(logFC<0)

DA.res.sig.up=DA_res|>
		dplyr::filter(FDR<FDR_CO & abs(logFC)>=lFC_CO)|>
        dplyr::filter(logFC>0)


DA_res$Geneid=rownames(DA_res)
DA.res.coords=DA_res|>
	left_join(count_table.sub[1:4],by="Geneid")|> #Geneid is PeakID here
	dplyr::rename(PeakID=Geneid)|>
	arrange(FDR)

#merge with annotations
DA_res_annot=left_join(DA.res.coords,all_annots_df,by="PeakID")
DA_res_annot=DA_res_annot|>arrange(FDR)

DA_res_annot_closest=left_join(DA.res.coords,peakAnno.closest,by="PeakID")
DA_res_annot_closest=DA_res_annot_closest|>arrange(FDR)


sub_resdir=file.path(res_subdir,"DA_all_peaks")
dir.create(sub_resdir,recursive=TRUE)

save.table(df=DA_res_annot, file=fname_annot_flanking,dir=sub_resdir)
save.table(df=DA_res_annot_closest, file=fname_annot_closest,dir=sub_resdir)



## ---- peaks_venn

# using eulerr package for consistency with other reports

# three lists
get_overlaps_euler_list_2=function(ls_lec=ls_lec,ls_bec=ls_bec,ls_all=ls_all){

  n_1=length(!is.na(unique(ls_lec)))
  n_2=length(!is.na(unique(ls_bec)))

  ls_overlap=setdiff(ls_all,c(ls_lec,ls_bec))

  n_1_2=length(unique(ls_overlap))

  set.seed(314)
  list_intersections = euler(c("LEC" = n_1, "BEC" = n_2,
                                 "LEC&BEC" = n_1_2), shape = "ellipse")

  return(list_intersections)

}

venn_DA_peaks.int=get_overlaps_euler_list_2(ls_lec=rownames(DA.res.sig.up),ls_bec=rownames(DA.res.sig.dn),ls_all=rownames(DA_res))

venn_DA_peaks=plot(venn_DA_peaks.int,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))


fname_DA_venn=paste("ATAC_peaks_Venn",DAstring,file_pref, sep=".")
save.plots(plot_ts=venn_DA_peaks,savedir=sub_resdir,pref=fname_DA_venn)


# no size cutoff

DA.res.sig.dn.0=DA_res|>
		dplyr::filter(FDR<FDR_CO)|>
        dplyr::filter(logFC<0)

DA.res.sig.up.0=DA_res|>
		dplyr::filter(FDR<FDR_CO)|>
        dplyr::filter(logFC>0)


venn_DA_peaks.int.0=get_overlaps_euler_list_2(ls_lec=rownames(DA.res.sig.up.0),ls_bec=rownames(DA.res.sig.dn.0),ls_all=rownames(DA_res))

venn_DA_peaks.0=plot(venn_DA_peaks.int.0,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))


fname_DA_venn=paste("ATAC_peaks_Venn","nolog2FC",DAstring,file_pref, sep=".")
save.plots(plot_ts=venn_DA_peaks.0,savedir=sub_resdir,pref=fname_DA_venn)


## ---- peaks_promoter

promoter_annot=c("Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)")

fname_annot_closest_prom=paste("ATAC_peaks_annot_closest_promoter.ensembl",DAstring,file_pref,"tsv", sep=".")


DA_res_annot_promoter=rbind(DA_res_annot[DA_res_annot$annotation.genes3%in%promoter_annot,], 
	DA_res_annot[DA_res_annot$annotation.genes5%in%promoter_annot,], 
	DA_res_annot[DA_res_annot$annotation.genesOverlap%in%promoter_annot,] )
DA_res_annot_promoter=DA_res_annot_promoter|>arrange(FDR)


# closest peak in promoter
DA_res_annot.closest.promoter=DA_res_annot_closest[DA_res_annot_closest$annotation.closest%in%promoter_annot,]
DA_res_annot.closest.promoter=DA_res_annot.closest.promoter|>arrange(FDR)

save.table(df=DA_res_annot.closest.promoter, file=fname_annot_closest_prom,dir=sub_resdir)


## ---- peaks_highest_count

n_top_peaks=500

select_exprs=order(rowMeans(logreads.TMM.combat),decreasing=TRUE)[1:n_top_peaks]

DA_res_annot.top500=DA_res_annot_closest[DA_res_annot_closest$PeakID%in%rownames(logreads.TMM.combat[select_exprs,]),]

## save table
fname_string1=paste("ATAC_peaks_annot",DAstring,sep=".")
fname_string2=paste("ATAC_peaks_heatmap",DAstring,sep=".")

fname_subset="top500peaks"

subset_resdir=file.path(res_subdir,fname_subset)
dir.create(subset_resdir)

fname_top500peaks=paste(fname_string1,fname_subset,file_pref,"tsv",sep=".")
save.table(df=DA_res_annot.top500, file=fname_top500peaks,dir=subset_resdir)

fname_top500peaks_heatmap=paste(fname_string2,fname_subset,file_pref,sep=".")

caption_top500peaks_heatmap="ATAC-seq signal of top 500 peaks, by highest signal."

peaks_for_hm2=DA_res_annot.top500$PeakID
reads.TMM.log.hm2=logreads.TMM.combat[rownames(logreads.TMM.combat)%in%peaks_for_hm2,]


hm2=pheatmap(reads.TMM.log.hm2, 
	cluster_rows=TRUE, 
	scale="row",
	column_title="top 500 peaks by signal", 
	heatmap_legend_param=list(title="Z-score"),
	fontsize_row = 7,
	show_rownames = FALSE,
	color=colour_scale_hm_div)


save.plots(plot_ts=hm2,savedir=subset_resdir,pref=fname_top500peaks_heatmap)

## ---- top100_promoters

atacseq_res=DA_res_annot_closest
colnames(atacseq_res)=gsub(".closest","",colnames(atacseq_res))

DA_res_annot.closest.promoter_top100=atacseq_res|>
			dplyr::select(logFC,logCPM,annotation,geneId,SYMBOL,PeakID)|>
			dplyr::filter(stringr::str_detect(annotation, "Promoter"))|>
			dplyr::arrange(desc(abs(logFC)))|>
			slice_head(n=100)

## save table

fname_string1=paste("ATAC_peaks_annot",DAstring,sep=".")
fname_string2=paste("ATAC_peaks_heatmap",DAstring,sep=".")
fname_subset="top100promoters"

subset_resdir=file.path(res_subdir,fname_subset)
dir.create(subset_resdir)

fname_top100prom=paste(fname_string1,fname_subset,file_pref,"tsv",sep=".")
save.table(df=DA_res_annot.closest.promoter_top100, file=fname_top100prom,dir=subset_resdir)


fname_top100prom_heatmap=paste(fname_string2,fname_subset,file_pref,sep=".")

caption_top100prom_heatmap="ATAC-seq signal of top 100 promoter peaks, by absolute fold change."

peaks_for_hm=DA_res_annot.closest.promoter_top100$PeakID
reads.TMM.log.hm=logreads.TMM.combat[rownames(logreads.TMM.combat)%in%peaks_for_hm,]


hm1.1=pheatmap(reads.TMM.log.hm, 
	cluster_rows=TRUE,  
	scale="row",
	column_title="top100promoters", 
	heatmap_legend_param=list(title="Z score"),
	fontsize_row = 7,
	show_rownames = FALSE,
	color=colour_scale_hm_div )


# rotate the cols dendrogram
hc=hclust(dist(t(reads.TMM.log.hm)))
col_dendro=as.dendrogram(hc)
col_dendro_rotated=dendextend::rotate(col_dendro, c("pp_R2","pp_R1","pp_R3","R2","R1","R3"))

reads_scaled=t(scale(t(reads.TMM.log.hm)))

hm1=ComplexHeatmap::Heatmap(reads_scaled, 
	cluster_rows=TRUE,  
	cluster_columns=col_dendro_rotated,
	column_title="top100promoters", 
	heatmap_legend_param=list(title="Z score"),
	show_row_names = FALSE,
	col=colour_scale_hm_div)

save.plots(plot_ts=hm1,savedir=subset_resdir,pref=fname_top100prom_heatmap)



## ---- bindetect_volcano

drer_hsa_ortholog.pth=file.path(annotdir,"Drer_hsa.orthology-alliance.tsv")

hsa_drer_ortho=read.delim(drer_hsa_ortholog.pth, sep="\t", header=FALSE, quote = "")
hsa_drer_ortho=hsa_drer_ortho[-1,]
hsa_drer_ortho=hsa_drer_ortho[c(1,2,5,6)]
colnames(hsa_drer_ortho)=c("HGNC-id","name","zfin-id","Symbol")

bindetect_dir=file.path(atac_datadir,"BINDdetect")
bindetect_res=file.path(bindetect_dir,"bindetect_results.txt")
lec_bec_tffprnt=read.delim(bindetect_res, sep="\t", header=TRUE, quote = "")
colnames(lec_bec_tffprnt)=gsub("_merged_replicates_footprints","", colnames(lec_bec_tffprnt))
lec_bec_tffprnt$name=toupper(lec_bec_tffprnt$name)

lec_bec_tffprnt_annot=left_join(lec_bec_tffprnt,hsa_drer_ortho, by="name", relationship = "many-to-many")



df_volcano=lec_bec_tffprnt|>
	dplyr::select(output_prefix,name,motif_id,LEC_BEC_change,LEC_BEC_pvalue,LEC_BEC_highlighted)|>
	dplyr::mutate(minuslog10pval=-log10(LEC_BEC_pvalue))|>
	dplyr::mutate(plot_colour=case_when(LEC_BEC_highlighted=="True" & LEC_BEC_change>0 ~ "up in LEC",
							LEC_BEC_highlighted=="True" & LEC_BEC_change<0 ~ "up in BEC",
							TRUE ~ "na"))

mycolour2.1 = c('na'="gray80", 'up in BEC' = BEC_dark,'up in LEC'= LEC_dark)

volcano_diff=ggplot(df_volcano ,aes(x=LEC_BEC_change,y=minuslog10pval, color=plot_colour)) +
	geom_point(alpha=0.8) +
	scale_color_manual(values=mycolour2.1, name = "Footprint score FC") +
  theme_bw() + theme(legend.position = "right")+ theme(aspect.ratio = 1)+
  ylab("-log10(pval)") + xlab("Footprint score fold change LEC vs. BEC") 


sub_resdir=file.path(res_subdir,"TF_footprinting_summary")
dir.create(sub_resdir)

volcano_bindif_pref="TF_footprinting_LECvsBEC_volcano"

save.plots(plot_ts=volcano_diff,savedir=sub_resdir,pref=volcano_bindif_pref)


## ---- mafba_targets


#              output_prefix name           motif_id
# 810           Mafb_MA0117.3 MAFB           MA0117.3
# 880           MAFB_MA0117.2 MAFB           MA0117.2
# 881           MAFB_MA0117.1 MAFB           MA0117.1
# 882 MAFB_MAFB.H12CORE.1.P.B MAFB MAFB.H12CORE.1.P.B
# 883         MAFB_YGVTGASTCA MAFB         YGVTGASTCA
# 884         MAFB_NGMTGASTCA MAFB         NGMTGASTCA


prom_annot=c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)")


# get bound sites


get_bound=function(TF_name=TF_name,BINDetect_dir=BINDetect_dir,atacseq_res=atacseq_res){
	
	bound_targets=list()

	TF_beds=file.path(BINDetect_dir,TF_name,"beds")

	smpl=c("LEC","BEC")

	for (smpl.i in smpl){

		bed_bound=file.path(TF_beds,paste0(TF_name,"_",smpl.i,"_merged_replicates_footprints_bound.bed"))
		bed_unbound=file.path(TF_beds,paste0(TF_name,"_",smpl.i,"_merged_replicates_footprints_unbound.bed"))

		bound=read.delim(bed_bound,sep="\t", header=FALSE)
		colnames(bound)=c("TFBS_chr.fp","TFBS_start.fp","TFBS_end.fp","TFBS_motif.fp","TFBS_score.fp","strand.fp","chr_peak.fp","start_peak.fp","end_peak.fp","PeakID","bed1","bed2","ensgeneid.fp","symbol.fp","footprint_score.fp")

		bound_annot_smpl.i=left_join(bound,atacseq_res, by="PeakID")

		bound_targets[[smpl.i]]=bound_annot_smpl.i
	}

	return(bound_targets)
}




get_TF_venn=function(targets=targets,f_pref=f_pref,outdir=outdir){

  targets_lec=targets[["LEC"]]
  targets_lec.prom=targets_lec|>dplyr::filter(annotation%in%prom_annot)
  targets_lec.ls=unique(targets_lec.prom$SYMBOL[!is.na(targets_lec.prom$SYMBOL)])

  targets_bec=targets[["BEC"]]
  targets_bec.prom=targets_bec|>dplyr::filter(annotation%in%prom_annot)
  targets_bec.ls=unique(targets_bec.prom$SYMBOL[!is.na(targets_bec.prom$SYMBOL)])

  targets_lec_bec=intersect(targets_lec.ls,targets_bec.ls)

  n_overlap=length(targets_lec_bec)
  n_lec=length(setdiff(targets_lec.ls,targets_lec_bec))
  n_bec=length(setdiff(targets_bec.ls,targets_lec_bec))

  targets_lec.df=targets_lec.prom|>
  	dplyr::filter(SYMBOL%in%targets_lec.ls)|>
  	dplyr::select(PeakID,Chr,Start,End,logFC,FDR,TFBS_start.fp,TFBS_end.fp,TFBS_motif.fp,annotation,geneId,SYMBOL)

  targets_bec.df=targets_bec.prom|>
  	dplyr::filter(SYMBOL%in%targets_bec.ls)|>
  	dplyr::select(PeakID,Chr,Start,End,logFC,FDR,TFBS_start.fp,TFBS_end.fp,TFBS_motif.fp,annotation,geneId,SYMBOL)

  targets_lec_bec.df.l=targets_lec.prom|>
  	dplyr::filter(SYMBOL%in%targets_lec_bec)|>
  	dplyr::select(PeakID,Chr,Start,End,logFC,FDR,TFBS_start.fp,TFBS_end.fp,TFBS_motif.fp,annotation,geneId,SYMBOL)|>
  	dplyr::mutate(target="LEC")

  targets_lec_bec.df.b=targets_bec.prom|>
  	dplyr::filter(SYMBOL%in%targets_lec_bec)|>
  	dplyr::select(PeakID,Chr,Start,End,logFC,FDR,TFBS_start.fp,TFBS_end.fp,TFBS_motif.fp,annotation,geneId,SYMBOL)|>
  	dplyr::mutate(target="BEC")

  targets_lec_bec.df=rbind(targets_lec_bec.df.l,targets_lec_bec.df.b)

  f_pref.l=paste(f_pref,"LEC","tsv",sep=".")
  save.table(df=targets_lec.df,file=f_pref.l,dir=outdir)

  f_pref.b=paste(f_pref,"BEC","tsv",sep=".")
  save.table(df=targets_bec.df,file=f_pref.b,dir=outdir)

  f_pref.lb=paste(f_pref,"LEC_BEC","tsv",sep=".")
  save.table(df=targets_lec_bec.df,file=f_pref.lb,dir=outdir)

  set.seed(314)
  list_intersections = euler(c("LEC" = n_lec, "BEC" = n_bec,
                                 "LEC&BEC" = n_overlap), shape = "ellipse")

  return(list_intersections)

}


sub_resdir=file.path(res_subdir,"TF_footprinting_Mafb")
dir.create(sub_resdir)


# endothelial motif

TF_name="MAFB_NGMTGASTCA"

sub_resdir.TF=file.path(sub_resdir,TF_name)
dir.create(sub_resdir.TF)
table_fname_pref=paste("Mafb_targets",TF_name,sep=".")

targets_NGM=get_bound(TF_name=TF_name,BINDetect_dir=bindetect_dir,atacseq_res=atacseq_res)
venns_NGM=get_TF_venn(targets=targets_NGM,f_pref=table_fname_pref,outdir=sub_resdir.TF)
venn_mafba_endothelial=plot(venns_NGM,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))


venn_mafba_endo.fname=paste("Mafba_TFfootprnt_Venn","promoter",TF_name,file_pref, sep=".")
save.plots(plot_ts=venn_mafba_endothelial,savedir=sub_resdir.TF,pref=venn_mafba_endo.fname)



# combo motifs Jaspar + Hocomoco

TF_name="Mafb_MA0117.3"
targets_ma0117.3=get_bound(TF_name=TF_name,BINDetect_dir=bindetect_dir,atacseq_res=atacseq_res)

TF_name="MAFB_MA0117.2"
targets_ma0117.2=get_bound(TF_name=TF_name,BINDetect_dir=bindetect_dir,atacseq_res=atacseq_res)

TF_name="MAFB_MA0117.1"
targets_ma0117.1=get_bound(TF_name=TF_name,BINDetect_dir=bindetect_dir,atacseq_res=atacseq_res)

TF_name="MAFB_MAFB.H12CORE.1.P.B"
targets_H12CORE=get_bound(TF_name=TF_name,BINDetect_dir=bindetect_dir,atacseq_res=atacseq_res)

mafb_db.lec=rbind(targets_ma0117.3$LEC, targets_ma0117.2$LEC, targets_ma0117.1$LEC, targets_H12CORE$LEC)
mafb_db.bec=rbind(targets_ma0117.3$BEC, targets_ma0117.2$BEC, targets_ma0117.1$BEC, targets_H12CORE$BEC)

targets_mafb_db=list()
targets_mafb_db[["LEC"]]=mafb_db.lec
targets_mafb_db[["BEC"]]=mafb_db.bec


TF_name="MAFB_Jaspar_Hocomoco"

sub_resdir.TF=file.path(sub_resdir,TF_name)
dir.create(sub_resdir.TF)
table_fname_pref=paste("Mafb_targets",TF_name,sep=".")


venns_combo=get_TF_venn(targets=targets_mafb_db,f_pref=table_fname_pref,outdir=sub_resdir.TF)
venn_mafba_combo=plot(venns_combo,  quantities = TRUE, fill_alpha = 0.5, 
  fills = c(LEC_light, BEC_light, "#80808060"),
  edges=c(LEC_dark,BEC_dark,"#808080"))


venn_mafba_combo.fname=paste("Mafba_TFfootprnt_Venn","promoter",TF_name,file_pref, sep=".")
save.plots(plot_ts=venn_mafba_combo,savedir=sub_resdir.TF,pref=venn_mafba_combo.fname)



## ---- foo



## ---- foo
