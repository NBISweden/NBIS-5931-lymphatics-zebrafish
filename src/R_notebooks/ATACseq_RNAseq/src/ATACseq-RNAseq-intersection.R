
## ---- prep-env

library(renv)

##renv::restore()
renv::activate()


library(sva)
library(edgeR)

library(tidyverse)
library(tibble)
library(dplyr)
library(DT)
library(htmltools)
library(magrittr)

library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggupset)
library(ggimage)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(eulerr)

library(knitr)
library(kableExtra)

##renv::snapshot()


tab_n_sv=0
tab_n=0 # table labels
fig_n=0
fig_n_cap=0

# for file names, consistent for all notebooks
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

mycolour2 = c('na'="gray80", 'up in BEC' = "#762A83",'up in LEC'= "#109630")


## ---- cutoffs

rnaseq_log2FC_CO=0.3
rnaseq_pct_expr_CO=0.3
rnaseq_padj_CO=0.1

atac_lFC_CO=0.5
#atac_lFC_CO=0
atac_FDR_CO=0.05



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

resdir=file.path(projrootdir,"results",paste0("ATACseq_RNAseq_",file_pref))
dir.create(resdir, recursive=TRUE)


## ---- annots




## ---- data-rnaseq

rnaseq_datadir=file.path(datadir,"RNAseq","RNA_seq_Marleen_29x2024")

DE_d5=file.path(rnaseq_datadir,"fullmarkers_5dpf.csv")
#DE_d7=file.path(rnaseq_datadir,"fullmarkers_7dpf.csv")

DE_rnaseq_d5=read.csv(DE_d5)
#DE_rnaseq_d7=read.csv(DE_d7)

#pct.1 is fraction in cluster 1 (LEC)

sc_markers_dpf5=file.path(datadir,"RNAseq","TOP50_DE_genes_5dpf_LEC_BEC.txt")
df_markers_dpf5=read.table(sc_markers_dpf5, sep="\t", header=TRUE, blank.lines.skip=TRUE)

genes_lec_dpf5=df_markers_dpf5[df_markers_dpf5$avg_log2FC>0,]
genes_lec_dpf5.ls=genes_lec_dpf5$names

genes_bec_dpf5=df_markers_dpf5[df_markers_dpf5$avg_log2FC<0,]
genes_bec_dpf5.ls=genes_bec_dpf5$names

rnaseq_ave_expr_5=file.path(datadir,"RNAseq","average_expression_5days.csv")
#rnaseq_fc_5=file.path(datadir,"RNAseq","240202_fold_change_5days.csv")




## ---- data-atacseq

atac_datadir=file.path(datadir,"ATACseq")
cnt_table.pth=file.path(atac_datadir,"consensus_merged_macs2broad.tsv")

cnt_table=read.table(cnt_table.pth, sep="\t", header=TRUE, blank.lines.skip=TRUE)
rownames(cnt_table)=cnt_table$Geneid
colnames(cnt_table)=c(colnames(cnt_table)[1:6],gsub(".filtered.sorted.bam","",colnames(cnt_table)[7:12]))
colnames(cnt_table)=c(colnames(cnt_table)[1:6],gsub("P17204_","",colnames(cnt_table)[7:12]))

cnt_table=cnt_table[!cnt_table$Chr=="chrM",]

reads.peak=cnt_table[,c(7:12)]

groups <- factor(c(rep("LEC",3),rep("BEC",3)))
sorting=c(1,2,3,1,2,3)
sorting=factor(sorting)
design <- model.matrix(~sorting+groups)
rownames(design)=colnames(reads.peak)

reads.dge = DGEList(counts=reads.peak, group=groups)
reads.dge = calcNormFactors(reads.dge)
keep = filterByExpr(reads.dge)
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



atac_resdir=file.path(projrootdir,"results",paste0("ATACseq_",file_pref),
        "consensus_intrareplicate_merged","DA_all_peaks")
atac_res_DA_pth=file.path(atac_resdir,paste("ATAC_peaks_annot_closest.ensembl.DA_groups",file_pref,"tsv",sep="."))

atacseq_res=read.delim(atac_res_DA_pth, sep="\t", header=TRUE, quote = "")

#atacseq_res$SYMBOL.closest=dplyr::na_if(atacseq_res$SYMBOL.closest, "<NA>")

atacseq_res.lec=atacseq_res|>
        dplyr::filter(FDR<atac_FDR_CO & abs(logFC)>=atac_lFC_CO)|>
        dplyr::filter(logFC>0)

atacseq_res.bec=atacseq_res|>
        dplyr::filter(FDR<atac_FDR_CO & abs(logFC)>=atac_lFC_CO)|>
        dplyr::filter(logFC<0)



## ---- atacseq-volcano

df_volcano=atacseq_res|>
        dplyr::select(logFC,logCPM,annotation.closest,geneId.closest,SYMBOL.closest,PeakID,FDR)|>
        dplyr::mutate(minuslog10FDR=-log10(FDR))|>
        dplyr::mutate(LEC_5=case_when(SYMBOL.closest%in%genes_lec_dpf5.ls ~ "up in LEC",
                TRUE ~ "na"))|>
        dplyr::mutate(BEC_5=case_when(SYMBOL.closest%in%genes_bec_dpf5.ls ~ "up in BEC",

                TRUE ~ "na"))
volcano5_lec=ggplot(df_volcano|>arrange(LEC_5) ,aes(x=logFC,y=minuslog10FDR, color=LEC_5)) +
        geom_point(alpha=0.8) +
        scale_color_manual(values=mycolour2, name = "DE genes at 5 dpf") +
        theme_bw() + 
        theme(legend.position = "right")+ 
        theme(aspect.ratio = 1)+
        ylab("-log10(FDR)") + xlab("log2 Fold Change") 

volcano5_bec=ggplot(df_volcano|>arrange(BEC_5) ,aes(x=logFC,y=minuslog10FDR, color=BEC_5)) +
        geom_point(alpha=0.8) +
        scale_color_manual(values=mycolour2, name = "DE genes at 5 dpf") +
        theme_bw() + 
        theme(legend.position = "right")+ 
        theme(aspect.ratio = 1)+
        ylab("-log10(FDR)") + xlab("log2 Fold Change") 


sub_resdir=file.path(resdir,"ATAC_volcano_DEGs")
dir.create(sub_resdir)

fname_volcano_LEC=paste("ATAC_peaks_volcano","LEC_DEGs",file_pref, sep=".")
save.plots(plot_ts=volcano5_lec,savedir=sub_resdir,pref=fname_volcano_LEC)

fname_volcano_BEC=paste("ATAC_peaks_volcano","BEC_DEGs",file_pref, sep=".")
save.plots(plot_ts=volcano5_bec,savedir=sub_resdir,pref=fname_volcano_BEC)


## ---- atacseq-hm-DEGs

caption_dpf5_heatmap="ATAC-seq signal\nLEC and BEC gene promoters\n5 dpf."

genes_dpf5.ls=c(genes_bec_dpf5.ls,genes_lec_dpf5.ls)

prom_peaks_DEG_5dpf=atacseq_res|>
        dplyr::select(logFC,logCPM,annotation.closest,geneId.closest,SYMBOL.closest,PeakID,FDR)|>
        dplyr::filter(SYMBOL.closest%in%genes_dpf5.ls )|>
        dplyr::filter(stringr::str_detect(annotation.closest, "Promoter"))

peaks_for_hm=prom_peaks_DEG_5dpf$PeakID

reads.TMM.log.5dpf=logreads.TMM.combat[rownames(logreads.TMM.combat)%in%peaks_for_hm,]
reads_scaled=t(scale(t(reads.TMM.log.5dpf)))

## annotation with RNA-seq log2FC


# add rnaseq data

rnaseq_expr_5_res=read.delim(rnaseq_ave_expr_5, sep=",", header=TRUE, quote = "", stringsAsFactors=FALSE)
rnaseq_expr_5_res=rnaseq_expr_5_res[,c(1,2,3)]
colnames(rnaseq_expr_5_res)=c("SYMBOL","clust0","clust1")
rnaseq_expr_5_res=as.data.frame(sapply(rnaseq_expr_5_res, function(x) gsub("\"", "", x)))
rnaseq_expr_5_res=tibble(rnaseq_expr_5_res)
rnaseq_expr_5_res$clust0=log2(as.numeric(rnaseq_expr_5_res$clust0)+1)
rnaseq_expr_5_res$clust1=log2(as.numeric(rnaseq_expr_5_res$clust1)+1)

annot_rows=prom_peaks_DEG_5dpf|>
        dplyr::select(PeakID,SYMBOL.closest)|>
        dplyr::rename(SYMBOL=SYMBOL.closest)|>
        dplyr::mutate(tissue=case_when(SYMBOL%in%genes_bec_dpf5.ls ~ "BEC",
                TRUE ~ "LEC"))

annot_rows=left_join(annot_rows,rnaseq_expr_5_res,by="SYMBOL")|>
        dplyr::rename(RNAseq_expr_LEC=clust0,RNAseq_expr_BEC=clust1)

rownames(annot_rows)=annot_rows$PeakID

annot_rows=annot_rows|>
        dplyr::select(tissue,RNAseq_expr_LEC,RNAseq_expr_BEC)|>
        dplyr::arrange(desc(tissue))


# order matrix same as annot
matrix_ord=reads.TMM.log.5dpf[match(rownames(annot_rows),rownames(reads.TMM.log.5dpf)),]

# > range(annot_rows$RNAseq_expr_LEC)
# [1]  3.192984 12.155844
# > range(annot_rows$RNAseq_expr_BEC)
# [1]  0.7916356 12.7691448

col_fun1=colorRamp2(c(0,7,13), c("#440144FF","#2A788EFF","#FDE725FF"))

row_annotations=rowAnnotation(
        Tissue=annot_rows$tissue, 
        RNAseq_expr_LEC=annot_rows$RNAseq_expr_LEC,
        RNAseq_expr_BEC=annot_rows$RNAseq_expr_BEC,
        col=list(Tissue=c(BEC="#762A83",LEC="#109630"), 
                RNAseq_expr_LEC=col_fun1,RNAseq_expr_BEC=col_fun1))

# rotate the dendrogram
hc = hclust(dist(t(reads.TMM.log.5dpf)))
col_dendro=as.dendrogram(hc)
col_dendro_rotated = dendextend::rotate(col_dendro, c("pp_R2","pp_R1","pp_R3","R2","R1","R3"))


hm.da_degs=ComplexHeatmap::Heatmap(reads_scaled,
        cluster_rows=TRUE, show_row_dend = FALSE,
        show_row_names = FALSE,
        heatmap_legend_param=list(title="Z score"),
        cluster_columns=col_dendro_rotated,
        row_split=annot_rows$tissue,     
        col=colour_scale_hm_div,
        column_title=caption_dpf5_heatmap,
        left_annotation=row_annotations
        )

sub_resdir=file.path(resdir,"ATAC_heatmap_DEGs")
dir.create(sub_resdir)

fname_heatmap_DEGs=paste("ATAC_peaks_heatmap","top_DEGs",file_pref, sep=".")
save.plots(plot_ts=hm.da_degs,savedir=sub_resdir,pref=fname_heatmap_DEGs)


## ---- atac-peak-gene-summary

# peaks
peaks_sig=atacseq_res|>
        dplyr::filter(FDR<atac_FDR_CO & abs(logFC)>=atac_lFC_CO)

peaks_sig.prom=peaks_sig|>
        dplyr::filter(stringr::str_detect(annotation.closest, "Promoter"))

peaks_LEC=peaks_sig|>
        dplyr::filter(logFC>0)

peaks_LEC.prom=peaks_sig.prom|>
        dplyr::filter(logFC>0)

peaks_BEC=peaks_sig|>
        dplyr::filter(logFC<0)

peaks_BEC.prom=peaks_sig.prom|>
        dplyr::filter(logFC<0)


#genes
DA_genes.ensembl=unique(peaks_sig$geneId.closest[!is.na(peaks_sig$geneId.closest)])
DA_genes.symbol=unique(peaks_sig$SYMBOL.closest[!is.na(peaks_sig$SYMBOL.closest)])

DA_genes.prom.ensembl=unique(peaks_sig.prom$geneId.closest[!is.na(peaks_sig.prom$geneId.closest)])
DA_genes.prom.symbol=unique(peaks_sig.prom$SYMBOL.closest[!is.na(peaks_sig.prom$SYMBOL.closest)])

DA_genes_LEC.ensembl=unique(peaks_LEC$geneId.closest[!is.na(peaks_LEC$geneId.closest)])
DA_genes_LEC.symbol=unique(peaks_LEC$SYMBOL.closest[!is.na(peaks_LEC$SYMBOL.closest)])

DA_genes_LEC.prom.ensembl=unique(peaks_LEC.prom$geneId.closest[!is.na(peaks_LEC.prom$geneId.closest)])
DA_genes_LEC.prom.symbol=unique(peaks_LEC.prom$SYMBOL.closest[!is.na(peaks_LEC.prom$SYMBOL.closest)])

DA_genes_BEC.ensembl=unique(peaks_BEC$geneId.closest[!is.na(peaks_BEC$geneId.closest)])
DA_genes_BEC.symbol=unique(peaks_BEC$SYMBOL.closest[!is.na(peaks_BEC$SYMBOL.closest)])

DA_genes_BEC.prom.ensembl=unique(peaks_BEC.prom$geneId.closest[!is.na(peaks_BEC.prom$geneId.closest)])
DA_genes_BEC.prom.symbol=unique(peaks_BEC.prom$SYMBOL.closest[!is.na(peaks_BEC.prom$SYMBOL.closest)])


peaks_stats=as.data.frame(rbind(
        c("all peaks",nrow(peaks_sig),length(DA_genes.symbol),length(DA_genes.ensembl)),
        c("promoter peaks",nrow(peaks_sig.prom),length(DA_genes.prom.symbol),length(DA_genes.prom.ensembl)),
        c("LEC DA peaks",nrow(peaks_LEC),length(DA_genes_LEC.symbol),length(DA_genes_LEC.ensembl)),
        c("LEC DA promoter peaks",nrow(peaks_LEC.prom),length(DA_genes_LEC.prom.symbol),length(DA_genes_LEC.prom.ensembl)),
        c("BEC DA peaks",nrow(peaks_BEC),length(DA_genes_BEC.symbol),length(DA_genes_BEC.ensembl)),
        c("BEC DA promoter peaks",nrow(peaks_BEC.prom),length(DA_genes_BEC.prom.symbol),length(DA_genes_BEC.prom.ensembl))
        ))

colnames(peaks_stats)=c("peak subset","DA","genes, SYMBOL","genes, Ensembl")


## ---- rnaseq-gene-summary

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




## ---- atac-rna-intersection


gex_lec.atac_lec=unique(intersect(peaks_LEC$SYMBOL.closest,rnaseq_fc.d5.lec$X))
gex_bec.atac_bec=unique(intersect(peaks_BEC$SYMBOL.closest,rnaseq_fc.d5.bec$X))
gex_lec.atac_bec=unique(intersect(peaks_BEC$SYMBOL.closest,rnaseq_fc.d5.lec$X))
gex_bec.atac_lec=unique(intersect(peaks_LEC$SYMBOL.closest,rnaseq_fc.d5.bec$X))



table_gene_counts=rbind(c("RNA-seq up LEC",length(gex_lec.atac_lec),length(gex_lec.atac_bec),length(rnaseq_fc.d5.lec$X)),
      c("RNA-seq up BEC",length(gex_bec.atac_lec),length(gex_bec.atac_bec),length(rnaseq_fc.d5.bec$X)),
      c("Total genes with DA peaks",length(unique(peaks_LEC$SYMBOL.closest)),length(unique(peaks_BEC$SYMBOL.closest)),"na")
  )

colnames(table_gene_counts)=c("Assay / Tissue","ATAC-seq DA LEC","ATAC-seq DA BEC", "Total DE genes")



table_gene_fractions=as.data.frame(rbind(c("RNA-seq up LEC",length(gex_lec.atac_lec)/length(rnaseq_fc.d5.lec$X),length(gex_lec.atac_bec)/length(rnaseq_fc.d5.lec$X) ),
      c("RNA-seq up BEC",length(gex_bec.atac_lec)/length(rnaseq_fc.d5.bec$X),length(gex_bec.atac_bec)/length(rnaseq_fc.d5.bec$X))
  ))

colnames(table_gene_fractions)=c("Assay / Tissue","ATAC-seq DA LEC","ATAC-seq DA BEC")
table_gene_fractions[,2]=as.numeric(table_gene_fractions[,2])
table_gene_fractions[,3]=as.numeric(table_gene_fractions[,3])

table_gene_fractions=as.data.frame(table_gene_fractions)%>%mutate_if(is.numeric, format, digits=4,nsmall = 0)



## ---- list-DEGs-topDARs

## Figure 3B, genes on the heatmap

peaks_LEC.rnaseq_LEC.top50=peaks_LEC|>
        dplyr::filter(SYMBOL.closest%in%rnaseq_fc.d5.lec$X)|>
        dplyr::filter(logFC>1 & logCPM>2)|>
        group_by(SYMBOL.closest)|>
        filter(logFC == max(logFC))|>
        ungroup()|>
        dplyr::arrange(desc(logFC))|>
        slice_head(n=50)

peaks_LEC.rnaseq_LEC.top50=as.data.frame(peaks_LEC.rnaseq_LEC.top50)

sub_resdir=file.path(resdir,"DEGs.ATAC_DA_list")
dir.create(sub_resdir)

f_pref=paste("DEGs.ATAC_DA_list",file_pref,"tsv",sep=".")
save.table(df=peaks_LEC.rnaseq_LEC.top50,file=f_pref,dir=sub_resdir)

genes.peaks_LEC.rnaseq_LEC.top50=peaks_LEC.rnaseq_LEC.top50$SYMBOL.closest
