#!/usr/bin/env Rscript

#rackham.uppmax.uu.se

#interactive -A uppmax2025-2-45 -t 2:00:00 -p core -n 10

#singularity exec /proj/snic2019-30-13/nobackup/nbis5931/tools/TADCompare/agatasm-tadcompare_r443.img bash

#Rscript TADCompare_consensus.R

## OBS!
## OnTAD cat bed files had to be manually edited
# /proj/snic2019-30-13/nobackup/nbis5931/analysis2024/TAD/OnTAD/OnTAD_24i2025_proc/LEC/Lsize_7_penalty_0.075_maxsize_300/LEC_OnTAD_24i2025.Lsize_7_penalty_0.075_maxsize_300.bed
# to rm these config lines, one per each chr
#track name="OnTAD 16" description="OnTAD 16" visibility=2 itemRgb="On"
#track name="OnTAD 24" description="OnTAD 24" visibility=2 itemRgb="On"



library(dplyr)
library(ggplot2)
library(HiCcompare)
library(TADCompare)

options(scipen=999)


wrkdir=""
setwd(wrkdir)

datadir=""

sp_BEC_101 = readRDS(file.path(datadir,"sparse_101_10kb_BEC.rds"))
sp_BEC_103 = readRDS(file.path(datadir,"sparse_103_10kb_BEC.rds"))
sp_LEC_102 = readRDS(file.path(datadir,"sparse_102_10kb_LEC.rds"))
sp_LEC_104 = readRDS(file.path(datadir,"sparse_104_10kb_LEC.rds"))



#for file names
today=""

########################
########################
#### read in TADs

read_ontad<-function(path=path){
	tads_df=read.csv(path, sep="\t", header=FALSE)


	tads_df=tads_df %>%
    dplyr::mutate(TAD_level = case_when(V9 == '56,108,176' ~ '1',
                           V9 == '127,201,127' ~ '2',
                           V9 == '190,174,212' ~ '3',
                           V9 == '253,192,134' ~ '4',
                           V9 == '255,0,0' ~ '5',                           
                           TRUE ~ '0')) %>%
    dplyr::filter(TAD_level!='0')


	colnames(tads_df)[1:9]=c("chr","start","end","tad_bed_4","tad_bed_5","tad_bed_6","tad_bed_7","tad_bed_8","tad_bed_9")

	tads_df$chr=as.character(tads_df$chr)

	return(tads_df)
}


### megamaps TADs

ontad_resdir=""
settings_pref="Lsize_7_penalty_0.075_maxsize_300"

taddir_lec=file.path(ontad_resdir,"LEC",settings_pref)
taddir_bec=file.path(ontad_resdir,"BEC",settings_pref)

lec_tad=file.path(taddir_lec,paste0("LEC","_OnTAD_24i2025.",settings_pref,".sorted.ed.bed") )
bec_tad=file.path(taddir_bec,paste0("BEC","_OnTAD_24i2025.",settings_pref,".sorted.ed.bed") )

TADs_lec=read_ontad(path=lec_tad)
TADs_bec=read_ontad(path=bec_tad)

TADs_lec.4=TADs_lec[,c(1:4)]
TADs_bec.4=TADs_bec[,c(1:4)]

#################################
#################################
#### TADcompare


get_diffTADs<-function(mat1=mat1,mat2=mat2,tad_list=tad_list, z_cutoff=z_cutoff){

	res_list=list()

	cmp_TAD_res_allchr=list()
	cmp_TAD_res_allchr_BS=list()

	for (i in names(mat1) ) {

		TADs.1.i=tad_list[[1]]|>dplyr::filter(chr==i)
		TADs.2.i=tad_list[[2]]|>dplyr::filter(chr==i)
		TAD_lst.i=list(TADs.1.i,TADs.2.i)

		TAD_cmp=TADCompare(mat1[[i]], mat2[[i]], resolution = 10000, pre_tads = TAD_lst.i, window_size=mat_bins, z_thresh=z_cutoff)

		cmp_TAD_res=TAD_cmp$TAD_Frame
		cmp_TAD_res=cmp_TAD_res|>
			dplyr::mutate(chr=i, .after=Boundary)

		cmp_TAD_res_bs=TAD_cmp$Boundary_Scores
		cmp_TAD_res_bs=cmp_TAD_res_bs|>
			dplyr::mutate(chr=i, .after=Boundary)

		cmp_TAD_res_allchr[[i]]=cmp_TAD_res
		cmp_TAD_res_allchr_BS[[i]]=cmp_TAD_res_bs

	}

	cmp_TAD_res_allchr=bind_rows(cmp_TAD_res_allchr, .id = "df_label")
	cmp_TAD_res_allchr=cmp_TAD_res_allchr[,!names(cmp_TAD_res_allchr)%in%c("df_label")]

	cmp_TAD_res_allchr_BS=bind_rows(cmp_TAD_res_allchr_BS, .id = "df_label")
	cmp_TAD_res_allchr_BS=cmp_TAD_res_allchr_BS[,!names(cmp_TAD_res_allchr_BS)%in%c("df_label")]


	res_list[["diff"]]=cmp_TAD_res_allchr
	res_list[["all_bins"]]=cmp_TAD_res_allchr_BS
	
	return(res_list)

}


######################
#### megamap TADs, sample matrices

mat_bins=50
z_cutoff=2.5

#### LEC

smpl_pref="LEC_mega_TADs.102_vs_104"

TAD_lst= list(TADs_lec.4, TADs_lec.4)

tadcompare_res=get_diffTADs(mat1=sp_LEC_102,mat2=sp_LEC_104,tad_list=TAD_lst,z_cutoff=z_cutoff)

saveRDS(tadcompare_res[[1]], file = paste0("TADcompare.ontad.diff_zCO2_5.",smpl_pref,".",settings_pref,".windows_",mat_bins,".",today,".rds"))
saveRDS(tadcompare_res[[2]], file = paste0("TADcompare.ontad.BoundaryScore_allpos_zCO2_5.",smpl_pref,".",settings_pref,".windows_",mat_bins,".",today,".rds"))



#### BEC

smpl_pref="BEC_mega_TADs.101_vs_103"

TAD_lst= list(TADs_bec.4, TADs_bec.4)

tadcompare_res=get_diffTADs(mat1=sp_BEC_101,mat2=sp_BEC_103,tad_list=TAD_lst)

saveRDS(tadcompare_res[[1]], file = paste0("TADcompare.ontad.diff_zCO2_5.",smpl_pref,".",settings_pref,".windows_",mat_bins,".",today,".rds"))
saveRDS(tadcompare_res[[2]], file = paste0("TADcompare.ontad.BoundaryScore_allpos_zCO2_5.",smpl_pref,".",settings_pref,".windows_",mat_bins,".",today,".rds"))


