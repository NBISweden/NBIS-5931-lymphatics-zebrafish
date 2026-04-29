#!/usr/bin/env Rscript


#pelle.uppmax.uu.se

#interactive -A uppmax2025-2-469 -t 4:00:00 -n 10

# mount /proj/snic2019-30-13/nobackup/nbis5931

#singularity exec --pwd $(pwd) --bind /proj/snic2019-30-13/nobackup/nbis5931 /proj/snic2019-30-13/nobackup/nbis5931/tools/TADCompare/agatasm-tadcompare_r443.img bash

#Rscript TADCompare_LECvsBEC.R

# compare megamap TADs on tissue matrices to get tissue differential TADs

# OnTAD TAD boundaries scored using several window sizes for contact summarisation
# appropriate for TADs of different length
# prev used TAD boundaries from
# Lsize_7_penalty_0.075_maxsize_300


## OBS!
## OnTAD cat files were manually edited
# /proj/snic2019-30-13/nobackup/nbis5931/analysis2024/TAD/OnTAD/OnTAD_24i2025_proc/LEC/Lsize_7_penalty_0.075_maxsize_300/LEC_OnTAD_24i2025.Lsize_7_penalty_0.075_maxsize_300.bed
# to rm these, one per each chr
#track name="OnTAD 16" description="OnTAD 16" visibility=2 itemRgb="On"
#track name="OnTAD 24" description="OnTAD 24" visibility=2 itemRgb="On"

# compare megamap TADs on tissue matrices to get tissue differential TADs



library(dplyr)
library(ggplot2)
library(HiCcompare)
library(TADCompare)

options(scipen=999)


wrkdir=""
dir.create(wrkdir,recursive=TRUE)

datadir=""

spLEC_mega = readRDS(file.path(datadir,"sparse_mega_10kb_LEC.rds"))
spBEC_mega = readRDS(file.path(datadir,"sparse_mega_10kb_BEC.rds"))

#for file names
today=""



########################
########################
#### read in TADs

read_ontad<-function(path=path){
	tads_df=read.csv(path, sep="\t", header=FALSE)


	tads_df=tads_df |>
    dplyr::mutate(TAD_level = case_when(V9 == '56,108,176' ~ '1',
                           V9 == '127,201,127' ~ '2',
                           V9 == '190,174,212' ~ '3',
                           V9 == '253,192,134' ~ '4',
                           V9 == '255,0,0' ~ '5',                           
                           TRUE ~ '0')) |>
    dplyr::filter(TAD_level!='0')


	colnames(tads_df)[1:9]=c("chr","start","end","tad_bed_4","tad_bed_5","tad_bed_6","tad_bed_7","tad_bed_8","tad_bed_9")

	tads_df$chr=as.character(tads_df$chr)

	return(tads_df)
}


###################
### sample TADs


### megamaps TADs

ontad_resdir=""
settings_pref="Lsize_7_penalty_0.075_maxsize_300"


taddir_lec=file.path(ontad_resdir,"LEC",settings_pref)
taddir_bec=file.path(ontad_resdir,"BEC",settings_pref)

lec_tad=file.path(taddir_lec,paste0("LEC","_OnTAD_24i2025.",settings_pref,".sorted.ed.bed") )
bec_tad=file.path(taddir_bec,paste0("BEC","_OnTAD_24i2025.",settings_pref,".sorted.ed.bed") )

TADs_lec=read_ontad(path=lec_tad)
TADs_bec=read_ontad(path=bec_tad)




#################################
#################################
#### TADcompare


get_diffTADs=function(mat1=mat1,mat2=mat2,tad_list=tad_list,mat_bins=mat_bins, z_co=z_co){

	res_list=list()

	cmp_TAD_res_allchr=list()
	cmp_TAD_res_allchr_BS=list()

	for (i in names(mat1) ) {

		TADs.1.i=tad_list[[1]]%>%dplyr::filter(chr==i)
		TADs.2.i=tad_list[[2]]%>%dplyr::filter(chr==i)
		TAD_lst.i=list(TADs.1.i,TADs.2.i)

		TAD_cmp=TADCompare(mat1[[i]], mat2[[i]], resolution = 10000, pre_tads = TAD_lst.i, window_size=mat_bins, z_thresh=z_co)

		cmp_TAD_res=TAD_cmp$TAD_Frame
		cmp_TAD_res=cmp_TAD_res%>%dplyr::mutate(chr=i, .after=Boundary)

		cmp_TAD_res_bs=TAD_cmp$Boundary_Scores
		cmp_TAD_res_bs=cmp_TAD_res_bs%>%dplyr::mutate(chr=i, .after=Boundary)

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
######################

# window sizes in bins (10kb bin)
window_sizes=c(20,35,50)

# for stat collection (stdout)
stats_header=c("smpl_pref","win_len","diff","nodiff")
write(stats_header,file=file.path(wrkdir,paste0("tadcompare_",today,".stats.txt")),sep="\t",append=TRUE)



######################
#### megamap TADs, megamap matrices

smpl_pref="LEC_BEC_mega_TADs.LEC_vs_BEC"


TADs_lec.4=TADs_lec[,c(1:4)]
TADs_bec.4=TADs_bec[,c(1:4)]

TAD_lst= list(TADs_lec.4, TADs_bec.4)


for (i in window_sizes){

	mat_bins=i

	tadcompare_res=get_diffTADs(mat1=spLEC_mega,mat2=spBEC_mega,tad_list=TAD_lst,mat_bins=mat_bins,z_co=1.8)

	outdir=file.path(wrkdir,mat_bins)
	dir.create(outdir,recursive=TRUE)
	smpl_suf=paste0(smpl_pref,".windows_",mat_bins,".",settings_pref,".",today)
	fname1=paste("TADcompare.ontad.diff",smpl_suf,"rds",sep=".")
	fname2=paste("TADcompare.ontad.BoundaryScore_allpos",smpl_suf,"rds",sep=".")

	saveRDS(tadcompare_res[[1]], file=file.path(outdir,fname1) )
	saveRDS(tadcompare_res[[2]], file=file.path(outdir,fname2) )

	# for stats saved to a file

	#diff
	diftads=tadcompare_res[[1]]|>dplyr::filter(Differential!="Non-Differential")
	nodiftads=tadcompare_res[[1]]|>dplyr::filter(Differential=="Non-Differential")


	stats_line=c(smpl_pref,mat_bins,nrow(diftads),nrow(nodiftads))
	write(stats_line,file=file.path(wrkdir,paste0("tadcompare_",today,".stats.txt")),sep="\t",append=TRUE)
	
}




