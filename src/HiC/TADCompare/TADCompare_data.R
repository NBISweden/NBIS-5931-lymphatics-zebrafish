# format and save data for TADcompare

#interactive -A uppmax2025-2-234 -t 2:00:00 -p core -n 10
#module load bioinfo-tools R_packages/4.3.1


library(dplyr)
#library(ggplot2)
library(HiCcompare)


wrkdir=""
setwd(wrkdir)

datadir=""

mat101=file.path(datadir,"101.chr/101.chr.BEC_P23106_2.10k_10000.cool")
mat102=file.path(datadir,"102.chr/102.chr.LEC_P23106_2.10k_10000.cool")
mat103=file.path(datadir,"103.chr/103.chr.BEC_P23106_2.10k_10000.cool")
mat104=file.path(datadir,"104.chr/104.chr.LEC_P23106_2.10k_10000.cool")



read_in_cool<-function(cool_pth=cool_pth){
	bedpe_lst=cooler2bedpe(path=cool_pth)
	bedpe_lst_cis=bedpe_lst$cis

	bedpe_df=bind_rows(bedpe_lst_cis, .id="df_label")
	bedpe2sparse=bedpe_df[,!names(bedpe_df)%in%c("df_label")]

	sparse_mat=cooler2sparse(bedpe2sparse)

	return(sparse_mat)
}


sp101=read_in_cool(cool_pth=mat101)
saveRDS(sp101, file="sparse_101_10kb_BEC.rds")

sp102=read_in_cool(cool_pth=mat102)
saveRDS(sp102, file="sparse_102_10kb_LEC.rds")

sp103=read_in_cool(cool_pth=mat103)
saveRDS(sp103, file="sparse_103_10kb_BEC.rds")

sp104=read_in_cool(cool_pth=mat104)
saveRDS(sp104, file="sparse_104_10kb_LEC.rds")

