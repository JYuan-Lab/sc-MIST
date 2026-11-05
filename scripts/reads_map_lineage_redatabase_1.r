library(tidyr)
library(dplyr)
library(data.table)
library(readxl)

argv=commandArgs(T)
path1=argv[1]
micro_out=fread(path1,header=F,sep='\t')
micro_out$taxid=unlist( lapply(micro_out$V3, function(x) { return( unlist( strsplit( x, "\\(taxid " ) )[ 2 ] ); } ) )
micro_out$taxid=gsub(')','',micro_out$taxid)


linea=read.csv(argv[2],header=T,sep='\t',check.names = F)
linea=linea %>% subset( taxid %in% micro_out$taxid & (genus != "" | species != ""))
linea=linea[,c('taxid','kindom','g_taxid','s_taxid','species')]
linea$taxid=as.character(linea$taxid)

database_score=read_xlsx(argv[3],sheet = 1)
database_score=database_score %>% subset( database_score$...17 >= argv[5] )

micro_out_1=micro_out %>% subset( taxid %in% linea$taxid)
micro_out_linea=merge(micro_out_1[,c('V2','taxid')],linea,by='taxid')

colnames(micro_out_linea)=c('readtaxid','readid','kindom',	'gtaxid',	'staxid','species')

#micro_out_linea = micro_out_linea %>% subset( is.na(staxid) | species %in% database_score$species)
micro_out_linea_viru <- micro_out_linea %>% subset(  !(kindom == "k__Viruses" & is.na(gtaxid) & is.na(staxid))    )
micro_out_linea_bac = micro_out_linea %>% subset( kindom == "k__Bacteria"  ) %>%
        subset( is.na(staxid) | species %in% database_score$species)

micro_read_gs=rbind(micro_out_linea_viru,micro_out_linea_bac)


write.table(micro_read_gs,argv[4],row.names=F,col.names=T,sep='\t',quote=F)
