argv=commandArgs(T)

library(dplyr)
library(tidyverse)
library(readr)


#cell.lines = read.delim('/hdd/yanhuan/micro/software/SAHMI/cell.lines.txt', header = T,sep=' ') %>% tibble() 
#cell.lines=cell.lines %>% 
#    filter (reads>0 & uniq >0)
#qtile = 0.90
#q_df = cell.lines %>%
#group_by(name,rank)%>%
#summarize(CLrpmm = 10^quantile(log10(rpmm),qtile, na.rm = T),
#.groups = 'keep')


#q_df=read.delim('/hdd/yanhuan/micro/software/SAHMI/cell_line_quantile.txt', header = T,sep='\t') %>% tibble()
#q_df=q_df[,c('taxid','rank','threshold_0.90_exclude0')]
#colnames(q_df)=c('taxid','rank','CLrpmm')  

q_df=read.delim('/hdd/yanhuan/micro/software/SAHMI/CCLE_quantile_mhc_rrna.txt', header = T,sep='\t') 
#q_df %>% subset( taxid == 10407)
#t(q_df %>% subset( taxid %in% c(210,10407)) )
q_df=q_df[,c('taxid','rank','threshold_0.9_exclude0')]
colnames(q_df)=c('taxid','rank','CLrpmm') 
q_df$CLrpmm=ifelse(q_df$taxid == 10407,10.507194,q_df$CLrpmm) # hbv 保留为 threshold_0.95_add0
#q_df$taxid=as.character(q_df$taxid)
#rownames(q_df)=q_df$taxid

kr=read.csv(argv[2],header=T,sep='\t') 
kr=kr %>% subset(rank %in% c('G','S'))
#kr$taxid=as.character(kr$taxid)

sprintf("kraken rpmm: %d",nrow(kr))
bh=read.csv(argv[1],header=T,sep='\t') 
bh$name<- NULL
#bh$taxid=as.character(bh$taxid)

kr_sub <- kr %>%
  select(taxid, name, rpmm) %>%
  distinct(taxid, .keep_all = TRUE) 
bh_sub <- bh %>%
  distinct(taxid, .keep_all = TRUE) 
q_df_sub <- q_df %>%
  distinct(taxid, .keep_all = TRUE) 

merge_bn_qdf_kr  <- full_join(kr_sub, bh_sub, by = "taxid") %>%
  left_join(q_df_sub, by = "taxid")
merge_bn_qdf_kr$CLrpmm[is.na(merge_bn_qdf_kr$CLrpmm)]= 0
hits_df <- merge_bn_qdf_kr

sample_level=read.csv(argv[6],header=T,sep='\t')
#hits_df=hits_df[ which(hits_df$taxid %in% sample_level$taxid), ]
s1=left_join(hits_df,sample_level, by = c('taxid' ))
s1$name.x<-NULL
s1$rank.x<-NULL
s1$name=s1$name.y
s1$rank=s1$rank.y
s1$name.y<-NULL
s1$rank.y<-NULL


kr_report=read.csv(argv[7],header=F,sep='\t')
kr_1=unique(data.frame(kr_report$V6,kr_report$V7,kr_report$V8))
colnames(kr_1)=c('rank','taxid','name')
s1$rank=ifelse(s1$taxid %in% kr_1$taxid, kr_1[match(s1$taxid,kr_1$taxid),]$rank,NA)
s1$name=ifelse(s1$taxid %in% kr_1$taxid, trimws(kr_1[match(s1$taxid,kr_1$taxid),]$name),NA)
s1=s1 %>% subset(taxid != 9606 & taxid !=9605)
s1$mreads=ifelse(s1$taxid %in% kr$taxid, kr[match(s1$taxid,kr$taxid),]$reads,NA)

sprintf("cell line quantile test(exclude zero counts,0.90): %d",nrow(hits_df))
write.table(s1,argv[4],row.names=F,col.names=T,quote=F,sep='\t')
write.table(s1$taxid,argv[5],row.names=F,col.names=F,quote=F,sep='\t')
