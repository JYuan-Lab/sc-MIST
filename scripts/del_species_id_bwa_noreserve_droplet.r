# 筛选出通过gnues bwa但species 没有通过的reads对应的spcies
# 例如 species 29382没有通过bwa，但其genus 1279通过了，但结果中不应出现29382这类结果

library(tidyr)
library(dplyr)
library(readxl)
library(parallel)
library(data.table)
library(optparse)

option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--fa2"), action="store", help = "path to fastq 2 file"),
  make_option(c("--work_path"), action="store", help = "output directory path"),
  make_option(c("--out_file"), action="store", help = "output directory path"),
  make_option(c("--condaenv"), action="store", help = "conda env path --- env/bin"),
  make_option(c("--reference_fa"), action="store", help = "path to each reference fa"),
  make_option(c("--read_linea"), action="store", help = "path to kraken output"),
  make_option(c("--bwa_filter"), action="store", help = "path to filter bwa result .sh")
    
)
opt = parse_args(OptionParser(option_list = option_list))

Sys.setenv(PATH = paste(opt$condaenv, Sys.getenv("PATH"), sep=":"))
#bwa_fa=opt$reference_fa
sample=opt$sample
#fq_path=opt$fq2
#pwdpath=opt$work_path
#outdir=opt$out_file


linea=read.csv(opt$read_linea,header=T,sep='\t')



#micro_out=fread(opt$micro_out,header=F,sep='\t')
#micro_out$taxid=unlist( lapply(micro_out$V3, function(x) { return( unlist( strsplit( x, "\\(taxid " ) )[ 2 ] ); } #) )
#micro_out$taxid=gsub(')','',micro_out$taxid)
#micro_out = micro_out %>% subset( taxid %in% as.character( linea$V7 ))


each_taxid_bwa<- function( ta){
    #ta='457429'
    #bwa_ta_path=paste0(opt$reference_fa,'/',ta,'.fa')
    if (ta %in% linea$gtaxid) {
        linea_1=linea %>% subset(  gtaxid  == ta) } else {
        linea_1=linea %>% subset(  staxid  == ta)
    }
    path_ta_srrid=paste0(opt$work_path,'/',opt$out_file,'/',sample,'_',ta,'srr_id')
    write.table(linea_1$readid,path_ta_srrid,col.names = F,row.names = F,quote=F)

    

    sub_fa_2= paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '_reads.sub_2.fa')
    str=paste0( 'seqkit grep -j 10 -f ',path_ta_srrid, '  ', opt$fa2 ,'>',sub_fa_2)
    system(str)
    sub_fa_rename_2= paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '_reads.sub_rename_2.fa')
    str=paste0( 'seqkit replace -p " .*" -w 0 ', sub_fa_2 ,'> ', sub_fa_rename_2 ) 
    system(str)
    # fa2fastq
    fa2fq_2= paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '_reads.sub_rename_2.fq')
    str=paste0( 'seqtk seq -F "I"  ', sub_fa_rename_2 ,'> ',fa2fq_2)
    system(str)
    
    # bwa
    bwa_sam=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '.sam')
    str1=paste0('/hdd/yanhuan/envs/tum/bin/bwa mem -t 4 ',paste0(opt$reference_fa,'/',ta,'.fa'),' ',' ',fa2fq_2,'> ', bwa_sam)
    system(str1)
#
    bwa_bam=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '.bam')
    str3=paste0('/hdd/yanhuan/envs/tum/bin/samtools collate ',bwa_sam, '  -o ', bwa_bam )
    system(str3)

    # get reads id 
    filter_id<- paste0(opt$work_path,'/',opt$out_file, '/',sample, '_',ta,'_map_reads_id')
    str2=paste0('bash ', opt$bwa_filter,' ',bwa_bam, '> ', filter_id)
    system(str2)

    if (file.info(filter_id)$size == 0) {
        del_id <- ta  # 空 data frame
    } else {
        del_id <- 'no'
    }
    system(paste0('rm ',sub_fa_2)) 
    system(paste0('rm ',path_ta_srrid))
    system(paste0('rm ', bwa_sam))
    system(paste0('rm ', sub_fa_rename_2))
    system(paste0('rm ', fa2fq_2))
    system(paste0('rm ', filter_id))
    system(paste0('rm ', bwa_bam))

  #  del_reads=setdiff(linea_1$readid,reads_id)
    return(del_id) }

species_id=do.call( rbind, mclapply( unique( c(linea$staxid[!is.na(linea$staxid)]) ), function(x) {each_taxid_bwa(x) },mc.cores=10))
print(species_id)
#species_reads_id=do.call( rbind, lapply( c(1485), function(x) {each_taxid_bwa(x) }))

path3=paste0( opt$work_path,'/',opt$out_file,'/',sample,'_bwa_black_species_id')
print(path3)
write.table(data.frame(taxid=species_id) %>% subset( taxid != 'no' ),path3,row.names = F,col.names = F,sep='\t',quote=F)