library(tidyr)
library(dplyr)
library(readxl)
library(parallel)
library(data.table)
library(optparse)

option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--fa1"), action="store", help = "path to fastq 1 file"),
  make_option(c("--fa2"), action="store", help = "path to fastq 2 file"),
  make_option(c("--work_path"), action="store", help = "output directory path"),
  make_option(c("--out_file"), action="store", help = "output directory path"),
  make_option(c("--condaenv"), action="store", help = "conda env path --- env/bin"),
  make_option(c("--reference_fa"), action="store", help = "path to each reference fa"),
  make_option(c("--read_linea"), action="store", help = "path to kraken output"),
  #make_option(c("--micro_out"), action="store", help = "path to step3 microbiome out "),
  make_option(c("--sfq1"), action="store", help = "path to species fastq 1 file"),
  make_option(c("--sfq2"), action="store", help = "path to speceis fastq 2 file"),
  make_option(c("--bwa_filter"), action="store", help = "path to speceis fastq 2 file"),
  make_option(c("--bc"), action="store", default=16, help = "nucleutide length of cell barcodes"),
  make_option(c("--umi"), action="store", default=10, help = "nucleutide length of umis")
    
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

    # step2 get species fasta
    sub_fa= paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '_reads.sub.fa')
    str=paste0( 'seqkit grep -j 10 -f ',path_ta_srrid, '  ', opt$fa2 ,'>',sub_fa)
    system(str)
    sub_fa_rename= paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '_reads.sub_rename.fa')
    str=paste0( 'seqkit replace -p " .*" -w 0 ', sub_fa ,'> ', sub_fa_rename ) 
    system(str)
    # fa2fastq
    fa2fq= paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '_reads.sub_rename.fq')
    str=paste0( 'seqtk seq -F "I"  ', sub_fa_rename ,'> ',fa2fq)
    system(str)
    # bwa
    bwa_sam=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_',ta, '.sam')
    str1=paste0('/hdd/yanhuan/envs/tum/bin/bwa mem -t 4 ',paste0(opt$reference_fa,'/',ta,'.fa'),' ',fa2fq,'> ', bwa_sam)
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
        reads_id <- data.frame()  # 空 data frame
    } else {
        reads_id <- read.csv(filter_id, header=FALSE, sep="\t")
    }
    system(paste0('rm ',sub_fa))
    system(paste0('rm ',path_ta_srrid))
#    system(paste0('rm ', bwa_sam))
    system(paste0('rm ', sub_fa_rename))
    system(paste0('rm ', fa2fq))
    system(paste0('rm ', filter_id))
    system(paste0('rm ', bwa_bam))

#    del_reads=setdiff(linea_1$readid,reads_id)
    
    
    return(reads_id) }
species_reads_id=do.call( rbind, mclapply( unique(c(linea$gtaxid[!is.na(linea$gtaxid)],linea$staxid[!is.na(linea$staxid)])), function(x) {each_taxid_bwa(x) },mc.cores=10))

#species_reads_id=do.call( rbind, lapply( c(1485), function(x) {each_taxid_bwa(x) }))

path3=paste0( opt$work_path,'/',opt$out_file,'/',sample,'_bwa_reserve_srrid')
print(path3)
write.table(species_reads_id,path3,row.names = F,col.names = F,sep='\t',quote=F)


# get all reserve species reads
f1=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_reserve_1.fa')
f2=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_reserve_2.fa')
str_fa1=paste0('seqkit grep -j 10 -f ',path3,' ', opt$fa1, '> ', f1)
str_fa2=paste0('seqkit grep -j 10 -f ',path3,' ', opt$fa2, '> ', f2)
system(str_fa1)
system(str_fa2)

f1_re=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_reserve_rename_1.fa')
f2_re=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_reserve_rename_2.fa')
str_fa1=paste0( 'seqkit replace -p " .*" -w 0 ', f1 ,'> ', f1_re )
str_fa2=paste0( 'seqkit replace -p " .*" -w 0 ', f2 ,'> ', f2_re )
system(str_fa1)
system(str_fa2)

#bcumi=opt$bc+opt$umi
#print(bcumi)
#f1_re_1=paste0(opt$work_path,'/',opt$out_file, '/',sample,'_reserve_rename_1_filter.fa')
#str_fa1_ren=paste0('seqkit subseq -r 1:',bcumi, ' ', f1_re, '> ',f1_re_1)
#str_fa1_ren
#system(str_fa1_ren)

# fa2fq
str_fq1=paste0( 'seqtk seq -F "I"  ', f1_re ,'> ',opt$sfq1)
str_fq2=paste0( 'seqtk seq -F "I"  ', f2_re ,'> ',opt$sfq2)
system(str_fq1)
system(str_fq2)

system(paste0('rm ',f1))
system(paste0('rm ',f2))
system(paste0('rm ',f1_re))
#system(paste0('rm ',f1_re_1))
system(paste0('rm ',f2_re))

print('finish')
