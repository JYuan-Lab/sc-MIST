library(tidyr)
library(dplyr)
library(readxl)
library(parallel)
library(data.table)
library(optparse)

option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--fq1"), action="store", help = "path to fastq 1 file"),
  make_option(c("--fq2"), action="store", help = "path to fastq 2 file"),
  make_option(c("--work_path"), action="store", help = "output directory path"),
  make_option(c("--out_file"), action="store", help = "output directory path"),
  make_option(c("--condaenv"), action="store", help = "conda env path --- env/bin"),
  make_option(c("--reference_fa"), action="store", help = "path to HLA reference fa"),
  make_option(c("--sfq1"), action="store", help = "path to species fastq 1 file"),
  make_option(c("--sfq2"), action="store", help = "path to speceis fastq 2 file")
)
opt = parse_args(OptionParser(option_list = option_list))
Sys.setenv(PATH = paste(opt$condaenv, Sys.getenv("PATH"), sep=":"))

bwa_sam=paste0(opt$work_path,'/',opt$out_file, '/',opt$sample,'_mhc.sam')
str1=paste0('/hdd/yanhuan/envs/tum/bin/bwa mem -t 4 ',opt$reference_fa,' ',opt$fq2,'> ', bwa_sam)
system(str1)
bwa_bam=paste0(opt$work_path,'/',opt$out_file, '/',opt$sample,'_mhc.bam')
str3=paste0('/hdd/yanhuan/envs/tum/bin/samtools collate ',bwa_sam, '  -o ', bwa_bam )
system(str3)

bwa_unmappe_bam=paste0(opt$work_path,'/',opt$out_file, '/',opt$sample,'_unmapped_mhc_read_id')
str4=paste0('/hdd/yanhuan/envs/tum/bin/samtools view -f 4 ',bwa_bam,' | cut -f1 | sort -u > ',bwa_unmappe_bam )
system(str4)

get_fq_1=paste0( 'seqkit grep -j 10 -f ',bwa_unmappe_bam,' ', opt$fq1,' > ', opt$sfq1)
system(get_fq_1)

get_fq_2=paste0( 'seqkit grep -j 10 -f ',bwa_unmappe_bam,' ', opt$fq2,' > ', opt$sfq2)
system(get_fq_2)