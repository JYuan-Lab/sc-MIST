library(dplyr)
library(tidyverse)
argv=args=commandArgs(T)

script=argv[1]
kraken_report=argv[2]
sckmer=argv[3]
barcode_hit=argv[4]
barcode_tax=argv[5]
kraken_report_all=argv[6]
source(script)

# report = read.delim('output/kraken2/P1/SCAF2963_3_Live.kraken.report.txt', header = F)
report = read.delim(argv[2], header = F)
report$V8 = trimws(report$V8)
report[report$V8 %in% c('Homo sapiens', 'Bacteria', 'Fungi', 'Viruses'), ]

# sckmer data
kmer_data = read.table(argv[3], header = T,sep=' ')
head(kmer_data)

#c = kmer_data %>% 
#  subset(kmer > 1) %>%
#  group_by(taxid) %>%
#  mutate(nn = n()) %>%
#  subset(nn > 3) %>% 
#  group_by(taxid) %>%
#  summarize(r =  cor.test(kmer, uniq, method = 'spearman')$estimate,
#            p =   cor.test(kmer, uniq, method = 'spearman')$p.value,
#            .groups='keep') %>%
#  mutate(p = p.adjust(p))

c= kmer_data %>%
  filter(kmer > 1) %>%
  group_by(taxid) %>%
  summarise(
    nn = n(),
    r = if (nn > 3 && sum(is.finite(kmer) & is.finite(uniq)) >= 2)
      suppressWarnings(as.numeric(cor.test(kmer, uniq, method = "spearman")$estimate))
    else NA_real_,
    p = if (nn > 3 && sum(is.finite(kmer) & is.finite(uniq)) >= 2)
      suppressWarnings(cor.test(kmer, uniq, method = "spearman")$p.value)
    else NA_real_,
    .groups = "keep"
  ) %>%
  mutate(p = p.adjust(p) )
c$nn<- NULL

print(c)
c$name = report$V8[match(c$taxid, report$V7)] # add taxa names 


#c=c[which(c$p < 0.05),]


x = read.delim(kraken_report_all, header = F)
x$V8 = trimws(x$V8)
n_microbiome_reads = sum(x$V2[x$V8 %in% c('Bacteria', 'Fungi', 'Viruses')])


k1 <- read_kraken_reports_v2(argv[2],n_microbiome_reads)
write.table(c, file=argv[4], sep="\t")
write.table(k1, file=argv[5], sep="\t")
