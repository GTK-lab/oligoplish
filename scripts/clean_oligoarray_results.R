## This generates the oligopaints bed files without the added primers

library(here)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
project <- "controls"


projectdir <- project

target_ranges <- rtracklayer::import(sprintf("%s/targets.bed",projectdir))
oa_ranges <- rtracklayer::import(sprintf("%s/oa_targets.bed",projectdir))
target_ranges
oa_ranges
oligos <- read_tsv(sprintf('%s/oligos.txt',projectdir),
                   col_types = cols(),
                   col_names <- c("name","pos","len",
                             "deltaG","deltaH","deltaS","Tm","targets","seq")) %>%
    filter(!str_detect(targets,";"))

oligo_table <- inner_join(as_data_frame(as.data.frame(oa_ranges)),oligos,by="name") %>%
    mutate(gene=str_extract(name,"^[A-Za-z0-9]+")) %>%
    group_by(gene) %>%
    mutate(name=sprintf("%s_o%02d",gene,row_number()),
           start=start+pos-1,
           end=start+len)  %>% ungroup()

rtracklayer::export(oligo_table,sprintf("%s/oligos.bed",projectdir))

dna <- Biostrings::DNAStringSet(oligo_table$seq)
names(dna) <- oligo_table$name
writeXStringSet(dna,sprintf("%s/oligopaints.fa",projectdir))


for (g in unique(oligo_table$gene)) {
    print(g)
    filter(oligo_table,gene==g) %>%
        dplyr::select(seqnames,start,end,name=seq,score=Tm) %>%
        rtracklayer::export(sprintf("%s/%s_oligopaints.bed",projectdir,g))
    }
