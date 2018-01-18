## This is a script to create files to use as input for OligoArray2.1
## parameters

## give it a project name for the set
project <- "controls"
## a vector of symbols
syms <- c("ACTB","GAPDH")
## how far upstream of the TSS to target
upstream <- 2000
## how far downstream of the TSS to target
downstream <- 1000
# This file is all hg38 with refseq ids, downloaded from ensembl
hg38_tss_fn <- 'data/hg38-tss.tsv'

## Don't change below this line
library(tidyverse)
library(here)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(GenomicFeatures)

# read in ENSEMBL hg38 transcript tss locations with refseq mrnas
TXDB <- read_tsv(here(hg38_tss_fn),skip=1,
                 col_names=c('gene_id','ts_id','tss','geneid','refseq','symbol','seqname','strand')) %>%
    filter(!is.na(geneid)) %>%
    mutate(strand=fct_recode(factor(strand),"+"="1","-"="-1"),
           start=tss,end=tss)


get_target_ranges <- function(symbols,up=2000,down=1000,txdb=TXDB) {
    tss <- dplyr::filter(txdb,symbol %in% syms) %>%
        as("GRanges") %>% promoters(up,down)
    tssl <- split(tss,tss$symbol)
    tssl <- lapply(tssl,function(x) as(union(x,x),"data.frame"))
    tssl  <- do.call("rbind",tssl) 
    tssl$name <- rownames(tssl)
    tssl
    }

ranges <- get_target_ranges(syms,upstream,downstream) %>%
    rowwise() %>%
    do(tibble(chr=sprintf("chr%s",.$seqnames),
              start=seq(from=.$start,to=.$end,by=950),
              end=.$end,
              name=sprintf(fmt="%s_%d",.$name,seq_along(start)))) %>%
    mutate(end=min(start+999,end)) %>% as("GRanges")

seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,ranges)
names(seqs) <- ranges$name

dir.create(here(sprintf("results/%s",project)))
writeXStringSet(seqs,here(sprintf("results/%s/oligoarray.fa",project)))
write_tsv(as(ranges,"data.frame"),here(sprintf("results/%s/targets.tsv",project)))
write_tsv(as(ranges,"data.frame") %>% dplyr::select(seqnames,start,end,name),
          here(sprintf("results/%s/targets.bed",project)),col_names=FALSE)
