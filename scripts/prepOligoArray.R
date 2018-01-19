#!/usr/bin/env Rscript

## This is a script to create files to use as input for OligoArray2.1
## parameters
cat("setting up\n")
library(argparse)
suppressPackageStartupMessages(library(tidyverse,quietly=TRUE,warn.conflicts = FALSE))
suppressPackageStartupMessages(library(here,quietly=TRUE))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38,quietly=TRUE,warn.conflicts = FALSE))
suppressPackageStartupMessages(library(GenomicRanges,quietly=TRUE,warn.conflicts = FALSE))
suppressPackageStartupMessages(library(GenomicFeatures,quietly=TRUE,warn.conflicts = FALSE))
## give it a project name for the set
parser <- ArgumentParser(description="Prepare promoters to run OligoArray2.1")
parser$add_argument("-u","--up",dest="upstream",default=2000,
                    help="length of upstream sequence, in bp (2000)")
parser$add_argument("-d","--down",dest="downstream",default=1000,
                    help="length of downstream sequence, in bp (1000)")
parser$add_argument("--project",dest="project",default=".",
                    help="directory where data will go")
args <- parser$parse_args()
upstream <- args$upstream
downstream <- args$downstream
project <- args$project
cat("ready to go\n")
f <- file("stdin")
syms <- readLines(f)
cat(sprintf("generating sequences for %d gene symbols\n",length(syms)))

                                        # This file is all hg38 with refseq ids, downloaded from ensembl
hg38_tss_fn <- 'data/hg38-tss.tsv'

## Don't change below this line

                                        # read in ENSEMBL hg38 transcript tss locations with refseq mrnas
txdb <- read_tsv(here(hg38_tss_fn),skip=1,col_types=cols(),
                 col_names=c('gene_id','ts_id','tss','geneid','refseq','symbol','seqname','strand')) %>%
    filter(!is.na(geneid),!str_detect(seqname,"[_.]")) %>%
    mutate(strand=fct_recode(factor(strand),"+"="1","-"="-1"),
           start=tss,end=tss)

for (s in syms) {
    if (!(s %in% txdb$symbol)) {
        cat(sprintf("Symbol '%s' not found.. skipping\n",s))
    }
}
outdir <- sprintf("results/%s",project)
cat(sprintf("results stored in %s\n",outdir))
dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

tss <- dplyr::filter(txdb,symbol %in% syms) %>%
    as("GRanges") %>% promoters(upstream,downstream)
tssl <- split(tss,tss$symbol)
tssl <- lapply(tssl,function(x) as(union(x,x),"data.frame"))
target_ranges <- do.call("rbind",tssl) 
target_ranges$name <- rownames(target_ranges)
rtracklayer::export(target_ranges,sprintf("%s/targets.bed",outdir))

oa_ranges <- rowwise(target_ranges) %>%
    do(tibble(chr=sprintf("chr%s",.$seqnames),
              start=seq(from=.$start,to=.$end,by=950),
              end=.$end,
              name=sprintf(fmt="%s_%d",.$name,seq_along(start)))) %>%
    mutate(end=min(start+999,end)) %>% as("GRanges")

seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,oa_ranges)
names(seqs) <- oa_ranges$name

rtracklayer::export(oa_ranges,sprintf("%s/oa_targets.bed",outdir))
writeXStringSet(seqs,sprintf("%s/oa_targets.fa",outdir))
## write_tsv(as(oa_ranges,"data.frame"),sprintf("%s/targets.tsv",outdir))
## write_tsv(as(oa_ranges,"data.frame") %>% dplyr::select(seqnames,start,end,name),
##           sprintf("%s/targets.bed",outdir),col_names=FALSE)
