# This has to be run with the "ranges" object still existing from the previous script.

library(here)
library(tidyverse)
oligos <- read_tsv(here(sprintf('results/%s/oligos.txt',project)),
              col_names <- c("name","pos","len",
                             "deltaG","deltaH","deltaS","Tm","targets","seq")) %>%
    filter(!str_detect(targets,";"))

oligo_table <- inner_join(as_data_frame(as.data.frame(ranges)),oligos,by="name")
oligo_bed <- oligo_table %>% transmute(chr=seqnames,start=start+pos-1,end=start+len,name=name)

write_tsv(oligo_bed,here(sprintf("results/%s/oligos.bed",project)),col_names=FALSE)
