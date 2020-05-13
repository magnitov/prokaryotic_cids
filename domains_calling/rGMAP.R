rm(list = ls(all.names = TRUE))
library(rGMAP)
set.seed(1)

our.data <- read.table(file = "in_folder/filename_rep2_raw.tab", row.names = NULL, 
                        sep = "\t", quote = "", stringsAsFactors = FALSE, header = F) 
resolt = 3000

RGMAP.CID <- function (map.cid, resolution){
  pre.results = rGMAP(map.cid, resl = resolution, dom_order = 1, min_d = 10, 
                      min_dp = 4, Max_dp =100,  fcthr = 0.75)
  first.column <- c('chr1')
  final.results <- pre.results$tads
  final.results <- cbind(first.column, final.results)
  final.results$start <- final.results$start - resolution/2
  final.results$end <- final.results$end - resolution/2
  write.table(final.results, "out_folder/filename_rep2", quote=F, 
              row.names=F, col.names=F, sep = "\t")
}

RGMAP.CID(our.data, resolt)
