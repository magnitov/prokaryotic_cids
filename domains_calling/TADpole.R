rm(list = ls(all.names = TRUE))
library(TADpole)

our.data <- read.table(file = "in_folder/filename_rep2.tab", row.names = NULL, 
                       sep = "\t", quote = "", stringsAsFactors = FALSE, header = F) 

bed <- read.table(file = "in_folder/filename.bed", row.names = NULL, 
                  sep = "\t", quote = "", stringsAsFactors = FALSE, header = F)
resolt = 15000

tadpole.CID <- function(map, bed, resolution){
  kr <- c(tail(bed, n=1))
  last.CID <- as.numeric(kr$V3)
  write.table(map, paste0("./tmp"), quote=F, row.names=F, col.names=F,sep = "\t")
  mat <- "./tmp"
  tadpole <- TADpole(mat, chr = "chr1", start = 0, end = last.CID, resol = resolution)
  i=0
  for (cluster in tadpole$clusters){
    i= i + 1
    tadpole.results <- data.frame(Reduce(cbind, cluster))
    tadpole.results$init <- (tadpole.results$init-1)*resolution
    tadpole.results$V2 <- tadpole.results$V2*resolution
    tadpole.results <- tadpole.results[which((tadpole.results$V2 - tadpole.results$init)>=4*resolution),-4]
    tadpole.results$init <- format(tadpole.results$init, scientific = F)
    tadpole.results$V2 <- format (tadpole.results$V2, scientific = F)
    new <- c('chr1')
    tadpole.results <- cbind(new, tadpole.results)
    write.table(tadpole.results, paste0("output_folder/filename_rep2_level_", i), quote=F, 
                row.names=F, col.names=F, sep = "\t")}
}

tadpole.CID(our.data, bed, resolt)
