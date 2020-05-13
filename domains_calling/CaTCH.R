rm(list = ls(all.names = TRUE))
library(CaTCH)

#3-columns matrix
rep.file <- read.table(file = "in_folder/filename_rep2", row.names = NULL, 
                        sep = "\t", quote = "", stringsAsFactors = FALSE, header = F)
resolt = 3000

CATCH.CIDS <- function (map, resolution){
  map$V4 <- as.numeric(format(map$V4, scientific = F))
  map$V4 <- map$V4*10000
  write.table(map,paste0("./tmp"),quote=F, row.names=F, col.names=F,sep = " ")
  fileinput<-"./tmp"
  resul.map <- domain.call(fileinput)
  resul.map <- resul.map$clusters
  size.f <- length(which(resul.map$RI < 0.0011))
  final.format <- data.frame(first = 'chr1',start = numeric(size.f), end = numeric(size.f))
  final.format$start <- (resul.map$start[which(resul.map$RI < 0.0011)]-1)*resolution
  final.format$end <- resul.map$end[which(resul.map$RI < 0.0011)]*resolution
  final.format <- final.format[which((final.format$end - final.format$start)>=4*resolution),-4]
  write.table(final.format, paste0("out_folder/filename_rep2"), quote=F, 
              row.names=F, col.names=F, sep = "\t")
  
}

CATCH.CIDS(rep.file, resolt)

