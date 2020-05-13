rm(list = ls(all.names = TRUE))
library(chromoR)

chromo.map <- read.table(file = "in_folder/filename_rep2_raw.tab", row.names = NULL, 
                         sep = "\t", quote = "", stringsAsFactors = FALSE, header = F) 

chromo.bed <- read.table(file = "in_folder/filename.bed", row.names = NULL, 
                         sep = "\t", quote = "", stringsAsFactors = FALSE, header = F)
resolt = 10000

CIM.CID <- function (map, bed, resolution){
  chr.map <-  Reduce(rbind, map)
  bed$V3 <- bed$V3-1
  bed[1,2] <- 1
  indices = which(bed$V1 == "chr1")
  i1 = indices[1]
  iN = indices[length(indices)]
  p = rowSums(chr.map[i1:iN, i1:iN]) - diag(chr.map[i1:iN, i1:iN])
  res = segmentCIM(p)
  vec.cid <- c(0,res$changePoints)
  size.f <- length(vec.cid)
  final.format <- data.frame(first = 'chr1', start = numeric(size.f-1), end = numeric(size.f-1))
  final.format$start <- vec.cid[1:size.f-1]*resolution
  final.format$end <- vec.cid[2:size.f]*resolution
  final.format <- final.format[which((final.format$end - final.format$start)>=4*resolution),-4]
  write.table(final.format, paste0("out_folder/filename_rep2"), quote=F, 
              row.names=F, col.names=F, sep = "\t")
}

CIM.CID(chromo.map, chromo.bed, resolt)



