rm(list = ls(all.names = TRUE))
library(HiCseg)
our.data <- read.table(file = "/hdd20tb/share/Ecoli_20M_5kb_rep2.tab", row.names = NULL, 
                       sep = "\t", quote = "", stringsAsFactors = FALSE, header = F) 
#!data should be the type of double
n.rows <- nrow(our.data)
our.data <- unlist(our.data)
k.range <- 2:100
resolution <- 5000

lapply(k.range, function(k) {
  our.res <- data.frame(HiCseg_linkC_R(n.rows, k, "G", our.data, "D"))
  n <- nrow(our.res) + 1
  col.elements <- our.res$t_hat*resolution
  bed <- data.frame(chr = 'chr1',
                    from = c(0, col.elements[-length(col.elements)]),
                    to = c(col.elements))
  bed <- bed[which((bed$to-bed$from) >= 4*resolution),] 
  bed <- format(bed, scientific = FALSE)
  write.table(bed,paste0("/hdd20tb/share/Ecoli_20M_5kb_rep2_", k, ".bed"),quote=F, 
              row.names=F, col.names=F,sep = "\t")
})