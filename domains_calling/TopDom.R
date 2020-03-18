library(TopDom)
library(data.table)

rm(list = ls(all.names = TRUE))

folder <- c("analysis/VERONIKA/COVERAGE")

param <- '5M'
name.list <- c(paste0("Ecoli_10kb_", param, "_rep1"),
               paste0("Ecoli_10kb_", param, "_rep2"))
resolution <- 10000
windows.size <- 2:20

lapply(name.list, function(rep) {
  tsv <- read.table(file.path(folder, paste0(rep, ".tab")), 
                    row.names = NULL, sep = "\t", quote = "", stringsAsFactors = FALSE, header = F)
  bed <- read.table(file.path(folder, "ec1_10kb_windows.bed"),
                    row.names = NULL, sep = "\t", quot = "", stringsAsFactors = FALSE, header = F)
  sum <- cbind(bed, tsv)
  sum <- format(sum, scientific = FALSE)
  path <- file.path(folder, "TopDom", paste0(rep, "_names"))
  write.table(sum, path,
              quote=F, row.names=F, col.names=F, sep="\t")
  lapply(windows.size, function(i) {
    topdom.res <- TopDom(data = path, window.size = i)
    topdom.res.bed <- topdom.res$bed
    topdom.res.bed <- topdom.res.bed[which(topdom.res.bed$name == "domain" & 
                   (topdom.res.bed$chromEnd - topdom.res.bed$chromStart) >= 4*resolution), -4]
    write.table(topdom.res.bed, file.path(folder, "TopDom", 
                                          paste0("Ecoli_10kb_", param, "_rep2_ws", i, ".bed")),
                quote=F, row.names=F, col.names=F,sep = "\t")
  })
})






