rm(list = ls(all.names = TRUE))
library(SpectralTAD)

path.to.input.raw.map <- c('../Data_for_SpectalTAD/COVERAGE/Ecoli_10kb_5M_rep1_raw.tab')
path.to.input.bed <- c('../Data_for_SpectalTAD/COVERAGE/ec1_10kb_windows.bed')
resolution <- 10000
path.to.output <- c('../Data_for_SpectalTAD/COVERAGE/SpectralTAD/Ecoli_10kb_5M_rep1_res')

SpectralTadResults <- function(path.to.input.raw.map, path.to.input.bed, resolution) {
  tsv <- read.table(path.to.input.raw.map, row.names = NULL, 
                  sep = "\t", quote = "", stringsAsFactors = FALSE, header = F)
  bed <- read.table(path.to.input.bed, row.names = NULL, 
                  sep = "\t", quot = "", stringsAsFactors = FALSE, header = F)
  sum <- cbind.data.frame(bed, tsv)
  colnames(sum) = NULL
  results <- SpectralTAD(sum, chr='1', resolution = resolution, min_size = 4)
  results <- as.data.frame(results)
  results <- results[which((results$Level_1.end - results$Level_1.start)>= 4*resolution),-4]
  results$Level_1.chr <- 'chr1'
  results <- format(results, scientific = FALSE)
  write.table(results, path.to.output, quote=F, row.names=F, col.names=F,sep = "\t")
}

SpectralTadResults(path.to.input.raw.map, path.to.input.bed, resolution)

  

