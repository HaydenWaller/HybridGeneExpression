BiocManager::install("pcaExplorer", dependencies = TRUE, force = TRUE)

library(pcaExplorer)
library(tidyverse)

countmatrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/merged.isoform.counts.matrix_all")
coldata<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/all_samples_pcaexplorer.txt")

pcaExplorer(countmatrix = countmatrix, coldata = coldata)
