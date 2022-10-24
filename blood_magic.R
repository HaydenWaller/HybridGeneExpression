library(tidyverse)

#isolate scOrthos
groups <- read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/Orthogroups.tsv")
one_one<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/Orthogroups_SingleCopyOrthologues.txt")

groups_one <- groups %>% filter(Orthogroup %in% one_one$Orthogroup)
write_tsv(groups_one,"C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/Orthogroups_SingleCopyOrthologues.tsv")


# merge matrices
fast.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/fast.isoform.counts.matrix")
slow.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/slow.isoform.counts.matrix")
lk.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/lk.isoform.counts.matrix")
lp.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/lp.isoform.counts.matrix")

lk.lp.matrix<-merge(lk.matrix, lp.matrix, by = "Orthogroup")
lk.lp.fast.matrix<-merge(lk.lp.matrix, fast.matrix, by = "Orthogroup")
lk.lp.fast.slow.matrix<-merge(lk.lp.fast.matrix, slow.matrix, by = "Orthogroup")

write_tsv(lk.lp.fast.slow.matrix, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/merged.isoform.counts.matrix")


#Merge TMM EXPR matrices
fast.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/fast.isoform.TMM.EXPR.matrix")
slow.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/slow.isoform.TMM.EXPR.matrix")
lk.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/lk.isoform.TMM.EXPR.matrix")
lp.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/lp.isoform.TMM.EXPR.matrix")

lk.lp.matrix<-merge(lk.matrix, lp.matrix, by = "Orthogroup")
lk.lp.fast.matrix<-merge(lk.lp.matrix, fast.matrix, by = "Orthogroup")
lk.lp.fast.slow.matrix<-merge(lk.lp.fast.matrix, slow.matrix, by = "Orthogroup")

write_tsv(lk.lp.fast.slow.matrix, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/matrices/merged.isoform.TMM.EXPR.matrix_all")
