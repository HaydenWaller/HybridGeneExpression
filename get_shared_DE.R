library(tidyverse)

#get shared on QTL
lklp<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/lk_vs_lp.annotated.results.recip_uniprot.DEonly.txt")
fastslow<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/fast_vs_slow.annotated.results.recip_uniprot.DEonly.txt")
QTL<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/all_peaks-CIs_Chapter1.txt")
GMAP<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/laupala_GMAP_filtered_withOrthos")

shared<-merge(lklp, fastslow, by="Orthogroup")
shared_mapped<-merge(shared, GMAP, by="Orthogroup")
shared_QTL<-merge(shared_mapped, QTL, by="scaffold")

write_tsv(shared_QTL, file="C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/shared_on_QTL.recip")

#get lklp unique on QTL

lklp_mapped<-merge(lklp, GMAP, by ="Orthogroup")
lklp_QTL<-merge(lklp_mapped, QTL, by = "scaffold")

write_tsv(lklp_QTL, file="C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/lk_vs_lp.DE_on_QTL.recips")

#get fastslow unique on QTL
fastslow_mapped<-merge(fastslow, GMAP, by ="Orthogroup")
fastslow_QTL<-merge(fastslow_mapped, QTL, by = "scaffold")

write_tsv(fastslow_QTL, file="C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/fast_vs_slow.DE_on_QTL.recips")
