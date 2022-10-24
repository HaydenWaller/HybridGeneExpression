library(tidyverse)
library(dplyr)

laupala_gmap<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/laupala_GMAP", sep="\t", skip=2, header=FALSE)
laupala_gmap_filtered<-laupala_gmap[laupala_gmap$V3 == "gene",]
write_tsv(laupala_gmap_filtered, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/laupala_GMAP_filtered")

#Step 2 - Clean and prepare filtered GMAP file along with all of the DE results files in Excel

#Step 3 - Merge cleaned GMAP file with Orthogroups file

orthos<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lk_transcripts_and_orthos.txt")
gmap_filtered<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/laupala_GMAP_filtered_cleaned.txt")
gmap_filtered_withOrthos<-merge(orthos, gmap_filtered, by="Transcript")


write_tsv(gmap_filtered_withOrthos, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/laupala_GMAP_filtered_withOrthos")


#Step 4 - Merge Step 3 output with full DE results

lk_lp<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lk_vs_lp.annotated.results")
lk_lp_GMAP<-merge(lk_lp, gmap_filtered_withOrthos, by = "Orthogroup")

lk_slow<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lk_vs_slow.annotated.results")
lk_slow_GMAP<-merge(lk_slow, gmap_filtered_withOrthos, by = "Orthogroup")

fast_slow<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fast_vs_slow.annotated.results")
fast_slow_GMAP<-merge(fast_slow, gmap_filtered_withOrthos, by = "Orthogroup")

fast_lk<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fast_vs_lk.annotated.results")
fast_lk_GMAP<-merge(fast_lk, gmap_filtered_withOrthos, by = "Orthogroup")

fast_lp<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fast_vs_lp.annotated.results")
fast_lp_GMAP<-merge(fast_lp, gmap_filtered_withOrthos, by = "Orthogroup")

lp_slow<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lp_vs_slow.annotated.results")
lp_slow_GMAP<-merge(lp_slow, gmap_filtered_withOrthos, by = "Orthogroup")


write_tsv(lk_lp_GMAP, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lk_vs_lp.annotated.GMAP.results")
write_tsv(lk_slow_GMAP, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lk_vs_slow.annotated.GMAP.results")
write_tsv(fast_slow_GMAP, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fast_vs_slow.annotated.GMAP.results")
write_tsv(fast_lk_GMAP, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fast_vs_lk.annotated.GMAP.results")
write_tsv(fast_lp_GMAP, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fast_vs_lp.annotated.GMAP.results")
write_tsv(lp_slow_GMAP, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lp_vs_slow.annotated.GMAP.results")

#Step 5 - Merge Step 4 output with peak/CIs

all_peaks<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/RNA-seq/reanalysis_project/all_peak-CIs_Chapter1.tsv")
fastslow_GMAP<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fast_vs_slow.annotated.GMAP.results")
lklp_GMAP<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lk_vs_lp.annotated.GMAP.results")

lklp_candidates<-merge(all_peaks, lklp_GMAP, by = "Scaffold")
fastslow_candidates<-merge(all_peaks, fastslow_GMAP, by = "Scaffold")

write_tsv(lklp_candidates, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/lklp_candidates")
write_tsv(fastslow_candidates, file = "C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/hybridz/annotations/fastslow_candidates")



#Step 6 - Filter Step 5 outputs to include only annotations with reciptrocal blast hits

#make recips file
doop<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/5_species_RNAseq/no_tantalis/annotations/lk_scOrthoProts.pep_blastp.outfmt6")
pood<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/5_species_RNAseq/no_tantalis/annotations/uniprot_blastp.outfmt6")

poopdood<-merge(doop, pood, by = "transcript", all = TRUE)


write_tsv(poopdood, file="C:/Users/hayde/Documents/GRAD_SKEWL/5_species_RNAseq/no_tantalis/annotations/recips_uniprot")

##clean file above, save as "annots_with_reciprocal_hits.uniprot.txt"

cleaned_recips_file<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/annots_with_reciprocal_hits.uniprot.txt")
lk_OGs<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/lk_OGs.txt")
recips_with_OGs<-merge(cleaned_recips_file, lk_OGs, by = "transcript", all = TRUE)

write_tsv(recips_with_OGs, file="C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/recips_with_OGs.uniprot")


#Merge after cleaning file above again (remove all columns except OG and annot)

fast_x_slow<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/DE/merged.isoform.counts.matrix_all.FAST_vs_SLOW.DESeq2.DE_results")
lk_x_lp<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/DE/merged.isoform.counts.matrix_all.KOHALENSIS_vs_PARANIGRA.DESeq2.DE_results")

fast_vs_slow.annotated.results<-merge(fast_x_slow, recips_with_OGs, by="Orthogroup", all=TRUE)
lk_vs_lp.annotated.results<-merge(lk_x_lp, recips_with_OGs, by="Orthogroup", all=TRUE)


write_tsv(fast_vs_slow.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/fast_vs_slow.annotated.results.recip_uniprot")
write_tsv(lk_vs_lp.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/recips/lk_vs_lp.annotated.results.recip_uniprot")
