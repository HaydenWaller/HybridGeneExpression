##Inheritence Plot

inheritance_laupala<- read.table("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/inheritance_table.txt", header = TRUE, stringsAsFactors = FALSE)

##slow hybrids
total_genes <- nrow(inheritance_laupala)
con <- inheritance_laupala[which(inheritance_laupala$padj_kvs > 0.001 & inheritance_laupala$padj_pvs > 0.001 & inheritance_laupala$padj_kvp > 0.001),]

add <- inheritance_laupala[which(inheritance_laupala$padj_kvs < 0.001 & inheritance_laupala$padj_pvs < 0.001 & inheritance_laupala$padj_kvp < 0.001),]
add1 <- add[which(add$log2FoldChange_kvs < -0.2 & add$log2FoldChange_pvs > 0.2),]
add2 <- add[which(add$log2FoldChange_kvs > 0.2 & add$log2FoldChange_pvs < -0.2),]
add <- rbind(add1, add2)

p1_dom <- inheritance_laupala[which(inheritance_laupala$padj_kvs > 0.001 & inheritance_laupala$padj_pvs < 0.001),]
p2_dom <- inheritance_laupala[which(inheritance_laupala$padj_kvs < 0.001 & inheritance_laupala$padj_pvs > 0.001),]
under_dom <- inheritance_laupala[which(inheritance_laupala$padj_kvp < 0.001 & inheritance_laupala$log2FoldChange_kvs > 0.2 & inheritance_laupala$log2FoldChange_pvs > 0.2),]
over_dom <- inheritance_laupala[which(inheritance_laupala$padj_kvp < 0.001 & inheritance_laupala$log2FoldChange_kvs < -0.2 & inheritance_laupala$log2FoldChange_pvs < -0.2),]
transgressive <- rbind(over_dom, under_dom)
dom <- rbind(p1_dom, p2_dom)

prop_con      <-paste(round((nrow(con)/total_genes * 100), digits = 2),"%", sep = "") 
prop_add      <-paste(round((nrow(add)/total_genes * 100), digits = 2),"%", sep = "")
prop_p1_dom   <-paste(round((nrow(p1_dom)/total_genes * 100), digits = 2),"%", sep = "")      
prop_p2_dom   <-paste(round((nrow(p2_dom)/total_genes * 100), digits = 2),"%", sep = "")      
prop_overdom  <-paste(round((nrow(over_dom)/total_genes * 100), digits = 2),"%", sep = "")         
prop_underdom <-paste(round((nrow(under_dom)/total_genes * 100), digits = 2),"%", sep = "")
prop_dom <- paste(round((nrow(dom)/total_genes * 100), digits = 2),"%", sep = "") 
prop_transgressive <- paste(round((nrow(transgressive)/total_genes * 100), digits = 2),"%", sep = "") 

inheritance <- rbind(con, add, dom, transgressive)
head(inheritance)


my_cols<-c("orange red","black","royal blue","blue violet","gray")
cols_in <- my_cols[inheritance$inheritance]

plot(inheritance$log2FoldChange_kvs, inheritance$log2FoldChange_pvs, pch =16, col = "white", 
      ylab = "log2 fold change mh",
      xlab = "log2 fold change ph",cex.axis=1.3, ylim = c(-10,10), xlim = c(-8,8))
  
  points(con$log2FoldChange_kvs,con$log2FoldChange_pvs, col = my_cols[2],pch = 1, cex = 0.5)
  points(transgressive$log2FoldChange_kvs,transgressive$log2FoldChange_pvs, col = my_cols[3],pch = 16, cex = 0.5)
  points(add$log2FoldChange_kvs,add$log2FoldChange_pvs, col = my_cols[4],pch = 16, cex = 0.5)
  points(dom$log2FoldChange_kvs,dom$log2FoldChange_pvs, col = my_cols[1],pch = 16, cex = 0.5)
  
  abline(v=0, col = "black", lty = 3, lwd = 1.8)
  abline(h=0, col = "black", lty = 3, lwd = 1.8)


write.table(add, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/slow_additive_genes")
write.table(transgressive, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/slow_transgressive_genes")


##fast hybrids
total_genes <- nrow(inheritance_laupala)
con <- inheritance_laupala[which(inheritance_laupala$padj_fvk > 0.001 & inheritance_laupala$padj_fvs > 0.001 & inheritance_laupala$padj_kvp > 0.001),]
  
add <- inheritance_laupala[which(inheritance_laupala$padj_fvk < 0.001 & inheritance_laupala$padj_fvp < 0.001 & inheritance_laupala$padj_kvp < 0.001),]
add1 <- add[which(add$log2FoldChange_fvk > -0.2 & add$log2FoldChange_fvp < 0.2),]
add2 <- add[which(add$log2FoldChange_fvk < 0.2 & add$log2FoldChange_fvp > -0.2),]
add <- rbind(add1, add2)
  
p1_dom <- inheritance_laupala[which(inheritance_laupala$padj_fvk < 0.001 & inheritance_laupala$padj_fvp > 0.001),]
p2_dom <- inheritance_laupala[which(inheritance_laupala$padj_fvk > 0.001 & inheritance_laupala$padj_fvp < 0.001),]
under_dom <- inheritance_laupala[which(inheritance_laupala$padj_kvp < 0.001 & inheritance_laupala$log2FoldChange_fvk < -0.2 & inheritance_laupala$log2FoldChange_fvp < -0.2),]
over_dom <- inheritance_laupala[which(inheritance_laupala$padj_kvp < 0.001 & inheritance_laupala$log2FoldChange_fvk > 0.2 & inheritance_laupala$log2FoldChange_fvp > 0.2),]
transgressive <- rbind(over_dom, under_dom)
dom <- rbind(p1_dom, p2_dom)
  
prop_con      <-paste(round((nrow(con)/total_genes * 100), digits = 2),"%", sep = "") 
prop_add      <-paste(round((nrow(add)/total_genes * 100), digits = 2),"%", sep = "")
prop_p1_dom   <-paste(round((nrow(p1_dom)/total_genes * 100), digits = 2),"%", sep = "")      
prop_p2_dom   <-paste(round((nrow(p2_dom)/total_genes * 100), digits = 2),"%", sep = "")      
prop_overdom  <-paste(round((nrow(over_dom)/total_genes * 100), digits = 2),"%", sep = "")         
prop_underdom <-paste(round((nrow(under_dom)/total_genes * 100), digits = 2),"%", sep = "")
prop_dom <- paste(round((nrow(dom)/total_genes * 100), digits = 2),"%", sep = "") 
prop_transgressive <- paste(round((nrow(transgressive)/total_genes * 100), digits = 2),"%", sep = "") 
  
inheritance <- rbind(con, add, dom, transgressive)
head(inheritance)
  
  
my_cols<-c("orange red","black","royal blue","blue violet","gray")
cols_in <- my_cols[inheritance$inheritance]
  
plot(inheritance$log2FoldChange_fvk, inheritance$log2FoldChange_fvp, pch =16, col = "white", 
    ylab = "log2 fold change mh",
    xlab = "log2 fold change ph",cex.axis=1.3, ylim = c(-10,10), xlim = c(-8,8))
  
points(con$log2FoldChange_fvk,con$log2FoldChange_fvp, col = my_cols[2],pch = 1, cex = 0.5)
points(transgressive$log2FoldChange_fvk,transgressive$log2FoldChange_fvp, col = my_cols[3],pch = 16, cex = 0.5)
points(add$log2FoldChange_fvk,add$log2FoldChange_fvp, col = my_cols[4],pch = 16, cex = 0.5)
points(dom$log2FoldChange_fvk,dom$log2FoldChange_fvp, col = my_cols[1],pch = 16, cex = 0.5)
  
abline(v=0, col = "black", lty = 3, lwd = 1.8)
abline(h=0, col = "black", lty = 3, lwd = 1.8)

write.table(add, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/fast_additive_genes")
write.table(transgressive, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/fast_transgressive_genes")

##Get additive Genes

slow_add_OGs<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/slow_additive_OGs.txt", header = TRUE, sep = "\t")
fast_add_OGs<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/fast_additive_OGs.txt", header = TRUE, sep = "\t")
lklp<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/lk_vs_lp.annotated.GMAP.results", header = TRUE, sep = "\t")

slow_add_genes<-merge(slow_add_OGs, lklp, by = "Orthogroup")
fast_add_genes<-merge(fast_add_OGs, lklp, by = "Orthogroup")

write.table(slow_add_genes, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/slow_additive_annots")
write.table(fast_add_genes, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/fast_additive_annots")

QTL_scaffs<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/all_peaks-CIs_Chapter1.txt", header = TRUE, sep = "\t")

slow_add_on_QTL<-merge(slow_add_genes, QTL_scaffs, by = "Scaffold")
fast_add_on_QTL<-merge(fast_add_genes, QTL_scaffs, by = "Scaffold")

##Get transgressive Genes

slow_trans_OGs<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/slow_transgressive_OGs.txt", header = TRUE, sep = "\t")
fast_trans_OGs<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/fast_transgressive_OGs.txt", header = TRUE, sep = "\t")
lklp<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/annotations/lk_vs_lp.annotated.GMAP.results", header = TRUE, sep = "\t")

slow_trans_genes<-merge(slow_trans_OGs, lklp, by = "Orthogroup")
fast_trans_genes<-merge(fast_trans_OGs, lklp, by = "Orthogroup")

write.table(slow_trans_genes, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/slow_transgressive_annots")
write.table(fast_trans_genes, file = "C:/Users/hayde/Documents/GRAD_SKEWL/hybridz/inheritence/fast_transgressive_annots")

QTL_scaffs<-read.delim("C:/Users/hayde/Documents/GRAD_SKEWL/reanalysis_project/all_peaks-CIs_Chapter1.txt", header = TRUE, sep = "\t")

slow_trans_on_QTL<-merge(slow_trans_genes, QTL_scaffs, by = "Scaffold")
fast_trans_on_QTL<-merge(fast_trans_genes, QTL_scaffs, by = "Scaffold")

