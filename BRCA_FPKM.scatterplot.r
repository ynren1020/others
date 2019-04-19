##################################2019-04-19########################################################
##/home/tywang/Projects/PanCanAtlas/Gene_Expression/HTSEQ/protein_coding_genes/FPKM/BRCA.FPKM.txt###
##FABP4vsDNMT3A,vsDNMT3B,DNMT1######################################################################
####################################################################################################

library(tidyverse)
library(ggpubr)

brca<-read.delim("BRCA.FPKM.txt",stringsAsFactors = FALSE)
dim(brca) #19703  1062
genes<-c("FABP4","DNMT1","DNMT3A","DNMT3B")
brcasub<-brca[brca$GeneID%in%genes,]
brcasubT<-as.data.frame(t(brcasub))
brcasubT<-apply(brcasubT,2,as.character)
colnames(brcasubT)<-brcasubT[1,]
brcasubT<-brcasubT[-1,]
brcasubT<-as.data.frame(apply(brcasubT,2,as.numeric))

##plot##
p1<-ggscatter(brcasubT, x = "FABP4", y = "DNMT1",
              xlab = "FABP4 (FPKM)",
              ylab = "DNMT1 (FPKM)"
          #add = "reg.line",                                 # Add regression line
          #conf.int = TRUE,                                  # Add confidence interval
          #add.params = list(color = "blue",
          #                  fill = "lightgray")
)+
  stat_cor(method = "pearson", label.x = 3, label.y = 20)  # Add correlation coefficient



