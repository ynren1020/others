####Jan9,2019########
####Coverage plot####

library(reshape)
library(ggplot2)
library(dplyr)

CASC9.bed<-cbind(c("chr8","chr8"),c())

nilo1_IN_CASC9=read.table("nilo1_IN_CASC9.bam.sort.coverage", sep="\t", header=F)
nilo1_IN_CASC9=rename(nilo1_IN_CASC9,chr1=V1) # renames the header
nilo1_IN_CASC9=rename(nilo1_IN_CASC9,location=V2)
nilo1_IN_CASC9=rename(nilo1_IN_CASC9,depth=V3)
nilo1_IN_CASC9$sample<-rep("IN",nrow(nilo1_IN_CASC9))



p<-ggplot(nilo1_IN_CASC9, aes(x=V2, y=V3)) +
  geom_line(colour="red", size=1)
  #scale_y_continuous(trans = scales::log10_trans(), breaks = scales::trans_breaks("log10", function(x) 10^x))