########2019-01-25########
##nilo parental release###
##exitron statistics######

##packages##
library(dplyr)
library(tidyr)
library(ggpubr)

##read data##
nilo1<-read.delim("nilo1_IN_sorted.junction.exitron",header=FALSE)
nilo2<-read.delim("nilo2_IN_sorted.junction.exitron",header=FALSE)
parental1<-read.delim("parental1_IN_sorted.junction.exitron",header=FALSE)
parental2<-read.delim("parental2_IN_sorted.junction.exitron",header=FALSE)
release1<-read.delim("release1_IN_sorted.junction.exitron",header=FALSE)
release2<-read.delim("release2_IN_sorted.junction.exitron",header=FALSE)

exitron<-dplyr::data_frame(nilo1=nrow(nilo1),nilo2=nrow(nilo2),parental1=nrow(parental1),parental2=nrow(parental2),release1=nrow(release1),release2=nrow(release2))
write.csv(exitron,"k562.exitronburden.csv",sep=",",col.names = TRUE,row.names = FALSE,quote=FALSE)
exitronT<-as.data.frame(t(exitron))
exitronT$sample<-rownames(exitronT)
exitronT<-exitronT%>%rename(Exitron=V1)
exitronT$group<-c(rep("nilo",2),rep("parental",2),rep("release",2))
exitronT$group<-factor(exitronT$group,levels = c("parental","nilo","release"))
##Barplot##
p<-ggbarplot(exitronT, x = "sample", y = "Exitron",
          fill = "group",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = TRUE,     # Don't sort inside each group
          x.text.angle = 90           # Rotate vertically x axis texts
)+theme(legend.position = "none")

p1<-ggbarplot(exitronT, x = "group", y = "Exitron",
             fill = "group", 
             add = "mean_se",# change fill color by cyl
             color = "black",            # Set bar border colors to white
             palette = "jco",            # jco journal color palett. see ?ggpar
             sort.val = "desc",          # Sort the value in dscending order
             sort.by.groups = TRUE,     # Don't sort inside each group
             x.text.angle = 90           # Rotate vertically x axis texts
)+stat_compare_means(ref.group = "parental", label = "p.signif",
                     label.y = c(240, 250))+theme(legend.position = "none")
ggsave("K562.exitronburden.barplot.pdf",width=6,height=6,dpi = 300)


library(gridExtra)
pdf("K562.exitronburden.bar_2.pdf")
p2<-grid.arrange(p, p1,ncol=2) 
dev.off()


  
