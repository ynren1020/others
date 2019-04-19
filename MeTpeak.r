#############################################Jan8,2019#############################################################
####Methlyation peaks in nilo1_IP.bam,nilo2_IP.bam,parental1_IP.bam,parental2_IP.bam(sorted,get from run_IP.sh)####
####Ref:Shujun Liu's Cell research paper:A dynamic N6-methyladenosine methylome regulates intrinsic and acquired### 
#########################################resistance to tyrosine kinase inhibitors##################################
####################################################################################################################

##packages##
library("devtools")
source("https://bioconductor.org/biocLite.R")
install_github("compgenomics/MeTPeak")
library(MeTPeak)

##modify gtf data align with the example gtf file##
hg38gtf<-read.delim("hg38_testsort.gtf",header = FALSE)



# in the real case, change the gtf to what you need
gtf <- file.path('hg38_testsort.gtf')

ip1 <- file.path('nilo1_IP_sorted.bam')
ip2 <- file.path('nilo2_IP_sorted.bam')
#ip3 <- system.file('IP3.bam',package='MeTPeak')
input1 <- file.path('nilo1_IN_sorted.bam')
input2 <- file.path('nilo2_IN_sorted.bam')
#input3 <- system.file('Input3.bam',package='MeTPeak')

IP_BAM <- c(ip1,ip2)
INPUT_BAM <- c(input1,input2)

K562_2vs1<-metpeak(GENE_ANNO_GTF=gtf,IP_BAM = IP_BAM,INPUT_BAM = INPUT_BAM,
                   EXPERIMENT_NAME="K562")

write.table(K562_2vs1,"K562_2vs1.txt",sep="\t")