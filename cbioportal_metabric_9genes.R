cbioMetabric_NonTNBC_gene_exp=read.table("/Users/zhaoz/Dropbox/SurvivalAnalyses_Jupyter/cbioportal_metabric_NonTNBC_GeneExp_9genes.txt",sep=",",header = T,stringsAsFactors = F)
cbioMetabric_TNBC_gene_exp=read.table("/Users/zhaoz/Dropbox/SurvivalAnalyses_Jupyter/cbioportal_metabric_TNBC_GeneExp_9genes.txt",sep=",",header = T,stringsAsFactors = F)
library(tidyr)
cbio_data_long2 <- gather(cbioMetabric_TNBC_gene_exp[,-1], GeneID, Gene_Expression, FANCL:GSK3B, factor_key=TRUE)
cbio_data_long2$METABRIC<-"TNBC"

cbio_data_long1 <- gather(cbioMetabric_NonTNBC_gene_exp[,-1], GeneID, Gene_Expression, FANCL:GSK3B, factor_key=TRUE)
cbio_data_long1$METABRIC<-"Non-TNBC"
library(ggpubr)
library(ggplot2)
cbio_dataM<-rbind(cbio_data_long1,cbio_data_long2)
cbio_dataM$METABRIC<-factor(cbio_dataM$METABRIC,levels(cbio_dataM$METABRIC)<-c("Non-TNBC", "TNBC"))


a<-ggplot(cbio_dataM,aes(x=GeneID,y=Gene_Expression,fill=METABRIC))+
  scale_fill_manual(values=c("grey","red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_boxplot(outlier.colour = "white", outlier.alpha = 0.1,outlier.shape = 1)+
  theme(axis.text=element_text(size=18,family="Times"),axis.text.x = element_text(angle=90, hjust=1),axis.title=element_text(size=24, family="Times"))+
  theme(plot.title = element_text(size=24,family="Times"))+
  theme(legend.title = element_text(colour="black", size=20,family="Times"))+
  theme(legend.text = element_text(colour="black", size = 16,family="Times"))+
  coord_cartesian(ylim = c(0,16))+
  labs(x='Genes',y='Gene expression')+
  stat_compare_means(method = "wilcox.test",size=6,label.y = c(14, 15,14, 15, 14,13, 14, 15),label='p.signif', paired = FALSE,symnum.args = list(cutpoints = c(0, 0.0001,  0.01, 1),symbols = c("**", "*", "ns")))
  #,symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", " "))
  
  #stat_compare_means(method = "wilcox.test",size=4,label.y = c(13, 15),label='p.signif', paired = FALSE)+
  ggsave("cBioportal_metabric_9genes_GE_comp.pdf") #boxplot
  #ggsave("cBioportal_metabric_2genes_GE_comp.pdf") #boxplot
a

