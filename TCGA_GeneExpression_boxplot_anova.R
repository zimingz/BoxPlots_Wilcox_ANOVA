tcga_Nontnbc_gene_exp=read.table("/Users/zhaoz/Dropbox/SurvivalAnalyses_Jupyter/TCGA_nonTNBC_9genes_976pts.csv",sep=",",header = T,stringsAsFactors = F)
tcga_tnbc_gene_exp=read.table("/Users/zhaoz/Dropbox/SurvivalAnalyses_Jupyter/TCGA_TNBC_9genes_115pts.csv",sep=",",header = T,stringsAsFactors = F)
tcga_normal=read.table("/Users/zhaoz/Dropbox/SurvivalAnalyses_Jupyter/Normal_9genes.csv",sep=",",header = T,stringsAsFactors = F)
library(tidyr)
#data2 <- gather(tcga_tnbc_gene_exp[,-1], GeneID, Gene_Expression, CCNE1:TPX2, factor_key=TRUE)
data2 <- gather(tcga_tnbc_gene_exp[,-1], GeneID, Gene_Expression, FANCL:GSK3B, factor_key=TRUE)
data2$Gene_Expression <- log2(data2$Gene_Expression+1)
data2$TCGA<-"TNBC"

#data1 <- gather(tcga_Nontnbc_gene_exp[,-1], GeneID, Gene_Expression, CCNE1:TPX2, factor_key=TRUE)
data1 <- gather(tcga_Nontnbc_gene_exp[,-1], GeneID, Gene_Expression, FANCL:GSK3B, factor_key=TRUE)
data1$Gene_Expression <- log2(data1$Gene_Expression+1)
data1$TCGA<-"Non_TNBC"

#data0 <- gather(tcga_normal[,-1], GeneID, Gene_Expression, CCNE1:TPX2, factor_key=TRUE)
data0 <- gather(tcga_normal[,-1], GeneID, Gene_Expression, FANCL:GSK3B, factor_key=TRUE)
data0$Gene_Expression <- log2(data0$Gene_Expression+1)
data0$TCGA<-"Normal"

data<-rbind(data1,data0,data2)
data$TCGA<-factor(data$TCGA,levels(data$TCGA)<-c("Normal", "Non_TNBC", "TNBC"))
library(ggpubr)

library(ggplot2)

positions <- c("Normal", "Non_TNBC", "TNBC")
ggplot(data,aes(GeneID,Gene_Expression,fill=TCGA))+
  scale_fill_manual(values=c("white","grey", "red"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_boxplot(outlier.colour = "white", outlier.alpha = 0.1,outlier.shape = 1)+
  theme(axis.text=element_text(size=18,family="Times"),axis.text.x = element_text(angle=90, hjust=1),axis.title=element_text(size=24, family="Times"))+
  theme(plot.title = element_text(size=24,family="Times"))+
  theme(legend.title = element_text(colour="black", size=20,family="Times"))+
  theme(legend.text = element_text(colour="black", size = 16,family="Times"))+
  coord_cartesian(ylim = c(0, 16))+
  #labs(x='Genes',y='Gene expression',title='ANOVA TEST')+
  labs(x='Genes',y='Gene expression')+
  #stat_compare_means(method = "anova",size=4,label.y = c(13, 15),label='p.signif', paired = FALSE)+
  #ggsave("tcga_GeneExpression_Comp_BoxPlot_anova_2genes_withNormal_log2.pdf") #boxplot
  #stat_compare_means(method = "anova",size=8,label.y = c(14, 15,14, 15, 14,13, 14, 15),label='p.signif', paired = FALSE)+
  stat_compare_means(method = "anova",size=6,label.y = c(14, 15,14, 15, 14,13, 14, 15),label='p.signif', paired = FALSE,symnum.args = list(cutpoints = c(0, 0.0001,  0.01, 1),symbols = c("**", "*", "ns")))
  ggsave("tcga_GeneExpression_Comp_BoxPlot_anova_9genes_withNormal_log2.pdf") #boxplot

#write.csv(data, file="TCGA_2Genes_withNormal_log2.csv", sep = "\t")
write.csv(data, file="TCGA_9Genes_withNormal_log2.csv", sep = "\t")



