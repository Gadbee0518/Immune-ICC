setwd("E:/GSEA_IMMU/CHOL/CODE/FIVE1")

install.packages("Seurat")
install.packages("outliers")
install.packages("pbmcapply")
install.packages("doFuture")
install.packages("rio")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("singscore")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("ggvenn")
BiocManager::install("cowplot")
BiocManager::install("paletteer")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scCATCH")
BiocManager::install("tidyverse")
BiocManager::install("mindr")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scater")
BiocManager::install("pheatmap")
BiocManager::install("psych")
BiocManager::install("monocle")
BiocManager::install("openxlsx")

library(Seurat)
library(outliers)
library(pbmcapply)
library(doFuture)
library(singscore)
library(GSVA)
library(GSEABase)
library(limma)
library(tidyverse)
library(Matrix)
library(ggsci)#改变配色
library(dplyr)
library(magrittr)
library(patchwork)
library(ggplot2)
library(mindr)

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggrepel)
##做四种亚型的差异基因图
load("./prediff_scCholdata.RData")

B14_DIFF<-FindMarkers(CHOL_data,ident.1 = 14,min.pct = 0.25,
                  logfc.threshold = 0,pseudocount.use=0.01)
B15_DIFF<-FindMarkers(CHOL_data,ident.1 = 15,min.pct = 0.25,
                 logfc.threshold = 0,pseudocount.use=0.01)
T1_DIFF<-FindMarkers(CHOL_data,ident.1 = 1,min.pct = 0.25,
                 logfc.threshold = 0,pseudocount.use=0.01)
M712_DIFF<-FindMarkers(CHOL_data,ident.1 = c(7,12),min.pct = 0.25,
                   logfc.threshold = 0,pseudocount.use=0.01)

FM<-T1_DIFF
FM<-FM %>% mutate(Difference= pct.1 - pct.2)
FM$group=0
FM$gene=rownames(FM)
FM$Difference<-FM$Difference*100
for(i in 1:nrow(FM)){
  if(FM$avg_log2FC[i]>=0.5 & FM$Difference[i]>=25 ){
    FM$group[i]='up'
  }else if(FM$avg_log2FC[i]<=(-0.5) & FM$Difference[i]<=(-25) ){
    FM$group[i]='down'
  }else{
    FM$group[i]='nodiff'
  }
}
pdf(file="./DOC_PIC/pipT1_diff.pdf",width=7,height=7)
ggplot(FM,aes(x=Difference,y=avg_log2FC))+
  geom_point(size=2.2,aes(color=group))+
  labs(title="MacrophagesLp subtype")+
  scale_color_manual(values = c('#0000ff','#9ac0cd','#ff0800'))+
  # geom_label_repel(data = subset(FM,group!='nodiff'),
  #                  aes(label=gene),
  #                  segment.size=0.3,size=3)+
  theme_classic()+
  theme(axis.text = element_text(size = 17,face = "bold"),
        axis.title = element_text(size = 20,face='bold'),
        axis.line = element_line(linetype = 1,size = 1),
        plot.title=element_text(hjust=0.5,size=22,face = "bold"),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 17))+
  geom_vline(xintercept = 0,linetype=2,size=1.2,color='#006400')+
  geom_hline(yintercept = 0,linetype=2,size=1.2,color='#006400')
#,arrow = arrow(length = unit(0.3, 'cm'))
dev.off()
# newmarkers<-FindAllMarkers(CHOL_data)


B14<-signewmarkers[signewmarkers$cluster==14,]
ego_B14<-enrichGO(gene = intersect(gene1310$V1,B14_gene[-which(B14_gene%in%c(T1_gene,B15_gene,M712_gene))]),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
write.csv(ego_B14,"./DOC_PIC/ego_B14.csv",row.names = T)

B15<-signewmarkers[signewmarkers$cluster==15,]
ego_B15<-enrichGO(gene = intersect(gene1310$V1,B15_gene[-which(B15_gene%in%c(B14_gene,T1_gene,M712_gene))]),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
write.csv(ego_B15,"./DOC_PIC/ego_B15.csv",row.names = T)

T1<-signewmarkers[signewmarkers$cluster==1,]
ego_T1<-enrichGO(gene = intersect(gene1310$V1,T1_gene[-which(T1_gene%in%c(B14_gene,B15_gene,M712_gene))]),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
write.csv(ego_T1,"./DOC_PIC/ego_T1.csv",row.names = T)

M712<-signewmarkers[signewmarkers$cluster==c(7,12),]
ego_M712<-enrichGO(gene = intersect(gene1310$V1,M712_gene[-which(M712_gene%in%c(B14_gene,B15_gene,T1_gene))]),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
iyuewrite.csv(ego_M712,"./DOC_PIC/ego_M712.csv",row.names = T)

newfoursub<-read.csv("./new four sub results.csv",check.names = F)
newfoursub$Description<-factor(newfoursub$Description,levels = unique(newfoursub$Description) )
for (i in 1:length(foursub$Group)) {
  if(foursub$Group[i]=="B14"){
    foursub$Group[i]<-"BLp"
  }
  if(foursub$Group[i]=="B15"){
    foursub$Group[i]<-"BHn"
  }
  if(foursub$Group[i]=="M712"){
    foursub$Group[i]<-"MacrophagesLp"
  }
  if(foursub$Group[i]=="T1"){
    foursub$Group[i]<-"THn"
  }
}
library(forcats)
pdf(file="./DOC_PIC/newfoursub.pdf",width=8,height=8)
ggplot(newfoursub,aes(Group,Description))+
  geom_point(aes(color=p.adjust,size=GeneRatio))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        legend.text = element_text(size=10,face = "bold"),
        legend.title = element_text(size = 11,face = "bold"))+
  scale_color_gradient(low = '#fcc200',high = '#50c878')+
  labs(x=NULL,y=NULL)
dev.off()
##寻找差异表达的特征
logFCfilter=0.5
adjPvalFilter=0.05
CHOL_data.markers <- FindAllMarkers(object = CHOL_data,
                                    only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)

scMarkers=CHOL_data.markers[(abs(as.numeric(as.vector(CHOL_data.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(CHOL_data.markers$p_val_adj))<adjPvalFilter),]
write.table(scMarkers,file="./DOC_PIC/markers.xls",sep="\t",row.names=F,quote=F)

library(ggvenn)
library(ggtext)
library(eulerr)
library(VennDiagram)
library(rio)

gene1310<-read.table("./DOC_PIC/1310diffGenenames.txt")
B14_gene<-rownames(B14_DIFF)[B14_DIFF$group!='nodiff']
B15_gene<-rownames(B15_DIFF)[B15_DIFF$group!='nodiff']
T1_gene<-rownames(T1_DIFF)[T1_DIFF$group!='nodiff']
M712_gene<-rownames(M712_DIFF)[M712_DIFF$group!='nodiff']

list(A=gene1310$V1,B=c(T1_gene,B14_gene,B15_gene,M712_gene)) %>% 
  ggvenn(show_percentage = T,show_elements = F,label_sep = ",",
         digits = 1,stroke_color = "white",
         fill_color = c("#E41A1C", "#1E90FF"),
         set_name_color = c("#E41A1C", "#1E90FF"))
intersect(gene1310$V1,rownames(sig))
