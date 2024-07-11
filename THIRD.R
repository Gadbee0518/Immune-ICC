setwd("E:\\GSEA_IMMU\\CHOL\\CODE\\THREE") #设置工作目录
##整理metadata文件
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("plot3D")
library("tidyverse")
library("DESeq2")
library("rjson")
library("MEGENA")
library("qvalue")
library("fdrtool")
library("ggplot2")
library("grid")
library("ggrepel")
library("ggpubr")
library("WGCNA")
library("data.table")
library("stringr")
library("pca3d")
library("rgl")
library("FactoMineR")
library("factoextra")
library("extrafont")
library("EnhancedVolcano")
library("scales")
library("ggsci")


######读取数据并且调整Tumor和Normal的顺序
options(stringsAsFactors = F)
content<-fromJSON(file="E:/GSEA_IMMU/CHOL/DATA/metadata.cart.2022-06-29.json")
metadata<-data.frame(t(sapply(content, function(x){
  id<-x$associated_entities[[1]]$entity_submitter_id
  file_name<-x$file_name
  all<-cbind(id,file_name)
})))
rownames(metadata)<-metadata[,2]

##生成矩阵
dir<-"E:/GSEA_IMMU/CHOL/DATA/gdc_download_20220629_224748.527789/"
samples<-list.files(dir)
samplesdir<-paste0(dir,samples)
mat<-do.call(cbind,lapply(samplesdir,function(x){
  rt<-data.table::fread(x,data.table=F)
  rownames(rt)<-rt[,1]
  rt<-rt[,4]
}))

##矩阵行名列名的替换
rt<-data.table::fread(samplesdir[1],data.table = F)
colnames(mat)<-sapply(strsplit(samplesdir,'/'), '[',6)
rownames(mat)<-rt$gene_id
colnames(mat)<-metadata$X1[match(colnames(mat),metadata$X2)]
cholCount<-cbind(rt$gene_name,mat)
colnames(cholCount)[1]<-"gene_name"
cholCount<-cholCount[-c(1:4),]
write.table(cholCount,file="./DOC_PIC/cholCount.txt",quote=F,sep="\t",col.names=T)

dealCount<-cholCount
rownames(dealCount)<-cholCount[,1]
dealCount<-dealCount[,-1]
for (i in 1:length(colnames(dealCount))) {
  colnames(dealCount)[i]<-substr(colnames(dealCount)[i],1,16)
}
dealCount=avereps(dealCount)
dealCount_dim=list(rownames(dealCount),colnames(dealCount))
dealCount=matrix(as.numeric(as.matrix(dealCount)),nrow=nrow(dealCount),dimnames=dealCount_dim)
dealCount=dealCount[rowMeans(dealCount)>0,]
write.table(dealCount,file="./DOC_PIC/dealCount.txt",sep="\t",quote=F,row.names = T)
group=sapply(strsplit(colnames(dealCount),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
dealCount44<-cbind(dealCount[,group==0],dealCount[,group==1])
write.table(dealCount44,file = "./DOC_PIC/dealCount44.txt",sep = "\t",quote = F,row.names = T)
cluster.Immunity<-read.table("./DOC_PIC/cluster.Immunity.txt",sep = "",header = T,check.names = F)
HLgroup<-cluster.Immunity[-c(10:28),c(1,3)]
HLgroup<-data.frame(row.names = HLgroup$ID,HLgroup$ImmuneGroup)
colnames(HLgroup)<-"ImmuneGroup"
write.table(HLgroup,file = "./DOC_PIC/HLgroup.txt",sep = "\t",quote = F,row.names = T,col.names = T)
Count16<-dealCount44[,HLgroup$ID]
Count16=Count16[rowMeans(Count16)>0,]
write.table(Count16,file = "./DOC_PIC/Count16.txt",sep = "\t",quote = F,row.names = T)

######DESeq2做差异分析
HLgroup$ImmuneGroup=factor(HLgroup$ImmuneGroup, levels=c("Immunity_L","Immunity_H"))
dds<-DESeqDataSetFromMatrix(Count16, HLgroup, design= ~ ImmuneGroup)
dds<-dds[rowSums(counts(dds)) > 1, ] 
dds<-DESeq(dds)  ####log2(High/Low)

#####标准化
vsd1<-assay(vst(dds, blind = FALSE))#####标准化
write.table(vsd1,file="./DOC_PIC/VSD1.txt",sep = "\t",quote = F,row.names = T)
vsd2 <- assay(varianceStabilizingTransformation(dds),blind = TRUE)
write.table(vsd2,file="./DOC_PIC/VSD2.txt",sep = "\t",quote = F,row.names = T)

res<-results(dds)
res<-res[order(res$pvalue),]
# head(res)
# summary(res)
res$pvalue
write.csv(res,file="./DOC_PIC/All_results.csv",sep = "\t",quote = F,row.names = T)
UP<-subset(res, padj < 0.05 & log2FoldChange > 1)
write.csv(UP,file="./DOC_PIC/UP.csv",sep = "\t",quote = F,row.names = T)
DOWN<-subset(res, padj < 0.05 & log2FoldChange < -1)
write.csv(DOWN,file="./DOC_PIC/DOWN.csv",sep = "\t",quote = F,row.names = T)
###筛选出1815个差异基因留下四个感兴趣的模块1310个基因
Diff<-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(Diff,file="./DOC_PIC/Diff.csv",sep = "\t",quote = F,row.names = T)
blackgenes<-read.table("./DOC_PIC/WGCNA/blackwgcnaGenes.txt")
bluegenes<-read.table("./DOC_PIC/WGCNA/bluewgcnaGenes.txt")
turquoisegenes<-read.table("./DOC_PIC/WGCNA/turquoisewgcnaGenes.txt")
yellowgenes<-read.table("./DOC_PIC/WGCNA/yellowwgcnaGenes.txt")
fourmodulegenes<-rbind(blackgenes,bluegenes,turquoisegenes,yellowgenes)
deletegenes<-Diff@rownames[!(Diff@rownames%in%fourmodulegenes$V1)]
######火山图
# vals<-c('IGLV1.36',
#         'IGLV3.27',
#         'PLK5',
#         'FMO2',
#         'PLAAT5',
#         'ADAMTS8',
#         'IGLV3.12',
#         'TMPRSS6')
rr<-res
rr@listData$pvalue=rr@listData$padj
rr<-rr[!(rr@rownames%in%deletegenes),]

gp<-ifelse(rr$log2FoldChange<(-2)&rr$padj<0.05,'#0000ff',#下调
  ifelse(rr$log2FoldChange>2&rr$padj<0.05,'#ff0800',#上调
  ifelse(abs(rr$log2FoldChange)<2&rr$padj<0.05,'#d891ef',#中上
  ifelse(abs(rr$log2FoldChange)>2&rr$padj>0.05,'#F18D00','#676767'))))#下两边  中下

gp[is.na(gp)]<-'#676767'
names(gp)[gp=='#ff0800']<-'Up'
names(gp)[gp=='#d891ef']<-'p-value'
names(gp)[gp=='#F18D00']<-'Log2FC'
names(gp)[gp=='#0000ff']<-'Down'
names(gp)[gp=='#676767']<-'NS'

write.csv(rr,"./tt.csv")
tt2<-read.csv("./tt2.csv")

t2<-ifelse(tt2$log2FoldChange<(-2)&tt2$pvalue<0.05,'#0000ff',#下调
           ifelse(tt2$log2FoldChange>2&tt2$pvalue<0.05,'#ff0800',#上调
                  ifelse(abs(tt2$log2FoldChange)<2&tt2$pvalue<0.05,'#d891ef',#中上
                         ifelse(abs(tt2$log2FoldChange)>2&tt2$pvalue>0.05,'#F18D00','#676767'))))#下两边  中下

t2[is.na(t2)]<-'#676767'
names(t2)[t2=='#ff0800']<-'Up'
names(t2)[t2=='#d891ef']<-'p-value'
names(t2)[t2=='#F18D00']<-'Log2FC'
names(t2)[t2=='#0000ff']<-'Down'
names(t2)[t2=='#676767']<-'NS'



pdf("./DOC_PIC/Volcano.pdf",width = 7,height = 10)
EnhancedVolcano(rr,
                lab = rownames(rr),
                # lab = "",
                pCutoff=0.05,#y轴阈值线(水平)
                FCcutoff=2,#x轴阈值线（垂直）
                # pointSize=c(ifelse(rownames(res) %in% vals,6,4)),
                pointSize=3,#散点大小
                labSize=4,#标签大小
                boxedLabels=T,#是否在框中绘制标签
                drawConnectors=T,#是否通过连线将标签连接到对应的点上
                widthConnectors=0.2,#连线的宽度
                endsConnectors= "last",#连线绘制箭头的方向，可选first、both、last
                # selectLab=c(vals),#使用selectLab参数选定所关注的标签
                # colCustom=gp,#用group覆盖默认配色方案
                colAlpha=0.6,#调整透明度
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim=c(-10, 10),#限制X轴范围
                ylim=c(0,40),#限制Y轴范围
                cutoffLineType='longdash',#阈值线类型，可选“blank”、“solid”、“dashed”、“dotted”、“dotdash”、“longdash”和“twodash”
                cutoffLineCol='#458b00',#阈值线颜色
                cutoffLineWidth=1,#阈值线粗细
                colCustom=gp#用group覆盖默认配色方案
                )+theme(axis.line.x = element_line(size = 1),
                        axis.title = element_text(size=20,face = "bold"),
                        axis.line.y = element_line(size = 1),
                        legend.text = element_text(size = 20))
dev.off()

WGCNA_GENES<-rownames(vsd1[rowSums(vsd1)>=44,])
WGCNA_DAT2<-vsd2[rownames(vsd2)%in%WGCNA_GENES,]
write.table(WGCNA_DAT2,file = "./DOC_PIC/WGCNA_DAT2_16.txt",sep = "\t",quote = F,row.names = T)
WGCNA_DAT1<-vsd1[rownames(vsd1)%in%WGCNA_GENES,]
write.table(WGCNA_DAT1,file = "./DOC_PIC/WGCNA_DAT1_16.txt",sep = "\t",quote = F,row.names = T)

diff_name<-c(upgenes,downgenes)
Tpm16<-read.table("./DOC_PIC/TPM44.txt",sep = "\t",row.names = 1,header = T,check.names = F)
Tpm16<-Tpm16[rownames(Tpm16)%in%diff_name,]
Tpm16<-Tpm16[,rownames(HLgroup)]
write.table(Tpm16,file="./DOC_PIC/DIFF_TPM16.txt",sep = "\t",quote = F,row.names = T,col.names = T)

allgeneTpm16<-read.table("./DOC_PIC/TPM44.txt",sep = "\t",row.names = 1,header = T,check.names = F)
allgeneTpm16<-allgeneTpm16[,rownames(HLgroup)]
ALL1310<-allgeneTpm16[rownames(allgeneTpm16)%in%fourmodulegenes$V1,]
write.table(rownames(ALL1310),file = "./DOC_PIC/1310diffGenenames.txt",sep = "\t",quote = F,row.names = F)
#####处理一下数据，选取35例患者（免疫低中高），一会做WGCNA与免疫低中高结合
# dealCount35<-dealCount[,group==0]
# dealCount35<-dealCount35[rowMeans(dealCount35)>0,]
# write.table(dealCount35,file = "./DOC_PIC/dealCount35.txt",sep = "\t",quote = F,row.names = T)
# ALL1815<-allgeneTpm16[rowMeans(ALL1815)>0,]
ALL1310.pca<-prcomp(t(ALL1310),scale. = TRUE)
ALL1310.pca_sum=summary(ALL1310.pca)
fviz_eig(ALL1310.pca, addlabels = T, ylim = c(0, 100))


pca3d(ALL1310.pca, 
      components = 1:3, 
      group = HLgroup$ImmuneGroup, 
      show.centroids=TRUE, 
      show.group.labels =F,
      col = hl,
      radius = 1.1) 
snapshotPCA3d(file="good1.png")


allgeneTpm16<-allgeneTpm16[rowMeans(allgeneTpm16)>0,]
allgeneTpm16.pca<-prcomp(t(allgeneTpm16),scale. = TRUE)
allgeneTpm16.pca_sum=summary(allgeneTpm16.pca)
pdf(file="./DOC_PIC/badpcabarPlot.pdf",height = 5,width=7)
fviz_eig(allgeneTpm16.pca,addlabels = T, ylim = c(0, 50))+
  theme(axis.title = element_text(size = 15,face = "bold"),
        axis.text = element_text(size = 15,face = "bold"))
dev.off()


pca3d(allgeneTpm16.pca, 
      components = 1:3, 
      group = HLgroup$ImmuneGroup, 
      show.centroids=TRUE, 
      show.group.labels =F,
      col = hl,
      radius = 1.1)
snapshotPCA3d(file="bad1.png")
# L_color<-rep("#e4d00a",9)
# names(L_color)<-cluster.Immunity$ID[1:9]
# H_color<-rep("blue",7)
# names(H_color)<-cluster.Immunity$ID[29:35]
# hl<-c(L_color,H_color)
# plot3d(allgeneTpm16.pca$x[,1:3],
#        xlab = "PC1",ylab = "PC2",zlab = "PC3",
#        col = hl,
#        type = "s",
#        size = 1,
#        lwd = 2,box=T)


Tpm16.pca<-prcomp(t(Tpm16),scale. = TRUE)
write.table(chol16.pca$rotation,file="./DOC_PIC/PC_gene.txt",quote=F,sep="\t")
write.table(predict(chol16.pca),file="./DOC_PIC/PC_sample.txt",quote=F,sep="\t")
Tpm16.pca_sum=summary(Tpm16.pca)
#输出PC比重
write.table(chol16.pca_sum$importance,file="./DOC_PIC/PC_importance.txt",quote=F,sep="\t")
# 
pdf(file="./DOC_PIC/pcabarPlot.pdf",height = 5,width=7)
fviz_eig(Tpm16.pca, addlabels = T, ylim = c(0, 50))+
  theme(axis.title = element_text(size = 15,face = "bold"),
        axis.text = element_text(size = 15,face = "bold"))
dev.off()

pca3d(Tpm16.pca, components = 1:3, group = HLgroup$ImmuneGroup, show.centroids=TRUE, show.group.labels =F)
# plot3d(Tpm16.pca$x[,1:3],
#        xlab = "PC1",ylab = "PC2",zlab = "PC3",
#        col = hl,
#        type = "s",
#        size = 1,
#        lwd = 10,box=T)
# grid3d(c("x","y","z"))
# rgl.postscript("good1.pdf", fmt = "pdf", drawText = F)


fviz_pca_ind(Tpm16.pca,
              geom.ind = "point", # show points only (nbut not "text")
              col.ind = HLgroup$ImmuneGroup, # color by groups
              palette = "jco", 
              addEllipses = TRUE, # Concentration ellipses
              legend.title = "Groups")

fviz_pca_ind(Tpm16.pca,
              geom.ind = "point", # show points only (or "text")
              col.ind = HLgroup$ImmuneGroup, # color by groups
              palette = "lancet",
              addEllipses = FALSE, # Concentration ellipses
              legend.title = "group",
              ggtheme = theme_minimal(base_size=18,base_family =fontfam),title="PCA(all phase)")


###GO通路富集
up1310<-as.vector(rownames(UP)[rownames(UP)%in%rownames(ALL1310)])
down1310<-as.vector(rownames(DOWN)[rownames(DOWN)%in%rownames(ALL1310)])
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("dplyr")
library("stringr")

#GO上调富集分析
entrezIDs <- as.character(mget(up1310, org.Hs.egSYMBOL2EG, ifnotfound=NA))    #找出基因对应的id
up_ID=cbind(up1310,entrezID=entrezIDs)
write.table(up_ID,file="./DOC_PIC/1310up_mark2id.txt",sep="\t",quote=F,row.names=F)    #输出结果
up_ID=read.table("./DOC_PIC/1310up_mark2id.txt",sep="\t",header=T,check.names=F)
up_ID<-up_ID[is.na(up_ID[,"entrezID"])==F,]#去除基因id为NA的基因
no.na.genes=up_ID$entrezID
upGO <- enrichGO(gene = no.na.genes,
               OrgDb = org.Hs.eg.db,
               pvalueCutoff =0.05,
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(upGO,file="./DOC_PIC/1310upGO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果
GOUP30<-read.csv("./DOC_PIC/GO_KEGG/GOUP30.csv",header = T)
GOUP30$ONTOLOGY<-factor(GOUP30$ONTOLOGY,levels = c("BP","CC","MF"))
GOUP30<-GOUP30%>%
  arrange(ONTOLOGY,Count)
GOUP30$Description<-factor(GOUP30$Description,levels = GOUP30$Description)
pdf("./DOC_PIC/1310upGOplot.pdf",width = 7.5,height = 6)
ggplot(GOUP30, aes(x = Description, y = Count)) +
  geom_bar(stat = "identity", aes(fill = ONTOLOGY), alpha = 2,width = 0.85) +
  scale_fill_manual(values=c("#006b3c","#8b4500","#5a4fcf"))+
  facet_grid(ONTOLOGY ~ ., scales = "free", space = "free", margins = F) + 
  coord_flip() +
  theme_light() +
  theme(axis.text.x = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        title = element_text(size = 14,face = "bold")) +
  labs(y = "Number of Genes", x = "") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black",fill = "#4682b4"),
        strip.text = element_text(size=12,colour = "white",face = "bold")) +
  labs(title = "GO-up barplot")
dev.off()
#GO下调富集分析
entrezIDs <- as.character(mget(down1310, org.Hs.egSYMBOL2EG, ifnotfound=NA))    #找出基因对应的id
down_ID=cbind(down1310,entrezID=entrezIDs)
write.table(down_ID,file="./DOC_PIC/1310down_mark2id.txt",sep="\t",quote=F,row.names=F)    #输出结果
down_ID=read.table("./DOC_PIC/1310down_mark2id.txt",sep="\t",header=T,check.names=F)
down_ID<-down_ID[is.na(down_ID[,"entrezID"])==F,]#去除基因id为NA的基因
no.na.genes=down_ID$entrezID
downGO <- enrichGO(gene = no.na.genes,
                 OrgDb = org.Hs.eg.db,
                 pvalueCutoff =0.1,
                 qvalueCutoff = 0.1,
                 ont="all",
                 readable =T)
write.table(downGO,file="./DOC_PIC/1310downGO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果
GODOWN30<-read.csv("./DOC_PIC/GO_KEGG/GODOWN30.csv",header = T)
GODOWN30$ONTOLOGY<-factor(GODOWN30$ONTOLOGY,levels = c("BP","CC","MF"))
GODOWN30<-GODOWN30%>%
  arrange(ONTOLOGY,Count)
GODOWN30$Description<-factor(GODOWN30$Description,levels = GODOWN30$Description)
pdf("./DOC_PIC/1310downGOplot.pdf",width = 9,height = 5)
ggplot(GODOWN30, aes(x = Description, y = Count)) +
  geom_bar(stat = "identity", aes(fill = ONTOLOGY), alpha = 2,width = 0.85) +
  scale_fill_manual(values=c("#00a86b","#ffb90f","#b57edc"))+
  facet_grid(ONTOLOGY ~ ., scales = "free", space = "free", margins = F) + 
  coord_flip() +
  theme_light() +
  theme(axis.text.x = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        title = element_text(size = 14,face = "bold"),
        panel.border = element_rect(size = 1)) +
  labs(y = "Number of Genes",x="") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "black",fill = "#63b8ff"),
        strip.text = element_text(size=12,colour = "black",face = "bold")) +
  labs(title = "GO-down barplot")
dev.off()

####KEGG上调富集分析
up_kegg <- enrichKEGG(gene = up_ID$entrezID, organism = "hsa", pvalueCutoff =0.1, qvalueCutoff =0.1)   #富集分析
UPKEGG=as.data.frame(up_kegg)
UPKEGG$geneID=as.character(sapply(UPKEGG$geneID,function(x)paste(up_ID$up1310[match(strsplit(x,"/")[[1]],as.character(up_ID$entrezID))],collapse="/")))
write.table(UPKEGG,file="./DOC_PIC/1310UPKEGG.txt",sep="\t",quote=F,row.names = F)
###KEGG下调富集分析
down_kegg <- enrichKEGG(gene = down_ID$entrezID, organism = "hsa", pvalueCutoff =0.5, qvalueCutoff =0.5)   #富集分析
DOWNKEGG=as.data.frame(down_kegg)
DOWNKEGG$geneID=as.character(sapply(DOWNKEGG$geneID,function(x)paste(down_ID$down1310[match(strsplit(x,"/")[[1]],as.character(down_ID$entrezID))],collapse="/")))
write.table(DOWNKEGG,file="./DOC_PIC/1310DOWNKEGG.txt",sep="\t",quote=F,row.names = F)

KEGG30<-read.csv("./DOC_PIC/GO_KEGG/KEGG30.csv",header = T)
KEGG30$LogP<-(log10(KEGG30$p.adjust))
KEGG30$LogP[1:15]<-abs(KEGG30$LogP[1:15])
KEGG30$Group<-c(rep("UP",15),rep("DOWN",11))
pdf("./DOC_PIC/1310UPvsDOWN_KEGGplot.pdf",width = 9,height = 6)
ggplot(KEGG30, aes(x = reorder(Description,LogP), LogP, fill=Group)) + 
  geom_bar(stat = 'identity',alpha = 1,width = 0.9) + 
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(),
        panel.border = element_rect(colour = "#b3b3b3"),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        title = element_text(size = 13,face = "bold"))+
  theme(panel.border = element_rect(size = 1.1))+
  scale_y_continuous(breaks=seq(-4, 24, 2))+
  labs(x = "", y="LogP",title = "KEGG UPvsDOWN plot")+
  scale_fill_manual(values = c("#548b54","#000080"))#设置颜色
dev.off()


# upgenes=as.vector(rownames(UP))
# entrezIDs <- as.character(mget(upgenes, org.Hs.egSYMBOL2EG, ifnotfound=NA))    #找出基因对应的id
# up_ID=cbind(upgenes,entrezID=entrezIDs)
# write.table(up_ID,file="./DOC_PIC/up_mark2id.txt",sep="\t",quote=F,row.names=F)    #输出结果
# #GO上调富集分析
# up_ID=read.table("./DOC_PIC/up_mark2id.txt",sep="\t",header=T,check.names=F)
# up_ID<-up_ID[is.na(up_ID[,"entrezID"])==F,]#去除基因id为NA的基因
# no.na.genes=up_ID$entrezID
# upGO <- enrichGO(gene = no.na.genes,
#                OrgDb = org.Hs.eg.db, 
#                pvalueCutoff =0.01, 
#                qvalueCutoff = 0.01,
#                ont="all",
#                readable =T)
# write.table(upGO,file="./DOC_PIC/upGO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果
# pdf("./DOC_PIC/upGOplot.pdf",width = 25,height = 25)
# plot1<-barplot(upGO, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# plot2<-dotplot(upGO,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# CombinePlots(plots = list(plot1,plot2))
# dev.off()
# ####KEGG上调富集分析
# up_kegg <- enrichKEGG(gene = up_ID$entrezID, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #富集分析
# UPKEGG=as.data.frame(up_kegg)
# UPKEGG$geneID=as.character(sapply(UPKEGG$geneID,function(x)paste(up_ID$upgenes[match(strsplit(x,"/")[[1]],as.character(up_ID$entrezID))],collapse="/")))
# write.table(UPKEGG,file="./DOC_PIC/UPKEGG.txt",sep="\t",quote=F,row.names = F)                          #保存富集结果
# pdf(file="./DOC_PIC/upKEGGplot.pdf",width = 25,height = 25)
# plot1<-barplot(up_kegg, drop = TRUE, showCategory = 30)
# plot2<-dotplot(up_kegg, showCategory = 30,orderBy = "GeneRatio")
# CombinePlots(plots = list(plot1,plot2))
# dev.off()
# 
# downgenes=as.vector(rownames(DOWN))
# entrezIDs <- as.character(mget(downgenes, org.Hs.egSYMBOL2EG, ifnotfound=NA))    #找出基因对应的id
# down_ID=cbind(downgenes,entrezID=entrezIDs)
# write.table(down_ID,file="./DOC_PIC/down_mark2id.txt",sep="\t",quote=F,row.names=F)    #输出结果
# #GO下调富集分析
# down_ID=read.table("./DOC_PIC/down_mark2id.txt",sep="\t",header=T,check.names=F)
# down_ID<-down_ID[is.na(down_ID[,"entrezID"])==F,]#去除基因id为NA的基因
# no.na.genes=down_ID$entrezID
# downGO <- enrichGO(gene = no.na.genes,
#                  OrgDb = org.Hs.eg.db, 
#                  pvalueCutoff =0.05, 
#                  qvalueCutoff = 0.05,
#                  ont="all",
#                  readable =T)
# write.table(downGO,file="./DOC_PIC/downGO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果
# pdf("./DOC_PIC/downGOplot.pdf",width = 25,height = 25)
# plot1<-barplot(downGO, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# plot2<-dotplot(downGO,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
# CombinePlots(plots = list(plot1,plot2))
# dev.off()
# 
# down_kegg <- enrichKEGG(gene = down_ID$entrezID, organism = "hsa", pvalueCutoff =0.2, qvalueCutoff =0.2)   #富集分析
# DOWNKEGG=as.data.frame(down_kegg)
# DOWNKEGG$geneID=as.character(sapply(DOWNKEGG$geneID,function(x)paste(down_ID$downgenes[match(strsplit(x,"/")[[1]],as.character(down_ID$entrezID))],collapse="/")))
# write.table(DOWNKEGG,file="./DOC_PIC/DOWNKEGG.txt",sep="\t",quote=F,row.names = F)                          #保存富集结果
# pdf(file="./DOC_PIC/downKEGGplot.pdf",width = 25,height = 15)
# plot1<-barplot(down_kegg, drop = TRUE, showCategory = 30)
# plot2<-dotplot(down_kegg, showCategory = 30,orderBy = "GeneRatio")
# CombinePlots(plots = list(plot1,plot2))
# dev.off()
# 
# #GO富集上下调
# library(readxl)
# library(dplyr)
# library(ggplot2)
# 
# # 整理数据格式
# upGO_data <- read.table("./DOC_PIC/upGO.txt",sep="\t",header=T,check.names=F)
# downGO_data <- read.table("./DOC_PIC/downGO.txt",sep="\t",header=T,check.names=F)
# 
# upGO_data <- upGO_data[,c(3,8,10)]
# upGO_data$change <- c("up")
# downGO_data <- downGO_data[,c(3,8,10)]
# downGO_data$change <- c("down")
# 
# upGO_vs_downGO <- rbind(upGO_data[1:15,], downGO_data[1:15,])
# # 改变数值的正负号，方便绘图时展示在左右两边，更加美观
# data <- upGO_vs_downGO %>% 
#   mutate(Qvalue = ifelse(change == "up", -log10(qvalue), log10(qvalue))) %>% 
#   arrange(change,Qvalue) 
# data <- data %>% 
#   mutate(Count = ifelse(change == "up", Count, -Count))
# data$Description = factor(data$Description, levels = unique(data$Description),ordered = T)
# head(data)
# # 设置一下坐标的范围
# tmp <- with(data, labeling::extended(-range(Qvalue)[2], range(Qvalue)[2], m = 7))
# lm <- tmp[c(1,length(tmp))]
# # 绘图
# pdf("./DOC_PIC/upGO_vs_downGO.pdf",width = 30,height = 10)
# ggplot(data, aes(x=Description, y=Qvalue)) +
#   geom_segment( aes(x=Description, xend=Description, y=0, yend=Qvalue, color=change), size=5, alpha=0.9) +
#   theme_light() +
#   theme(
#     panel.border = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     axis.title.x = element_text(size = 15),
#     # legend.text = element_text(size = 20),
#     # legend.position = c(0.9, 0.5),
#     # 不显示图例
#     legend.position = "none"
#   ) +
#   xlab("") +
#   ylab("-log10(Qvalue)")+
#   ylim(lm)+
#   scale_y_continuous(breaks = tmp,
#                      labels = abs(tmp),
#                      # 双Y轴
#                      sec.axis = sec_axis(~./0.3,
#                                          name = 'Number of genes',
#                                          breaks = seq(-100,100,10),
#                                          labels = abs(seq(-100,100,10))))+
#   # 翻转一下
#   coord_flip() + 
#   # 绘制折线图
#   geom_line(aes(Description, Count*0.3, group = 1), size=0.8, color='#71ad46') + geom_point(aes(Description, Count*0.3, group = 1), color='#71ad46')
# 
# dev.off()
# 
# ###KEGG富集上下调
# upKEGG_data <- read.table("./DOC_PIC/UPKEGG.txt",sep="\t",header=T,check.names=F)
# downKEGG_data <- read.table("./DOC_PIC/DOWNKEGG.txt",sep="\t",header=T,check.names=F)
# 
# upKEGG_data <- upKEGG_data[,c(2,7,9)]
# upKEGG_data$change <- c("up")
# downKEGG_data <- downKEGG_data[,c(2,7,9)]
# downKEGG_data$change <- c("down")
# 
# upKEGG_vs_downKEGG <- rbind(upKEGG_data[1:15,], downKEGG_data[1:15,])
# # 改变数值的正负号，方便绘图时展示在左右两边，更加美观
# data <- upKEGG_vs_downKEGG %>% 
#   mutate(Qvalue = ifelse(change == "up", -log10(qvalue), log10(qvalue))) %>% 
#   arrange(change,Qvalue) 
# data <- data %>% 
#   mutate(Count = ifelse(change == "up", Count, -Count))
# data$Description = factor(data$Description, levels = unique(data$Description),ordered = T)
# head(data)
# # 设置一下坐标的范围
# tmp <- with(data, labeling::extended(-range(Qvalue)[2], range(Qvalue)[2], m = 7))
# lm <- tmp[c(1,length(tmp))]
# # 绘图
# pdf("./DOC_PIC/upKEGG_vs_downKEGG.pdf",width = 30,height = 10)
# ggplot(data, aes(x=Description, y=Qvalue)) +
#   geom_segment( aes(x=Description, xend=Description, y=0, yend=Qvalue, color=change), size=5, alpha=0.9) +
#   theme_light() +
#   theme(
#     panel.border = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     axis.title.x = element_text(size = 15),
#     # legend.text = element_text(size = 20),
#     # legend.position = c(0.9, 0.5),
#     # 不显示图例
#     legend.position = "none"
#   ) +
#   xlab("") +
#   ylab("-log10(Qvalue)")+
#   ylim(lm)+
#   scale_y_continuous(breaks = tmp,
#                      labels = abs(tmp),
#                      # 双Y轴
#                      sec.axis = sec_axis(~./0.3,
#                                          name = 'Number of genes',
#                                          breaks = seq(-100,100,10),
#                                          labels = abs(seq(-100,100,10))))+
#   # 翻转一下
#   coord_flip() + 
#   # 绘制折线图
#   geom_line(aes(Description, Count*0.3, group = 1), size=0.8, color='#71ad46') + geom_point(aes(Description, Count*0.3, group = 1), color='#71ad46')
# 
# dev.off()
# Generate some example plot

