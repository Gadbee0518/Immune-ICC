setwd("E:\\GSEA_IMMU\\CHOL\\CODE\\ONE") #设置工作目录
##整理metadata文件
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("tidyr")
library("tidyverse")
library("factoextra")
library("rjson")
library("factoextra")
library("gplots")
library("data.table")
library("dplyr")

library("ggplot2")
library("GenVisR")
library("descriptr")
library("vioplot")
library("pheatmap")
library("RColorBrewer")
library("tidyr")
# update.packages("tidyr")
library("viridis")
library("corrplot")
source("CIBERSORT.R")
require("data.table")
library("GSVA")
library("limma")
library("GSEABase")
library("methods")
library("ggpubr")
library("stringr")
library("survival")
library("survminer")

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
  rt<-rt[,7]
}))

##矩阵行名列名的替换
rt<-data.table::fread(samplesdir[1],data.table = F)
colnames(mat)<-sapply(strsplit(samplesdir,'/'), '[',6)
rownames(mat)<-rt$gene_id
colnames(mat)<-metadata$X1[match(colnames(mat),metadata$X2)]
CHOL<-cbind(rt$gene_name,mat)
colnames(CHOL)[1]<-"gene_name"
CHOL<-CHOL[-c(1:4),]
write.table(CHOL,file="./DOC_PIC/CHOL.txt",quote=F,sep="\t",col.names=T)

dealChol<-CHOL
rownames(dealChol)<-CHOL[,1]
dealChol<-dealChol[,-1]
for (i in 1:length(colnames(dealChol))) {
  colnames(dealChol)[i]<-substr(colnames(dealChol)[i],1,16)
}
dealChol=avereps(dealChol)
dealChol_dim=list(rownames(dealChol),colnames(dealChol))
dealChol=matrix(as.numeric(as.matrix(dealChol)),nrow=nrow(dealChol),dimnames=dealChol_dim)
dealChol=dealChol[rowMeans(dealChol)>0,]
write.table(dealChol,file="./DOC_PIC/dealChol.txt",sep="\t",quote=F,row.names = T)  

chol.44<-read.table("./DOC_PIC/dealChol.txt",header = T,sep = '\t',check.names = F)
group=sapply(strsplit(colnames(chol.44),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
chol.44<-cbind(chol.44[,group==0],chol.44[,group==1])
write.table(chol.44,file="./DOC_PIC/chol44add.txt",sep="\t",quote=F,row.names = T)
####ssGSEA分析#######
dimnames=list(rownames(chol.44),colnames(chol.44))
chol.44=matrix(as.numeric(as.matrix(chol.44)),nrow=nrow(chol.44),dimnames=dimnames)
geneSet=getGmt("le_immune.gmt",geneIdType=SymbolIdentifier())
ssgseaScore=gsva(chol.44, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
write.table(ssgseaScore,file="./DOC_PIC/ssgseaScore.txt",sep="\t",quote=F,col.names = T,row.names = T)
##B79/CD4Tconv121/CD8T68/CD8Tex72/DC79/Endothelial3/Fibroblasts12/Mast51/Mono+Macro131/Myofibroblasts6/
##Neutrophils60/NK78/pDC81/Plasma49/TMKI6754/Treg64
library(sparcl) #引用包
#删掉正常样品
group=sapply(strsplit(colnames(ssgseaScore),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
data=ssgseaScore[,group==0]
hc = hclust(dist(t(data)),method = "ward.D")#欧式距离，离差平方和法
HCTREE=cutree(hc,3)#将聚类结果分成3组
#输出聚类表格
write.table(HCTREE,file="./DOC_PIC/HCTREE.txt",sep="\t",quote=F,col.names=F,row.names = T)
#输出聚类图片
pdf(file="./DOC_PIC/35treeCluster.pdf",width=11,height=11)
fviz_dend(hc, k = 3, # 聚类的类别数目为2
          cex = 1.2, # 数据标签的字体大小
          k_colors = c("#e4d00a", "blue", "#458b74"),
          color_labels_by_k = T, # 数据标签也根据颜色设定
          rect_border = c("#fbec5d", "blue", "#9acd32"),
          rect = TRUE, # 使用不同颜色的矩形框标定类别
          rect_fill = TRUE,
          labels_track_height=5,
          show_labels = T)
dev.off()

####绘制热图
library(pheatmap)#三组字体的颜色("#e4d00a", "blue", "#458b74")
#背景颜色("#fbec5d", "blue", "#9acd32")
Type=read.table("./DOC_PIC/HCTREE.txt",sep="\t",check.names=F,header=F)
Type[,2]=paste0("Cluster",Type[,2])
Type=Type[order(Type[,2]),]
rt=ssgseaScore[,as.vector(Type[,1])]

cluster=as.data.frame(Type[,2])
row.names(cluster)=Type[,1]
colnames(cluster)="Cluster"

an_color<-list(Cluster=c(Cluster1="#fdee00",Cluster2="#9acd32",Cluster3="blue"))
pdf("./DOC_PIC/heatmap.pdf",height=4,width=8)
pheatmap(rt, annotation=cluster,
         annotation_colors = an_color,
         color = colorRampPalette(c("#ffff31","blue"))(100),
         cluster_rows = T,
         cluster_cols =F,
         fontsize=8,
         fontsize_row=7,
         scale="row",
         show_colnames=T,
         fontsize_col=5)
dev.off()
#定义每个cluster分组信息，需要修改
Cluster1="Immunity_L"
Cluster2="Immunity_M"
Cluster3="Immunity_H"
#输出分组结果
a=c()
a[Type[,2]=="Cluster1"]=Cluster1
a[Type[,2]=="Cluster2"]=Cluster2
a[Type[,2]=="Cluster3"]=Cluster3
clusterOut=cbind(Type,a)
colnames(clusterOut)<-c("ID","Cluster","ImmuneGroup")
row.names(clusterOut)<-seq(1,35,1)
write.table(clusterOut,file="./DOC_PIC/cluster.Immunity.txt",sep="\t",quote=F,col.names=T,row.names=F)

#####肿瘤微环境
library(limma)
library(estimate)
#删除正常，只保留肿瘤样品,输出整理后的矩阵文件
write.table(chol.44[,1:35],file="./DOC_PIC/35tumorsymbol.txt",sep="\t",quote=F,row.names = T,col.names=T)
#运行estimate包
filterCommonGenes(input.f = "./DOC_PIC/35tumorsymbol.txt", 
                  output.f = "./DOC_PIC/commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "./DOC_PIC/commonGenes.gct",
              output.ds="./DOC_PIC/estimateScore.gct",
              platform = "illumina")

#输出每个样品的打分
estimateScores=read.table("./DOC_PIC/estimateScore.gct",skip = 2,header = T,check.names = F)
rownames(estimateScores)=estimateScores[,1]
estimateScores=t(estimateScores[,3:ncol(estimateScores)])
rownames(estimateScores)=gsub("\\.","\\-",rownames(estimateScores))
write.table(estimateScores,file="./DOC_PIC/estimateScores.txt",sep="\t",quote=F,col.names=T,row.names = T)

######绘制肿瘤微环境的热图
#读取免疫分组文件
rownames(clusterOut)<-clusterOut$ID
clusterOut<-clusterOut[,-1]
cluster=cbind(clusterOut,estimateScores)
write.table(cluster,file="./DOC_PIC/35cluster.txt",sep="\t",quote=F,col.names=T,row.names = T)
#绘制热图
ann_colors<-list(
  # Cluster=c(Cluster1="#fdee00",Cluster2="green",Cluster3="blue"),
  ImmuneGroup=c(Immunity_L="#e4d00a",Immunity_M="#458b74",Immunity_H="blue"),
  StromalScore=colorRampPalette(c("Tan1","Tan4"))(1000),
  ImmuneScore=colorRampPalette(c("Cyan1","Cyan4"))(1000),
  ESTIMATEScore=colorRampPalette(c("Plum1","Plum4"))(1000)
)
pdf("./DOC_PIC/35estimateLMH_heatmap.pdf",height=10,width=12)
#三组字体的颜色("#e4d00a", "blue", "#458b74")

cluster1<-cluster[,-1]
rownames(ssgseaScore)[9]<-"Macrophages"
pm<-pheatmap(ssgseaScore[,rownames(cluster1)],
         annotation = cluster1,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#e4d00a","blue"))(100),
         cluster_cols =F,
         fontsize=9,
         cellheight = 10,
         cellwidth = 10,
         fontsize_row=9,
         scale="row",
         show_colnames=T,
         fontsize_col=10)
dev.off()
############免疫分组内部相关性
library(GGally)
# cluster$Cluster<-factor(cluster$Cluster)
cluster1$ImmuneGroup<-factor(cluster1$ImmuneGroup)
pdf("./DOC_PIC/pipLMH_ggally.pdf",height=10,width=10)
ggpairs(cluster1[,2:4],columns = 1:3,
        upper = list(continuous = wrap("cor", size = 5)),
        ggplot2::aes(color=cluster1$ImmuneGroup))+
  scale_color_manual(values = c("blue","#e4d00a","#458b74"))+
  scale_fill_manual(values = c("blue","#e4d00a","#458b74"))+
  theme(axis.text = element_text(colour = "black",size = 10,face = "bold"),
        strip.text = element_text(colour = "black",size = 15,face = "bold"))
dev.off()

##########绘制小提琴图免疫分组与肿瘤微环境的相关性
#三组字体的颜色("#e4d00a", "blue", "#458b74")
clustercopy=cluster
clustercopy<-clustercopy[cluster$Cluster!="Cluster2",]
dc=c("Immunity_L","Immunity_H")
immunescore_pvalue<-t.test(ImmuneScore ~ ImmuneGroup,data = clustercopy)

pdf(file="./DOC_PIC/LHvioplot.pdf",width=6,height=6)
ggviolin(clustercopy, x="ImmuneGroup",
         legend="none",
         y="ImmuneScore",
         fill = "ImmuneGroup",
         palette = c("#e4d00a","blue"),
         add = "boxplot", add.params = list(fill = "white")) +
  # geom_signif(comparisons = list(c("Immunity_L","Immunity_H")),#设置需要比较的组
  #        test = t.test, ##计算方法
  #        textsize = 6,
  #        y_position = 5000,#图中横线位置设置
  #        tip_length = c(0.02,0.02),#横线下方的竖线设置
  #        size=0.8,color="black")+
  stat_compare_means(method = "t.test", label = "p.format",y_position = 5000,
                     label.x.npc ="center", size = 6,comparisons = list(c("Immunity_L","Immunity_H")))+
  theme(axis.line.x=element_line(color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1),
        axis.title.x = element_text(size = 15,face = "bold"),
        axis.line.y=element_line(color="black",size=1),
        axis.ticks.y=element_line(color="black",size=1),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text = element_text(colour = "black",size = 15,face = "bold"))+ylim(-2000,6000)
dev.off()

###########提取HLA相关基因
HLA<-chol.44[which(grepl("HLA-",row.names(chol.44))==TRUE),]
HLA<-HLA[,rownames(clustercopy)]
HLA_data=data.frame()
for(i in 1:nrow(HLA)){
  HLA_data=rbind(HLA_data,cbind(expression=log2(HLA[i,]+1),gene=row.names(HLA)[i],ImmuneGroup=as.vector(clustercopy[,2])))
}

write.table(HLA_data,file="./DOC_PIC/HLA16S_39G.txt",sep="\t",row.names=F,quote=F,col.names = T)

#绘制箱型图
hla=read.table("./DOC_PIC/HLA16S_39G.txt",sep="\t",header=T,check.names=F)       #读取箱线图输入文件
hla$ImmuneGroup=factor(hla$ImmuneGroup, levels=c("Immunity_L","Immunity_H"))
p=ggboxplot(hla, x="gene", y="expression", 
            color = "ImmuneGroup",
            orientation = "horizontal",
            ylab="Gene expression",
            xlab="",
            palette = c("#e4d00a","blue"))
# p=p+rotate_x_text(0)

pdf(file="./DOC_PIC/HLA_boxplot.pdf",width=6,height=10)#输出图片文件
p+stat_compare_means(aes(group=ImmuneGroup),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*","#")),label = "p.signif",
                     label.y = 15,label.x.npc = "center",size=4)+
  theme(axis.line.x=element_line(color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1),
        axis.title.x = element_text(size = 12,face = "bold"),
        axis.line.y=element_line(color="black",size=1),
        axis.ticks.y=element_line(color="black",size=1),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text = element_text(colour = "black",size = 12,face = "bold"),
        legend.key.size = unit(0.3,"inches"),
        legend.text = element_text(size = 12,face = "bold"),
        title = element_text(size = 12,face = "bold"))
dev.off()


LH_ImmuneGenes<-data.frame(t(data[,clinical_LH16$ID]))
LH_ImmuneGenes<-cbind(LH_ImmuneGenes,clinical_LH16$ImmuneGroup)
colnames(LH_ImmuneGenes)[17]<-"ImmuneGroup"
colnames(LH_ImmuneGenes)[9]<-"Macrophages"
outTab=data.frame()
for(i in colnames(LH_ImmuneGenes)[1:16]){
  rt1=LH_ImmuneGenes[,c(i,"ImmuneGroup")]
  colnames(rt1)=c("expression","ImmuneGroup")
  outTab=rbind(outTab,cbind(rt1,ImmuneCell=i))
}


pdf(file="./DOC_PIC/LH_ImmuneGenes_cor_boxplot.pdf",width=9.6,height=7.2) 
ggplot(data=outTab, aes(x=ImmuneCell, y=expression, fill=ImmuneGroup)) +
  scale_fill_manual(values = c("blue","#e4d00a"))+##   L"#9966cc","#56a0d3"/H"#eead0e","#a2cd5a"
  stat_boxplot(size=0.75, # 线条的粗细
               width=0.5, # 误差棒的宽度，和下边箱型图的宽度设置一样宽，不然会错开
               linetype="solid", # 线条的类型，solid是实线，dashed是虚线，dotted是点线，blank是没有线
               position=position_dodge(.7), #同一个X对应的不同颜色的组别别之间的间距
               color="black" #线条颜色同一设置成黑色，不根据分组来设置
  )+
  stat_compare_means(aes(group = ImmuneGroup),
                     label = "p",
                     method = "t.test",
                     label.y = c(0.92,0.87,1,1,
                                 0.83,1.2,1.1,0.8,
                                 0.65,1.1,0.75,0.9,
                                 0.8,0.8,1,1),
                     label.x.npc = "center",size=4)+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.line.x=element_line(color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1),
        axis.title.x = element_blank(),
        axis.line.y=element_line(color="black",size=1),
        axis.ticks.y=element_line(color="black",size=1),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(colour = "black",size = 12,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 12,face = "bold",angle = 45,hjust = 1),
        legend.key.size = unit(0.3,"inches"),
        legend.text = element_text(size = 12,face = "bold"),
        title = element_text(size = 12,face = "bold")
        ) + ylim(0,1.2)
dev.off()


##################生存分析
###############免疫分组和生存相关分析
####提取临床数据  1-Dead   0-Alive
CLINICAL<-read.csv("E:/GSEA_IMMU/CHOL/DATA/TCGA-CHOL.GDC_phenotype.csv",header = T)
for (i in 1:nrow(CLINICAL)) {
  if(CLINICAL$STATUS[i]=="Dead"){
    CLINICAL$STATUS[i]<-1
  }else{
    CLINICAL$STATUS[i]<-0
  }
}
##提取整体44的临床数据  0-Alive   1-Dead
clinical.44<-CLINICAL[match(colnames(chol.44),CLINICAL$ID),]
rownames(clinical.44)<-seq(1,nrow(clinical.44),1)
clinical.44$STATUS<-as.numeric(clinical.44$STATUS)
clinical.44$TIME<-clinical.44$TIME/30
clinical.44$STAGE[which(sapply(strsplit(clinical.44$STAGE,""),"[",8)=="v")]<-"stage iv"
write.table(clinical.44,file="./DOC_PIC/clinical44.txt",sep="\t",quote=F,col.names = T)

clinical_LH16<-clinical.44[match(rownames(clustercopy),clinical.44$ID),]
rownames(clinical_LH16)<-seq(1,16,1)
clinical_LH16<-cbind(clinical_LH16,clustercopy)
clinical_LH16$TIME<-clinical_LH16$TIME*30
write.table(clinical_LH16,file = "./DOC_PIC/clinical_LH16.txt",sep="\t",quote=F,col.names = T)
# fit<-survfit(Surv(TIME,STATUS) ~ ImmuneGroup, data = clinical_LH16)
# pdf("./DOC_PIC/16ImmuneLH_SUR.pdf",height=12,width=8)
# ggsurvplot(fit,                     # survfit object with calculated statistics.
#            pval = TRUE,             # show p-value of log-rank test.
#            conf.int = TRUE,         # show confidence intervals for 
#            # point estimaes of survival curves.
#            fun = "event",
#            conf.int.style = "step",  # customize style of confidence intervals
#            xlab = "Time in months",   # customize X axis label.
#            break.time.by = 10,     # break X axis in time intervals by 200.
#            ggtheme = theme_light(), # customize plot and risk table with a theme.
#            # risk.table.col="strata",
#            risk.table = "abs_pct",  # absolute number and percentage at risk.
#            risk.table.y.text.col = T,# colour risk table text annotations.
#            risk.table.y.text = F,# show bars instead of names in text annotations
#            # in legend of risk table.
#            ncensor.plot = T,      # plot the number of censored subjects at time t
#            surv.median.line = "hv" # add the median survival pointer.
#            #legend.labs =    c("Male", "Female"),    # change legend labels.
#            #palette =     c("#E7B800", "#2E9FDF") # custom color palettes.
#            
# )
# dev.off()

#查看五年生存率
summary(fit)

clinical_LMH35<-clinical.44[match(rownames(cluster),clinical.44$ID),]
rownames(clinical_LMH35)<-seq(1,35,1)
clinical_LMH35<-cbind(clinical_LMH35,cluster)
fit<-survfit(Surv(TIME,STATUS) ~ ImmuneGroup, data = clinical_LMH35)
pdf("./DOC_PIC/35ImmuneLMH_SUR.pdf",height=12,width=8)
ggsurvplot(fit,                     # survfit object with calculated statistics.
           pval = F,             # show p-value of log-rank test.
           conf.int = TRUE,         # show confidence intervals for 
           # point estimaes of survival curves.
           # fun = "event",
           conf.int.style = "step",  # customize style of confidence intervals
           xlab = "Time in months",   # customize X axis label.
           break.time.by = 10,     # break X axis in time intervals by 200.
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           risk.table.col="strata",
           risk.table = "abs_pct",  # absolute number and percentage at risk.
           risk.table.y.text.col = T,# colour risk table text annotations.
           risk.table.y.text = F,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = T,      # plot the number of censored subjects at time t
           surv.median.line = "hv",  # add the median survival pointer.
           # add.all = T,
           palette = "hue"# custom color palettes.
)
dev.off()

clinical_LH16<-clinical_LMH35[-c(10:28),]
rownames(clinical_LH16)<-seq(1,nrow(clinical_LH16),1)
fit16<-survfit(Surv(TIME,STATUS) ~ ImmuneGroup, data = clinical_LH16)
pdf("./DOC_PIC/16ImmuneLH_SUR.pdf",height=12,width=8)
ggsurvplot(fit16,                     # survfit object with calculated statistics.
           pval = F,             # show p-value of log-rank test.
           conf.int = TRUE,         # show confidence intervals for 
           # point estimaes of survival curves.
           # fun = "event",
           conf.int.style = "step",  # customize style of confidence intervals
           xlab = "Time in months",   # customize X axis label.
           break.time.by = 10,     # break X axis in time intervals by 200.
           ggtheme = theme_light(), # customize plot and risk table with a theme.
           risk.table.col="strata",
           risk.table = "abs_pct",  # absolute number and percentage at risk.
           risk.table.y.text.col = T,# colour risk table text annotations.
           risk.table.y.text = F,# show bars instead of names in text annotations
           # in legend of risk table.
           ncensor.plot = T,      # plot the number of censored subjects at time t
           surv.median.line = "hv",  # add the median survival pointer.
           # add.all = T,
           palette = "hue"# custom color palettes.
)
dev.off()

N<-chol.44[,36:44]
L<-chol.44[,rownames(cluster[as.character(cluster$Cluster)%in%"Cluster1",])]
H<-chol.44[,rownames(cluster[as.character(cluster$Cluster)%in%"Cluster3",])]
save(N,L,H,file = "./NLH.RData")
