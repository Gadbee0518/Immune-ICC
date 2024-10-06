setwd("E:\\GSEA_IMMU\\CHOL\\CODE\\TWO") #设置工作目录
##整理metadata文件
install.packages("cli")
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("cli")
library("tidyverse")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
source("CIBERSORT.R")
source("split_violin_ggplot.R")
library("ggpubr")
library("stringr")
library("reshape2")
library("ggplot2")
library("Rmisc")
library("ggprism")


dealChol<-read.table("./DOC_PIC/dealChol.txt",check.names = F,sep = "",header = T)
CIBERSORTresult<-CIBERSORT("LM22.txt", "./DOC_PIC/dealChol.txt", perm = 10000, QN = T) #perm置换次数=1000，QN分位数归一化=TRUE
cibersort_pvalue<-CIBERSORTresult[which(CIBERSORTresult[,23]<=0.05),]
# write.csv(cibersort_pvalue,"./DOC_PIC/cib_p.csv",row.names = T,col.names = T,sep = "\t",quote = F)
cib_p<-cibersort_pvalue[,1:17]
cib_p$stage_one<-apply(cib_p[,cli_T$ID[which(cli_T$STAGE=="stage i")]],1,mean)
cib_p$stage_two<-apply(cib_p[,cli_T$ID[which(cli_T$STAGE=="stage ii")]],1,mean)
cib_p$stage_three<-cib_p[,cli_T$ID[which(cli_T$STAGE=="stage iii")]]
cib_p$stage_four<-cib_p[,cli_T$ID[which(cli_T$STAGE=="stage iv")]]
cib_p<-cib_p[,18:21]
cib_p<-data.frame(t(cib_p),check.names = F)
cib_p$STAGE<-rownames(cib_p)
rownames(cib_p)<-seq(1,4,1)
outTab=data.frame()
for(i in colnames(cib_p[,1:(ncol(cib_p)-1)])){
  rt1=cib_p[,c(i,"STAGE")]
  colnames(rt1)=c("Proportion","STAGE")
  outTab=rbind(outTab,cbind(rt1,ImmuneCell=i))
}
cib_p<-outTab
rownames(cib_p)<-seq(1,nrow(cib_p),1)
cib_p<-cib_p[-c(21:24),]
# write.csv(cib_p,"./DOC_PIC/cib_p.csv",row.names = T,col.names = T,sep = "\t",quote = F)
# cib_p<-read.csv("./DOC_PIC/cib_p.csv")
cib_p$STAGE<-factor(cib_p$STAGE,levels = c("stage_one","stage_two","stage_three","stage_four"))
sj<-cib_p[c(1:8,17:20),]
js<-cib_p[c(13:16,25:28,33:36),]
j<-cib_p[c(9:12,21:24,29:32),]
j$STAGE<-as.character(j$STAGE)
for(i in 1:length(j$STAGE)){
  if(j$STAGE[i]=="stage_one"){
    j$STAGE[i]<-"stage I"
  }
  if(j$STAGE[i]=="stage_two"){
    j$STAGE[i]<-"stage II"
  }
  if(j$STAGE[i]=="stage_three"){
    j$STAGE[i]<-"stage III"
  }
  if(j$STAGE[i]=="stage_four"){
    j$STAGE[i]<-"stage IV"
  }
}
pdf("./DOC_PIC/j.pdf",width = 6,height = 4)
ggplot(j,aes(x=STAGE,y=Proportion))+
  geom_line(aes(group=ImmuneCell,color=ImmuneCell),stat = "identity",size=1.5)+
  scale_color_manual(values = c("#eaa221","#00ee00","#cd0000"))+#j"#eaa221","#00ee00","#cd0000"
  geom_point(aes(color=ImmuneCell,shape=ImmuneCell),size=4)+#js"#ffff00","#ff5349","#8b0000"
  # scale_shape_manual(values = 0:9)+#sj "#00f5ff","#cd6600","#ff8c00"
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black",size = 1,lineend = 1),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.key = element_blank())
dev.off()

ssgseaScore<-read.csv("./DOC_PIC/ssgseaScore.txt",check.names = F,row.names = 1,header = T,sep = "\t")
ssgseaScore$stage_one<-apply(ssgseaScore[,cli_T$ID[which(cli_T$STAGE=="stage i")]],1,median)
ssgseaScore$stage_two<-apply(ssgseaScore[,cli_T$ID[which(cli_T$STAGE=="stage ii")]],1,median)
ssgseaScore$stage_three<-ssgseaScore[,cli_T$ID[which(cli_T$STAGE=="stage iii")]]
ssgseaScore$stage_four<-ssgseaScore[,cli_T$ID[which(cli_T$STAGE=="stage iv")]]
ssgseaScore<-ssgseaScore[,45:48]
ssgseaScore<-data.frame(t(ssgseaScore),check.names = F)
ssgseaScore$Macrophages<-ssgseaScore[,9]
colnames(ssgseaScore)[9]<-"Monocytes"
ssgseaScore<-ssgseaScore[,c(1:3,9,16,17)]
ssgseaScore$STAGE<-rownames(ssgseaScore)
rownames(ssgseaScore)<-seq(1,4,1)
colnames(ssgseaScore)[2]<-"T cells CD4"
colnames(ssgseaScore)[3]<-"T cells CD8"
colnames(ssgseaScore)[5]<-"Tregs"
outTab=data.frame()
for(i in colnames(ssgseaScore[,1:(ncol(ssgseaScore)-1)])){
  rt1=ssgseaScore[,c(i,"STAGE")]
  colnames(rt1)=c("pathwayScore","STAGE")
  outTab=rbind(outTab,cbind(rt1,ImmuneCell=i))
}
for(i in 1:length(outTab$STAGE)){
  if(outTab$STAGE[i]=="stage_one"){
    outTab$STAGE[i]<-"stage I"
  }
  if(outTab$STAGE[i]=="stage_two"){
    outTab$STAGE[i]<-"stage II"
  }
  if(outTab$STAGE[i]=="stage_three"){
    outTab$STAGE[i]<-"stage III"
  }
  if(outTab$STAGE[i]=="stage_four"){
    outTab$STAGE[i]<-"stage IV"
  }
}
ssgseaScore<-outTab
ssgseaScore$STAGE<-factor(ssgseaScore$STAGE,levels = c("stage I","stage II","stage III","stage IV"))
rownames(ssgseaScore)<-seq(1,nrow(ssgseaScore),1)
sjbar<-ssgseaScore[c(1:4,17:20),]
sjbar$ImmuneCell[c(1:4)]<-"B cells"
jsbar<-ssgseaScore[c(5:8,21:24),]
jbar<-ssgseaScore[c(9:16,21:24),]
pdf("./DOC_PIC/jsbar.pdf",width = 6,height = 2)
ggplot(jsbar,aes(x=STAGE,y=pathwayScore,fill=ImmuneCell))+
  geom_bar(stat = "identity",position = "dodge",width = 0.5,
           color="#3b3c36",size=0.7)+
  scale_fill_manual(values = c("#d99058","#ba55d3"))+
  #jbar "#d99058","#00ee00","#cd0000"
  #jsbar "#d99058","#ba55d3"
  #sjbar "#77b5fe","#ff8c00"
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black",size = 1,lineend = 1),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"))
dev.off()

oneDiff<-read.csv("./DOC_PIC/one_results.csv",check.names = F,header = T,row.names = 1)
twoDiff<-read.csv("./DOC_PIC/two_results.csv",check.names = F,header = T,row.names = 1)
threeDiff<-read.csv("./DOC_PIC/three_results.csv",check.names = F,header = T,row.names = 1)
fourDiff<-read.csv("./DOC_PIC/four_results.csv",check.names = F,header = T,row.names = 1)
allDiff<-read.csv("./DOC_PIC/All_results.csv",check.names = F,header = T,row.names = 1)

B<-read.table("./DOC_PIC/B.txt")
CD4Tconv<-read.table("./DOC_PIC/CD4Tconv.txt")
CD8T<-read.table("./DOC_PIC/CD8T.txt")
Treg<-read.table("./DOC_PIC/Treg.txt")
Macrophages<-read.table("./DOC_PIC/MonoMacro.txt")



library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(reshape2)
library(cli)
library(ggridges)
library(ggsci)
devtools::install_github("nicolash2/gggsea")
library("gggsea")
data<-fourDiff[rownames(fourDiff)%in%Treg$V1,]
data$SYMBOL<-rownames(data)
rownames(data)<-seq(1,nrow(data),1)
head(data)
gene = data$SYMBOL
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene = dplyr::distinct(gene,SYMBOL,.keep_all=T)
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(log2FoldChange))
geneList = data_all_sort$log2FoldChange
names(geneList) <- data_all_sort$ENTREZID 

gsea1 <- gseKEGG(geneList, organism = "hsa", pvalueCutoff = 0.05)
gsea1<- setReadable(gsea1, OrgDb=org.Hs.eg.db, keyType = 'ENTREZID')
gsea2<-gseGO(geneList, OrgDb=org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05)
gsea2<-setReadable(gsea2, OrgDb=org.Hs.eg.db, keyType = 'ENTREZID')
# dotplot(gsea)
# ridgeplot(gsea,label_format = 100)
# write.table(gsea1@result,"./DOC_PIC/GSEA/TWO_MonoMacro.txt",sep = "\t",col.names = T,row.names = T,quote = F)
pdf("./DOC_PIC/GSEA/FOUR_Tregs.pdf",width = 7,height = 6)
gseaplot2(gsea1,2,
          base_size = 15,
          rel_heights = c(0.6, 0.1, 0.3),
          pvalue_table = T,ES_geom = "line",
          color = "#ff8c00")
dev.off()
#"#77b5fe",B
#"#cd0000",CD8T
#"#ff8c00",#Tregs
#"#d99058"Macrophages
#"#ba55d3"CD4T

####tumor和normal免疫细胞丰度比较
group=sapply(strsplit(rownames(cibersort_pvalue),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
gp0<-data.frame(Group = group[which(group==0)],row.names = rownames(cibersort_pvalue[group==0,]))
gp0$Group<-"Tumor"
gp1<-data.frame(Group = group[which(group==1)],row.names = rownames(cibersort_pvalue[group==1,]))
gp1$Group<-"Normal"
gp<-rbind(gp0,gp1)
cibersort_pvalue<-as.data.frame(t(cibersort_pvalue[,-c(23:25)]))
cibersort_pvalue<-cibersort_pvalue[,rownames(gp)]
####热图17(T)+5(N)
rownames(cibersort_pvalue)[9]<-"Tregs"
gp_color<-list(
  Group=c(Tumor="red",Normal="#00cc99")
)
pdf("./DOC_PIC/22TN_heatmap.pdf",height=6,width=9)
pheatmap(cibersort_pvalue, 
         annotation_col=gp,
         annotation_colors = gp_color,
         color = colorRampPalette(colors = c("#e4d00a","grey","blue"))(100),
         cluster_cols = F,
         cluster_rows = F,
         fontsize=9,
         fontsize_row=9,
         scale="column",
         show_colnames=T,
         fontsize_col=10,
         drop_levels = T,
         cutree_cols = 2,gaps_col = c(17,22),
         cellwidth = 10, cellheight = 10,
         display_numbers = F)
dev.off()

###########排序箱式图
# mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat <- t(cibersort_pvalue) %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

a=dat%>%
  group_by(Cell_type)%>%
  summarise(m = median(Proportion))%>%
  arrange(desc(m))%>%
  pull(Cell_type)

dat$Cell_type = factor(dat$Cell_type,levels = a)

pdf("./DOC_PIC/22TN_orderboxplot.pdf",height=10,width=14)
ggplot(dat,aes(Cell_type,Proportion,fill = Cell_type)) +
  geom_boxplot(outlier.shape = 21,color = "black",size=0.7) +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold")) +
  scale_fill_manual(values = mycolor)
dev.off()

##17Tumor的排序箱型图
group=sapply(strsplit(dat$Sample,"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)

dat_T <- t(cibersort_pvalue[,1:17]) %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  gather(key = Cell_type,value = Proportion,-Sample)
a_T = dat_T %>%
  group_by(Cell_type) %>%
  summarise(m = median(Proportion)) %>%
  arrange(desc(m)) %>%
  pull(Cell_type)
dat_T$Cell_type = factor(dat_T$Cell_type,levels = a_T)
color17<-c("#8b0000",#T cells CD4 memory resting
           "#cd0000",#T cells CD8
           "#ff8c00",#Tregs
           "#ff5349",#Macrophages M2
           "#eaa221",#Macrophages M1
           "#ffff00",#Macrophages M0
           "#cdcd00",#NK cells activated
           "#cd6600",#B cells naive
           "#8b8b00",#Plasma cells
           "#008b00",#T cells follicular helper
           "#00a550",#Dendritic cells resting
           "#adff2f",#Mast cells resting
           "#bfefff",#Neutrophils
           "#00ee00",#Monocytes
           "#00f5ff",#B cells memory
           "#63b8ff",#Dendritic cells activated
           "#4f94cd",#Eosinophils
           "#4166f5",#Mast cells activated
           "#3f00ff",#NK cells resting
           "#e0b0ff",#T cells CD4 memory activated
           "#acace6",#T cells CD4 naive
           "#9370db")#T cells gamma delta
pdf("./DOC_PIC/17T_orderboxplot.pdf",height=6,width=7)
ggplot(dat_T,aes(Cell_type,Proportion,fill = Cell_type)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold")) +
  scale_fill_manual(values = color17)
dev.off()


dat_N <- t(cibersort_pvalue[,18:22]) %>% as.data.frame() %>%
  rownames_to_column("Sample") %>%
  gather(key = Cell_type,value = Proportion,-Sample)
a_N = dat_N %>%
  group_by(Cell_type) %>%
  summarise(m = median(Proportion)) %>%
  arrange(desc(m)) %>%
  pull(Cell_type)
dat_N$Cell_type = factor(dat_N$Cell_type,levels = a_N)
color5<-c("#ff5349",#Macrophages M2
          "#8b0000",#T cells CD4 memory resting
          "#cd0000",#T cells CD8
          "#cd6600",#B cells naive
          "#eaa221",#Macrophages M1
          "#8b8b00",#Plasma cells
          "#cdcd00",#NK cells activated
          "#ffff00",#Macrophages M0
          "#00ee00",#Monocytes
          "#adff2f",#Mast cells resting
          "#ff8c00",#Tregs
          "#00a550",#Dendritic cells resting
          "#008b00",#T cells follicular helper
          "#9370db",#T cells gamma delta
          "#bfefff",#Neutrophils
          "#00f5ff",#B cells memory
          "#63b8ff",#Dendritic cells activated
          "#4f94cd",#Eosinophils
          "#4166f5",#Mast cells activated
          "#3f00ff",#NK cells resting
          "#e0b0ff",#T cells CD4 memory activated
          "#acace6")#T cells CD4 naive

pdf("./DOC_PIC/5N_orderboxplot.pdf",height=6,width=7)
ggplot(dat_N,aes(Cell_type,Proportion,fill = Cell_type)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold")) +
  scale_fill_manual(values = color5)
dev.off()

######直方图
# mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
mycolor<-c("#8b0000",#T cells CD4 memory resting
           "#cd0000",#T cells CD8
           "#ff5349",#Macrophages M2
           "#cd6600",#B cells naive
           "#eaa221",#Macrophages M1
           "#ff8c00",#Tregs
           "#ffff00",#Macrophages M0
           "#cdcd00",#NK cells activated
           "#8b8b00",#Plasma cells
           "#008b00",#T cells follicular helper
           "#00a550",#Dendritic cells resting
           "#00ee00",#Monocytes
           "#adff2f",#Mast cells resting
           "#bfefff",#Neutrophils
           "#00f5ff",#B cells memory
           "#63b8ff",#Dendritic cells activated
           "#4f94cd",#Eosinophils
           "#4166f5",#Mast cells activated
           "#3f00ff",#NK cells resting
           "#e0b0ff",#T cells CD4 memory activated
           "#acace6",#T cells CD4 naive
           "#9370db")#T cells gamma delta
# dat$Cell_type<-factor(dat$Cell_type,levels = rownames(cibersort_pvalue))


dat$Sample<-factor(dat$Sample,levels = rownames(gp))
pdf("./DOC_PIC/22TN_barplot.pdf",height=8,width=11)
ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity",color="#3b3c36",size=0.15,width = 0.8) +
  labs(fill = "Cell Type",x = "",y = "Estimated Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=270,vjust = 0.55,size = 10,face = "bold"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12,face = "bold"),
        axis.text.y = element_text(size = 10,face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 12,face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = mycolor)
dev.off()      
#######添加分组显著性17(T)+5(N)
gp_color<-list(
  Group=c(Tumor="red",Normal="#00cc99")
)
dat$Group = ifelse(as.numeric(str_sub(dat$Sample,14,15))<10,"Tumor","Normal")
pdf("./DOC_PIC/22TN_comparedboxplot.pdf",height=8,width=8)
ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black",position = position_dodge(.9)) + 
  theme_bw() + 
  labs(x = "", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"))+
  scale_fill_manual(values = c(Tumor="red",Normal="#00cc99"))+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "wilcox.test")
dev.off()

#############免疫分组和免疫细胞的相关性
LMH_Immune<-read.table("./DOC_PIC/35cluster.txt",sep="\t",header = T)
LMH<-LMH_Immune[,1:2]
LMH=LMH_Immune[colnames(cibersort_pvalue[,1:17]),]
cib35<-CIBERSORTresult[match(rownames(LMH),rownames(CIBERSORTresult)),1:22]
ciberLMH=cbind(cib35,LMH)

pb<-ciberLMH[-c(10:28),-23]
pvaluetable=data.frame()
for(i in colnames(pb[,1:(ncol(pb)-1)])){
  rt1=pb[,c(i,"ImmuneGroup")]
  colnames(rt1)=c("expression","ImmuneGroup")
  temp<-kruskal.test(expression ~ ImmuneGroup, data = rt1)
  pValue=temp$p.value
  if(is.na(pValue)==TRUE){
    pValue=0
  }
  if(pValue<=0.05&pValue!=0){
    pvaluetable=rbind(pvaluetable,cbind(pValue,ImmuneCell=i))
    print(i)
    print(pValue)
    print("end")
  }
}
write.table(pvaluetable,file="./DOC_PIC/LH_ImmuneCell_pvalue.txt",sep="\t",row.names=T,quote=F,col.names = T)

outTab=data.frame()
for(i in colnames(pb[,1:(ncol(pb)-1)])){
  rt1=pb[,c(i,"ImmuneGroup")]
  colnames(rt1)=c("expression","ImmuneGroup")
  outTab=rbind(outTab,cbind(rt1,ImmuneCell=i))
}
OT<-outTab[outTab$ImmuneCell %in% pvaluetable$ImmuneCell,]

OT_summary<-summarySE(OT,measurevar = "expression",groupvars = c("ImmuneGroup","ImmuneCell"))

pdf(file="./DOC_PIC/LH_ImmuneCell_cor_boxplot.pdf",width=6,height=7) 
ggplot(data=OT, aes(x=ImmuneCell, y=expression, fill=ImmuneGroup)) +
  scale_fill_manual(values = c("blue","#e4d00a"))+##   L"#9966cc","#56a0d3"/H"#eead0e","#a2cd5a"
  stat_boxplot(size=0.75, # 线条的粗细
               width=0.4, # 误差棒的宽度，和下边箱型图的宽度设置一样宽，不然会错开
               linetype="solid", # 线条的类型，solid是实线，dashed是虚线，dotted是点线，blank是没有线
               position=position_dodge(.7), #同一个X对应的不同颜色的组别别之间的间距
               color="black" #线条颜色同一设置成黑色，不根据分组来设置
  )+
  stat_compare_means(aes(group = ImmuneGroup),
                     label = "p.format",
                     method = "wilcox.test",
                     label.y = 0.48,
                     size=5)+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.title.x = element_blank(),
        axis.line.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 1),
        axis.text.x = element_text(size = 12,face = "bold",hjust = 1,angle = 45),
        axis.line.y=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 1),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.3,"inches"))+ylim(0,0.5)
dev.off()

#################临床
CLINICAL44<-read.table("./DOC_PIC/clinical44.txt",sep="\t",header = T)
cli_T<-CLINICAL44[CLINICAL44$ID%in%colnames(cibersort_pvalue)[1:17],]
rownames(cli_T)<-seq(1,nrow(cli_T),1)
cli_N<-CLINICAL44[CLINICAL44$ID%in%colnames(cibersort_pvalue)[18:22],]
rownames(cli_N)<-seq(1,nrow(cli_N),1)
library("DESeq2")
dealCount<-read.table("./DOC_PIC/dealCount.txt",sep="\t",header = T,check.names = F)
dealCount<-dealCount[,colnames(cibersort_pvalue)]
stageOne<-dealCount[,c(cli_T$ID[which(cli_T$STAGE=="stage i")],colnames(dealCount)[18:22])]
stageOne<-stageOne[rowMeans(stageOne)>0,]
stageTwo<-dealCount[,c(cli_T$ID[which(cli_T$STAGE=="stage ii")],colnames(dealCount)[18:22])]
stageTwo<-stageTwo[rowMeans(stageTwo)>0,]
stageThree<-dealCount[,c(cli_T$ID[which(cli_T$STAGE=="stage iii")],colnames(dealCount)[18:22])]
stageThree<-stageThree[rowMeans(stageThree)>0,]
stageFour<-dealCount[,c(cli_T$ID[which(cli_T$STAGE=="stage iv")],colnames(dealCount)[18:22])]
stageFour<-stageFour[rowMeans(stageFour)>0,]
######DESeq2做差异分析
condition<-factor(c(rep("T",1),rep("N",5)))
dds<-DESeqDataSetFromMatrix(stageFour, DataFrame(condition), ~ condition)
dds<-dds[rowSums(counts(dds)) > 1, ] 
dds<-DESeq(dds)  ####log2(High/Low)
res<-results(dds)
res<-res[order(res$pvalue),]
write.csv(res,file="./DOC_PIC/four_results.csv",sep = "\t",quote = F,row.names = T)

##############################################################
###########新添加的内容，强调T细胞加强，NK细胞弱化############
##############################################################
propertion<-CIBERSORTresult[which(CIBERSORTresult[,23]<=0.05),]
propertion<-as.data.frame(propertion[,1:22])
propertion$`B cells naive`<-propertion$`B cells naive`+propertion$`B cells memory`
propertion<-propertion[,-2]
colnames(propertion)[1]<-"B cells proportion"
propertion$`T cells CD8`<-propertion$`T cells CD8`+propertion$`T cells CD4 naive`+
  propertion$`T cells CD4 memory resting`+propertion$`T cells CD4 memory activated`+
  propertion$`T cells follicular helper`+propertion$`T cells regulatory (Tregs)`+
  propertion$`T cells gamma delta`
propertion<-propertion[,-c(4:9)]
colnames(propertion)[3]<-"T cells proportion"
propertion$`NK cells resting`<-propertion$`NK cells resting`+propertion$`NK cells activated`
propertion<-propertion[,-5]
colnames(propertion)[4]<-"NK cells proportion"
propertion$`Macrophages M0`<-propertion$`Macrophages M0`+propertion$`Macrophages M1`+propertion$`Macrophages M2`
propertion<-propertion[,-c(7,8)]
colnames(propertion)[6]<-"Macrophages proportioin"

library("GSVA")
library("limma")
library("GSEABase")
library("methods")
library("ggpubr")
library("ggplot2")

gene_le_immuneSet=getGmt("./le_immune.gmt",geneIdType=SymbolIdentifier())
gene_10_immuneSet=getGmt("./10immune.gmt",geneIdType=SymbolIdentifier())
tpmbulkchol<-as.matrix(dealChol)
immune_ssgsea_Score=gsva(tpmbulkchol, gene_10_immuneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
immune_ssgsea_Score<-immune_ssgsea_Score[,rownames(propertion)]
rownames(immune_ssgsea_Score)<-c("T cells immuneScore",
                                 "B cells immuneScore",
                                 "NK cells immuneScore",
                                 "dendritic cells immuneScore",
                                 "Plasma cells immuneScore",
                                 "Macrophages immuneScore",
                                 "Eosinophils immuneScore",
                                 "Mast cells immuneScore",
                                 "Monocytes immuneScore",
                                 "Neutrophils immuneScore")

pro_immune<-cbind(propertion,t(immune_ssgsea_Score))
group=sapply(strsplit(rownames(pro_immune),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
pro_immune<-rbind(pro_immune[group==0,],pro_immune[group==1,])
pro_immune$type<-c(rep("Tumor",17),rep("Control",5))


df<-pro_immune
df<-df[1:17,c(4,15)]
df<-df[18:22,c(4,15)]

model_nk_t<-lm(`NK cells proportion`~`NK cells immuneScore`-1,data=df)
summary(model_nk_t)

model_nk_n<-lm(`NK cells proportion`~`NK cells immuneScore`-1,data=df)
summary(model_nk_n)

pdf("./DOC_PIC/NK_lineplot.pdf",height=4,width=6)
ggplot(pro_immune,aes(x=`NK cells immuneScore`,y=`NK cells proportion`,color=type,fill=type))+
  geom_point(size=4,shape=21)+
  geom_smooth(method = "lm",formula=y~x,fullrange=F,size=1.5,se=T)+
  scale_color_manual(values = c("#006b3c","#3f00ff"))+  
  scale_fill_manual(values = c("#429d75","#c4b4f3"))+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.3,"inches"))
dev.off()


df<-pro_immune
df<-df[1:17,c(3,13)]
df<-df[18:22,c(3,13)]

model_t_t<-lm(`T cells proportion`~`T cells immuneScore`-1,data=df)
summary(model_t_t)
model_t_n<-lm(`T cells proportion`~`T cells immuneScore`-1,data=df)
summary(model_t_n)


pdf("./DOC_PIC/T_lineplot.pdf",height=4,width=6)
ggplot(pro_immune,aes(x=`T cells immuneScore`,y=`T cells proportion`,color=type,fill=type))+
  geom_point(size=4,shape=21)+
  geom_smooth(method = "lm",formula=y~x,fullrange=F,size=1.5,se=T)+
  scale_color_manual(values = c("#006b3c","#ff8c00"))+  
  scale_fill_manual(values = c("#429d75","#dbaf79"))+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.3,"inches"))
dev.off()


df<-pro_immune
df<-df[1:17,c(1,14)]
df<-df[18:22,c(1,14)]

model_b_t<-lm(`B cells proportion`~`B cells immuneScore`-1,data=df)
summary(model_b_t)
model_b_n<-lm(`B cells proportion`~ `B cells immuneScore`-1,data=df)
summary(model_b_n)


pdf("./DOC_PIC/B_lineplot.pdf",height=4,width=6)
ggplot(pro_immune,aes(x=`B cells immuneScore`,y=`B cells proportion`,color=type,fill=type))+
  geom_point(size=4,shape=21)+
  geom_smooth(method = "lm",formula=y~x,fullrange=F,size=1.5,se=T)+
  scale_color_manual(values = c("#006b3c","#318ce7"))+  
  scale_fill_manual(values = c("#429d75","#89cff0"))+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.3,"inches"))
dev.off()

df<-pro_immune
df<-df[1:17,c(6,18)]
df<-df[18:22,c(6,18)]

model_m_t<-lm(`Macrophages proportion`~`Macrophages immuneScore`-1,data=df)
summary(model_m_t)
model_m_n<-lm(`Macrophages proportion`~`Macrophages immuneScore`-1,data=df)
summary(model_m_n)

pdf("./DOC_PIC/Macrophages_lineplot.pdf",height=4,width=6)
ggplot(pro_immune,aes(x=`Macrophages immuneScore`,y=`Macrophages proportion`,color=type,fill=type))+
  geom_point(size=4,shape=21)+
  geom_smooth(method = "lm",formula=y~x,fullrange=F,size=1.5,se=T)+
  scale_color_manual(values = c("#006b3c","#ff5349"))+  
  scale_fill_manual(values = c("#429d75","#f4b2aa"))+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 12,face = "bold"),
        legend.text = element_text(size = 12,face = "bold"),
        legend.title = element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.3,"inches"))
dev.off()
