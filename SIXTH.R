setwd("E:/GSEA_IMMU/CHOL0406/CODE/SIX")

B_cells<-read.csv("./DOC_PIC/CELLTYPE/B_cells.csv",header = T,check.names = F)
colnames(B_cells)[2]<-"celltype"
DC_cells<-read.csv("./DOC_PIC/CELLTYPE/DC.csv",header = T,check.names = F)
colnames(DC_cells)[2]<-"celltype"
Endothelial_cells<-read.csv("./DOC_PIC/CELLTYPE/Endothelial_cells.csv",header = T,check.names = F)
colnames(Endothelial_cells)[2]<-"celltype"
Epithelial_cells<-read.csv("./DOC_PIC/CELLTYPE/Epithelial_cells.csv",header = T,check.names = F)
colnames(Epithelial_cells)[2]<-"celltype"
Fibroblasts_cells<-read.csv("./DOC_PIC/CELLTYPE/Fibroblasts.csv",header = T,check.names = F)
colnames(Fibroblasts_cells)[2]<-"celltype"
Hepatocytes_cells<-read.csv("./DOC_PIC/CELLTYPE/Hepatocytes.csv",header = T,check.names = F)
colnames(Hepatocytes_cells)[2]<-"celltype"
Macrophage_cells<-read.csv("./DOC_PIC/CELLTYPE/Macrophage.csv",header = T,check.names = F)
colnames(Macrophage_cells)[2]<-"celltype"
Monocyte_cells<-read.csv("./DOC_PIC/CELLTYPE/Monocyte.csv",header = T,check.names = F)
colnames(Monocyte_cells)[2]<-"celltype"
Neutrophils_cells<-read.csv("./DOC_PIC/CELLTYPE/Neutrophils.csv",header = T,check.names = F)
colnames(Neutrophils_cells)[2]<-"celltype"
NK_cells<-read.csv("./DOC_PIC/CELLTYPE/NK_cell.csv",header = T,check.names = F)
colnames(NK_cells)[2]<-"celltype"
Smooth_muscle_cells<-read.csv("./DOC_PIC/CELLTYPE/Smooth_muscle_cells.csv",header = T,check.names = F)
colnames(Smooth_muscle_cells)[2]<-"celltype"
T_cells<-read.csv("./DOC_PIC/CELLTYPE/T_cells.csv",header = T,check.names = F)
colnames(T_cells)[2]<-"celltype"
Tissue_stem_cells<-read.csv("./DOC_PIC/CELLTYPE/Tissue_stem_cells.csv",header = T,check.names = F)
colnames(Tissue_stem_cells)[2]<-"celltype"
######################################################################################
######################################################################################
H_pos<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_H_pos.csv",row.names = 1,header = T,check.names = F)
H_neg<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_H_neg.csv",row.names = 1,header = T,check.names = F)
L_pos<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_L_pos.csv",row.names = 1,header = T,check.names = F)
L_neg<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_L_neg.csv",row.names = 1,header = T,check.names = F)
ssgseaScore<-read.table("./DOC_PIC/ssgseaScore.txt",header = T,row.names = 1,check.names = F)
ssgseaScore_H<-ssgseaScore[,rownames(H_neg)]
ssgseaScore_L<-ssgseaScore[,rownames(L_neg)]

HpEpithelial<-H_pos[,intersect(Epithelial_cells$Barcode,colnames(H_pos))]#33.965%上皮细胞占比
# HpFibroblasts<-H_pos[,intersect(Fibroblasts_cells$Barcode,colnames(H_pos))]#2.624%纤维母细胞占比
HpHepatocytes<-H_pos[,intersect(Hepatocytes_cells$Barcode,colnames(H_pos))]#57.726%肝细胞占比
# HpSmooth<-H_pos[,intersect(Smooth_muscle_cells$Barcode,colnames(H_pos))]#3.353%平滑肌细胞占比

HnT<-H_neg[,intersect(T_cells$Barcode,colnames(H_neg))]#24.554%T细胞占比
HnB<-H_neg[,intersect(B_cells$Barcode,colnames(H_neg))]#68.75%B细胞占比
# HnNK<-H_neg[,intersect(NK_cells$Barcode,colnames(H_neg))]#5.357%NK细胞占比

LpB<-L_pos[,intersect(B_cells$Barcode,colnames(L_pos))]#7.463%B细胞占比
LpMacrophage<-L_pos[,intersect(Macrophage_cells$Barcode,colnames(L_pos))]#23.507%巨噬细胞占比
LpMonocyte<-L_pos[,intersect(Monocyte_cells$Barcode,colnames(L_pos))]#3.358%单核白细胞占比
LpHepatocytes<-L_pos[,intersect(Hepatocytes_cells$Barcode,colnames(L_pos))]#63.06%肝细胞占比

# LnT<-L_neg[,intersect(T_cells$Barcode,colnames(L_neg))]#1.651%T细胞占比
# LnEndothelial<-L_neg[,intersect(Endothelial_cells$Barcode,colnames(L_neg))]#0.661%内皮细胞占比
LnEpithelial<-L_neg[,intersect(Epithelial_cells$Barcode,colnames(L_neg))]#72.853%上皮细胞占比
LnHepatocytes<-L_neg[,intersect(Hepatocytes_cells$Barcode,colnames(L_neg))]#24.371%肝细胞占比

write.csv(HnB,"./DOC_PIC/HnB.csv",col.names = T)
write.csv(LpB,"./DOC_PIC/LpB.csv",col.names = T)
write.csv(HnT,"./DOC_PIC/HnT.csv",col.names = T)
write.csv(LpMonocyte,"./DOC_PIC/LpMonocyte.csv",col.names = T)
write.csv(LpMacrophage,"./DOC_PIC/LpMacrophage.csv",col.names = T)

rownames(ssgseaScore_L)[9]<-"Macrophages"
rownames(ssgseaScore_H)[9]<-"Macrophages"
###############################################################################################
##############重要回归图
#############################################################################################
library(tidyverse)
library(corrplot)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("performance")
library(SignifReg)
library(broom)
library(performance)
#royalblue-neg   indianred1-pos
df<-data.frame(BHn=rowMeans(HnB),HnB)
ssg<-cbind(t(ssgseaScore_H[1,]),df)[,1:2]
colnames(ssg)[1]<-"ImmuneScore"


ssg %>%
  select_if(is.numeric) %>%
  cor() %>%
  corrplot(method = 'color',order = 'AOE',addCoef.col = 'grey')

pdf(file="./DOC_PIC/RESULT/TPM_NEW/HnB.pdf",width=6,height=4)
ssg %>%
  # select(B,where(is.numeric),-`AAAGCAATCAACACTG-1`) %>%
  #"#77b5fe",B
  #"#cd0000",CD8T
  #"#ff8c00",#Tregs
  #"#d99058"Macrophages
  #"#ba55d3"CD4T
  pivot_longer(cols = -ImmuneScore,
               names_to = 'variable',
               values_to = 'value') %>%
  ggplot(aes(x=value,y=ImmuneScore)) +
  geom_point(size=6,color="#77b5fe",shape=20) + #b385ff-NK/#f07f4d-B/#e936a7-CD4Tconv/#2243b6-Monocyte
                                                  #7-CD8T/3-CD8Tex/1-TMKI67/2-Treg/#03c03c-Fibroblasts/
                                                  #00ff00-Myofibroblasts/#9400d3-Endothelial/#00bfff-Macrophage
  geom_smooth(method = 'lm',color="royalblue",size=2) +
  facet_wrap(~variable, scales = 'free')+
  theme(axis.title = element_text(size = 16,face = "bold"),
        axis.text = element_text(size = 14,face = "bold"),
        strip.text = element_text(size = 20,face = "bold"))
dev.off()

pdf(file="./DOC_PIC/RESULT/TPM_NEW/HnT.pdf",width=6,height=4)

ssg %>%
  # select(B,where(is.numeric),-`AAAGCAATCAACACTG-1`) %>%
  pivot_longer(cols = -ImmuneScore,
               names_to = 'variable',
               values_to = 'value') %>%
  ggplot(aes(x=value,y=ImmuneScore)) +
  geom_point(size=6,color="#cd0000",shape=20) +
  geom_smooth(method = 'lm',formula = y~poly(x,2),color="royalblue",size=2) +
  facet_wrap(~variable,scales = "free")+
  theme(axis.title = element_text(size = 16,face = "bold"),
        axis.text = element_text(size = 14,face = "bold"),
        strip.text = element_text(size = 20,face = "bold"))
dev.off()
# llmm<-lm(formula = MoMa~poly(mean,3),data=ssg)


model1 = lm(Macrophages~Macrophages_Cor,data=ssg)
model2 = lm(Macrophages~poly(Macrophages_Cor,4),data=ssg)
model3 = lm(Macrophages~poly(Macrophages_Cor,5),data=ssg)
model4 = step(model3,direction = 'both',trace = 0)
model5 = SignifReg(model3,criterion = 'p-value',alpha = 0.05,direction = 'both')


####模型对比
sum_table = data.frame(Model = c('model1','model2','model3','model4','model5'),
                       adj.R2 = c(summary(model1)$adj.r.squared,
                                  summary(model2)$adj.r.squared,
                                  summary(model3)$adj.r.squared,
                                  summary(model4)$adj.r.squared,
                                  summary(model5)$adj.r.squared),
                       AIC = AIC(model1,model2,model3,model4,model5)[,2],
                       BIC = BIC(model1,model2,model3,model4,model5)[,2])

tidy(model1)
tidy(model2)
tidy(model3)
tidy(model4)
tidy(model5)

check_model(model1)
check_model(model2)
check_model(model3)
check_model(model4)
check_model(model5)

#################################################################################
##################重要箱型图
#################################################################################
library(ggpubr)
library(ggplot2)
#
TPM_H_neg<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_H_neg.csv",row.names = 1,header = T,check.names = F)
TPM_H_pos<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_H_pos.csv",row.names = 1,header = T,check.names = F)
TPM_L_neg<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_L_neg.csv",row.names = 1,header = T,check.names = F)
TPM_L_pos<-read.csv("./DOC_PIC/RESULT/TPM_NEW/TPM_L_pos.csv",row.names = 1,header = T,check.names = F)
#
tabPos=data.frame()
for(i in colnames(TPM_H_pos)){
  rtpos=TPM_H_pos[,i]
  tabPos=rbind(tabPos,cbind(rtpos,Group="Poor survival cells",Immune="High_Immune"))
}
colnames(tabPos)[1]<-"corScore"

tabNeg=data.frame()
for(j in colnames(TPM_H_neg)){
  rtneg=TPM_H_neg[,j]
  tabNeg=rbind(tabNeg,cbind(rtneg,Group="Good survival cells",Immune="High_Immune"))
}
colnames(tabNeg)[1]<-"corScore"

TPMH<-rbind(tabPos,tabNeg)
write.csv(TPMH,"./DOC_PIC/RESULT/TPM_NEW/TPMH.csv",col.names = T)

TPMH$Group=factor(TPMH$Group, levels=c("Poor survival cells","Good survival cells"))
TPMH$corScore<-as.numeric(TPMH$corScore)
H_pvalue<-wilcox.test(corScore ~ Group,data = TPMH)

pdf(file="./DOC_PIC/RESULT/TPM_NEW/Hpos_vs_Hneg_boxplot.pdf",width=5,height=7)#输出图片文件
# par(mar=c(1,1,1,1))
ggplot(data=TPMH, aes(x=Group, y=corScore, fill=Group)) +
  scale_fill_manual(values = c("indianred1","royalblue"))+##   L"#9966cc","#56a0d3"/H"#eead0e","#a2cd5a"/"#0048ba","#b0bf1a"
  geom_boxplot(size=0.8, # 线条的粗细
               width=0.2, # 误差棒的宽度，和下边箱型图的宽度设置一样宽，不然会错开
               linetype="solid", # 线条的类型，solid是实线，dashed是虚线，dotted是点线，blank是没有线
               position=position_dodge(.3), #同一个X对应的不同颜色的组别别之间的间距
               color="black"#线条颜色同一设置成黑色，不根据分组来设置
  )+
  geom_signif(comparisons = list(c("Poor survival cells","Good survival cells")),#设置需要比较的组
         test = wilcox.test, ##计算方法
         textsize = 8,
         y_position = 0.65,#图中横线位置设置
         tip_length = c(0.02,0.02),#横线下方的竖线设置
         size=0.75,color="black",
         map_signif_level = T)+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13,face = "bold",angle = 45,hjust = 1),
        axis.text.y = element_text(size = 13,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.line.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 2),
        axis.line.y=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 2),
        legend.position = "none",
        plot.title=element_text(hjust=0.5,size=15,face = "bold"))+
  labs(title="High Immunoactivity Group")+
  guides(fill=guide_legend( keywidth=0.5, keyheight=0.10, default.unit="inch") )
dev.off()

TPM_HvsL<-rbind(TPMH,TPML)
TPM_HvsL$Immune=factor(TPM_HvsL$Immune, levels=c("High_Immune","Low_Immune"))
TPM_HvsL$corScore<-as.numeric(TPM_HvsL$corScore)
HvsL_pvalue<-wilcox.test(corScore ~ Immune,data = TPM_HvsL)
#################################################################################
#################################################################################
tabPos=data.frame()
for(i in colnames(TPM_L_pos)){
  rtpos=TPM_L_pos[,i]
  tabPos=rbind(tabPos,cbind(rtpos,Group="Poor survival cells",Immune="High_Immune"))
}
colnames(tabPos)[1]<-"corScore"

tabNeg=data.frame()
for(j in colnames(TPM_L_neg)){
  rtneg=TPM_L_neg[,j]
  tabNeg=rbind(tabNeg,cbind(rtneg,Group="Good survival cells",Immune="High_Immune"))
}
colnames(tabNeg)[1]<-"corScore"

TPML<-rbind(tabPos,tabNeg)
write.csv(TPML,"./DOC_PIC/RESULT/TPM_NEW/TPML.csv",col.names = T)

TPML$Group=factor(TPML$Group, levels=c("Poor survival cells","Good survival cells"))
TPML$corScore<-as.numeric(TPML$corScore)
L_pvalue<-wilcox.test(corScore ~ Group,data = TPML)

pdf(file="./DOC_PIC/RESULT/TPM_NEW/Lpos_vs_Lneg_boxplot.pdf",width=5,height=7)#输出图片文件
# par(mar=c(1,1,1,1))
ggplot(data=TPML, aes(x=Group, y=corScore, fill=Group)) +
  scale_fill_manual(values = c("indianred1","royalblue"))+##   L"#9966cc","#56a0d3"/H"#eead0e","#a2cd5a"/"#0048ba","#b0bf1a"
  geom_boxplot(size=0.8, # 线条的粗细
               width=0.2, # 误差棒的宽度，和下边箱型图的宽度设置一样宽，不然会错开
               linetype="solid", # 线条的类型，solid是实线，dashed是虚线，dotted是点线，blank是没有线
               position=position_dodge(.3), #同一个X对应的不同颜色的组别别之间的间距
               color="black"#线条颜色同一设置成黑色，不根据分组来设置
  )+
  geom_signif(comparisons = list(c("Poor survival cells","Good survival cells")),#设置需要比较的组
              test = wilcox.test, ##计算方法
              textsize = 8,
              y_position = 0.65,#图中横线位置设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.75,color="black",
              map_signif_level = T)+
  theme_classic()+ #设置图背景，我常用theme_classic()，或者theme_test，或者theme_minimal
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13,face = "bold",angle = 45,hjust = 1),
        axis.text.y = element_text(size = 13,face = "bold"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.line.x=element_line(color="black",size=1,lineend = 1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 2),
        axis.line.y=element_line(color="black",size=1,lineend = 1),
        axis.ticks.y=element_line(color="black",size=1,lineend = 2),
        legend.position = "none",
        plot.title=element_text(hjust=0.5,size=15,face = "bold"))+
  labs(title="Low Immunoactivity Group")
dev.off()

source("CIBERSORT.R")
mergeSTADCIBERSORT<-CIBERSORT("./test/LM22.txt", "./test/mergeSTAD.txt", perm = 1000, QN = T)
write.csv(mergeSTADCIBERSORT,"./test/mergeSTADCIBERSORT.csv",row.names = T,col.names = T,sep = "\t",quote = F)






