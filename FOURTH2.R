setwd("E:\\GSEA_IMMU\\CHOL\\CODE\\FOUR2") #设置工作目录
##整理metadata文件
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggstatsplot")
library("tidyverse")
library("ggplot2")
library("grid")
library("ggrepel")
library("ggpubr")
library("WGCNA")
library("stringr")

######WGCNA
######读取数据
chol16<-read.table("./DOC_PIC/WGCNA_DAT1_16.txt",row.names = 1,header = T,sep = "\t",check.names = F)
HLgroup<-read.table("./DOC_PIC/HLgroup.txt",row.names = 1,header = T,sep = "\t",check.names = F)
clinical16<-read.table("./DOC_PIC/clinical_LH16.txt",row.names = 1,header = T,sep = "\t",check.names = F)
UP<-read.csv("./DOC_PIC/UP.csv",row.names = 1,header = T)
DOWN<-read.csv("./DOC_PIC/DOWN.csv",row.names = 1,header = T)
diffnames<-c(rownames(UP),rownames(DOWN))
chol16<-chol16[rownames(chol16)%in%diffnames,]
######样品聚类，删除离散的样品,结果发现没有离散样本，那就验证一下就好了，不作为论文中出现的图
datExpr0<-as.data.frame(t(chol16))
sampleTree<-hclust(dist(datExpr0),method = "ward.D")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#####聚类树热图，也是没什么明显的效果，不用这个结果图了.删除了“TCGA-W5-AA2I-01A”离群样本
clinicalColors = numbers2colors(clinical16, signed = FALSE)
pdf(file="2_sample_heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree, clinicalColors,
                    groupLabels = names(clinical16),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

######软阈值
powers = c(1:20)       #幂指数范围1:20
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 10)
pdf(file="./DOC_PIC/scale_independence.pdf",width=10,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#获得 TOM 矩阵
adjacency <- adjacency(datExpr0, power = 9)
tomSim <- TOMsimilarity(adjacency)
rownames(tomSim) <- rownames(adjacency)
colnames(tomSim) <- colnames(adjacency)
#输出 TOM 矩阵
write.table(tomSim, './DOC_PIC/TOMsimilarity.txt', sep = '\t', col.names = T,row.names = T, quote = FALSE)

#层次聚类树，使用中值的非权重成对组法的平均聚合聚类
#TOM 相异度 = 1 – TOM 相似度
tomDis  <- 1-tomSim
geneTree <- hclust(as.dist(tomDis), method = 'average')
plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
     labels = FALSE, hang = 0.04)

#使用动态剪切树挖掘模块
minModuleSize <- 30  #模块基因数目
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tomDis,
                             deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)

table(dynamicMods)

#模块颜色指代
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="./DOC_PIC/Gene_dendrogram_and_Module_Colors.pdf",width=12,height=9)
plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
                    dendroLabels = FALSE, addGuide = TRUE, hang = 0.001, guideHang = 1,
                    main = 'Gene dendrogram and module colors')
dev.off()

###############################需要在服務器
#基因表达聚类树和共表达拓扑热图
plotSim <- -(1-tomSim)
#plot_sim <- log(tom_sim)
diag(plotSim) <- NA
pdf(file="./DOC_PIC/NetworkHeatmapPlot.pdf",width=6,height=6)
TOMplot(plotSim, geneTree, dynamicColors,
        main = 'Network heatmap plot, selected genes')
dev.off()


#计算基因表达矩阵中模块的特征基因（第一主成分）
MEList <- moduleEigengenes(datExpr0, colors = dynamicColors)
MEs <- MEList$eigengenes
head(MEs)[1:6]
#输出模块特征基因矩阵
write.table(MEs, './DOC_PIC/moduleEigengenes.txt', sep = '\t', row.names = T,col.names = T, quote = FALSE)

#通过模块特征基因计算模块间相关性，表征模块间相似度
ME_cor <- cor(MEs)
ME_cor[1:6,1:6]

#绘制聚类树观察
METree <- hclust(as.dist(1-ME_cor), method = 'average')
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')

# #探索性分析，观察模块间的相似性
# #height 值可代表模块间的相异度，并确定一个合适的阈值作为剪切高度
# #以便为低相异度（高相似度）的模块合并提供依据
# abline(h = 0.1, col = 'blue')
# abline(h = 0.05, col = 'red')
# 
# #模块特征基因聚类树热图
# plotEigengeneNetworks(MEs, '', cex.lab = 0.8, xLabelsAngle= 90,
#                       marDendro = c(0, 3, 1, 3), marHeatmap = c(3, 4, 1, 2))
# 
# #相似模块合并，以 0.25 作为合并阈值（剪切高度），在此高度下的模块将合并
# #近似理解为相关程度高于 0.75 的模块将合并到一起
# merge_module <- mergeCloseModules(datExpr0, dynamicColors, cutHeight = 0.15, verbose = 3)
# mergedColors <- merge_module$colors
# table(mergedColors)
# 
# #基因表达和模块聚类树
# pdf(file="./DOC_PIC/Gene_dendrogram_and_Merged_dynamic.pdf",width=12,height=9)
# plotDendroAndColors(geneTree, mergedColors, 'Merged dynamic',
#                     dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05)
# dev.off()
# 
# #使用上一步新组合的共表达模块的结果
# newMEs <- merge_module$newMEs

#患者基因共表达模块和临床表型的相关性分析
moduleClinicalCor <- cor(MEs, clinical16, use = 'p')
moduleClinicalCor[1:5,1:5]  #相关矩阵

#相关系数的 p 值矩阵
moduleClinicalPvalue <- corPvalueStudent(moduleClinicalCor, nrow(MEs))

#输出相关系数矩阵或 p 值矩阵
write.table(moduleClinicalCor, 'moduleTraitCor.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(moduleClinicalPvalue, 'moduleTraitPvalue.txt', sep = '\t', col.names = NA, quote = FALSE)

#相关图绘制
textMatrix <- paste(signif(moduleClinicalCor, 2), '\n(', signif(moduleClinicalPvalue, 1), ')', sep = '')
dim(textMatrix) <- dim(moduleClinicalCor)

pdf(file="./DOC_PIC/Module&trait_relationships.pdf",width=8,height=7)
par(mar=c(5,7,5,1))
labeledHeatmap(Matrix = moduleClinicalCor,
               xLabels = names(clinical16),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = T,
               # colors = blueWhiteRed(100),
               colors = colorRampPalette(c("blue","white","red"))(100),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = "Module-trait relationships")
dev.off()


######保存各个模块的数据
# moduleLables=net$colors
# moduleColors=labels2colors(net$colors)
# MEs=net$MEs
# geneTree<-net$dendrograms[[1]]
# geneTree1<-net$dendrograms
# save(MEs,moduleLables,moduleColors,geneTree,file = "./DOC_PIC/networkConstruction.RData")


############感兴趣性状的模块的具体基因分析
nGenes<-ncol(datExpr0)
nSamples<-nrow(datExpr0)
IG<-as.data.frame(clinical16$ImmuneGroup)
names(IG)<-"ImmuneGroup"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr0是每个基因在每个样本的表达量
gMMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("gMM", modNames, sep="")
names(gMMPvalue) = paste("p.gMM", modNames, sep="")

geneClinicalSignificance = as.data.frame(cor(datExpr0, IG, use = "p"))
gCSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneClinicalSignificance), nSamples))
names(geneClinicalSignificance) = paste("gCS.", names(IG), sep="")
names(gCSPvalue) = paste("p.CS.", names(IG), sep="")

module = "blue"
column = match(module, modNames)
moduleGenes = dynamicColors==module
pdf(file=paste("./DOC_PIC/",module,"Module_membership_vs_gene_significance.pdf",sep = ""),
    width=6,height=6)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneClinicalSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,pch = 10,
                   bg=0,abline = T,abline.color = "#458b00",abline.lty = 1,lwd=1,
                   lmFnc = lm)
dev.off()

# names(datExpr0)[dynamicColors==module]
geneInfo<-data.frame(geneName=names(datExpr0),
                     moduleColor=dynamicColors,
                     geneClinicalSignificance,
                     gCSPvalue)
modOrder<-order(-abs(cor(MEs,IG,use = "p")))
for (mod in 1:ncol(geneModuleMembership)) {
  colNames<-names(geneInfo)
  geneInfo<-data.frame(geneInfo,geneModuleMembership[,modOrder[mod]],gMMPvalue[,modOrder[mod]])
  names(geneInfo)<-c(colNames,paste("gMM.",modNames[modOrder[mod]],sep = ""),
                     paste("p.gMM",modNames[modOrder[mod]],sep = ""))
  print(mod)
}
geneOrder<-order(geneInfo$moduleColor,-abs(geneInfo$gCS.ImmuneGroup))
geneInfo<-geneInfo[geneOrder,]
write.csv(geneInfo,"./DOC_PIC/geneInfo.csv")

########提取指定模块的基因名主要是关心具体某个模块内部的基因
module = "turquoise"
# Select module probes
probes = colnames(datExpr0) ## 我们例子里面的probe就是基因
modProbes = probes[dynamicColors==module];
table(modProbes)
write(modProbes,paste("./DOC_PIC/",module,"wgcnaGenes.txt",sep = ""))
