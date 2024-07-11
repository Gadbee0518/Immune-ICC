setwd("E:/GSEA_IMMU/CHOL/CODE/FOUR3")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("survivalROC")
library(survivalROC)
library(survival)
library(survminer)
library(pROC)
library(ggplot2)


diff_TPM16<-read.table("./DOC_PIC/DIFF_TPM16.txt",header = T,row.names = 1,sep = "\t",check.names = F)
clinical16<-read.table("./DOC_PIC/clinical_LH16.txt",header = T,row.names = 1,sep = "\t",check.names = F)
diff_TPM16<-diff_TPM16[,rownames(clinical16)]
clinical16<-clinical16[,c(3,6)]
clinical16$TIME<-clinical16$TIME/30


########################################################################
#black
##########################################################################
blackgenes<-read.table("./DOC_PIC/blackwgcnaGenes.txt")
blackBulk<-diff_TPM16[blackgenes$V1,]
blackBulk<-data.frame(t(blackBulk))
black<-cbind(clinical16,blackBulk)
# remove(blackBulk,blackgenes)
blackSingleCox=data.frame()
for(i in colnames(black[,3:ncol(black)])){
  cox <- coxph(Surv(black$TIME, black$STATUS) ~ black[,i], data = black)
  coxSummary = summary(cox)
  blackSingleCox=rbind(blackSingleCox,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                          z=coxSummary$coefficients[,"z"],
                                          pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}
write.csv(blackSingleCox,file="./DOC_PIC/blackSingleCox.csv",sep="\t",row.names=F,quote=F)
blackSecletGenes<-read.table("./DOC_PIC/blackSecletGenes.txt")
black<-cbind(black[,c(1,2)],black[,blackSecletGenes$V1])

blackMcox <- coxph(Surv(black$TIME, black$STATUS) ~ ., data = black)
blackMcox <- step(blackMcox,direction = "both")
blackRiskScore <- predict(blackMcox,type="risk",newdata=black)
blackRisk=as.vector(ifelse(blackRiskScore>median(blackRiskScore),"high","low"))
black_risk<-cbind(black[,1:2],blackRiskScore,blackRisk)
write.table(black_risk,file="./DOC_PIC/black_risk.txt",sep="\t",quote=F,row.names=T,col.names = T)

black_diff=survdiff(Surv(black_risk$TIME, black_risk$STATUS) ~ blackRisk,data = black_risk)
black_pValue=1-pchisq(black_diff$chisq,df=1)
black_pValue=round(black_pValue,5)
black_fit <- survfit(Surv(black_risk$TIME, black_risk$STATUS) ~ blackRisk, data = black_risk)
summary(black_fit)    #查看五年生存率
pdf(file="./DOC_PIC/black_survival.pdf",width = 6,height = 6)
par(mar=c(5,5,5,1))
plot(black_fit, lty = 2:2,col=c("red","black"),
     xlab="time (month)",
     ylab="black module genes survival rate",
     main=paste("survival curve (p=", black_pValue ,")",sep=""),
     mark.time=T,
     lwd = 2,
     cex.main=2, cex.lab=1.5, cex.axis=1.5
     )
legend("topright", c("high risk", "low risk"), 
       lty = 2:2, 
       col=c("red","black"),
       text.font = 1,
       cex = 1.2)
dev.off()

blackExpected <- predict(blackMcox,type="expected",newdata=black)
blackroc<-cbind(black_risk,blackExpected)
pdf(file="./DOC_PIC/black_ROC.pdf",width=6,height = 6)
par(mar=c(5,5,5,1))
roc=survivalROC(Stime=blackroc$TIME, 
                status=blackroc$STATUS, marker = blackroc$blackExpected, 
                predict.time =60, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='black', 
     xlab="False positive rate", ylab="True positive rate",
     main="black module ROC curve",
     lwd = 4, cex.main=1.5, cex.lab=2, cex.axis=1.5)+
  grid(col = "grey", lty = 2,lwd = 2)+
  abline(0,1,col="#006b3c",lty=5,lwd=3)+
  legend(x=0.5,y=0.1, paste0("AUC = ",round(roc$AUC,5)),
       box.lty = 0,
       text.font = 3,text.col = "black",
       adj = c(0.25,-0.5),cex = 2,bg="transparent")
dev.off()

########################################################################
#blue
##########################################################################
bluegenes<-read.table("./DOC_PIC/bluewgcnaGenes.txt")
blueBulk<-diff_TPM16[bluegenes$V1,]
blueBulk<-data.frame(t(blueBulk))
blue<-cbind(clinical16,blueBulk)
# remove(blueBulk,bluegenes)

blueSingleCox=data.frame()
for(i in colnames(blue[,3:ncol(blue)])){
  cox <- coxph(Surv(blue$TIME, blue$STATUS) ~ blue[,i], data = blue)
  coxSummary = summary(cox)
  blueSingleCox=rbind(blueSingleCox,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                                    z=coxSummary$coefficients[,"z"],
                                                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}
write.csv(blueSingleCox,file="./DOC_PIC/blueSingleCox.csv",sep="\t",row.names=F,quote=F)
blueSecletGenes<-read.table("./DOC_PIC/blueSecletGenes.txt")
blue<-cbind(blue[,c(1,2)],blue[,blueSecletGenes$V1])

blueMcox <- coxph(Surv(blue$TIME, blue$STATUS) ~ ., data = blue)
blueMcox <- step(blueMcox,direction = "both")
blueRiskScore <- predict(blueMcox,type="risk",newdata=blue)
blueRisk=as.vector(ifelse(blueRiskScore>median(blueRiskScore),"high","low"))
blue_risk<-cbind(blue[,1:2],blueRiskScore,blueRisk)
write.table(blue_risk,file="./DOC_PIC/blue_risk.txt",sep="\t",quote=F,row.names=T,col.names = T)

blue_diff=survdiff(Surv(blue_risk$TIME, blue_risk$STATUS) ~ blueRisk,data = blue_risk)
blue_pValue=1-pchisq(blue_diff$chisq,df=1)
blue_pValue=round(blue_pValue,5)
blue_fit <- survfit(Surv(blue_risk$TIME, blue_risk$STATUS) ~ blueRisk, data = blue_risk)
summary(blue_fit)    #查看五年生存率
pdf(file="./DOC_PIC/blue_survival.pdf",width = 6,height = 6)
par(mar=c(5,5,5,1))
plot(blue_fit, lty = 2:2,col=c("red","blue"),xlab="time (month)",ylab="blue module genes survival rate",
     main=paste("survival curve (p=", blue_pValue ,")",sep=""),mark.time=T,
     lwd = 2,
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
legend("topright", 
       c("high risk", "low risk"), 
       lty = 2:2, 
       col=c("red","blue"),
       text.font = 1,
       cex = 1.2)
dev.off()

blueExpected <- predict(blueMcox,type="expected",newdata=blue)
blueroc<-cbind(blue_risk,blueExpected)
pdf(file="./DOC_PIC/blue_ROC.pdf",width=6,height = 6)
par(mar=c(5,5,5,1))
roc=survivalROC(Stime=blueroc$TIME, 
                status=blueroc$STATUS, marker = blueroc$blueExpected, 
                predict.time =60, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='blue', 
     xlab="False positive rate", ylab="True positive rate",
     main="blue module ROC curve",
     lwd = 4, cex.main=1.5, cex.lab=2, cex.axis=1.5)+
  grid(col = "grey", lty = 2,lwd = 2)+
  abline(0,1,col="#006b3c",lty=5,lwd=3)+
  legend(x=0.5,y=0.1, paste0("AUC = ",round(roc$AUC,5)),
       box.lty = 0,
       text.font = 3,text.col = "blue",
       adj = c(0.25,-0.5),cex = 2,bg="transparent")
dev.off()

########################################################################
#turquoise
##########################################################################
turquoisegenes<-read.table("./DOC_PIC/turquoisewgcnaGenes.txt")
turquoiseBulk<-diff_TPM16[turquoisegenes$V1,]
turquoiseBulk<-data.frame(t(turquoiseBulk))
turquoise<-cbind(clinical16,turquoiseBulk)
# remove(turquoiseBulk,turquoisegenes)

turquoiseSingleCox=data.frame()
for(i in colnames(turquoise[,3:ncol(turquoise)])){
  cox <- coxph(Surv(turquoise$TIME, turquoise$STATUS) ~ turquoise[,i], data = turquoise)
  coxSummary = summary(cox)
  turquoiseSingleCox=rbind(turquoiseSingleCox,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                            z=coxSummary$coefficients[,"z"],
                            pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}
write.csv(turquoiseSingleCox,file="./DOC_PIC/turquoiseSingleCox.csv",sep="\t",row.names=F,quote=F)

turquoiseSecletGenes<-read.table("./DOC_PIC/turquoiseSecletGenes.txt")
turquoise<-cbind(turquoise[,c(1,2)],turquoise[,turquoiseSecletGenes$V1])
turquoiseMcox <- coxph(Surv(turquoise$TIME, turquoise$STATUS) ~ ., data = turquoise)
turquoiseMcox <- step(turquoiseMcox,direction = "both")
turquoiseRiskScore <- predict(turquoiseMcox,type="risk",newdata=turquoise)
turquoiseRisk=as.vector(ifelse(turquoiseRiskScore>median(turquoiseRiskScore),"high","low"))
turquoise_risk<-cbind(turquoise[,1:2],turquoiseRiskScore,turquoiseRisk)
write.table(turquoise_risk,file="./DOC_PIC/turquoise_risk.txt",sep="\t",quote=F,row.names=T,col.names = T)

turquoise_diff=survdiff(Surv(turquoise_risk$TIME, turquoise_risk$STATUS) ~ turquoiseRisk,data = turquoise_risk)
turquoise_pValue=1-pchisq(turquoise_diff$chisq,df=1)
turquoise_pValue=round(turquoise_pValue,5)
turquoise_fit <- survfit(Surv(turquoise_risk$TIME, turquoise_risk$STATUS) ~ turquoiseRisk, data = turquoise_risk)
summary(turquoise_fit)    #查看五年生存率
pdf(file="./DOC_PIC/turquoise_survival.pdf",width = 6,height = 6)
par(mar=c(5,5,5,1))
plot(turquoise_fit, lty = 2:2,col=c("red","turquoise"),xlab="time (month)",ylab="turquoise module genes survival rate",
     main=paste("survival curve (p=", turquoise_pValue ,")",sep=""),mark.time=T,
     lwd = 2,
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
legend("topright", 
       c("high risk", "low risk"), 
       lty = 2:2, 
       col=c("red","turquoise"),
       text.font = 1,
       cex = 1.2)
dev.off()

turquoiseExpected <- predict(turquoiseMcox,type="expected",newdata=turquoise)
turquoiseroc<-cbind(turquoise_risk,turquoiseExpected)
pdf(file="./DOC_PIC/turquoise_ROC.pdf",width=6,height = 6)
par(mar=c(5,5,5,1))
roc=survivalROC(Stime=turquoiseroc$TIME, 
                status=turquoiseroc$STATUS, marker = turquoiseroc$turquoiseExpected, 
                predict.time =60, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='turquoise', 
     xlab="False positive rate", ylab="True positive rate",
     main="turquoise module ROC curve",
     lwd = 4, cex.main=1.5, cex.lab=2, cex.axis=1.5)+
  grid(col = "grey", lty = 2,lwd = 2)+
  abline(0,1,col="#006b3c",lty=5,lwd=3)+
  legend(x=0.5,y=0.1, paste0("AUC = ",round(roc$AUC,5)),
       box.lty = 0,
       text.font = 3,text.col = "turquoise",
       adj = c(0.25,-0.5),cex = 2,bg="transparent")
dev.off()

########################################################################
#yellow
##########################################################################
yellowgenes<-read.table("./DOC_PIC/yellowwgcnaGenes.txt")
yellowBulk<-diff_TPM16[yellowgenes$V1,]
yellowBulk<-data.frame(t(yellowBulk))
yellow<-cbind(clinical16,yellowBulk)
# remove(yellowBulk,yellowgenes)

yellowSingleCox=data.frame()
for(i in colnames(yellow[,3:ncol(yellow)])){
  cox <- coxph(Surv(yellow$TIME, yellow$STATUS) ~ yellow[,i], data = yellow)
  coxSummary = summary(cox)
  yellowSingleCox=rbind(yellowSingleCox,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                                    z=coxSummary$coefficients[,"z"],
                                                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}
write.csv(yellowSingleCox,file="./DOC_PIC/yellowSingleCox.csv",sep="\t",row.names=F,quote=F)

yellowSecletGenes<-read.table("./DOC_PIC/yellowSecletGenes.txt")
yellow<-cbind(yellow[,c(1,2)],yellow[,yellowSecletGenes$V1])

yellowMcox <- coxph(Surv(yellow$TIME, yellow$STATUS) ~ ., data = yellow)
yellowMcox <- step(yellowMcox,direction = "both")
yellowRiskScore <- predict(yellowMcox,type="risk",newdata=yellow)
yellowRisk=as.vector(ifelse(yellowRiskScore>median(yellowRiskScore),"high","low"))
yellow_risk<-cbind(yellow[,1:2],yellowRiskScore,yellowRisk)
write.table(yellow_risk,file="./DOC_PIC/yellow_risk.txt",sep="\t",quote=F,row.names=T,col.names = T)

yellow_diff=survdiff(Surv(yellow_risk$TIME, yellow_risk$STATUS) ~ yellowRisk,data = yellow_risk)
yellow_pValue=1-pchisq(yellow_diff$chisq,df=1)
yellow_pValue=round(yellow_pValue,5)
yellow_fit <- survfit(Surv(yellow_risk$TIME, yellow_risk$STATUS) ~ yellowRisk, data = yellow_risk)
summary(yellow_fit)    #查看五年生存率
pdf(file="./DOC_PIC/yellow_survival.pdf",width = 6,height = 6)
par(mar=c(5,5,5,1))
plot(yellow_fit, lty = 2:2,col=c("red","#e4d00a"),xlab="time (month)",ylab="yellow module genes survival rate",
     main=paste("survival curve (p=", yellow_pValue ,")",sep=""),mark.time=T,
     lwd = 2,
     cex.main=2, cex.lab=1.5, cex.axis=1.5)
legend("topright", 
       c("high risk", "low risk"), 
       lty = 2:2, 
       col=c("red","#e4d00a"),
       text.font = 1,
       cex = 1.2)
dev.off()

yellowExpected <- predict(yellowMcox,type="expected",newdata=yellow)
yellowroc<-cbind(yellow_risk,yellowExpected)
pdf(file="./DOC_PIC/yellow_ROC.pdf",width=6,height = 6)
par(mar=c(5,5,5,1))
roc=survivalROC(Stime=yellowroc$TIME, status=yellowroc$STATUS, marker = yellowroc$yellowExpected,
                predict.time =60, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='#e4d00a',
     xlab="False positive rate", ylab="True positive rate",
     main="yellow module ROC curve",
     lwd = 4, cex.main=1.5, cex.lab=2, cex.axis=1.5)+
  grid(col = "grey", lty = 2,lwd = 2)+
  abline(0,1,col="#006b3c",lty=5,lwd=3)+
  legend(x=0.5,y=0.1, paste0("AUC = ",round(roc$AUC,5)),
       box.lty = 0,
       text.font = 3,text.col = "#e4d00a",
       adj = c(0.25,-0.5),cex = 2,bg="transparent")
dev.off()


bulkgenes1310<-as.data.frame(c(blackgenes$V1,bluegenes$V1,turquoisegenes$V1,yellowgenes$V1))
colnames(bulkgenes1310)<-"BulkGene"
write.table(bulkgenes1310,file = "./DOC_PIC/1310bulkgenes.txt",col.names = F,row.names = F)
