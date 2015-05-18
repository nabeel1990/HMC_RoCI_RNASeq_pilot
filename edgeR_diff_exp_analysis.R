##  Load appropriate libraries
library(Biobase)
library(limma)
library(pamr)
library(gplots)
require(limma)
require(edgeR)
library(xlsx)



##	load the rna-seq data
ReadCounts=read.table("featureCounts_STAR.txt",header=T,row.names=1,sep="\t")
ReadCounts[1:3,]
dim(ReadCounts)
logReadCounts=log2(ReadCounts+.5)


#Non-specific filtering 
mads = apply(ReadCounts,1,mad)
qkeep = quantile(mads,.9,na.rm=TRUE)
ok = which(mads > qkeep)
ReadCounts2 = ReadCounts[ok,]
logReadCounts2=log2(ReadCounts2+.5)
plotMDS(logReadCounts2)


#MA plots based on counts
par(ask=F,mfrow=c(2,2))
tmp=maPlot(ReadCounts2[,1],ReadCounts2[,2])
lines(lowess(tmp[["A"]],tmp[["M"]]),col=2)
abline(h=0,col="grey")


lib.size=apply(ReadCounts2,2,sum)
pCounts=t(ReadCounts2)/lib.size
pCounts=t(pCounts)
tmp=maPlot(pCounts[,1],pCounts[,2])
lines(lowess(tmp[["A"]],tmp[["M"]]),col=2)
abline(h=0,col="grey")

par(mfrow=c(1,1),ask=F)
dist=as.dist(1-cor(logReadCounts2))
plot(hclust(dist))

plotMDS(logReadCounts2)


par(mfrow=c(1,1))
#Convert data into matrix

TotalReads<- as.matrix(ReadCounts2)
rownames(TotalReads)=rownames(ReadCounts2)
colnames(TotalReads)=colnames(ReadCounts2)
write.xlsx(TotalReads, file = "Edge_counts_filtered.xls")

#Specify the treatment levels
treatments<-c("sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","ARDS","ARDS","ARDS","ARDS","ARDS","ARDS","ARDS")
treatments=factor(treatments)

#Estimate the porportion of mapped reads occupied by the highest expressed genes
maxp=function(v) max(v)/sum(v) 
apply(TotalReads,2,maxp)
#The most highly expressed gene in each sample takes up a significant proprotion of reads.

#Estimate the cummulative proportion of mapped reads occupied by 10 highly expressing genes 
totRead=apply(TotalReads,1,sum)
cumsum(sort(totRead,decreasing=TRUE)/sum(totRead))[1:10]
bigGenes=which(totRead>= sort(totRead,decreasing=TRUE)[10])
lib.size=colSums(TotalReads)
apply(TotalReads[bigGenes,],2,sum)/lib.size
#this number is even more.In some samples more than 58% of reads are in top 10 expressed genes while for SI12 it is 68%
#Hence we need to normalize the values in order to detect significantly differentially expressed genes


#Calculate overall normalization factors
normFact=calcNormFactors(TotalReads,method="TMM")
normFact

#Now we start the edgeR analysis
#Create edgeR object 
d=DGEList(counts=TotalReads, group=treatments,genes=rownames(TotalReads))
d=calcNormFactors(d,method="TMM")
d=estimateCommonDisp(d)
d$common.dispersion

#We now estimate the tagwise dispersion(different dispersion for each individual gene) 
d=estimateTagwiseDisp(d)

#PART-I of edgeR analysis
#Pairwise comparison using Fisher's exact test
exact_results=exactTest(d,pair=c("sepsis","ARDS"))
head(exact_results$table)
topTags(exact_results)
exact.p=exact_results$table[,3]
hist(exact.p,main="Sepsis vs ARDS",xlab="p-values",nclass=50)

# Ajust using FDR
de1 <- decideTestsDGE(exact_results, p = 0.05, adjust = "fdr")
# Get the number of genes differentially expressed
summary(de1)

DE.unpaired <- as.data.frame(topTags(exact_results, n = 1)) # all DE genes
Results.unpaired <- as.data.frame(topTags(exact_results, n=2511)) # all results
write.xlsx(Results.unpaired, file = "EdgeR_Fisher's exact test.xls")

#MA plot
detags <- rownames(d)[as.logical(de1)]
plotSmear(exact_results, de.tags = detags, cex = 0.5, main = "MA Plot edgeR Exact Test")
abline(h = c(-1, 1), col = "blue")


png("VolcanoPlot_ARDS_v_sepsis_FDR_vs_FC_Fisher_Exact_test_EdgeR.png",    # create PNG for the heat map
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
# Volcano plot using FDR values
plot(-log10(Results.unpaired$FDR) ~ Results.unpaired$logFC, xlab = "log2 fold change",
     ylab = "-log10 FDR", main = "ARDS vs. Sepsis edgeR Exact test")
points(-log10(Results.unpaired$FDR[row.names(Results.unpaired) %in%
                                     row.names(DE.unpaired)]) ~ Results.unpaired$logFC[row.names(Results.unpaired)
                                                                                       %in% row.names(DE.unpaired)], col = "red", pch = 16)
dev.off()

png("VolcanoPlot_ARDS_v_sepsis_PValues_vs_FC_Fisher_exact_Test_EdgeR.png",    # create PNG for the heat map
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
#Volcano plot using p-values
plot(-log10(Results.unpaired$PValue) ~ Results.unpaired$logFC, xlab = "log2 fold change",
     ylab = "-log10 P-value", main = "ARDS vs. Sepsis edgeR Exact test")
points(-log10(Results.unpaired$PValue[row.names(Results.unpaired) %in%
                                     row.names(DE.unpaired)]) ~ Results.unpaired$logFC[row.names(Results.unpaired)

                                                                                       %in% row.names(DE.unpaired)], col = "red", pch = 16)
dev.off()
#Part-II of edgeR analysis
#GLM model
design=model.matrix(~0+d$samples$group)
colnames(design)=levels(d$samples$group)

D=d
D=estimateGLMCommonDisp(D,design)
D$common.dispersion
D=estimateGLMTagwiseDisp(D,design)
#comparison of Ordinary and GLM dispersion 
plot(d$tag,D$tag,xlab="ordinary dispersion",ylab="GLM dispersion")

#Fititng the data in the model design 
fitLung=glmFit(D,design)
colnames(design)
ARDSvssepsis=c(1,-1)
lrt.ARDSvssepsis=glmLRT(fitLung,contrast=ARDSvssepsis)
head(lrt.ARDSvssepsis$table)
Lungs.p=lrt.ARDSvssepsis$table[,4]
hist(Lungs.p,main="Sepsis vs ARDS",xlab="p-values",nclass=50)

topTags(lrt.ARDSvssepsis)
de2 <- decideTestsDGE(lrt.ARDSvssepsis, adjust.method = "fdr")
summary(de2)

DE.paired <- as.data.frame(topTags(lrt.ARDSvssepsis, n = 1)) # all DE Genes
Results.paired <- as.data.frame(topTags(lrt.ARDSvssepsis, n = 2511)) # all results
# Number of DE genes with p < 0.05 (uncorrected):
dim(Results.paired[Results.paired$PValue < 0.05, ])[1]

write.xlsx(Results.paired, file = "EdgeR_GLM_model.xls")

detags2 <- rownames(D)[as.logical(de2)]
plotSmear(lrt.ARDSvssepsis, de.tags = detags2)
abline(h = c(-1, 1), col = "blue")

# Select Genes with log2FC > 1 (2-fold change):
BigFC <- DE.paired[abs(DE.paired$logFC) > 1, ]
# Select Genes with log2CPM > 5 (32 CPM):
BiglogCPM <- BigFC[BigFC$logCPM > 5, ]

detags3 <- BiglogCPM$genes
plotSmear(lrt.ARDSvssepsis, de.tags = detags3, cex = 0.4, main = "MA Plot paired edgeR")
abline(h = c(-1, 1), col = "blue")
legend("topright", legend = "DE genes (FC > 2, CPM > 32)", pch = 16, col = "red")

png("VolcanoPlot_ARDS_v_sepsis_PValues_vs_FC_EdgeR.png",    # create PNG for the heat map
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot(-log10(Results.paired$PValue) ~ Results.paired$logFC, xlab = "log2 fold change",
     ylab = "-log10 P-value", main = "ARDS vs. Sepsis edgeR analysis")
points(-log10(Results.paired$PValue[Results.paired$genes %in% BiglogCPM$genes]) ~
         Results.paired$logFC[Results.paired$genes %in% BiglogCPM$genes],
       col = "red", pch = 16)
dev.off()


#Part-III of edgeR analysis
#Adjust for phenotypes
dx<-c("sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","sepsis","ARDS","ARDS","ARDS","ARDS","ARDS","ARDS","ARDS")
sepsis<-c("septic shock","septic shock","severe sepsis","septic shock","severe sepsis","severe sepsis","severe sepsis","septic shock","septic shock","septic shock","severe sepsis","septic shock","severe sepsis","septic shock","severe sepsis","septic shock","severe sepsis","septic shock","severe sepsis")
apache<-c(15,27,19,14,19,20,9,9,31,15,19,13,18,14,9,11,25,27,25)
phenos<-cbind(dx,sepsis)
phenos<-cbind(phenos,apache)
demos<-read.table("demos.txt",header=TRUE,sep="\t")
phenos<-cbind(demos,phenos)

#GLM model
design_adj=model.matrix(~0+d$samples$group+as.numeric(phenos$Age)+factor(phenos$Race,levels=c("White","Black","Hispanic"))+factor(phenos$Gender,levels=c("Male","Female")))
colnames(design_adj)<-c("ARDS","sepsis","age","black","hispanic","female")

D=d
D=estimateGLMCommonDisp(D,design_adj)
D$common.dispersion
D=estimateGLMTagwiseDisp(D,design_adj)
#comparison of Ordinary and GLM dispersion 
plot(d$tag,D$tag,xlab="ordinary dispersion",ylab="GLM dispersion")

#Fititng the data in the model design 
fitLung_adj=glmFit(D,design_adj)
colnames(design_adj)
ARDSvssepsis_adj=c(1,-1, 0, 0, 0, 0)
lrt.ARDSvssepsis_adj=glmLRT(fitLung_adj,contrast=ARDSvssepsis_adj)
head(lrt.ARDSvssepsis_adj$table)
Lungs_adj.p=lrt.ARDSvssepsis_adj$table[,4]
hist(Lungs_adj.p,main="Sepsis vs ARDS",xlab="p-values",nclass=50)

topTags(lrt.ARDSvssepsis_adj)
de2_adj <- decideTestsDGE(lrt.ARDSvssepsis_adj, adjust.method = "fdr")
summary(de2_adj)

DE.paired <- as.data.frame(topTags(lrt.ARDSvssepsis_adj, n = 3)) # all DE Genes
Results.paired <- as.data.frame(topTags(lrt.ARDSvssepsis_adj, n = 2511)) # all results
# Number of DE genes with p < 0.05 (uncorrected):
dim(Results.paired[Results.paired$PValue < 0.05, ])[1]

write.xlsx(Results.paired, file = "EdgeR_GLM_model_adjusted_race_gender_age.xls")

detags2 <- rownames(D)[as.logical(de2_adj)]
plotSmear(lrt.ARDSvssepsis_adj, de.tags = detags2)
abline(h = c(-1, 1), col = "blue")

# Select Genes with log2FC > 1 (2-fold change):
BigFC <- DE.paired[abs(DE.paired$logFC) > 1, ]
# Select Genes with log2CPM > 5 (32 CPM):
BiglogCPM <- BigFC[BigFC$logCPM > 5, ]

detags3 <- BiglogCPM$genes
plotSmear(lrt.ARDSvssepsis_adj, de.tags = detags3, cex = 0.4, main = "MA Plot paired edgeR")
abline(h = c(-1, 1), col = "blue")
legend("topright", legend = "DE genes (FC > 2, CPM > 32)", pch = 16, col = "red")


x <- topTags(lrt.ARDSvssepsis_adj, n = 3)
idx<-as.numeric(rownames(x))
pdf("boxplot_ARDS_v_sepsis_adj.pdf")
for (i in 1:3)
{
  boxplot(TotalReads[idx[i],1:12],TotalReads[idx[i],13:19],
          main=names[i], names=c("Sepsis","ARDS"),
          cex.main=0.7,sub=paste("FDR =",signif(DE.paired$FDR[i],3)))
}
dev.off()


png("VolcanoPlot_ARDS_v_sepsis_FDR_vs_FC_adjusted_EdgeR.png",    # create PNG for the heat map
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
# Volcano plot using FDR values
plot(-log10(Results.paired$FDR) ~ Results.paired$logFC, xlab = "log2 fold change",
     ylab = "-log10 FDR", main = "ARDS vs. Sepsis edgeR Adjusted for Race Age and Gender")
points(-log10(Results.paired$FDR[row.names(Results.paired) %in%
                                     row.names(DE.paired)]) ~ Results.paired$logFC[row.names(Results.paired)
                                                                                       %in% row.names(DE.paired)], col = "red", pch = 16)

dev.off()


png("VolcanoPlot_ARDS_v_sepsis_PValues_vs_FC_adjusted_EdgeR.png",    # create PNG for the heat map
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot(-log10(Results.paired$PValue) ~ Results.paired$logFC, xlab = "log2 fold change",
     ylab = "-log10 P-value", main = "ARDS vs. Sepsis edgeR analysis Adjusted for Race Age and Gender")
points(-log10(Results.paired$PValue[Results.paired$genes %in% BiglogCPM$genes]) ~
         Results.paired$logFC[Results.paired$genes %in% BiglogCPM$genes],
       col = "red", pch = 16)
dev.off()
x <- topTags(lrt.ARDSvssepsis_adj, n = 3)

idx<-as.numeric(rownames(x))
y <- TotalReads[idx,]
heatmap(y, Rowv=NA, Colv= NA)
