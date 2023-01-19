setwd('/Users/maria/Desktop/Kenkel_lab/Anthopleura/Nov2020_AeleHeatStress/2bRAD/bams/')

library(WGCNA)
library(sparcl)

# reading list of bam files = order of samples in IBS matrix
bams=read.table("bam_list",header=F)[,1]
#bams=read.table("bam_clean",header=F)[,1] # without 13-O-4 bc low coverage
bams=sub("\\_host.bam","",bams,perl=T)

# reading IBS matrix based on SNPs with allele frequency >= 0.05:
ma = as.matrix(read.table("AeleResult.ibsMat"))
#ma = as.matrix(read.table("AeleResult2.ibsMat")) # without 13-O-4 bc low coverage
#ma = as.matrix(read.table("AeleResult5.ibsMat")) # without restricting non-HWE sites bc clonality
dimnames(ma)=list(bams,bams)

# plotting hierarchical clustering tree
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.6)

#known technical replicates
clones=bams[grep('*_TR',bams)]

#first look at tech reps
col_clone=rep('black',nrow(ma))
col_clone[bams %in% clones] = 'red'
par(cex=1)
ColorDendrogram(hclust(as.dist(ma),"ave"), y = col_clone, branchlength=0.05, labels = bams,ylab='(1-IBS)')
abline(h=0.22,col="black",lty=3)

# looks like 7 genotypes, not 12 like planned
# 10 and 11 same genet
# 1, 2, 3, 7, & 8 same genet

# PROBLEM -- one rep of 11 clustered with 13, ugggghhh I feel like this means we need to genotype all samples :(
# Carly thought this issue might change with settings but it looks like it stays the same

# looking at results for 11 and 13 run 3
setwd('/Users/maria/Desktop/Kenkel_lab/Anthopleura/Nov2020_AeleHeatStress/2bRAD/run3/')

bams=read.table("bam_list",header=F)[,1]
names=read.csv("lig_to_sample.csv")
names=names[order(names$bam,bams),]

samples=names$sample
# color tech reps
#known technical replicates
clones=c('11-I-6','13-I-6')
#first look at tech reps
col_clone=rep('black',nrow(ma))
col_clone[samples %in% clones] = 'red'

ma = as.matrix(read.table("AeleResult.ibsMat")) # without restricting non-HWE sites bc clonality
dimnames(ma)=list(samples,samples)

# plotting hierarchical clustering tree
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.6)

par(cex=1)
ColorDendrogram(hclust(as.dist(ma),"ave"), y = col_clone, branchlength=0.05, labels = samples,ylab='(1-IBS)')
abline(h=0.32,col="black",lty=3)

