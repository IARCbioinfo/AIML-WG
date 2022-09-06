#load libraries
library(DT)
library(SummarizedExperiment)
library(tidyverse)

## load data
load("TCGA-MESOTranscriptome_ProfilingFri_Jul_29_11:39:27_2022.RData")
MESO.Rnaseq.SE = data
load("TCGA-BRCATranscriptome_ProfilingFri_Jul_29_12:11:29_2022.RData")
BRCA.Rnaseq.SE = data

## Get expression in TPM
MESO.Rnaseq.expr.tpm = assay(MESO.Rnaseq.SE,"tpm_unstrand")
BRCA.Rnaseq.expr.tpm = assay(BRCA.Rnaseq.SE,"tpm_unstrand")

# Analysis
## PCA of log(tpm+1) expression of 5000 most variable genes, color by subtype
library(ade4)
pca.MESO = dudi.pca(t(log(MESO.Rnaseq.expr.tpm[order(apply(MESO.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)))
pca.BRCA = dudi.pca(t(log(BRCA.Rnaseq.expr.tpm[order(apply(BRCA.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)))

## UMAP
library(umap)
umap.MESO = umap(t(log(MESO.Rnaseq.expr.tpm[order(apply(MESO.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)))
umap.BRCA = umap(t(log(BRCA.Rnaseq.expr.tpm[order(apply(BRCA.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)))

# plot
library(patchwork)
## create tibble objects for plotting
embeddings.MESO.tib = tibble(PCA=pca.MESO$li, UMAP = umap.MESO$layout, data=as.data.frame(colData(MESO.Rnaseq.SE) ) )
embeddings.BRCA.tib = tibble(PCA=pca.BRCA$li, UMAP = umap.BRCA$layout, data=as.data.frame(colData(BRCA.Rnaseq.SE) ) )

## plot first 2 axes
ggMESO_PCA  = ggplot(embeddings.MESO.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=data$primary_diagnosis)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
ggMESO_UMAP = ggplot(embeddings.MESO.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=data$primary_diagnosis)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")

ggBRCA_PCA  = ggplot(embeddings.BRCA.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=data$paper_BRCA_Subtype_PAM50)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
ggBRCA_UMAP = ggplot(embeddings.BRCA.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=data$paper_BRCA_Subtype_PAM50)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")

((ggMESO_PCA+guides(color="none"))+ ggMESO_UMAP)/((ggBRCA_PCA+guides(color="none"))+ ggBRCA_UMAP)


