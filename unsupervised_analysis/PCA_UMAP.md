---
title: "PCA_UMAP"
author: "N. Alcala"
date: "9/6/2022"
output: 
  html_document: 
    keep_md: yes
---


# R Code to compare embeddings on TCGA expression data

## load libraries 

```r
library(tidyverse)
library(SummarizedExperiment)
library(patchwork)
library(umap)
library(ade4)
library(keras)
library(cluster)
library(readxl)
```

## Load expression data from TCGA

```r
load("AIML-WG/data/TCGA-MESOTranscriptome_ProfilingFri_Jul_29_11:39:27_2022.RData")
MESO.Rnaseq.SE = data
# clean histology
MESO.Rnaseq.SE$primary_diagnosis[MESO.Rnaseq.SE$primary_diagnosis=="Epithelioid mesothelioma, malignant"] = "MME"
MESO.Rnaseq.SE$primary_diagnosis[MESO.Rnaseq.SE$primary_diagnosis=="Mesothelioma, biphasic, malignant"] = "MMB"
MESO.Rnaseq.SE$primary_diagnosis[MESO.Rnaseq.SE$primary_diagnosis%in%c("Mesothelioma, malignant","Fibrous mesothelioma, malignant")] = NA
table(MESO.Rnaseq.SE$primary_diagnosis)
```

```
## 
## MMB MME 
##  22  58
```

```r
load("AIML-WG/data/TCGA-BRCATranscriptome_ProfilingFri_Jul_29_12:11:29_2022.RData")
BRCA.Rnaseq.SE = data
```
### load additional clinical data from TCGA papers

```r
MESO_S1 = read_xlsx("AIML-WG/data/TCGA_TableS1.xlsx",sheet=2)
# match barcodes 
MESO_S1$barcode = MESO.Rnaseq.SE$barcode[sapply(MESO_S1$TCGA_barcode, function(x) grep(x = MESO.Rnaseq.SE$barcode,pattern = x) )]
MESOorder = sapply(MESO.Rnaseq.SE$barcode,function(x){res=which(MESO_S1$barcode==x);if(length(res)==0){res=NA};return(res)})
# remove samples not reported in study
MESO.Rnaseq.SE = MESO.Rnaseq.SE[,!is.na(MESOorder)]
MESO_S1 = MESO_S1[MESOorder[!is.na(MESOorder)],]
```

### make survival groups

```r
## MESO
MESO.Rnaseq.SE$survival_group = cut(MESO.Rnaseq.SE$days_to_death,breaks = quantile(MESO.Rnaseq.SE$days_to_death,c(0,0.33,0.67,1),na.rm=T))
MESO.Rnaseq.SE$survival_group[which(as.numeric(MESO.Rnaseq.SE$survival_group)<3 & MESO.Rnaseq.SE$vital_status=="Alive")] = NA

## BRCA
BRCA.Rnaseq.SE$paper_days_to_death = as.numeric(BRCA.Rnaseq.SE$paper_days_to_death)
```

```
## Warning: NAs introduced by coercion
```

```r
BRCA.Rnaseq.SE$survival_group = cut(BRCA.Rnaseq.SE$paper_days_to_death,breaks = quantile(BRCA.Rnaseq.SE$paper_days_to_death,c(0,0.33,0.67,1),na.rm=T))
BRCA.Rnaseq.SE$survival_group[which(as.numeric(BRCA.Rnaseq.SE$survival_group)<3 & BRCA.Rnaseq.SE$paper_vital_status=="Alive")] = NA
```


### Get expression in TPM

```r
MESO.Rnaseq.expr.tpm = assay(MESO.Rnaseq.SE,"tpm_unstrand")
BRCA.Rnaseq.expr.tpm = assay(BRCA.Rnaseq.SE,"tpm_unstrand")

MESO.Rnaseq.expr.tpm.logmat = t(log(MESO.Rnaseq.expr.tpm[order(apply(MESO.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1))
BRCA.Rnaseq.expr.tpm.logmat = t(log(BRCA.Rnaseq.expr.tpm[order(apply(BRCA.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1))
```

## Analysis
### PCA of log(tpm+1) expression of 5000 most variable genes

```r
pca.MESO = dudi.pca(MESO.Rnaseq.expr.tpm.logmat,scannf = F,nf=10)
pca.BRCA = dudi.pca(BRCA.Rnaseq.expr.tpm.logmat,scannf = F,nf=10)
```

### UMAP of log(tpm+1) expression of 5000 most variable genes

```r
umap.MESO = umap(MESO.Rnaseq.expr.tpm.logmat)
umap.BRCA = umap(BRCA.Rnaseq.expr.tpm.logmat)
```

### Autoencoder of of log(tpm+1) expression of 5000 most variable genes

```r
# set models
modelMESO <- keras_model_sequential()
modelBRCA <- keras_model_sequential()
modelMESO %>%
  layer_dense(units = 50, activation = "linear", input_shape = 5000) %>%
  layer_dense(units = 2, activation = "linear", name = "bottleneck") %>%
  layer_dense(units = 50, activation = "linear") %>%
  layer_dense(units = 5000)
modelBRCA %>%
  layer_dense(units = 50, activation = "linear", input_shape = 5000) %>%
  layer_dense(units = 2, activation = "linear", name = "bottleneck") %>%
  layer_dense(units = 50, activation = "linear") %>%
  layer_dense(units = 5000)
# view model layers
summary(modelMESO)
```

```
## Model: "sequential"
## ________________________________________________________________________________
## Layer (type)                        Output Shape                    Param #     
## ================================================================================
## dense_2 (Dense)                     (None, 50)                      250050      
## ________________________________________________________________________________
## bottleneck (Dense)                  (None, 2)                       102         
## ________________________________________________________________________________
## dense_1 (Dense)                     (None, 50)                      150         
## ________________________________________________________________________________
## dense (Dense)                       (None, 5000)                    255000      
## ================================================================================
## Total params: 505,302
## Trainable params: 505,302
## Non-trainable params: 0
## ________________________________________________________________________________
```

```r
# compile model
modelMESO %>% compile(loss = "mean_squared_error", optimizer = "rmsprop")
modelBRCA %>% compile(loss = "mean_squared_error", optimizer = "rmsprop")

# fit model
modelMESO %>% fit(x = MESO.Rnaseq.expr.tpm.logmat, y = MESO.Rnaseq.expr.tpm.logmat, epochs = 200,verbose = 1)
modelBRCA %>% fit(x = BRCA.Rnaseq.expr.tpm.logmat, y = BRCA.Rnaseq.expr.tpm.logmat, epochs = 200,verbose = 1)

# evaluate the performance of the model
mse.aeMESO <- evaluate(modelMESO, MESO.Rnaseq.expr.tpm.logmat, MESO.Rnaseq.expr.tpm.logmat)
mse.aeBRCA <- evaluate(modelBRCA, BRCA.Rnaseq.expr.tpm.logmat, BRCA.Rnaseq.expr.tpm.logmat)
mse.aeMESO
```

```
##    loss 
## 1.61935
```

```r
mse.aeBRCA
```

```
##      loss 
## 0.5589278
```

```r
intermediate_layer_modelMESO <- keras_model(inputs = modelMESO$input, outputs = get_layer(modelMESO, "bottleneck")$output)
intermediate_layer_modelBRCA <- keras_model(inputs = modelBRCA$input, outputs = get_layer(modelBRCA, "bottleneck")$output)
AE.MESO <- predict(intermediate_layer_modelMESO, MESO.Rnaseq.expr.tpm.logmat)
AE.BRCA <- predict(intermediate_layer_modelBRCA, BRCA.Rnaseq.expr.tpm.logmat)
```

### plot
We create tibble objects for plotting

```r
embeddings.MESO.tib = tibble(PCA=pca.MESO$li, UMAP = umap.MESO$layout, AE= AE.MESO, data=as.data.frame(colData(MESO.Rnaseq.SE) ) )
embeddings.BRCA.tib = tibble(PCA=pca.BRCA$li, UMAP = umap.BRCA$layout, AE=AE.BRCA, data=as.data.frame(colData(BRCA.Rnaseq.SE) ) )
```

#### Histopathology
Plot first 2 PCA axes. For MESO, we use the histopathological types as colors; for BRCA, we use the molecular subtypes

```r
ggMESO_PCA  = ggplot(embeddings.MESO.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=MESO.Rnaseq.SE$primary_diagnosis)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
ggBRCA_PCA  = ggplot(embeddings.BRCA.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=BRCA.Rnaseq.SE$paper_BRCA_Pathology)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
```

We plot the UMAP

```r
ggMESO_UMAP = ggplot(embeddings.MESO.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=MESO.Rnaseq.SE$primary_diagnosis)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
ggBRCA_UMAP = ggplot(embeddings.BRCA.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=BRCA.Rnaseq.SE$paper_BRCA_Pathology)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
```

We plot the autoencoders

```r
ggMESO_AE = ggplot(embeddings.MESO.tib, aes(x=AE[,1],y=AE[,2] ,col=MESO.Rnaseq.SE$primary_diagnosis)) + geom_point() + theme_bw() + xlab("AE1") + ylab("AE2")
ggBRCA_AE = ggplot(embeddings.BRCA.tib, aes(x=AE[,1],y=AE[,2] ,col=BRCA.Rnaseq.SE$paper_BRCA_Pathology)) + geom_point() + theme_bw() + xlab("AE1") + ylab("AE2")
```

We arrange the plots with patchwork

```r
((ggMESO_PCA+guides(color="none"))+ ggMESO_UMAP+guides(color="none")+ (ggMESO_AE + guides(col=guide_legend(title="Histology"))) )/((ggBRCA_PCA+guides(color="none"))+ ggBRCA_UMAP+guides(color="none")+ (ggBRCA_AE + guides(col=guide_legend(title="Histology"))) )
```

![](PCA_UMAP_files/figure-html/allplots-1.png)<!-- -->

#### Molecular groups

```r
ggMESO_PCA  = ggplot(embeddings.MESO.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=MESO_S1$iCluster_k4.5types)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
ggBRCA_PCA  = ggplot(embeddings.BRCA.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
```

We plot the UMAP

```r
ggMESO_UMAP = ggplot(embeddings.MESO.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=MESO_S1$iCluster_k4.5types)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
ggBRCA_UMAP = ggplot(embeddings.BRCA.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
```

We plot the autoencoders

```r
ggMESO_AE = ggplot(embeddings.MESO.tib, aes(x=AE[,1],y=AE[,2] ,col=MESO_S1$iCluster_k4.5types)) + geom_point() + theme_bw() + xlab("AE1") + ylab("AE2")
ggBRCA_AE = ggplot(embeddings.BRCA.tib, aes(x=AE[,1],y=AE[,2] ,col=BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50)) + geom_point() + theme_bw() + xlab("AE1") + ylab("AE2")
```

We arrange the plots with patchwork

```r
((ggMESO_PCA+guides(color="none"))+ ggMESO_UMAP+guides(color="none")+ (ggMESO_AE + guides(col=guide_legend(title="Molecular group"))) )/((ggBRCA_PCA+guides(color="none"))+ ggBRCA_UMAP+guides(color="none")+ (ggBRCA_AE + guides(col=guide_legend(title="Molecular group"))) )
```

![](PCA_UMAP_files/figure-html/allplotsM-1.png)<!-- -->

#### Clinical groups

```r
ggMESO_PCA  = ggplot(embeddings.MESO.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=MESO.Rnaseq.SE$survival_group)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
ggBRCA_PCA  = ggplot(embeddings.BRCA.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=BRCA.Rnaseq.SE$survival_group)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
```

We plot the UMAP

```r
ggMESO_UMAP = ggplot(embeddings.MESO.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=MESO.Rnaseq.SE$survival_group)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
ggBRCA_UMAP = ggplot(embeddings.BRCA.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=BRCA.Rnaseq.SE$survival_group)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
```

We plot the autoencoders

```r
ggMESO_AE = ggplot(embeddings.MESO.tib, aes(x=AE[,1],y=AE[,2] ,col=MESO.Rnaseq.SE$survival_group)) + geom_point() + theme_bw() + xlab("AE1") + ylab("AE2")
ggBRCA_AE = ggplot(embeddings.BRCA.tib, aes(x=AE[,1],y=AE[,2] ,col=BRCA.Rnaseq.SE$survival_group)) + geom_point() + theme_bw() + xlab("AE1") + ylab("AE2")
```

We arrange the plots with patchwork

```r
((ggMESO_PCA+guides(color="none"))+ ggMESO_UMAP+guides(color="none")+ (ggMESO_AE + guides(col=guide_legend(title="Survival group"))) )/((ggBRCA_PCA+guides(color="none"))+ ggBRCA_UMAP+guides(color="none")+ (ggBRCA_AE + guides(col=guide_legend(title="Survival group"))) )
```

![](PCA_UMAP_files/figure-html/allplotsC-1.png)<!-- -->

## Compute metrics
### For histopathological types

```r
## MESO
silMESO.PCA  = silhouette(as.numeric(as.factor(MESO.Rnaseq.SE$primary_diagnosis))[!is.na(MESO.Rnaseq.SE$primary_diagnosis)],
           dist = dist(embeddings.MESO.tib$PCA[!is.na(MESO.Rnaseq.SE$primary_diagnosis),]))
silMESO.UMAP = silhouette(as.numeric(as.factor(MESO.Rnaseq.SE$primary_diagnosis))[!is.na(MESO.Rnaseq.SE$primary_diagnosis)],
           dist = dist(embeddings.MESO.tib$UMAP[!is.na(MESO.Rnaseq.SE$primary_diagnosis),]))
silMESO.AE   = silhouette(as.numeric(as.factor(MESO.Rnaseq.SE$primary_diagnosis))[!is.na(MESO.Rnaseq.SE$primary_diagnosis)],
           dist = dist(embeddings.MESO.tib$AE[!is.na(MESO.Rnaseq.SE$primary_diagnosis),]))
## BRCA
silBRCA.PCA  = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$paper_BRCA_Pathology))[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Pathology)],
           dist = dist(embeddings.BRCA.tib$PCA[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Pathology),]))
silBRCA.UMAP = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$paper_BRCA_Pathology))[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Pathology)],
           dist = dist(embeddings.BRCA.tib$UMAP[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Pathology),]))
silBRCA.AE   = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$paper_BRCA_Pathology))[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Pathology)],
           dist = dist(embeddings.BRCA.tib$AE[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Pathology),]))
```

### For molecular groups

```r
## MESO
silMESO.PCA.mol  = silhouette(as.numeric(as.factor(MESO_S1$iCluster_k4.5types)), dist = dist(embeddings.MESO.tib$PCA))
silMESO.UMAP.mol = silhouette(as.numeric(as.factor(MESO_S1$iCluster_k4.5types)), dist = dist(embeddings.MESO.tib$UMAP))
silMESO.AE.mol   = silhouette(as.numeric(as.factor(MESO_S1$iCluster_k4.5types)), dist = dist(embeddings.MESO.tib$AE))
## BRCA
silBRCA.PCA.mol  = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50))[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50)],
           dist = dist(embeddings.BRCA.tib$PCA[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50),]))
silBRCA.UMAP.mol = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50))[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50)],
           dist = dist(embeddings.BRCA.tib$UMAP[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50),]))
silBRCA.AE.mol   = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50))[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50)],
           dist = dist(embeddings.BRCA.tib$AE[!is.na(BRCA.Rnaseq.SE$paper_BRCA_Subtype_PAM50),]))
```

### For clinical variables

```r
## MESO
silMESO.PCA.clin  = silhouette(as.numeric(MESO.Rnaseq.SE$survival_group)[!is.na(MESO.Rnaseq.SE$survival_group)], dist = dist(embeddings.MESO.tib$PCA[!is.na(MESO.Rnaseq.SE$survival_group),]))
silMESO.UMAP.clin = silhouette(as.numeric(MESO.Rnaseq.SE$survival_group)[!is.na(MESO.Rnaseq.SE$survival_group)], dist = dist(embeddings.MESO.tib$UMAP[!is.na(MESO.Rnaseq.SE$survival_group),]))
silMESO.AE.clin   = silhouette(as.numeric(MESO.Rnaseq.SE$survival_group)[!is.na(MESO.Rnaseq.SE$survival_group)], dist = dist(embeddings.MESO.tib$AE[!is.na(MESO.Rnaseq.SE$survival_group),]))

## BRCA
silBRCA.PCA.clin  = silhouette(as.numeric(BRCA.Rnaseq.SE$survival_group)[!is.na(BRCA.Rnaseq.SE$survival_group)],
           dist = dist(embeddings.BRCA.tib$PCA[!is.na(BRCA.Rnaseq.SE$survival_group),]))
silBRCA.UMAP.clin = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$survival_group))[!is.na(BRCA.Rnaseq.SE$survival_group)],
           dist = dist(embeddings.BRCA.tib$UMAP[!is.na(BRCA.Rnaseq.SE$survival_group),]))
silBRCA.AE.clin   = silhouette(as.numeric(as.factor(BRCA.Rnaseq.SE$survival_group))[!is.na(BRCA.Rnaseq.SE$survival_group)],
           dist = dist(embeddings.BRCA.tib$AE[!is.na(BRCA.Rnaseq.SE$survival_group),]))
```


### Plot metrics

```r
sill = tibble(embedding=rep(c("PCA","UMAP","AE"),2*3),variable=rep(c("Histological types","Molecular groups","Survival groups"),each=3*2), 
              Cohort = rep(rep(c("MESO","BRCA"),each=3),3),
              silhouette=c(mean(silMESO.PCA[,3]),mean(silMESO.UMAP[,3]),mean(silMESO.AE[,3]),
                           mean(silBRCA.PCA[,3]),mean(silBRCA.UMAP[,3]),mean(silBRCA.AE[,3]),
                           mean(silMESO.PCA.mol[,3]),mean(silMESO.UMAP.mol[,3]),mean(silMESO.AE.mol[,3]),
                           mean(silBRCA.PCA.mol[,3]),mean(silBRCA.UMAP.mol[,3]),mean(silBRCA.AE.mol[,3]),
                           mean(silMESO.PCA.clin[,3]),mean(silMESO.UMAP.clin[,3]),mean(silMESO.AE.clin[,3]),
                           mean(silBRCA.PCA.clin[,3]),mean(silBRCA.UMAP.clin[,3]),mean(silBRCA.AE.clin[,3])
                           ))

ggplot(sill, aes(x=embedding,y=silhouette,fill=Cohort) ) + geom_bar(stat="identity",beside=T,position=position_dodge()) + 
  theme_bw() +ylab("Mean silhouette width") + facet_grid(.~variable)
```

```
## Warning: Ignoring unknown parameters: beside
```

![](PCA_UMAP_files/figure-html/silplot-1.png)<!-- -->

## Session Info 

```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] readxl_1.4.1                cluster_2.1.4              
##  [3] keras_2.9.0                 ade4_1.7-18                
##  [5] umap_0.2.8.0                patchwork_1.1.1            
##  [7] SummarizedExperiment_1.24.0 Biobase_2.54.0             
##  [9] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
## [11] IRanges_2.28.0              S4Vectors_0.32.4           
## [13] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
## [15] matrixStats_0.62.0          forcats_0.5.1              
## [17] stringr_1.4.1               dplyr_1.0.10               
## [19] purrr_0.3.4                 readr_2.1.2                
## [21] tidyr_1.2.1                 tibble_3.1.8               
## [23] ggplot2_3.3.6               tidyverse_1.3.1            
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7           fs_1.5.2               lubridate_1.8.0       
##  [4] httr_1.4.4             tools_4.1.2            backports_1.4.1       
##  [7] bslib_0.4.0            utf8_1.2.2             R6_2.5.1              
## [10] DBI_1.1.3              colorspace_2.0-3       withr_2.5.0           
## [13] tidyselect_1.1.2       compiler_4.1.2         cli_3.4.1             
## [16] rvest_1.0.3            xml2_1.3.3             DelayedArray_0.20.0   
## [19] labeling_0.4.2         sass_0.4.2             scales_1.2.1          
## [22] askpass_1.1            tfruns_1.5.1           digest_0.6.29         
## [25] rmarkdown_2.16         XVector_0.34.0         base64enc_0.1-3       
## [28] pkgconfig_2.0.3        htmltools_0.5.3        highr_0.9             
## [31] dbplyr_2.2.1           fastmap_1.1.0          rlang_1.0.6           
## [34] rstudioapi_0.14        farver_2.1.1           jquerylib_0.1.4       
## [37] generics_0.1.3         jsonlite_1.8.0         tensorflow_2.9.0      
## [40] RCurl_1.98-1.8         magrittr_2.0.3.9000    GenomeInfoDbData_1.2.7
## [43] Matrix_1.5-1           Rcpp_1.0.9             munsell_0.5.0         
## [46] fansi_1.0.3            reticulate_1.26        lifecycle_1.0.2       
## [49] whisker_0.4            stringi_1.7.8          yaml_2.3.5            
## [52] MASS_7.3-58.1          zlibbioc_1.40.0        grid_4.1.2            
## [55] crayon_1.5.1           lattice_0.20-45        haven_2.5.1           
## [58] hms_1.1.2              zeallot_0.1.0          knitr_1.38            
## [61] pillar_1.8.1           reprex_2.0.2           glue_1.6.2            
## [64] evaluate_0.15          modelr_0.1.9           vctrs_0.4.1           
## [67] png_0.1-7              tzdb_0.3.0             cellranger_1.1.0      
## [70] gtable_0.3.1           openssl_2.0.3          assertthat_0.2.1      
## [73] cachem_1.0.6           xfun_0.33              broom_1.0.1           
## [76] RSpectra_0.16-1        ellipsis_0.3.2
```
