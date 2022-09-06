---
title: "PCA_UMAP"
author: "N. Alcala"
date: "9/6/2022"
output: 
  html_document: 
    keep_md: yes
---


# R Code to run PCA and UMAP on TCGA expression data

## load libraries 

```r
library(tidyverse)
library(SummarizedExperiment)
library(patchwork)
library(umap)
library(ade4)
```

## Load expression data from TCGA

```r
load("data/TCGA-MESOTranscriptome_ProfilingFri_Jul_29_11:39:27_2022.RData")
MESO.Rnaseq.SE = data
load("data/TCGA-BRCATranscriptome_ProfilingFri_Jul_29_12:11:29_2022.RData")
BRCA.Rnaseq.SE = data
```

### Get expression in TPM

```r
MESO.Rnaseq.expr.tpm = assay(MESO.Rnaseq.SE,"tpm_unstrand")
BRCA.Rnaseq.expr.tpm = assay(BRCA.Rnaseq.SE,"tpm_unstrand")
```

## Analysis
### PCA of log(tpm+1) expression of 5000 most variable genes

```r
pca.MESO = dudi.pca(t(log(MESO.Rnaseq.expr.tpm[order(apply(MESO.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)),scannf = F,nf=10)
pca.BRCA = dudi.pca(t(log(BRCA.Rnaseq.expr.tpm[order(apply(BRCA.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)),scannf = F,nf=10)
```

### UMAP of log(tpm+1) expression of 5000 most variable genes

```r
umap.MESO = umap(t(log(MESO.Rnaseq.expr.tpm[order(apply(MESO.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)))
umap.BRCA = umap(t(log(BRCA.Rnaseq.expr.tpm[order(apply(BRCA.Rnaseq.expr.tpm,1,var),decreasing = T)[1:5000],]+1)))
```

### plot
We create tibble objects for plotting

```r
embeddings.MESO.tib = tibble(PCA=pca.MESO$li, UMAP = umap.MESO$layout, data=as.data.frame(colData(MESO.Rnaseq.SE) ) )
embeddings.BRCA.tib = tibble(PCA=pca.BRCA$li, UMAP = umap.BRCA$layout, data=as.data.frame(colData(BRCA.Rnaseq.SE) ) )
```

Plot first 2 PCA axes. For MESO, we use the histopathological types as colors; for BRCA, we use the molecular subtypes

```r
ggMESO_PCA  = ggplot(embeddings.MESO.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=data$primary_diagnosis)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
ggBRCA_PCA  = ggplot(embeddings.BRCA.tib, aes(x=PCA$Axis1,y=PCA$Axis2 ,col=data$paper_BRCA_Subtype_PAM50)) + geom_point() + theme_bw() + xlab("PC1") + ylab("PC2")
```

We plot the UMAP

```r
ggMESO_UMAP = ggplot(embeddings.MESO.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=data$primary_diagnosis)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
ggBRCA_UMAP = ggplot(embeddings.BRCA.tib, aes(x=UMAP[,1],y=UMAP[,2] ,col=data$paper_BRCA_Subtype_PAM50)) + geom_point() + theme_bw() + xlab("UMAP1") + ylab("UMAP2")
```

We arrange the plots with patchwork

```r
((ggMESO_PCA+guides(color="none"))+ ggMESO_UMAP)/((ggBRCA_PCA+guides(color="none"))+ ggBRCA_UMAP)
```

![](PCA_UMAP_files/figure-html/allplots-1.png)<!-- -->

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
##  [1] ade4_1.7-18                 umap_0.2.8.0               
##  [3] patchwork_1.1.1             SummarizedExperiment_1.24.0
##  [5] Biobase_2.54.0              GenomicRanges_1.46.1       
##  [7] GenomeInfoDb_1.30.1         IRanges_2.28.0             
##  [9] S4Vectors_0.32.4            BiocGenerics_0.40.0        
## [11] MatrixGenerics_1.6.0        matrixStats_0.62.0         
## [13] forcats_0.5.1               stringr_1.4.0              
## [15] dplyr_1.0.9                 purrr_0.3.4                
## [17] readr_2.1.2                 tidyr_1.2.0                
## [19] tibble_3.1.7                ggplot2_3.3.6              
## [21] tidyverse_1.3.1            
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7           fs_1.5.2               lubridate_1.8.0       
##  [4] httr_1.4.3             tools_4.1.2            backports_1.4.1       
##  [7] bslib_0.3.1            utf8_1.2.2             R6_2.5.1              
## [10] DBI_1.1.3              colorspace_2.0-3       withr_2.5.0           
## [13] tidyselect_1.1.2       compiler_4.1.2         cli_3.3.0             
## [16] rvest_1.0.2            xml2_1.3.3             DelayedArray_0.20.0   
## [19] labeling_0.4.2         sass_0.4.1             scales_1.2.0          
## [22] askpass_1.1            digest_0.6.29          rmarkdown_2.14        
## [25] XVector_0.34.0         pkgconfig_2.0.3        htmltools_0.5.2       
## [28] highr_0.9              dbplyr_2.2.1           fastmap_1.1.0         
## [31] rlang_1.0.4            readxl_1.4.0           rstudioapi_0.13       
## [34] farver_2.1.0           jquerylib_0.1.4        generics_0.1.2        
## [37] jsonlite_1.8.0         RCurl_1.98-1.7         magrittr_2.0.3.9000   
## [40] GenomeInfoDbData_1.2.7 Matrix_1.4-1           Rcpp_1.0.8.3          
## [43] munsell_0.5.0          fansi_1.0.3            reticulate_1.25       
## [46] lifecycle_1.0.1        stringi_1.7.6          yaml_2.3.5            
## [49] MASS_7.3-57            zlibbioc_1.40.0        grid_4.1.2            
## [52] crayon_1.5.1           lattice_0.20-45        haven_2.5.0           
## [55] hms_1.1.1              knitr_1.38             pillar_1.7.0          
## [58] reprex_2.0.1           glue_1.6.2             evaluate_0.15         
## [61] modelr_0.1.8           vctrs_0.4.1            png_0.1-7             
## [64] tzdb_0.3.0             cellranger_1.1.0       gtable_0.3.0          
## [67] openssl_2.0.2          assertthat_0.2.1       xfun_0.31             
## [70] broom_0.8.0            RSpectra_0.16-1        ellipsis_0.3.2
```
