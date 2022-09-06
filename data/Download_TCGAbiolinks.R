# download with these links to solve bug in GDCprepare for RNAseq:
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
library(TCGAbiolinksGUI.data)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

# Getting the data ready
## Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(
  project = "TCGA-MESO", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

query.BRCA <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

## Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)
GDCdownload(query.BRCA)

## Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
## rsem.genes.results as values
MESO.Rnaseq.SE <- GDCprepare(query,save = TRUE)
BRCA.Rnaseq.SE <- GDCprepare(query.BRCA,save = TRUE)

## alternatively, if already prepared and saved, load from disk
load("TCGA-MESOTranscriptome_ProfilingFri_Jul_29_11:39:27_2022.RData")
MESO.Rnaseq.SE = data
load("TCGA-BRCATranscriptome_ProfilingFri_Jul_29_12:11:29_2022.RData")
BRCA.Rnaseq.SE = data

## Get expression in TPM
MESO.Rnaseq.expr.tpm = assay(MESO.Rnaseq.SE,"tpm_unstrand")
BRCA.Rnaseq.expr.tpm = assay(BRCA.Rnaseq.SE,"tpm_unstrand")

## access clinical data
colData(MESO.Rnaseq.SE)
colData(BRCA.Rnaseq.SE)
