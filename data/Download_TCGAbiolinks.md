---
title: "Download_TCGAbiolinks"
author: "N. Alcala"
date: "9/6/2022"
output: 
  html_document: 
    keep_md: yes
---


# R Code to download expression data from TCGA with TCGAbiolinks

## Install TCGAbiolink
Note that due to a bug in  GDCprepare for RNAseq, it is for now necessary to install the package with these links:

```r
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
```

## load libraries 

```r
library(TCGAbiolinksGUI.data)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
```

## Getting the data ready
### Query platform Illumina HiSeq with a list of barcode 
We first build queries, for example for the malignant pleural mesothelioma (MESO) cohort and the breast cancer cohort (BRCA)

```r
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
```

We then download the queries

```r
GDCdownload(query,directory = "data/")
GDCdownload(query.BRCA,directory = "data/")
```

### Prepare expression matrix 
Matrices contain geneID in the rows and samples (barcode) in the columns, and rsem.genes.results as values

```r
MESO.Rnaseq.SE <- GDCprepare(query,save = TRUE,directory = "data/")
BRCA.Rnaseq.SE <- GDCprepare(query.BRCA,save = TRUE,directory = "data/")
```

### alternatively, if already prepared and saved, load from disk

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

### Access clinical data

```r
colData(MESO.Rnaseq.SE)
```

```
## DataFrame with 87 rows and 62 columns
##                                             barcode      patient
##                                         <character>  <character>
## TCGA-TS-A7P1-01A-41R-A40A-07 TCGA-TS-A7P1-01A-41R.. TCGA-TS-A7P1
## TCGA-3U-A98D-01A-12R-A40A-07 TCGA-3U-A98D-01A-12R.. TCGA-3U-A98D
## TCGA-UD-AABY-01A-12R-A40A-07 TCGA-UD-AABY-01A-12R.. TCGA-UD-AABY
## TCGA-3U-A98E-01A-11R-A40A-07 TCGA-3U-A98E-01A-11R.. TCGA-3U-A98E
## TCGA-LK-A4O2-01A-11R-A34F-07 TCGA-LK-A4O2-01A-11R.. TCGA-LK-A4O2
## ...                                             ...          ...
## TCGA-TS-A7PB-01A-11R-A34F-07 TCGA-TS-A7PB-01A-11R.. TCGA-TS-A7PB
## TCGA-UD-AAC5-01A-11R-A40A-07 TCGA-UD-AAC5-01A-11R.. TCGA-UD-AAC5
## TCGA-SH-A9CT-01A-11R-A40A-07 TCGA-SH-A9CT-01A-11R.. TCGA-SH-A9CT
## TCGA-ZN-A9VO-01A-11R-A40A-07 TCGA-ZN-A9VO-01A-11R.. TCGA-ZN-A9VO
## TCGA-SC-A6LR-01A-11R-A34F-07 TCGA-SC-A6LR-01A-11R.. TCGA-SC-A6LR
##                                        sample shortLetterCode
##                                   <character>     <character>
## TCGA-TS-A7P1-01A-41R-A40A-07 TCGA-TS-A7P1-01A              TP
## TCGA-3U-A98D-01A-12R-A40A-07 TCGA-3U-A98D-01A              TP
## TCGA-UD-AABY-01A-12R-A40A-07 TCGA-UD-AABY-01A              TP
## TCGA-3U-A98E-01A-11R-A40A-07 TCGA-3U-A98E-01A              TP
## TCGA-LK-A4O2-01A-11R-A34F-07 TCGA-LK-A4O2-01A              TP
## ...                                       ...             ...
## TCGA-TS-A7PB-01A-11R-A34F-07 TCGA-TS-A7PB-01A              TP
## TCGA-UD-AAC5-01A-11R-A40A-07 TCGA-UD-AAC5-01A              TP
## TCGA-SH-A9CT-01A-11R-A40A-07 TCGA-SH-A9CT-01A              TP
## TCGA-ZN-A9VO-01A-11R-A40A-07 TCGA-ZN-A9VO-01A              TP
## TCGA-SC-A6LR-01A-11R-A34F-07 TCGA-SC-A6LR-01A              TP
##                                       definition sample_submitter_id
##                                      <character>         <character>
## TCGA-TS-A7P1-01A-41R-A40A-07 Primary solid Tumor    TCGA-TS-A7P1-01A
## TCGA-3U-A98D-01A-12R-A40A-07 Primary solid Tumor    TCGA-3U-A98D-01A
## TCGA-UD-AABY-01A-12R-A40A-07 Primary solid Tumor    TCGA-UD-AABY-01A
## TCGA-3U-A98E-01A-11R-A40A-07 Primary solid Tumor    TCGA-3U-A98E-01A
## TCGA-LK-A4O2-01A-11R-A34F-07 Primary solid Tumor    TCGA-LK-A4O2-01A
## ...                                          ...                 ...
## TCGA-TS-A7PB-01A-11R-A34F-07 Primary solid Tumor    TCGA-TS-A7PB-01A
## TCGA-UD-AAC5-01A-11R-A40A-07 Primary solid Tumor    TCGA-UD-AAC5-01A
## TCGA-SH-A9CT-01A-11R-A40A-07 Primary solid Tumor    TCGA-SH-A9CT-01A
## TCGA-ZN-A9VO-01A-11R-A40A-07 Primary solid Tumor    TCGA-ZN-A9VO-01A
## TCGA-SC-A6LR-01A-11R-A34F-07 Primary solid Tumor    TCGA-SC-A6LR-01A
##                              sample_type_id tumor_descriptor
##                                 <character>      <character>
## TCGA-TS-A7P1-01A-41R-A40A-07             01     Not Reported
## TCGA-3U-A98D-01A-12R-A40A-07             01     Not Reported
## TCGA-UD-AABY-01A-12R-A40A-07             01     Not Reported
## TCGA-3U-A98E-01A-11R-A40A-07             01     Not Reported
## TCGA-LK-A4O2-01A-11R-A34F-07             01     Not Reported
## ...                                     ...              ...
## TCGA-TS-A7PB-01A-11R-A34F-07             01     Not Reported
## TCGA-UD-AAC5-01A-11R-A40A-07             01     Not Reported
## TCGA-SH-A9CT-01A-11R-A40A-07             01     Not Reported
## TCGA-ZN-A9VO-01A-11R-A40A-07             01     Not Reported
## TCGA-SC-A6LR-01A-11R-A34F-07             01     Not Reported
##                                           sample_id   sample_type  composition
##                                         <character>   <character>  <character>
## TCGA-TS-A7P1-01A-41R-A40A-07 56956a36-3890-4558-a.. Primary Tumor Not Reported
## TCGA-3U-A98D-01A-12R-A40A-07 51b1915d-b527-4248-8.. Primary Tumor Not Reported
## TCGA-UD-AABY-01A-12R-A40A-07 9b4479ab-b8f8-4e9a-8.. Primary Tumor Not Reported
## TCGA-3U-A98E-01A-11R-A40A-07 432c1fa8-a151-4ea2-8.. Primary Tumor Not Reported
## TCGA-LK-A4O2-01A-11R-A34F-07 e3fdf51e-9057-4197-9.. Primary Tumor Not Reported
## ...                                             ...           ...          ...
## TCGA-TS-A7PB-01A-11R-A34F-07 383a44cb-be6b-4e3a-9.. Primary Tumor Not Reported
## TCGA-UD-AAC5-01A-11R-A40A-07 d644008a-8c2c-4d6a-b.. Primary Tumor Not Reported
## TCGA-SH-A9CT-01A-11R-A40A-07 d1fc0218-4d71-4ab4-a.. Primary Tumor Not Reported
## TCGA-ZN-A9VO-01A-11R-A40A-07 d39da5a5-2c0d-40e1-b.. Primary Tumor Not Reported
## TCGA-SC-A6LR-01A-11R-A34F-07 cecee06a-8161-4b3b-a.. Primary Tumor Not Reported
##                              days_to_collection       state initial_weight
##                                       <integer> <character>      <numeric>
## TCGA-TS-A7P1-01A-41R-A40A-07               1161    released            570
## TCGA-3U-A98D-01A-12R-A40A-07                628    released            230
## TCGA-UD-AABY-01A-12R-A40A-07               5129    released            190
## TCGA-3U-A98E-01A-11R-A40A-07               1321    released            390
## TCGA-LK-A4O2-01A-11R-A34F-07                630    released            310
## ...                                         ...         ...            ...
## TCGA-TS-A7PB-01A-11R-A34F-07                869    released            800
## TCGA-UD-AAC5-01A-11R-A40A-07               2910    released             60
## TCGA-SH-A9CT-01A-11R-A40A-07                187    released            570
## TCGA-ZN-A9VO-01A-11R-A40A-07                765    released             90
## TCGA-SC-A6LR-01A-11R-A34F-07                443    released            480
##                               pathology_report_uuid submitter_id oct_embedded
##                                         <character>  <character>  <character>
## TCGA-TS-A7P1-01A-41R-A40A-07 0FAC8AE9-5E06-44F7-9.. TCGA-TS-A7P1         true
## TCGA-3U-A98D-01A-12R-A40A-07 8F55DF05-CE3D-42FF-9.. TCGA-3U-A98D         true
## TCGA-UD-AABY-01A-12R-A40A-07 BED6A020-507B-41BD-8.. TCGA-UD-AABY         true
## TCGA-3U-A98E-01A-11R-A40A-07 11FB2902-12C0-4F01-A.. TCGA-3U-A98E         true
## TCGA-LK-A4O2-01A-11R-A34F-07 D4751900-6378-445D-B.. TCGA-LK-A4O2         true
## ...                                             ...          ...          ...
## TCGA-TS-A7PB-01A-11R-A34F-07 10F97F63-EE6F-448D-9.. TCGA-TS-A7PB         true
## TCGA-UD-AAC5-01A-11R-A40A-07 5BB63340-2A48-4ACB-A.. TCGA-UD-AAC5         true
## TCGA-SH-A9CT-01A-11R-A40A-07 D7747749-3C50-4563-A.. TCGA-SH-A9CT         true
## TCGA-ZN-A9VO-01A-11R-A40A-07 D10B7DE3-1AE6-42A3-B.. TCGA-ZN-A9VO         true
## TCGA-SC-A6LR-01A-11R-A34F-07 E1F77E03-1351-48F3-9.. TCGA-SC-A6LR         true
##                                is_ffpe  tissue_type synchronous_malignancy
##                              <logical>  <character>            <character>
## TCGA-TS-A7P1-01A-41R-A40A-07     FALSE Not Reported                     No
## TCGA-3U-A98D-01A-12R-A40A-07     FALSE Not Reported                     No
## TCGA-UD-AABY-01A-12R-A40A-07     FALSE Not Reported                     No
## TCGA-3U-A98E-01A-11R-A40A-07     FALSE Not Reported           Not Reported
## TCGA-LK-A4O2-01A-11R-A34F-07     FALSE Not Reported                     No
## ...                                ...          ...                    ...
## TCGA-TS-A7PB-01A-11R-A34F-07     FALSE Not Reported                     No
## TCGA-UD-AAC5-01A-11R-A40A-07     FALSE Not Reported                     No
## TCGA-SH-A9CT-01A-11R-A40A-07     FALSE Not Reported                     No
## TCGA-ZN-A9VO-01A-11R-A40A-07     FALSE Not Reported                     No
## TCGA-SC-A6LR-01A-11R-A34F-07     FALSE Not Reported                     No
##                              ajcc_pathologic_stage days_to_diagnosis
##                                        <character>         <integer>
## TCGA-TS-A7P1-01A-41R-A40A-07               Stage I                 0
## TCGA-3U-A98D-01A-12R-A40A-07             Stage III                 0
## TCGA-UD-AABY-01A-12R-A40A-07             Stage III                 0
## TCGA-3U-A98E-01A-11R-A40A-07              Stage IV                 0
## TCGA-LK-A4O2-01A-11R-A34F-07             Stage III                 0
## ...                                            ...               ...
## TCGA-TS-A7PB-01A-11R-A34F-07             Stage III                 0
## TCGA-UD-AAC5-01A-11R-A40A-07               Stage I                 0
## TCGA-SH-A9CT-01A-11R-A40A-07              Stage II                 0
## TCGA-ZN-A9VO-01A-11R-A40A-07              Stage IV                 0
## TCGA-SC-A6LR-01A-11R-A34F-07             Stage III                 0
##                                                                                                                                 treatments
##                                                                                                                                     <list>
## TCGA-TS-A7P1-01A-41R-A40A-07 NA:2019-07-31T22:07:41...:360dcca6-f2d1-575d-b..:...,NA:2019-07-31T22:07:41...:39ddc85e-012f-53e4-a..:...,...
## TCGA-3U-A98D-01A-12R-A40A-07 NA:2019-07-31T22:08:36...:418b05df-ef07-581a-a..:...,NA:2019-07-31T22:08:36...:d0f85001-3636-57e8-9..:...,...
## TCGA-UD-AABY-01A-12R-A40A-07 NA:2019-07-31T22:08:15...:37bf541d-61c4-5646-b..:...,NA:2019-07-31T22:08:15...:8106affc-430c-5c20-b..:...,...
## TCGA-3U-A98E-01A-11R-A40A-07 NA:2019-07-31T22:09:09...:28503648-d450-5739-a..:...,NA:2019-07-31T22:09:09...:b54f6125-585a-5d55-a..:...,...
## TCGA-LK-A4O2-01A-11R-A34F-07 NA:2019-07-31T22:07:59...:12dd3210-92f4-53fc-9..:...,NA:2019-07-31T22:07:59...:77d564cd-da97-5cbd-8..:...,...
## ...                                                                                                                                    ...
## TCGA-TS-A7PB-01A-11R-A34F-07 NA:2019-07-31T16:10:18...:3515ac50-0449-5457-a..:...,NA:2019-07-31T16:10:18...:b0d915d5-853c-5c86-8..:...,...
## TCGA-UD-AAC5-01A-11R-A40A-07                                         NA:NA:5c6a9a88-7214-538d-a..:...,NA:NA:d877bab9-d40e-5740-9..:...,...
## TCGA-SH-A9CT-01A-11R-A40A-07 NA:2019-07-31T22:09:52...:719c49bf-347f-5c26-a..:...,NA:2019-07-31T22:09:52...:a2b61c93-f7d3-53cd-9..:...,...
## TCGA-ZN-A9VO-01A-11R-A40A-07 NA:2019-07-31T16:31:37...:9be96e8d-4806-5ea9-9..:...,NA:2019-07-31T16:31:37...:c38837db-2d99-54b8-b..:...,...
## TCGA-SC-A6LR-01A-11R-A34F-07 NA:2019-07-31T17:43:07...:679a8723-0b10-55e8-b..:...,NA:2019-07-31T17:43:07...:bb04be49-b03d-5f2f-b..:...,...
##                              last_known_disease_status
##                                            <character>
## TCGA-TS-A7P1-01A-41R-A40A-07              not reported
## TCGA-3U-A98D-01A-12R-A40A-07              not reported
## TCGA-UD-AABY-01A-12R-A40A-07              not reported
## TCGA-3U-A98E-01A-11R-A40A-07              not reported
## TCGA-LK-A4O2-01A-11R-A34F-07              not reported
## ...                                                ...
## TCGA-TS-A7PB-01A-11R-A34F-07              not reported
## TCGA-UD-AAC5-01A-11R-A40A-07              not reported
## TCGA-SH-A9CT-01A-11R-A40A-07              not reported
## TCGA-ZN-A9VO-01A-11R-A40A-07              not reported
## TCGA-SC-A6LR-01A-11R-A34F-07              not reported
##                              tissue_or_organ_of_origin days_to_last_follow_up
##                                            <character>              <integer>
## TCGA-TS-A7P1-01A-41R-A40A-07               Pleura, NOS                     NA
## TCGA-3U-A98D-01A-12R-A40A-07               Pleura, NOS                    125
## TCGA-UD-AABY-01A-12R-A40A-07               Pleura, NOS                     NA
## TCGA-3U-A98E-01A-11R-A40A-07               Pleura, NOS                   1916
## TCGA-LK-A4O2-01A-11R-A34F-07               Pleura, NOS                     NA
## ...                                                ...                    ...
## TCGA-TS-A7PB-01A-11R-A34F-07               Pleura, NOS                     NA
## TCGA-UD-AAC5-01A-11R-A40A-07               Pleura, NOS                     NA
## TCGA-SH-A9CT-01A-11R-A40A-07               Pleura, NOS                    256
## TCGA-ZN-A9VO-01A-11R-A40A-07               Pleura, NOS                    261
## TCGA-SC-A6LR-01A-11R-A34F-07               Pleura, NOS                   1168
##                              age_at_diagnosis      primary_diagnosis
##                                     <integer>            <character>
## TCGA-TS-A7P1-01A-41R-A40A-07            19191 Epithelioid mesothel..
## TCGA-3U-A98D-01A-12R-A40A-07            29500 Mesothelioma, biphas..
## TCGA-UD-AABY-01A-12R-A40A-07            18315 Epithelioid mesothel..
## TCGA-3U-A98E-01A-11R-A40A-07            26264 Epithelioid mesothel..
## TCGA-LK-A4O2-01A-11R-A34F-07            20757 Epithelioid mesothel..
## ...                                       ...                    ...
## TCGA-TS-A7PB-01A-11R-A34F-07            22822 Mesothelioma, biphas..
## TCGA-UD-AAC5-01A-11R-A40A-07            24314 Epithelioid mesothel..
## TCGA-SH-A9CT-01A-11R-A40A-07            24528 Epithelioid mesothel..
## TCGA-ZN-A9VO-01A-11R-A40A-07            13243 Epithelioid mesothel..
## TCGA-SC-A6LR-01A-11R-A34F-07            26773 Epithelioid mesothel..
##                              prior_malignancy year_of_diagnosis prior_treatment
##                                   <character>         <integer>     <character>
## TCGA-TS-A7P1-01A-41R-A40A-07               no              2010              No
## TCGA-3U-A98D-01A-12R-A40A-07               no              2012              No
## TCGA-UD-AABY-01A-12R-A40A-07               no              1999              No
## TCGA-3U-A98E-01A-11R-A40A-07              yes              2010              No
## TCGA-LK-A4O2-01A-11R-A34F-07               no              2011              No
## ...                                       ...               ...             ...
## TCGA-TS-A7PB-01A-11R-A34F-07               no              2011              No
## TCGA-UD-AAC5-01A-11R-A40A-07               no              2005             Yes
## TCGA-SH-A9CT-01A-11R-A40A-07               no              2013              No
## TCGA-ZN-A9VO-01A-11R-A40A-07               no              2011              No
## TCGA-SC-A6LR-01A-11R-A34F-07               no              2012              No
##                              ajcc_staging_system_edition ajcc_pathologic_t
##                                              <character>       <character>
## TCGA-TS-A7P1-01A-41R-A40A-07                         6th                T1
## TCGA-3U-A98D-01A-12R-A40A-07                         7th                T2
## TCGA-UD-AABY-01A-12R-A40A-07                         5th                T3
## TCGA-3U-A98E-01A-11R-A40A-07                         7th                T4
## TCGA-LK-A4O2-01A-11R-A34F-07                         7th                T3
## ...                                                  ...               ...
## TCGA-TS-A7PB-01A-11R-A34F-07                         6th                T3
## TCGA-UD-AAC5-01A-11R-A40A-07                         6th                T1
## TCGA-SH-A9CT-01A-11R-A40A-07                         7th                T2
## TCGA-ZN-A9VO-01A-11R-A40A-07                         7th                T4
## TCGA-SC-A6LR-01A-11R-A34F-07                         7th               T1b
##                               morphology ajcc_pathologic_n ajcc_pathologic_m
##                              <character>       <character>       <character>
## TCGA-TS-A7P1-01A-41R-A40A-07      9052/3                N0                MX
## TCGA-3U-A98D-01A-12R-A40A-07      9053/3                N2                M0
## TCGA-UD-AABY-01A-12R-A40A-07      9052/3                N0                M0
## TCGA-3U-A98E-01A-11R-A40A-07      9052/3                N2                M0
## TCGA-LK-A4O2-01A-11R-A34F-07      9052/3                N0                M0
## ...                                  ...               ...               ...
## TCGA-TS-A7PB-01A-11R-A34F-07      9053/3                N2                MX
## TCGA-UD-AAC5-01A-11R-A40A-07      9052/3                N0                M0
## TCGA-SH-A9CT-01A-11R-A40A-07      9052/3                N0                M0
## TCGA-ZN-A9VO-01A-11R-A40A-07      9052/3                N2                M1
## TCGA-SC-A6LR-01A-11R-A34F-07      9052/3                N2                M0
##                              classification_of_tumor           diagnosis_id
##                                          <character>            <character>
## TCGA-TS-A7P1-01A-41R-A40A-07            not reported 4f603cc3-3d1c-57bf-8..
## TCGA-3U-A98D-01A-12R-A40A-07            not reported 4c06d5fe-4318-5e12-9..
## TCGA-UD-AABY-01A-12R-A40A-07            not reported c3cd5299-18cf-5b3f-b..
## TCGA-3U-A98E-01A-11R-A40A-07            not reported 300ab384-85d4-5d91-b..
## TCGA-LK-A4O2-01A-11R-A34F-07            not reported a6bd6e1c-3c8f-50f2-b..
## ...                                              ...                    ...
## TCGA-TS-A7PB-01A-11R-A34F-07            not reported f408d64d-5d16-5e81-a..
## TCGA-UD-AAC5-01A-11R-A40A-07            not reported d5f361a5-aa90-53fd-b..
## TCGA-SH-A9CT-01A-11R-A40A-07            not reported 7354609e-52b6-5cd2-a..
## TCGA-ZN-A9VO-01A-11R-A40A-07            not reported 772e8ca4-6034-5e23-b..
## TCGA-SC-A6LR-01A-11R-A34F-07            not reported a7ff271e-a068-5a93-a..
##                              icd_10_code site_of_resection_or_biopsy
##                              <character>                 <character>
## TCGA-TS-A7P1-01A-41R-A40A-07       C45.0                 Pleura, NOS
## TCGA-3U-A98D-01A-12R-A40A-07       C45.0                 Pleura, NOS
## TCGA-UD-AABY-01A-12R-A40A-07       C45.0                 Pleura, NOS
## TCGA-3U-A98E-01A-11R-A40A-07       C45.0                 Pleura, NOS
## TCGA-LK-A4O2-01A-11R-A34F-07       C45.0                 Pleura, NOS
## ...                                  ...                         ...
## TCGA-TS-A7PB-01A-11R-A34F-07       C45.0                 Pleura, NOS
## TCGA-UD-AAC5-01A-11R-A40A-07       C45.0                 Pleura, NOS
## TCGA-SH-A9CT-01A-11R-A40A-07       C45.0                 Pleura, NOS
## TCGA-ZN-A9VO-01A-11R-A40A-07       C45.0                 Pleura, NOS
## TCGA-SC-A6LR-01A-11R-A34F-07       C45.0                 Pleura, NOS
##                               tumor_grade progression_or_recurrence
##                               <character>               <character>
## TCGA-TS-A7P1-01A-41R-A40A-07 not reported              not reported
## TCGA-3U-A98D-01A-12R-A40A-07 not reported              not reported
## TCGA-UD-AABY-01A-12R-A40A-07 not reported              not reported
## TCGA-3U-A98E-01A-11R-A40A-07 not reported              not reported
## TCGA-LK-A4O2-01A-11R-A34F-07 not reported              not reported
## ...                                   ...                       ...
## TCGA-TS-A7PB-01A-11R-A34F-07 not reported              not reported
## TCGA-UD-AAC5-01A-11R-A40A-07 not reported              not reported
## TCGA-SH-A9CT-01A-11R-A40A-07 not reported              not reported
## TCGA-ZN-A9VO-01A-11R-A40A-07 not reported              not reported
## TCGA-SC-A6LR-01A-11R-A34F-07 not reported              not reported
##                              alcohol_history            exposure_id        race
##                                  <character>            <character> <character>
## TCGA-TS-A7P1-01A-41R-A40A-07    Not Reported 41088f57-1a9a-526d-8..       white
## TCGA-3U-A98D-01A-12R-A40A-07    Not Reported 7bac06ec-d8e1-5dea-9..       white
## TCGA-UD-AABY-01A-12R-A40A-07    Not Reported 04f1447d-c71a-570d-b..       white
## TCGA-3U-A98E-01A-11R-A40A-07    Not Reported 41dec49d-c305-54c1-8..       white
## TCGA-LK-A4O2-01A-11R-A34F-07    Not Reported 4a9c1158-67f4-522b-8..       white
## ...                                      ...                    ...         ...
## TCGA-TS-A7PB-01A-11R-A34F-07    Not Reported c1fb9976-e25c-5134-b..       white
## TCGA-UD-AAC5-01A-11R-A40A-07    Not Reported bfc419bc-004b-51bd-8..       white
## TCGA-SH-A9CT-01A-11R-A40A-07    Not Reported bccb1aae-3daa-5391-8..       white
## TCGA-ZN-A9VO-01A-11R-A40A-07    Not Reported 0d921ec5-d3cf-5f9e-a..       white
## TCGA-SC-A6LR-01A-11R-A34F-07    Not Reported b8230a1a-48e1-5cd0-b..       white
##                                   gender              ethnicity vital_status
##                              <character>            <character>  <character>
## TCGA-TS-A7P1-01A-41R-A40A-07        male           not reported         Dead
## TCGA-3U-A98D-01A-12R-A40A-07        male           not reported         Dead
## TCGA-UD-AABY-01A-12R-A40A-07        male not hispanic or latino         Dead
## TCGA-3U-A98E-01A-11R-A40A-07      female not hispanic or latino        Alive
## TCGA-LK-A4O2-01A-11R-A34F-07      female not hispanic or latino         Dead
## ...                                  ...                    ...          ...
## TCGA-TS-A7PB-01A-11R-A34F-07        male not hispanic or latino         Dead
## TCGA-UD-AAC5-01A-11R-A40A-07        male not hispanic or latino         Dead
## TCGA-SH-A9CT-01A-11R-A40A-07        male not hispanic or latino         Dead
## TCGA-ZN-A9VO-01A-11R-A40A-07        male not hispanic or latino        Alive
## TCGA-SC-A6LR-01A-11R-A34F-07        male not hispanic or latino        Alive
##                              age_at_index days_to_birth year_of_birth
##                                 <integer>     <integer>     <integer>
## TCGA-TS-A7P1-01A-41R-A40A-07           52        -19191          1958
## TCGA-3U-A98D-01A-12R-A40A-07           80        -29500          1932
## TCGA-UD-AABY-01A-12R-A40A-07           50        -18315          1949
## TCGA-3U-A98E-01A-11R-A40A-07           71        -26264          1939
## TCGA-LK-A4O2-01A-11R-A34F-07           56        -20757          1955
## ...                                   ...           ...           ...
## TCGA-TS-A7PB-01A-11R-A34F-07           62        -22822          1949
## TCGA-UD-AAC5-01A-11R-A40A-07           66        -24314          1939
## TCGA-SH-A9CT-01A-11R-A40A-07           67        -24528          1946
## TCGA-ZN-A9VO-01A-11R-A40A-07           36        -13243          1975
## TCGA-SC-A6LR-01A-11R-A34F-07           73        -26773          1939
##                                      demographic_id days_to_death year_of_death
##                                         <character>     <integer>     <integer>
## TCGA-TS-A7P1-01A-41R-A40A-07 75710032-edb5-5dc7-9..           563          2011
## TCGA-3U-A98D-01A-12R-A40A-07 fcb6364b-5aeb-5c67-9..           142            NA
## TCGA-UD-AABY-01A-12R-A40A-07 a5a63a50-bcb8-5b5b-9..           843          2001
## TCGA-3U-A98E-01A-11R-A40A-07 c8b9ff26-8639-50f7-8..            NA            NA
## TCGA-LK-A4O2-01A-11R-A34F-07 1c773221-4d6f-54cb-8..            57          2011
## ...                                             ...           ...           ...
## TCGA-TS-A7PB-01A-11R-A34F-07 f7be51ca-2b04-5edd-8..           982          2013
## TCGA-UD-AAC5-01A-11R-A40A-07 a57d246d-224c-5983-a..           253          2005
## TCGA-SH-A9CT-01A-11R-A40A-07 49ec427e-b0ad-5082-9..           449            NA
## TCGA-ZN-A9VO-01A-11R-A40A-07 0a4e6011-32d6-587f-a..            NA            NA
## TCGA-SC-A6LR-01A-11R-A34F-07 cf8b8ef8-dd39-55cd-8..            NA            NA
##                              bcr_patient_barcode
##                                      <character>
## TCGA-TS-A7P1-01A-41R-A40A-07    TCGA-TS-A7P1-01A
## TCGA-3U-A98D-01A-12R-A40A-07    TCGA-3U-A98D-01A
## TCGA-UD-AABY-01A-12R-A40A-07    TCGA-UD-AABY-01A
## TCGA-3U-A98E-01A-11R-A40A-07    TCGA-3U-A98E-01A
## TCGA-LK-A4O2-01A-11R-A34F-07    TCGA-LK-A4O2-01A
## ...                                          ...
## TCGA-TS-A7PB-01A-11R-A34F-07    TCGA-TS-A7PB-01A
## TCGA-UD-AAC5-01A-11R-A40A-07    TCGA-UD-AAC5-01A
## TCGA-SH-A9CT-01A-11R-A40A-07    TCGA-SH-A9CT-01A
## TCGA-ZN-A9VO-01A-11R-A40A-07    TCGA-ZN-A9VO-01A
## TCGA-SC-A6LR-01A-11R-A34F-07    TCGA-SC-A6LR-01A
##                                                          primary_site
##                                                                <list>
## TCGA-TS-A7P1-01A-41R-A40A-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-3U-A98D-01A-12R-A40A-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-UD-AABY-01A-12R-A40A-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-3U-A98E-01A-11R-A40A-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-LK-A4O2-01A-11R-A34F-07 Heart, mediastinum, ..,Bronchus and lung
## ...                                                               ...
## TCGA-TS-A7PB-01A-11R-A34F-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-UD-AAC5-01A-11R-A40A-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-SH-A9CT-01A-11R-A40A-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-ZN-A9VO-01A-11R-A40A-07 Heart, mediastinum, ..,Bronchus and lung
## TCGA-SC-A6LR-01A-11R-A34F-07 Heart, mediastinum, ..,Bronchus and lung
##                               project_id          disease_type         name
##                              <character>                <list>  <character>
## TCGA-TS-A7P1-01A-41R-A40A-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-3U-A98D-01A-12R-A40A-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-UD-AABY-01A-12R-A40A-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-3U-A98E-01A-11R-A40A-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-LK-A4O2-01A-11R-A34F-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## ...                                  ...                   ...          ...
## TCGA-TS-A7PB-01A-11R-A34F-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-UD-AAC5-01A-11R-A40A-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-SH-A9CT-01A-11R-A40A-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-ZN-A9VO-01A-11R-A40A-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
## TCGA-SC-A6LR-01A-11R-A34F-07   TCGA-MESO Mesothelial Neoplasms Mesothelioma
##                              releasable  released      sample.aux
##                               <logical> <logical>     <character>
## TCGA-TS-A7P1-01A-41R-A40A-07       TRUE      TRUE TCGA-TS-A7P1-01
## TCGA-3U-A98D-01A-12R-A40A-07       TRUE      TRUE TCGA-3U-A98D-01
## TCGA-UD-AABY-01A-12R-A40A-07       TRUE      TRUE TCGA-UD-AABY-01
## TCGA-3U-A98E-01A-11R-A40A-07       TRUE      TRUE TCGA-3U-A98E-01
## TCGA-LK-A4O2-01A-11R-A34F-07       TRUE      TRUE TCGA-LK-A4O2-01
## ...                                 ...       ...             ...
## TCGA-TS-A7PB-01A-11R-A34F-07       TRUE      TRUE TCGA-TS-A7PB-01
## TCGA-UD-AAC5-01A-11R-A40A-07       TRUE      TRUE TCGA-UD-AAC5-01
## TCGA-SH-A9CT-01A-11R-A40A-07       TRUE      TRUE TCGA-SH-A9CT-01
## TCGA-ZN-A9VO-01A-11R-A40A-07       TRUE      TRUE TCGA-ZN-A9VO-01
## TCGA-SC-A6LR-01A-11R-A34F-07       TRUE      TRUE TCGA-SC-A6LR-01
```

```r
colData(BRCA.Rnaseq.SE)
```

```
## DataFrame with 1226 rows and 85 columns
##                                             barcode      patient
##                                         <character>  <character>
## TCGA-AC-A8OP-01A-11R-A36F-07 TCGA-AC-A8OP-01A-11R.. TCGA-AC-A8OP
## TCGA-D8-A1XU-01A-11R-A14M-07 TCGA-D8-A1XU-01A-11R.. TCGA-D8-A1XU
## TCGA-BH-A18L-01A-32R-A12D-07 TCGA-BH-A18L-01A-32R.. TCGA-BH-A18L
## TCGA-B6-A0IK-01A-12R-A056-07 TCGA-B6-A0IK-01A-12R.. TCGA-B6-A0IK
## TCGA-BH-A18L-11A-42R-A12D-07 TCGA-BH-A18L-11A-42R.. TCGA-BH-A18L
## ...                                             ...          ...
## TCGA-GM-A2DN-01A-11R-A180-07 TCGA-GM-A2DN-01A-11R.. TCGA-GM-A2DN
## TCGA-B6-A0I1-01A-11R-A21T-07 TCGA-B6-A0I1-01A-11R.. TCGA-B6-A0I1
## TCGA-EW-A1IW-01A-11R-A13Q-07 TCGA-EW-A1IW-01A-11R.. TCGA-EW-A1IW
## TCGA-BH-A42V-01A-11R-A24H-07 TCGA-BH-A42V-01A-11R.. TCGA-BH-A42V
## TCGA-AR-A251-01A-12R-A169-07 TCGA-AR-A251-01A-12R.. TCGA-AR-A251
##                                        sample shortLetterCode
##                                   <character>     <character>
## TCGA-AC-A8OP-01A-11R-A36F-07 TCGA-AC-A8OP-01A              TP
## TCGA-D8-A1XU-01A-11R-A14M-07 TCGA-D8-A1XU-01A              TP
## TCGA-BH-A18L-01A-32R-A12D-07 TCGA-BH-A18L-01A              TP
## TCGA-B6-A0IK-01A-12R-A056-07 TCGA-B6-A0IK-01A              TP
## TCGA-BH-A18L-11A-42R-A12D-07 TCGA-BH-A18L-11A              NT
## ...                                       ...             ...
## TCGA-GM-A2DN-01A-11R-A180-07 TCGA-GM-A2DN-01A              TP
## TCGA-B6-A0I1-01A-11R-A21T-07 TCGA-B6-A0I1-01A              TP
## TCGA-EW-A1IW-01A-11R-A13Q-07 TCGA-EW-A1IW-01A              TP
## TCGA-BH-A42V-01A-11R-A24H-07 TCGA-BH-A42V-01A              TP
## TCGA-AR-A251-01A-12R-A169-07 TCGA-AR-A251-01A              TP
##                                       definition sample_submitter_id
##                                      <character>         <character>
## TCGA-AC-A8OP-01A-11R-A36F-07 Primary solid Tumor    TCGA-AC-A8OP-01A
## TCGA-D8-A1XU-01A-11R-A14M-07 Primary solid Tumor    TCGA-D8-A1XU-01A
## TCGA-BH-A18L-01A-32R-A12D-07 Primary solid Tumor    TCGA-BH-A18L-01A
## TCGA-B6-A0IK-01A-12R-A056-07 Primary solid Tumor    TCGA-B6-A0IK-01A
## TCGA-BH-A18L-11A-42R-A12D-07 Solid Tissue Normal    TCGA-BH-A18L-11A
## ...                                          ...                 ...
## TCGA-GM-A2DN-01A-11R-A180-07 Primary solid Tumor    TCGA-GM-A2DN-01A
## TCGA-B6-A0I1-01A-11R-A21T-07 Primary solid Tumor    TCGA-B6-A0I1-01A
## TCGA-EW-A1IW-01A-11R-A13Q-07 Primary solid Tumor    TCGA-EW-A1IW-01A
## TCGA-BH-A42V-01A-11R-A24H-07 Primary solid Tumor    TCGA-BH-A42V-01A
## TCGA-AR-A251-01A-12R-A169-07 Primary solid Tumor    TCGA-AR-A251-01A
##                              sample_type_id              sample_id
##                                 <character>            <character>
## TCGA-AC-A8OP-01A-11R-A36F-07             01 c32e8a89-bc6b-4e8a-8..
## TCGA-D8-A1XU-01A-11R-A14M-07             01 4fbe5e67-d18a-48d5-b..
## TCGA-BH-A18L-01A-32R-A12D-07             01 796e565c-3841-425a-b..
## TCGA-B6-A0IK-01A-12R-A056-07             01 357e1eca-0f5c-49d6-b..
## TCGA-BH-A18L-11A-42R-A12D-07             11 67f32715-fea8-4b4e-9..
## ...                                     ...                    ...
## TCGA-GM-A2DN-01A-11R-A180-07             01 e46b29b9-c697-4739-b..
## TCGA-B6-A0I1-01A-11R-A21T-07             01 2a5f4696-6e02-4fe9-b..
## TCGA-EW-A1IW-01A-11R-A13Q-07             01 bfd30c90-9b30-4e7b-9..
## TCGA-BH-A42V-01A-11R-A24H-07             01 2cf478d5-5400-4724-b..
## TCGA-AR-A251-01A-12R-A169-07             01 8cf391f9-d743-408f-8..
##                                      sample_type days_to_collection       state
##                                      <character>          <integer> <character>
## TCGA-AC-A8OP-01A-11R-A36F-07       Primary Tumor                921    released
## TCGA-D8-A1XU-01A-11R-A14M-07       Primary Tumor                102    released
## TCGA-BH-A18L-01A-32R-A12D-07       Primary Tumor               2669    released
## TCGA-B6-A0IK-01A-12R-A056-07       Primary Tumor               5724    released
## TCGA-BH-A18L-11A-42R-A12D-07 Solid Tissue Normal               2669    released
## ...                                          ...                ...         ...
## TCGA-GM-A2DN-01A-11R-A180-07       Primary Tumor               2436    released
## TCGA-B6-A0I1-01A-11R-A21T-07       Primary Tumor               4248    released
## TCGA-EW-A1IW-01A-11R-A13Q-07       Primary Tumor                213    released
## TCGA-BH-A42V-01A-11R-A24H-07       Primary Tumor                197    released
## TCGA-AR-A251-01A-12R-A169-07       Primary Tumor               1611    released
##                              initial_weight  pathology_report_uuid submitter_id
##                                   <numeric>            <character>  <character>
## TCGA-AC-A8OP-01A-11R-A36F-07             80 F4F5C477-30BB-41EE-B.. TCGA-AC-A8OP
## TCGA-D8-A1XU-01A-11R-A14M-07            210 845F8FCF-CF3C-4CEF-B.. TCGA-D8-A1XU
## TCGA-BH-A18L-01A-32R-A12D-07            260 E4DF4771-074D-4870-A.. TCGA-BH-A18L
## TCGA-B6-A0IK-01A-12R-A056-07            140 3A38A97C-2CBB-4802-9.. TCGA-B6-A0IK
## TCGA-BH-A18L-11A-42R-A12D-07            530                     NA TCGA-BH-A18L
## ...                                     ...                    ...          ...
## TCGA-GM-A2DN-01A-11R-A180-07            210 8829274B-EB82-4155-A.. TCGA-GM-A2DN
## TCGA-B6-A0I1-01A-11R-A21T-07            100 1E498731-FE74-493E-A.. TCGA-B6-A0I1
## TCGA-EW-A1IW-01A-11R-A13Q-07             90 D4194F52-3E15-4105-9.. TCGA-EW-A1IW
## TCGA-BH-A42V-01A-11R-A24H-07             60 5BB5CEFD-CC07-4451-A.. TCGA-BH-A42V
## TCGA-AR-A251-01A-12R-A169-07            450 F5B17108-8904-48E7-9.. TCGA-AR-A251
##                              oct_embedded   is_ffpe  tissue_type
##                               <character> <logical>  <character>
## TCGA-AC-A8OP-01A-11R-A36F-07        false     FALSE Not Reported
## TCGA-D8-A1XU-01A-11R-A14M-07        false     FALSE Not Reported
## TCGA-BH-A18L-01A-32R-A12D-07         true     FALSE Not Reported
## TCGA-B6-A0IK-01A-12R-A056-07         true     FALSE Not Reported
## TCGA-BH-A18L-11A-42R-A12D-07         true     FALSE Not Reported
## ...                                   ...       ...          ...
## TCGA-GM-A2DN-01A-11R-A180-07        false     FALSE Not Reported
## TCGA-B6-A0I1-01A-11R-A21T-07         true     FALSE Not Reported
## TCGA-EW-A1IW-01A-11R-A13Q-07         true     FALSE Not Reported
## TCGA-BH-A42V-01A-11R-A24H-07         true     FALSE Not Reported
## TCGA-AR-A251-01A-12R-A169-07         true     FALSE Not Reported
##                              synchronous_malignancy ajcc_pathologic_stage
##                                         <character>           <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                     No              Stage IA
## TCGA-D8-A1XU-01A-11R-A14M-07                     No              Stage IA
## TCGA-BH-A18L-01A-32R-A12D-07                     No            Stage IIIA
## TCGA-B6-A0IK-01A-12R-A056-07                     No            Stage IIIB
## TCGA-BH-A18L-11A-42R-A12D-07                     No            Stage IIIA
## ...                                             ...                   ...
## TCGA-GM-A2DN-01A-11R-A180-07                     No             Stage IIA
## TCGA-B6-A0I1-01A-11R-A21T-07                     No             Stage IIA
## TCGA-EW-A1IW-01A-11R-A13Q-07                     No             Stage IIB
## TCGA-BH-A42V-01A-11R-A24H-07                     No              Stage IB
## TCGA-AR-A251-01A-12R-A169-07                     No            Stage IIIA
##                              days_to_diagnosis
##                                      <integer>
## TCGA-AC-A8OP-01A-11R-A36F-07                 0
## TCGA-D8-A1XU-01A-11R-A14M-07                 0
## TCGA-BH-A18L-01A-32R-A12D-07                 0
## TCGA-B6-A0IK-01A-12R-A056-07                 0
## TCGA-BH-A18L-11A-42R-A12D-07                 0
## ...                                        ...
## TCGA-GM-A2DN-01A-11R-A180-07                 0
## TCGA-B6-A0I1-01A-11R-A21T-07                 0
## TCGA-EW-A1IW-01A-11R-A13Q-07                 0
## TCGA-BH-A42V-01A-11R-A24H-07                 0
## TCGA-AR-A251-01A-12R-A169-07                 0
##                                                                                                                                 treatments
##                                                                                                                                     <list>
## TCGA-AC-A8OP-01A-11R-A36F-07                                         NA:NA:1fb7b3f7-a087-5f9b-9..:...,NA:NA:2c052d2c-abe7-5b35-a..:...,...
## TCGA-D8-A1XU-01A-11R-A14M-07 NA:2019-07-31T21:36:15...:19097760-fda9-5a37-a..:...,NA:2019-07-31T21:36:15...:4e239d8d-2557-536f-a..:...,...
## TCGA-BH-A18L-01A-32R-A12D-07 NA:2019-07-31T21:51:11...:46f0e5da-3fd4-5dea-b..:...,NA:2019-07-31T21:51:11...:b02f574e-fdcc-5066-a..:...,...
## TCGA-B6-A0IK-01A-12R-A056-07                                         NA:NA:dc370e08-20f4-5af4-9..:...,NA:NA:f6a2399a-b1f6-5d42-b..:...,...
## TCGA-BH-A18L-11A-42R-A12D-07 NA:2019-07-31T21:51:11...:46f0e5da-3fd4-5dea-b..:...,NA:2019-07-31T21:51:11...:b02f574e-fdcc-5066-a..:...,...
## ...                                                                                                                                    ...
## TCGA-GM-A2DN-01A-11R-A180-07                                         NA:NA:258fc301-61b6-5cf3-9..:...,NA:NA:ebc0bc09-8b4d-5fb7-8..:...,...
## TCGA-B6-A0I1-01A-11R-A21T-07                                         NA:NA:64394cd4-ef14-5d11-8..:...,NA:NA:c47ad272-22b4-5352-b..:...,...
## TCGA-EW-A1IW-01A-11R-A13Q-07                                         NA:NA:32cc83cc-7d68-51b5-9..:...,NA:NA:f4fcec2f-71a4-50c0-9..:...,...
## TCGA-BH-A42V-01A-11R-A24H-07                                         NA:NA:04dd058c-b48a-5a07-8..:...,NA:NA:49948f50-d824-59c6-8..:...,...
## TCGA-AR-A251-01A-12R-A169-07                                         NA:NA:3e2be726-7e39-5bd0-9..:...,NA:NA:7369d53e-5224-5a0b-b..:...,...
##                              last_known_disease_status
##                                            <character>
## TCGA-AC-A8OP-01A-11R-A36F-07              not reported
## TCGA-D8-A1XU-01A-11R-A14M-07              not reported
## TCGA-BH-A18L-01A-32R-A12D-07              not reported
## TCGA-B6-A0IK-01A-12R-A056-07              not reported
## TCGA-BH-A18L-11A-42R-A12D-07              not reported
## ...                                                ...
## TCGA-GM-A2DN-01A-11R-A180-07              not reported
## TCGA-B6-A0I1-01A-11R-A21T-07              not reported
## TCGA-EW-A1IW-01A-11R-A13Q-07              not reported
## TCGA-BH-A42V-01A-11R-A24H-07              not reported
## TCGA-AR-A251-01A-12R-A169-07              not reported
##                              tissue_or_organ_of_origin days_to_last_follow_up
##                                            <character>              <integer>
## TCGA-AC-A8OP-01A-11R-A36F-07               Breast, NOS                    614
## TCGA-D8-A1XU-01A-11R-A14M-07               Breast, NOS                    395
## TCGA-BH-A18L-01A-32R-A12D-07               Breast, NOS                     NA
## TCGA-B6-A0IK-01A-12R-A056-07               Breast, NOS                     NA
## TCGA-BH-A18L-11A-42R-A12D-07               Breast, NOS                     NA
## ...                                                ...                    ...
## TCGA-GM-A2DN-01A-11R-A180-07               Breast, NOS                   3091
## TCGA-B6-A0I1-01A-11R-A21T-07               Breast, NOS                     NA
## TCGA-EW-A1IW-01A-11R-A13Q-07               Breast, NOS                    371
## TCGA-BH-A42V-01A-11R-A24H-07               Breast, NOS                    635
## TCGA-AR-A251-01A-12R-A169-07               Breast, NOS                   3030
##                              age_at_diagnosis      primary_diagnosis
##                                     <integer>            <character>
## TCGA-AC-A8OP-01A-11R-A36F-07            26550 Infiltrating duct ca..
## TCGA-D8-A1XU-01A-11R-A14M-07            20715 Infiltrating duct ca..
## TCGA-BH-A18L-01A-32R-A12D-07            18519 Infiltrating duct ca..
## TCGA-B6-A0IK-01A-12R-A056-07            23150 Infiltrating duct ca..
## TCGA-BH-A18L-11A-42R-A12D-07            18519 Infiltrating duct ca..
## ...                                       ...                    ...
## TCGA-GM-A2DN-01A-11R-A180-07            21288 Infiltrating duct ca..
## TCGA-B6-A0I1-01A-11R-A21T-07            26902 Infiltrating duct ca..
## TCGA-EW-A1IW-01A-11R-A13Q-07            29529 Lobular carcinoma, NOS
## TCGA-BH-A42V-01A-11R-A24H-07            15158 Infiltrating duct ca..
## TCGA-AR-A251-01A-12R-A169-07            18771 Infiltrating duct ca..
##                              prior_malignancy year_of_diagnosis prior_treatment
##                                   <character>         <integer>     <character>
## TCGA-AC-A8OP-01A-11R-A36F-07               no              2011              No
## TCGA-D8-A1XU-01A-11R-A14M-07               no              2010              No
## TCGA-BH-A18L-01A-32R-A12D-07               no              2003              No
## TCGA-B6-A0IK-01A-12R-A056-07               no              1994              No
## TCGA-BH-A18L-11A-42R-A12D-07               no              2003              No
## ...                                       ...               ...             ...
## TCGA-GM-A2DN-01A-11R-A180-07               no              2004              No
## TCGA-B6-A0I1-01A-11R-A21T-07               no              1998              No
## TCGA-EW-A1IW-01A-11R-A13Q-07               no              2010              No
## TCGA-BH-A42V-01A-11R-A24H-07               no              2012              No
## TCGA-AR-A251-01A-12R-A169-07               no              2006              No
##                              ajcc_staging_system_edition ajcc_pathologic_t
##                                              <character>       <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                         7th               T1c
## TCGA-D8-A1XU-01A-11R-A14M-07                         7th               T1c
## TCGA-BH-A18L-01A-32R-A12D-07                         6th                T3
## TCGA-B6-A0IK-01A-12R-A056-07                         4th                T4
## TCGA-BH-A18L-11A-42R-A12D-07                         6th                T3
## ...                                                  ...               ...
## TCGA-GM-A2DN-01A-11R-A180-07                         6th                T2
## TCGA-B6-A0I1-01A-11R-A21T-07                         5th                T2
## TCGA-EW-A1IW-01A-11R-A13Q-07                         7th                T2
## TCGA-BH-A42V-01A-11R-A24H-07                         7th               T1c
## TCGA-AR-A251-01A-12R-A169-07                         6th                T2
##                               morphology ajcc_pathologic_n ajcc_pathologic_m
##                              <character>       <character>       <character>
## TCGA-AC-A8OP-01A-11R-A36F-07      8500/3           N0 (i-)                MX
## TCGA-D8-A1XU-01A-11R-A14M-07      8500/3                N0                M0
## TCGA-BH-A18L-01A-32R-A12D-07      8500/3              N1mi                M0
## TCGA-B6-A0IK-01A-12R-A056-07      8500/3                N1                M0
## TCGA-BH-A18L-11A-42R-A12D-07      8500/3              N1mi                M0
## ...                                  ...               ...               ...
## TCGA-GM-A2DN-01A-11R-A180-07      8500/3           N0 (i-)                M0
## TCGA-B6-A0I1-01A-11R-A21T-07      8500/3                N0                M0
## TCGA-EW-A1IW-01A-11R-A13Q-07      8520/3               N1a                MX
## TCGA-BH-A42V-01A-11R-A24H-07      8500/3              N1mi                M0
## TCGA-AR-A251-01A-12R-A169-07      8500/3                N2                M0
##                              classification_of_tumor           diagnosis_id
##                                          <character>            <character>
## TCGA-AC-A8OP-01A-11R-A36F-07            not reported 7cfc3db9-9892-5ba9-b..
## TCGA-D8-A1XU-01A-11R-A14M-07            not reported 535e0da6-6128-5eb3-9..
## TCGA-BH-A18L-01A-32R-A12D-07            not reported 82eef841-3a64-5fb0-9..
## TCGA-B6-A0IK-01A-12R-A056-07            not reported f12376bc-6e8a-5ac0-a..
## TCGA-BH-A18L-11A-42R-A12D-07            not reported 82eef841-3a64-5fb0-9..
## ...                                              ...                    ...
## TCGA-GM-A2DN-01A-11R-A180-07            not reported 57c58553-d1d1-5e0d-a..
## TCGA-B6-A0I1-01A-11R-A21T-07            not reported b5d064de-7697-57e5-a..
## TCGA-EW-A1IW-01A-11R-A13Q-07            not reported 8cafa6ba-2732-5c82-b..
## TCGA-BH-A42V-01A-11R-A24H-07            not reported 7eb6a599-2ee2-5137-9..
## TCGA-AR-A251-01A-12R-A169-07            not reported 591c58d8-d324-528a-9..
##                              icd_10_code site_of_resection_or_biopsy
##                              <character>                 <character>
## TCGA-AC-A8OP-01A-11R-A36F-07       C50.9                 Breast, NOS
## TCGA-D8-A1XU-01A-11R-A14M-07       C50.9                 Breast, NOS
## TCGA-BH-A18L-01A-32R-A12D-07       C50.9                 Breast, NOS
## TCGA-B6-A0IK-01A-12R-A056-07       C50.9                 Breast, NOS
## TCGA-BH-A18L-11A-42R-A12D-07       C50.9                 Breast, NOS
## ...                                  ...                         ...
## TCGA-GM-A2DN-01A-11R-A180-07       C50.9                 Breast, NOS
## TCGA-B6-A0I1-01A-11R-A21T-07       C50.9                 Breast, NOS
## TCGA-EW-A1IW-01A-11R-A13Q-07       C50.9                 Breast, NOS
## TCGA-BH-A42V-01A-11R-A24H-07       C50.9                 Breast, NOS
## TCGA-AR-A251-01A-12R-A169-07       C50.9                 Breast, NOS
##                               tumor_grade progression_or_recurrence
##                               <character>               <character>
## TCGA-AC-A8OP-01A-11R-A36F-07 not reported              not reported
## TCGA-D8-A1XU-01A-11R-A14M-07 not reported              not reported
## TCGA-BH-A18L-01A-32R-A12D-07 not reported              not reported
## TCGA-B6-A0IK-01A-12R-A056-07 not reported              not reported
## TCGA-BH-A18L-11A-42R-A12D-07 not reported              not reported
## ...                                   ...                       ...
## TCGA-GM-A2DN-01A-11R-A180-07 not reported              not reported
## TCGA-B6-A0I1-01A-11R-A21T-07 not reported              not reported
## TCGA-EW-A1IW-01A-11R-A13Q-07 not reported              not reported
## TCGA-BH-A42V-01A-11R-A24H-07 not reported              not reported
## TCGA-AR-A251-01A-12R-A169-07 not reported              not reported
##                              alcohol_history            exposure_id
##                                  <character>            <character>
## TCGA-AC-A8OP-01A-11R-A36F-07    Not Reported 8829454d-9ece-5771-a..
## TCGA-D8-A1XU-01A-11R-A14M-07    Not Reported 5ec9e2d6-e38e-5a21-9..
## TCGA-BH-A18L-01A-32R-A12D-07    Not Reported 2dbc9779-131f-53fb-9..
## TCGA-B6-A0IK-01A-12R-A056-07    Not Reported d5ef80f3-8ae3-5574-a..
## TCGA-BH-A18L-11A-42R-A12D-07    Not Reported 2dbc9779-131f-53fb-9..
## ...                                      ...                    ...
## TCGA-GM-A2DN-01A-11R-A180-07    Not Reported 5edac88e-386b-5835-b..
## TCGA-B6-A0I1-01A-11R-A21T-07    Not Reported 652b7ab6-af74-57e0-8..
## TCGA-EW-A1IW-01A-11R-A13Q-07    Not Reported 124c6c50-fbbf-5610-a..
## TCGA-BH-A42V-01A-11R-A24H-07    Not Reported 5d950ae3-9d04-570f-b..
## TCGA-AR-A251-01A-12R-A169-07    Not Reported 017eef2e-5c00-5dd7-b..
##                                                race      gender
##                                         <character> <character>
## TCGA-AC-A8OP-01A-11R-A36F-07 black or african ame..      female
## TCGA-D8-A1XU-01A-11R-A14M-07                  white      female
## TCGA-BH-A18L-01A-32R-A12D-07                  white      female
## TCGA-B6-A0IK-01A-12R-A056-07                  white      female
## TCGA-BH-A18L-11A-42R-A12D-07                  white      female
## ...                                             ...         ...
## TCGA-GM-A2DN-01A-11R-A180-07                  white      female
## TCGA-B6-A0I1-01A-11R-A21T-07 black or african ame..      female
## TCGA-EW-A1IW-01A-11R-A13Q-07                  white      female
## TCGA-BH-A42V-01A-11R-A24H-07 black or african ame..      female
## TCGA-AR-A251-01A-12R-A169-07                  white      female
##                                           ethnicity vital_status age_at_index
##                                         <character>  <character>    <integer>
## TCGA-AC-A8OP-01A-11R-A36F-07 not hispanic or latino        Alive           72
## TCGA-D8-A1XU-01A-11R-A14M-07 not hispanic or latino        Alive           56
## TCGA-BH-A18L-01A-32R-A12D-07 not hispanic or latino         Dead           50
## TCGA-B6-A0IK-01A-12R-A056-07 not hispanic or latino         Dead           63
## TCGA-BH-A18L-11A-42R-A12D-07 not hispanic or latino         Dead           50
## ...                                             ...          ...          ...
## TCGA-GM-A2DN-01A-11R-A180-07 not hispanic or latino        Alive           58
## TCGA-B6-A0I1-01A-11R-A21T-07 not hispanic or latino         Dead           73
## TCGA-EW-A1IW-01A-11R-A13Q-07     hispanic or latino        Alive           80
## TCGA-BH-A42V-01A-11R-A24H-07 not hispanic or latino        Alive           41
## TCGA-AR-A251-01A-12R-A169-07 not hispanic or latino        Alive           51
##                              days_to_birth year_of_birth         demographic_id
##                                  <integer>     <integer>            <character>
## TCGA-AC-A8OP-01A-11R-A36F-07        -26550          1939 e93f6e03-699d-51ff-8..
## TCGA-D8-A1XU-01A-11R-A14M-07        -20715          1954 2d23b71d-cca7-5db8-a..
## TCGA-BH-A18L-01A-32R-A12D-07        -18519          1953 6ee59a01-876d-5f82-9..
## TCGA-B6-A0IK-01A-12R-A056-07        -23150          1931 9a206ac3-aff1-50c1-b..
## TCGA-BH-A18L-11A-42R-A12D-07        -18519          1953 6ee59a01-876d-5f82-9..
## ...                                    ...           ...                    ...
## TCGA-GM-A2DN-01A-11R-A180-07        -21288          1946 91a06fb2-6332-5200-9..
## TCGA-B6-A0I1-01A-11R-A21T-07        -26902          1925 9f1d5c3f-f926-5594-a..
## TCGA-EW-A1IW-01A-11R-A13Q-07        -29529          1930 2cd46bce-530d-593c-9..
## TCGA-BH-A42V-01A-11R-A24H-07        -15158          1971 f3b2aad9-add6-5d43-a..
## TCGA-AR-A251-01A-12R-A169-07        -18771          1955 9b1440be-8ed5-5c5e-b..
##                              year_of_death days_to_death bcr_patient_barcode
##                                  <integer>     <integer>         <character>
## TCGA-AC-A8OP-01A-11R-A36F-07            NA            NA    TCGA-AC-A8OP-01A
## TCGA-D8-A1XU-01A-11R-A14M-07            NA            NA    TCGA-D8-A1XU-01A
## TCGA-BH-A18L-01A-32R-A12D-07          2005           811    TCGA-BH-A18L-01A
## TCGA-B6-A0IK-01A-12R-A056-07          1995           571    TCGA-B6-A0IK-01A
## TCGA-BH-A18L-11A-42R-A12D-07          2005           811    TCGA-BH-A18L-11A
## ...                                    ...           ...                 ...
## TCGA-GM-A2DN-01A-11R-A180-07            NA            NA    TCGA-GM-A2DN-01A
## TCGA-B6-A0I1-01A-11R-A21T-07          2004          2361    TCGA-B6-A0I1-01A
## TCGA-EW-A1IW-01A-11R-A13Q-07            NA            NA    TCGA-EW-A1IW-01A
## TCGA-BH-A42V-01A-11R-A24H-07            NA            NA    TCGA-BH-A42V-01A
## TCGA-AR-A251-01A-12R-A169-07            NA            NA    TCGA-AR-A251-01A
##                              primary_site  project_id
##                                    <list> <character>
## TCGA-AC-A8OP-01A-11R-A36F-07       Breast   TCGA-BRCA
## TCGA-D8-A1XU-01A-11R-A14M-07       Breast   TCGA-BRCA
## TCGA-BH-A18L-01A-32R-A12D-07       Breast   TCGA-BRCA
## TCGA-B6-A0IK-01A-12R-A056-07       Breast   TCGA-BRCA
## TCGA-BH-A18L-11A-42R-A12D-07       Breast   TCGA-BRCA
## ...                                   ...         ...
## TCGA-GM-A2DN-01A-11R-A180-07       Breast   TCGA-BRCA
## TCGA-B6-A0I1-01A-11R-A21T-07       Breast   TCGA-BRCA
## TCGA-EW-A1IW-01A-11R-A13Q-07       Breast   TCGA-BRCA
## TCGA-BH-A42V-01A-11R-A24H-07       Breast   TCGA-BRCA
## TCGA-AR-A251-01A-12R-A169-07       Breast   TCGA-BRCA
##                                                                                        disease_type
##                                                                                              <list>
## TCGA-AC-A8OP-01A-11R-A36F-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-D8-A1XU-01A-11R-A14M-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-BH-A18L-01A-32R-A12D-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-B6-A0IK-01A-12R-A056-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-BH-A18L-11A-42R-A12D-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## ...                                                                                             ...
## TCGA-GM-A2DN-01A-11R-A180-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-B6-A0I1-01A-11R-A21T-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-EW-A1IW-01A-11R-A13Q-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-BH-A42V-01A-11R-A24H-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
## TCGA-AR-A251-01A-12R-A169-07 Basal Cell Neoplasms,Complex Epithelial N..,Adenomas and Adenoca..,...
##                                                name releasable  released
##                                         <character>  <logical> <logical>
## TCGA-AC-A8OP-01A-11R-A36F-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-D8-A1XU-01A-11R-A14M-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-BH-A18L-01A-32R-A12D-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-B6-A0IK-01A-12R-A056-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-BH-A18L-11A-42R-A12D-07 Breast Invasive Carc..       TRUE      TRUE
## ...                                             ...        ...       ...
## TCGA-GM-A2DN-01A-11R-A180-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-B6-A0I1-01A-11R-A21T-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-EW-A1IW-01A-11R-A13Q-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-BH-A42V-01A-11R-A24H-07 Breast Invasive Carc..       TRUE      TRUE
## TCGA-AR-A251-01A-12R-A169-07 Breast Invasive Carc..       TRUE      TRUE
##                              preservation_method days_to_sample_procurement
##                                      <character>                  <integer>
## TCGA-AC-A8OP-01A-11R-A36F-07                  NA                         NA
## TCGA-D8-A1XU-01A-11R-A14M-07                  NA                         NA
## TCGA-BH-A18L-01A-32R-A12D-07                  NA                         NA
## TCGA-B6-A0IK-01A-12R-A056-07                  NA                         NA
## TCGA-BH-A18L-11A-42R-A12D-07                  NA                         NA
## ...                                          ...                        ...
## TCGA-GM-A2DN-01A-11R-A180-07                  NA                         NA
## TCGA-B6-A0I1-01A-11R-A21T-07                  NA                         NA
## TCGA-EW-A1IW-01A-11R-A13Q-07                  NA                         NA
## TCGA-BH-A42V-01A-11R-A24H-07                  NA                         NA
## TCGA-AR-A251-01A-12R-A169-07                  NA                         NA
##                              paper_patient paper_Tumor.Type
##                                <character>      <character>
## TCGA-AC-A8OP-01A-11R-A36F-07  TCGA-AC-A8OP             BRCA
## TCGA-D8-A1XU-01A-11R-A14M-07  TCGA-D8-A1XU             BRCA
## TCGA-BH-A18L-01A-32R-A12D-07  TCGA-BH-A18L             BRCA
## TCGA-B6-A0IK-01A-12R-A056-07  TCGA-B6-A0IK             BRCA
## TCGA-BH-A18L-11A-42R-A12D-07            NA               NA
## ...                                    ...              ...
## TCGA-GM-A2DN-01A-11R-A180-07  TCGA-GM-A2DN             BRCA
## TCGA-B6-A0I1-01A-11R-A21T-07  TCGA-B6-A0I1             BRCA
## TCGA-EW-A1IW-01A-11R-A13Q-07  TCGA-EW-A1IW             BRCA
## TCGA-BH-A42V-01A-11R-A24H-07  TCGA-BH-A42V             BRCA
## TCGA-AR-A251-01A-12R-A169-07  TCGA-AR-A251             BRCA
##                              paper_Included_in_previous_marker_papers
##                                                           <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                                       NO
## TCGA-D8-A1XU-01A-11R-A14M-07                                      YES
## TCGA-BH-A18L-01A-32R-A12D-07                                      YES
## TCGA-B6-A0IK-01A-12R-A056-07                                      YES
## TCGA-BH-A18L-11A-42R-A12D-07                                       NA
## ...                                                               ...
## TCGA-GM-A2DN-01A-11R-A180-07                                      YES
## TCGA-B6-A0I1-01A-11R-A21T-07                                       NO
## TCGA-EW-A1IW-01A-11R-A13Q-07                                      YES
## TCGA-BH-A42V-01A-11R-A24H-07                                       NO
## TCGA-AR-A251-01A-12R-A169-07                                      YES
##                              paper_vital_status paper_days_to_birth
##                                     <character>         <character>
## TCGA-AC-A8OP-01A-11R-A36F-07              Alive              -26550
## TCGA-D8-A1XU-01A-11R-A14M-07              Alive              -20715
## TCGA-BH-A18L-01A-32R-A12D-07               Dead              -18519
## TCGA-B6-A0IK-01A-12R-A056-07               Dead              -23150
## TCGA-BH-A18L-11A-42R-A12D-07                 NA                  NA
## ...                                         ...                 ...
## TCGA-GM-A2DN-01A-11R-A180-07              Alive              -21288
## TCGA-B6-A0I1-01A-11R-A21T-07               Dead              -26902
## TCGA-EW-A1IW-01A-11R-A13Q-07              Alive              -29529
## TCGA-BH-A42V-01A-11R-A24H-07              Alive              -15158
## TCGA-AR-A251-01A-12R-A169-07              Alive              -18771
##                              paper_days_to_death paper_days_to_last_followup
##                                      <character>                 <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                  NA                         614
## TCGA-D8-A1XU-01A-11R-A14M-07                  NA                         395
## TCGA-BH-A18L-01A-32R-A12D-07                 811                          NA
## TCGA-B6-A0IK-01A-12R-A056-07                 571                          NA
## TCGA-BH-A18L-11A-42R-A12D-07                  NA                          NA
## ...                                          ...                         ...
## TCGA-GM-A2DN-01A-11R-A180-07                  NA                        3091
## TCGA-B6-A0I1-01A-11R-A21T-07                2361                          NA
## TCGA-EW-A1IW-01A-11R-A13Q-07                  NA                         371
## TCGA-BH-A42V-01A-11R-A24H-07                  NA                         635
## TCGA-AR-A251-01A-12R-A169-07                  NA                        3030
##                              paper_age_at_initial_pathologic_diagnosis
##                                                              <numeric>
## TCGA-AC-A8OP-01A-11R-A36F-07                                        72
## TCGA-D8-A1XU-01A-11R-A14M-07                                        56
## TCGA-BH-A18L-01A-32R-A12D-07                                        50
## TCGA-B6-A0IK-01A-12R-A056-07                                        63
## TCGA-BH-A18L-11A-42R-A12D-07                                        NA
## ...                                                                ...
## TCGA-GM-A2DN-01A-11R-A180-07                                        58
## TCGA-B6-A0I1-01A-11R-A21T-07                                        73
## TCGA-EW-A1IW-01A-11R-A13Q-07                                        80
## TCGA-BH-A42V-01A-11R-A24H-07                                        41
## TCGA-AR-A251-01A-12R-A169-07                                        51
##                              paper_pathologic_stage paper_Tumor_Grade
##                                         <character>       <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                Stage_I                NA
## TCGA-D8-A1XU-01A-11R-A14M-07                Stage_I                NA
## TCGA-BH-A18L-01A-32R-A12D-07              Stage_III                NA
## TCGA-B6-A0IK-01A-12R-A056-07              Stage_III                NA
## TCGA-BH-A18L-11A-42R-A12D-07                     NA                NA
## ...                                             ...               ...
## TCGA-GM-A2DN-01A-11R-A180-07               Stage_II                NA
## TCGA-B6-A0I1-01A-11R-A21T-07               Stage_II                NA
## TCGA-EW-A1IW-01A-11R-A13Q-07               Stage_II                NA
## TCGA-BH-A42V-01A-11R-A24H-07                Stage_I                NA
## TCGA-AR-A251-01A-12R-A169-07              Stage_III                NA
##                              paper_BRCA_Pathology paper_BRCA_Subtype_PAM50
##                                       <character>              <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                   NA                     LumA
## TCGA-D8-A1XU-01A-11R-A14M-07                  IDC                     LumA
## TCGA-BH-A18L-01A-32R-A12D-07                Mixed                     LumB
## TCGA-B6-A0IK-01A-12R-A056-07                  IDC                     Her2
## TCGA-BH-A18L-11A-42R-A12D-07                   NA                       NA
## ...                                           ...                      ...
## TCGA-GM-A2DN-01A-11R-A180-07                  IDC                     LumA
## TCGA-B6-A0I1-01A-11R-A21T-07                   NA                    Basal
## TCGA-EW-A1IW-01A-11R-A13Q-07                  ILC                     LumA
## TCGA-BH-A42V-01A-11R-A24H-07                  IDC                     LumA
## TCGA-AR-A251-01A-12R-A169-07                   NA                    Basal
##                              paper_MSI_status paper_HPV_Status
##                                   <character>      <character>
## TCGA-AC-A8OP-01A-11R-A36F-07               NA               NA
## TCGA-D8-A1XU-01A-11R-A14M-07               NA               NA
## TCGA-BH-A18L-01A-32R-A12D-07               NA               NA
## TCGA-B6-A0IK-01A-12R-A056-07               NA               NA
## TCGA-BH-A18L-11A-42R-A12D-07               NA               NA
## ...                                       ...              ...
## TCGA-GM-A2DN-01A-11R-A180-07               NA               NA
## TCGA-B6-A0I1-01A-11R-A21T-07               NA               NA
## TCGA-EW-A1IW-01A-11R-A13Q-07               NA               NA
## TCGA-BH-A42V-01A-11R-A24H-07               NA               NA
## TCGA-AR-A251-01A-12R-A169-07               NA               NA
##                              paper_tobacco_smoking_history paper_CNV Clusters
##                                                <character>        <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                            NA                 C4
## TCGA-D8-A1XU-01A-11R-A14M-07                            NA                 C4
## TCGA-BH-A18L-01A-32R-A12D-07                            NA                 C6
## TCGA-B6-A0IK-01A-12R-A056-07                            NA                 C3
## TCGA-BH-A18L-11A-42R-A12D-07                            NA                 NA
## ...                                                    ...                ...
## TCGA-GM-A2DN-01A-11R-A180-07                            NA                 C1
## TCGA-B6-A0I1-01A-11R-A21T-07                            NA                 C4
## TCGA-EW-A1IW-01A-11R-A13Q-07                            NA                 C1
## TCGA-BH-A42V-01A-11R-A24H-07                            NA                 C5
## TCGA-AR-A251-01A-12R-A169-07                            NA                 C4
##                              paper_Mutation Clusters
##                                          <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                      C1
## TCGA-D8-A1XU-01A-11R-A14M-07                      C4
## TCGA-BH-A18L-01A-32R-A12D-07                      C1
## TCGA-B6-A0IK-01A-12R-A056-07                      C9
## TCGA-BH-A18L-11A-42R-A12D-07                      NA
## ...                                              ...
## TCGA-GM-A2DN-01A-11R-A180-07                      C4
## TCGA-B6-A0I1-01A-11R-A21T-07                      C4
## TCGA-EW-A1IW-01A-11R-A13Q-07                      C3
## TCGA-BH-A42V-01A-11R-A24H-07                      C4
## TCGA-AR-A251-01A-12R-A169-07                      C1
##                              paper_DNA.Methylation Clusters paper_mRNA Clusters
##                                                 <character>         <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                             C2                  C1
## TCGA-D8-A1XU-01A-11R-A14M-07                             C2                  C1
## TCGA-BH-A18L-01A-32R-A12D-07                             C1                  C1
## TCGA-B6-A0IK-01A-12R-A056-07                             C2                  C2
## TCGA-BH-A18L-11A-42R-A12D-07                             NA                  NA
## ...                                                     ...                 ...
## TCGA-GM-A2DN-01A-11R-A180-07                             C1                  C2
## TCGA-B6-A0I1-01A-11R-A21T-07                             C4                  C4
## TCGA-EW-A1IW-01A-11R-A13Q-07                             C2                  C1
## TCGA-BH-A42V-01A-11R-A24H-07                             C1                  C2
## TCGA-AR-A251-01A-12R-A169-07                             C4                  C4
##                              paper_miRNA Clusters paper_lncRNA Clusters
##                                       <character>           <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                   C2                    NA
## TCGA-D8-A1XU-01A-11R-A14M-07                   C3                    C6
## TCGA-BH-A18L-01A-32R-A12D-07                   C3                    C2
## TCGA-B6-A0IK-01A-12R-A056-07                   C3                    C1
## TCGA-BH-A18L-11A-42R-A12D-07                   NA                    NA
## ...                                           ...                   ...
## TCGA-GM-A2DN-01A-11R-A180-07                   C2                    C2
## TCGA-B6-A0I1-01A-11R-A21T-07                   C7                    NA
## TCGA-EW-A1IW-01A-11R-A13Q-07                   C3                    NA
## TCGA-BH-A42V-01A-11R-A24H-07                   C2                    NA
## TCGA-AR-A251-01A-12R-A169-07                   C7                    C6
##                              paper_Protein Clusters paper_PARADIGM Clusters
##                                         <character>             <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                     NA                      C6
## TCGA-D8-A1XU-01A-11R-A14M-07                     C1                      C5
## TCGA-BH-A18L-01A-32R-A12D-07                     C1                      C5
## TCGA-B6-A0IK-01A-12R-A056-07                     C2                      C4
## TCGA-BH-A18L-11A-42R-A12D-07                     NA                      NA
## ...                                             ...                     ...
## TCGA-GM-A2DN-01A-11R-A180-07                     C1                      C6
## TCGA-B6-A0I1-01A-11R-A21T-07                     NA                      C2
## TCGA-EW-A1IW-01A-11R-A13Q-07                     NA                      C6
## TCGA-BH-A42V-01A-11R-A24H-07                     C1                      C6
## TCGA-AR-A251-01A-12R-A169-07                     C2                      C2
##                              paper_Pan-Gyn Clusters
##                                         <character>
## TCGA-AC-A8OP-01A-11R-A36F-07                     NA
## TCGA-D8-A1XU-01A-11R-A14M-07                     C1
## TCGA-BH-A18L-01A-32R-A12D-07                     C5
## TCGA-B6-A0IK-01A-12R-A056-07                     C5
## TCGA-BH-A18L-11A-42R-A12D-07                     NA
## ...                                             ...
## TCGA-GM-A2DN-01A-11R-A180-07                     C1
## TCGA-B6-A0I1-01A-11R-A21T-07                     NA
## TCGA-EW-A1IW-01A-11R-A13Q-07                     NA
## TCGA-BH-A42V-01A-11R-A24H-07                     C1
## TCGA-AR-A251-01A-12R-A169-07                     C4
```

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
##  [1] SummarizedExperiment_1.24.0 Biobase_2.54.0             
##  [3] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
##  [5] IRanges_2.28.0              S4Vectors_0.32.4           
##  [7] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
##  [9] matrixStats_0.62.0          DT_0.23                    
## [11] dplyr_1.0.9                 TCGAbiolinks_2.25.2        
## [13] TCGAbiolinksGUI.data_1.15.1
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7           bit64_4.0.5            filelock_1.0.2        
##  [4] progress_1.2.2         httr_1.4.3             tools_4.1.2           
##  [7] bslib_0.3.1            utf8_1.2.2             R6_2.5.1              
## [10] DBI_1.1.3              colorspace_2.0-3       tidyselect_1.1.2      
## [13] prettyunits_1.1.1      bit_4.0.4              curl_4.3.2            
## [16] compiler_4.1.2         cli_3.3.0              rvest_1.0.2           
## [19] xml2_1.3.3             DelayedArray_0.20.0    sass_0.4.1            
## [22] scales_1.2.0           readr_2.1.2            rappdirs_0.3.3        
## [25] stringr_1.4.0          digest_0.6.29          rmarkdown_2.14        
## [28] XVector_0.34.0         pkgconfig_2.0.3        htmltools_0.5.2       
## [31] dbplyr_2.2.1           fastmap_1.1.0          htmlwidgets_1.5.4     
## [34] rlang_1.0.4            rstudioapi_0.13        RSQLite_2.2.14        
## [37] jquerylib_0.1.4        generics_0.1.2         jsonlite_1.8.0        
## [40] RCurl_1.98-1.7         magrittr_2.0.3.9000    GenomeInfoDbData_1.2.7
## [43] Matrix_1.4-1           Rcpp_1.0.8.3           munsell_0.5.0         
## [46] fansi_1.0.3            lifecycle_1.0.1        stringi_1.7.6         
## [49] yaml_2.3.5             zlibbioc_1.40.0        plyr_1.8.6            
## [52] BiocFileCache_2.2.1    grid_4.1.2             blob_1.2.3            
## [55] crayon_1.5.1           lattice_0.20-45        Biostrings_2.62.0     
## [58] hms_1.1.1              KEGGREST_1.34.0        knitr_1.38            
## [61] pillar_1.7.0           biomaRt_2.50.3         XML_3.99-0.9          
## [64] glue_1.6.2             evaluate_0.15          downloader_0.4        
## [67] data.table_1.14.2      png_0.1-7              vctrs_0.4.1           
## [70] tzdb_0.3.0             gtable_0.3.0           purrr_0.3.4           
## [73] tidyr_1.2.0            assertthat_0.2.1       cachem_1.0.6          
## [76] ggplot2_3.3.6          xfun_0.31              tibble_3.1.7          
## [79] AnnotationDbi_1.56.2   memoise_2.0.1          ellipsis_0.3.2
```
