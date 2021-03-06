---
title: "Analysing sex"
author: Emeric Texeraud
date: '2022-03-03'
categories: 
tags:
  - rmarkdown
subtitle: '' 
summary: 'An Rmarkdown vignette showing how to anaylse sex with single-cell data'
authors: []
lastmod: '2020-05-28T20:45:57-07:00'
featured: no
disable_jquery: no
image:
  caption: 'UMAP.png'
  focal_point: ''
  preview_only: no
projects: []
---


```r
# Load packages
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(gridExtra)
})

seur <- readRDS("/home/william/Téléchargements/p12_new.rds")
```

# Sex analysis

To assign cells as either male or female, we first need to find the male genes. We found the list of Y genes on [this website.]() We keep the Y genes that are in our dataset.

```r
Y = read.table("/home/william/Téléchargements/YGenes.csv", header = TRUE)
Ygenes = as.character(Y$Ygenes)
Ygenes = colnames(FetchData(seur, vars = Ygenes))
```

We then find the mean expression of all Y genes across all cells, and the expression of Xist, marker of female cells. We assign as male every cell that an expression of any of the male gene, female the cells that do not express male genes, but express Xist, and unknown the cells that express none of those genes.

```r
DefaultAssay(seur) <- "SCT"
male_expression <- rowMeans(FetchData(seur, vars = Ygenes, slot = "data"))
female_expression <- rowMeans(FetchData(seur, vars = "Xist", slot = "data"))
seur@meta.data$male_expression <- male_expression
seur@meta.data$female_expression <- female_expression
seur@meta.data$sex <- ifelse(male_expression > 0, "male",
                             ifelse(female_expression > 0, "female", "unknown"))

FeaturePlot(seur, features = c("male_expression", "female_expression"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" width="768" />

```r
DimPlot(seur, group.by = "sex")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-3-1.png" width="672" />
