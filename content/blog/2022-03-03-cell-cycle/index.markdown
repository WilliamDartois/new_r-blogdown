---
title: "Analysing cell cycle"
author: Emeric Texeraud
date: '2022-03-03'
categories: 
tags:
  - rmarkdown
subtitle: '' 
summary: 'An Rmarkdown vignette showing how to anaylse cell cycle phases with single-cell data'
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

# Cell cycle analysis

To find our cycling cells, we first need to get the genes involved in the G2M and S phases of the cell cycle. For that, we will get the human cycling genes from seurat, and transform them to mouse format (uppercase first letter).

```r
s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
```

We then use the `CellcycleScoring` function from Seurat to score the S and G2M phases, and assign them to cells. What it does is scoring the S phase using S genes, then scoring the G2M phase using G2M genes. Finally, it will assign the cells to S phase if they have a high score for S genes, or G2M for a high G2M score.

```r
seur <- CellCycleScoring(seur, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
```

We can then see which cells have a high S score or high G2M score.

```r
FeaturePlot(seur, features = c("S.Score", "G2M.Score"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-1-1.png" width="768" />
We can also see which cells are assigned S or G2M score.

```r
DimPlot(seur, group.by = "Phase")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-2-1.png" width="672" />
