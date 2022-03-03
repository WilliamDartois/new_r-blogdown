---
title: "Platform Rpackage tutorial"
author: Emeric Texersaud
date: '2022-03-03'
categories: 
tags:
  - rmarkdown
subtitle: '' 
summary: 'An Rmarkdown vignette showing how to use our own Rpackage.'
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
suppressMessages({
  library(LxCxUtils)
  library(Seurat)
  library(SeuratData)
})
```

The `LxCxUtils` package mainly uses Seurat objects to get the expression data for the plots, so we first need to setup a Seurat object. We will use the pbmc 3k dataset from the `SeuratData` package.

```r
suppressMessages({
  data("pbmc3k")
  pbmc3k
})
#> An object of class Seurat 
#> 13714 features across 2700 samples within 1 assay 
#> Active assay: RNA (13714 features, 0 variable features)
```




```r
pbmc3k <- NormalizeData(pbmc3k, verbose = FALSE)
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pbmc3k <- ScaleData(pbmc3k, verbose = FALSE)
pbmc3k <- RunPCA(pbmc3k, verbose = FALSE)
```


The first visualisation function from `LxCxUtils` is the `PlotGeneRank`. This function will rank the genes by their expression in each identity of the dataset, and we can even highlight some genes of interest.

```r
#highlight no genes
plotGenesRank(pbmc3k, genes_highlight = c())
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
#highlight some genes
plotGenesRank(pbmc3k, genes_highlight = c("FTL", "LYZ", "MALAT1"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-5-2.png" width="672" />

And we can separate the plots by any cluster column from the metadata.

```r
#highlight some genes and separate by clusters
plotGenesRank(pbmc3k, genes_highlight = c("FTL", "LYZ", "MALAT1"), group.by = "seurat_annotations")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="768" />

Another visualization function from the package, which can be used for the QC, is the `QC_barplot`. It will plot the distribution of cells along the QC column from the metadata, as an histogram.

```r
#Grouping by default : orig.ident
QC_barplot(pbmc3k, features = c("nFeature_RNA"))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="672" />

We can color the histogram by another column from the metadata

```r
#Coloring by clusters
QC_barplot(pbmc3k, features = c("nFeature_RNA"), group.by = "seurat_annotations")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" width="672" />
And even separate the histograms, with one for each cluster of this column.

```r
#Coloring and grouping by clusters
QC_barplot(pbmc3k, features = c("nFeature_RNA"), group.by = "seurat_annotations", wrap = TRUE)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-9-1.png" width="960" />

Another feature of the package is to compute and plot the PCA variance explained. For this we use the function `plotVarianceExplained` or `plotCumulativeVarExplained`, which respectively plot the variance for each PC, and plot the cumulative variance for each PC, that is adding the variance of all previous PC to each PC.

```r
plotVarianceExplained(pbmc3k)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" width="672" />

```r
plotCumulativeVarExplained(pbmc3k)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-2.png" width="672" />

One of the main plotting function from this package is the `seuratBoxPlot`, which aims to complement the `VlnPlot` from Seurat with a boxplot representation.

```r
#Box plot
p1 <- seuratBoxPlot(pbmc3k, "FTL", group.by = "seurat_annotations", pt.size = 0.1) + ggtitle("Boxplot of FTL")
# seuratBoxPlot(pbmc3k, "FTL", group.by = "seurat_annotations", pt.size = 0.1)
#VlnPlot from Seurat
p2 <- VlnPlot(pbmc3k, "FTL", group.by = "seurat_annotations", pt.size = 0.1) + ggtitle("Violin plot of FTL")
gridExtra::grid.arrange(p1,p2,ncol = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-1.png" width="1152" />

To choose the best dimension for UMAP, the package has a function that computes all the dimension you ask, and give a plot with the umap for each of these dimensions. This function is called `multiUMAP`.  

```r
multiUMAP(pbmc3k, dims = 10:20, out = TRUE, save = FALSE, verbose = FALSE)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-12-1.png" width="1152" />



