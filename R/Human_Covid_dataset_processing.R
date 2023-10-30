library(Seurat)
library(ggplot2)
source("~/Dropbox/RProjects/scrnaseq/plotting_functions.R")

load("Covid_vs_healthy_fibroblasts_Teichman.RData")
table(covid.fibroblasts$Study)
table(covid.fibroblasts$Disease)

DimPlot(covid.fibroblasts, group.by = "Study", label = TRUE, split.by = "Disease")

################################################################
### Subset the Teichman healthy samples and perform analysis ###
################################################################

teichman.fibroblasts <- subset(covid.fibroblasts, Study == "Teichmann")
dim(teichman.fibroblasts)
table(teichman.fibroblasts$Sample)

## Keep samples with > 100 fibroblasts
samples.keep <- names(table(teichman.fibroblasts$Sample))[which(table(teichman.fibroblasts$Sample) > 100)]
samples.keep

teichman.fibroblasts <- subset(teichman.fibroblasts, Sample %in% samples.keep)
dim(teichman.fibroblasts)

dataset.list <- Seurat::SplitObject(teichman.fibroblasts, split.by = "Sample")

## Run log-normalization
dataset.list <- lapply(X = dataset.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
  return(x)
})

int.features <- SelectIntegrationFeatures(object.list = dataset.list, nfeatures = 2000)

teichman.anchors <- FindIntegrationAnchors(object.list = dataset.list, 
                                           anchor.features = int.features,
                                           normalization.method = "LogNormalize")
teichman.combined.seurat <- IntegrateData(anchorset = teichman.anchors)

remove(teichman.anchors)
remove(dataset.list)
remove(teichman.fibroblasts)

DefaultAssay(teichman.combined.seurat) <- "integrated"

teichman.combined.seurat <- ScaleData(teichman.combined.seurat)
teichman.combined.seurat <- RunPCA(teichman.combined.seurat, verbose = FALSE)
teichman.combined.seurat <- RunUMAP(teichman.combined.seurat, reduction = "pca", dims = 1:40)

teichman.combined.seurat <- FindNeighbors(teichman.combined.seurat, reduction = "pca", dims = 1:40)
teichman.combined.seurat <- FindClusters(teichman.combined.seurat, resolution = seq(0.2, 0.6, 0.1))

Idents(teichman.combined.seurat) <- teichman.combined.seurat$integrated_snn_res.0.4

DimPlot(teichman.combined.seurat, label = TRUE)
DimPlot(teichman.combined.seurat, label = TRUE, split.by = "Sample") + facet_wrap(~Sample)

DefaultAssay(teichman.combined.seurat) <- "RNA"
teichman.combined.seurat <- NormalizeData(teichman.combined.seurat)

### Check cluster 6
res.table <- FindMarkers(teichman.combined.seurat, ident.1 = "9", only.pos = TRUE, logfc.threshold = 0.5)

## Glial cell markers
PlotExpressionUMAP(teichman.combined.seurat, c("CHL1", "PTPRZ1", "CADM2", "ANK3"))

## Immune cells
PlotExpressionUMAP(teichman.combined.seurat, c("PTPRC", "CD79A", "CD3D", "CSF3R"))

## Mural cells
PlotExpressionUMAP(teichman.combined.seurat, c("RGS5", "PDGFRB", "VTN", "KCNJ8"))

## EC markers
PlotExpressionUMAP(teichman.combined.seurat, c("KDR", "PECAM1", "CDH5", "EMCN"))

## 
PlotExpressionUMAP(teichman.combined.seurat, c("SCARA5", "CD248", "CD55", "GFPT2"))
PlotExpressionUMAP(teichman.combined.seurat, c("POSTN", "MEOX1", "CILP", "COL8A1"))

DimPlot(teichman.combined.seurat, label = TRUE)

### Subset high-confidence fibroblasts
teichman.combined.seurat <- subset(teichman.combined.seurat, ident = c("0", "1", "2", "3", "5", "6", "10"))
DimPlot(teichman.combined.seurat, label = TRUE)

teichman.combined.seurat <- subset( teichman.combined.seurat, UMAP_1 < 5 & UMAP_2 < 5 )
DimPlot(teichman.combined.seurat, label = TRUE)

### Save the fibroblast subset Seurat object
save(teichman.combined.seurat, file = "Teichman_healthy_fibroblasts_integrated_PC40.RData")


###############################################################
### Subset the Ellinor healthy samples and perform analysis ###
###############################################################

ellinor.fibroblasts <- subset(covid.fibroblasts, Study == "Ellinor")
dim(ellinor.fibroblasts)
table(ellinor.fibroblasts$Sample)

## Keep samples with > 100 fibroblasts
samples.keep <- names(table(ellinor.fibroblasts$Sample))[which(table(ellinor.fibroblasts$Sample) > 100)]
samples.keep

ellinor.fibroblasts <- subset(ellinor.fibroblasts, Sample %in% samples.keep)
dim(ellinor.fibroblasts)

dataset.list <- Seurat::SplitObject(ellinor.fibroblasts, split.by = "Sample")

## Run log-normalization
dataset.list <- lapply(X = dataset.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
  return(x)
})

ellinor.features <- SelectIntegrationFeatures(object.list = dataset.list, nfeatures = 2000)

ellinor.anchors <- FindIntegrationAnchors(object.list = dataset.list, 
                                          anchor.features = int.features,
                                          normalization.method = "LogNormalize")
ellinor.combined.seurat <- IntegrateData(anchorset = ellinor.anchors)

remove(ellinor.anchors)
remove(dataset.list)
remove(ellinor.fibroblasts)

DefaultAssay(ellinor.combined.seurat) <- "integrated"

ellinor.combined.seurat <- ScaleData(ellinor.combined.seurat)
ellinor.combined.seurat <- RunPCA(ellinor.combined.seurat, verbose = FALSE)
ellinor.combined.seurat <- RunUMAP(ellinor.combined.seurat, reduction = "pca", dims = 1:40)

ellinor.combined.seurat <- FindNeighbors(ellinor.combined.seurat, reduction = "pca", dims = 1:40)
ellinor.combined.seurat <- FindClusters(ellinor.combined.seurat, resolution = seq(0.2, 0.6, 0.1))

Idents(ellinor.combined.seurat) <- ellinor.combined.seurat$integrated_snn_res.0.4

DimPlot(ellinor.combined.seurat, label = TRUE)
DimPlot(ellinor.combined.seurat, label = TRUE, split.by = "Sample") + facet_wrap(~Sample)

DefaultAssay(ellinor.combined.seurat) <- "RNA"
ellinor.combined.seurat <- NormalizeData(ellinor.combined.seurat)

### Check cluster 6
res.table <- FindMarkers(ellinor.combined.seurat, ident.1 = "9", only.pos = TRUE, logfc.threshold = 0.5)

### Appears to correspond to Neural cells 
VlnPlot(ellinor.combined.seurat, c("CHL1", "PTPRZ1", "PLP1", "SOX10", "CADM2", "ANK3"))
PlotExpressionUMAP(ellinor.combined.seurat, c("CHL1", "PTPRZ1", "PLP1", "SOX10", "CADM2", "ANK3"))
PlotExpressionUMAP(ellinor.combined.seurat, c("PTPRC", "KDR", "TNNT2", "RGS5"))

### Subset high-confidence fibroblasts
ellinor.combined.seurat <- subset(ellinor.combined.seurat, ident = c("0", "1", "6", "7"))

DimPlot(ellinor.combined.seurat, label = TRUE)
ellinor.combined.seurat <- subset( ellinor.combined.seurat, UMAP_1 < 4 & UMAP_2 > -3 )
DimPlot(ellinor.combined.seurat, label = TRUE)

### Save the fibroblast subset Seurat object
save(ellinor.combined.seurat, file = "Ellinor_healthy_fibroblasts_integrated_PC40.RData")


###############################################
### Subset the Broad COVID-infected samples ###
###############################################

broad.fibroblasts <- subset(covid.fibroblasts, Study == "broad")
dim(broad.fibroblasts)
table(broad.fibroblasts$Sample)

## Keep samples with > 100 fibroblasts
samples.keep <- names(table(broad.fibroblasts$Sample))[which(table(broad.fibroblasts$Sample) > 100)]
samples.keep

broad.fibroblasts <- subset(broad.fibroblasts, Sample %in% samples.keep)
dim(broad.fibroblasts)

dataset.list <- Seurat::SplitObject(broad.fibroblasts, split.by = "Sample")

## Run log-normalization
dataset.list <- lapply(X = dataset.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
  return(x)
})

int.features <- SelectIntegrationFeatures(object.list = dataset.list, nfeatures = 2000)

broad.anchors <- FindIntegrationAnchors(object.list = dataset.list, 
                                        anchor.features = int.features,
                                        normalization.method = "LogNormalize")
broad.combined.seurat <- IntegrateData(anchorset = broad.anchors)

remove(broad.anchors)
remove(dataset.list)
remove(broad.fibroblasts)

DefaultAssay(broad.combined.seurat) <- "integrated"

broad.combined.seurat <- ScaleData(broad.combined.seurat)
broad.combined.seurat <- RunPCA(broad.combined.seurat, verbose = FALSE)
broad.combined.seurat <- RunUMAP(broad.combined.seurat, reduction = "pca", dims = 1:40)

broad.combined.seurat <- FindNeighbors(broad.combined.seurat, reduction = "pca", dims = 1:40)
broad.combined.seurat <- FindClusters(broad.combined.seurat, resolution = seq(0.2, 0.6, 0.1))

Idents(broad.combined.seurat) <- broad.combined.seurat$integrated_snn_res.0.4

DimPlot(broad.combined.seurat, label = TRUE)
DimPlot(broad.combined.seurat, label = TRUE, split.by = "Sample") + facet_wrap( ~Sample )

DefaultAssay(broad.combined.seurat) <- "RNA"
broad.combined.seurat <- NormalizeData(broad.combined.seurat)


## Glial cell markers
PlotExpressionUMAP(broad.combined.seurat, c("CHL1", "PTPRZ1", "CADM2", "ANK3"))

## Immune cells
PlotExpressionUMAP(broad.combined.seurat, c("PTPRC", "CD79A", "CD3D", "CSF3R"))

## Mural cells
PlotExpressionUMAP(broad.combined.seurat, c("RGS5", "PDGFRB", "VTN", "KCNJ8"))

## EC markers
PlotExpressionUMAP(broad.combined.seurat, c("KDR", "PECAM1", "CDH5", "EMCN"))

## CM markers
PlotExpressionUMAP(broad.combined.seurat, c("TNNC1", "MYH6", "FABP3", "TNNT2"))

DimPlot(broad.combined.seurat, label = TRUE)

## Subset high-confidence fibroblasts
## Cluster 7 = Glial and Cluster 5 = ECs
broad.combined.seurat <- subset(broad.combined.seurat, ident = c("0", "1", "2", "3", "4", "6", "8"))
DimPlot(broad.combined.seurat, label = TRUE)

broad.combined.seurat <- subset( broad.combined.seurat, UMAP_2 < 6 )
DimPlot(broad.combined.seurat, label = TRUE)

save(broad.combined.seurat, file = "Broad_Covid_fibroblasts_integrated_PC40.RData")



##############################################
### Subset the Izar COVID-infected samples ###
##############################################

izar.fibroblasts <- subset(covid.fibroblasts, Study == "Izar")
dim(izar.fibroblasts)
table(izar.fibroblasts$Sample)

## Keep samples with > 100 fibroblasts
samples.keep <- names(table(izar.fibroblasts$Sample))[which(table(izar.fibroblasts$Sample) > 100)]
samples.keep

izar.fibroblasts <- subset(izar.fibroblasts, Sample %in% samples.keep)
dim(izar.fibroblasts)

dataset.list <- Seurat::SplitObject(izar.fibroblasts, split.by = "Sample")

## Run log-normalization
dataset.list <- lapply(X = dataset.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
  return(x)
})

int.features <- SelectIntegrationFeatures(object.list = dataset.list, nfeatures = 2000)

izar.anchors <- FindIntegrationAnchors(object.list = dataset.list, 
                                       anchor.features = int.features,
                                       normalization.method = "LogNormalize")
izar.combined.seurat <- IntegrateData(anchorset = izar.anchors)

remove(izar.anchors)
remove(dataset.list)

DefaultAssay(izar.combined.seurat) <- "integrated"

izar.combined.seurat <- ScaleData(izar.combined.seurat)
izar.combined.seurat <- RunPCA(izar.combined.seurat, verbose = FALSE)
izar.combined.seurat <- RunUMAP(izar.combined.seurat, reduction = "pca", dims = 1:40)

izar.combined.seurat <- FindNeighbors(izar.combined.seurat, reduction = "pca", dims = 1:40)
izar.combined.seurat <- FindClusters(izar.combined.seurat, resolution = seq(0.2, 0.6, 0.1))

Idents(izar.combined.seurat) <- izar.combined.seurat$integrated_snn_res.0.4

DimPlot(izar.combined.seurat, label = TRUE)
DimPlot(izar.combined.seurat, label = TRUE, split.by = "Sample") + facet_wrap( ~Sample )

DefaultAssay(izar.combined.seurat) <- "RNA"
izar.combined.seurat <- NormalizeData(izar.combined.seurat)

## Glial cell markers
PlotExpressionUMAP(izar.combined.seurat, c("CHL1", "PTPRZ1", "CADM2", "ANK3"))

## Immune cells
PlotExpressionUMAP(izar.combined.seurat, c("PTPRC", "CD79A", "CD3D", "CSF3R"))

## Mural cells
PlotExpressionUMAP(izar.combined.seurat, c("RGS5", "PDGFRB", "VTN", "KCNJ8"))

## EC markers
PlotExpressionUMAP(izar.combined.seurat, c("KDR", "PECAM1", "CDH5", "EMCN"))

## CM markers
PlotExpressionUMAP(izar.combined.seurat, c("TNNC1", "MYH6", "FABP3", "TNNT2"))

DimPlot(izar.combined.seurat, label = TRUE)

### Subset high-confidence fibroblasts
izar.combined.seurat <- subset(izar.combined.seurat, ident = c( "0", "1", "2", "3", "4", "5", "6" ))

DimPlot(izar.combined.seurat, label = TRUE)

save(izar.combined.seurat, file = "Izar_Covid_fibroblasts_integrated_PC40.RData")

