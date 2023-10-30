library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)


################################
### Broad COVID fibroblasts  ###
################################

load("Heart_integrated_fibroblast_singlets_Seurat.RData")
DimPlot(fibroblasts.integrated, label = TRUE)
DefaultAssay(fibroblasts.integrated) <- "integrated"

load("Broad_Covid_fibroblasts_integrated_PC40.RData")
DimPlot(broad.combined.seurat, label = TRUE)
DefaultAssay(broad.combined.seurat) <- "RNA"

## Label transfer from healthy heart integrated fibroblasts
fibro.anchors <- FindTransferAnchors(reference = fibroblasts.integrated, 
                                     query = broad.combined.seurat, 
                                     dims = 1:40)
predictions <- TransferData(anchorset = fibro.anchors, 
                            refdata = Idents(fibroblasts.integrated), 
                            dims = 1:40)
colnames(predictions) <- paste0("Porello_", colnames(predictions))

broad.combined.seurat <- AddMetaData(broad.combined.seurat, metadata = predictions)

broad.combined.seurat@meta.data$Porello_predicted.id <- factor(broad.combined.seurat@meta.data$Porello_predicted.id,
                                                               levels = as.character(seq(1, 5)))
DimPlot(broad.combined.seurat, group.by = "Porello_predicted.id", label = TRUE)

Idents(broad.combined.seurat) <- broad.combined.seurat@meta.data$Porello_predicted.id
table(Idents(broad.combined.seurat))

DimPlot(broad.combined.seurat, label = TRUE, reduction = "umap")

FeaturePlot(broad.combined.seurat, "Porello_prediction.score.2")

PlotExpressionUMAP(broad.combined.seurat, c("POSTN", "MEOX1", "CILP", "COL8A1"))

doBoxPlot(broad.combined.seurat, c("POSTN", "MEOX1", "CILP", "COL8A1"))
doBoxPlot(broad.combined.seurat, c("SCARA5", "CD248", "CD55", "GFPT2"))

save( broad.combined.seurat, file = "Broad_Covid_fibroblasts_integrated_PC40.RData" )

###############################
### Izhar COVID fibroblasts ###
###############################

load("Heart_integrated_fibroblast_singlets_Seurat.RData")
DimPlot(fibroblasts.integrated, label = TRUE)
DefaultAssay(fibroblasts.integrated) <- "integrated"

load("Izar_Covid_fibroblasts_integrated_PC40.RData")
DimPlot(izar.combined.seurat, label = TRUE)
DefaultAssay(izar.combined.seurat) <- "RNA"


fibro.anchors <- FindTransferAnchors(reference = fibroblasts.integrated, 
                                     query = izar.combined.seurat, 
                                     dims = 1:40)
predictions <- TransferData(anchorset = fibro.anchors, 
                            refdata = Idents(fibroblasts.integrated), 
                            dims = 1:40)
colnames(predictions) <- paste0("Porello_", colnames(predictions))

izar.combined.seurat <- AddMetaData(izar.combined.seurat, metadata = predictions)

izar.combined.seurat@meta.data$Porello_predicted.id <- factor(izar.combined.seurat@meta.data$Porello_predicted.id,
                                                              levels = as.character(seq(1, 5)))
DimPlot(izar.combined.seurat, group.by = "Porello_predicted.id", label = TRUE)

Idents(izar.combined.seurat) <- izar.combined.seurat@meta.data$Porello_predicted.id
table(Idents(izar.combined.seurat))

DimPlot(izar.combined.seurat, label = TRUE, reduction = "umap")

FeaturePlot(izar.combined.seurat, "Porello_prediction.score.2")

PlotExpressionUMAP(izar.combined.seurat, c("POSTN", "MEOX1", "CILP", "COL8A1"))

doBoxPlot(izar.combined.seurat, c("POSTN", "MEOX1", "CILP", "COL8A1"))
doBoxPlot(izar.combined.seurat, c("SCARA5", "CD248", "CD55", "GFPT2"))

save( izar.combined.seurat, file = "Izar_Covid_fibroblasts_integrated_PC40.RData" )

#####################################
### Teichman healthy fibroblasts  ###
#####################################

load("../Porello_human_heart/Heart_integrated_fibroblast_singlets_Seurat.RData")
DimPlot(fibroblasts.integrated, label = TRUE)
DefaultAssay(fibroblasts.integrated) <- "integrated"

load("Teichman_healthy_fibroblasts_integrated_PC40.RData")
DimPlot(teichman.combined.seurat, label = TRUE)
DefaultAssay(teichman.combined.seurat) <- "RNA"


fibro.anchors <- FindTransferAnchors(reference = fibroblasts.integrated, 
                                     query = teichman.combined.seurat, 
                                     dims = 1:40)
predictions <- TransferData(anchorset = fibro.anchors, 
                            refdata = Idents(fibroblasts.integrated), 
                            dims = 1:40)
colnames(predictions) <- paste0("Porello_", colnames(predictions))

teichman.combined.seurat <- AddMetaData(teichman.combined.seurat, metadata = predictions)

teichman.combined.seurat@meta.data$Porello_predicted.id <- factor(teichman.combined.seurat@meta.data$Porello_predicted.id,
                                                                  levels = as.character(seq(1, 5)))
DimPlot(teichman.combined.seurat, group.by = "Porello_predicted.id", label = TRUE)

Idents(teichman.combined.seurat) <- teichman.combined.seurat@meta.data$Porello_predicted.id
table(Idents(teichman.combined.seurat))

DimPlot(teichman.combined.seurat, label = TRUE, reduction = "umap")

save( teichman.combined.seurat, file = "Teichman_healthy_fibroblasts_integrated_PC40.RData" )


####################################
### Ellinor healthy fibroblasts  ###
####################################

load("/Heart_integrated_fibroblast_singlets_Seurat.RData")
DimPlot(fibroblasts.integrated, label = TRUE)
DefaultAssay(fibroblasts.integrated) <- "integrated"

load("Ellinor_healthy_fibroblasts_integrated_PC40.RData")
DimPlot(ellinor.combined.seurat, label = TRUE)
DefaultAssay(ellinor.combined.seurat) <- "RNA"


fibro.anchors <- FindTransferAnchors(reference = fibroblasts.integrated, 
                                     query = ellinor.combined.seurat, 
                                     dims = 1:40)
predictions <- TransferData(anchorset = fibro.anchors, 
                            refdata = Idents(fibroblasts.integrated), 
                            dims = 1:40)
colnames(predictions) <- paste0("Porello_", colnames(predictions))

ellinor.combined.seurat <- AddMetaData(ellinor.combined.seurat, metadata = predictions)

ellinor.combined.seurat@meta.data$Porello_predicted.id <- factor(ellinor.combined.seurat@meta.data$Porello_predicted.id,
                                                                 levels = as.character(seq(1, 5)))
DimPlot(ellinor.combined.seurat, group.by = "Porello_predicted.id", label = TRUE)

Idents(ellinor.combined.seurat) <- ellinor.combined.seurat@meta.data$Porello_predicted.id
table(Idents(ellinor.combined.seurat))

#DimPlot(teichman.combined.seurat, label = TRUE, reduction = "umap")

PlotExpressionUMAP(ellinor.combined.seurat, c("POSTN", "MEOX1", "CILP", "COL8A1"))
PlotExpressionUMAP(ellinor.combined.seurat, c("SCARA5", "CD248", "CD55", "GFPT2"))

doBoxPlot(ellinor.combined.seurat, c("POSTN", "MEOX1", "CILP", "COL8A1"))
doBoxPlot(ellinor.combined.seurat, c("SCARA5", "CD248", "CD55", "GFPT2"))

save( ellinor.combined.seurat, file = "Ellinor_healthy_fibroblasts_integrated_PC40.RData" )


######################################################################################
### Build a combined Seurat object and evaluate differences in cluster proportions ###
######################################################################################

## Create a combined object of all fibroblasts incorporating the reference set
load("Heart_integrated_fibroblast_singlets_Seurat.RData")
fibroblasts.integrated@meta.data$Study <- rep("Porrello", ncol(fibroblasts.integrated))
fibroblasts.integrated@meta.data$Sample <- fibroblasts.integrated$orig.ident
fibroblasts.integrated@meta.data$Disease <- rep("healthy", ncol(fibroblasts.integrated))
fibroblasts.integrated@meta.data$Porello_predicted.id <- Idents(fibroblasts.integrated)


load("Ellinor_healthy_fibroblasts_integrated_PC40.RData")
load("Izar_Covid_fibroblasts_integrated_PC40.RData")
load("Broad_Covid_fibroblasts_integrated_PC40.RData")
load("Teichman_healthy_fibroblasts_integrated_PC40.RData")

all.fibroblasts.combined <- merge(x = fibroblasts.integrated, 
                                  y = c(teichman.combined.seurat, broad.combined.seurat, ellinor.combined.seurat, izar.combined.seurat),
                                  add.cell.ids = c("Porrello", "Teichman", "Ellinor", "Broad", "Izar"))

groups <- all.fibroblasts.combined@meta.data[, c("Porello_predicted.id", "Sample", "Disease", "Study")]
groups$Sample <- stringr::str_replace_all(groups$Sample, "_", "-")
table(groups$Sample)
table(groups$Study)
table(groups$Porello_predicted.id)


### Run differential proportion analysis using Propeller
library(speckle)
library(limma)

cluster.table <- all.fibroblasts.combined@meta.data[, c("Porello_predicted.id", "Sample")]
cluster.table %>%
  dplyr::group_by(Sample, Porello_predicted.id, .drop = FALSE) %>%
  summarise(n=n()) %>% 
  group_by(Sample) %>% 
  mutate(nGroup = sum(n), Percentage = n / sum(n) * 100) %>%
  as.data.frame() -> cluster.summary

all.fibroblasts.combined@meta.data[, c("Sample", "Disease")] %>%
  distinct(Sample, .keep_all = TRUE) -> sample.table
rownames(sample.table) <- sample.table$Sample


### Proportion comparisons for all conditions
cluster.table <- all.fibroblasts.combined@meta.data[, c("Porello_predicted.id", "Sample", "Disease")]

prop.res <- propeller(clusters = cluster.table$Porello_predicted.id, sample = cluster.table$Sample, 
                      group = cluster.table$Disease)

write.csv(prop.res, file = "results/propeller_results.csv")

pl <- plotCellTypeProps(clusters = cluster.table$Porello_predicted.id, sample = cluster.table$Sample)
ggData <- pl$data
sample.table.extended <- sample.table[ggData$Sample, ]

ggData$Condition <- sample.table$Disease
ggData$Condition <- plyr::mapvalues(ggData$Condition, "healthy", to = "Healthy")

ggplot(ggData, aes(x=Condition, y=Proportions*100, fill=Condition)) + 
  geom_point(aes(colour = Condition)) +
  facet_wrap(~Clusters, nrow=2) + xlab("Population") + ylab("% of total cells") + 
  ggtitle("Covid labels proportions") +
  theme_classic(base_size = 16) +
  theme(legend.position="none", strip.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

ggplot(ggData, aes(x=Condition, y=Proportions*100, fill=Condition)) + 
  geom_boxplot(colour = "black") + geom_jitter() +
  facet_wrap(~Clusters, nrow=1) + xlab("Population") + ylab("% of total cells") + 
  ggtitle("Covid vs healthy cluster proportions") +
  theme_classic(base_size = 16) +
  theme(legend.position="none", strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))

ggplot(ggData, aes(x=Condition, y=Proportions*100, fill=Condition)) + 
  geom_violin(colour = "black") + geom_jitter() +
  facet_wrap(~Clusters, nrow=1) + xlab("Population") + ylab("% of total cells") + 
  ggtitle("Covid vs healthy cluster proportions") +
  theme_classic(base_size = 16) +
  theme(legend.position="none", strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))



