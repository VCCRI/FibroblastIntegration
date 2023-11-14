library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)
library(cFIT)

load("Forte_GFPD3_GFPD7_CTHRC1_Hif1a_varGenes4000_r15_cFIT.RData")
var.genes <- rownames(int.out$W)

###############################
### Forte et al fibroblasts ###
###############################
load("Forte/fibroblasts_no_norm_singlets_NegBinomReg_nCount_dims1-30.RData")
DimPlot(fibroblasts.seurat, label = TRUE)
forte.fibroblasts <- fibroblasts.seurat
remove(fibroblasts.seurat)

forte.fibroblasts$ClusteringIdent <- Idents(forte.fibroblasts)
forte.cell.idents <- Idents(forte.fibroblasts)
forte.genes.expressed <- rownames(forte.fibroblasts)
forte.conditions <- forte.fibroblasts@meta.data$Condition

##############################
###### Farbehi et al. data ###
##############################
load("GFP_day3/GFP_day3_Seurat.RData")
DimPlot(genes.seurat, label = TRUE)
gfp.day3 = genes.seurat
remove(genes.seurat)

gfp.day3$ClusteringIdent <- Idents(gfp.day3)

gfp.day3.cell.idents <- Idents(gfp.day3)
gfp.day3.genes.expressed <- rownames(gfp.day3)
gfp.day3.conditions <- gfp.day3@meta.data$Condition

##############################
### Read in GFP-day 7 data ###
##############################
load("GFP_day7/GFP_day7_Seurat.RData")
DimPlot(genes.seurat, label = TRUE)
gfp.day7 = genes.seurat
remove(genes.seurat)


gfp.day7$ClusteringIdent <- Idents(gfp.day7)
gfp.day7.cell.idents <- Idents(gfp.day7)
gfp.day7.genes.expressed <- rownames(gfp.day7)
gfp.day7.conditions <- gfp.day7@meta.data$Condition

######################################
### Hif1a WT sham vs MI day 3 data ###
######################################
load("TdTom_day3/tdTom_day3_Seurat.RData")
DimPlot(genes.seurat, label = TRUE)
tdTom.day3 = genes.seurat
remove(genes.seurat)

tdTom.day3$ClusteringIdent <- Idents(tdTom.day3)

tdTom.day3.cell.idents <- Idents(tdTom.day3)
tdTom.day3.genes.expressed <- rownames(tdTom.day3)
tdTom.day3.conditions <- tdTom.day3@meta.data$Condition

#####################################
### Read in CTHRC1 MI-day 14 data ###
#####################################
load("CTHRC1/fibroblasts_MI_day14_Seurat.RData")
DimPlot(mid14.fibroblasts, label = TRUE)

mid14.fibroblasts$ClusterID <- mid14.fibroblasts$FibroblastClusters
table(mid14.fibroblasts$ClusterID)
Idents(mid14.fibroblasts) <- mid14.fibroblasts$ClusterID

mid14.fibroblasts$ClusteringIdent <- Idents(mid14.fibroblasts)
mid14.fibroblasts.cell.idents <- Idents(mid14.fibroblasts)
mid14.fibroblasts.genes.expressed <- rownames(mid14.fibroblasts)
mid14.fibroblasts.conditions <- rep("MI-day 14", ncol(mid14.fibroblasts))

#####################################
### Read in CTHRC1 MI-day 30 data ###
#####################################
load("CTHRC1/fibroblasts_MI_day30_Seurat.RData")
DimPlot(mid30.fibroblasts, label = TRUE, group.by = "FibroblastClusters")

table(mid30.fibroblasts$FibroblastClusters)
Idents(mid30.fibroblasts) <- mid30.fibroblasts$FibroblastClusters

mid30.fibroblasts$ClusteringIdent <- Idents(mid30.fibroblasts)
mid30.fibroblasts.cell.idents <- Idents(mid30.fibroblasts)
mid30.fibroblasts.genes.expressed <- rownames(mid30.fibroblasts)
mid30.fibroblasts.conditions <- rep("MI-day 30", ncol(mid30.fibroblasts))

########################################################################
### Regenerate the Seurat objects with variable genes to save memory ###
########################################################################

#overlapping.genes <- Reduce(union, list(gfp.day3.genes.expressed, gfp.day7.genes.expressed,
#                                        forte.genes.expressed, tdTom.day3.genes.expressed, 
#                                        mid14.fibroblasts.genes.expressed, mid30.fibroblasts.genes.expressed))
overlapping.genes <- Reduce(intersect, list(gfp.day3.genes.expressed, gfp.day7.genes.expressed,
                                            forte.genes.expressed, tdTom.day3.genes.expressed, 
                                            mid14.fibroblasts.genes.expressed, mid30.fibroblasts.genes.expressed))
length(overlapping.genes)

### GFP-day 3 data
count.data <- Read10X("GFP_day3/GFP_day3_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(gfp.day3), ]
dim(count.data)
gfp.day3 <- CreateSeuratObject(counts = count.data, project = "GFP-day3")
gfp.day3@meta.data$Condition <- gfp.day3.conditions
gfp.day3@meta.data$ClusterID <- gfp.day3.cell.idents
gfp.day3@meta.data$Experiment <- rep("GFPD3", ncol(gfp.day3))

### GFP-day 7 data
count.data <- Read10X("GFP_day7/GFP_day7_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(gfp.day7), ]
dim(count.data)
gfp.day7 <- CreateSeuratObject(counts = count.data, project = "GFP-day7")
gfp.day7@meta.data$Condition <- gfp.day7.conditions
gfp.day7@meta.data$ClusterID <- gfp.day7.cell.idents
gfp.day7@meta.data$Experiment <- rep("GFPD7", ncol(gfp.day7))

### TdTom day 3 data
count.data <- Read10X("TdTom_day3/TdTom_day3_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(tdTom.day3), ]
dim(count.data)
tdTom.day3 <- CreateSeuratObject(counts = count.data, project = "tdTom-day3")
tdTom.day3@meta.data$Condition <- tdTom.day3.conditions
tdTom.day3@meta.data$ClusterID <- tdTom.day3.cell.idents
tdTom.day3@meta.data$Experiment <- rep("tdTomD3", ncol(tdTom.day3))

### Forte data
count.data <- Read10X("Forte/MI_timecourse_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(forte.fibroblasts), ]
dim(count.data)
forte.fibroblasts <- CreateSeuratObject(counts = count.data, project = "FORTE")
forte.fibroblasts@meta.data$Condition <- forte.conditions
forte.fibroblasts@meta.data$ClusterID <- forte.cell.idents
forte.fibroblasts@meta.data$Experiment <- rep("FORTE", ncol(forte.fibroblasts))

### CTHRC1 MI-day 14 data
count.data <- Read10X("CTHRC1/MI_day14_CR6/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(mid14.fibroblasts), ]
dim(count.data)
mid14.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day14")
mid14.fibroblasts$Condition <- mid14.fibroblasts.conditions
mid14.fibroblasts@meta.data$ClusterID <- mid14.fibroblasts.cell.idents
mid14.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID14", ncol(mid14.fibroblasts))

### CTHRC1 MI-day 30 data
count.data <- Read10X("CTHRC1/MI_day30_CR6/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(mid30.fibroblasts), ]
dim(count.data)
mid30.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day30")
mid30.fibroblasts$Condition <- mid30.fibroblasts.conditions
mid30.fibroblasts@meta.data$ClusterID <- mid30.fibroblasts.cell.idents
mid30.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID30", ncol(mid30.fibroblasts))

fibroblasts.combined <- merge(x = gfp.day3, y = c(gfp.day7, forte.fibroblasts, tdTom.day3,
                                                  mid14.fibroblasts, mid30.fibroblasts),
                              add.cell.ids = c("GFPD3", "GFPD7", "FORTE", "tdTomD3",
                                               "CTHRC1D14", "CTHRC1D30"))
table(fibroblasts.combined$ClusterID)

fibroblasts.combined$ClusterID <- plyr::mapvalues(fibroblasts.combined$ClusterID,
                                                  from = c("MYO-1", "MYO-2", "MYO-3", "Cyc",
                                                           "MFC-1", "MFC-2", "MFC-3"),
                                                  to = c("MYO", "MYO", "MYO", "F-Cyc",
                                                         "MFC", "MFC", "MFC"))
table(fibroblasts.combined$ClusterID)
Idents(fibroblasts.combined) <- fibroblasts.combined$ClusterID

table(fibroblasts.combined$Condition)
fibroblasts.combined$Condition <- plyr::mapvalues(fibroblasts.combined$Condition,
                                                  from = "Sham", to = "Healthy")
table(fibroblasts.combined$Condition)

#clusters.keep <- setdiff(unique(fibroblasts.combined$ClusterID), c("EC", "EPI", "MAC"))
#clusters.keep

#fibroblasts.combined <- subset(fibroblasts.combined, idents = clusters.keep)
#table(fibroblasts.combined$ClusterID)

remove(gfp.day3)
remove(gfp.day7)
remove(forte.fibroblasts)
remove(tdTom.day3)
remove(mid14.fibroblasts)
remove(mid30.fibroblasts)

fibroblasts.combined <- NormalizeData(fibroblasts.combined)

#########################
### Run cFIT analysis ###
#########################

# extract the raw counts and metadata for data sets from 5 technologies
data.list = split_dataset_by_batch(X=t(as.matrix(fibroblasts.combined@assays$RNA@counts)),
                                   batch = fibroblasts.combined@meta.data$Experiment,
                                   labels = fibroblasts.combined@meta.data$ClusterID,
                                   metadata = fibroblasts.combined@meta.data,
                                   dataset.name = 'fibroblast:')

exprs.list = preprocess_for_integration(data.list$X.list, var.genes, scale.factor=10^4, scale=T, center=F)

######################################################
### Read in the integrated data and run UMAP ###
######################################################

load("MI_datasets_cFIT_integration_Seurat_intersect_genes.RData")

exprs.int = do.call(rbind, int.out$H.list) %*% t(int.out$W)

# ncell-by-r low dimensiional representation
Hnorm = do.call(rbind, int.out$H.list) %*% diag(colSums(int.out$W))

exprs.int = t(exprs.int)
exprs.int <- exprs.int[, colnames(fibroblasts.combined)]

cfit.assay <- Seurat::CreateAssayObject(data = exprs.int)
fibroblasts.combined@assays$cFIT <- cfit.assay

colnames(Hnorm) <- paste0("HN_", seq(1:ncol(Hnorm)))
hnorm.object <- CreateDimReducObject(embeddings = Hnorm, loadings = int.out$W,
                                     key = "HN_", assay = "cFIT")
fibroblasts.combined$cFITnorm <- hnorm.object

DefaultAssay(fibroblasts.combined) <- "cFIT"

fibroblasts.combined <- RunUMAP(fibroblasts.combined, 
                                n.neighbors = 200, 
                                min.dist = 0.5,
                                metric = "cosine",
                                reduction = "cFITnorm", 
                                dims = 1:15)


Idents(fibroblasts.combined) <- fibroblasts.combined$ClusterID

col.set <- c("#fde800", "#fb8500", "#00b600", "#00896e", "#f60000", "#d647bf",
             "#3bfff2", "#0053c8", "#0099cc", "#a7b5c9", "pink")
pops <- c("F-SH", "F-SL", "F-Trans", "F-WntX", "F-Act", "IR", "F-CI", "F-Cyc", "MYO", "MFC", "F-IFNS")
names(col.set) <- pops 
pl1 <- DimPlot(fibroblasts.combined, label = TRUE, cols = col.set)
pl1

ggData <- pl1$data
ggData$Condition <- fibroblasts.combined$Condition
ggData$Condition <- factor(ggData$Condition, levels = c("Healthy", "MI-day 1", "MI-day 3", "MI-day 5",
                                                        "MI-day 7", "MI-day 14", "MI-day 28", "MI-day 30"))
ggData$Condition <- plyr::mapvalues(ggData$Condition, 
                                    from = c( "MI-day 28", "MI-day 30" ),
                                    to = c("MI-day 28/30", "MI-day 28/30"))

ggData2 <- ggData
ggData2$Condition <- rep("Combined", nrow(ggData2))
ggData2 <- rbind(ggData2, ggData)
ggData2$Condition <- factor(ggData2$Condition, levels = c("Combined","Healthy", "MI-day 1", "MI-day 3", "MI-day 5",
                                                          "MI-day 7", "MI-day 14", "MI-day 28/30"))
pl1 <- ggplot(ggData2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 16) +
  facet_wrap( ~Condition, ncol = 4 ) +  
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(strip.background = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
pl1 

pl2 <- ggplot(ggData, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 16) +
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4))) 

ggData %>%
  dplyr::group_by(ident) %>%
  summarize(UMAP_1 = median(x = UMAP_1), UMAP_2 = median(x = UMAP_2)) -> centers

pl2 <- pl2 + geom_text(data = centers, mapping = aes(label = ident), size = 4.5, colour="black", fontface="bold") +
  theme(text = element_text(size = 14, family="Helvetica")) + xlab("UMAP 1") + ylab("UMAP 2")
pl2 

pl1 <- pl1 + theme(legend.position = "none") 
pl2 <- pl2 + theme(legend.position = "none") + ggtitle("Integrated fibroblasts")

plot_grid(pl2, pl1, ncol=2, align = "h", rel_widths = c( 0.35, 0.65 ))



#####################################################
### Update the cluster IDs based on KNN matching ###
#####################################################

Idents(fibroblasts.combined) <- fibroblasts.combined$ClusterID

populations <- names(table(Idents(fibroblasts.combined)))

current.idents <- Idents(fibroblasts.combined)
new.idents <- as.character(current.idents)
names(new.idents) <- names(current.idents)

k = 25
verbose = T

## Get KNN graph
knn.out = FNN::get.knn(do.call(rbind, int.out$H.list), k = k)
knn.index = knn.out$nn.index

## Identify the nearest neighbor for each cell in the cluster
labels <- Idents(fibroblasts.combined)
set.seed(1)
for (this.cell in populations) {
  print(this.cell)
  
  cl.labels <- which(labels == this.cell)
  pred.class <- unlist(lapply(cl.labels, function(x) {
    ind = knn.index[x, ]
    this.nn <- factor(labels)[ind[!is.na(ind)]]
    nn.norm <- table(this.nn)
    this.pred <- names(nn.norm)[which(nn.norm == max(nn.norm))]
    if (length(this.pred) > 1) {
      ## In rare cases of more than one that one cell type matching, 
      ## randomly select one of the best matches
      this.pred <- sample(this.pred, size = 1)
    }
    return(this.pred)
  }))
  print(table(pred.class))
  print(length(cl.labels))
  print(length(pred.class))
  
  new.idents[cl.labels] <- pred.class
  
}
table(current.idents)

table(new.idents)

Idents(fibroblasts.combined) <- new.idents

fibroblasts.combined$Clusters_updated_knn_basic_k25 <- new.idents
Idents(fibroblasts.combined) <- fibroblasts.combined$Clusters_updated_knn_basic_k25


### Redo the UMAP with the updated cluster IDs
ggData <- pl1$data
ggData$Condition <- fibroblasts.combined$Condition
ggData$Condition <- factor(ggData$Condition, levels = c("Healthy", "MI-day 1", "MI-day 3", "MI-day 5",
                                                        "MI-day 7", "MI-day 14", "MI-day 28", "MI-day 30"))
ggData$Condition <- plyr::mapvalues(ggData$Condition, 
                                    from = c( "MI-day 28", "MI-day 30" ),
                                    to = c("MI-day 28/30", "MI-day 28/30"))

ggData2 <- ggData
ggData2$Condition <- rep("Combined", nrow(ggData2))
ggData2 <- rbind(ggData2, ggData)
ggData2$Condition <- factor(ggData2$Condition, levels = c("Combined","Healthy", "MI-day 1", "MI-day 3", "MI-day 5",
                                                          "MI-day 7", "MI-day 14", "MI-day 28/30"))
pl1 <- ggplot(ggData2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 16) +
  facet_wrap( ~Condition, ncol = 4 ) +  
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(strip.background = element_blank()) + xlab("UMAP 1") + ylab("UMAP 2")
pl1 

pl2 <- ggplot(ggData, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 16) +
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4))) 

ggData %>%
  dplyr::group_by(ident) %>%
  summarize(UMAP_1 = median(x = UMAP_1), UMAP_2 = median(x = UMAP_2)) -> centers

pl2 <- pl2 + geom_text(data = centers, mapping = aes(label = ident), size = 4.5, colour="black", fontface="bold") +
  theme(text = element_text(size = 14, family="Helvetica")) + xlab("UMAP 1") + ylab("UMAP 2")
pl2 

pl1 <- pl1 + theme(legend.position = "none") 
pl2 <- pl2 + theme(legend.position = "none") + ggtitle("Integrated fibroblasts")

plot_grid(pl2, pl1, ncol=2, align = "h", rel_widths = c( 0.35, 0.65 ))

