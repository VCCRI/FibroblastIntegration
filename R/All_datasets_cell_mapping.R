library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)
library(cFIT)

load("All_datasets_IR_2samples_cFIT_varGenes4000_r15_cFITSketched05.RData")
var.genes <- rownames(int.out$W)

###############################
### Forte et al fibroblasts ###
###############################
load("Forte/fibroblasts_MI_timecourse.RData.RData")
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


####################################
### Read in AngII treatment data ###
####################################
load("AngII_Pinto/fibroblasts_AngII.RData")
DimPlot(fibroblasts.seurat, label = TRUE)
angII.fibroblasts <- fibroblasts.seurat
remove(fibroblasts.seurat)

angII.fibroblasts$ClusteringIdent <- Idents(angII.fibroblasts)
angII.fibroblasts.cell.idents <- Idents(angII.fibroblasts)
angII.fibroblasts.genes.expressed <- rownames(angII.fibroblasts)
angII.fibroblasts.conditions <- angII.fibroblasts$Condition


########################
### Read in TAC data ###
########################
load("TAC/fibroblasts_TAC.RData")
DimPlot(fibroblasts.seurat, label = TRUE)
tac.fibroblasts <- fibroblasts.seurat
remove(fibroblasts.seurat)

tac.fibroblasts$ClusteringIdent <- Idents(tac.fibroblasts)
tac.fibroblasts.cell.idents <- Idents(tac.fibroblasts)
tac.fibroblasts.genes.expressed <- rownames(tac.fibroblasts)
tac.fibroblasts.conditions <- tac.fibroblasts$Condition

#######################
### Read in IR data ###
#######################
load("ischaemia_repurfusion/IR_integrated_fibroblasts.RData")
DimPlot(fibroblasts.seurat, label = TRUE)
ir.fibroblasts <- fibroblasts.seurat
remove(fibroblasts.seurat)

ir.fibroblasts <- subset(ir.fibroblasts, Experiment %in% c("IR_1", "IR_3", "Sham_1", "Sham_3"))
dim(ir.fibroblasts)

ir.fibroblasts$ClusteringIdent <- Idents(ir.fibroblasts)
ir.fibroblasts.cell.idents <- Idents(ir.fibroblasts)
ir.fibroblasts.genes.expressed <- rownames(ir.fibroblasts)
ir.fibroblasts.conditions <- ir.fibroblasts$Condition

########################################################################
### Regenerate the Seurat objects with variable genes to save memory ###
########################################################################

overlapping.genes <- Reduce(intersect, list(gfp.day3.genes.expressed, gfp.day7.genes.expressed,
                                            forte.genes.expressed, tdTom.day3.genes.expressed, 
                                            mid14.fibroblasts.genes.expressed, mid30.fibroblasts.genes.expressed,
                                            angII.fibroblasts.genes.expressed, tac.fibroblasts.genes.expressed,
                                            ir.fibroblasts.genes.expressed))
length(overlapping.genes)

### GFP-day 3 data
count.data <- Read10X("GFP_day3/GFP_day3_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(gfp.day3), ]
dim(count.data)
gfp.day3 <- CreateSeuratObject(counts = count.data, project = "GFP-day3")
gfp.day3@meta.data$ClusterID <- gfp.day3.cell.idents
gfp.day3@meta.data$Experiment <- rep("GFPD3", ncol(gfp.day3))
gfp.day3@meta.data$Condition <- gfp.day3.conditions

### GFP-day 7 data
count.data <- Read10X("GFP_day7/GFP_day7_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(gfp.day7), ]
dim(count.data)
gfp.day7 <- CreateSeuratObject(counts = count.data, project = "GFP-day7")
gfp.day7@meta.data$ClusterID <- gfp.day7.cell.idents
gfp.day7@meta.data$Experiment <- rep("GFPD7", ncol(gfp.day7))
gfp.day7@meta.data$Condition <- gfp.day7.conditions

### TdTom day 3 data
count.data <- Read10X("TdTom_day3/TdTom_day3_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(tdTom.day3), ]
dim(count.data)
tdTom.day3 <- CreateSeuratObject(counts = count.data, project = "tdTom-day3")
tdTom.day3@meta.data$ClusterID <- tdTom.day3.cell.idents
tdTom.day3@meta.data$Experiment <- rep("tdTomD3", ncol(tdTom.day3))
tdTom.day3@meta.data$Condition <- tdTom.day3.conditions

### Forte data
count.data <- Read10X("Forte/MI_timecourse_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(forte.fibroblasts), ]
dim(count.data)
forte.fibroblasts <- CreateSeuratObject(counts = count.data, project = "FORTE")
forte.fibroblasts@meta.data$ClusterID <- forte.cell.idents
forte.fibroblasts@meta.data$Experiment <- rep("FORTE", ncol(forte.fibroblasts))
forte.fibroblasts@meta.data$Condition <- forte.conditions

### CTHRC1 MI-day 14 data
count.data <- Read10X("CTHRC1/MI_day14_CR6/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(mid14.fibroblasts), ]
dim(count.data)
mid14.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day14")
mid14.fibroblasts@meta.data$ClusterID <- mid14.fibroblasts.cell.idents
mid14.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID14", ncol(mid14.fibroblasts))
mid14.fibroblasts@meta.data$Condition <- mid14.fibroblasts.conditions

### CTHRC1 MI-day 30 data
count.data <- Read10X("CTHRC1/MI_day30_CR6/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(mid30.fibroblasts), ]
dim(count.data)
mid30.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day30")
mid30.fibroblasts@meta.data$ClusterID <- mid30.fibroblasts.cell.idents
mid30.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID30", ncol(mid30.fibroblasts))
mid30.fibroblasts@meta.data$Condition <- mid30.fibroblasts.conditions

### AngII treatment data 
count.data <- Read10X("AngII_Pinto/AngII_treat_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(angII.fibroblasts), ]
dim(count.data)
angII.fibroblasts <- CreateSeuratObject(counts = count.data, project = "AngII")
angII.fibroblasts@meta.data$Condition <- angII.fibroblasts.conditions
angII.fibroblasts@meta.data$ClusterID <- angII.fibroblasts.cell.idents
angII.fibroblasts@meta.data$Experiment <- rep("AngII", ncol(angII.fibroblasts))
remove(count.data)

### TAC data
count.data <- Read10X("TAC/TAC_vs_sham_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(tac.fibroblasts), ]
dim(count.data)
tac.fibroblasts <- CreateSeuratObject(counts = count.data, project = "TAC")
tac.fibroblasts@meta.data$Condition <- tac.fibroblasts.conditions
tac.fibroblasts@meta.data$ClusterID <- tac.fibroblasts.cell.idents
tac.fibroblasts@meta.data$Experiment <- rep("TAC", ncol(tac.fibroblasts))
remove(count.data)

### IR data
count.data <- Read10X("ischaemia_repurfusion/IR_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(ir.fibroblasts), ]
dim(count.data)
ir.fibroblasts <- CreateSeuratObject(counts = count.data, project = "IR")
ir.fibroblasts@meta.data$Condition <- ir.fibroblasts.conditions
ir.fibroblasts@meta.data$ClusterID <- ir.fibroblasts.cell.idents
ir.fibroblasts@meta.data$Experiment <- rep("IR", ncol(ir.fibroblasts))
remove(count.data)

fibroblasts.combined <- merge(x = gfp.day3, y = c(gfp.day7, forte.fibroblasts, tdTom.day3,
                                                  mid14.fibroblasts, mid30.fibroblasts,
                                                  angII.fibroblasts, tac.fibroblasts,
                                                  ir.fibroblasts),
                              add.cell.ids = c("GFPD3", "GFPD7", "FORTE", "tdTomD3",
                                               "CTHRC1D14", "CTHRC1D30", "ANGII", "TAC", "IR"))
table(fibroblasts.combined$ClusterID)

fibroblasts.combined$ClusterID <- plyr::mapvalues(fibroblasts.combined$ClusterID,
                                                  from = c("MYO-1", "MYO-2", "MYO-3", "Cyc",
                                                           "MFC-1", "MFC-2", "MFC-3"),
                                                  to = c("MYO", "MYO", "MYO", "F-Cyc",
                                                         "MFC", "MFC", "MFC"))
table(fibroblasts.combined$ClusterID)
Idents(fibroblasts.combined) <- fibroblasts.combined$ClusterID

fibroblasts.combined$Condition <- plyr::mapvalues(fibroblasts.combined$Condition,
                                                  from = c("Sham", "MI-day 28", "MI-day 30"), 
                                                  to = c("Healthy", "MI-day 28/30", "MI-day 28/30"))
table(fibroblasts.combined$Condition)


remove(gfp.day3)
remove(gfp.day7)
remove(forte.fibroblasts)
remove(tdTom.day3)
remove(mid14.fibroblasts)
remove(mid30.fibroblasts)
remove(angII.fibroblasts)
remove(tac.fibroblasts)
remove(ir.fibroblasts)

#########################
### Run cFIT analysis ###
#########################

meta.data <- read.csv("all_datasets_metadata.csv", row.names = 1)
head(meta.data)
fibroblasts.combined@meta.data <- meta.data

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

fibroblasts.combined$KNN25_update_RF_500trees_prediction <- factor(fibroblasts.combined$KNN25_update_RF_500trees_prediction,
                                                                   levels = c("F-SH", "F-SL", "F-Trans", "F-WntX", "F-Act", "IR",
                                                                              "F-CI", "F-Cyc", "MYO", "MFC", "F-IFNS") )
Idents(fibroblasts.combined) <- fibroblasts.combined$KNN25_update_RF_500trees_prediction

col.set <- c("#fde800", "#fb8500", "#00b600", "#00896e", "#f60000", "#d647bf",
             "#00eefd", "#0053c8", "#0099cc", "#a7b5c9", "pink")
pops <- c("F-SH", "F-SL", "F-Trans", "F-WntX", "F-Act", "IR", "F-CI", "F-Cyc", "MYO", "MFC", "F-IFNS")
names(col.set) <- pops 
pl1 <- DimPlot(fibroblasts.combined, label = TRUE, cols = col.set, raster = FALSE)
pl1


pl1 <- DimPlot(fibroblasts.combined, label = TRUE, split.by = "Condition", ncol = 5, cols = col.set)
#pl1

ggData <- pl1$data
ggData$Condition <- fibroblasts.combined$Condition
ggData2 <- ggData
ggData2$Condition <- factor(ggData2$Condition, levels = c("Healthy", "MI-day 1", "MI-day 3", "MI-day 5", "IR",
                                                          "MI-day 7", "MI-day 14", "MI-day 28/30",  "AngII", "TAC"))
pl1 <- ggplot(ggData2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 15) +
  facet_wrap( ~Condition, ncol = 5 ) +  
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(strip.background = element_blank())
pl1

pl2 <- ggplot(ggData, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 15) +
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4)))

ggData %>%
  dplyr::group_by(ident) %>%
  summarize(UMAP_1 = median(x = UMAP_1), UMAP_2 = median(x = UMAP_2)) -> centers

pl2 <- pl2 + geom_text(data = centers, mapping = aes(label = ident), size = 4.5, colour="black", fontface="bold") +
  theme(text = element_text(size = 14, family="Helvetica")) + xlab("UMAP 1") + ylab("UMAP 2")
pl2

##############################################################################
### Update the MI cells to use the cell labels defined from the MI dataset ###
##############################################################################
load("All_datasets_cFIT_integration_Seurat.RData")
head(fibroblasts.combined@meta.data)

mi.meta.data <- read.csv("MI_fibroblasts_KNN25_meta_data.csv",
                         row.names = 1, stringsAsFactors = FALSE)

overlapping.cells <- intersect(colnames(fibroblasts.combined), rownames(mi.meta.data))

cluster.idents <- fibroblasts.combined$ClusterID
updated.idents <- cluster.idents

updated.idents[overlapping.cells] <- mi.meta.data[overlapping.cells, 'Clusters_updated_knn_basic_k25']

fibroblasts.combined$Clusters_MI_updated_KNN25 <- updated.idents
Idents(fibroblasts.combined) <- fibroblasts.combined$Clusters_MI_updated_KNN25

col.set <- c("#fde800", "#fb8500", "#00b600", "#00896e", "#f60000", "#d647bf",
             "#00eefd", "#0053c8", "#0099cc", "#a7b5c9", "pink")
pops <- c("F-SH", "F-SL", "F-Trans", "F-WntX", "F-Act", "IR", "F-CI", "F-Cyc", "MYO", "MFC", "F-IFNS")
names(col.set) <- pops 
pl1 <- DimPlot(fibroblasts.combined, label = TRUE, cols = col.set, raster = FALSE)
pl1

save(fibroblasts.combined, file = "All_datasets_cFIT_integration_Seurat.RData")

############################################################################
### Train a Random Forest classifier to detect fibroblast subpopulations ###
############################################################################
library(caret)
library(ranger)

table(fibroblasts.combined$Experiment)
fibroblasts.combined.mi <- subset(fibroblasts.combined, Experiment %in% c("GFPD3", "GFPD7", "FORTE", "tdTomD3",
                                                                          "CTHRC1_MID14", "CTHRC1_MID30"))

exp.data <- GetAssayData(fibroblasts.combined.mi, assay = "cFIT")

feature.matrix <- t(as.matrix(exp.data))
feature.matrix <- as.data.frame(feature.matrix)

class.labels <- Idents(fibroblasts.combined.mi)

set.seed(1)
folds <- createFolds(class.labels, k = 10, list = TRUE, returnTrain = FALSE)
#Iterate through the training/test folds

label.set <- class.labels
foldNum = 0
all.predictions.continuous <- c()
all.predictions <- c()
all.test.labels <- c()
for (thisFold in folds) {
  test <- feature.matrix[thisFold, ]
  train <- feature.matrix[-thisFold, ]
  
  train.labels <- label.set[-thisFold]
  
  rf.classifier <- ranger(x = train,
                          y = as.factor(train.labels),
                          num.trees = 500)
  print("Fibroblast classifier trained")
  
  #pred.test <- predict(rf.classifier, data.frame(test))
  pred.fibroblasts <- predict(rf.classifier, test)
  
  pred.test <- as.character(pred.fibroblasts$predictions)
  
  test.labels <- label.set[thisFold]
  
  all.predictions <- append(all.predictions, as.character(pred.test))
  
  all.test.labels <- c(all.test.labels, test.labels)
  foldNum = foldNum + 1
  print(paste0(foldNum, " folds completed"))
}
all.test.labels.update <- plyr::mapvalues(all.test.labels,
                                          from = seq(1,11), to  = levels(label.set))
cm <- confusionMatrix( data=as.factor(all.predictions), reference=as.factor(all.test.labels.update) )
cm
cm$byClass
class.accuracy.matrix <- cm$byClass
rownames(class.accuracy.matrix) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix))

## Plot the prediction accuracy derived from the confusion matrix
clusters <- unique(class.labels)
ggData <- data.frame(
  CellType = rep(rownames(class.accuracy.matrix), 4),
  Metric = c(rep("Precision", length(clusters)), rep("Recall", length(clusters)), 
             rep("Balanced Accuracy", length(clusters)), rep("Specificity", length(clusters))),
  Value = c(class.accuracy.matrix[, "Precision"], class.accuracy.matrix[, "Recall"],
            class.accuracy.matrix[, "Balanced Accuracy"], class.accuracy.matrix[, "Specificity"])
)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
pl
plot.label <- paste0("Forte_GFPD3_GFPD7_CTHRC1_Hif1a_varGenes4000_r15_cFIT_intExp_RF_500tree_CrossVal10_15x4.5inch.pdf")
ggsave(pl, filename = plot.label, width = 15, height = 4.5, units = "in", device = cairo_pdf)

file.label <- paste0("Forte_GFPD3_GFPD7_CTHRC1_Hif1a_varGenes4000_r15_cFIT_intExp_RF_500tree_CrossVal10_predictions_accuracy.csv")
write.csv(ggData, file = file.label, row.names = FALSE)

## Plot accuracy metrics with updated colour code
col.set <- c("#fde800", "#fb8500", "#00b600", "#00896e", "#f60000", "#d647bf",
             "#3bfff2", "#0053c8", "#0099cc", "#a7b5c9", "pink")
pops <- c("F-SH", "F-SL", "F-Trans", "F-WntX", "F-Act", "IR", "F-CI", "F-Cyc", "MYO", "MFC", "F-IFNS")
names(col.set) <- pops 
this.file <-paste0("Forte_GFPD3_GFPD7_CTHRC1_Hif1a_varGenes4000_r15_cFIT_intExp_RF_500tree_CrossVal10_predictions_accuracy.csv")
ggData <- read.csv(this.file)
pl <- ggplot(ggData, aes(x = reorder(CellType, -Value), y = Value, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap( ~Metric, nrow=1) + xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = col.set) + theme(legend.position = "none") +
  theme(strip.background = element_blank()) + labs(fill = "Cell type")
pl


############################################################################
### Train RF classifier on all data and make predictions to new datasets ###
############################################################################
library(ranger)

fibroblasts.combined.mi <- subset(fibroblasts.combined, Experiment %in% c("GFPD3", "GFPD7", "FORTE", "tdTomD3",
                                                                          "CTHRC1_MID14", "CTHRC1_MID30"))
exp.data <- GetAssayData(fibroblasts.combined.mi, assay = "cFIT")

feature.matrix <- t(as.matrix(exp.data))
feature.matrix <- as.data.frame(feature.matrix)

class.labels <- fibroblasts.combined.mi$Clusters_MI_updated_KNN25

set.seed(1)
start <- Sys.time()
rf.classifier.complete <- ranger(x = feature.matrix,
                                 y = as.factor(class.labels),
                                 num.trees = 500)
end <- Sys.time()
end - start

new.cell.idents <- fibroblasts.combined$Clusters_MI_updated_KNN25

### AngII cell predictions
cells.angII <- colnames(fibroblasts.combined)[which(fibroblasts.combined$Experiment == "AngII")]
angII.exp <- GetAssayData(fibroblasts.combined, assay = "cFIT")[, cells.angII]
angII.exp <- t(as.matrix(angII.exp)) %>% as.data.frame()

angII.rf.predictions <- predict(rf.classifier.complete, angII.exp)$predictions
table(angII.rf.predictions)
names(angII.rf.predictions) <- cells.angII

new.cell.idents[cells.angII] <- angII.rf.predictions

### TAC cell predictions
cells.tac <- colnames(fibroblasts.combined)[which(fibroblasts.combined$Experiment == "TAC")]
tac.exp <- GetAssayData(fibroblasts.combined, assay = "cFIT")[, cells.tac]
tac.exp <- t(as.matrix(tac.exp)) %>% as.data.frame()

tac.rf.predictions <- predict(rf.classifier.complete, tac.exp)$predictions
table(tac.rf.predictions)
names(tac.rf.predictions) <- cells.tac

new.cell.idents[cells.tac] <- tac.rf.predictions

### IR cell predictions
cells.ir <- colnames(fibroblasts.combined)[which(fibroblasts.combined$Experiment == "IR")]
ir.exp <- GetAssayData(fibroblasts.combined, assay = "cFIT")[, cells.ir]
ir.exp <- t(as.matrix(ir.exp)) %>% as.data.frame()

ir.rf.predictions <- predict(rf.classifier.complete, ir.exp)$predictions
table(ir.rf.predictions)
names(ir.rf.predictions) <- cells.ir

new.cell.idents[cells.ir] <- ir.rf.predictions

fibroblasts.combined$KNN25_update_RF_500trees_prediction <- new.cell.idents

### UMAP plot of updated clusters
Idents(fibroblasts.combined) <- fibroblasts.combined$KNN25_update_RF_500trees_prediction

### Plot the new identities
col.set <- c("#fde800", "#fb8500", "#00b600", "#00896e", "#f60000", "#d647bf",
             "#3bfff2", "#0053c8", "#0099cc", "#a7b5c9", "pink")
pops <- c("F-SH", "F-SL", "F-Trans", "F-WntX", "F-Act", "IR", "F-CI", "F-Cyc", "MYO", "MFC", "F-IFNS")
names(col.set) <- pops 
pl1 <- DimPlot(fibroblasts.combined, label = TRUE, cols = col.set, raster = FALSE)
pl1

pl1 <- DimPlot(fibroblasts.combined, label = TRUE, split.by = "Condition", ncol = 5, cols = col.set)

ggData <- pl1$data
ggData$Condition <- fibroblasts.combined$Condition
ggData2 <- ggData
ggData2$Condition <- factor(ggData2$Condition, levels = c("Healthy", "MI-day 1", "MI-day 3", "MI-day 5", "IR",
                                                          "MI-day 7", "MI-day 14", "MI-day 28/30",  "AngII", "TAC"))
pl1 <- ggplot(ggData2, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 15) +
  facet_wrap( ~Condition, ncol = 5 ) +  
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(strip.background = element_blank())
pl1

pl2 <- ggplot(ggData, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 15) +
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4)))

ggData %>%
  dplyr::group_by(ident) %>%
  summarize(UMAP_1 = median(x = UMAP_1), UMAP_2 = median(x = UMAP_2)) -> centers

pl2 <- pl2 + geom_text(data = centers, mapping = aes(label = ident), size = 4.5, colour="black", fontface="bold") +
  theme(text = element_text(size = 14, family="Helvetica")) + xlab("UMAP 1") + ylab("UMAP 2")
pl2

## Do a merged condition plot
pl3 <- DimPlot(fibroblasts.combined, label = TRUE, split.by = "Condition", ncol = 5, cols = col.set)

ggData <- pl3$data
ggData$Condition <- fibroblasts.combined$Condition
ggData3 <- ggData
ggData3$ConditionsMerged <- plyr::mapvalues(ggData3$Condition, 
                                            from = c("MI-day 1", "MI-day 3", "MI-day 5","MI-day 7",
                                                     "MI-day 14", "MI-day 28/30", "IR", "AngII", "TAC"),
                                            to = c("MI-day 1-7", "MI-day 1-7", "MI-day 1-7", "MI-day 1-7",
                                                   "MI-day 14-30", "MI-day 14-30", "I/R-day 5", "AngII-day 14",
                                                   "TAC-day 62"))

ggData3$ConditionsMerged <- factor(ggData3$ConditionsMerged, 
                                   levels = c("Healthy", "MI-day 1-7", "I/R-day 5", 
                                              "MI-day 14-30", "AngII-day 14", "TAC-day 62"))
table(ggData3$ConditionsMerged)
pl3 <- ggplot(ggData3, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(size = 0.1, aes(colour = ident)) + theme_classic(base_size = 15) +
  facet_wrap( ~ConditionsMerged, ncol = 3 ) +  
  scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") + 
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(strip.background = element_blank())
pl3

pl2 <- pl2 + theme(legend.position = "none")
pl2

plot_grid(pl2, pl3, ncol=2, align = "tb", rel_widths = c( 0.38, 0.62 ))