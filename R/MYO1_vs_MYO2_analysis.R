library(Seurat)
library(ranger)
library(caret)
library(ggplot2)
library(patchwork)
library(magrittr)
library(SingleCellExperiment)
library(heatmap3)
library(dplyr)

load("MI_datasets_cFIT_integration_Seurat_intersect_genes.RData")

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

#myo1.cells <- colnames(gfp.day7)[which(gfp.day7.cell.idents == "MYO-1")]
#myo1.cells <- paste0("GFPD7_", myo1.cells)

#myo2.cells <- colnames(gfp.day7)[which(gfp.day7.cell.idents == "MYO-2")]
#myo2.cells <- paste0("GFPD7_", myo2.cells)

gfp.day7.integrated <- subset(fibroblasts.combined, Experiment == "GFPD7")

gfp.day7.idents <- Idents(gfp.day7)
names(gfp.day7.idents) <- paste0("GFPD7_", names(gfp.day7.idents))

gfp.day7.idents <- gfp.day7.idents[colnames(gfp.day7.integrated)]

gfp.day7.integrated <- AddMetaData(gfp.day7.integrated, metadata = gfp.day7.idents, col.name = "GFPD7_Idents")
Idents(gfp.day7.integrated) <- gfp.day7.integrated$GFPD7_Idents

DimPlot(gfp.day7.integrated, split.by = "Condition", label = TRUE)

gfp.day7.myo.subsets <- subset(gfp.day7.integrated, idents = c("MYO-1", "MYO-2"))
DimPlot(gfp.day7.myo.subsets, label = TRUE)

########################################################################################
### Train RF classifier on MYO-1 vs MYO-2 cells and make predictions to new datasets ###
########################################################################################

### First do some cross-validation to check if we can classify accurately
### Train a Random Forest classifier to detect fibroblast subpopulations 

exp.data <- GetAssayData(gfp.day7.myo.subsets, assay = "cFIT")

feature.matrix <- t(as.matrix(exp.data))
feature.matrix <- as.data.frame(feature.matrix)

class.labels <- gfp.day7.myo.subsets$GFPD7_Idents
class.labels <- droplevels(class.labels)
table(class.labels)

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
                                          from = c(1,2), to  = levels(label.set))
cm <- confusionMatrix( data=as.factor(all.predictions), reference=as.factor(all.test.labels.update) )
cm
cm$byClass
class.accuracy.matrix <- cm$byClass
#rownames(class.accuracy.matrix) <- sub("Class: (.*)", "\\1", rownames(class.accuracy.matrix))

## Plot the prediction accuracy derived from the confusion matrix
clusters <- unique(class.labels)
ggData <- data.frame(
  Metric = c("Precision", "Recall", "Balanced Accuracy", "Specificity"),
  Value = c(class.accuracy.matrix["Precision"], class.accuracy.matrix["Recall"],
            class.accuracy.matrix["Balanced Accuracy"], class.accuracy.matrix["Specificity"])
)
pl <- ggplot(ggData, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("Cell lineage") +
  theme_classic(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(strip.background = element_blank()) + ggtitle("MYO-1 vs MYO-2 RF accuracy")
pl
#plot.label <- paste0("transfer_predictions/Forte_GFPD3_GFPD7_CTHRC1_Hif1a_varGenes4000_r15_cFIT_intExp_RF_500tree_CrossVal10_15x4.5inch.pdf")
#ggsave(pl, filename = plot.label, width = 15, height = 4.5, units = "in", device = cairo_pdf)


### Now train a RF classifier on all data and make predictions
exp.data <- GetAssayData(gfp.day7.myo.subsets, assay = "cFIT")

feature.matrix <- t(as.matrix(exp.data))
feature.matrix <- as.data.frame(feature.matrix)

class.labels <- gfp.day7.myo.subsets$GFPD7_Idents
class.labels <- droplevels(class.labels)
table(class.labels)

set.seed(1)
start <- Sys.time()
rf.classifier.complete <- ranger(x = feature.matrix,
                                 y = as.factor(class.labels),
                                 num.trees = 500)
end <- Sys.time()
end - start

### Subset for the MYO/MFC subsets at days 5 - 28/30
myo.integrated <- subset(fibroblasts.combined, idents = c("MYO", "MFC"))
myo.integrated <- subset(myo.integrated, Condition %in% c("MI-day 5", "MI-day 7", "MI-day 14", "MI-day 28", "MI-day 30"))
myo.integrated$Condition <- factor(myo.integrated$Condition, levels = c("MI-day 5", "MI-day 7", 
                                                                        "MI-day 14", "MI-day 28", "MI-day 30"))
myo.integrated$Condition <- plyr::mapvalues(myo.integrated$Condition, 
                                            from = c( "MI-day 28", "MI-day 30" ),
                                            to = c("MI-day 28/30", "MI-day 28/30"))
myo.integrated$Condition <- factor(myo.integrated$Condition, 
                                   levels = c("MI-day 5", "MI-day 7", "MI-day 14", "MI-day 28/30"))

DimPlot(myo.integrated, label = TRUE)

new.cell.idents <- myo.integrated$Clusters_updated_knn_basic_k25
new.cell.idents <- droplevels(new.cell.idents)
table(new.cell.idents)

### Make predictions for MYO-1 vs MYO-2 
myo.exp <- GetAssayData(myo.integrated, assay = "cFIT")
myo.exp <- t(myo.exp) %>% as.data.frame()

myo.rf.predictions <- predict(rf.classifier.complete, myo.exp)$predictions
table(myo.rf.predictions)
names(myo.rf.predictions) <- colnames(myo.exp)

new.cell.idents <- myo.rf.predictions

### Update the myofibroblast subset with RF predictions
myo.integrated$RF_predictions <- myo.rf.predictions
Idents(myo.integrated) <- myo.integrated$RF_predictions

DimPlot(myo.integrated, label = TRUE)

DimPlot(myo.integrated, label = TRUE, split.by = "Condition", ncol = 4)


## Do a bar plot of MYO-1 vs MYO-2 cells per 

cell.data <- data.frame(Condition = myo.integrated$Condition,
                        Prediction = myo.integrated$RF_predictions)

cell.data %>%
  group_by(Condition, Prediction) %>%
  summarise(n = n()) %>%
  group_by(Condition) %>% 
  mutate(nGroup = sum(n), Percentage = n / sum(n) * 100) -> ggData


ggplot(ggData, aes(x = Prediction, y = Percentage, fill = Prediction)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  facet_wrap( ~Condition, nrow = 1 ) + theme_classic(base_size = 16) +
  theme(strip.background = element_blank()) + scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  xlab("Prediction") + ggtitle("MYO-1 vs MYO-2 predictions")

## Plot the MYO-1 and MYO-2 predictions back on the main UMAP
all.cell.idents <- rep("Other", ncol(fibroblasts.combined))
names(all.cell.idents) <- colnames(fibroblasts.combined)
all.cell.idents[colnames(myo.integrated)] <- as.character(Idents(myo.integrated))
all.cell.idents <- as.factor(all.cell.idents)

fibroblasts.combined <- AddMetaData(fibroblasts.combined, all.cell.idents, col.name = "MYO_subtype_RF_predictions")

col.set <- RColorBrewer::brewer.pal(3, name = "Dark2")
col.set[3] <- "#d9d9d9"

pl <- DimPlot(fibroblasts.combined, group.by = "MYO_subtype_RF_predictions", cols = col.set, combine = FALSE)
pl <- pl[[1]]

pl$data$Condition <- fibroblasts.combined$Condition
pl$data <- subset(pl$data, Condition %in% c("MI-day 5", "MI-day 7", "MI-day 14", "MI-day 28", "MI-day 30"))

pl$data$Condition <- factor(pl$data$Condition, levels = c("MI-day 5", "MI-day 7", 
                                                          "MI-day 14", "MI-day 28", "MI-day 30"))
pl$data$Condition <- plyr::mapvalues(pl$data$Condition, 
                                     from = c( "MI-day 28", "MI-day 30" ),
                                     to = c("MI-day 28/30", "MI-day 28/30"))
pl$data$Condition <- factor(pl$data$Condition, levels = c("MI-day 5", "MI-day 7", "MI-day 14", "MI-day 28/30"))

pl + facet_wrap( ~Condition, nrow = 1 ) +
  theme(strip.background = element_blank()) +
  xlab("UMAP 1") + ylab("UMAP 2")



