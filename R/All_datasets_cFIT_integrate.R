library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(Matrix)
library(cowplot)
library(cFIT)

###############################
### Forte et al fibroblasts ###
###############################
load("Forte/fibroblasts_MI_timecourse.RData")
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

######################################################
### Create new Seurat objects of overlapping genes ###

overlapping.genes <- Reduce(union, list(gfp.day3.genes.expressed, gfp.day7.genes.expressed,
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

### GFP-day 7 data
count.data <- Read10X("GFP_day7/GFP_day7_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(gfp.day7), ]
dim(count.data)
gfp.day7 <- CreateSeuratObject(counts = count.data, project = "GFP-day7")
gfp.day7@meta.data$ClusterID <- gfp.day7.cell.idents
gfp.day7@meta.data$Experiment <- rep("GFPD7", ncol(gfp.day7))

### TdTom day 3 data
count.data <- Read10X("TdTom_day3/TdTom_day3_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(tdTom.day3), ]
dim(count.data)
tdTom.day3 <- CreateSeuratObject(counts = count.data, project = "tdTom-day3")
tdTom.day3@meta.data$ClusterID <- tdTom.day3.cell.idents
tdTom.day3@meta.data$Experiment <- rep("tdTomD3", ncol(tdTom.day3))

### Forte data
count.data <- Read10X("Forte/MI_timecourse_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(forte.fibroblasts), ]
dim(count.data)
forte.fibroblasts <- CreateSeuratObject(counts = count.data, project = "FORTE")
forte.fibroblasts@meta.data$ClusterID <- forte.cell.idents
forte.fibroblasts@meta.data$Experiment <- rep("FORTE", ncol(forte.fibroblasts))

### CTHRC1 MI-day 14 data
count.data <- Read10X("CTHRC1/MI_day14_CR6/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(mid14.fibroblasts), ]
dim(count.data)
mid14.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day14")
mid14.fibroblasts@meta.data$ClusterID <- mid14.fibroblasts.cell.idents
mid14.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID14", ncol(mid14.fibroblasts))

### CTHRC1 MI-day 30 data
count.data <- Read10X("CTHRC1/MI_day30_CR6/filtered_feature_bc_matrix")
count.data <- count.data[overlapping.genes, colnames(mid30.fibroblasts), ]
dim(count.data)
mid30.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day30")
mid30.fibroblasts@meta.data$ClusterID <- mid30.fibroblasts.cell.idents
mid30.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID30", ncol(mid30.fibroblasts))

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

### Merge the fibroblast data together
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



########################################
### Select the variable genes to use ###
########################################

data.list <- Seurat::SplitObject(fibroblasts.combined, split.by = "Experiment")
length(data.list)

data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 4000)
  return(x)
})

var.genes <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 4000)
length(var.genes)
### Regenerate the Seurat objects with variable genes to save memory

### GFP-day 3 data
count.data <- Read10X("GFP_day3/GFP_day3_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(gfp.day3), ]
dim(count.data)
gfp.day3 <- CreateSeuratObject(counts = count.data, project = "GFP-day3")
gfp.day3@meta.data$ClusterID <- gfp.day3.cell.idents
gfp.day3@meta.data$Experiment <- rep("GFPD3", ncol(gfp.day3))

### GFP-day 7 data
count.data <- Read10X("GFP_day7/GFP_day7_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(gfp.day7), ]
dim(count.data)
gfp.day7 <- CreateSeuratObject(counts = count.data, project = "GFP-day7")
gfp.day7@meta.data$ClusterID <- gfp.day7.cell.idents
gfp.day7@meta.data$Experiment <- rep("GFPD7", ncol(gfp.day7))

### TdTom day 3 data
count.data <- Read10X("TdTom_day3/TdTom_day3_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(tdTom.day3), ]
dim(count.data)
tdTom.day3 <- CreateSeuratObject(counts = count.data, project = "tdTom-day3")
tdTom.day3@meta.data$ClusterID <- tdTom.day3.cell.idents
tdTom.day3@meta.data$Experiment <- rep("tdTomD3", ncol(tdTom.day3))

### Forte data
count.data <- Read10X("Forte/MI_timecourse_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(forte.fibroblasts), ]
dim(count.data)
forte.fibroblasts <- CreateSeuratObject(counts = count.data, project = "FORTE")
forte.fibroblasts@meta.data$ClusterID <- forte.cell.idents
forte.fibroblasts@meta.data$Experiment <- rep("FORTE", ncol(forte.fibroblasts))

### CTHRC1 MI-day 14 data
count.data <- Read10X("CTHRC1/MI_day14_CR6/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(mid14.fibroblasts), ]
dim(count.data)
mid14.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day14")
mid14.fibroblasts@meta.data$ClusterID <- mid14.fibroblasts.cell.idents
mid14.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID14", ncol(mid14.fibroblasts))

### CTHRC1 MI-day 30 data
count.data <- Read10X("CTHRC1/MI_day30_CR6/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(mid30.fibroblasts), ]
dim(count.data)
mid30.fibroblasts <- CreateSeuratObject(counts = count.data, project = "CTHRC1_MI-day30")
mid30.fibroblasts@meta.data$ClusterID <- mid30.fibroblasts.cell.idents
mid30.fibroblasts@meta.data$Experiment <- rep("CTHRC1_MID30", ncol(mid30.fibroblasts))

### AngII treatment data 
count.data <- Read10X("AngII_Pinto/AngII_treat_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(angII.fibroblasts), ]
dim(count.data)
angII.fibroblasts <- CreateSeuratObject(counts = count.data, project = "AngII")
angII.fibroblasts@meta.data$Condition <- angII.fibroblasts.conditions
angII.fibroblasts@meta.data$ClusterID <- angII.fibroblasts.cell.idents
angII.fibroblasts@meta.data$Experiment <- rep("AngII", ncol(angII.fibroblasts))
remove(count.data)

### TAC data
count.data <- Read10X("TAC/TAC_vs_sham_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(tac.fibroblasts), ]
dim(count.data)
tac.fibroblasts <- CreateSeuratObject(counts = count.data, project = "TAC")
tac.fibroblasts@meta.data$Condition <- tac.fibroblasts.conditions
tac.fibroblasts@meta.data$ClusterID <- tac.fibroblasts.cell.idents
tac.fibroblasts@meta.data$Experiment <- rep("TAC", ncol(tac.fibroblasts))
remove(count.data)

### IR data
count.data <- Read10X("ischaemia_repurfusion/IR_CR6_no_norm/filtered_feature_bc_matrix")
count.data <- count.data[var.genes, colnames(ir.fibroblasts), ]
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


#fibroblasts.combined <- subset(fibroblasts.combined, idents = clusters.keep)
#table(fibroblasts.combined$ClusterID)

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

# extract the raw counts and metadata for data sets from 5 technologies
data.list = split_dataset_by_batch(X=t(as.matrix(fibroblasts.combined@assays$RNA@counts)),
                                   batch = fibroblasts.combined@meta.data$Experiment,
                                   labels = fibroblasts.combined@meta.data$ClusterID,
                                   metadata = fibroblasts.combined@meta.data,
                                   dataset.name = 'fibroblast:')

exprs.list = preprocess_for_integration(data.list$X.list, var.genes, scale.factor=10^4, scale=T, center=F)

save(exprs.list, file = "All_datasets_IR_2samples_cFIT_var4000_exprList.RData")

### Code for running cFIT analysis ###
### Use random sketching method
library(Seurat)
library(cFIT)
load("All_datasets_IR_2samples_cFIT_var4000_exprList.RData")

options(future.globals.maxSize= 3491289600)
int.out = CFITIntegrate_sketched(X.list=exprs.list,
                                 r=15, 
                                 max.niter=100,
                                 future.plan='multicore', 
                                 workers = 4,
                                 time.out = 60 * 4,
                                 seed=0, verbose=T)

save(int.out, file = "All_datasets_IR_2samples_cFIT_varGenes4000_r15_cFITSketched05.RData")

