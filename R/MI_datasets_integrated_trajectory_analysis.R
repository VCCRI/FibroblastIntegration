library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(SingleCellExperiment)

load("MI_datasets_cFIT_integration_Seurat_intersect_genes.RData")
Idents(fibroblasts.combined) <- fibroblasts.combined$Clusters_updated_knn_basic_k25

col.set <- c("#fde800", "#fb8500", "#00b600", "#00896e", "#f60000", "#d647bf",
             "#3bfff2", "#0053c8", "#0099cc", "#a7b5c9", "pink")

pops <- c("F-SH", "F-SL", "F-Trans", "F-WntX", "F-Act", "IR", "F-CI", "F-Cyc", "MYO", "MFC", "F-IFNS")
names(col.set) <- pops 

DimPlot(fibroblasts.combined, cols = col.set, label = TRUE)

#####################################
### Monocle 3 pseudotime analysis ###
#####################################

fibroblast.cds <- as.cell_data_set(fibroblasts.combined)

list_cluster <- fibroblasts.combined$Clusters_updated_knn_k50
names(list_cluster) <- colnames(fibroblasts.combined)

fibroblast.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

#fibroblast.cds <- cluster_cells(fibroblast.cds)
graph.control <- list(rann.k = 30, maxiter = 15)
fibroblast.cds <- learn_graph(fibroblast.cds, learn_graph_control = graph.control, verbose = TRUE)

plot_cells(fibroblast.cds, label_groups_by_cluster = FALSE, 
           label_leaves = FALSE, label_branch_points = FALSE)


### Run on subsets of the data for Healthy
fibroblasts.subset <- subset(fibroblasts.combined, Condition %in% c("Healthy"))
DimPlot(fibroblasts.subset)

fibroblast.cds <- as.cell_data_set(fibroblasts.subset)

list_cluster <- fibroblasts.subset$Clusters_updated_knn_basic_k25
names(list_cluster) <- colnames(fibroblasts.subset)

fibroblast.cds <- cluster_cells(fibroblast.cds)
fibroblast.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
graph.control <- list(rann.k = 50, maxiter = 15, minimal_branch_len = 10)
fibroblast.cds <- learn_graph(fibroblast.cds, learn_graph_control = graph.control, verbose = TRUE)

pl0 <- plot_cells(fibroblast.cds, label_groups_by_cluster = FALSE, 
                  label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE,
                  trajectory_graph_segment_size = 1.25, cell_size = 0.5)
pl0 <- pl0 + scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") +
  theme(legend.position = "none")
pl0


### Run on subsets of the data for MI-days 1-3
fibroblasts.subset <- subset(fibroblasts.combined, Condition %in% c("MI-day 1", "MI-day 3"))
DimPlot(fibroblasts.subset)

fibroblast.cds <- as.cell_data_set(fibroblasts.subset)

list_cluster <- fibroblasts.subset$Clusters_updated_knn_basic_k25
names(list_cluster) <- colnames(fibroblasts.subset)

fibroblast.cds <- cluster_cells(fibroblast.cds)
fibroblast.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
graph.control <- list(rann.k = 50, maxiter = 15, minimal_branch_len = 10)
fibroblast.cds <- learn_graph(fibroblast.cds, learn_graph_control = graph.control, verbose = TRUE)

pl1 <- plot_cells(fibroblast.cds, label_groups_by_cluster = FALSE, 
                  label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE,
                  trajectory_graph_segment_size = 1.25, cell_size = 0.5)
pl1 <- pl1 + scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") +
  theme(legend.position = "none")
pl1

### Run on subsets of the data for MI-days 3-5
fibroblasts.subset <- subset(fibroblasts.combined, Condition %in% c("MI-day 3", "MI-day 5"))
DimPlot(fibroblasts.subset)

fibroblast.cds <- as.cell_data_set(fibroblasts.subset)

list_cluster <- fibroblasts.subset$Clusters_updated_knn_basic_k25
names(list_cluster) <- colnames(fibroblasts.subset)

fibroblast.cds <- cluster_cells(fibroblast.cds)
fibroblast.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
graph.control <- list(rann.k = 50, maxiter = 15, minimal_branch_len = 10)
fibroblast.cds <- learn_graph(fibroblast.cds, learn_graph_control = graph.control, verbose = TRUE)

pl2 <- plot_cells(fibroblast.cds, label_groups_by_cluster = FALSE, 
                  label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE,
                  trajectory_graph_segment_size = 1.25, cell_size = 0.5)
pl2 <- pl2 + scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") +
  theme(legend.position = "none")
pl2


### Run on subsets of the data for MI-days 5-7
fibroblasts.subset <- subset(fibroblasts.combined, Condition %in% c("MI-day 5", "MI-day 7"))
DimPlot(fibroblasts.subset)

fibroblast.cds <- as.cell_data_set(fibroblasts.subset)

list_cluster <- fibroblasts.subset$Clusters_updated_knn_basic_k25
names(list_cluster) <- colnames(fibroblasts.subset)

fibroblast.cds <- cluster_cells(fibroblast.cds)
fibroblast.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
graph.control <- list(rann.k = 50, maxiter = 15, minimal_branch_len = 10)
fibroblast.cds <- learn_graph(fibroblast.cds, learn_graph_control = graph.control, verbose = TRUE)

pl3 <- plot_cells(fibroblast.cds, label_groups_by_cluster = FALSE, 
                  label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE,
                  trajectory_graph_segment_size = 1.25, cell_size = 0.5)
pl3 <- pl3 + scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") +
  theme(legend.position = "none")
pl3


### Run on subsets of the data for MI-days 7-14
fibroblasts.subset <- subset(fibroblasts.combined, Condition %in% c("MI-day 7", "MI-day 14"))
DimPlot(fibroblasts.subset)

fibroblast.cds <- as.cell_data_set(fibroblasts.subset)

list_cluster <- fibroblasts.subset$Clusters_updated_knn_basic_k25
names(list_cluster) <- colnames(fibroblasts.subset)

fibroblast.cds <- cluster_cells(fibroblast.cds)
fibroblast.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
graph.control <- list(rann.k = 50, maxiter = 15, minimal_branch_len = 10)
fibroblast.cds <- learn_graph(fibroblast.cds, learn_graph_control = graph.control, verbose = TRUE)

pl4 <- plot_cells(fibroblast.cds, label_groups_by_cluster = FALSE, 
                  label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE,
                  trajectory_graph_segment_size = 1.25, cell_size = 0.5)
pl4 <- pl4 + scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") +
  theme(legend.position = "none")
pl4


### Run on subsets of the data for MI-days 14-28/30
fibroblasts.subset <- subset(fibroblasts.combined, Condition %in% c("MI-day 14", "MI-day 28", "MI-day30"))
DimPlot(fibroblasts.subset)

fibroblast.cds <- as.cell_data_set(fibroblasts.subset)

list_cluster <- fibroblasts.subset$Clusters_updated_knn_basic_k25
names(list_cluster) <- colnames(fibroblasts.subset)

fibroblast.cds <- cluster_cells(fibroblast.cds)
fibroblast.cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
graph.control <- list(rann.k = 50, maxiter = 15, minimal_branch_len = 10)
fibroblast.cds <- learn_graph(fibroblast.cds, learn_graph_control = graph.control, verbose = TRUE)

pl5 <- plot_cells(fibroblast.cds, label_groups_by_cluster = FALSE, 
                  label_leaves = FALSE, label_branch_points = FALSE, label_cell_groups = FALSE,
                  trajectory_graph_segment_size = 1.25, cell_size = 0.5)
pl5 <- pl5 + scale_color_manual(values=col.set, breaks=pops, labels=pops, name="Cluster") +
  theme(legend.position = "none")
pl5


### Combine the plots with cowplot

plots.combined <- cowplot::plot_grid(pl0, pl1, pl2, pl3, pl4, pl5,
                                     align = "hv", ncol = 2, 
                                     labels = c("Healthy", "MI-days 1-3", "MI-days 3-5",
                                                "MI-days 5-7", "MI-days 7-14", "MI-days 14-28/30"))
plots.combined

