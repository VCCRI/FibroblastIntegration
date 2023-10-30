
library(Seurat)
library(ggplot2)
library(dplyr)
library(ViSEAGO)
library(biomaRt)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GO.db)

load("MI_datasets_cFIT_integration_Seurat_intersect_genes.RData")

Idents(fibroblasts.combined) <- fibroblasts.combined$Clusters_updated_knn_basic_k25
DefaultAssay(fibroblasts.combined) <- "RNA"

####################################################
### Do GO term test on subsets of MI time points ###
####################################################

## Define BiomaRt for gene ID mapping
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
head(datasets)

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

all.attributes = listAttributes(ensembl)
#uniprotswissprot

### Read in background genes
background <- rownames(fibroblasts.combined)

attributes=c( 'mgi_symbol','uniprotswissprot' )
bg.genes <- background
mart.list <- getBM(attributes=attributes, filters="mgi_symbol",values=bg.genes,
                   mart=mart, uniqueRows=T)

length(setdiff(bg.genes, mart.list$mgi_symbol))
mart.list <- subset(mart.list, uniprotswissprot != '')
head(mart.list)
bg.uniprot <- unique(mart.list$uniprotswissprot)
length(bg.uniprot)

### Define the ViSEAGO objects
Uniprot<-ViSEAGO::Uniprot2GO()
myGENE2GO<-ViSEAGO::annotate(
  "mouse",
  Uniprot
)

cells.1 <- rownames(subset(fibroblasts.combined@meta.data, Clusters_updated_knn_basic_k25 == "F-Cyc" &
                             Condition %in% c("MI-day 1", "MI-day 3")))
cells.2 <- rownames(subset(fibroblasts.combined@meta.data, 
                           Clusters_updated_knn_basic_k25 %in% c("F-CI", "IR", "F-Act") &
                             Condition %in% c("MI-day 1", "MI-day 3")))

deg.label <- "Days1,3_F-Cyc_vs_F-CI,IR,F-Act"
res.table <- FindMarkers(fibroblasts.combined, 
                         ident.1 = cells.1, 
                         ident.2 = cells.2, 
                         min.pct = 0.25,
                         logfc.threshold = 0.5, 
                         only.pos = TRUE, assay = "RNA", test.use = "MAST")

head(res.table)
res.table <- subset(res.table, p_val_adj < 1e-05 & avg_log2FC > 0.5)
nrow(res.table)
selection <- rownames(res.table)

mart.list <- getBM(attributes=attributes, filters="mgi_symbol",values=selection,
                   mart=mart, uniqueRows=T)

length(setdiff(selection, mart.list$mgi_symbol))
mart.list <- subset(mart.list, uniprotswissprot != '')
head(mart.list)
selection.uniprot <- unique(mart.list$uniprotswissprot)
length(selection.uniprot)

# create topGOdata for BP
BP<-ViSEAGO::create_topGOdata(
  geneSel=selection.uniprot,
  allGenes=bg.uniprot,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)

classic<-topGO::runTest(
  BP,
  algorithm ="classic",
  statistic = "fisher"
)

allRes <- GenTable(BP, classicFisher = classic, orderBy = "classicFisher", 
                   ranksOf = "classicFisher", topNodes = length(classic@score), numChar = 100)
all.scores <- classic@score
all.scores <- all.scores[allRes$GO.ID]
allRes$classicFisher <- all.scores
allRes$Fisher_p_adj <- p.adjust(allRes$classicFisher, method = "BH")
head(allRes, n = 15)

resSig <- subset(allRes, Fisher_p_adj < 0.05)
write.csv(resSig, paste0("Results_MI_integration/GO_terms/", deg.label, "_TopGO_classic_Fisher_padj05.csv"))

godata <- resSig[1:20, ]
go_plot <- ggplot(godata, aes(x=reorder(Term, -log10(Fisher_p_adj)), y=-log10(Fisher_p_adj))) + 
  geom_bar(stat="identity", fill = "grey", colour = "black") + coord_flip() + ylab("-Log10(adjusted p-value)") + 
  theme_bw(base_size=22) + xlab("") + ggtitle(deg.label) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(label=Significant), hjust=-0.15, size = 5) + 
  scale_y_continuous(limits = c(0, max(-log10(godata$Fisher_p_adj)) + 0.5))
print(go_plot)

