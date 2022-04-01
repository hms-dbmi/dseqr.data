### EXPLORATORY ANALYSIS OF SINGLE-CELL RNA-SEQ DATA ###
# Author: Eddie Cano-Gamez (ecg@sanger.ac.uk)
#
# This code performs normalization, regression of covariates and scaling of scRNA-seq data. Next, it performs dimensionality reduction and visualisation.

# Adapted to generate Seurat reference
# data downloaded from https://www.opentargets.org/projects/effectorness

# LOADING LIBRARIES
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(org.Hs.eg.db)

data_dir <- 'data-raw/azimuth/tcell_effectorness/data'

# recover ensemble gene ids
# genes <- read.table(file.path(data_dir, "NCOMMS-19-7936188_scRNAseq_genes.tsv"))[[1]]
# enids <- mapIds(org.Hs.eg.db,
#                 keys=genes,
#                 column="ENSEMBL",
#                 keytype="SYMBOL",
#                 multiVals="first")

# write.table(data.frame(enids, genes), file.path(data_dir, "features.tsv"), col.names = FALSE, row.names = FALSE, sep = '\t')
# R.utils::gzip(file.path(data_dir, "features.tsv"))


# LOADING DATA
## Loading expression matrix from a directory containing the matrix, genes and barcodes file
expressionMatrix <- Read10X(data_dir)

## Loading cell annotations
metadata <- read.table(file.path(data_dir, "NCOMMS-19-7936188_metadata.txt"), header=T, sep="\t", check.names = FALSE)

# CREATING SEURAT OBJECT
cells <- CreateSeuratObject(expressionMatrix, min.cells = 3, min.genes = 200, project = "Cytoimmgen_10X_optimization")

## Adding cell annotations as metadata
cells@meta.data <- metadata

# NORMALISING GENE EXPRESSION
cells <- NormalizeData(cells, normalization.method = "LogNormalize")

# REGRESSING COVARAITES
cells <- ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch.10X"))

# FINDING HIGHLY VARIABLE GENES
cells <- FindVariableFeatures(cells, mean.cutoff = c(0.1, 5), dispersion.cutoff = c(2, Inf))

# REDUCING DIMENSIONS
cells <- SeuratObject::RenameAssays(cells, 'RNA' = 'refAssay')
cells <- RunPCA(cells, npcs = 50, reduction.name = 'refDR')

# UMAP
cells <- RunUMAP(cells, dims=1:30, reduction = "refDR", reduction.name = "refUMAP", return.model = TRUE)

# VISUALISING DATA
DimPlot(cells, reduction = "refUMAP", group.by = "cell.type", pt.size = 0.5)
DimPlot(cells, reduction = "refUMAP", group.by = "donor.id", pt.size = 0.5)
DimPlot(cells, reduction = "refUMAP", group.by = "Phase", pt.size = 0.5)
DimPlot(cells, reduction = "refUMAP", group.by = "cytokine.condition", pt.size = 0.5)

cells$celltype.cytokines <- paste(cells$cell.type, cells$cytokine.condition, sep='.')
cells$celltype.stimulated <- paste(cells$cell.type, ifelse(cells$cytokine.condition == 'UNS', 'UNS', 'STIM'), sep='.')
DimPlot(cells, reduction = "refUMAP", group.by = "celltype.cytokines", pt.size = 0.5)


# compute first 50 neighbors in PCA space of reference and cache annoy index
cells <- FindNeighbors(
  object = cells,
  reduction = "refDR",
  dims = 1:50,
  graph.name = "refdr.annoy.neighbors",
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)


# save in same format as Azimuth references
ref <- SeuratObject::CreateAssayObject(data = cells[['refAssay']]@data)

map <- SeuratObject::CreateSeuratObject(ref, assay = 'refAssay')
map@reductions$refDR <- cells@reductions$refDR
map@reductions$refUMAP <- cells@reductions$refUMAP

map[['refdr.annoy.neighbors']] <- cells[['refdr.annoy.neighbors']]
map$celltype.cytokines <- cells$celltype.cytokines
map$celltype.stimulated <- cells$celltype.stimulated

qs::qsave(list(map = map), 'inst/extdata/human_differentiated_tcell.qs')

# save Seurat object
scdata <- cells
names(scdata@reductions) <- c('PCA', 'umap')
Idents(scdata) <- 'celltype.cytokines'
scdata <- SeuratObject::RenameAssays(scdata, 'refAssay' = 'RNA')
scdata$sample <- paste(scdata$donor.id, scdata$cytokine.condition, sep='.')

# necessary to preserve reference resolutions
scdata$predicted.celltype.cytokines <- scdata$celltype.cytokines
scdata$predicted.celltype.stimulated <- scdata$celltype.stimulated
scdata@misc$ref_name <- 'human_differentiated_tcell'
scdata@misc$resoln <- 'predicted.celltype.cytokines'

qs::qsave(scdata, file.path(data_dir, "differentiated_cd4_tcells.qs"))
