# downloaded from https://github.com/carmonalab/ProjecTILs/issues/23
library(Seurat)
scdata <- readRDS("data-raw/azimuth/til_tcells/ref_TILAtlas_mouse_wcounts.rds")

scdata <- RenameAssays(scdata, 'integrated' = 'refAssay')
scdata@reductions$refDR <- scdata@reductions$pca
DimPlot(scdata)

# re-run UMAP saving model
scdata@reductions$umap <- scdata@reductions$pca <- NULL
scdata <- RunUMAP(scdata, dims=1:15, seed.use = 123, n.neighbors = 30, min.dist = 0.3, reduction = "refDR", reduction.name = "refUMAP", return.model = TRUE)
DimPlot(scdata, reduction = "refUMAP")

scdata <- FindNeighbors(
    object = scdata,
    reduction = "refDR",
    dims = 1:50,
    graph.name = "refdr.annoy.neighbors",
    k.param = 50,
    cache.index = TRUE,
    return.neighbor = TRUE,
    l2.norm = TRUE
)

# save in same format as Azimuth references
ref <- SeuratObject::CreateAssayObject(data = scdata[['refAssay']]@data)

map <- SeuratObject::CreateSeuratObject(ref, assay = 'refAssay')
map@reductions$refDR <- scdata@reductions$refDR
map@reductions$refUMAP <- scdata@reductions$refUMAP

map[['refdr.annoy.neighbors']] <- scdata[['refdr.annoy.neighbors']]
map$celltype <- scdata$functional.cluster

qs::qsave(list(map = map), 'inst/extdata/mouse_til_tcells.qs')

# save Seurat object
names(scdata@reductions) <- c('PCA', 'UMAP')
scdata <- SeuratObject::RenameAssays(scdata, 'refAssay' = 'integrated')
scdata$sample <- scdata$Sample
scdata$sample[is.na(scdata$sample)] <- 'NA'

# necessary to preserve reference resolutions
scdata$predicted.celltype <- scdata$functional.cluster
scdata@misc$ref_name <- 'mouse_til_tcells'
scdata@misc$resoln <- 'predicted.celltype'


qs::qsave(scdata, "data-raw/azimuth/til_tcells/til_tcell_atlas.qs")
