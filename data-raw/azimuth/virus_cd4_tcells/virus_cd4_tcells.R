library(Seurat)

# downloaded from https://doi.org/10.6084/m9.figshare.16592693.v1
scdata <- readRDS("data-raw/azimuth/virus_cd4_tcells/ref_LCMV_CD4_mouse_release_v1.rds")
DimPlot(scdata)

DefaultAssay(scdata) <- 'integrated'
scdata <- ScaleData(scdata, features = scdata@assays$integrated@var.features)
scdata <- RunPCA(scdata, features = scdata@assays$integrated@var.features, seed.use = 42)
scdata <- RenameAssays(scdata, 'integrated' = 'refAssay')

scdata@reductions$refDR <- scdata@reductions$pca

# re-run UMAP saving model
scdata@reductions$umap <- scdata@reductions$pca <- NULL
scdata <- RunUMAP(scdata, dims=1:20, seed.use = 123, reduction = "refDR", reduction.name = "refUMAP", return.model = TRUE)
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

qs::qsave(list(map = map), 'inst/extdata/mouse_virus_cd4_tcells.qs')

# save Seurat object
scdata@reductions$ica <- NULL
names(scdata@reductions) <- c('PCA', 'UMAP')
scdata <- SeuratObject::RenameAssays(scdata, 'refAssay' = 'integrated')
scdata$sample <- scdata$Sample

# necessary to preserve reference resolutions
scdata$predicted.celltype <- scdata$functional.cluster
scdata@misc$ref_name <- 'mouse_virus_cd4_tcells'
scdata@misc$resoln <- 'predicted.celltype'

# no counts: use library-size corrected, de-logged counts
# - used to get mitochondrial percents during import
scdata[['RNA']]@counts <- expm1(scdata[['RNA']]@data)

qs::qsave(scdata, "data-raw/azimuth/virus_cd4_tcells/virus_cd4_tcell_atlas.qs")
