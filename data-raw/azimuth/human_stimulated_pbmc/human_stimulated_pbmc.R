library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)

# downloaded data from GSE181897
data_dir <- 'data-raw/azimuth/human_stimulated_pbmc/data'

# re-run PCA on scaled data (need loadings etc in reference)
scdata.scaled <- zellkonverter::readH5AD(file.path(data_dir, 'GSE181897_concat.4.h5ad'))
pca <- RunPCA(
    scdata.scaled@assays@data@listData$X,
    reduction.key = 'refdr_',
    reduction.name = 'refDR',
    assay = 'refAssay')

rm(scdata.scaled); gc()

scdata <- zellkonverter::readH5AD(file.path(data_dir, 'GSE181897_concat.4.raw.h5ad'))
sdata <- Seurat::CreateSeuratObject(counts = scdata@assays@data@listData$X)

scater::plotReducedDim(scdata, 'X_umap', colour_by = 'cond')

# adding cell annotations as metadata
sdata@meta.data <- as.data.frame(scdata@colData)

# rename
sdata <- SeuratObject::RenameAssays(sdata, 'RNA' = 'refAssay')

# add PCA (calculated above)
sdata@reductions$refDR <- pca
sdata <- RunUMAP(sdata,
                 dims = 1:30,
                 reduction = "refDR",
                 reduction.name = "refUMAP",
                 return.model = TRUE)

DimPlot(sdata, reduction = "refUMAP", group.by = 'cond', pt.size = .5, cols = scdata@metadata$cond_colors)
DimPlot(sdata, reduction = "refUMAP", group.by = 'ct1', pt.size = .5, cols = scdata@metadata$ct1_colors)
DimPlot(sdata, reduction = "refUMAP", group.by = 'ct2', pt.size = .5, cols = scdata@metadata$ct2_colors)
DimPlot(sdata, reduction = "refUMAP", group.by = 'ct3', pt.size = .5, cols = scdata@metadata$ct3_colors)

rm(scdata); gc()

# format annotation
# A = TNF-alpha
# B = IFN-beta
# C = control?
# G = IFN-gamma
# P = PMA/Ionomycin (TH1)
# R = Resiquimod
# 0 = ? (only 455 cells)
levels(sdata$cond)
sdata$condition <- sdata$cond
levels(sdata$condition) <- c('0', 'TNF-alpha', 'IFN-beta', 'Control', 'IFN-gamma', 'PMA/Ionomycin', 'Resiquimod')
sdata$celltype.l1 <- sdata$ct1
sdata$celltype.l2 <- sdata$ct2
sdata$celltype.l3 <- sdata$ct3

# compute first 50 neighbors in PCA space of reference and cache annoy index
sdata <- FindNeighbors(
    object = sdata,
    reduction = "refDR",
    dims = 1:50,
    graph.name = "refdr.annoy.neighbors",
    k.param = 50,
    cache.index = TRUE,
    return.neighbor = TRUE,
    l2.norm = TRUE
)


# save in same format as Azimuth references
ref <- SeuratObject::CreateAssayObject(data = sdata[['refAssay']]@data)

map <- SeuratObject::CreateSeuratObject(ref, assay = 'refAssay')
map@reductions$refDR <- sdata@reductions$refDR
map@reductions$refUMAP <- sdata@reductions$refUMAP

map[['refdr.annoy.neighbors']] <- sdata[['refdr.annoy.neighbors']]
map$condition <- sdata$condition
map$celltype.l1 <- sdata$celltype.l1
map$celltype.l2 <- sdata$celltype.l2
map$celltype.l3 <- sdata$celltype.l3

# also upload to S3
qs::qsave(list(map = map), 'inst/extdata/human_stimulated_pbmc.qs')

# save Seurat object for import
scdata <- sdata
names(scdata@reductions) <- c('PCA', 'umap')
Idents(scdata) <- 'condition'
scdata <- SeuratObject::RenameAssays(scdata, 'refAssay' = 'RNA')
scdata$sample <- paste0(scdata$condition, '.D', scdata$exp_id)

# necessary to preserve reference resolutions
scdata$predicted.condition <- scdata$condition
scdata$predicted.celltype.l1 <- scdata$celltype.l1
scdata$predicted.celltype.l2 <- scdata$celltype.l2
scdata$predicted.celltype.l3 <- scdata$celltype.l3
scdata@misc$ref_name <- 'human_stimulated_pbmc'
scdata@misc$resoln <- 'predicted.condition'

# seurat standard names
scdata$percent.mt <- scdata$percent_mito

qs::qsave(scdata, file.path(data_dir, "human_stimulated_pbmc.qs"))

