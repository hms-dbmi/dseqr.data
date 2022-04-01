library(Seurat)
library(tibble)
library(dplyr)
library(harmony)
library(symphony)
source('data-raw/symphony/utils_seurat.R')

load('scMuscle_mm10_slim_v1-1.RData')
sdata <- scMuscle.slim.seurat
rm(scMuscle.slim.seurat)

sdata <- RunPCA(sdata)
sdata <- RunHarmony.Seurat(sdata, 'sample')

# find n.pcs
n.pcs <- npcs(sdata, reduction = "harmony")
n.pcs == sdata@reductions[['umap_harmony']]@misc$n.pcs.used

# currently, Seurat does not let you cache the umap model for future mapping
# therefore, please use this custom function to learn a saveable UMAP model
sdata[['umap']] <- RunUMAP2(
    Embeddings(sdata, 'harmony')[, 1:n.pcs],
    assay='RNA',
    verbose=FALSE,
    umap.method='uwot',
    return.model=TRUE)

DimPlot(sdata, reduction = 'umap', group.by = 'harmony_res.1.2_IDs', shuffle = TRUE)

# store cell_type
sdata$cell_type <- sdata$harmony_res.1.2_IDs

# remove what won't use
sdata@meta.data <- sdata@meta.data[, 'cell_type', drop = FALSE]


ref <- buildReferenceFromSeurat(
    sdata,
    verbose = TRUE,
    save_umap = TRUE,
    save_uwot_path = '/home/alex/Documents/Batcave/zaklab/dseqr/inst/extdata/scmuscle_uwot_model')

ref$normalization_method <- 'log(CP10k+1)'
ref$save_uwot_path <- NULL

# save reference
qs::qsave(ref, 'inst/extdata/scmuscle_reference.qs')

# store clusters as idents
Idents(sdata) <- 'harmony_res.1.2_IDs'

# necessary to preserve reference resolutions
sdata$predicted.celltype <- sdata$harmony_res.1.2_IDs
sdata@misc$ref_name <- 'scmuscle'
sdata@misc$resoln <- 'predicted.celltype'

# save Seurat object
qs::qsave(sdata, 'data-raw/symphony/scmuscle.qs')

