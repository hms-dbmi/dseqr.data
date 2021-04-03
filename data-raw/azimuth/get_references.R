library(qs)
ref_name <- 'human_pbmc'
version <- 'v1.0.0'
ref_dir <- 'inst/extdata'

reference <- LoadReference(path = file.path("https://seurat.nygenome.org/azimuth/references",
                                            version,
                                            ref_name),
                           seconds=200L)

qsave(reference, file.path(ref_dir, paste0(ref_name, '.qs')))
