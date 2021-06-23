library(qs)
library(Azimuth)

get_reference <- function(ref_name, version = 'v1.0.0', ref_dir = 'inst/extdata') {
    reference <- LoadReference(path = file.path("https://seurat.nygenome.org/azimuth/references",
                                                version,
                                                ref_name),
                               seconds=200L)

    qsave(reference, file.path(ref_dir, paste0(ref_name, '.qs')))
    return(reference)
}

# human pbmc
# ref <- get_reference('human_pbmc')

#human lung
ref <- get_reference('human_lung')

# human/mouse motor cortex
ref <- get_reference('human_motorcortex')
ref <- get_reference('mouse_motorcortex')


# UPLOAD ADDED REFS TO AWS
# ADD TO dl_data files argument
