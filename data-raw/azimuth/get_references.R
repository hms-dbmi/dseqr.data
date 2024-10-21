library(qs)
library(Azimuth)
library(Seurat)

refs <- c(human_pbmc = 'pbmcref',
          human_lung_v2 = 'lungref',
          human_motorcortex = 'humancortexref',
          mouse_motorcortex = 'mousecortexref',
          human_bonemarrow = 'bonemarrowref',
          human_fetus = 'fetusref',
          human_liver = 'bonemarrowref')

options(timeout = 1000)

for (i in 1:length(refs)) {
    ref <- refs[i]
    SeuratData::InstallData(ref, force.reinstall = TRUE)
    ref <- SeuratData::LoadData(ref, type = 'azimuth')
    ref_path <- file.path('inst/extdata', paste0(names(refs)[i], '.qs'))
    qsave(ref, ref_path)
}

# UPLOAD ADDED REFS TO AWS
# ADD TO dl_data files argument
