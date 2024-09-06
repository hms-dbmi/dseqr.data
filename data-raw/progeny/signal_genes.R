# ligand-receptor database from multinichenetr
lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
lr_network = lr_network |>
    dplyr::rename(ligand = from, receptor = to) |>
    dplyr::distinct(ligand, receptor) |>
    dplyr::mutate(ligand = make.names(ligand), receptor = make.names(receptor))

# all ExpressionSets fro LINCS L1000
esets_dir <- '/mnt/12tb/Batcave/GEO/l1000/level1/6-limma/esets'
eset_fnames <- list.files(esets_dir)

# all experiments where a gene was knocked-down (sh), over-expressed (oe), or stimulated with ligand (lig)
gene_fnames <- grep('_sh_|_oe_|_lig_', eset_fnames, value = TRUE)

pert_genes <- gsub('^(.+?)_(oe|sh|lig)_.+?$', '\\1', gene_fnames)
pert_types <- gsub('^(.+?)_(oe|sh|lig)_.+?$', '\\2', gene_fnames)

# restrict to ligand or receptor (signal) perturbations
is.signal <- pert_genes %in% unlist(lr_network) | pert_types == 'lig'
signal.genes <- unique(pert_genes[is.signal])

saveRDS(signal.genes, 'data-raw/progeny/signal_genes.rds')
