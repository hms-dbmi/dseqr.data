# based on saezlab/footprints and fork jperales/footprints
library(dplyr)


# load dprimes from LINCS L1000
dprimes <- readRDS('/mnt/12tb/Batcave/GEO/l1000/level1/6-limma/l1000_es.rds')
deltap <- readRDS('/mnt/12tb/Batcave/GEO/l1000/level1/6-limma/l1000_es_deltap.rds')

common_genes <- intersect(row.names(dprimes), row.names(deltap))

# all experiments where a gene was knocked-down (sh), over-expressed (oe), or stimulated with ligand (lig)
pert_names1 <- grep('_sh_|_oe_|_lig_', colnames(dprimes), value = TRUE)
pert_names2 <- grep('_sh_|_oe_|_lig_', colnames(deltap), value = TRUE)

pert_genes <- gsub('^(.+?)_(oe|sh|lig)_.+?$', '\\1', pert_fnames)
pert_types <- gsub('^(.+?)_(oe|sh|lig)_.+?$', '\\2', pert_fnames)


get_record <- function(eset, eset_fname, pert_gene, pert_type) {
    pdata <- Biobase::pData(eset)
    is.ctl <- pdata$drug == 'ctl'

    effect <- ifelse(pert_type %in% c('lig', 'oe'), 'activating', 'inhibiting')

    # pathway is pert_gene
    # TODO: explore correlations between pert_gene and any paired ligand/receptor

    list(
        id = gsub('[.]rds$', '', eset_fname),
        control = row.names(pdata)[is.ctl],
        perturbed = row.names(pdata)[!is.ctl],
        pathway = pert_gene,
        effect = effect
    )
}

# get exprs and records for all signal perts
data <- list(
    expr = list(),
    records = list()
)

for (i in 1:length(pert_fnames)) {
    cat('Working on', i, 'of', length(pert_fnames), '...\n')
    pert_fname <- pert_fnames[i]
    pert_gene <- pert_genes[i]
    pert_type <- pert_types[i]

    eset_fpath <- file.path(esets_dir, pert_fname)
    eset <- readRDS(eset_fpath)

    emat <- Biobase::exprs(eset)
    rec <- get_record(eset, pert_fname, pert_gene, pert_type)

    n.ctrl <- length(rec$control)
    n.test <- length(rec$perturbed)
    n.tot <- n.ctrl + n.test

    if (n.ctrl > 0 & n.test > 0 & n.tot > 2) {
        data$expr[[rec$id]] <- emat
        data$records[[rec$id]] <- rec
    }
}

n.removed <- length(pert_fnames) - length(data$expr)
cat(n.removed, 'of', length(pert_fnames), "removed because no control/test sample or fewer than 3 samples total")

# get zscores for all signal perts
zscores <- data2zscores(data)

saveRDS(zscores, 'data-raw/progeny/zscores.rds')


