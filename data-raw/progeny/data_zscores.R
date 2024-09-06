# based on saezlab/footprints and fork jperales/footprints
library(dplyr)

#' Compute z-scores for one contrast
#'
#' @param rec   The experiment record (`data.frame`, from the yaml files)
#' @param emat  The expression matrix [genes x experiments]
expr2zscore = function(rec, emat) {
    message(rec$id)

    # get expression values from source name
    control = emat[,rec$control, drop=FALSE]
    perturbed = emat[,rec$perturbed, drop=FALSE]

    # build speed models
    mean_control= apply(control, 1, mean)
    sd_control = apply(control, 1, sd)
    mean_perturbed = apply(perturbed, 1, mean)
    logFC = mean_perturbed - mean_control
    model = loess(sd_control ~ mean_control)

    # NOTE: loess won't predict outside of its range => NA
    zscore <- logFC / predict(model, mean_perturbed)

    return(zscore)
}

#' Calculate Z-scores for all experiments
#'
#' @param data  A list with elements `records` and `expr`
#' @return      The `zscores` and `index` objects
data2zscores = function(data) {
    records = data$records
    expr = data$expr

    zscores = mapply(expr2zscore, rec=records, emat=expr, SIMPLIFY=FALSE) %>%
        narray::stack(along=2)

    idx_remove = c("control", "perturbed")
    sign_lookup = setNames(c(1,-1), c("activating", "inhibiting"))
    index = lapply(records, function(x) x[setdiff(names(x), idx_remove)]) %>%
        do.call(bind_rows, .) %>%
        mutate(sign = sapply(effect, function(x) sign_lookup[x]))

    stopifnot(colnames(zscores) == index$id)

    list(zscores=zscores, index=index)
}

# all ExpressionSets fro LINCS L1000
esets_dir <- '/mnt/12tb/Batcave/GEO/l1000/level1/6-limma/esets'
eset_fnames <- list.files(esets_dir)

# all experiments where a gene was knocked-down (sh), over-expressed (oe), or stimulated with ligand (lig)
pert_fnames <- grep('_sh_|_oe_|_lig_', eset_fnames, value = TRUE)

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


