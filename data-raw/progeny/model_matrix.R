# point of this file:
# - use the zscores to create a linear model
library(dplyr)
library(broom)
library(Biobase)

min_mask = function(x, N=2) {
    if (N > length(x))
        !is.na(x)
    else
        seq_along(x) %in% order(x, decreasing=FALSE, na.last=NA)[1:N]
}

# throws out cov.unscaled as too much memory use
summary.mlm <- function (object, ny = ncol(coef), names = TRUE, ...)
{
    coef <- coef(object)
    effects <- object$effects
    resid <- object$residuals
    fitted <- object$fitted.values
    value <- if (names) {
        ynames <- colnames(coef)
        if (is.null(ynames))
            ynames <- {
                lhs <- object$terms[[2L]]
                if (mode(lhs) == "call" && lhs[[1L]] == "cbind")
                    as.character(lhs)[-1L]
                else paste0("Y", seq_len(ny))
            }
        ind <- ynames == ""
        if (any(ind))
            ynames[ind] <- paste0("Y", seq_len(ny))[ind]
        setNames(vector("list", ny), paste("Response", ynames))
    }
    else vector("list", ny)
    cl <- oldClass(object)
    class(object) <- cl[match("mlm", cl):length(cl)][-1L]
    object$call$formula <- formula(object)
    for (i in seq(ny)) {
        object$coefficients <- setNames(coef[, i], rownames(coef))
        object$residuals <- resid[, i]
        object$fitted.values <- fitted[, i]
        object$effects <- effects[, i]
        object$call$formula[[2L]] <- object$terms[[2L]] <- as.name(if (names)
            ynames[i]
            else paste0("Y", i))
        value[[i]] <- summary(object, ...)

        # add to reduce memory use
        value[[i]]$cov.unscaled <- NULL
        gc()
    }
    class(value) <- "listof"
    value
}

# used instead of broom::tidy so that can use above summary.mlm
tidy.mlm <- function(x, s, conf.int = FALSE, conf.level = 0.95) {
    co <- stats::coef(s)
    nn <- c("estimate", "std.error", "statistic", "p.value")
    ret <- broom:::map_as_tidy_tibble(
        co,
        new_names = nn[1:ncol(co[[1]])],
        id_column = "response")
    ret$response <- stringr::str_replace(ret$response, "Response ", "")
    ret <- tibble::as_tibble(ret)
    if (conf.int) {
        CI <- tryCatch(stats::confint(x, level = conf.level),
                       error = function(x) {
                           NULL
                       })
        if (is.null(CI)) {
            CI <- confint_mlm(x, level = conf.level)
        }
        colnames(CI) <- c("conf.low", "conf.high")
        ret <- dplyr::bind_cols(ret, as_tibble(CI))
    }
    as_tibble(ret)
}

#' Fits a linear model on Z-scores
#'
#' @param zdata  A list with the zscore matrix and index object
#' @return       The coefficients matrix [gene x pathway]
zscore2model = function(zdata, hpc_args=NULL) {
    index = zdata$index

    # TODO: understand NAs
    zscores = t(zdata$zscores) * index$sign
    gene.nas = apply(is.na(zscores), 2, sum)
    gene.keep = colnames(zscores)[gene.nas < 10]
    zscores = zscores[, gene.keep]

    # fit model to pathway perturbations
    pathway = t(narray::mask(index$pathway)) + 0
    pathway = t(pathway)

    # TODO: data-driven linking of pathways (e.g. correlation)
    # pathway["PI3K",] = pathway["EGFR",] + pathway["PI3K",]
    # pathway["MAPK",] = pathway["EGFR",] + pathway["MAPK",]
    # pathway["NFkB",] = pathway["TNFa",] + pathway["NFkB",]

    mod <- lm(zscores ~ 0 + pathway)
    summary_list <- summary.mlm(mod)
    summary_tb <- tidy.mlm(mod, summary_list)
    summary_tb <- summary_tb |>
        transmute(gene = response,
                  pathway = sub("^pathway", "", term),
                  zscore = estimate,
                  p.value = p.value) |>
        mutate(adj.p = p.adjust(p.value, method="fdr"))


    #NOTE: Original shuttle lm replaced by stats::lm and broom::tidy() reformat
    #---
    #index = zdata$index
    #zscores = t(zdata$zscores) * index$sign

    # fit model to pathway perturbations
    #pathway = t(ar$mask(index$pathway)) + 0
    #pathway["EGFR",] = pathway["EGFR",] + pathway["MAPK",] + pathway["PI3K",]
    #pathway["TNFa",] = pathway["TNFa",] + pathway["NFkB",]

    #mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=30, atomic="pathway",
    #            hpc_args=hpc_args) %>%
    #    transmute(gene = zscores,
    #              pathway = sub("^pathway", "", term),
    #              zscore = estimate,
    #              p.value = p.value) %>%
    #    mutate(adj.p = p.adjust(p.value, method="fdr"))
    #---

    zfit = narray::construct(zscore ~ gene + pathway, data=summary_tb)
    pval = narray::construct(p.value ~ gene + pathway, data=summary_tb)

    # filter zfit to only include top 100 genes per pathway
    model = zfit
    model[apply(pval, 2, function(p) !min_mask(p, 100))] = 0

    list(assocs=summary_tb, model=model, zfit=zfit)
}

# load z-scores and convert probe --> symbol
zdata <- readRDS('data-raw/progeny/zscores.rds')

# setup mock eset to get map
es_probes <- row.names(zdata$zscores)
adat <- matrix(NA, length(es_probes), 2, dimnames=list(es_probes, c('s1', 's2')))
fdat <- AnnotatedDataFrame(data.frame(PROBE=es_probes, row.names = es_probes, stringsAsFactors = FALSE))
eset <- ExpressionSet(adat, featureData = fdat)

# get map
ensql <- '/mnt/12tb/Batcave/GEO/crossmeta/data-raw/entrezdt/ensql.sqlite'
annotation(eset) <- 'GPL96'

map <- fData(crossmeta::symbol_annot(eset, ensql = ensql))
map <- map[, c('PROBE', 'SYMBOL')]
map <- map[!is.na(map$SYMBOL), ]

# add symbol
zscores <- zdata$zscores
zscores <- zscores[map$PROBE, ] #expands 1:many
row.names(zscores) <- make.unique(map$SYMBOL) # one duplicated gene ARID4B

# restrict to signal genes (ligand or receptor)
signal.genes <- readRDS("data-raw/progeny/signal_genes.rds")


zdata$zscores <- zscores
model <- zscore2model(zdata)
saveRDS(model, 'data-raw/progeny/model.rds')
