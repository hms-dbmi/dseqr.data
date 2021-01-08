#' Get transcript to gene map.
#'
#' @param species Character vector indicating species. Genus and species should
#'   be space separated, not underscore. Default is \code{'Homo sapiens'}.
#' @param release EnsemblDB release. Should be same as used in
#'   \code{\link[rkal]{build_kallisto_index}}.
#' @param columns Character vector of columns from ensdb package to return or
#'   \code{'list'} to print available options.
#'
#' @return \code{data.frame} with columns \code{tx_id}, \code{gene_name}, and
#'   \code{entrezid}
#'
#' @keywords internal
#' @examples
#'
#' # tx2gene <- get_tx2gene()
get_tx2gene <- function(species = "Homo sapiens",
    release = "94",
    columns = c(
        "tx_id", "gene_name", "entrezid",
        "gene_id", "seq_name", "description"
    )) {

    # load EnsDb package
    ensdb_package <- get_ensdb_package(species, release)
    if (!require(ensdb_package, character.only = TRUE)) {
        build_ensdb(species, release)
        require(ensdb_package, character.only = TRUE)
    }

    # for printing available columns in EnsDb package
    if (columns[1] == "list") {
        print(ensembldb::listColumns(get(ensdb_package)))
        return(NULL)
    }

    # map from transcripts to genes
    tx2gene <- ensembldb::transcripts(get(ensdb_package),
        columns = columns, return.type = "data.frame"
    )
    tx2gene[tx2gene == ""] <- NA
    tx2gene <- tx2gene[!is.na(tx2gene$gene_name), ]

    if ("description" %in% colnames(tx2gene)) {
          tx2gene$description <- gsub(" \\[Source.+?\\]",
                                      "",
                                      tx2gene$description)
      }

    return(tx2gene)
}

#' Get ensembldb package name
#'
#' @inheritParams get_tx2gene
#'
#' @return Character vector with ensembldb package name. e.g.
#'   \code{'EnsDb.Hsapiens.v94'}.
#' @export
#' @examples
#'
#' get_ensdb_package(species = "Homo sapiens", release = "94")
get_ensdb_package <- function(species, release) {
    ensdb_species <- strsplit(species, " ")[[1]]
    ensdb_species[1] <- toupper(substr(ensdb_species[1], 1, 1))

    ensdb_package <- paste("EnsDb",
        paste0(ensdb_species, collapse = ""),
        paste0("v", release),
        sep = "."
    )
    return(ensdb_package)
}


#' Build and install ensembldb annotation package.
#'
#' @inheritParams get_tx2gene
#'
#' @return Called for side effects.
#' @keywords internal
#'
#' @examples
#'
#' # build ensembldb annotation package for human
#' # build_ensdb()
build_ensdb <- function(species = "Homo sapiens", release = "94") {

    # store ensembl databases in built package
    ensdb_dir <- "EnsDb"
    unlink("EnsDb", recursive = TRUE)
    dir.create(ensdb_dir)

    # format is genus_species in multiple other functions but not here
    species <- gsub("_", " ", species)

    # generate new ensembl database from specified release
    ah <- AnnotationHub::AnnotationHub()
    ahDb <- AnnotationHub::query(ah, pattern = c(species, "EnsDb", release))

    if (!length(ahDb)) {
          stop("Specified ensemble species/release not found in AnnotationHub.")
      }

    ahEdb <- ahDb[[1]]

    ensembldb::makeEnsembldbPackage(
        AnnotationDbi::dbfile(ensembldb::dbconn(ahEdb)),
        "0.0.1", "Alex Pickering <alexvpickering@gmail.com>",
        "Alex Pickering",
        ensdb_dir
    )

    # install new ensemble database
    ensdb_name <- list.files(ensdb_dir)
    ensdb_path <- file.path(ensdb_dir, ensdb_name)
    utils::install.packages(ensdb_path, repos = NULL)

    # remove source files
    unlink(ensdb_dir, recursive = TRUE)
}




#' Load transcript to gene map
#'
#' @inheritParams get_tx2gene
#'
#' @return \code{data.frame} with columns \code{tx_id}, \code{gene_name}, and
#'   \code{entrezid}
#' @export
#'
#' @examples
#'
#' tx2gene <- load_tx2gene("Homo sapiens", "94")
load_tx2gene <- function(species = "Homo sapiens", release = "94") {
    if (grepl("musculus", species)) {
        tx2gene <- readRDS(system.file("extdata",
            "tx2gene_mouse.rds",
            package = "drugseqr.data"
        ))
    } else if (grepl("sapiens", species)) {
        tx2gene <- readRDS(system.file("extdata",
            "tx2gene.rds",
            package = "drugseqr.data"
        ))
    } else {
        tx2gene <- get_tx2gene(species, release,
            columns = c("tx_id", "gene_name", "entrezid")
        )
    }

    return(tx2gene)
}
