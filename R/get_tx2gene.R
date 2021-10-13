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
                        release = '103',
                        with_hgnc = FALSE,
                        columns = c("tx_id", "gene_name", "entrezid",
                                    "gene_id", "seq_name", "description")
) {

  # use latest
  if (is.null(release)) {
    ah <- AnnotationHub::AnnotationHub()
    ahDb <- AnnotationHub::query(ah, pattern = c(species, "EnsDb"))
    last <- tail(ahDb$title, 1)
    release <- gsub('^Ensembl ([0-9]+).+?$', '\\1', last)
  }

  # load EnsDb package
  ensdb_package <- get_ensdb_package(species, release)
  if (!require(ensdb_package, character.only = TRUE)) {
    ensdb_package <- build_ensdb(species, release)
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

  # try to add human hgnc symbols
  if (with_hgnc) tx2gene <- add_hgnc(tx2gene, species)
  return(tx2gene)
}

get_biomart_ensdb_species <- function(species) {

  # a few cases that don't match below pattern
  switch(species,
         'Bos indicus x Bos taurus' = return('bihybrid'),
         'Gorilla gorilla gorilla' = return('ggorilla'),
         'Mus musculus musculus' = return('mmusculus'),
         'Mus musculus domesticus' = return('mmusculus')
  )


  ensdb_species <- strsplit(species, " ")[[1]]
  nparts <- length(ensdb_species)

  for (i in seq_len(nparts-1)) {
    ensdb_species[i] <- tolower(substr(ensdb_species[i], 1, 1))
  }

  ensdb_species <- paste0(ensdb_species, collapse = "")
  return(ensdb_species)
}

add_hgnc <- function(tx2gene, species) {

  ensdb_species <- get_biomart_ensdb_species(species)

  mart <- tryCatch(
    biomaRt::useEnsembl(
      biomart = 'genes',
      dataset = paste0(ensdb_species, '_gene_ensembl')),
    error = function(e) {
      message(e)
      return(NULL)
    })

  # return if couldn't find
  if (is.null(mart)) return(tx2gene)

  map <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"),
    filters = "ensembl_gene_id",
    values = tx2gene$gene_id,
    mart = mart)

  # make useful
  map <- map[map$hsapiens_homolog_ensembl_gene != '', ]
  map <- map[!duplicated(map$ensembl_gene_id), ]
  map <- map[map$ensembl_gene_id %in% tx2gene$gene_id, ]

  # add to tx2gene
  tx2gene <- dplyr::left_join(tx2gene, map, by = c("gene_id" = 'ensembl_gene_id'))
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
#' get_ensdb_package(species = "Canis lupus familiaris", release = "99")
#'
get_ensdb_package <- function(species, release) {
  ensdb_species <- .abbrevOrganismName(.organismName(species))

  ensdb_package <- paste0(
    "EnsDb.",
    .abbrevOrganismName(.organismName(species)),
    ".v",
    release)

  return(ensdb_package)
}

# from ensembldb:::.organismName
.organismName <- function(x) {
  substring(x, 1, 1) <- toupper(substring(x, 1, 1))
  return(x)
}

# from ensembldb:::.abbrevOrganismName
.abbrevOrganismName <- function(organism) {
  spc <- unlist(strsplit(organism, "_|[[:space:]]"))
  return(paste0(substr(spc[[1]], 1, 1), spc[[2]]))
}


#' Build and install ensembldb annotation package.
#'
#' @inheritParams get_tx2gene
#'
#' @return EnsDb package name.
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

  return(ensdb_name)
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
load_tx2gene <- function(species = "Homo sapiens", release = '103', with_hgnc = FALSE) {
  if (grepl("musculus", species)) {
    tx2gene <- readRDS(system.file("extdata",
                                   "tx2gene_mouse.rds",
                                   package = "dseqr.data"
    ))
  } else if (grepl("sapiens", species)) {
    tx2gene <- readRDS(system.file("extdata",
                                   "tx2gene.rds",
                                   package = "dseqr.data"
    ))
  } else {
    tx2gene <- get_tx2gene(species, release, with_hgnc)
  }

  return(tx2gene)
}
