#' Download ensembl transcriptome and build index for salmon quantification
#'
#' This index is used for bulk RNA seq quantification. See \code{\link{build_gencode_index}} for single cell RNA-seq equivalent.
#'
#' @param species The species. Default is \code{homo_sapiens.}
#' @param release ensembl release. Default is \code{94} (latest in release for AnnotationHub -
#'   needs to match with \code{\link{build_ensdb}}) and corresponds to Gencode release 29 for \code{\link{build_gencode_index}}.
#'
#' @return NULL
#' @export
#'
#' @examples
#' # build salmon index for humans
#' build_salmon_index()
#'
build_salmon_index <- function(species = 'homo_sapiens', release = '94') {

  salmon_version <- get_pkg_version('salmon')
  if (salmon_version == '0.14.0')
    stop('EnsemblDb indices not yet implemented for salmon 0.14.0')

  indices_dir <- system.file(package = 'drugseqr.data')
  indices_dir <- file.path(indices_dir, paste0('indices/salmon_', salmon_version), 'ensdb')

  if (!dir.exists(indices_dir))
    dir.create(indices_dir, recursive = TRUE)

  # construct ensembl url for transcriptome
  ensembl_species <- gsub(' ', '_', tolower(species))
  ensembl_release <- paste0('release-', release)
  ensembl_url <- paste0('ftp://ftp.ensembl.org/pub/', ensembl_release, '/fasta/', ensembl_species, '/cdna/')

  # get list of all files
  handle <- curl::new_handle(dirlistonly=TRUE)
  con <- curl::curl(ensembl_url, "r", handle)
  tbl <- utils::read.table(con, stringsAsFactors=FALSE)
  close(con)

  # get transcripts cdna.all file
  ensembl_all <- grep('cdna.all.fa.gz$', tbl[[1]], value = TRUE)
  ensembl_url <- paste0(ensembl_url, ensembl_all)

  work_dir <- getwd()
  setwd(indices_dir)
  curl::curl_download(ensembl_url, ensembl_all)

  # build index
  tryCatch(system2(command, args=c('index',
                                   '-t', ensembl_all,
                                   '-i', ensembl_species)),
           error = function(err) {err$message <- 'Is salmon installed and on the PATH?'; stop(err)})

  unlink(ensembl_all)
  setwd(work_dir)
}

#' Download ensembl transcriptome and build index for salmon quantification
#'
#' This index is used for single cell RNA-seq quantification. See \code{\link{build_ensdb_index}} for bulk RNA-seq equivalent.
#'
#' @param indices_dir directory to place indices in.
#' @param species The species. Default is \code{homo_sapiens.}
#' @param release gencode release. Default is \code{29} (matches ensembl release 94 for \code{\link{build_ensdb_index}})).
#'
#' @return NULL
#' @export
#'
#' @examples
#' # build salmon alevin index for humans
#' indices_dir <- 'data-raw/indices'
#' build_gencode_index(indices_dir)
#'
build_gencode_index <- function(indices_dir, species = 'human', release = '29', command = 'salmon') {

  salmon_version <- get_pkg_version(command)
  salmon_old <- salmon_lt_0.14.0(salmon_version)

  # newer salmon uses decoys in index
  if (salmon_old)
    build_gencode_index_no_decoys(indices_dir, species, release)
  else
    build_gencode_index_decoys(indices_dir, species, release)
}

build_gencode_index_decoys <- function(indices_dir, species, release) {

  indices_dir <- file.path(indices_dir, '0.14.0', 'gencode')
  if (!dir.exists(indices_dir)) dir.create(indices_dir, recursive = TRUE)

  if (species != 'human' | release != '29')
    stop('Only implemented for human gencode v29.')

  gencode_file <- 'human_GENCODEv29.tar.gz'
  gencode_url <- paste0('https://s3.us-east-2.amazonaws.com/drugseqr/', gencode_file)

  work_dir <- getwd()
  setwd(indices_dir)

  download.file(gencode_url, gencode_file)
  utils::untar(gencode_file)

  # build index
  tryCatch(system2('salmon', args=c('index',
                                    '-t', 'human_GENCODEv29/gentrome.fa',
                                    '--gencode',
                                    '-d', 'human_GENCODEv29/decoys.txt',
                                    '-i', species)),
           error = function(err) {err$message <- 'Is salmon installed and on the PATH?'; stop(err)})


  unlink(c(gencode_file, 'human_GENCODEv29'), recursive = TRUE)
  setwd(work_dir)
}


build_gencode_index_no_decoys <- function(indices_dir, species, release) {

  indices_dir <- file.path(indices_dir, '0.13.1', 'gencode')
  if (!dir.exists(indices_dir)) dir.create(indices_dir, recursive = TRUE)

  # construct ensembl url for protein coding transcriptome
  gencode_file <- paste0('gencode.v', release, '.pc_transcripts.fa.gz')
  gencode_url <- paste0('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_', species, '/release_', release, '/', gencode_file)

  work_dir <- getwd()
  setwd(indices_dir)
  curl::curl_download(gencode_url, gencode_file)

  # build index
  tryCatch(system2('salmon', args=c('index',
                                    '-t', gencode_file,
                                    '--gencode',
                                    '-i', species)),
           error = function(err) {err$message <- 'Is salmon installed and on the PATH?'; stop(err)})

  unlink(gencode_file)
  setwd(work_dir)
}
