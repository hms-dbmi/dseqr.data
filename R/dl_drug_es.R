#' Download drug effect size data.
#'
#' @param files Character vector of file names to download.
#' @param check Check that existing drug effect size data is loadable? Default is FALSE.
#'
#' @return NULL
#' @keywords internal
#'
dl_drug_es <- function(files = c('cmap_es_ind.rds', 'l1000_drugs_es.rds', 'l1000_genes_es.rds'), check = FALSE) {

  timeout <- options()$timeout
  options(timeout = 120)

  # make sure doesn't already exist
  dest_dir <- system.file(package = 'drugseqr.data', mustWork = TRUE)
  dest_dir <- file.path(dest_dir, 'extdata')
  dir.create(dest_dir, showWarnings = FALSE)

  can_load <- c()
  exist_files <- file.exists(file.path(dest_dir, files))
  exist_files <- files[exist_files]
  if (length(exist_files) & check) {
    message(paste(exist_files, collapse = ' and '), ' already exists.')

    # check that can load
    for (exist_file in exist_files) {
      message('Checking that ', exist_file, ' can be loaded.')

      # store file name if can load
      fname <- tryCatch({
        drug_es <- readRDS(file.path(dest_dir, exist_file))
        exist_file
      },
      error = function(err) {
        message("Couldn't load ", exist_file, '. Will download.')
        unlink(file.path(dest_dir, exist_file))
        return(NULL)
      })
      can_load <- c(can_load, fname)
    }
    exist_files <- can_load
  }

  need_files <- setdiff(files, can_load)
  if (!length(need_files)) return(NULL)

  for (need_file in need_files) {
    message('downloading: ', need_file)
    dl_url <- paste0('https://s3.us-east-2.amazonaws.com/drugseqr/', need_file)
    download.file(dl_url, file.path(dest_dir, need_file))
  }

  options(timeout = timeout)
}


#' Load drug effect size data
#'
#' Downloads requested file if not done so previously.
#'
#' @param file Character vector of drug effect size datasets to load. One of
#'   \code{'cmap_es_ind.rds'} (CMAP02),
#'   \code{'l1000_drugs_es.rds'} (L1000 compounds), or
#'   \code{'l1000_genes_es.rds'} (L1000 genetic perturbations).
#'
#' @return data.frame of expression values. Rows are genes, columns are
#'   perturbations.
#' @export
#'
#' @examples
#'
#' # dummy example (actual files are large)
#' load_drug_es('example.rds')
#'
load_drug_es <- function(file = c('cmap_es_ind.rds', 'l1000_drugs_es.rds', 'l1000_genes_es.rds')) {

  dest_dir <- system.file(package = 'drugseqr.data', mustWork = TRUE)
  fpath <- file.path(dest_dir, 'extdata', file[1])


  drug_es <- NULL

  while(is.null(drug_es)) {

    drug_es <- tryCatch({
      readRDS(fpath)
    },
    error = function(err) {
      message("Couldn't load ", file)
      unlink(fpath)
      return(NULL)
    })

    if (is.null(drug_es))
      tryCatch(dl_drug_es(file),
               error = function(err) message("Couldn't download", file))

  }

  return(drug_es)
}
