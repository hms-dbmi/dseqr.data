#' Download drug effect size and Azimuth reference data.
#'
#' Will download to R package extdata folder if the environmental variable
#' `DSEQR_DATA_PATH` is unset.
#'
#' @param files Character vector of file names to download.
#' @param check Check that existing data is loadable? Default
#'   is FALSE.
#'
#' @return Downloads data into package folder.
#' @export
#' @examples
#'
#' dl_data('example.qs')
#'
dl_data <- function(files = c("cmap_es_ind.qs",
                              "l1000_drugs_es.qs",
                              "l1000_genes_es.qs",
                              "human_pbmc.qs",
                              "human_lung.qs",
                              "human_motorcortex.qs",
                              "mouse_motorcortex.qs",
                              "zhang_reference.qs"),
    check = FALSE) {


    # default to downloading to package directory
    dest_dir <- Sys.getenv('DSEQR_DATA_PATH')
    if (dest_dir == "") {
        dest_dir <- system.file(package = "dseqr.data", mustWork = TRUE)
        dest_dir <- file.path(dest_dir, "extdata")
    }

    dir.create(dest_dir, showWarnings = FALSE)

    timeout <- options()$timeout
    options(timeout = 600)


    # make sure doesn't already exist
    can_load <- c()
    exist_files <- file.exists(file.path(dest_dir, files))
    exist_files <- files[exist_files]
    if (length(exist_files) & check) {
        message(paste(exist_files, collapse = " and "), " already exists.")

        # check that can load
        for (exist_file in exist_files) {
            message("Checking that ", exist_file, " can be loaded.")

            # store file name if can load
            fname <- tryCatch(
                {
                    drug_es <- qs::qread(file.path(dest_dir, exist_file))
                    exist_file
                },
                error = function(err) {
                    message("Couldn't load ", exist_file, ". Will download.")
                    unlink(file.path(dest_dir, exist_file))
                    return(NULL)
                }
            )
            can_load <- c(can_load, fname)
        }
        exist_files <- can_load
    }

    need_files <- setdiff(files, can_load)
    if (!length(need_files)) {
        return(NULL)
    }

    for (need_file in need_files) {
        message("downloading: ", need_file)
        dl_url <- paste0("https://s3.us-east-2.amazonaws.com/dseqr/",
                         need_file)
        utils::download.file(dl_url, file.path(dest_dir, need_file))
    }

    options(timeout = timeout)
}


#' Load drug effect size and Azimuth reference data
#'
#' Downloads requested file if not done so previously.
#'
#' @param file Character vector of drug effect size datasets to load. One of:
#'   * `'cmap_es_ind.qs'` - CMAP02
#'   * `'l1000_drugs_es.qs'` - L1000 compounds
#'   * `'l1000_genes_es.qs'` - L1000 genetic perturbations
#'   * `'human_pbmc.qs'` - Azimuth human PBMC reference
#'
#' @return data.frame of expression values. Rows are genes, columns are
#'   perturbations.
#' @export
#'
#' @examples
#'
#' # dummy example (actual files are large)
#' load_drug_es("example.qs")
load_data <- function(
    file = c("cmap_es_ind.qs", "l1000_drugs_es.qs", "l1000_genes_es.qs", "human_pbmc.qs")) {

    # default to loading from package directory
    dest_dir <- Sys.getenv('DSEQR_DATA_PATH')
    if (dest_dir == "") {
        dest_dir <- system.file(package = "dseqr.data", mustWork = TRUE)
        dest_dir <- file.path(dest_dir, "extdata")
    }

    fpath <- file.path(dest_dir, file[1])
    data <- NULL

    i <- 1
    while (is.null(data) && i < 5) {
        i <- i + 1
        data <- tryCatch(
            {
                qs::qread(fpath)
            },
            error = function(err) {
                message("Couldn't load ", file)
                unlink(fpath)
                return(NULL)
            }
        )

        if (is.null(data)) {
            tryCatch(dl_data(file, dest_dir = dest_dir),
                     error = function(err) message("Couldn't download", file)
            )
        }
    }

    return(data)
}
