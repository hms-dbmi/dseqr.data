#' Download drug effect size and Azimuth reference data.
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
                              "mouse_motorcortex.qs"),
    check = FALSE) {

    timeout <- options()$timeout
    options(timeout = 600)

    # make sure doesn't already exist
    dest_dir <- system.file(package = "dseqr.data", mustWork = TRUE)
    dest_dir <- file.path(dest_dir, "extdata")
    dir.create(dest_dir, showWarnings = FALSE)

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
                    drug_es <- readRDS(file.path(dest_dir, exist_file))
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
#' @param file Character vector of drug effect size datasets to load. One of
#'   \code{'cmap_es_ind.qs'} (CMAP02),
#'   \code{'l1000_drugs_es.qs'} (L1000 compounds),
#'   \code{'l1000_genes_es.qs'} (L1000 genetic perturbations), or
#'   \code{'human_pbmc.qs'} (Azimuth human PBMC reference).
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
    dest_dir <- system.file(package = "dseqr.data", mustWork = TRUE)
    fpath <- file.path(dest_dir, "extdata", file[1])


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
            tryCatch(dl_data(file),
                     error = function(err) message("Couldn't download", file)
            )
        }
    }

    return(data)
}
