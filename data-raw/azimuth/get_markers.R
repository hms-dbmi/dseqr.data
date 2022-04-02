# get azimuth specified markers
# downloaded from https://downgit.github.io/#/home?url=https://github.com/satijalab/azimuth_website/tree/master/static/csv

marker_files <- unzip('data-raw/azimuth/markers_2022-04-02.zip')

markers <- list()
for (i in seq_along(marker_files)) {
    marker_file <- marker_files[i]
    markers[[i]] <- read.csv(marker_file, row.names = 1, check.names = FALSE)
}

names(markers) <- tools::file_path_sans_ext(basename(marker_files))
qs::qsave(markers, 'inst/extdata/azimuth_markers.qs')

unlink('csv', recursive = TRUE)
