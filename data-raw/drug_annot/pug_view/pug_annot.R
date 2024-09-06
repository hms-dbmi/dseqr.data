source('data-raw/drug_annot/pug_view/parse_pug.R')

cids <- readRDS('data-raw/drug_annot/pug_view/cids.rds')

# load in previous
pug_annot_old <- readRDS('data-raw/drug_annot/pug_view/pug_annot.rds')

# setup new
pug_annot <- tibble::tibble(pubchem_cid = cids, pert_iname = names(cids), drugbank = NA_character_, wikipedia = NA_character_, gras = FALSE)

for (i in 1:length(cids)) {
  cat('Working on', i, 'of', length(cids), '\n')
  cid <- cids[i]

  # load pug view
  pug_file <- file.path('data-raw/drug_annot/pug_view/views', paste0(cid, '.json'))
  if (!file.exists(pug_file)) {
      pug_old <- pug_annot_old[pug_annot_old$pubchem_cid == cid, ]
      if (nrow(pug_old) == 0) next()

      pug_annot[i, 'drugbank'] <- pug_old$drugbank
      pug_annot[i, 'gras'] <- pug_old$gras
      pug_annot[i, 'wikipedia'] <- pug_old$wikipedia
      next()
  }

  pug_view <- rjson::fromJSON(file=pug_file)

  pug_annot[i, 'drugbank'] <- get_drugbank(pug_view)
  pug_annot[i, 'gras'] <- check_gras(pug_view)
  pug_annot[i, 'wikipedia'] <- get_wikipedia(pug_view)
}

sum(!is.na(pug_annot$drugbank))
# 1239

sum(pug_annot$gras)
# 7

sum(!is.na(pug_annot$wikipedia))
# 1704

saveRDS(pug_annot, 'pug_annot.rds')
