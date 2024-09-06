# get unique compound ids for 2-get_views.sh

cmap_pdata <- readRDS('inst/extdata/CMAP02_pdata.rds')
l1000_drugs_pdata <- readRDS('inst/extdata/L1000_drugs_pdata.rds')
lincs2_drugs_pdata <- readRDS('inst/extdata/lincs2_drugs_pdata.rds')

cids <- c(cmap_pdata$`Pubchem CID`, l1000_drugs_pdata$`Pubchem CID`, lincs2_drugs_pdata$`Pubchem CID`)

names(cids) <- c(cmap_pdata$title, l1000_drugs_pdata$title, lincs2_drugs_pdata$title)

names(cids) <- gsub('^([^_]+)_.+?$', '\\1', names(cids))
cids <- cids[!duplicated(cids) & !is.na(cids)]

cids_old <- readRDS('data-raw/drug_annot/pug_view/cids.rds')
cids_new <- setdiff(cids, cids_old)

# save new for getting views (slow)
write.table(cids_new, 'data-raw/drug_annot/pug_view/cids.csv', row.names = FALSE, col.names = FALSE, quote = FALSE)

# save old for history
saveRDS(cids_old, 'data-raw/drug_annot/pug_view/cids_old.rds')

# overwrite new
saveRDS(cids, 'data-raw/drug_annot/pug_view/cids.rds')
