
check_sider <- function(cid) {
  url <- paste0('http://sideeffects.embl.de/drugs/', cid, '/')
  is_url <- RCurl::url.exists(url)
  return(is_url)
}

sider_old <- readRDS('data-raw/drug_annot/pug_view/sider.rds')
cids <- readRDS('data-raw/drug_annot/pug_view/cids.rds')

sider <- rep(FALSE, length(cids))
names(sider) <- cids

for (i in 1:length(cids)) {
  cid <- cids[i]
  if (cid %in% names(sider_old)) {
      sider[cid] <- sider_old[cid]
      next()
  }

  sider[cid] <- check_sider(cid)
  cat('Working on', i, 'of', length(cids), '\n')
}

saveRDS(sider, 'data-raw/drug_annot/pug_view/sider.rds')

sum(sider)
# 477
