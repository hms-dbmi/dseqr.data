# this script sets up L1000_pdata.rds
# which are used to match Pubchem CIDS to signatures

library(Biobase)
data_dir <- file.path('inst', 'extdata')

# setup pdata for L1000 ------
cmap_new_path <- '/home/alex/Documents/Batcave/zaklab/cmap_new/data-raw/'
lincs2_es <- qs::qread(file.path(cmap_new_path, 'siamese/data/lincs2_embedding.qs'))
siginfo <- qs::qread(file.path(cmap_new_path, 'lincs_beta/processed/siginfo.qs'))

# save siamese model for embedding of query data
file.copy(file.path(cmap_new_path, 'siamese/models/siamese_model.pt'),
          file.path(data_dir, 'siamese_model.pt'))

# save gene info needed to run embedding model
geneinfo <- qs::qread(file.path(cmap_new_path, 'lincs_beta/processed/geneinfo.qs'))
keep.genes <- geneinfo$feature_space %in% c('best inferred', 'landmark')
geneinfo <- geneinfo[keep.genes, ]
geneinfo <- geneinfo[, c('gene_id', 'gene_symbol')]
qs::qsave(geneinfo, file.path(data_dir, 'siamese_geneinfo.qs'))

# make pert titles match previous lincs2 format
lincs2_perts <- gsub('_', '.', siginfo$cmap_name)
lincs2_types <- gsub('trt_', '-', siginfo$pert_type)
lincs2_types <- gsub('-cp', '', lincs2_types)
lincs2_cells <- gsub('_', '.', siginfo$cell_iname)
lincs2_doses <- gsub(' ', '', siginfo$pert_idose)
lincs2_doses <- gsub('/', '.', lincs2_doses)
lincs2_doses[lincs2_doses == ''] <- '-700-666.0'
lincs2_duration <- gsub(' ', '', siginfo$pert_itime)
lincs2_duration[lincs2_duration == ''] <- '?h'

cmpinfo <- qs::qread(file.path(cmap_new_path, 'lincs_beta/processed/compoundinfo.qs'))

# pretty up pubchem cids
cmpinfo$pubchem_cid[cmpinfo$pubchem_cid %in% c('0', '-666')] <- NA

# lincs2_genes_pdata has columns 'title' and 'Samples(n)'
# lincs2_drugs_pdata has columns 'title', 'Pubchem CID', and 'Samples(n)'

lincs2_title <- paste(lincs2_perts, lincs2_cells, lincs2_doses, lincs2_duration, sep='_')

# determine genetic perts
is.genetic <- siginfo$pert_type != 'trt_cp'

# join siginfo with cmpinfo to get pubchem
cmpinfo <- cmpinfo[, c('pert_id', 'pubchem_cid')]
cmpinfo <- unique(cmpinfo)

siginfo <- siginfo[, c('pert_id', 'nsample')]
siginfo <- dplyr::left_join(siginfo, cmpinfo, by = 'pert_id')

# construct pdatas
lincs2_pdata <- data.frame(
    title = lincs2_title,
    `Pubchem CID` = siginfo$pubchem_cid,
    `Samples(n)` = siginfo$nsample,
    check.names = FALSE
)

# rename lincs2_es cols as title and split off genetic perts
colnames(lincs2_es) <- lincs2_title

lincs2_genes_es <- lincs2_es[, is.genetic]
lincs2_drugs_es <- lincs2_es[, !is.genetic]

# pubchem cids don't exists for genes
lincs2_genes_pdata <- lincs2_pdata[is.genetic, -2]
lincs2_drugs_pdata <- lincs2_pdata[!is.genetic, ]

# save pdata and overwrite existing lincs2_es data with fixed names
# ALSO UPDATE S3 DATA for lincs2_es in dseqr BUCKET IF ANYTHING CHANGES
qs::qsave(lincs2_genes_es, file.path(data_dir, 'lincs2_genes_es.qs'))
qs::qsave(lincs2_drugs_es, file.path(data_dir, 'lincs2_drugs_es.qs'))
saveRDS(lincs2_genes_pdata, file.path(data_dir, 'lincs2_genes_pdata.rds'))
saveRDS(lincs2_drugs_pdata, file.path(data_dir, 'lincs2_drugs_pdata.rds'))
