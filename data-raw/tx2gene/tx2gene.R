# get tx2gene
library(rkal)
library(tidyr)
library(dplyr)
library(tibble)
library(dseqr.data)

tx2gene <- dseqr.data:::get_tx2gene(release = '94')
tx2gene_mouse <- dseqr.data:::get_tx2gene(species = "Mus musculus", release = "98")
tx2gene_mouse <- dseqr.data:::add_hgnc(tx2gene_mouse, 'Mus musculus', 98)

saveRDS(tx2gene_mouse, "inst/extdata/tx2gene_mouse.rds")

# get entrezid --> HGNC map used by cmap_es_ind and l1000_es
load("/mnt/12tb/Batcave/GEO/crossmeta/R/sysdata.rda")
hs$ENTREZID <- as.integer(hs$ENTREZID)

# unnest and join on entrezid
tx2gene_unnest <- unnest(tx2gene, entrezid)
tx2gene_unnest <- left_join(tx2gene_unnest, hs, by = c("entrezid" = "ENTREZID"))

summarise_symbol <- function(gnames, syms, species = "human") {
    if (species == "human") {
        # use SYMBOL_9606 if available (from CMAP data)
        syms <- unique(toupper(na.omit(syms)))
        if (length(syms) == 1) {
            return(syms)
        }

        # otherwise gene_name
        gnames <- unique(toupper(na.omit(gnames)))
        return(gnames)
    } else {
        unique(na.omit(gnames))
    }
}

# re-nest
tx2gene <- group_by(tx2gene_unnest, tx_id) %>%
    summarise(
        gene_name = summarise_symbol(gene_name, SYMBOL_9606, species = "human"),
        entrezid = entrezid[1],
        gene_id = unique(gene_id),
        seq_name = unique(seq_name),
        description = unique(description)
    )

# add homolog column to human tx2gene for equivalence with other species
tx2gene$hsapiens_homolog_ensembl_gene <- tx2gene$gene_id

# need tx_id, gene_name and entrezid for rkal::load_seq
# need gene_id for annotation
# need description for app
saveRDS(tx2gene, "inst/extdata/tx2gene.rds")


# check concordance with l1000_es/cmap_es ----
data_dir <- system.file("extdata", package = "dseqr.data")
cmap_es_ind <- qs::qread(file.path(data_dir, "cmap_es_ind.qs"))
l1000_es <- qs::qread(file.path(data_dir, "l1000_drugs_es.qs"))

sum(!row.names(cmap_es_ind) %in% tx2gene$gene_name)
# EnsDb89: 207
# EnsDb90: 203
# EnsDb91: 201
# EnsDb92: 220
# EnsDb94: 169
# EnsDb95: 197
sum(!row.names(l1000_es) %in% tx2gene$gene_name)
# EnsDb89: 4
# EnsDb90: 4
# EnsDb91: 4
# EnsDb92: 4
# EnsDb94: 3
# EnsDb95: 4
