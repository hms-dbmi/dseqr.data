test_that("tx2gene can be loaded", {
  tx2gene <- load_tx2gene()
  gene1 <- tx2gene[1,]$gene_name
  expect_equal(gene1, toupper(gene1))
})


test_that("tx2gene_mouse can be loaded", {
    tx2gene <- load_tx2gene('Mus musculus')
    gene1 <- tx2gene[1,]$gene_name
    expect_false(gene1 == toupper(gene1))
})
