test_that("download works", {
  expect_equal(load_drug_es('example.rds'), 'hello there!')
})
