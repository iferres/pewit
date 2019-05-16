context('splitClusters')

rda_file <- system.file('testdata', 'preclust_1.rda', package = 'pewit')
load(rda_file) #preclust_1
x <- pewit:::splitCluster(preclust_1, sep = '___')

test_that("splitCluster works", {
  expect_is(x, 'data.frame')
  d <- dim(x)
  expect_equal(d[1], 30)
  expect_equal(d[2], 2)
  expect_named(x, c('Gene', 'NODE'))
  colc <- vapply(x, class, NA_character_)
  expect_identical(colc, structure(c("character", "character"),
                                   .Names = c("Gene", "NODE")))
  nNodes <- length(unique(x$NODE))
  expect_equal(nNodes, 3)
  tNodes <- table(x$NODE)
  expect_equivalent(tNodes, as.table(c(10, 10, 10)))
})
