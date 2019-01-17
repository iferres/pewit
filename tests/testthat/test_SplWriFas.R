context('splitAndWriteFastas')

sq_file <- system.file('testdata', 'gene_SeqFastaAA.rds', package = 'pewit')
sq <- readRDS(sq_file)

#create fake multi sequence
atts <- attributes(sq)
sq <- paste0(sq, collapse = '')
attributes(sq) <- atts
fastas <- list(sq_1 = sq, sq_2 = sq, sq_3 = sq, sq_4 = sq)

test_that('splitAndWriteFastas works', {

  tmps <- pewit:::splitAndWriteFastas(fastas, 3)
  expect_length(tmps, 3)
  expect_is(tmps, 'character')
  expect_true(all(file.exists(tmps)))
  expect_true(all(file.info(tmps)$size %in% c(259, 518, 259))) #order depends
  # on names given by tempfile() function.

})

file.remove(tmps)
