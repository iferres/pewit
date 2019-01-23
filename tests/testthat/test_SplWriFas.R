context('splitAndWriteFastas')

sq <- structure(c("M", "K", "T", "E", "S", "M", "I", "L", "R", "Q",
                  "S", "L", "Q", "N", "R", "I", "V", "H", "W", "G", "V", "A", "F",
                  "S", "T", "F", "I", "L", "I", "A", "S", "G", "I", "F", "Q", "M",
                  "P", "V", "S", "K", "R", "Y", "M", "I", "N", "E", "L", "P", "L",
                  "M", "A", "W", "S", "G", "D", "Y", "H", "I", "S", "L", "M", "L",
                  "H", "Y", "V", "G", "A", "F", "G", "L", "I", "F", "F", "V", "A",
                  "F", "H", "L", "Y", "F", "H", "I", "A", "R", "A", "E", "F", "D",
                  "I", "F", "P", "K", "K", "G", "D", "G", "V", "K", "S", "L", "K",
                  "I", "I", "K", "A", "M", "L", "F", "G", "G", "Q", "E", "P", "K",
                  "S", "E", "K", "Y", "L", "P", "E", "Q", "R", "L", "A", "Y", "F",
                  "F", "I", "G", "L", "V", "L", "L", "L", "L", "I", "I", "T", "G",
                  "L", "I", "K", "T", "L", "K", "N", "L", "A", "G", "W", "N", "I",
                  "S", "D", "S", "L", "Y", "L", "W", "S", "A", "Q", "L", "H", "N",
                  "L", "G", "M", "F", "L", "I", "I", "L", "G", "I", "I", "G", "H",
                  "L", "A", "A", "F", "I", "F", "K", "A", "N", "R", "P", "L", "L",
                  "R", "A", "M", "F", "S", "G", "K", "V", "D", "S", "H", "Y", "I",
                  "V", "E", "R", "H", "S", "L", "W", "Q", "E", "G", "V", "K", "M",
                  "A", "R", "D", "K", "V", "A", "E", "N", "N", "L", "D", "S", "I",
                  "S", "C", "H", "I", "K", "P", "L", "G", "E", "V", "S", "N", "T",
                  "E", "S", "T", "Q", "N", "K", "D", "*"),
                name = "37372.27_00001",
                Annot = "hypothetical protein",
                class = "SeqFastaAA")

#create fake multi sequence
atts <- attributes(sq)
sq <- paste0(sq, collapse = '')
attributes(sq) <- atts
fastas <- list(sq_1 = sq, sq_2 = sq, sq_3 = sq, sq_4 = sq)
tmps <- pewit:::splitAndWriteFastas(fastas, 3)

test_that('splitAndWriteFastas works', {
  expect_length(tmps, 3)
  expect_is(tmps, 'character')
  expect_true(all(file.exists(tmps)))
  expect_true(all(file.info(tmps)$size %in% c(259, 518, 259))) #order depends
  # on names given by tempfile() function.

})

file.remove(tmps)
