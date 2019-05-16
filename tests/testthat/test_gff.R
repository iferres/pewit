context('gff')

example_data <- system.file('extdata', 'Hinfluenzae.tar.gz', package = 'pewit')
untar(tarfile = example_data, files = 'Hinfluenzae_2019.gff')
example_gff <- 'Hinfluenzae_2019.gff'


test_that("extracting table from gff works", {
  rl <- readLines(example_gff)
  x <- pewit:::extractGffTable(rl)
  d <- dim(x)
  d1 <- d[1]
  d2 <- d[2]
  colc <- vapply(x, class, NA_character_)
  expect_is(x, class = "data.frame")
  expect_equal(d1, 1975)
  expect_equal(d2, 10)
  expect_named(x, c("Contig", "ID", "LocusTag", "Gene",
                    "Product", "Type", "From", "To",
                    "Strand", "Phase"))
  expect_identical(colc, structure(c("character", "character", "character", "character",
                                     "character", "character", "integer", "integer",
                                     "character","character"),
                                   .Names = c("Contig", "ID", "LocusTag", "Gene",
                                              "Product", "Type", "From", "To",
                                              "Strand", "Phase")))
})

# require(Biostrings)
test_that("extracting sequences from gff works",{
  x <- pewit:::extractSeqsFromGff3(example_gff)
  meta <- mcols(x)
  d <- dim(meta)
  colc <- vapply(meta, class, NA_character_)
  expect_is(x, 'DNAStringSet')
  expect_length(x, 1916)
  expect_is(meta, 'DataFrame')
  expect_equal(d[1], 1916)
  expect_equal(d[2], 7)
  expect_named(meta, c("geneName", "organism", "product", "Contig",
                       "From", "To", "Strand"))
  expect_identical(colc, structure(c("character", "character", "character", "character",
                                     "integer", "integer", "character"),
                                   .Names = c("geneName", "organism", "product",
                                              "Contig", "From", "To", "Strand")))
})


file.remove(example_gff)
