context('gff')

tgz <- system.file('testdata', 'test_gff.tar.gz', package = 'pewit')

untar(tgz, files = 'test.gff', exdir = tempdir())
test_gff <- list.files(path = tempdir(),
                       pattern = '^test.gff$',
                       full.names = TRUE)

gff_table_df_file <- system.file('testdata', 'gff_table_df.rds', package = 'pewit')
gff_table_df <- readRDS(gff_table_df_file)


test_that("extracting table from gff works", {
  rl <- readLines(test_gff)
  tt <- pewit:::extractGffTable(rl)
  d <- dim(tt)
  d1 <- d[1]
  d2 <- d[2]
  expect_is(tt, class = "data.frame")
  expect_equal(d1, 236)
  expect_equal(d2, 10)
  expect_identical(tt, gff_table_df)
})


gene_SeqFastadna_file <- system.file('testdata',
                                     'gene_SeqFastadna.rds',
                                     package = 'pewit')
gene_SeqFastadna <- readRDS(gene_SeqFastadna_file)

gene_SeqFastaAA_file <- system.file('testdata',
                                    'gene_SeqFastaAA.rds',
                                    package = 'pewit')
gene_SeqFastaAA <- readRDS(gene_SeqFastaAA_file)


test_that("translate works", {
  tr <- pewit:::translate(gene_SeqFastadna, numcode = 11)
  expect_equivalent(tr, unclass(gene_SeqFastaAA))
})


SeqFastadna_file <- system.file('testdata',
                                'SeqFastadna.rds',
                                package = 'pewit')
SeqFastadna <- readRDS(SeqFastadna_file)


test_that("extracting sequences from gff table works",{
  cds <- gff_table_df[1, ]
  sq <- pewit:::getFfnFaa(SeqFastadna,
                          contig = cds$Contig,
                          strand = cds$Strand,
                          from = cds$From,
                          to = cds$To,
                          id = cds$ID,
                          product = cds$Product)

  expect_is(sq, 'list')
  expect_length(sq, 2)
  expect_is(sq[[1]], "SeqFastadna")
  expect_is(sq[[2]], "SeqFastaAA")
  expect_equivalent(sq[[1]], gene_SeqFastadna)
  expect_equivalent(sq[[2]], gene_SeqFastaAA)
})


file.remove(test_gff)
