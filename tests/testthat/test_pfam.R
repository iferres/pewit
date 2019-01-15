context('pfam')

example_file <- system.file('testdata', 'Pfam-A.hmm.dat.tar.gz', package = 'pewit')
untar(example_file, files = 'Pfam-A.hmm.dat', exdir = tempdir())
example_pfamdat <- list.files(path = tempdir(),
                              pattern = '^Pfam-A.hmm.dat$',
                              full.names = TRUE)
ref_pfamdat_file <- system.file('testdata', 'processed_pfamadat.rds', package = 'pewit')
ref_pfamdat <- readRDS(ref_pfamdat_file)

x <- pewit:::processPfam_A_Dat(example_pfamdat, n_threads = 1)
xcln <- c("ID", "TP", "CL")
d <- dim(x)

test_that('processPfam_A_Dat works', {
  expect_is(x, 'data.frame')
  expect_identical(d[1], 6L)
  expect_identical(d[2], 3L)
  expect_identical(colnames(x), xcln)
  expect_identical(x, ref_pfamdat)
})

file.remove(example_pfamdat)
