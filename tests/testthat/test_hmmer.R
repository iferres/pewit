context('Processing HMMER results')

example_domtblout_tgz <- system.file('testdata', 'example_domtblout.tar.gz', package = 'pewit')
untar(example_domtblout_tgz, files = 'example_domtblout.tab', exdir = tempdir())
example_domtblout_file <- list.files(path = tempdir(),
                                     pattern = '^example_domtblout.tab$',
                                     full.names = TRUE)

ref <- readRDS(system.file('testdata', 'processed_pfamadat2.rds', package = 'pewit'))
ccls <- structure(c("character", "character", "character", "numeric",
                    "numeric", "numeric", "numeric", "character",
                    "character", "character"),
                  .Names = c("Query", "Hit", "PfamID", "Evalue",
                             "Score", "Start", "End", "Description",
                             "Type", "Clan"))

test_that('outhmmsearch works', {

  x <- pewit:::outhmmsearch(example_domtblout_file, ref)
  expect_is(x, 'data.frame')
  expect_length(x, 10)
  expect_equal(dim(x), c(166, 10))
  expect_identical(sapply(x, class), ccls)

})


test_that('single hit grouping works', {
  x <- pewit:::outhmmsearch(example_domtblout_file, ref)
  y <- pewit:::grouping(x, singlehit = TRUE, type = 'Domain')

  expect_true(all(sapply(y, class)=='data.frame'))

  y_dims <- all(vapply(y, function(z){
    identical(dim(z), c(1L, 10L))
  }, FUN.VALUE = NA))
  expect_true(y_dims)

})


test_that('multiple hit grouping works',{
  x <- pewit:::outhmmsearch(example_domtblout_file, ref)
  y <- pewit:::grouping(x, singlehit = FALSE)

  expect_true(all(sapply(y, class)=='data.frame'))

  dims <- lapply(y, dim)
  expect_true(all(sapply(dims, '[', 1) > 1L))
  expect_true(all(sapply(dims, '[', 2) == 10L))

})

#fake overlapping domain matrix:
m <- matrix(data = c(T, T, F, T, T, F, F, F, T),
            nrow = 3, ncol = 3,
            dimnames = list(c('A', 'B', 'C'),
                            c('A', 'B', 'C')))

test_that('determining overlapping domains works',{
  x <- pewit:::determineOverlap(m)
  expect_is(x, 'list')
  expect_length(x, 2L)
  expect_length(x[[1]], 2L)
  expect_length(x[[2]], 1L)
  expect_true(all(unlist(x)%in%dimnames(m)[[1]]))
})


# sample_hit with four domains, two of them overlap and belong to the same clan.
sample_hit <- structure(list(Query = c("genomeA;gene_0001", "genomeA;gene_0001",
                                       "genomeA;gene_0001", "genomeA;gene_0001"),
                             Hit = c("Helicase_C", "ResIII", "UVR", "UvrB"),
                             PfamID = c("PF00271.30", "PF04851.14","PF02151.18", "PF12344.7"),
                             Evalue = c(9.9e-21, 6.1e-12, 1.3e-10,2.4e-24),
                             Score = c(66.9, 38.5, 33.5, 77.7),
                             Start = c(446, 20, 640, 557),
                             End = c(560, 102, 674, 598),
                             Description = c("-","-", "-", "-"),
                             Type = c("Family", "Family", "Family", "Family"),
                             Clan = c("CL0023", "CL0023", "No-Clan", "CL0023")),
                        .Names = c("Query","Hit", "PfamID", "Evalue", "Score",
                                   "Start", "End", "Description","Type", "Clan"),
                        row.names = c(76L, 121L, 158L, 161L),
                        class = "data.frame")

test_that('removing overlaping domains works',{
  x <- pewit:::rm_overlaping_clans(sample_hit, 'Family')
  expect_is(x, 'data.frame')
  expect_length(x, 10L)
  expect_equal(dim(x), c(3, 10))
  expect_true(all(colnames(x)%in%colnames(sample_hit)))
  expect_identical(sapply(x, class), sapply(sample_hit, class))
})




example_tblout_tgz <- system.file('testdata', 'example_tblout.tar.gz', package = 'pewit')
untar(example_tblout_tgz, exdir = tempdir())
example_tblout_file <- list.files(path = tempdir(),
                                  pattern = '^example_tblout.tab$',
                                  full.names = TRUE)

test_that('outphmmer works', {
  x <- pewit:::outphmmer(example_tblout_file)
  expect_is(x, 'data.frame')
  expect_equal(dim(x), c(3, 4))
  xcl <- structure(c("character", "character", "numeric", "numeric"),
                   .Names = c("Query","Hit", "Evalue", "Score"))
  expect_identical(sapply(x, class), xcl)
})




file.remove(example_domtblout_file)
file.remove(example_tblout_file)



