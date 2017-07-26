

realocateSingletons <- function(final.clusters,
                                panm,
                                fastas,
                                n_threads,
                                seed = numeric()){

  ln <- rowSums(panm)
  aog <- which(ln>1 & ln<ncol(panm))
  sog <- which(ln==1)
  ge <- unlist(final.clusters[sog])

  tge <- as.data.frame(do.call(rbind, strsplit(ge, ';')))
  sge <- split(tge, tge$V1)

  f <- seqinr::a()#[-1]

  mn <- parallel::mclapply(aog, function(x){
    fc <- final.clusters[[x]]
    if (length(fc)>5) {set.seed(seed); fc <- sample(fc, 5)}
    fa <- fastas[fc]
    ali <- align(rf = fa, type = 'AA', n_threads = 1, accu = TRUE)
    ap <- apply(ali, 1, function(x){paste0(x, collapse = '')})

    tmp <- tempfile()
    seqinr::write.fasta(sapply(ap, seqinr::s2c, simplify = FALSE),
                        names = names(ap),
                        file.out = tmp)
    hmmModel <- hmmBuild(ali = tmp, name = names(final.clusters)[x])
    file.remove(tmp)

    return(hmmModel)

    # ssp <- strsplit(unlist(fastas[final.clusters[[x]]]),'')
    # sap <-sapply(ssp, function(z){table(factor(z, levels = f))}, simplify = F)
    # ori <- do.call(rbind, sap)
    # re <- colMeans(ori)
    # d <- as.data.frame(as.matrix(vegan::vegdist(rbind(ori,re), method = 'bray')))
    # rmv <- grep('^re$',rownames(d))
    # m <- min(d$re[-rmv])
    # M <- max(d$re[-rmv])
    # attr(re, 'max') <- M
    # attr(re, 'min') <- m
    # re
  }, mc.cores = n_threads)
  mn <- unlist(mn)

  mdls <- paste0(tempdir(),'/accsModels.hmm')
  ct <- paste0('cat ', paste0(mn, collapse = ' '), ' > ', mdls)
  system(ct)
  file.remove(mn)
  press <- hmmPress(model = mdls)

  sqs <- unlist(final.clusters[sog])
  sing <- tempfile()
  seqinr::write.fasta(fastas[sqs], names = sqs, file.out = sing)

  tblout <- runHmmsearch(fasta = sing,
                         hmm = mdls,
                         pfam = FALSE,
                         n_threads = n_threads)
  file.remove(press)




  # for (i in 1:length(sge)){
  #   or <- as.character(sge[[i]]$V1[1])
  #   gor <- paste(or, sge[[i]]$V2, sep = ';')
  #   lp <- lapply(fastas[gor], function(x){
  #     table(factor(strsplit(x, '')[[1]], levels = f))
  #   })
  #
  #   nd <- sapply(lp, function(x){
  #     lapply(mn, function(y){
  #       d <- as.vector(vegan::vegdist(rbind(x,y)))
  #       ifelse(d>attr(y, 'min') & d<attr(y, 'max'), TRUE, FALSE)
  #     })
  #   }, simplify = TRUE)
  #
  #   ap <- apply(nd, 2, function(x){which(unlist(x))})
  #
  #
  #
  #   # ...
  # }

}




hmmBuild <- function(ali, name){
  tmp <- tempfile()
  hmmbuild <- paste('hmmbuild -o /dev/null --amino -n', name, tmp, ali)
  system(hmmbuild)
  tmp
}


hmmPress <- function(model){
  hmmpress <- paste('hmmpress', models)
  system(hmmpress)
  paste0(model, c('','.h3f', '.h3i', '.h3m', '.h3p'))
}
