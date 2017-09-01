
#' @name realocateSingletons
#' @title Refining Singletons
#' @description Performs the final refinement which consist in evaluate the
#' singletons against the accessory genes. This step is performed because it
#' has benn seen that previous steps tends to generate more singletons than
#' there really are, specially if working with many genomes (100+).
#' @param clusters A \code{list} of clusters.
#' @param panm A binary panmatrix.
#' @param fastas A \code{list} of sequences.
#' @param n_threads \code{integer} The number of cpus to use.
#' @param seed \code{integer}. The algorithm uses a random step. This parameter
#' is used to allow reproducible results.
#' @return A new \code{list} of clusters.
#' @author Ignacio Feres
realocateSingletons <- function(clusters,
                                panm,
                                fastas,
                                n_threads,
                                seed = numeric()){

  ln <- rowSums(panm)
  aog <- which(ln>1 & ln<ncol(panm))
  sog <- which(ln==1)

  mn <- parallel::mclapply(aog, function(x){


    # Takes up to five sequences from each accs cluster
    fc <- clusters[[x]]
    if (length(fc)>5) {set.seed(seed); fc <- sample(fc, 5)}
    fa <- fastas[fc]

    # Align
    ali <- align(rf = fa, type = 'AA', n_threads = 1, accu = TRUE)
    ap <- apply(ali, 1, function(x){paste0(x, collapse = '')})

    tmp <- tempfile()
    seqinr::write.fasta(sapply(ap, seqinr::s2c, simplify = FALSE),
                        names = names(ap),
                        file.out = tmp)

    #Build hmm model
    hmmModel <- hmmBuild(ali = tmp, name = names(clusters)[x])
    press <- hmmPress(hmmModel)

    # set threshold. hmmsearch model vs proteins used to build model,
    #  save 4/5 of the minimum threshold in tblout.
    se <- runHmmsearch(tmp, hmmModel, pfam = FALSE, n_threads = 1L)
    tbl <- readTblout(tblout = se)
    tr <- min(tbl$Score) * 4/5
    file.remove(se)
    file.remove(press[-1])

    # Output model filename and computed threshold.
    attr(hmmModel, 'threshold') <- tr
    return(hmmModel)


  }, mc.cores = n_threads)

  tr <- sapply(mn, attr, 'threshold')
  mn <- unlist(mn)

  # Concatenates all the generated models
  mdls <- paste0(tempdir(),'/accsModels.hmm')
  ct <- paste0('cat ', paste0(mn, collapse = ' '), ' > ', mdls)
  system(ct)
  file.remove(mn)
  press <- hmmPress(model = mdls)

  sqs <- unlist(clusters[sog])
  sing <- tempfile()
  seqinr::write.fasta(fastas[sqs], names = sqs, file.out = sing)

  # Search all singletons vs hmm models
  tblout <- runHmmsearch(fasta = sing,
                         hmm = mdls,
                         pfam = FALSE,
                         n_threads = n_threads)
  file.remove(press)
  sqs <- as.data.frame(cbind(sqs))
  m <- readTblout(tblout = tblout)[,-3]

  # Only conserve hits above the computed threshold.
  spl <- split(m, m$queryName)
  rm(m)
  filt <- lapply(names(spl), function(x){
    cu <- tr[x]
    ta <- spl[[x]]
    ta <- ta[which(ta$Score>=cu),]
    ta$singClu <- sapply(ta$targetName, function(y){
      rownames(sqs)[which(sqs$sqs==y)]
    })
    ta
  })
  m <- do.call(rbind, filt)

  spl <- split(m, m$targetName)
  en <- lapply(spl, function(x){
    x[which.max(x$Score),]
  })
  m <- do.call(rbind, en)

  clusters <- assignOrphans(m = m,
                            clusters = clusters,
                            panm = panm)
  return(clusters)

}



#' @name hmmBuild
#' @title hmmBuild
#' @description Wrapper function of \code{hmmbuild} (HMMER 3).
#' @param ali \code{character} A fasta file name.
#' @param name \code{character} The name of the model.
#' @return The name of the hmm model file.
#' @author Ignacio Ferres
hmmBuild <- function(ali, name){
  tmp <- tempfile()
  hmmbuild <- paste('hmmbuild -o /dev/null --amino -n', name, tmp, ali)
  system(hmmbuild)
  tmp
}

#' @name hmmPress
#' @title hmmPress
#' @description Wrapper function of \code{hmmpress} (HMMER 3).
#' @param model \code{character} The name of the hmm file.
#' @return The names of indexed files.
#' @author Ignacio Ferres
hmmPress <- function(model){
  hmmpress <- paste('hmmpress', model)
  system(hmmpress, ignore.stdout = TRUE)
  o <- paste0(model, c('','.h3f', '.h3i', '.h3m', '.h3p'))
  o
}


#' @name readTblout
#' @title Read tblout Formated Files
#' @description Read hmmsearch output.
#' @param tblout The file name.
#' @return A \code{data.frame}.
#' @author Ignacio Ferres
readTblout <- function(tblout) {
  rl <- readLines(tblout)
  rl <- rl[which(!grepl("^\\#", rl))]
  rl <- gsub("[ ]+", " ", rl)
  lst <- strsplit(rl, " ")

  targetName <- sapply(lst, function(x) {
    x[1]
  })
  queryName <- sapply(lst, function(x) {
    x[3]
  })
  evalue <- sapply(lst, function(x) {
    x[5]
  })
  score <- sapply(lst, function(x) {
    x[6]
  })

  hmmer.table <- data.frame(targetName = targetName,
                            queryName = queryName,
                            Evalue = as.numeric(evalue),
                            Score = as.numeric(score),
                            stringsAsFactors = F)

  return(hmmer.table)
}

#' @name assignOrphans
#' @title Assign Orphan Clusters to other clusters.
#' @description .
#' @param m A tblout \code{data.frame}.
#' @param clusters A \code{list} of clusters.
#' @param panm A \code{data.frame}. The panmatrix.
#' @return A \code{list} of refined clusters.
#' @author Ignacio Ferres
assignOrphans <- function(m,
                          clusters,
                          panm){
  for (i in 1:nrow(m)){
    ocl <- clusters[m$singClu[i]]
    at <- attr(ocl[[1]], 'paralogues')
    org <- strsplit(ocl[[1]][1], ';')[[1]][1]
    hav <- panm[ m$queryName[i], org ]

    if (hav == 1){
      atv <- attr(clusters[[m$queryName[i]]], 'paralogues')
      atv <- c(atv, ocl[[1]][1], at)
      attr(clusters[[m$queryName[i]]], 'paralogues') <- atv
      clusters[m$singClu[i]] <- NULL
    }else{
      ln <- length(clusters[[m$queryName[i]]])
      clusters[[m$queryName[i]]][ln + 1] <- ocl[[1]][1]
      atv <- attr(clusters[[m$queryName[i]]], 'paralogues')
      atv <- c(atv, at)
      attr(clusters[[m$queryName[i]]], 'paralogues') <- atv
      clusters[m$singClu[i]] <- NULL
      panm[ m$queryName[i], org ] <- 1L
    }


  }

  return(clusters)
}




