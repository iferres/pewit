#' @importFrom reshape2 melt
#' @importFrom S4Vectors elementNROWS mcols mcols<- List
fast_clust <- function(faas, verbose = TRUE){
  # Take largest as representative sequence
  faas <- faas[order(elementNROWS(faas), decreasing = TRUE)]
  minclus <- minhash_clust_k4(faas, verbose = verbose)
  mminclu <- melt(minclus)
  mminclu$uniques <- as.list(mcols(faas[mminclu$value])$X)
  fastclu <- lapply(split(mminclu$uniques, f = mminclu$L1), unlist)
  rm(minclus)
  rm(mminclu)
  names(fastclu) <- NULL
  reps <- sapply(fastclu, '[', 1) # Representative
  mcols(faas[reps])[[1]] <- List(fastclu) #
  faas <- faas[reps]
  return(faas)
}




#' @importFrom reshape2 melt
#' @importFrom textreuse minhash_generator
minhash_clust_k4 <- function(faas, n = 16L, cutoff = (n-1L)/(n+1L), verbose = TRUE){

  minhash <- minhash_generator(n = n)

  if (verbose){
    message('Computing minhashes.')
  }
  lp <- lapply(as.character(faas), function(x){
    y <- strsplit(x, '')[[1]]
    # k = 4
    minhash(paste0(y[c(T, F, F, F)],
                   y[c(F, T, F, F)] ,
                   y[c(F, F, T, F)],
                   y[c(F, F, F, T)] ))
  })

  if (verbose){
    message('Generating minhashes-to-sequences hash table.')
  }
  ints <- melt(lp)
  ee <- list2env(split(ints$L1, ints$value), hash = TRUE)
  rm(ints)

  cutoff <- cutoff
  res <- list()
  i <- 1L
  ln <- length(lp)
  left <- ln

  if (verbose){
    message('Clustering highly similar sequences using MinHash algorithm.')
  }
  while (left) {
    x <- lp[[1]]
    an <- unique(unlist(lapply(as.character(x), function(y) ee[[y]]), use.names = FALSE))
    an <- an[which(an %in% names(lp))]

    if (length(an)){
      mch <- an[-1L]
      hsh <- .subset(lp, mch)
      jdis <- jaccard_similarity_V(x, hsh)
      nms <- c(names(lp[1]), names(which(jdis>=cutoff)))
      res[[i]] <- nms
    }else{
      nms <- names(lp[1])
      res[[i]] <- nms
    }

    lp <- lp[!names(lp)%in%nms]
    i <- i + 1L
    left <- length(lp)

    if (verbose){
      pcnt <- round(left * 100 / ln)
      message(paste(pcnt, '%'), appendLF = FALSE)
      message(paste0(rep('\r', nchar(pcnt) + 3), collapse = ''),appendLF = FALSE)
    }

  }

  return(res)

}



# Optimization
jaccard_similarity <- function(a, b){
  length(intersect(a, b)) / length(unique.default(c(a, b)))
}

# Vectorizacion
jaccard_similarity_V <- Vectorize(jaccard_similarity, 'b', SIMPLIFY = FALSE)

