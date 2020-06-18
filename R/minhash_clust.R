# #' @importFrom reshape2 melt
#' @importFrom S4Vectors elementNROWS mcols mcols<- List
fast_clust <- function(faas, verbose = TRUE){
  # Take largest as representative sequence
  faas <- faas[order(elementNROWS(faas), decreasing = TRUE)]
  minclus <- minhash_clust_k4(faas, verbose = verbose)
  # mminclu <- melt(minclus)
  lns <- sapply(minclus, length)
  mminclu <- data.frame(L1 = rep(seq_along(minclus), lns))
  mminclu$value <- unlist(minclus)
  mminclu$uniques <- as.list(mcols(faas[mminclu$value])$X)
  fastclu <- lapply(split(mminclu$uniques, f = mminclu$L1), unlist)
  rm(minclus)
  rm(mminclu)
  names(fastclu) <- NULL
  reps <- sapply(fastclu, '[', 1) # Representative
  mcols(faas[reps])[["X"]] <- List(fastclu) #
  faas <- faas[reps]
  return(faas)
}



#' @importFrom parallel splitIndices
#' @importFrom digest digest
#' @importFrom reshape2 melt
#' @importFrom textreuse minhash_generator
#' @importFrom stats setNames
minhash_clust_k4 <- function(faas, n = 16L, cutoff = (n-1L)/(n+1L), verbose = TRUE){

  minhash <- minhash_generator(n = n)

  if (verbose){
    message('Computing 4-mers and minhashes.')
  }
  lp <- list2env(lapply(as.character(faas), function(x){
    minhash(compute_kmers(x, k=4L))
  }), hash = TRUE)

  bix <- splitIndices(n, 2L)

  bands <- lapply(lp, function(x) {
    vapply(bix, function(y){
      digest(.subset(x, y))
    }, FUN.VALUE = NA_character_)
  })

  # bands <- do.call(rbind, bands)
  ml <- melt(do.call(rbind, bands))
  dml <- dim(ml)
  mlx_1 <- seq.int(1L, dml[[1]]/2)
  mlx_2 <- seq.int((dml[[1]]/2) + 1L, dml[[1]])

  ee1 <- list2env(split(as.character(ml$Var1[mlx_1]), ml$value[mlx_1], drop = TRUE), hash = TRUE)
  ee2 <- list2env(split(as.character(ml$Var1[mlx_2]), ml$value[mlx_2], drop = TRUE), hash = TRUE)
  rm(list = c('ml', 'mlx_1', 'mlx_2'))

  res <- list()
  i <- 1L
  # left <- dim(bands)[[1]]
  left <- length(bands)
  oln <- left
  while (left>0){
    # bn <- dimnames(bands)[[1]]
    # x <- bands[1L, ]
    bn <- names(bands)
    x <- bands[[1]]

    candi <- unique(c(ee1[[x[1]]], ee2[[x[2]]]))
    candix <- match(candi, bn, nomatch = NA)
    nas <- is.na(candix)
    candi <- candi[!nas]
    candix <- candix[!nas]

    ln <- length(candi)
    if (ln>1){
      r1 <- .subset2(lp, candi[1])
      r2 <- lapply(setNames(candi[-1], candi[-1]), function(x) lp[[x]])
      jsim <- jaccard_similarity_V(r1, r2)
      hits <- c(candi[1], names(which(jsim>=cutoff)))
      rmix <- candix[match(hits, candi)]
      # bands <- bands[-rmix, , drop=FALSE]
      # bands <- bands[-rmix]
      bands <- .subset(bands, -rmix)
      res[[i]] <- hits
    }else{
      hits <- candi
      # bands <- bands[-1L, ,drop=FALSE]
      bands <- bands[-1]
      res[[i]] <- hits
    }

    # left <- dim(bands)[[1]]
    left <- length(bands)
    i <- i + 1L

    if (verbose){
      pcnt2 <- if (exists('pcnt')) pcnt else 101
      pcnt <- round(left * 100 / oln)
      if (pcnt<pcnt2){
        message(paste(pcnt, '%'), appendLF = FALSE)
        message(paste0(rep('\r', nchar(pcnt) + 6), collapse = ''), appendLF = FALSE)
      }
    }
  }

  return(res)

}



# Optimization
jaccard_similarity <- function(a, b){
  length(intersect(a, b)) / length(unique.default(c(a, b)))
}

# Vectorizacion
jaccard_similarity_V <- Vectorize(jaccard_similarity, 'b', SIMPLIFY = TRUE)



#' @importFrom parallel mclapply
#' @importFrom textreuse minhash_generator
#' @importFrom S4Vectors elementNROWS
#' @importFrom stats hclust cutree
#' @importFrom utils capture.output
clust_orgs <- function(faas_orgs, n_threads = 1, k = 4L, n = 512L, h = 0.3, verbose = TRUE){

  minhash <- minhash_generator(n = n)

  if (verbose) message("Computing whole proteome minhashes")

  lp <- mclapply(faas_orgs, function(x){
    ln <- elementNROWS(x)
    ln <- ln[-length(ln)]
    krm <- unlist(lapply(cumsum(ln), function(y) seq(y - (k - 2L) , y, 1)))
    cr <- paste(as.character(x), collapse = '')
    km <- compute_kmers(cr, k)[-krm]
    minhash(km)
  }, mc.preschedule = TRUE, mc.cores = n_threads)

  if (verbose) message("Computing whole proteome distance estimates")
  djaccard <- dist_jaccard(lp)

  co <- capture.output(summary(djaccard))
  if (verbose){
    message(co[1])
    message(co[2])
  }

  hc <- hclust(djaccard, method = "ward.D2")
  ct <- cutree(hc, h = h)
  spl <- split(names(ct), ct)
  if (verbose){
    lno <- length(faas_orgs)
    lnc <- length(spl)
    if (lnc==1){
      message("Great! All proteomes are closely related (estimated above 70% 4-mer shared)")
    }else{
      mssg <- paste("Found", lnc, "clusters (estimated above 70% 4-mer shared) from", lno, "organisms.")
      message(mssg)
    }
  }

  return(spl)
}

dist_jaccard <- function(x){
  n <- length(x)
  d <- vector('integer', length = (n * n/2) - n/2)
  y <- 1L
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      d[y] <- jaccard_similarity(x[[i]], x[[j]])
      y <- y + 1L
    }
  }
  d <- 1 - d
  attr(d, 'class') <- 'dist'
  attr(d, 'Labels') <- names(x)
  attr(d, 'Size') <- n
  attr(d, 'Diag') <- FALSE
  attr(d, 'Upper') <- FALSE
  d
}
