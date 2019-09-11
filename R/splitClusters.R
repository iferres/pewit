#' @importFrom parallel mclapply
#' @importFrom S4Vectors split
#' @importFrom reshape2 melt
splitPreClusters <- function(fastas, n_threads, sep, minhash_split = FALSE, verbose){

  preclusters <- split(fastas, mcols(fastas)$Arch)

  splited <- mclapply(preclusters,
                      splitCluster,
                      sep = sep,
                      minhash_split = minhash_split,
                      verbose = verbose,
                      mc.cores = n_threads)
  attr(splited, 'varname') <- 'Arch'

  splited_df <- melt(splited, id.vars=c('Gene', 'NODE'))

  splited_df$CLUSTER <- paste(splited_df$Arch, splited_df$NODE, sep = sep)

  splited_df$CLUSTER <- as.integer(as.factor(splited_df$CLUSTER))

  mcols(fastas[splited_df$Gene])$Cluster <- splited_df$CLUSTER

  if (verbose){
    lpre <- length(preclusters)
    lpos <- length(unique(splited_df$CLUSTER))
    mssg <- paste(lpre, 'pre-clusters splited into', lpos, 'final groups of orthologues.')
    message(mssg)
  }

  fastas
}


#' @import DECIPHER
#' @importFrom phangorn midpoint Descendants Ancestors
#' @importFrom ape bionjs cophenetic.phylo drop.tip
#' @importFrom reshape2 melt
#' @importFrom stats as.dist
splitCluster <- function(x, sep, minhash_split= FALSE, verbose = TRUE){

  ntips <- length(x)
  fc <- as.integer(as.factor(as.character(x)))
  anydup <- any(table(fc)>1)
  if (anydup){
    xl <- split(names(x), fc)
    x2 <- x
    x <- x[sapply(xl, '[', 1)]
  }

  mcls <- mcols(x)

  if (verbose){
    arch <- strsplit(mcls$Arch[1], ';')[[1]]
    mssg <- paste('   Resolving precluster tree: ',
                  paste(paste0('[', arch, ']'), collapse = ' '),
                  '(',ntips, 'tips )')
    message(mssg)
  }

  if (length(unique(mcls$organism)) == 1){
    # Future versions should be able to resolve cases where there are only one
    # organism but there's also a filogenetic structure, like:
    # ((A; A);(A; A));
    # For now, in that cases all sequences are assigned to the same NODE.
    df <- data.frame(Gene = names(x), NODE = 'NODE_1', row.names = NULL)

  }else if (length(x)>2){

    if (any(duplicated(mcls$organism))){

      # Align, compute distance, compute nj tree.
      if (!minhash_split){
        ali <- AlignTranslation(x, verbose = FALSE)
        dm <- DistanceMatrix(ali, verbose = FALSE)
        d <- as.dist(dm)
      }else{
        d <- minhash_dist(x, k = 16)
      }

      tree <- midpoint(bionjs(d))
      tiplab <- tree$tip.label
      trorgs <- sapply(strsplit(tiplab, sep), '[', 1)

      nodes_recpar <- c()
      for (i in seq_along(tiplab)){
        nd <- NULL
        anc <- Ancestors(tree, i, type = 'parent')
        dec <- Descendants(tree, node = anc)[[1]]
        while(all(trorgs[i]==trorgs[dec])){
          nd <- anc
          anc <- Ancestors(tree, anc, type = 'parent')
          dec <- Descendants(tree, node = anc)[[1]]
        }
        nodes_recpar <- c(nodes_recpar, nd)
      }
      nodes_recpar <- unique(nodes_recpar)

      par_list <- Descendants(tree, node = nodes_recpar)
      names(par_list) <- lapply(par_list, function(x) tiplab[x[1]])
      par_list <- lapply(par_list, function(x) tiplab[x[-1]])
      paralogs <- unlist(par_list)
      tree2 <- tree
      for (i in seq_along(paralogs)){
        tree2 <- drop.tip(tree2, tip = paralogs[i])
      }

      # Compute descendants tips for every node
      alldes <- Descendants(tree2, type = 'tips')
      names(alldes) <- paste('NODE', seq_along(alldes), sep = '_')

      # Which nodes have repeated organisms?
      rep_orgs <- which(vapply(alldes, function(y){
        orgs <- vapply(strsplit(tree2$tip.label[y], sep), '[', 1,
                       FUN.VALUE = NA_character_)
        any(duplicated(orgs))
      }, FUN.VALUE = NA))

      # If more than 1 clade (length rep_orgs > 0)...
      if (length(rep_orgs)){
        # Remove nodes with repeated organisms (they are not orthologues)
        ogs <- alldes[-rep_orgs]

        # Compute containment matrix (if each node is contained in another node)
        in_mat <- sapply(ogs, function(y){
          sapply(ogs, function(z){
            all(y%in%z)
          }, USE.NAMES = TRUE)
        }, USE.NAMES = TRUE)

        # Sum. Those which sum == 1 are those which only contain themselves. Those
        # are true orthologues. Extract tip labels.
        true_ogs <- lapply(ogs[colSums(in_mat)==1], function(y){
          tree2$tip.label[y]
        })
        attr(true_ogs, 'varname') <- 'NODE'

        # ... else (if only a single clan exists)..
      }else{

        true_ogs <- list(NODE_1 = tree2$tip.label)
        attr(true_ogs, 'varname') <- 'NODE'

      }

      # Melt list into long-formatted data.frame
      df <- melt(true_ogs, value.name = 'Gene')
      df[] <- lapply(df, as.character)

      # If recent paralogues exists, add them to result
      if (length(paralogs)){

        dfpar <- lapply(names(par_list), function(x){
          NODE <- df$NODE[which(df$Gene == x)]
          data.frame(Gene = par_list[[x]], NODE = NODE)
        })

        dfpar <- do.call(rbind, dfpar)
        df <- rbind(df, dfpar)
      }

      #else (aren't duplicated organisms, so they are a single orthologous group)
    }else{
      df <- data.frame(Gene = names(x), NODE = 'NODE_1', row.names = NULL)
    }


    # else (length of precluster is 1 or 2)...
  }else{
    df <- data.frame(Gene = names(x), NODE = 'NODE_1', row.names = NULL)
  }

  if (anydup){
    xl <- xl[which(sapply(xl, length)>1)]
    dfdup <- lapply(xl, function(x){
      NODE <- df$NODE[which(df$Gene == x[1])]
      data.frame(Gene = x[-1], NODE = NODE)
    })
    dfdup <- do.call(rbind, dfdup)
    df <- rbind(df, dfdup)
  }

  # attr(splited, 'varname') <- 'Arch'
  return(df)
}

# #' @importFrom ape subtrees
# recent_par_nodes <- function(phy, sep = '___'){
#   subtrees <- subtrees(phy)
#   names(subtrees) <- as.character(sapply(subtrees, '[[', 'name'))
#   rcnt <- vapply(subtrees, function(y){
#     orgs <- sapply(strsplit(y$tip.label, sep), '[', 1)
#     uorgs <- unique(orgs)
#     if (length(orgs) > 1 & length(uorgs) == 1){
#       TRUE
#     }else{
#       FALSE
#     }
#   }, FUN.VALUE = NA)
#   as.integer(names(which(rcnt)))
# }
#
#
#
#
# tip2node_len <- function(phy, tip, node){
#   edge <- phy$edge
#   elen <- phy$edge.length
#   elen[elen<0] <- abs(elen[elen<0])
#   wtip <- which(phy$tip.label==tip)
#   len <- c()
#   t1 <- wtip
#   rw <- which(edge[, 2]==t1)[1]
#   t2 <- edge[rw, 1]
#   while(t2!=node){
#     len <- c(len, elen[rw])
#     #next...
#     edge <- edge[-rw, ]
#     elen <- elen[-rw]
#     t1 <- t2
#     rw <- which(edge[, 2]==t1)[1]
#     t2 <- edge[rw, 1]
#   }
#   sum(len)
# }
#
# Vtip2node_len <- Vectorize(tip2node_len, vectorize.args = 'tip')


#' @importFrom textreuse minhash_generator
minhash_dist <- function(x, k){

  xc <- as.character(x)
  # dup <- any(duplicated(x))
  # if (dup){
  #   xcn <- split(names(xc), as.integer(as.factor(xc)))
  #   xc <- xc[sapply(xcn, '[', 1)]
  # }

  minhash_fun <- minhash_generator(200)
  mhs <- lapply(xc, compute_minhash, k = k, minhash_fun = minhash_fun)
  n <- length(xc)

  if (n>70){

    d <- lsh_candidates(mhs)
    y <- 1L
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        if (.subset2(d, y) == 0L){
          d[y] <- jaccard_dissimilarity(.subset2(mhs,i), .subset2(mhs,j))
        }
        y <- y + 1L
      }
    }

  }else{

    d <- vector('numeric', length = (n * n-1L) / 2)
    y <- 1L
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        d[y] <- jaccard_dissimilarity(.subset2(mhs,i), .subset2(mhs,j))
        y <- y + 1L
      }
    }

  }

  attr(d, 'class') <- 'dist'
  attr(d, 'Labels') <- names(mhs)
  attr(d, 'Size') <- n
  attr(d, 'Diag') <- FALSE
  attr(d, 'Upper') <- FALSE

  # if (dup){
  #   m <- as.matrix(d)
  #   ln <- unname(sapply(xcn, length, USE.NAMES = F))
  #   du <- which(ln>1, useNames = F)
  #
  #   d <- as.dist(m)
  # }

  return(d)
}

compute_minhash <- function(x, minhash_fun = NULL, k =16){
  kmers <- compute_kmers(x, k)
  minhash_fun(kmers)
}


#' @importFrom parallel splitIndices
#' @importFrom digest digest
lsh_candidates <- function(mhs){
  n <- 200L
  ln <- length(mhs)
  # Determine sensitivity
  prop_diff <- length(unique(unlist(mhs))) / (n * ln)
  b <- ifelse(prop_diff > .5, 50L, 40L)
  bix <- splitIndices(n, b)
  Mhs <- t(do.call(rbind, mhs))
  bands <- lapply(bix, function(x) Mhs[x, ])
  cna <- dimnames(Mhs)[[2]]
  mres <- matrix(1L,
                 nrow = ln,
                 ncol = ln,
                 dimnames = list(cna, cna))
  lapply(bands, function(y){
    ap <- apply(y, 2, digest)
    un <- unname(split(cna, ap))
    lapply(un, function(x) {
      mres[x, x] <<- 0L
    })
  })

  mres[upper.tri(mres, diag = T)] <- NA
  as.dist(mres)
}

jaccard_dissimilarity <- function(a, b){
  1 - length(intersect2(a, b)) / length(unique.default(c(a, b)))
}


intersect2 <- function(x, y){
  m <- match(x, y, 0L)
  unique.default(.subset(y, m))
}
