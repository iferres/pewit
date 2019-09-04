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
#' @importFrom phangorn midpoint Descendants
#' @importFrom ape bionjs cophenetic.phylo drop.tip
#' @importFrom reshape2 melt
#' @importFrom stats as.dist
splitCluster <- function(x, sep, minhash_split= FALSE, verbose = TRUE){

  mcls <- mcols(x)

  if (verbose){
    arch <- strsplit(mcls$Arch[1], ';')[[1]]
    mssg <- paste('   Resolving structure: ', paste(paste0('[', arch, ']'), collapse = ' '))
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
      }else{
        dm <- minhash_dist(x, k = 16, s = 1)
      }

      tree <- midpoint(bionjs(as.dist(dm)))

      # Determine recent paralogues, prune
      wh_rec <- recent_par_nodes(phy = tree, sep = sep)
      tree2 <- tree
      paralogs <- c()
      if (length(wh_rec)){
        while(length(wh_rec)>0){
          # Compute descendants tips for every node
          alldes <- Descendants(tree2, type = 'tips')
          names(alldes) <- paste('NODE', seq_along(alldes), sep = '_')
          node <- wh_rec[1]
          tps <- tree2$tip.label[alldes[[node]]]
          pars <- Vtip2node_len(tree2, tps, node)
          pmin <- which.min(pars)
          par_rm <- names(pars[-pmin])
          paralogs <- c(paralogs, par_rm)
          tree2 <- drop.tip(tree2, par_rm)
          wh_rec <- recent_par_nodes(phy = tree2, sep = sep)
        }
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
        # Cophenetic distances between tips
        coph <- cophenetic.phylo(tree)
        # Only keep non-paralogue columns
        coph <- coph[ , as.character(df$Gene), drop = FALSE]
        glist <- colnames(coph)
        # Look for min distance between paralog and all the non-paralog. Return
        # cluster name of closer tip.
        reloc <- vapply(paralogs, function(y){
          idx <- which.min(coph[y, , drop = TRUE])
          df[which(df$Gene == names(idx)), 2, drop = TRUE]
        }, FUN.VALUE = NA_character_, USE.NAMES = TRUE)
        dfpar <- data.frame(Gene = names(reloc), NODE = reloc, row.names = NULL)
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

  # attr(splited, 'varname') <- 'Arch'
  return(df)
}

#' @importFrom ape subtrees
recent_par_nodes <- function(phy, sep = '___'){
  subtrees <- subtrees(phy)
  names(subtrees) <- as.character(sapply(subtrees, '[[', 'name'))
  rcnt <- vapply(subtrees, function(y){
    orgs <- sapply(strsplit(y$tip.label, sep), '[', 1)
    uorgs <- unique(orgs)
    if (length(orgs) > 1 & length(uorgs) == 1){
      TRUE
    }else{
      FALSE
    }
  }, FUN.VALUE = NA)
  as.integer(names(which(rcnt)))
}




tip2node_len <- function(phy, tip, node){
  edge <- phy$edge
  elen <- phy$edge.length
  elen[elen<0] <- abs(elen[elen<0])
  wtip <- which(phy$tip.label==tip)
  len <- c()
  t1 <- wtip
  rw <- which(edge[, 2]==t1)[1]
  t2 <- edge[rw, 1]
  while(t2!=node){
    len <- c(len, elen[rw])
    #next...
    edge <- edge[-rw, ]
    elen <- elen[-rw]
    t1 <- t2
    rw <- which(edge[, 2]==t1)[1]
    t2 <- edge[rw, 1]
  }
  sum(len)
}

Vtip2node_len <- Vectorize(tip2node_len, vectorize.args = 'tip')


#' @importFrom textreuse minhash_generator
#' @importFrom S4Vectors elementNROWS
minhash_dist <- function(x, k, s){

  minhash <- minhash_generator(200)
  mhs <- mapply(compute_minhash,
                x,
                length = elementNROWS(x),
                MoreArgs = list(minhash_fun = minhash,
                                k = k,
                                s = s,
                                type = type),
                SIMPLIFY = F)

  n <- length(x)
  d <- vector('numeric', length = (n * n/2) - n/2)
  y <- 1L
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      d[y] <- jaccard_dissimilarity(mhs[[i]], mhs[[j]])
      y <- y + 1L
    }
  }

  attr(d, 'class') <- 'dist'
  attr(d, 'Labels') <- names(mhs)
  attr(d, 'Size') <- n
  attr(d, 'Diag') <- FALSE
  attr(d, 'Upper') <- FALSE
  d
}

#' @importFrom Biostrings DNAStringSet AAStringSet
compute_minhash <- function(x, minhash_fun = NULL,length, k, s){
  dss <- DNAStringSet(x, seq(1L, length - k + 1L, s), seq(k, length, s))
  minhash_fun(as.character(dss))
}


jaccard_dissimilarity <- function(a, b){
  1 - length(intersect(a, b)) / length(unique.default(c(a, b)))
}

