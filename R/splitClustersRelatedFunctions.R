#' @name splitPreClusters
#' @title Split Pre Clusters
#' @description Takes coarse clusters and split them using a tree prunning
#' algorithm.
#' @param fastas Amino acid sequences
#' @param pre.clusters A \code{list} of clusters.
#' @param accuAli \code{logical}. Use accurate alignment?
#' @param n_threads \code{integer} The number of cpus to use.
#' @return A curated \code{list} of clusters.
#' @author Ignacio Ferres
splitPreClusters <- function(fastas,
                             pre.clusters,
                             accuAli = FALSE,
                             n_threads = 1L) {

  ind.withparalogues <- which(sapply(pre.clusters, function(x) {
    checkIfParalogues(x)
  }))

  cat("Splitting pre-clusters..")
  splitedClusters <- mclapply(ind.withparalogues, function(x) {
    splitClusters(clstr = fastas[pre.clusters[[x]]], accuAli = accuAli)
  }, mc.cores = n_threads, mc.preschedule = FALSE)

  cat(" DONE!\n")


  # Merge clusters (Clusters splited in previous step with clusters which doesn't
  # contained paralogues')
  cat("Merging clusters.. ")
  clusters <- c(unlist(splitedClusters, recursive = F),
                pre.clusters[-ind.withparalogues])

  # No domain and no phmmer hit (names):
  hu <- names(fastas)[which(!names(fastas) %in% unlist(clusters))]


  # Merge clusters (Clusters merged in previous step with orphans sequences which doesn't
  # have any pfam domain or similarity with any other sequence[so coudn't be captured by
  # phmmer thus neither by mcl]).
  clusters <- c(clusters, as.list(hu))

  if (!is.null(names(clusters))) {
    pfamstr <- strsplit(names(clusters), ";")

    for (i in 1:length(clusters)) {
      attr(clusters[[i]], "pfamStr") <- pfamstr[[i]]
    }
  }

  # Set cluster names
  names(clusters) <- setClusterNames(clusters = clusters)

  # Clusters with recent paralogues:
  rcnt <- which(sapply(clusters, function(x) {
    any(duplicated(do.call(rbind, strsplit(x, ";"))[, 1]))
  }))
  # Use the first one as representative and pass the other to attr(,'paralogues')
  for (i in rcnt) {
    repr <- clusters[[i]][1]
    toattr <- clusters[[i]][-1]
    clusters[[i]] <- repr
    attr(clusters[[i]], "paralogues") <- toattr
  }
  cat("DONE!\n")

  return(clusters)
}


#' @name splitClusters
#' @title Split clusters
#' @description Takes coarse protein groups clustered by ethier domain
#' structure or by phmmer comparison plus MCL, and split them in true
#' orthologues groups.
#' @param clstr A \code{list} of \code{SeqFastaAA} objects.
#' @param accuAli \code{logical}. Perform accurate alignment? If number of
#' sequences is grater than 400, it is set to FALSE.
#' @details A multiple alignment is performed by mafft aligner. The alignment
#' is trimmed and a p-distance matrix is calculated. A neighbour joining (NJ)
#' tree is inferred, and midpoint rooted.
#'
#' Recent paralogues are detected by identifying those nodes which only contain
#' tips from the same organism, and just the one with the smaller distance to
#' its neighbours is kept.
#'
#' Maximum subtrees which contain at most one tip of each organism are
#' considered true orthologues.
#' @return A \code{list} of true orthologues protein clusters.
#' @author Ignacio Ferres
#' @importFrom phangorn midpoint
#' @importFrom ape njs subtrees
splitClusters <- function(clstr,
                          accuAli = FALSE){
  #Align
  if (length(clstr)>400){
    accuAli <- FALSE
  }
  align(rf = clstr,type = 'AA',accu = accuAli) -> al
  #Trim alignment
  round(nrow(al)*0.8) -> trtr
  which(apply(al,2,function(x){length(which(x=='-'))})>=trtr) -> trimal
  if (length(trimal)!=0){
    al[,-trimal] -> al
  }

  #Calculate p distances
  pdist.aa(al) -> d
  #Infer NJ tree
  try(phangorn::midpoint(ape::njs(d)),silent = T) -> tree
  if(class(tree)=='try-error'){
    ape::njs(d) -> tree
  }


  #Extract true orthologues
  ape::subtrees(tree) -> l
  whichNodesHaveRecentParalogues(subtrees = l) -> w
  if(!is.null(w)){
    #..and then remove them leaving just one (the closer one).
    rmvRecentParalogues(njtree = tree,w = w,d = d) -> tree.trim
    ape::subtrees(tree.trim$tree) -> l
  }
  whichSubtreesAreTrueOrthologues(subtrees = l) -> tru.or

  maxTrueOGs(tru.or = tru.or,l = l) -> tru.or2
  sapply(l[which(tru.or2)],function(x){x$tip.label},simplify = F) -> fin

  if(!is.null(w)){
    #Some singletones may have out of 'fin' because they doesn't form a
    # subtree, so..
    sapply(names(tree.trim$paralogs),function(x){any(grep(x,fin))}) -> sap
    if(length(which(!sap))>0){
      as.list(names(which(!sap))) -> sng
      append(fin,sng) -> fin
    }
    # attr(fin,'paralogues') <- tree.trim$paralogs
    sapply(names(tree.trim$paralogs),function(x){grep(x,fin)}) -> where
    for (i in 1:length(fin)){
      if (any(i==where)){
        which(names(tree.trim$paralogs)%in%fin[[i]]) -> prl
        attr(fin[[i]],'paralogues') <- as.vector(unlist(tree.trim$paralogs[prl]))
      }
    }

  }
  # cat('.')
  fin
}


#' @name align
#' @title Multiple alignment with mafft.
#' @description Takes a list of \code{SeqFastaAA} objects from \link{seqinr}
#' package and performs a multiple alignment with mafft aligner. Outputs a
#' \code{matrix} with columns as positions in the alignment and rows as
#' sequences.
#' @param rf A \code{list} of \code{SeqFastaAA} objects.
#' @param type \code{character}. One of 'AA' or 'DNA'.
#' @param n_threads \code{integer}. The number of threads to use.
#' @param accu \code{logical}. Accurate alignment?
#' @param mxit1000 \code{logical}. Set \code{mafft} parameter
#' \code{--maxiterate 1000} .
#' @return A \code{matrix} with columns as positions in the alignment and rows
#' as sequences.
#' @author Ignacio Ferres
#' @references Katoh, Misawa, Kuma, Miyata (2002). MAFFT: a novel method for
#' rapid multiple sequence alignment based on fast Fourier transform. \emph{Nucleic
#' Acids Res} \strong{30}:3059-3066.
#' @importFrom seqinr write.fasta
align<-function(rf,type='AA',n_threads=1L,accu=TRUE,mxit1000=TRUE){
  tempfile(fileext = ".faa")->tmp
  write.fasta(rf,names = names(rf),file.out = tmp)

  type <- match.arg(type,c('AA','DNA'))

  if(type=='AA'){

    if(accu){
      args <- c("--quiet --amino --bl 30","--globalpair")
      if(mxit1000){
        args <- c(args,"--maxiterate 1000")
      }

    }else{
      if(length(rf)<2000){
        args <- c("--quiet --amino --bl 30")
      }else{
        args <- c('--quiet --amino --bl 30 --retree 1 --maxiterate 0')
      }

    }

  }else{
    if(accu){
      args <- c("--quiet --nuc","--globalpair")
      if(mxit1000){
        args <- c(args,"--maxiterate 1000")
      }

    }else{
      if(length(rf)<2000){
        args <- c("--quiet --nuc")
      }else{
        args <- c('--quiet --nuc --retree 1 --maxiterate 0')
      }
    }
  }


  if (n_threads>1){
    args <- c(args,"--thread",n_threads)
  }
  args <- c(args,tmp)


  #RUN MAFFT
  if (length(rf) <= 100){

    system2(command = "mafft",
            args = args,
            stdout = T)->rl
  }else{

    tmpo <-tempfile()
    paste0('mafft ',paste0(args,collapse = ' '),' > ',tmpo) -> run
    system(run)
    readLines(tmpo) -> rl
    file.remove(tmpo)

  }

  #Transform alignment in matrix
  grep(">",rl)->gp
  cbind(gp+1,c(gp[-1]-1,length(rl)))->coorsq
  strsplit(apply(coorsq,1,function(x){
    paste(toupper(rl[x[1]:x[2]]),collapse = '')
  }),split = '')->sqs
  names(sqs)<-sub(">","",rl[gp])
  do.call(rbind,sqs)->m
  file.remove(tmp)
  #apply(m,1:2,charToRaw)->m
  return(m)
}

#' @name pdist.aa
#' @title Calculate p distance from a multiple alignment
#' @description Calculate p distance from a multiple alignment.
#' @param x \code{matrix}. A multiple alignment as passed by
#' \link{align}.
#' @param scaled \code{logical}. Return scaled by alignment length?
#' @return A \code{dist} object.
#' @note Taken and adapted from \code{\link[ape]{dist.aa}}, from ape package by
#' Paradis et al.
#' @author Ignacio Ferres
pdist.aa <- function (x, scaled = TRUE) {
  if (any(x%in%letters)){
    x <- toupper(x)
  }
  x <- apply(x,1:2,charToRaw)
  n <- nrow(x)
  d <- numeric(n * (n - 1)/2)
  X <- charToRaw("-")
  k <- 0L
  for (i in 1:(n - 1)) {
    a <- x[i, ]
    for (j in (i + 1):n) {
      b <- x[j, ]
      del <- a == X & b == X
      p <- length(b <- b[!del])
      tmp <- sum(a[!del] != b)
      k <- k + 1L
      d[k] <- if (!scaled){
        tmp
      }else{ tmp/p }
    }
  }
  attr(d, "Size") <- n
  attr(d, "Labels") <- rownames(x)
  attr(d, "Diag") <- attr(d, "Upper") <- FALSE
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}

#' @name whichNodeHaveRecentParalogues
#' @title Which nodes have recent paralogues
#' @description Detect which nodes contain only tips from the same organism.
#' @param subtrees A \code{list} of trees of class "phylo".
#' @return A \code{vector} of nodes that contain recen paralogues.
#' @author Ignacio Ferres
whichNodesHaveRecentParalogues <- function(subtrees){
  unlist(lapply(subtrees,function(x){
    sapply(x$tip.label,function(y){strsplit(y,';')[[1]][1]}) -> s
    if(length(unique(s))==1 & length(s)>1){
      x$name
    }
  }))
}


#' @name rmvRecentParalogues
#' @title Remove recent paralogues
#' @description If recent paralogues are detected, just the one with the closer
#' distance to the tips which are descendants of the immediate previous node,
#' is kept.
#' @param njtree A tree of class \code{phylo}.
#' @param w A \code{vector} of nodes which contain recent paralogues, as passed
#' by \link{whichNodesHaveRecentParalogues}.
#' @param d A \code{dist} object with the pairwise p-distance between proteins
#' as passed by \link{pdist.aa}.
#' @return A \code{list} with two elements: A tree of class \code{phylo} with
#' paralogues removed, and a \code{list} of paralogues removed.
#' @author Ignacio Ferres
#' @importFrom phangorn Descendants Ancestors
#' @importFrom ape drop.tip
rmvRecentParalogues <- function(njtree,w,d){
  #Removes subtrees and nodes contained in other subtrees
  sapply(w,function(x){
    phangorn::Descendants(njtree,x,type='tips')
  }) -> tps
  if(length(w)>1){
    rmv <- c()
    for (i in 1:length(tps)){
      any(sapply(tps[-i],function(x){all(tps[[i]]%in%x)})) -> rr
      rmv <- c(rmv,rr)
    }

    if(any(rmv==TRUE)){
      w[-which(rmv)] -> w
      tps[-which(rmv)] -> tps
    }
  }


  as.matrix(d) -> d2

  #Removes recent paralogues until just one is kept in each case
  keep <- c()
  paralogs <- list()
  rem <- c()
  for (i in 1:length(w)){
    njtree$tip.label[tps[[i]]] -> gns
    #Ancestral node
    phangorn::Ancestors(njtree,w[i],'parent') -> nd
    #Descendantes of the ancestral node
    phangorn::Descendants(njtree,nd,'tips') -> dsc
    #Genes of the closest node (descendants from the parent node)
    njtree$tip.label[dsc[[1]][-which(dsc[[1]]%in%tps[[i]])]] -> bros
    #Kept the one with the smaller distance between it and their brothers
    which.min(apply(d2[gns,bros,drop=F],1,sum)) -> k
    paralogs[[i]] <- gns[-k]
    keep <- c(keep,k)
    rem <- c(rem,gns[-k])
  }
  #list of paralogues. The chosen one is keept as element list names
  names(paralogs) <- names(keep)
  #Descard the rest of the tips (recent paralogues)
  ape::drop.tip(njtree,rem) -> njtree
  list(njtree,paralogs) -> out
  names(out) <- c('tree','paralogs')
  return(out)
}


#' @name whichSubtreesAreTrueOrthologues
#' @title Which subtrees are true orthologues
#' @description Which subtrees are true orthologues, i.e. which nodes contain
#' at most only one tip of each organism.
#' @param subtrees A \code{list} of trees of class \code{phylo}.
#' @return A \code{vector} of \code{logical} values indicating if subtrees
#' contain at most one tip of each organism.
#' @author Ignacio Ferres
whichSubtreesAreTrueOrthologues<-function(subtrees){
  res<-c()
  for (i in 1:length(subtrees)){
    sapply(subtrees[[i]]$tip.label,function(x){strsplit(x,';')[[1]][1]}) -> p
    table(p) -> t
    if(all(t==1)){
      res<-c(res,TRUE)
    }else{
      res<-c(res,FALSE)
    }
  }
  res
}


#' @name maxTrueOGs
#' @title Maximum subtrees which contains only true orthologues
#' @description Detect maximum subtrees which contains only true orthologues.
#' @param tru.or A \code{vector} of logical values.
#' @param l A \code{list} of trees of class \code{phylo}.
#' @return A logical \code{vector} indicating if true orthologue subtree is
#' maximum or not.
#' @author Ignacio Ferres
maxTrueOGs <- function(tru.or,l){
  res2<-tru.or
  for(i in 1:length(tru.or)){
    if(tru.or[i]==FALSE){
      next
    }else{
      l[[i]]$tip.label -> tps
      res2[i] <- FALSE
      sapply(l[which(res2)],function(x){x$tip.label},simplify = F)->s
      any(sapply(s,function(x){all(tps%in%x)})) -> r
      if(!r){
        res2[i] <- TRUE
      }
    }
  }
  res2
}
