#' @name buildPanMatrix
#' @title Build the pan-matrix
#' @description Builds the pan-matrix from the list of orthologue clusters.
#' @param pangenome A \code{list} of protein clusters or a \code{pangenome}
#' object as returned by \link{pangenome}.
#' @param type One of: "none" (no panmatrix), "binary" (presence/absence
#' matrix) (DEFAULT), "nparalog" (number of paralogues for each organism and
#' orthologous group), "representative" (just one sequence of each organism, in
#' case it has more than one paralogue), or "allgenes" (write all sequences,
#' paralogues included).
#' @details The function detects if \code{pangenome} argument is a list
#' (as used internally by \code{pangenome()} function, or returned  by
#' \link{getClusterGenes} function with params \code{annot=F} and
#' \code{sequence=F}) or a \code{pangenome} object as returned by the cited
#' function (the most probable case for the user).
#' @return A pan-genome \code{data.frame}. Each organism is a column, and
#' each orthologue cluster is a row. If \code{binary=T}, the values (0 or 1)
#  indicates absence or presence of the protein gene in the organism,
#' respectively. Otherwise each cell contain the number of paralogues of each
#' orthologue group for each organism.
#' @author Ignacio Ferres
#' @export
buildPanMatrix <- function(pangenome, type='binary'){

  if(class(pangenome)=='pangenome' & type=='binary'){
    return(pangenome$panmatrix)
  }else if(class(pangenome)=='pangenome' & all(type!=c('none','binary'))){
    pangenome$clusters -> final.clusters
  }else if(class(pangenome)=='list'){
    pangenome -> final.clusters
  }else{
    stop('"pangenome" must be an object of class "pangenome" or "list".')
  }

  type <- match.arg(type,
                    choices = c('binary',
                                'nparalog',
                                'representative',
                                'allgenes'))

  sort(unique(sapply(unlist(final.clusters),function(x){
    strsplit(x,';')[[1]][1]
  }))) -> cls
  m <- matrix(nrow = length(final.clusters),
              ncol = length(cls))
  rownames(m) <- names(final.clusters)
  colnames(m) <- cls

  if(type=='binary'){
      for (i in seq_along(final.clusters)){
        as.numeric(cls%in%sapply(final.clusters[[i]],function(x){
          strsplit(x,';')[[1]][1]
        })) -> m[i,]
      }

  }else if (type=='nparalog'){

    for (i in seq_along(final.clusters)){
      c(final.clusters[[i]],attr(final.clusters[[i]],'paralogues')) -> alls
      factor(do.call(rbind,strsplit(alls,';'))[,1],levels = cls) -> this
      as.vector(table(this)) -> m[i,]
    }

  }else if(type=='representative'){

    for (i in seq_along(final.clusters)){
      final.clusters[[i]] -> alls
      factor(do.call(rbind,strsplit(alls,';'))[,1,drop=F],levels = cls) -> this
      which(table(this)>0) -> hav
      sapply(names(hav),function(x){
        grep(paste0(x,';'),alls,value = T)
        },simplify = T) -> gens
      sapply(gens,function(x){
        do.call(rbind,strsplit(x,';'))[,2]
        },simplify = T) -> gens
      m[i,hav] <- as.vector(unlist(gens))
    }

  }else if(type=='allgenes'){

    for (i in seq_along(final.clusters)){
      c(final.clusters[[i]],attr(final.clusters[[i]],'paralogues')) -> alls
      factor(do.call(rbind,strsplit(alls,';'))[,1,drop=F],levels = cls) -> this
      which(table(this)>0) -> hav
      sapply(names(hav),function(x){
        grep(paste0(x,';'),alls,value = T)
        },simplify = F) -> gens
      sapply(gens,function(x){
        paste0(do.call(rbind,strsplit(x,';'))[,2,drop=F],collapse = ';')
        },simplify = T) -> gens
      m[i,hav] <- as.vector(gens)
    }


  }

  return(as.data.frame(m))
}


#' @name getClusterGenes
#' @title Extract the genes of a set of clusters
#' @author Ignacio Ferres
#' @description Extract a list of genes from a set of clusters in the
#' \code{pangenome} object.
#' @param x A \code{pangenome} object.
#' @param clustNames A \code{vector} of cluster names.
#' @param annot \code{logical} Should original gene annotation be retrieved?
# #' @param sequence \code{logical} Should sequences as \code{SeqFastaAA} objects
#' of \link{seqinr} package be retrieved?
#' @details Something
#' @return A \code{list} of ethier internal gene names, original annotation
#' names, or the sequences as class \code{SeqFastaAA} from the clusters names
#' passed in \code{clustNames} argument.
#' @export
getClusterGenes <- function(x,clustNames=c(),annot=F){
  if(class(x)!='pangenome') {
    stop('x must be an object of class "pangenome" (PANDORA package).')
  }

  if(is.null(clustNames)){
    clustNames <- rownames(x$panmatrix)
  }

  unlist(lapply(clustNames,function(y){
    x$clusters[y]
  }),recursive = F) -> a
  if (!annot) {
    attr(a,'Annot') <- FALSE
    a
  }else{
    lapply(a,function(y){
      sapply(y,function(z){
        attr(x$cds[z][[1]],'Annot')
      })
    }) -> a
    attr(a,'Annot') <- TRUE
    a
  }
}


#' @name getCoreClusters
#' @title Extract core cluster names
#' @author Ignacio Ferres
#' @description Extract core cluster names at a certain level of permissiveness
#' @param x A \code{pangenome} object.
#' @param level \code{numeric} A number between 0.9 and 1.
#' @details Extracts the core cluster names at a certain level of
#' permissiveness.
#' \code{level} from 1 to 0.99 is the strict definition of core-genome, which
#' is the set of genes shared by the 100%-99% of the organisms.
#'
#' A \code{level}=0.95 extracts the softcore clusters, i.e. those clusters that
#' contain genes shared by a 95% of the organisms.
#'
#' This level can be adjusted to up to 0.9, which is only recomended in case of
#' working with very distant species or in case of bad sequence qualities so a
#' relatively big proportion of genes (proteins) may be absent or truncated in
#' the original genome sequence.
#' @return A \code{vector} of the core genes cluster names at the specified
#' level permissiveness.
#' @export
getCoreClusters <- function(x,level=1L){
  if(class(x)!='pangenome') {
    stop('x must be an object of class "pangenome" (PANDORA package).')
  }
  if(level>1|level<0.9){
    stop('level must be a number between 0.9 and 1. Recomended: >=0.95 .')
  }

  names(which(rowSums(x$panmatrix)>=(ncol(x$panmatrix)*level)))
}


# Deprecated
# #' @name getPanMatrix
# #' @title Extract pan-matrix from pangenome object
# #' @author Ignacio Ferres
# #' @description Extracts the binary pan-matrix from a \code{pangenome} object.
# #' @param x A \code{pangenome} object.
# #' @details A presence/absence \code{data.frame} is extracted from the
# #' \code{pangenome} object. Rows are cluster names and columns organism names.
# #' @return A binary \code{data.frame} with rows as cluster names and columns
# #' as organism names.
# #' @export
# getPanMatrix <- function(x){
#   if(class(x)!='pangenome') {
#     stop('x must be an object of class "pangenome" (PANDORA package).')
#   }
#   x$panmatrix
# }


#' #' @name plotRarefaction
#' #' @title Plot pangenome rarefaction curves
#' #' @author Ignacio Ferres
#' #' @description Plot core-genome and pan-genome rarefaction curves.
#' #' @param x A \code{pangenome} object.
#' #' @param nsamp \code{integer} The number of random samples with no replace of
#' #' the lallaa...
#' #' @details Both the number of shared genes and the total number of genes as
#' #' a function of the number of organisms sequencially added are plotted. For
#' #' each new genome added, a sample of \code{nsamp} (\code{default} 10) genomes
#' #' are evaluated with no replace.
#' #'
#' #' A nice \code{ggplot2} plot is drawn. A smooth shadow representing the mean
#' #' plus or minus a constant times the standard deviation is also ploted.
#' #' Future versions will allow more customization.
#' #'
#' #' Also 2 matrices are invisibly returned, the first for the core rarefaction
#' #' curve, and the second for the pangenome rarefaction curve. Each cell is the
#' #' count of either core or pan genes for the ith sample (rows) of jth organism
#' #' (columns) added.
#' #' @return A \code{list} of 2 \code{nsamp}*# of organisms \code{matrix} is
#' #' returned.
#' #' @importFrom ggplot2 ggplot aes xlab ylab geom_point stat_summary
#' #' @importFrom reshape2 melt
#' #' @export
#' plotRarefaction <- function(x,nsamp=10){
#'
#'   seq(1,ncol(x$panmatrix),1)->br
#'
#'   corev<-matrix(nrow = nsamp,ncol = length(br))->panev
#'   rownames(corev)<-paste("sample",1:nsamp,sep = "")
#'   colnames(corev)<-br
#'   rownames(panev)<-paste("sample",1:nsamp,sep = "")
#'   colnames(panev)<-br
#'
#'   for (b in 1:length(br)){
#'     for (i in 1:nsamp){
#'
#'       x$panmatrix[,as.vector(sample(colnames(x$panmatrix),br[b],replace = F)),drop=F]->mu
#'       if(length(which(rowSums(mu)==0))>0){
#'         mu[-which(rowSums(mu)==0),,drop=F]->mu
#'       }
#'       nrow(mu)->panev[i,b]
#'       length(which(apply(mu,1,function(x){all(x>0)})))->corev[i,b]
#'     }
#'   }
#'   #require(reshape2)
#'   melt(corev)->mcorev
#'   "Core-genome"->mcorev$set
#'   melt(panev)->mpanev
#'   "Pan-genome"->mpanev$set
#'   rbind(mcorev,mpanev)->dat
#'
#'   ggplot(data=dat,aes(x=Var2,y=value,color=set)) +
#'     xlab("# of genomes") + ylab("# of genes") +
#'     geom_point(alpha=0.5) +
#'     stat_summary(fun.data="mean_sdl",geom="smooth",
#'                  alpha=0.25)
#'
#'   invisible(list(corev,panev))
#' }
#'

#' @name writeFastaClusters
#' @title Write Cluster Sequences in Fasta Format
#' @description Write cluster's cds in fasta format.
#' @param x A \code{pangenome} object.
#' @param clustNames The names of the clusters (OGXXXX...) to be written.
#' @param ffns The path to 'all.ffn'.
#' @param outdir Where to output the files.
#' @param paralogues \code{logical}. If paralogues should be also written (if
#' present).
#' @details One file per cluster is written. In each file, one sequence from
#' each organism is written at the top. If \code{paralogues=TRUE}, and some of
#' the sequences do contain paralogues, those sequences are written at the end
#' of the file.
#' @return Fasta formated files.
#' @author Ignacio Ferres
#' @importFrom seqinr read.fasta write.fasta
#' @export
writeFastaClusters <- function(x,
                               clustNames=c(),
                               ffns='',
                               outdir='',
                               paralogues=T){

  if(class(x)!='pangenome') {
    stop('x must be an object of class "pangenome" (PANDORA package).')
  }

  if(is.null(clustNames)){
    clustNames <- names(x$clusters)
  }

  normalizePath(ffns) -> ffndir
  seqinr::read.fasta(paste0(ffndir,'all.ffn'),
                     seqtype = 'DNA',
                     as.string = T) -> seqs

  paste0(normalizePath(outdir),'/') -> outdir

  for (i in clustNames){
    x$clusters[clustNames[i]] -> a

    if(paralogues){
      c(a,attr(x$clusters[clustNames[i]],'paralogues')) -> a
    }

    seqinr::write.fasta(seqs[a],
                        names = a,
                        file.out = paste0(outdir,clustNames[i],'.fasta'))

  }
}


#' @name coreAlign
#' @title Align Core-Genes
#' @description Align core-genes and outputs a concatenation of the aligned
#' genes. Core-genes are specified at a certain level of permissiveness.
#' @param x A \code{pangenome} object.
#' @param ffns A vector of nucleotide fasta files.
#' @param level \code{numeric} between 1 and 0.9.
#' @param n_threads \code{integer}. The number of threads to use.
#' @param file.out \code{character}. The name of the output file.
#' @details This function performs an alignment of each orthologous group and
#' concatenates them into a "super-gene". This is specially suitable to perform
#' subsecuent phylogenetic analysis.
#' @return A fasta file with the concatenated core-genes alignment.
#' @author Ignacio Ferres
#' @importFrom seqinr write.fasta read.fasta
#' @export
coreAlign <- function(x,
                      ffns=character(),
                      level=1,
                      n_threads=1L,
                      file.out=''){

  if(class(x)!='pangenome'){
    stop('"x" must be an object of class "pangenome".')
  }
  if(file.out==''){
    stop('You must provide a output file name.')
  }

  if(length(which(rowSums(x$panmatrix)>=round(level * ncol(x$panmatrix))))==0){
    stop('No core-genes for the pangenome at the specified level.')
  }

  cat('Getting core-genes.\n')
  getCoreClusters(x,level=level) -> clusters
  x$clusters[clusters] -> sid


  mclapply(ffns,function(x){
    read.fasta(x,seqtype = 'DNA',as.string = T)
  },mc.cores = n_threads) -> seqs
  unlist(seqs,recursive = F) -> seqs

  lapply(sid,function(y){seqs[y]}) -> sqs

  colnames(x$panmatrix) -> orgs

  #Align.
  #Save in list the alignments as matrices.
  cat('Aligning.\n')
  al <- list()
  pb<- txtProgressBar(min = 0,max = length(sqs),style = 3)
  for (i in 1:length(sqs)){
    align(rf = sqs[[i]],n_threads = n_threads,mxit1000 = TRUE) -> a

    #If some organism/s doesn't have a core gene...
    if(nrow(a)!=length(orgs)){
      orgs[which(!orgs%in%do.call(rbind,strsplit(rownames(a),';'))[,1])] -> nw
      matrix('-',nrow = length(nw),ncol = ncol(a)) -> nwm
      rownames(nwm) <- nw
      rbind(a,nwm) -> a
    }

    al[[i]] <- a
    setTxtProgressBar(pb,i)
  }
  close(pb)

  #Make sure alignment are concatenated in the correct order.
  cat('Concatenating.\n')
  lapply(orgs,function(y){
    sapply(al,function(z){z[grep(paste0(y,';'),rownames(z),fixed=T),]})
  }) -> coral
  names(coral) <- orgs

  #Concatenate
  lapply(coral,function(y){unlist(y,use.names = F)}) -> o

  #Write output
  paste('Writing output at',file.out,'..') -> p
  cat(p)
  write.fasta(o,names = names(o),file.out = file.out)
  cat(' DONE!\n')
}













#' @name coreAlign2
#' @title Align Core-Genes
#' @description Align core-genes and outputs a concatenation of the aligned
#' genes. Core-genes are specified at a certain level of permissiveness.
#' @param x A \code{pangenome} object.
#' @param ffns A vector of nucleotide fasta files.
#' @param level \code{numeric} between 1 and 0.9.
#' @param n_threads \code{integer}. The number of threads to use.
#' @param file.out \code{character}. The name of the output file.
#' @details This function performs an alignment of each orthologous group and
#' concatenates them into a "super-gene". This is specially suitable to perform
#' subsecuent phylogenetic analysis.
#' @return A fasta file with the concatenated core-genes alignment.
#' @author Ignacio Ferres
#' @importFrom seqinr write.fasta read.fasta
#' @export
coreAlign2 <- function(x,
                      ffns=character(),
                      level=1,
                      n_threads=1L,
                      file.out=''){

  if(class(x)!='pangenome'){
    stop('"x" must be an object of class "pangenome".')
  }
  if(file.out==''){
    stop('You must provide a output file name.')
  }

  if(length(which(rowSums(x$panmatrix)>=round(level * ncol(x$panmatrix))))==0){
    stop('No core-genes for the pangenome at the specified level.')
  }

  cat('Getting core-genes.\n')
  getCoreClusters(x,level=level) -> clusters
  x$clusters[clusters] -> sid


  mclapply(ffns,function(x){
    read.fasta(x,seqtype = 'DNA',as.string = T)
  },mc.cores = n_threads) -> seqs
  unlist(seqs,recursive = F) -> seqs

  lapply(sid,function(y){seqs[y]}) -> seqs

  colnames(x$panmatrix) -> orgs

  #Align.
  #Save in list the alignments as matrices.
  cat('Aligning.\n')
  al <- list()
  pb<- txtProgressBar(min = 0,max = length(sqs),style = 3)
  for (i in 1:length(seqs)){
    align(rf = seqs[[i]],n_threads = n_threads,mxit1000 = TRUE) -> a

    #If some organism/s doesn't have a core gene...
    if(nrow(a)!=length(orgs)){
      orgs[which(!orgs%in%do.call(rbind,strsplit(rownames(a),';'))[,1])] -> nw
      matrix('-',nrow = length(nw),ncol = ncol(a)) -> nwm
      rownames(nwm) <- paste0(nw,';')
      rbind(a,nwm) -> a
    }


    a[sapply(paste0(orgs,';'),grep,rownames(a)),] -> a
    tmp <- tempfile(fileext = '.ali')
    as.alignment(a) -> a
    write.fasta(sequences = as.list(a$seq),names = a$nam,file.out = tmp)
    al[[i]] <- tmp
    setTxtProgressBar(pb,i)
  }
  close(pb)


  unlist(al) -> al

  mclapply(orgs,function(x){

    li <- list()
    for (i in 1:length(al)){
      read.alignment(al[i],format = 'fasta') -> ral
      ral$seq[[grep(paste0(x,';'),ral$nam)]] -> li[[i]]
    }
    he <- paste0('>',x,'\n')
    writeLines(c(he,paste0(unlist(li),collapse = '')))

  })


  #Concatenate

  #Write output
  paste('Writing output at',file.out,'..') -> p
  cat(p)
  write.fasta(o,names = names(o),file.out = file.out)
  cat(' DONE!\n')
}




