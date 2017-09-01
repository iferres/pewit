# pangenome() function and methods associated to generic functions to handle
# pangenome objects.

#' @name pangenome
#' @title Compute the pangenome of a set of related organisms.
#' @author Ignacio Ferres and Gregorio Iraola
#'
#' @description Computes the pangenome of a set of related organisms. Takes a
#' list of gff3 files, a path to Pfam HMM models and data associated to those
#' models.
#'
#' It first uses hmmsearch (HMMER 3) to look for Pfam domains in the aminoacid
#' sequences and, where it finds them, creates a domain profile. Then clusterize
#' sequences with the same domain structure in protein families.
#'
#' Those proteins without any domain identified are compared with phmmer
#' (HMMER 3) and clustered with mcl into protein families.
#'
#' Protein families obtained from both two previous steps are then splited into
#' true orthologues by a gene tree-prunning algorithm, creating a pangenome
#' profile of the group of organisms.
#' @param gffs A \code{vector} of gff3 files as retrieved by prokka (Seeman,
#' 2014).
#' @param hmmPfam \code{character} with the path to Pfam-A.hmm file. If any of
#' \code{hmmPfam} or \code{datPfam} is left \code{NULL}, then the pfam step is
#' skipped.
#' @param datPfam \code{character} with the path to Pfam-A.hmm.dat file. If any
#'  of \code{hmmPfam} or \code{datPfam} is left \code{NULL}, then the pfam step
#'  is skipped.
#' @param n_threads \code{integer}. The number of threads to use.
#' @param dir_out Name of output directory to be created.
#' @param writeFfns \code{logical}. If should write ffn sequences. If \code{alignCore}
#' is set to \code{TRUE} (see below), this option is set to TRUE.
#' @param writeFastas \code{logical}. Write fasta files with gene sequences for
#' each cluster?
#' @param pmOutfileType The type of pan-matrix you want to be written in
#' \code{out.dir} folder. One of: "none" (no panmatrix), "binary" (presence/
#' absence matrix) (DEFAULT), "nparalog" (number of paralogues for each organism
#' and orthologous group), "representative" (just one sequence of each organism,
#' in case it has more than one paralogue), or "allgenes" (write all sequences,
#' paralogues included).
#' @param alignCore \code{logical}. Align core genes and concatenate them into
#' a "super-gene" alignment? Suitable for phylogenetic analysis.
#' @param accuAli \code{logical}. Perform an accurate alignment? If \code{TRUE}
#' the alignment stage could take a while depending both on the number of
#' genomes and on the core size.
#' @param coreLevel \code{numeric}. A number between 1-0.9 which determines at
#' what proportion of presence should clusters be considered as part of the
#' core-genome.
#' @param seed \code{integer}. A seed to allow reproducible analyses.
#' @details A scan against Pfam-A database is performed and only those hits of
#' PF class 'domain' or 'family' are considered for further analysis.
#' If two or more domains of the same clan overlap, then the one with the
#' smaller e-value is kept. For those proteins which have Pfam hits, a profile
#' of each protein is determined by the consecutive sequence of Pfam domains,
#' and those proteins with the same domain structure are then clustered
#' together to form protein families.
#'
#' Phmmer is then used to compare all those proteins with no Pfam domains
#' identified and the Markov Clustering algorithm (mcl) clusters them into
#' coarse groups (i.e. protein families).
#'
#' As this package is designed to handle groups of organisms related at a
#' greater level than species (although it works perfectly fine with close
#' related clades), a paralogue split based on conserved neighbourhood was not
#' considered because this relationships get lose at distantly related
#' organisms. Instead of that it relies on the concept of orthology itself to
#' diferenciate between paralogues and true orthologues, this is, a gene tree
#' based inference.
#'
#' With each protein familiy groups which contain paralogues, an alignment
#' followed by all-vs-all distance calculation is performed. A
#' neighbour-joining tree is then generated and midpoint rooted.
#' If recent paralogues are detected (this is, nodes which all its descendants
#' belong to the same organism), then all but one tips are removed from the
#' gene-tree and saved as recent paralogues. Just the one sequence with the
#' closer distance to its neighbours is kept.
#' Finally, the gene-tree is splited in as many sub-trees as necessary to have
#' the minimum set of sub-trees with just one gene of each organism.
#'
#' In that way true-orthologue clusters are retrieved and saved in an object of
#' class "pangenome".
#' @return An object of class \code{pangenome} which basically consist on a
#' \code{list} with the first element being the list of clusters, the second
#' element the pan-matrix as a binary \code{data.frame} with columns being the
#' organisms and rows the cluster names, and a series of attributes.
# #' @references
#' @note External dependencies are HMMER3, MCL and MAFFT, without them running
#' the function will stop with an error.
#' @importFrom parallel mclapply splitIndices
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach '%dopar%'
#' @importFrom utils write.table
#' @export
pangenome<-function(gffs=c(),
                    hmmPfam=NULL,
                    datPfam=NULL,
                    n_threads=1L,
                    dir_out='out',
                    writeFfns=FALSE,
                    writeFastas=FALSE,
                    pmOutfileType='representative',
                    alignCore=TRUE,
                    accuAli=FALSE,
                    coreLevel=1,
                    seed = 5){

  if (Sys.which("hmmsearch")==""){
    stop("\n\tHMMER (v.3) is not installed. (Couldn't find 'hmmsearch' in $PATH)
         \tPlease install it before re-running pangenome().\n\n")
  }

  if (Sys.which("mcl")==""){
    stop("\n\tMCL is not installed. (Couldn't find 'mcl' in $PATH)
         \tPlease install it before re-running pangenome().\n\n")
  }

  if(Sys.which("mafft")==""){
    stop("\n\tMAFFT aligner is not installed. (Couldn't find 'mafft' in $PATH)
         \tPlease install it before re-running pangenome().\n\n")
  }

  if(any(grepl('.gff$',gffs)==F)){
    stop("Files ('gffs') must have '.gff' extension.")
  }

  if(dir.exists(dir_out)){
    stop(paste0('Directory ',dir_out,'/ already exists.'))
  }

  if(any(file.info(gffs)$isdir)){
    stop("At least one of the 'gffs' is a directory.")
  }


  #Create output directory
  dir.create(dir_out)
  outdir <- paste0(normalizePath(dir_out),'/')

  #Run on exit
  on.exit(runOnExit(outdir))

  #Type of panmatrix arg
  pmOutfileType <- match.arg(pmOutfileType,
                             c('binary',
                               'nparalog',
                               'representative',
                               'allgenes',
                               'none'))




  cat('Extracting and translating gene sequences from gff files..')

  wffns <- ifelse(writeFfns,'dna','none')
  fastas <- mclapply(gffs,function(x){
    extractSeqsFromGff3(infile = x,
                        in.path = outdir,
                        keep = 'aa',
                        write.in.path = ifelse(alignCore, 'dna', wffns))
  }, mc.cores = n_threads, mc.preschedule = FALSE)

  fastas <- unlist(fastas,recursive = F)
  cat(' DONE!\n')




  if (!any(sapply(list(hmmPfam,datPfam), is.null))){

    clus <- domainSearch(fastas = fastas,
                         hmmPfam = hmmPfam,
                         datPfam = datPfam,
                         n_threads = n_threads)
    clu_dom <- clus[[1]]
    clu_fam <- clus[[2]]
    tout <- clus[[3]]

  }else{

    cat('Either hmmPfam or/and datPfam arguments were left NULL. Pfam search will be skipped.\n')
    clu_dom <- NULL
    clu_fam <- NULL
    tout <- data.frame(Domain = character(),
                       Family = character(),
                       Motif = character(),
                       Repeat = character())

  }




  # PHMMER + MCL
  cat('Clustering orphan sequences (phmmer + mcl) ..')
  orphans <- clusterOrphans(fastas = fastas,
                            tout = tout,
                            n_threads = n_threads)
  cat(' DONE!\n')




  #Split paralogs
  pre.clusters <- c(clu_dom, clu_fam, orphans)
  clusters <- splitPreClusters(fastas = fastas,
                               pre.clusters = pre.clusters,
                               accuAli = accuAli,
                               n_threads = n_threads)



  #Pan-matrix (presence/absence)
  cat('Computing provisory binary pan-matrix..')
  panm <- buildPanMatrix(pangenome = clusters,
                         type='binary')
  si1 <- length(which(rowSums(panm)==1))
  cat(paste0(' there are currently ',si1,' singletons.\n'))




  cat('Refining..\n')
  cat('        ..realocating misassigned singletones..\n')
  clusters <- realocateSingletons(clusters = clusters,
                                 panm = panm,
                                 fastas = fastas,
                                 n_threads = n_threads,
                                 seed = seed)
  names(clusters) <- setClusterNames(clusters = clusters)
  cat('        ..computing binary pan-matrix..\n')
  panm <- buildPanMatrix(pangenome = clusters,
                         type='binary')
  si2 <- length(which(rowSums(panm)==1))
  asgn <- si1 - si2
  cat(paste0(' ..DONE! ',asgn,' singletons realocated.\n'))




  cat('Preparing output..\n')
  #panmatrix.tab
  cat('        ..writing panmatrix:')
  write.table(panm,file = paste0(outdir,'panmatrix.tab'),sep = '\t',quote = F)
  cat(paste0(' DONE! Saved at ',outdir,'panmatrix.tab\n'))

  #clusters.txt
  cat('        ..writing clusters:')
  writeClusters(outdir = outdir,
                final.clusters = clusters,
                filename = 'clusters.txt')
  cat(paste0(' DONE! Saved at ',outdir,'clusters.txt\n'))

  #paralogues.txt
  cat('        ..writing paralogues:')
  writeParalogues(outdir = outdir,final.clusters = clusters)
  cat(paste0(' DONE! Saved at ',outdir,'paralogues.txt\n'))

  #Out
  ncds <- length(fastas)
  norgs <- ncol(panm)
  nclust <- nrow(panm)
  ncore100 <- length(which(rowSums(panm)==ncol(panm)))
  ncore95 <- length(which(rowSums(panm)>=round(ncol(panm)*0.95)))
  nsingles <- length(which(rowSums(panm)==1))
  naccs <- length(which(rowSums(panm)<round(ncol(panm)*0.95)))-nsingles
  out <- list(clusters,
              panm)

  names(out) <- c('clusters','panmatrix')
  attr(out,'ncds') <- ncds
  attr(out,'norgs') <- norgs
  attr(out,'nclust') <- nclust
  attr(out, 'ncore100') <- ncore100
  attr(out, 'ncore95') <- ncore95
  attr(out, 'nsingles') <- nsingles
  attr(out, 'naccs') <- naccs
  attr(out,'output') <- outdir
  attr(out,'package') <- 'pewit'
  class(out) <- c('pangenome')

  cat('        ..saving pangenome object:')
  saveRDS(out,file = paste0(outdir,'pangenome.rds'))
  cat(paste0(' DONE! Saved at ',outdir,'pangenome.rds\n'))




  ## OTHER OPTIONS

  #Do you want to write fasta files for each cluster?
  if(writeFastas){
    cat('          ..writing clusters:')
    paste0(outdir,'clusters/') -> cludir
    dir.create(cludir)
    writeFastaClusters(x=out,
                       clustNames = c(),
                       ffns = list.files(path = outdir,
                                         pattern = 'ffn$',
                                         full.names = T),
                       outdir=cludir,
                       type = 'all',
                       n_threads = n_threads)
    cat(paste0(' DONE!'),
        attr(out,'nclust'),
        'clusters written at',
        cludir,'\n')
  }

  #Type of pan-matrix to be written on out directory (not in pangenome object)
  if (pmOutfileType!='none'){
    if(pmOutfileType=='binary'){
      write.table(panm,
                  file = paste0(outdir,'panmatrix.tab'),
                  quote = F,
                  na = '-',
                  sep = '\t')
    }else if(pmOutfileType%in%c('nparalog',
                                'representative',
                                'allgenes')){

      panm <- buildPanMatrix(out, type = pmOutfileType)
      write.table(panm,
                  file = paste0(outdir,'panmatrix.tab'),
                  quote = F,
                  na = '-',
                  sep = '\t')
    }
  }

  #Align core genome with mafft aligner?
  if(alignCore){
    if(length(which(rowSums(out$panmatrix)>=(coreLevel * ncol(out$panmatrix))))>0){
      list.files(path = outdir,pattern = '.ffn$',full.names = T) -> ffns
      coreAlign(x = out,
                ffns = ffns,
                level = coreLevel,
                accu = accuAli,
                n_threads = n_threads,
                file.out = paste0(outdir,'coreAlignment_',coreLevel,'.fasta'))
    }else{
      cat('No core-genes for the pangenome at the specified level. Try obtaining a suitable
number of "soft" core genes lowering the level through "getCoreClusters()" function.',fill=T)
      cat('Then pass the desire level to "coreAlign()" function')
    }

  }

  cat('\nFINNISH!\nThanks for using PEWIT.\n')
  out

}



#Print generic
#' @name print.pangenome
#' @title Print Generic for Pangenome Objects
#' @description Generic function to print "pangenome" class objects.
#' @param x An object of class \code{pangenome}.
#' @param ... Unused.
#' @return \code{character}. A compact summary of a pangenome object.
#' @author Ignacio Ferres
#' @export
print.pangenome <- function(x, ...){

  titl <- 'Object of class "pangenome" (pewit package).\n\n'

  norgs <- attr(x, 'norgs')
  ncds <- attr(x, 'ncds')
  nclust <- attr(x, 'nclust')
  ncore100 <- attr(x, 'ncore100')
  ncore95 <- attr(x, 'ncore95')
  nsingles <- attr(x, 'nsingles')
  naccs <- attr(x, 'naccs')

  out <- c(norgs,
           ncds,
           nclust,
           ncore100,
           ncore95,
           naccs,
           nsingles)

  li <- c('Number of organisms:               ',
          'Number of cds:                     ',
          'Number of clusters:                ',
          'Number of core-genes (100%):       ',
          'Number of soft core-genes (>=95%): ',
          'Number of accessory genes(<95%):   ',
          'Number of singletones:             ')

  pr <- paste('\t', li, out, '\n')

  cat(titl, pr, '\n')

}


#' @name summary.pangenome
#' @title Summarizing pangenome object
#' @author Ignacio Ferres
#' @description Summary method for class \code{pangenome}
#' @param object An object of class \code{pangenome}
#' @param ... Unused
#' @return A \code{matrix} with softcore, accessory and singletons numbers for
#' each genome.
#' @export
summary.pangenome <- function(object, ...){

  socore <- object$panmatrix[which(rowSums(object$panmatrix)>=(ncol(object$panmatrix)*0.95)),]
  nscore <- colSums(socore)
  ssing <- object$panmatrix[which(rowSums(object$panmatrix)==1),]
  nsing <- colSums(ssing)
  accs <- object$panmatrix[which(rowSums(object$panmatrix)<(ncol(object$panmatrix)*0.95)),]
  naccs <- colSums(accs) - nsing

  o <- data.frame(row.names=colnames(object$panmatrix))
  o$Softcore <- nscore
  o$Accessory <- naccs
  o$Singletones <- nsing

  o <- as.matrix(o)
  o <- t(o)

  return(o)
}


#' @name plot.pangenome
#' @title Plot a \code{pangenome} object.
#' @description Generic function to plot a \code{pangenome} object.
#' @param x \code{pangenome} object.
#' @param hclustfun Function used to compute the hierarchical clustering.
#' @param distfun Function used to compute the distance (dissimilarity).
#' @param edgePar A list of parameters to be passed to plot.dendrogram().
#' @param margins A \code{numeric} vector of length 2 containing the margins
#' for column and row names, repectively.
#' @param col A \code{character} vector with the colors of the presence/absence
#' plot. The first color indicate absence, and the second indicates presence.
#' @param raster \code{logical} See \code{useRaster} parameter in \code{?image}.
#' @param labRow A \code{character} vector of length equal to the number of
#' organisms in the pangenome object. By default it uses the original file
#' names without the extension (.gff).
#' @param cexRow \code{numeric} The size of row labels.
#' @param rowSideColors A \code{character} vector of length equals the number
#' of organisms present in the pangenome object. The \code{default} is \code{NULL},
#' so the rowSideColor column is not plotted.
#' @param colSideColors A \code{character} vector of length 3 specifying the
#' colors to use on the column bar. The colors correspond, in order, to
#' singletons, accessory and soft-core clusters.
#' @param colSideLabels A \code{character} vector of length 3 to specify the
#' labels to use to indicate Singletons, Accessory genes, and Soft-core genes.
#' @param cexKey \code{numeric} The size of the key labels.
#' @param ... Unused
#' @details A heatmap-like image is plotted attached with a dendrogram. The
#' grid shows the presence or absence of a gene (column), for each
#' organism (rows). The columns are ordered so the core genes appear at the
#' left, the accessory genes at the middle, and the singletones at the right
#' of the grid.
#'
#' The left dendrogram is computed based on a dissimilarity measure between
#' the vectors of presence/absence of 2 genomes at a time.
#'
#' The top bar indicates the genes (columns) assigned as a part of the soft-
#' core genome, the accessory genome, and the singletones.
#' @importFrom graphics image plot plot.new axis par layout rect
#' @importFrom grDevices dev.hold dev.flush
#' @import stats
#' @author Ignacio Ferres
#' @export
plot.pangenome <- function(x,
                           hclustfun = hclust,
                           distfun = dist,
                           edgePar = list(lwd = 1.5),
                           margins = c(5, 7),
                           col = c('white','#E64B35FF'),
                           raster = TRUE,
                           labRow = NULL,
                           cexRow = 0.5,
                           rowSideColors = NULL,
                           colSideColors = c('#4DBBD5FF','#00A087FF','#3C5488FF'),
                           colSideLabels = c('Softcore genome','Accessory','Singletones'),
                           cexKey = 1,
                           ...){

  if (class(x)!="pangenome"){
    stop('x is not of class "pangenome".')
  }

  ncore95 <- attr(x, 'ncore95')
  nsingles <- attr(x, 'nsingles')
  naccs <- attr(x, 'naccs')


  ### Distance and clustering ###
  dd <- distfun(t(x$panmatrix))
  hc <- hclustfun(dd)

  mm <- x$panmatrix[order(rowSums(x$panmatrix), decreasing = T), ]
  mm <- as.matrix(mm)
  mm <- mm[,hc$order]

  dhc <- as.dendrogram(hc)

  nr <- nrow(mm)
  nc <- ncol(mm)

  ### RowSideColors ###
  if (!is.null(rowSideColors)){
    if (!is.character(rowSideColors) || length(rowSideColors) != nc){
      stop("'rowSideColors' must be a character vector of length equal to the
           number of organisms in the pangenome")
    }else{
      rowSideColors <- rowSideColors[hc$order]
      mat <- rbind(c(5,5,4), c(3,2,1))
      wlo <- c(0.9, 0.1, 3)
      # wlo <- c(0.8, 3)
      # mim <- 0.2
    }
  }else{
    mat <- rbind(c(4,3), c(2,1))
    wlo <- c(1, 3)
    #   wlo <- c(1,3)
    #   mim <- 0L
  }


  ### Start device functions ###
  op <- par(no.readonly = T)
  on.exit(dev.flush())
  on.exit(par(op), add = TRUE)

  dev.hold()

  ### Set layout ###


  layout(mat = mat,
         widths = wlo,
         heights = c(1, 3),
         respect = F)

  ### Image ###
  par(mar = c(margins[1],
              0,
              0.5,
              margins[2]))
  image(x = 1:nr,
        y = 1:nc,
        z = mm,
        xlim = 0.5 + c(0, nr),
        ylim = 0.5 + c(0, nc),
        frame.plot=F,
        col = col,
        yaxt = 'n',
        xaxs= 'r',
        ylab = '',
        xlab = '',
        useRaster = raster,
        xpd = T)

  if(is.null(labRow)){
    labRow <- colnames(mm)
  }

  axis(4,
       1:nc,
       labels = labRow,
       las = 2,
       line = -1,
       tick = 0,
       cex.axis = cexRow)

  if (!is.null(rowSideColors)){
    ### rowSideColors ###
    par(mar = c(margins[1], 0, 0.5, 0))
    image(rbind(1L:nc), col = rowSideColors, axes = F)
  }

  ### Dendrogram ###
  par(mar = c(margins[1], 0, 0.5, 0), yaxs= 'i')
  plot(dhc,
       horiz = TRUE,
       # frame.plot = F,
       leaflab = 'n',
       edgePar = edgePar,
       axes=F)

  ### Column colors ###
  par(mar = c(0, 0, 0, margins[2]))
  plot.new()
  #position
  k1 <- c((ncore95+naccs)/nr, 0, (ncore95+naccs+nsingles)/nr, 0.2)
  k2 <- c(ncore95/nr, 0, (ncore95+naccs)/nr, 0.2)
  k3 <- c(0, 0, ncore95/nr, 0.2)
  #draw
  rect(k1[1],k1[2],k1[3],k1[4], col = colSideColors[1], border = NA)
  rect(k2[1],k2[2],k2[3],k2[4], col = colSideColors[2], border = NA)
  rect(k3[1],k3[2],k3[3],k3[4], col = colSideColors[3], border = NA)

  ### Key ###
  par(mar=c(1,0,2,7))
  plot(NA,xlim = c(0,1),ylim = c(0,1),xlab = '',ylab = '',axes = F)
  #draw
  rect(4/5, 2/3, 5/5, 3/3, col = colSideColors[1], border = NA)
  rect(4/5, 1/3, 5/5, 2/3, col = colSideColors[2], border = NA)
  rect(4/5, 0/3, 5/5, 1/3, col = colSideColors[3], border = NA)

  axis(4,
       at = (1:3+0:2)/6,
       labels = colSideLabels,
       cex.axis = cexKey,
       las = 2,
       line = -1,
       tick = 0)

}
