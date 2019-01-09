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
    final.clusters <- pangenome$clusters
  }else if(class(pangenome)=='list'){
    final.clusters <- pangenome
  }else{
    stop('"pangenome" must be an object of class "pangenome" or "list".')
  }

  type <- match.arg(type,
                    choices = c('binary',
                                'nparalog',
                                'representative',
                                'allgenes'))

  cls <- sort(unique(sapply(unlist(final.clusters),function(x){
    strsplit(x,';')[[1]][1]
  })))
  m <- matrix(nrow = length(final.clusters),
              ncol = length(cls))
  rownames(m) <- names(final.clusters)
  colnames(m) <- cls

  if(type=='binary'){
      for (i in seq_along(final.clusters)){
        m[i, ] <- as.numeric(cls%in%sapply(final.clusters[[i]],function(x){
          strsplit(x,';')[[1]][1]
        }))
      }

  }else if (type=='nparalog'){

    for (i in seq_along(final.clusters)){
      alls <- c(final.clusters[[i]],attr(final.clusters[[i]],'paralogues'))
      this <- factor(do.call(rbind,strsplit(alls,';'))[,1],levels = cls)
      m[i, ] <- as.vector(table(this))
    }

  }else if(type=='representative'){

    for (i in seq_along(final.clusters)){
      alls <- final.clusters[[i]]
      this <- factor(do.call(rbind,strsplit(alls,';'))[,1,drop=F],levels = cls)
      hav <- which(table(this)>0)
      gens <- sapply(names(hav),function(x){
        grep(paste0(x,';'),alls,value = T)
        },simplify = T)
      gens <- sapply(gens,function(x){
        do.call(rbind,strsplit(x,';'))[,2]
        },simplify = T)
      m[i,hav] <- as.vector(unlist(gens))
    }

  }else if(type=='allgenes'){

    for (i in seq_along(final.clusters)){
      alls <- c(final.clusters[[i]],attr(final.clusters[[i]],'paralogues'))
      this <- factor(do.call(rbind,strsplit(alls,';'))[,1,drop=F],levels = cls)
      hav <- which(table(this)>0)
      gens <- sapply(names(hav),function(x){
        grep(paste0(x,';'),alls,value = T)
        },simplify = F)
      gens <- sapply(gens,function(x){
        paste0(do.call(rbind,strsplit(x,';'))[,2,drop=F],collapse = ';')
        },simplify = T)
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
# #' of \link{seqinr} package be retrieved?
# #' @details Something
#' @return A \code{list} of ethier internal gene names, original annotation
#' names, or the sequences as class \code{SeqFastaAA} from the clusters names
#' passed in \code{clustNames} argument.
#' @export
getClusterGenes <- function(x,clustNames=c(),annot=F){
  if(class(x)!='pangenome') {
    stop('x must be an object of class "pangenome" (pewit package).')
  }

  if(is.null(clustNames)){
    clustNames <- rownames(x$panmatrix)
  }

   a <- unlist(lapply(clustNames,function(y){
    x$clusters[y]
  }),recursive = F)
  if (!annot) {
    attr(a,'Annot') <- FALSE
    a
  }else{
    a <- lapply(a,function(y){
      sapply(y,function(z){
        attr(x$cds[z][[1]],'Annot')
      })
    })
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
    stop('x must be an object of class "pangenome" (pewit package).')
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
#     stop('x must be an object of class "pangenome" (pewit package).')
#   }
#   x$panmatrix
# }


#' @name plotRarefaction
#' @title Plot Pangenome Rarefaction Curves
#' @description Plot tipical pangenome and coregenome curves.
#' @param x An object of class \code{pangenome}.
#' @param nsamp \code{integer} The number of genomes to sample on each stage.
#' @param plot \code{logical} If \code{TRUE}, then a plot is produced. If not,
#' just the result matrices are returned invisibly.
#' @param xlab \code{character} X axis label.
#' @param ylab \code{character} Y axis label.
#' @param y.mar A \code{vector} of length 2, indicating a proportion of the
#' minimum and maximum values to take as bottom and top margins, respectively.
#' @param las \code{numeric} Style of axis labels. See \code{?par}.
#' @param pg.col Pangenome points and curve color.
#' @param cg.col Coregenome points and curve color.
#' @param pt.pch An \code{integer} specifying a symbol. See \code{?par}.
#' @param pt.cex \code{numeric} The point size.
#' @param pt.alpha A factor modifying the opacity alpha of the points.
#' @param li.lwd \code{integer} The line width.
#' @param shadow \code{logical} Plot a ribbon showing a function of points
#' quantiles (quantiles to plot specifyied by \code{shadow.quant}).
#' @param shadow.quant A \code{vector} of length 2 in {0, 1, 2, 3, 4, 5},
#' specifying the quantiles to plot.
#' @param shadow.col The ribbon color.
#' @param leg \code{logical} Plot a legend?
#' @param legend A \code{vector} of length 2 with the legend names to plot.
#' @param leg.pos The position of the legend.
#' @param leg.cex \code{numeric} The legend size.
#' @param leg.pt.cex \code{numeric} The legend point size.
#' @param leg.pch An \code{integer} specifying the legend symbols.
#' @param leg.bty The legend box type. See \code{?legend}.
#' @param ... Arguments to be passed to \code{plot()}
#' @details Both the number of shared genes and the total number of genes as
#' a function of the number of organisms sequencially added are plotted. For
#' each new genome added, a sample of \code{nsamp} (\code{default} 10) genomes
#' are evaluated with no replace.
#'
#' A scatter plot is drawn.
#' Future versions will allow more customization.
#'
#' Also 2 matrices are invisibly returned, the first for the core rarefaction
#' curve, and the second for the pangenome rarefaction curve. Each cell is the
#' count of either core or pan genes for the ith sample (rows) of jth organism
#' (columns) added.
#' @return A \code{list} of 2 \code{nsamp}*# of organisms \code{matrix} is
#' returned.
#' @importFrom graphics par plot box axis polygon points lines
#' @importFrom stats quantile predict loess
#' @importFrom grDevices adjustcolor
#' @export
plotRarefaction <- function(x,
                            nsamp = 10,
                            plot = TRUE,
                            xlab = "# of genomes",
                            ylab = "# of genes",
                            y.mar = c(0.85, 1.15),
                            las = 0,
                            pg.col = "#E64B35FF",
                            cg.col = "#4DBBD5FF",
                            pt.pch = 20,
                            pt.cex = 1,
                            pt.alpha = 0.4,
                            li.lwd = 3,
                            shadow = TRUE,
                            shadow.quant=c(2,4),
                            shadow.col = "grey",
                            leg = TRUE,
                            legend = c("Pangenome", "Coregenome"),
                            leg.pos = "topleft",
                            leg.cex = par("cex"),
                            leg.pt.cex = par("cex"),
                            leg.pch = 22,
                            leg.bty = 'n',
                            ...){

  if (class(x) != "pangenome"){
    stop('Object is not of class "pangenome" (pewit).')
  }

  n <- ncol(x$panmatrix)
  br <- seq(1, n, 1)

  corev <- matrix(nrow = nsamp, ncol = n) -> panev
  rownames(corev) <- paste("sample", 1:nsamp, sep = "")
  colnames(corev) <- br
  rownames(panev) <- paste("sample",1:nsamp,sep = "")
  colnames(panev) <- br

  for (b in 1:n){
    for (i in 1:nsamp){

      mu <- x$panmatrix[,as.vector(sample(colnames(x$panmatrix),br[b],replace = F)),drop=F]
      if(length(which(rowSums(mu)==0))>0){
        mu <- mu[-which(rowSums(mu)==0),,drop=F]
      }
      panev[i, b] <- nrow(mu)
      corev[i, b] <- length(which(apply(mu,1,function(x){all(x>0)})))
    }
  }

  vl <- unlist(c(corev,panev))
  ymar <- c(min(vl)*y.mar[1], max(vl)*y.mar[2])

  if(plot){

    op <- par(no.readonly = TRUE)
    on.exit(op)

    dev.hold()
    plot(NA,
         type='n',
         xlim = c(1, ncol(x$panmatrix)),
         ylim = ymar,
         axes = FALSE,
         xlab = xlab,
         ylab = ylab,
         ...)

    box()
    axis(1,at = 1:n, las = las)
    axis(2, las = las)

    if (shadow){

      pg.quant <- apply(panev, 2, quantile)
      pg.mi <- predict(loess(pg.quant[shadow.quant[1],] ~ c(1:n)))
      pg.ma <- predict(loess(pg.quant[shadow.quant[2],] ~ c(1:n)))
      # pg.mi[n] <- pg.quant[1,n]

      cg.quant <- apply(corev, 2, quantile)
      cg.mi <- predict(loess(cg.quant[shadow.quant[1],] ~ c(1:n)))
      cg.ma <- predict(loess(cg.quant[shadow.quant[2],] ~ c(1:n)))
      # cg.mi[n] <- cg.quant[1,n]

      polygon(x = c(1:n, (n-1):1),
              y = c(pg.mi[-1], rev(pg.ma)),
              col = adjustcolor(shadow.col, alpha.f = 0.3),
              border = F)
      polygon(x = c(1:n, (n-1):1),
              y = c(cg.mi, rev(cg.ma)[-1]),
              col = adjustcolor(shadow.col, alpha.f = 0.3),
              border = F)
    }


    apply(panev,
          1,
          points,
          pch = pt.pch,
          col = adjustcolor(pg.col, alpha.f = pt.alpha),
          cex = pt.cex)

    apply(corev,
          1,
          points,
          pch = pt.pch,
          col = adjustcolor(cg.col, alpha.f = pt.alpha),
          cex = pt.cex)

    p.me <- loess(apply(panev, 2, mean) ~ c(1:n))
    c.me <- loess(apply(corev, 2, mean) ~ c(1:n))

    lines(predict(p.me), col = pg.col, lwd = li.lwd)
    lines(predict(c.me), col = cg.col, lwd = li.lwd)

    if (leg){
      legend(x = leg.pos,
             legend = legend,
             cex = leg.cex,
             pch = leg.pch,
             pt.cex = leg.pt.cex,
             pt.bg = c(pg.col, cg.col),
             bty = leg.bty,
             x.intersp = 0.7,
             y.intersp = 0.7)
    }

    dev.flush()
  }

  o <- list(corev, panev)
  # class(o) <- "Rarefaction"

  invisible(o)
}


#' @name writeFastaClusters
#' @title Write Cluster Sequences in Fasta Format
#' @description Write cluster's cds in fasta format.
#' @param x A \code{pangenome} object.
#' @param clustNames The names of the clusters (OGXXXX...) to be written.
#' @param ffns The fasta files (as a vector of paths, as returned by
#' \code{list.files} with options \code{pattern = 'ffn$'} and
#' \code{full.names = TRUE}).
#' @param outdir Where to output the files.
#' @param type \code{character}. One of "\code{default}" (one sequence per
#' genome on each fasta file), "\code{all}" (including paralogues), or
#' "\code{representative}" (one sequence per fasta file, a representative of
#' the cluster).
#' @param n_threads The number of cpus to use. Used when reading fasta files.
# #' @details One file per cluster is written. In each file, one sequence from
# #' each organism is written at the top. If \code{paralogues=TRUE}, and some of
# #' the sequences do contain paralogues, those sequences are written at the end
# #' of the file.
#' @return Fasta formated files.
#' @author Ignacio Ferres
#' @importFrom seqinr read.fasta write.fasta
#' @importFrom parallel mclapply
#' @export
writeFastaClusters <- function(x,
                               clustNames=c(),
                               ffns=c(),
                               outdir='',
                               type='default',
                               # paralogues=T,
                               n_threads=1){

  if(class(x)!='pangenome') {
    stop('x must be an object of class "pangenome" (pewit package).')
  }

  type <- match.arg(type, c('default','all','representative'))

  if(is.null(clustNames)){
    clustNames <- names(x$clusters)
  }

  ffndir <- normalizePath(ffns)
  # seqinr::read.fasta(paste0(ffndir,'all.ffn'),
  #                    seqtype = 'DNA',
  #                    as.string = T) -> seqs
  cat('Reading fasta files..')
  seqs <- parallel::mclapply(ffndir,function(x){
    seqinr::read.fasta(x,seqtype = 'DNA',as.string = T)
  },mc.cores = n_threads)
  seqs <- unlist(seqs,recursive = F)
  cat(' DONE!\n')

  out.dir <- paste0(normalizePath(outdir),'/')

  cat('Writing..')
  for (i in seq_along(clustNames)){
    a <- x$clusters[[clustNames[i]]]

    # if(paralogues){
    #   c(a,attr(x$clusters[clustNames[i]],'paralogues')) -> a
    # }

    if (type=='default'){
      seqinr::write.fasta(seqs[a],
                          names = a,
                          file.out = paste0(outdir,clustNames[i],'.fasta'))
    }else if(type=='all'){
      a <- c(a,attr(x$clusters[clustNames[i]],'paraligues'))
      seqinr::write.fasta(seqs[a],
                          names = a,
                          file.out = paste0(outdir,clustNames[i],'.fasta'))
    }else if(type=='representative'){
      a <- a[1]
      seqinr::write.fasta(seqs[a],
                          names = a,
                          file.out = paste0(outdir,clustNames[i],'.fasta'))
    }
  }
  cat(' DONE!\n')
}


#' @name coreAlign
#' @title Align Core-Genes
#' @description Align core-genes and outputs a concatenation of the aligned
#' genes. Core-genes are specified at a certain level of permissiveness.
#' @param x A \code{pangenome} object.
#' @param ffns A vector of nucleotide fasta files.
#' @param level \code{numeric} between 1 and 0.9.
#' @param accu \code{logical}. Accurate alignment?
#' @param n_threads \code{integer}. The number of threads to use.
#' @param file.out \code{character}. The name of the output file.
#' @details This function performs an alignment of each orthologous group and
#' concatenates them into a "super-gene". This is specially suitable to perform
#' subsecuent phylogenetic analysis.
#' @return A fasta file with the concatenated core-genes alignment.
#' @author Ignacio Ferres
#' @importFrom seqinr write.fasta read.fasta read.alignment
#' @importFrom ape as.alignment
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
coreAlign <- function(x,
                      ffns=c(),
                      level=1,
                      accu=TRUE,
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
  clusters <- getCoreClusters(x,level=level)
  sid <- x$clusters[clusters]


  seqs <- mclapply(ffns,function(x){
    seqinr::read.fasta(x,seqtype = 'DNA',as.string = T)
  },mc.cores = n_threads)
  seqs <- unlist(seqs,recursive = F)

  seqs <- mclapply(sid,function(y){
    seqs[y]
  },mc.cores = n_threads)

  orgs <- colnames(x$panmatrix)

  #Align.
  cat('Aligning.\n')
  al <- list()
  pb <- txtProgressBar(min = 0,max = length(seqs),style = 3)
  for (i in 1:length(seqs)){
    if(accu){
      a <- align(rf = seqs[[i]],
            type = 'DNA',
            n_threads = n_threads,
            accu = TRUE,
            mxit1000 = TRUE)
    }else{
      a <- align(rf = seqs[[i]],
            type = 'DNA',
            n_threads = n_threads,
            accu=FALSE)
    }


    #If some organism/s doesn't have a core gene...
    if(nrow(a)!=length(orgs)){
      nw <- orgs[which(!orgs%in%do.call(rbind,strsplit(rownames(a),';'))[,1])]
      nwm <- matrix('-',nrow = length(nw),ncol = ncol(a))
      rownames(nwm) <- paste0(nw,';')
      a <- rbind(a,nwm)
    }

    #Order
    a <- a[sapply(paste0('^',orgs,';'),grep,rownames(a)),]
    #Write
    tmp <- tempfile(fileext = '.ali')
    a <- ape::as.alignment(a)
    seqinr::write.fasta(sequences = as.list(a$seq),names = a$nam,file.out = tmp)
    al[[i]] <- tmp
    setTxtProgressBar(pb,i)
  }
  close(pb)


  al <- unlist(al)


  #Concatenate (horizontal)
  cat('Merging..')
  supergenes <- mclapply(orgs,function(x){

    li <- list()
    for (i in 1:length(al)){
      ral <- seqinr::read.alignment(al[i],format = 'fasta')
      li[[i]] <- ral$seq[[grep(paste0('^',x,';'),ral$nam)]]
    }
    he <- paste0('>',x)
    tmp <- tempfile(fileext = '.supergene')
    writeLines(c(he,paste0(unlist(li),collapse = '')),
               con = tmp,
               sep = '\n')
    tmp

  },mc.cores = n_threads)
  supergenes <- unlist(supergenes)
  cat(' DONE!\n')
  file.remove(al)

  #Concatenate (vertical). Output.
  p <- paste('Writing output at',file.out,'..')
  cat(p)
  for (i in 1:length(supergenes)){
    cat(readLines(supergenes[i]),
        file = file.out,
        sep = '\n',
        append = T)
  }
  cat(' DONE!\n')
  file.remove(supergenes)

}




