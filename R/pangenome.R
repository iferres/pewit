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
#' @param hmmPfam \code{character} with the path to Pfam-A.hmm file.
#' @param datPfam \code{character} with the path to Pfam-A.hmm.dat file.
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
                    hmmPfam=character(),
                    datPfam=character(),
                    n_threads=1L,
                    dir_out='out',
                    writeFfns=FALSE,
                    writeFastas=FALSE,
                    pmOutfileType='representative',
                    alignCore=TRUE,
                    accuAli=FALSE,
                    coreLevel=1){

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


  #Create output directory
  dir.create(dir_out)
  paste0(normalizePath(dir_out),'/') -> outdir

  #Run on exit
  on.exit(runOnExit(outdir))

  #Type of panmatrix arg
  pmOutfileType <- match.arg(pmOutfileType,
                             c('binary',
                               'nparalog',
                               'representative',
                               'allgenes',
                               'none'))

  #Process Pfam-A.dat
  cat('\n\nProcessing Pfam-A.hmm.dat..')
  ref<-processPfam_A_Dat(datPfam = datPfam,
                         n_threads = n_threads)
  cat(' DONE!\n')

  #deprecated:
#   #Format fasta headers
#   fastas<-formatFastaHeaders(prots = prots,n_threads = n_threads)
#   fastas<-unlist(lapply(fastas,function(x){x[[1]]}),recursive = F)
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

  cat('Preparing HMMSEARCH.\n')
  # Split indices and write fastas to distribute among threads with hmmscan
  temps<-splitAndWriteFastas(fastas = fastas,
                             n_threads = n_threads)

  #Deprecated:
  # #Run hmmscan (HMMER)
  # cat('Running HMMSCAN against Pfam-A database (this can take a while)..')
  # registerDoParallel(cores = n_threads)
  # hmm.temps<-foreach(i=seq_along(temps),.inorder = F) %dopar% {
  #   tempfile(pattern = "tmpo",tmpdir = tempdir(),fileext = ".tab")->dmblout
  #   paste("hmmscan -o /dev/null --domtblout ",dmblout,
  #         " --noali --cut_ga --cpu 0 ",
  #         hmmPfam," ",temps[i],sep = "")->pfm
  #   system(pfm)
  #   dmblout
  # }
  # cat(' DONE!\n')

  #Index hmm if not yet
  if(any(!file.exists(paste0(hmmPfam,c('.h3f','.h3i','.h3m','.h3p'))))){
    cat('Preparing Pfam-A.hmm files for hmmscan search.\n')
    paste0('hmmpress ',hmmPfam) -> hmmpress
    system(hmmpress)
  }

  #Run hmmsearch (HMMER)
  cat('Running HMMSEARCH against Pfam-A database (this can take a while)..')
  mclapply(temps, function(x){

    runHmmsearch(fasta = x,
                 pfam = hmmPfam,
                 n_threads = 0L)

  }, mc.cores = n_threads) -> hmm.temps
  cat(' DONE!\n')
  # file.remove(temps)

  #Load hmmscan output and process
  cat('Processing hmmsearch output (resolving overlapping Pfam hits and building protein families profiles)..')
  unlist(hmm.temps)->pout

  registerDoParallel(cores = n_threads)
  tout<-foreach(i=seq_along(pout),.combine =rbind,.inorder = T)%dopar%{
    processHmmsearch(pout=pout[i],ref=ref)
  }
  cat(' DONE!\n')
  file.remove(pout)

  #Clustering by domain structure only
  cat('Clustering sequences per domain structure..')
  split(x = rownames(tout),
        f = tout$Domain)[-1]->clu_dom
  cat(' DONE!\n')

  #Cluster sequences without domains asigned, by family
  cat('Clustering sequences per family..')
  split(x = rownames(tout[which(tout$Domain==''),]),
        f = tout$Family[which(tout$Domain=='')])[-1]->clu_fam
  cat(' DONE!\n')

  # PHMMER + MCL
  cat('Clustering orphan sequences (phmmer + mcl) ..')
  clusterOrphans(fastas = fastas,
                 tout = tout,
                 n_threads = n_threads) -> orphans
  cat(' DONE!\n')

  #Split paralogs
  pre.clusters <- c(clu_dom, clu_fam, orphans)

  which(sapply(pre.clusters,function(x){
    checkIfParalogues(x)
    })) -> ind.withparalogues

  cat('Splitting pre-clusters..')
  mclapply(ind.withparalogues,function(x){
    # splitClusters(clstr = lapply(fastas[pre.clusters[[x]]],memDecompress,'gzip',TRUE))
    splitClusters(clstr = fastas[pre.clusters[[x]]], accuAli = accuAli)
  },mc.cores=n_threads,mc.preschedule = FALSE) -> splitedClusters

  # registerDoParallel(cores = n_threads)
  # splitedClusters<-foreach(i=ind.withparalogues,
  #                          .inorder = T,
  #                          .options.multicore=list(preschedule=F))%dopar%{
  #                            splitClusters(clstr = fastas[pre.clusters[[i]]])
  #                          }

  # A patch for avoiding cds prediction errors in gff (strange 'proteins' that
  #  gives errors when are aligned with mafft).
  # which(sapply(splitedClusters,class)!='list') -> errors
  # if (length(errors)>0){
  #   pre.clusters[ind.withparalogues[errors]] -> werr
  #   names(werr) <- paste0('cluster_with_error(s)_',seq_len(length(werr)))
  #   writeClusters(outdir = outdir,
  #                 final.clusters = werr,
  #                 filename='err.txt')
  #   splitedClusters[-errors] -> splitedClusters
  #   pre.clusters[-ind.withparalogues[errors]] -> pre.clusters
  #   ind.withparalogues[-errors] -> ind.withparalogues
  #   fastas[-sapply(unlist(werr),function(x){which(names(fastas)==x)})] -> fastas
  # }

  cat(' DONE!\n')

  cat('Preparing output..\n')

  #Merge clusters (Clusters splited in previous step with clusters which
  # doesn't contained paralogues')
  cat('          ..merging clusters.. ')
  c(unlist(splitedClusters,recursive = F),
    pre.clusters[-ind.withparalogues]) -> final.clusters

  #No domain and no phmmer hit (names):
  names(fastas)[which(!names(fastas)%in%unlist(final.clusters))] -> hu


  #Merge clusters (Clusters merged in previous step with orphans sequences
  # which doesn't have any pfam domain or similarity with any other sequence[so
  # coudn't be captured by phmmer thus neither by mcl]).
  c(final.clusters,
    as.list(hu)
  ) -> final.clusters

  strsplit(names(final.clusters),';') -> pfamstr

  #Set cluster names
  names(final.clusters) <- setClusterNames(final.clusters = final.clusters)
  for (i in 1:length(final.clusters)){
    attr(final.clusters[[i]],'pfamStr') <- pfamstr[[i]]
  }

  #Clusters with recent paralogues:
  which(sapply(final.clusters,function(x){
    any(duplicated(do.call(rbind,strsplit(x,';'))[,1]))
    })) -> rcnt
  #Use the first one as representative and pass the other to attr(,'paralogues')
  for(i in rcnt){
    final.clusters[[i]][1] -> repr
    final.clusters[[i]][-1] -> toattr
    final.clusters[[i]] <- repr
    attr(final.clusters[[i]],'paralogues') <- toattr
  }

  #Pan-matrix (presence/absence)
  cat('          ..computing binary pan-matrix:')
  buildPanMatrix(pangenome = final.clusters,
                 type='binary') -> panm
  cat(' DONE!\n')
  #write.table(panm,file = paste0(outdir,'panmatrix.tab'),sep = '\t',quote = F)
  #cat(paste0(' DONE, saved at ',outdir,'panmatrix.tab\n'))

  #clusters.txt
  cat('          ..writing clusters:')
  writeClusters(outdir=outdir,
                final.clusters=final.clusters,
                filename='clusters.txt')
  cat(paste0(' DONE, saved at ',outdir,'clusters.txt\n'))

  #paralogues.txt
  cat('          ..writing paralogues:')
  writeParalogues(outdir=outdir,final.clusters=final.clusters)
  cat(paste0(' DONE, saved at ',outdir,'paralogues.txt\n'))

  #Out
  ncds <- length(fastas)
  norgs <- ncol(panm)
  nclust <- nrow(panm)

  list(final.clusters,
       panm) -> out
       # ffns,
       # pfamstr) -> out
  names(out) <- c('clusters','panmatrix')#,'pfam_structure')
  attr(out,'ncds') <- ncds
  attr(out,'norgs') <- norgs
  attr(out,'nclust') <- nclust
  attr(out,'output') <- outdir
  attr(out,'package') <- 'pewit'

  class(out) <- c('pangenome')

  saveRDS(out,file = paste0(outdir,'pangenome.rds'))

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

      buildPanMatrix(out,type=pmOutfileType) -> panm
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

  cat('\nFINNISH!\n\nThanks for using PEWIT.\n')
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
  titl <- 'Object of class "pangenome"\n\n'
  txt <- paste('Pangenome of',
               attr(x,'norgs'),
               'organisms with',
               attr(x,'ncds'),
               'genes.\n There are a total of',
               attr(x,'nclust'),
               'gene clusters identified.')
  writeLines(paste(titl,txt))
}


#' @name summary.pangenome
#' @title Summarizing pangenome object
#' @author Ignacio Ferres
#' @description Summary method for class \code{pangenome}
#' @param object An object of class \code{pangenome}
#' @param ... Unused
#' @return A \code{table} with: Number of organisms, Number of genes, Number of
#' clusters, Number of core-genes present at 100% of organisms, Number of soft-
#' core genes present in at least 95% of organisms, Number of accessory genes
#' present at less than 95% of organisms but that there are not singletones,
#' and Number of singletones.
#' @export
summary.pangenome <- function(object, ...){

  norgs <- ncol(object$panmatrix)
  ngenes <- attr(object,'ncds')
  nclust <- nrow(object$panmatrix)
  ncore100 <- length(which(rowSums(object$panmatrix)==ncol(object$panmatrix)))
  ncore95 <- length(which(rowSums(object$panmatrix)>=round(ncol(object$panmatrix)*0.95)))
  nsingles <- length(which(rowSums(object$panmatrix)==1))
  naccs <- length(which(rowSums(object$panmatrix)<round(ncol(object$panmatrix)*0.95)))-nsingles

  out <- c(norgs,
           ngenes,
           nclust,
           ncore100,
           ncore95,
           naccs,
           nsingles)

  names(out) <- c('Number of organisms',
                  'Number of genes',
                  'Number of clusters',
                  'Number of core-genes (100%)',
                  'Number of soft core-genes (>=95%)',
                  'Number of accessory genes(<95%)',
                  'Number of singletones')

  as.table(out, ...) -> out

  attr(out,'package') <- 'pewit'
  out
}


#' @name plot.pangenome
#' @title Plot a \code{pangenome} object.
#' @description Generic function to plot a \code{pangenome} object.
#' @param x \code{pangenome} object.
#' @param horiz \code{logical}. Should barplot be drawn horizontal?
#' @param border \code{logical}. Should borders be drawn?
#' @param las \code{numeric}. see \code{\link[graphics]{par}}.
#' @param xlim \code{vector}. limits for the x axis.
#' @param legend.text \code{logical}. see \code{\link[graphics]{barplot}}.
#' @param args.legend \code{list} of additional arguments to pass to
#'  \code{\link[graphics]{legend}}.
#' @param mar \code{vector}. see \code{\link[graphics]{par}}.
#' @param plot \code{logical} If FALSE, nothing is plotted.
#' @param ... arguments to be passed to \code{\link[graphics]{barplot}}.
#' @details A stacked barplot is generated. Each bar corresponds to each
#' organisms, and the bars show the proportion of soft-core, accessory and
#' singletons.
#' @return a \code{martrix} with the count of softcore, accessory and singleton
#' genes for each organism used to draw the barplot.
#' @importFrom graphics barplot par
#' @export
plot.pangenome <- function(x,
                           horiz=T,
                           border=T,
                           las=2,
                           xlim=c(0,pretty(max(colSums(x$panmatrix)))[2]),
                           legend.text=T,
                           args.legend=list(x='topright',bg=NULL,cex=0.5),
                           mar=c(4,7,3,2),
                           plot=T,
                           ...){

  socore <- x$panmatrix[which(rowSums(x$panmatrix)>=(ncol(x$panmatrix)*0.95)),]
  nscore <- colSums(socore)
  ssing <- x$panmatrix[which(rowSums(x$panmatrix)==1),]
  nsing <- colSums(ssing)
  accs <- x$panmatrix[which(rowSums(x$panmatrix)<(ncol(x$panmatrix)*0.95)),]
  naccs <- colSums(accs) - nsing

  o <- data.frame(row.names=colnames(x$panmatrix))
  o$Softcore <- nscore
  o$Accessory <- naccs
  o$Singletones <- nsing

  if (plot){
    par(mar=mar)
    barplot(t(as.matrix(o)),
            horiz = horiz,
            border = border,
            las=las,
            xlim = xlim,
            legend.text = legend.text,
            args.legend = args.legend,
            mar=mar,
            ...)
  }

  invisible(t(as.matrix(o)))
}



