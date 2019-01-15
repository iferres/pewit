# Internal functions.

#' @name runOnExit
#' @title Removes everything if an error occur
#' @description Removes the output folder if an error occurs after finnishing.
#' @param outdir The normalized path of the output directory
#' @return Warning message if exits before finnishing all the process.
#' @author Ignacio Ferres
runOnExit <- function(outdir) {
  fi <- paste0(outdir, "pangenome.rds")
  if (!file.exists(fi)) {
    warning("Something gone wrong: removing output directory.")
    unlink(outdir, recursive = TRUE)
  }
}


#' @name domainSearch
#' @title Compute All Pfam Related Functions
#' @description Computes all Pfam-A related functions.
#' @param fastas A \code{list} of amino acid sequences.
#' @param hmmPfam \code{character} The path to the hmm file.
#' @param datPfam \code{character} The path to the dat file.
#' @param n_threads \code{integer} The number of cpus to use.
#' @return a \code{list} with the domain and family clustering.
#' @importFrom parallel mclapply
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach '%dopar%'
#' @author Ignacio Ferres
domainSearch <- function(fastas, hmmPfam, datPfam, n_threads = 1L) {

  # Process Pfam-A.dat
  cat("Processing Pfam-A.hmm.dat..")
  ref <- processPfam_A_Dat(datPfam = datPfam, n_threads = n_threads)
  cat(" DONE!\n")

  cat("Preparing HMMSEARCH.\n")
  # Split indices and write fastas to distribute among threads with hmmscan
  temps <- splitAndWriteFastas(fastas = fastas, n_threads = n_threads)

  # Deprecated: #Run hmmscan (HMMER) cat('Running HMMSCAN against Pfam-A database
  # (this can take a while)..') registerDoParallel(cores = n_threads)
  # hmm.temps<-foreach(i=seq_along(temps),.inorder = F) %dopar% { tempfile(pattern
  # = 'tmpo',tmpdir = tempdir(),fileext = '.tab')->dmblout paste('hmmscan -o
  # /dev/null --domtblout ',dmblout, ' --noali --cut_ga --cpu 0 ', hmmPfam,'
  # ',temps[i],sep = '')->pfm system(pfm) dmblout } cat(' DONE!\n')

  # Index hmm if not yet
  if (any(!file.exists(paste0(hmmPfam, c(".h3f", ".h3i", ".h3m", ".h3p"))))) {
    cat("Preparing Pfam-A.hmm files for hmmscan search.\n")
    hmmpress <- paste0("hmmpress ", hmmPfam)
    system(hmmpress)
  }

  # Run hmmsearch (HMMER)
  cat("Running HMMSEARCH against Pfam-A database (this can take a while)..")
  hmm.temps <- mclapply(temps, function(x) {

    runHmmsearch(fasta = x, hmm = hmmPfam, pfam = TRUE, n_threads = 0L)

  }, mc.cores = n_threads)
  cat(" DONE!\n")
  # file.remove(temps)

  # Load hmmscan output and process
  cat("Processing hmmsearch output (resolving overlapping Pfam hits and building protein families profiles)..")
  pout <- unlist(hmm.temps)

  registerDoParallel(cores = n_threads)
  tout <- foreach(i = seq_along(pout), .combine = rbind, .inorder = T) %dopar%
  {
    processHmmsearch(pout = pout[i], ref = ref)
  }
  cat(" DONE!\n")
  file.remove(pout)

  # Clustering by domain structure only
  cat("Clustering sequences per domain structure..")
  clu_dom <- split(x = rownames(tout), f = tout$Domain)[-1]
  cat(" DONE!\n")

  # Cluster sequences without domains asigned, by family
  cat("Clustering sequences per family..")
  clu_fam <- split(x = rownames(tout[which(tout$Domain == ""), ]), f = tout$Family[which(tout$Domain ==
                                                                                           "")])[-1]
  cat(" DONE!\n")

  o <- list(clu_dom, clu_fam, tout)
  return(o)
}

# Distributes all proteins in many files as n_threads set in order to optimize
# the computing power in the following step (hmmscan).
#' @name splitAndWriteFastas
#' @title Split sequences to distribute among threads
#' @description Split sequences to distribute among threads.
#' @param fastas A \code{list} of aminoacid sequences in \code{SeqFastaAA}
#' class.
#' @param n_threads \code{integer} The number of threads to use.
#' @details Distributes all proteins in many files as n_threads
#' set in order to optimize the computing power in the following step(hmmscan).
#' @return A \code{vector} of temporary file names.
#' @author Ignacio Ferres
#' @importFrom parallel splitIndices
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach '%dopar%'
#' @importFrom seqinr write.fasta
splitAndWriteFastas <- function(fastas, n_threads) {

  spi <- splitIndices(length(fastas), n_threads)
  names(spi) <- seq_along(spi)
  tdir <- tempdir()
  registerDoParallel(cores = n_threads)
  temps <- foreach(i = seq_along(spi)) %dopar% {
    tpo <- tempfile(pattern = "tmpo", tmpdir = tdir, fileext = ".fasta")
    # lapply(fastas[spi[[i]]],function(x){ memDecompress(x,type = 'gzip',asChar = T)
    # }) -> sq names(sq) -> nsq write.fasta(sq, names = nsq, file.out = tpo)
    write.fasta(fastas[spi[[i]]], names = names(fastas[spi[[i]]]), file.out = tpo)
    tpo
  }

  temps <- unlist(temps)
  return(temps)
}


#' @name processPfam_A_Dat
#' @title Process Pfam-A.dat file
#' @description Process Pfam-A.dat file.
#' @param datPfam \code{character}. The path to Pfam-A.dat file.
#' @param n_threads \code{integer}. The number of threads to use.
#' @return A \code{data.frame} with information of each Pfam entry.
#' @author Ignacio Ferres
#' @importFrom parallel mclapply
processPfam_A_Dat <- function(datPfam, n_threads) {
  rl <- readLines(datPfam, skipNul = T)
  pr <- grep("STOCKHOLM 1.0", rl) + 1
  fn <- grep("//", rl) - 1
  pa <- cbind(pr, fn)
  li <- apply(pa, 1, function(x) {
    rl[x[1]:x[2]]
  })
  ref <- do.call(rbind, mclapply(li, function(y) {
    g1 <- grep('^#=GF AC', y, value = TRUE)
    g2 <- grep('^#=GF TP', y, value = TRUE)
    g3 <- grep('^#=GF CL', y)
    c(g1, g2, ifelse(length(g3)!=0, y[g3], 'No-Clan'))
  }, mc.cores = n_threads))
  colnames(ref) <- c("ID", "TP", "CL")
  ref <- sub('^#=\\w{2}[ ]\\w{2}[ ]{3}', '', ref)
  ref <- data.frame(ref, stringsAsFactors = F)
  return(ref)
}

#' # Deprecated
#' #@name outhmmscan
#' #@title Process hmmscan output to make it 'R'eadable.
#' #@description Process hmmscan output to make it readable by R.
#' #@param pouti \code{character}. hmmscan output temporary file.
#' #@param ref \code{data.frame}. Information about each Pfam-A entry.
#' #@return A \code{data.frame} with the hmmscan output plus information about
#' Pfam-A hits.
#' #@note Taken and adapted from \code{micropan} package (Lars Snipen and
#' Kristian Hovde Liland).
#' #@author Ignacio Ferres
# outhmmscan<-function(pouti,ref){ readLines(pouti)->rl
# rl[which(!grepl('^\\#',rl))]->rl gsub('[ ]+',' ',rl)->rl strsplit(rl,'
# ')->lst hit<-sapply(lst,function(x){x[1]}) pfmID<-sapply(lst,function(x){x[2]})
# query<-sapply(lst,function(x){x[4]})
# eval<-as.numeric(sapply(lst,function(x){x[13]}))
# score<-as.numeric(sapply(lst,function(x){x[14]})) st<-as.numeric(sapply(lst,
# function(x){x[18]})) en<-as.numeric(sapply(lst,function(x){x[19]}))
# desc<-sapply(lst,function(x){paste(x[23:length(x)],collapse = ' ')})
# hmmer.table<-data.frame(Query=query,Hit=hit,PfamID=pfmID,Evalue=eval,Score=score,
# Start=st,End=en,Description=desc,stringsAsFactors = F)
# hmmer.table<-hmmer.table[-which(hmmer.table$Evalue>0.1),]
# sub('\\.\\d+','',ref$ID) -> codpf sub('\\.\\d+','',hmmer.table$PfamID)
# -> codhit sapply(codhit,function(x){which(codpf==x)})->ind
# ref$TP[ind]->hmmer.table$Type ref$CL[ind]->hmmer.table$Clan return(hmmer.table)
# }

#' @name checkIfParalogues
#' @title Check if cluster contains paralogues
#' @description Check if a cluster contains paralogues.
#' @param p A vector of proteins
#' @return \code{logical}, if contains paralogues (i.e. proteins of the same
#' organism).
#' @author Ignacio Ferres
checkIfParalogues <- function(p) {
  st <- sapply(p, function(y) {
    strsplit(y, ";")[[1]][1]
  })
  if (any(duplicated(st))) {
    TRUE
    if (length(unique(st)) == 1) {
      FALSE
    } else {
      TRUE
    }
  } else {
    FALSE
  }
}


#' @name setClusterNames
#' @title Set cluster names
#' @description Creates othologue cluster names.
#' @param clusters A \code{list} of protein clusters.
#' @return A \code{vector} of names for the clusters
#' @author Ignacio Ferres
setClusterNames <- function(clusters) {
  np <- nchar(as.character(length(clusters)))
  np <- paste("%0", np, "d", sep = "")
  npnam <- paste("OG", sprintf(np, 1:length(clusters)), sep = "")
  npnam
}


#' @name writeClusters
#' @title Write Clusters
#' @description Write clusters in a fasta-like format.
#' @param outdir Output directory
#' @param final.clusters A list of clusters
#' @param filename File name.
#' @return A file 'clusters.txt' in \code{outdir}.
#' @author Ignacio Ferres
writeClusters <- function(outdir, final.clusters, filename) {
  sink(paste0(outdir, filename))
  for (i in 1:length(final.clusters)) {
    pfm <- attr(final.clusters[[i]], "pfamStr")
    if (length(pfm) > 0) {
      pfm <- paste("[", paste(pfm, collapse = ";"), "]", sep = "")
    } else {
      pfm <- ""
    }
    cat(paste0(">", names(final.clusters[i]), " ", pfm, "\n"))
    cat(final.clusters[[i]], sep = " ", fill = F)
    cat("\n")
  }
  sink()
}


#' @name writeParalogues
#' @title Write Paralogues
#' @description Write clusters which contain paralogues in a fasta-like format.
#' Only clusters which do contain paralogues are written; the genes that
#' were already written in 'clusters.txt' (representatives) are not re-written.
#' @param outdir Output directory
#' @param final.clusters A list of clusters
#' @return A file 'paralogues.txt' in \code{outdir}.
#' @author Ignacio Ferres
writeParalogues <- function(outdir, final.clusters) {
  sink(paste0(outdir, "paralogues.txt"))
  for (i in 1:length(final.clusters)) {
    pfm <- attr(final.clusters[[i]], "pfamStr")
    if (length(pfm) > 0) {
      pfm <- paste("[", paste(pfm, collapse = ";"), "]", sep = "")
      if (!is.null(attributes(final.clusters[[i]])$paralogues)) {
        cat(paste0(">", names(final.clusters[i]), " ", pfm, "\n"))
        cat(attr(final.clusters[[i]], "paralogues"), sep = " ", fill = F)
        cat("\n")
      }
    } else {
      pfm <- ""
      if (!is.null(attributes(final.clusters[[i]])$paralogues)) {
        cat(paste0(">", names(final.clusters[i]), " ", pfm, "\n"))
        cat(attr(final.clusters[[i]], "paralogues"), sep = " ", fill = F)
        cat("\n")
      }
    }

  }
  sink()
}
