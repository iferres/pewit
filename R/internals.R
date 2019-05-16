# @name domainSearch
# @title Compute All Pfam Related Functions
# @description Computes all Pfam-A related functions.
# @param faas A \code{list} of amino acid sequences.
# @param hmmPfam \code{character} The path to the hmm file.
# @param datPfam \code{character} The path to the dat file.
# @param n_threads \code{integer} The number of cpus to use.
# @return a \code{list} with the domain and family clustering.
# author Ignacio Ferres
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
domainSearch <- function(faas,
                         hmm_pfam,
                         dat_pfam,
                         n_threads = 1L,
                         sep = '___',
                         verbose = TRUE) {

  # Process Pfam-A.dat
  if (verbose) message("   Processing Pfam-A.hmm.dat.")
  ref <- processPfam_A_Dat(datPfam = dat_pfam)

  if (verbose) message("   Preparing HMMSEARCH.")
  # Split indices and write fastas to distribute among threads with hmmscan
  temps <- splitAndWriteFastas(faas = faas, n_threads = n_threads)

  # Deprecated: #Run hmmscan (HMMER) cat('Running HMMSCAN against Pfam-A database
  # (this can take a while)..') registerDoParallel(cores = n_threads)
  # hmm.temps<-foreach(i=seq_along(temps),.inorder = F) %dopar% { tempfile(pattern
  # = 'tmpo',tmpdir = tempdir(),fileext = '.tab')->dmblout paste('hmmscan -o
  # /dev/null --domtblout ',dmblout, ' --noali --cut_ga --cpu 0 ', hmmPfam,'
  # ',temps[i],sep = '')->pfm system(pfm) dmblout } cat(' DONE!\n')

  # Index hmm if not yet
  if (any(!file.exists(paste0(hmm_pfam, c(".h3f", ".h3i", ".h3m", ".h3p"))))) {
    if (verbose) message("  Preparing Pfam-A.hmm files for hmmseach.")
    hmmpress <- paste0("hmmpress ", hmm_pfam)
    system(hmmpress)
  }

  # Run hmmsearch (HMMER)
  if (verbose) message("   Running HMMSEARCH against Pfam-A database (this can take a while).")
  hmm.temps <- mclapply(temps, function(x) {

    runHmmsearch(fasta = x, hmm = hmm_pfam, pfam = TRUE, n_threads = 0L)

  }, mc.cores = n_threads)
  # file.remove(temps)

  # Load hmmscan output and process
  if (verbose) message("   Processing hmmsearch output (resolving overlapping Pfam hits and building protein families profiles)")
  pout <- unlist(hmm.temps)

  tout <- mclapply(pout, processHmmsearch, ref = ref, mc.cores = n_threads)
  tout <- do.call(rbind, tout)
  file.remove(pout)

  dm <- which(tout$Domain!='')
  rw <- rownames(tout[dm, ])
  mcols(faas[rw])$arch <- tout$Domain[dm]

  fm <- which(tout$Domain=='' & tout$Family!='')
  rw <- rownames(tout[fm, ])
  mcols(faas[rw])$arch <- tout$Family[fm]

  if (verbose){
    lfaas <- length(faas)
    march <- mcols(faas)$arch
    lfawd <- length(which(!is.na(march)))
    alld <- length(unlist(strsplit(march, ';')))
    mssg1 <- paste0(alld, ' domains distributed among ', lfawd, ' proteins (out of ',lfaas,').')
    message(mssg1)
    diff <- length(unique(march))
    mssg2 <- paste(diff, 'distinct domain architectures found.')
    message(mssg2)
    mgene <- mcols(faas)$X
    mgene <- mgene[!is.na(march)]
    mgene[] <- lapply(mgene, function(x) sapply(strsplit(x, sep), '[', 1))
    tab <- table(unlist(mgene))
    nam <- names(tab)
    nch <- nchar(nam)
    upt <- nch[which.max(nch)] + 3L
    lna <- sapply(nch, function(x) paste0(rep('-', upt - x), collapse = ''))
    arr <- paste(lna, '>', sep = '')
    for (i in seq_along(tab)){
      mssg <- paste('', nam[i], arr[i], tab[i], 'proteins with domain architecture.')
      message(mssg)
    }
  }

  faas
}

# Distributes all proteins in many files as n_threads set in order to optimize
# the computing power in the following step (hmmscan).
# @name splitAndWriteFastas
# title Split sequences to distribute among threads
# description Split sequences to distribute among threads.
# param faas A \code{list} of aminoacid sequences in \code{SeqFastaAA}
# class.
# param n_threads \code{integer} The number of threads to use.
# details Distributes all proteins in many files as n_threads
# set in order to optimize the computing power in the following step(hmmscan).
# return A \code{vector} of temporary file names.
# author Ignacio Ferres
#' @importFrom parallel splitIndices
#' @importFrom Biostrings writeXStringSet
splitAndWriteFastas <- function(faas, n_threads) {
  spi <- splitIndices(length(faas), n_threads)
  names(spi) <- seq_along(spi)
  temps <- lapply(seq_along(spi), function(i){
    tpo <- tempfile(fileext = '_pewit.fasta')
    writeXStringSet(faas[spi[[i]]], filepath = tpo, format = 'fasta')
    tpo
  })
  temps <- unlist(temps)
  return(temps)
}


# @name processPfam_A_Dat
# title Process Pfam-A.dat file
# description Process Pfam-A.dat file.
# param datPfam \code{character}. The path to Pfam-A.dat file.
# return A \code{data.frame} with information of each Pfam entry.
# author Ignacio Ferres
#' @importFrom parallel mclapply
processPfam_A_Dat <- function(datPfam) {
  rl <- readLines(datPfam, skipNul = T)
  pr <- grep("STOCKHOLM 1.0", rl) + 1
  fn <- grep("//", rl) - 1
  pa <- cbind(pr, fn)
  li <- apply(pa, 1, function(x) {
    rl[x[1]:x[2]]
  })
  ref <- do.call(rbind, lapply(li, function(y) {
    g1 <- grep('^#=GF AC', y, value = TRUE)
    g2 <- grep('^#=GF TP', y, value = TRUE)
    g3 <- grep('^#=GF CL', y)
    c(g1, g2, ifelse(length(g3)!=0, y[g3], 'No-Clan'))
  }))
  colnames(ref) <- c("ID", "TP", "CL")
  ref <- sub('^#=\\w{2}[ ]\\w{2}[ ]{3}', '', ref)
  ref <- data.frame(ref, stringsAsFactors = F)
  return(ref)
}




