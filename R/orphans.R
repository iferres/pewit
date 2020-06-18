#' @importFrom utils write.table
#' @importFrom Biostrings writeXStringSet
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
clusterOrphans <- function(faas, n_threads, sep = '___', verbose = TRUE) {

  # Retrieve all orphan sequences from previous steps.
  orph.n <- names(faas)[is.na(mcols(faas)$arch)]

  orps <- orph.n

  count <- 0

  if (verbose){
    mssg <- paste('   Clustering', length(orps), 'sequences with no Pfam domains.')
    message(mssg)
  }
  while (length(orps) > 0) {

    tcont <- table(sapply(strsplit(orps, sep), '[', 1))
    so <- sort(tcont, decreasing = T)

    if (length(so) >= 3) {
      tk <- 3
    } else {
      tk <- length(so)
    }

    ma <- names(so[1:tk])
    fi <- unlist(sapply(ma, function(x) {
      grep(paste0(x, sep), orps, fixed = T)
    }, USE.NAMES = F))

    if (length(fi) >= n_threads) {
      parts <- n_threads
    } else {
      parts <- length(fi)
    }

    spin <- splitIndices(length(fi), parts)
    temps <- unlist(lapply(spin, function(x) {
      fixx <- orps[fi[x]]
      tmp1 <- tempfile(fileext = '_pewit.fasta')
      faa <- faas[fixx]
      writeXStringSet(faa, filepath = tmp1)
      tmp1
    }))

    tmp2 <- tempfile(fileext = '_pewit.fasta')
    writeXStringSet(faas[orps], filepath = tmp2)

    abc <- mclapply(seq_along(spin), function(x){
      phmm_tmp <- tempfile(fileext = '_pewit.tblout')
      system2(command = 'phmmer',
              args = c("-o /dev/null",
                       paste("--tblout", phmm_tmp),
                       "--cpu 0 --mx BLOSUM45",
                       temps[x], tmp2))
      file.remove(temps[x])
      file.remove(tmp2)
      outphmmer(pouti = phmm_tmp)
    }, mc.cores = n_threads)

    abc <- do.call(rbind, abc)

    # res[[length(res) + 1]] <- abc

    # MCL
    tmp <- tempfile(fileext = '_pewit.tab')
    write.table(abc[, c(1, 2, 4)],
                file = tmp,
                sep = "\t",
                quote = F,
                row.names = F,
                col.names = F)

    mclust <- strsplit(runMCL(abc = tmp, neg.log10 = F, infl = 4), split = "\t")
    # Append singletones without phmmer hits:
    mclust <- c(mclust, as.list(orps[fi][!orps[fi]%in%unlist(mclust)]))

    ncount <- count + length(mclust)
    names(mclust) <- paste('NOARCH', (count+1):ncount, sep = '_')
    count <- ncount

    mel <- melt(mclust)
    mcols(faas[mel$value])$arch <- mel$L1

    w <- which(orps %in% mel$value)
    orps <- orps[-w]

    if (verbose){
      mssg <- paste('   ', length(orps), 'orphans left.')
      message(mssg)
    }
  }

  faas
}



# name outphmmer
# title Process phmmer output to make ir 'R'eadable.
# description Process phmmer output to make it readable by R. Taken and
# adapted from \code{micropan} package (Lars Snipen and Kristian Hovde
# Liland).
# param pouti \code{character}. phmmer output temporary file.
# return A \code{data.frame} with the phmmer output.
# author Ignacio Ferres
outphmmer <- function(pouti) {
  rl <- readLines(pouti)
  file.remove(pouti)
  rl <- rl[which(!grepl("^\\#", rl))]
  rl <- gsub("[ ]+", " ", rl)
  lst <- strsplit(rl, " ")

  hit <- sapply(lst, function(x) {
    x[1]
  })
  query <- sapply(lst, function(x) {
    x[3]
  })
  eval <- as.numeric(sapply(lst, function(x) {
    x[5]
  }))
  score <- as.numeric(sapply(lst, function(x) {
    x[6]
  }))

  hmmer.table <- data.frame(Query = query, Hit = hit, Evalue = eval, Score = score,
                            stringsAsFactors = F)
  return(hmmer.table)
}


# name runMCL
# title Run MCL
# description Cluster protein sequences by running the Markov Clustering
# algorithm (MCL) over the phmmer comparison scores.
# param abc \code{character}. Temporary files of abc formatted phmmer
# comparisons.
# param neg.log10 \code{logical}. Should the negative of the log10 value be
# used as similarity measure?
# param infl \code{integer}. Inflacion value.
# return MCL output.
# author Ignacio Ferres
runMCL <- function(abc, neg.log10 = TRUE, infl = 6) {

  if (neg.log10) {
    arg <- paste("--abc --abc-neg-log10 -discard-loops y -abc-tf 'ceil(300)' -I",
                 infl)
  } else {
    arg <- paste("--abc -discard-loops y -I", infl)
  }

  # deprecated Run MCL system2(command = 'mcl', args = c(abc,arg, '-q x -o -'),
  # stdout = TRUE, stderr = FALSE)->rl rl

  # Run MCL
  tmpmcl <- tempfile()
  system(paste0("mcl ", abc, " ", arg, " -q x -o ", tmpmcl))
  rl <- readLines(tmpmcl)
  file.remove(tmpmcl)
  file.remove(abc)
  rl
}
