#' @name runHmmsearch
#' @title Run Hmmsearch (HMMER 3)
#' @description Takes a fasta file and a Hidden Markov Model profile and
#' performs a search of the former over the latter.
#' @param fasta A protein fasta file.
#' @param hmm A hmm file. Must be pressed (see hmmpress from HMMER manual).
#' @param pfam \code{logical}. If hmm file is the Pfam-A.hmm file or not
#' (custom hmm models).
#' @param n_threads An \code{integer}. The number of cores to use.
#' @return The path to a temporary file where the hmmsearch output is placed.
#' @author Ignacio Ferres
runHmmsearch <- function(fasta, hmm, pfam = FALSE, n_threads = 1L) {
  # run hmmsearch
  blout <- tempfile(pattern = "tmpo", fileext = ".tab")
  pfm <- paste0("hmmsearch -o /dev/null", ifelse(pfam, " --domtblout ", " --tblout "),
                blout, ifelse(pfam, " --noali --cut_ga --cpu ", " --noali -E 1e-10 --cpu "),
                as.character(n_threads), " ", hmm, " ", fasta)
  system(pfm)
  file.remove(fasta)
  blout

}

#' @name outhmmsearch
#' @title Process hmmsearch output to make it 'R'eadable.
#' @description Process hmmsearch output to make it readable by R.
#' @param pouti \code{character}. hmmsearch output temporary file.
#' @param ref \code{data.frame}. Information about each Pfam-A entry.
#' @return A \code{data.frame} with the hmmsearch output plus information about
#' Pfam-A hits.
#' @note Taken and adapted from \code{micropan} package (Lars Snipen and
#' Kristian Hovde Liland).
#' @author Ignacio Ferres
outhmmsearch <- function(pouti, ref) {
  rl <- readLines(pouti)
  rl <- rl[which(!grepl("^\\#", rl))]
  rl <- gsub("[ ]+", " ", rl)
  lst <- strsplit(rl, " ")

  query <- sapply(lst, function(x) {
    x[1]
  })
  hit <- sapply(lst, function(x) {
    x[4]
  })
  pfmID <- sapply(lst, function(x) {
    x[5]
  })
  eval <- as.numeric(sapply(lst, function(x) {
    x[13]
  }))
  score <- as.numeric(sapply(lst, function(x) {
    x[14]
  }))
  st <- as.numeric(sapply(lst, function(x) {
    x[18]
  }))
  en <- as.numeric(sapply(lst, function(x) {
    x[19]
  }))
  desc <- sapply(lst, function(x) {
    paste(x[23:length(x)], collapse = " ")
  })

  hmmer.table <- data.frame(Query = query, Hit = hit, PfamID = pfmID, Evalue = eval,
                            Score = score, Start = st, End = en, Description = desc, stringsAsFactors = F)

  wh <- which(hmmer.table$Evalue > 0.1)
  if (length(wh) > 0) {
    hmmer.table <- hmmer.table[-which(hmmer.table$Evalue > 0.1), ]
  }

  codpf <- sub("\\.\\d+", "", ref$ID)
  codhit <- sub("\\.\\d+", "", hmmer.table$PfamID)
  ind <- lapply(codhit, function(x) {
    which(codpf == x)
  })
  ind <- unlist(ind)
  hmmer.table$Type <- ref$TP[ind]
  hmmer.table$Clan <- ref$CL[ind]
  return(hmmer.table)
}


#' @name processHmmsearch
#' @title Main hmmsearch output processing function
#' @description Outputs a data.frame with the domain structure of each sequence
#' (rows). It has 6 columns, each one with the semi-colon separated found
#' sequences of Pfam (PF) entries of the following 6 categories: 'Domain',
#' 'Family', 'Motif', 'Repeat','Disordered' or 'Coiled-coil'.
#' @param pout \code{character}. Temporary file name (hmmsearch output).
#' @param ref \code{data.frame}. Information about each Pfam-A entry.
#' @return A \code{data.frame} with the domain structure of each sequence.
#' @author Ignacio Ferres
processHmmsearch <- function(pout, ref) {
  # writeLines('Loading hmmsearch output..') outhmmscan(pouti = pout,ref = ref)->t
  t <- outhmmsearch(pout, ref)


  #### Keep single hit ####

  # Single hit (Domain)
  dom.sing <- grouping(t = t, singlehit = TRUE, type = "Domain")
  # Single hit (Family)
  fam.sing <- grouping(t = t, singlehit = TRUE, type = "Family")
  # Single hit (Motif)
  mot.sing <- grouping(t = t, singlehit = TRUE, type = "Motif")
  # Single hit (Repeat)
  rep.sing <- grouping(t = t, singlehit = TRUE, type = "Repeat")
  # Single hit (Disordered)
  dis.sing <- grouping(t = t, singlehit = TRUE, type = "Disordered")
  # Single hit (Coiled-coil)
  cco.sing <- grouping(t = t, singlehit = TRUE, type = "Coiled-coil")



  #### Group by query those which have more than one hit ####

  gping <- grouping(t = t, singlehit = FALSE)


  #### Remove same-clan overlapping hits #### cat('Removing same-clan overlapping
  #### domains..')

  # Domain
  dom.novlp <- lapply(gping, function(x) {
    rm_overlaping_clans(tping = x, type = "Domain")
  })
  # Family
  fam.novlp <- lapply(gping, function(x) {
    rm_overlaping_clans(tping = x, type = "Family")
  })
  # Motif
  mot.novlp <- lapply(gping, function(x) {
    rm_overlaping_clans(tping = x, type = "Motif")
  })
  # Repeat
  rep.novlp <- lapply(gping, function(x) {
    rm_overlaping_clans(tping = x, type = "Repeat")
  })
  # Disordered
  dis.novlp <- lapply(gping, function(x) {
    rm_overlaping_clans(tping = x, type = "Disordered")
  })
  # Coiled-coil
  cco.novlp <- lapply(gping, function(x) {
    rm_overlaping_clans(tping = x, type = "Coiled-coil")
  })
  # cat(' DONE!')

  #### Merge ####

  # List of lists
  lpm <- list(dom.sing, fam.sing, mot.sing, rep.sing, dis.sing, cco.sing, dom.novlp,
              fam.novlp, mot.novlp, rep.novlp, dis.novlp, cco.novlp)
  # Index of non-empthy lists
  ind.lpm <- which(unlist(lapply(lpm, function(x) {
    !is.null(unlist(x))
  })))

  # Filtering done!
  t <- do.call(rbind, unlist(lpm[ind.lpm], recursive = F))
  t <- t[order(t$Query), ]
  rownames(t) <- seq_len(nrow(t))


  #### Build table (data.frame) with domain structure of each secuence ####

  tout <- as.data.frame(matrix(nrow = length(unique(t$Query)), ncol = 6), row.names = unique(as.character(t$Query)))
  colnames(tout) <- c("Domain", "Family", "Motif", "Repeat", "Disordered", "Coiled-coil")

  tdom <- t[which(t$Type == "Domain"), ]
  tfam <- t[which(t$Type == "Family"), ]
  tmot <- t[which(t$Type == "Motif"), ]
  trep <- t[which(t$Type == "Repeat"), ]
  tdis <- t[which(t$Type == "Disordered"), ]
  tcoi <- t[which(t$Type == "Coiled-coil"), ]

  # Dom
  tout$Domain <- do.call(rbind, lapply(unique(t$Query), function(x) {
    paste(tdom$Hit[grep(x, tdom$Query, fixed = T)], collapse = ";")
  }))[, 1]
  # Fam
  tout$Family <- do.call(rbind, lapply(unique(t$Query), function(x) {
    paste(tfam$Hit[grep(x, tfam$Query, fixed = T)], collapse = ";")
  }))[, 1]
  # Mot
  tout$Motif <- do.call(rbind, lapply(unique(t$Query), function(x) {
    paste(tmot$Hit[grep(x, tmot$Query, fixed = T)], collapse = ";")
  }))[, 1]
  # Rep
  tout$Repeat <- do.call(rbind, lapply(unique(t$Query), function(x) {
    paste(trep$Hit[grep(x, trep$Query, fixed = T)], collapse = ";")
  }))[, 1]
  # Dis
  tout$Disordered <- do.call(rbind, lapply(unique(t$Query), function(x) {
    paste(tdis$Hit[grep(x, tdis$Query, fixed = T)], collapse = ";")
  }))[, 1]
  # Coi
  tout$`Coiled-coil` <- do.call(rbind, lapply(unique(t$Query), function(x) {
    paste(tcoi$Hit[grep(x, tcoi$Query, fixed = T)], collapse = ";")
  }))[, 1]



  return(tout)
}


#' @name grouping
#' @title Split hmmscan output by query.
#' @description Split hmmscan output by query.
#' @param t \code{data.frame}. hmmscan output passed by \link{outhmmsearch}.
#' @param singlehit \code{logical}. Has this query just one hit?
#' @param type \code{character}. One of 'Domain', 'Family', 'Motif', 'Repeat',
#' 'Disordered' or 'Coiled-coil'.
#' @return A \code{list} of \code{data.frame} with hits by query.
#' @author Ignacio Ferres
grouping <- function(t, singlehit = FALSE, type = "Domain") {
  if (!singlehit) {
    x <- split.data.frame(t, t$Query)
    y <- unlist(lapply(x, function(x) {
      nrow(x) > 1
    }))
    return(x[which(y)])
  } else {
    x <- split.data.frame(t, t$Query)
    y <- which(unlist(lapply(x, function(z) {
      nrow(z) == 1
    })))
    j <- which(unlist(lapply(x[y], function(z) {
      z$Type == type
    })))
    return(x[y[j]])
  }
}


#' @name rm_overlaping_clans
#' @title Removes same-clan overlapping domains
#' @description Removes those overlapping domains (which belong to the same
#' clan) with higher e-value. In a case like (a[b)c] with
#' \code{which.min(evalue)=='a'}, the algorithm will keep only domain 'a', even
#' if 'a' doesn't overlaps with 'c'. Will also remove overlapping domains with
#' no clan assigned, following the same above criteria as if all no clan
#' assigned domains were from the same clan.
#' @param tping A \code{data.frame} with all hits for each sequence.
#' @param type \code{character}. One of 'Domain', 'Family', 'Motif', 'Repeat',
#' 'Disordered' or 'Coiled-coil'.
#' @return a \code{data.frame} of non-overlapping domains.
#' @author Ignacio Ferres
rm_overlaping_clans <- function(tping, type) {

  # Just consider one type tpingo<-tping[which(tping$Type==type),,drop=F]
  tpingo <- as.data.frame(tping)[which(tping["Type"] == type), ]

  if (dim(tpingo)[1] != 0) {
    # Subset by clan
    lping <- split.data.frame(tpingo, as.factor(tpingo$Clan))

    for (i in 1:length(lping)) {

      if (dim(lping[[i]])[1] > 1) {

        # Order by start (this is in order to allow determineOverlap() to work)
        lping[[i]] <- lping[[i]][order(lping[[i]]$Start), ]

        # Create list with sequence positions
        inc <- unlist(apply(lping[[i]], 1, function(p) {
          list(seq(as.numeric(p[6]), as.numeric(p[7]), 1))
        }), recursive = F)

        # Logical matrix of (non-)overlaping PF -determines shared positions-.
        m <- matrix(unlist(lapply(inc, function(y) {
          lapply(inc, function(p) {
            any(p %in% y)
          })
        })), nrow = length(inc), dimnames = list(rownames(lping[[i]]), rownames(lping[[i]])))

        # Output list of overlapping domains (row names)
        ovlp <- determineOverlap(m = m)


        # Process those with same-clan overlapping domains
        multiclan <- which(sapply(ovlp, length) > 1)
        if (length(multiclan) >= 1) {
          dej <- list()
          for (k in seq_along(multiclan)) {
            #columna de evalue
            indx <- which.min(lping[[i]][ovlp[[multiclan[k]]],
                                    4])
            # Remove the rest
            dej[[k]] <- lping[[i]][ovlp[[multiclan[k]]][indx], , drop = F]
          }
          dej <- do.call(rbind, dej)
          lping[[i]] <- do.call(rbind, list(lping[[i]][unlist(ovlp[-multiclan]),
                                                       ], dej))
        }


      }

    }
    lping <- do.call(rbind, lping)

  } else {
    lping <- NULL
  }
  return(lping)
}


#' @name determineOverlap
#' @title Determine which Pfam domains overlaps.
#' @description Takes an overlappling-domain logical matrix \code{m} and
#' through a 'trace-forward' precedure from \code{m[1,1]} to \code{m[ ,ncol(m)]}
#' outputs a list with overlappling-domains row names.
#' @param m A logical \code{matrix} indicating if 'x' domain overlaps with 'y'.
#' @return A \code{list} of overlapping domains.
#' @author Ignacio Ferres
determineOverlap <- function(m) {
  st <- c(1, 1)
  li <- list()
  # start from m[1,1]
  li[[1]] <- colnames(m)[st[2]]
  # while 'st' doesn't reach the right edge of the matrix...
  while (st[2] != ncol(m)) {
    # If next (right) cell is TRUE..
    if (m[st[1], st[2] + 1] == T) {
      # ..then move to next cell.
      st[2] <- st[2] + 1
      li[[length(li)]][length(li[[length(li)]]) + 1] <- colnames(m)[st[2]]
      # If next cell is FALSE but bottom cell is TRUE...
    } else if (m[st[1] + 1, st[2]] == T) {
      # ..then move to bottom cell.
      st[1] <- st[1] + 1
      # If both next and bottom cells are FALSE...
    } else {
      # ..then move to the bottom-right (diagonal) cell.
      st[1] <- st[1] + 1
      st[2] <- st[2] + 1
      li[length(li) + 1] <- colnames(m)[st[2]]
    }
  }
  return(li)
}
