#' @name extractSeqsFromGff3
#' @title Extract Sequences from Gff3 file
#' @description Extract Sequences from Gff3 file and outputs a list of element,
#' one per each CDS. Each element contain the DNA gene sequence and the AA
#' protein sequence.
#' @param infile \code{character}. The gff3 filename.
#' @return A \code{DNAStringSet}.
#' @author Gregorio Iraola and Ignacio Ferres.
#' @importFrom Biostrings DNAStringSet subseq reverseComplement
#' @importFrom S4Vectors mcols<-
extractSeqsFromGff3 <- function(infile) {

  rl <- readLines(infile)
  if (rl[1] != "##gff-version 3") {
    stop("One or more gff files seems not to be gff3 (version 3).")
  }

  if (grepl("#", infile)) {
    warning("Conflictive character \"#\", will be substituted by \"_\".")
    infile <- gsub("#", "_", infile)
  }

  # Save sequence in SeqFastadna object
  if (!any(grepl("##FASTA", rl))) {
    stop("One or more gff files do(es) not contain the fasta genome sequences.")
  }
  fasta.ini <- grep("##FASTA", rl) + 1
  fasta.end <- length(rl)
  fasta <- rl[fasta.ini:fasta.end]
  fheaders <- grep(">", fasta, fixed = T)
  t <- matrix(nrow = length(fheaders), ncol = 2)
  t[, 1] <- fheaders + 1
  t[, 2] <- c((fheaders - 1)[-1], length(fasta))
  ap <- apply(t, 1, function(x) {
    paste(fasta[x[1]:x[2]], collapse = "")
  })
  names(ap) <- gsub(">", "", fasta[fheaders])
  fnas <- DNAStringSet(ap)

  gfftable <- extractGffTable(rl = rl)

  cds <- gfftable[which(gfftable$Type == "CDS"), ]

  # Patch to avoid wrong formatted cds in gff3 file (predictions with artemis
  # return a different format than Prodigal or Aragorn, so they are discarded.
  # Usually there are very few (if any) proteins predicted with this method,
  # nothing to worry about).
  sqreg <- rl[which(grepl("^\\#\\#sequence-region", rl))]
  contigs <- do.call(rbind, strsplit(sqreg, " "))[, 2]
  bad <- which(!cds$Contig %in% contigs)
  if (length(bad) > 0) {
    cds <- cds[-bad, ]
  }

  fin <- subseq(fnas[cds$Contig], start = cds$From, end = cds$To)
  fin[which(cds$Strand=='-')] <- reverseComplement(fin[which(cds$Strand=='-')])
  nam <- cds$LocusTag

  if (any(grepl("#", nam))) {
    warning("Conflictive character \"#\", will be substituted by \"_\".")
    nam <- gsub("#", "_", nam)
  }
  names(fin) <- nam

  mcols(fin) <- DataFrame(geneName = cds$Gene,
                          organism = sub('[.]gff', '', basename(infile)),
                          product = cds$Product,
                          Contig = cds$Contig,
                          From = cds$From,
                          To = cds$To,
                          Strand = cds$Strand)


  return(fin)

}


#' @name extractGffTable
#' @title Extract the gff3 table and make it 'R'eadable.
#' @description Read a gff3 file and transforms the table in a \code{data.frame}.
#' @param rl \code{character}. A vector of character strings as passed by
#' \code{readLines()}, reading the gff3 file.
#' @return A \code{data.frame}.
#' @author Ignacio Ferres
extractGffTable <- function(rl) {
  w <- which(grepl("^\\#\\#", rl))
  upto <- rev(w)[1] - 1
  from <- rev(w)[2] + 1
  o <- rl[from:upto]

  lst <- strsplit(o, "\t")

  contig <- sapply(lst, function(x) {
    x[1]
  })
  type <- sapply(lst, function(x) {
    x[3]
  })
  from <- sapply(lst, function(x) {
    x[4]
  })
  to <- sapply(lst, function(x) {
    x[5]
  })
  strand <- sapply(lst, function(x) {
    x[7]
  })
  phase <- sapply(lst, function(x) {
    x[8]
  })
  attrib <- sapply(lst, function(x) {
    x[9]
  })

  metadata <- strsplit(attrib, ";")

  id <- sapply(metadata, function(x) {
    gp <- grep("ID=", x, value = T)
    if (length(gp) > 0) {
      sub("ID=", "", gp)
    } else {
      ""
    }
  })
  locustag <- sapply(metadata, function(x) {
    gp <- grep("locus_tag=", x, value = T)
    if (length(gp) > 0) {
      sub("locus_tag=", "", gp)
    } else {
      ""
    }
  })
  gene <- sapply(metadata, function(x) {
    gp <- grep("gene=", x, value = T)
    if (length(gp) > 0) {
      sub("gene=", "", gp)
    } else {
      ""
    }
  })
  product <- sapply(metadata, function(x) {
    gp <- grep("product=", x, value = T)
    if (length(gp) > 0) {
      sub("product=", "", gp)
    } else {
      ""
    }
  })

  out <- data.frame(Contig = contig,
                    ID = id,
                    LocusTag = locustag,
                    Gene = gene,
                    Product = product,
                    Type = type,
                    From = as.integer(from),
                    To = as.integer(to),
                    Strand = strand,
                    Phase = phase,
                    stringsAsFactors = F)

  out
}
