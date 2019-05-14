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
#' @param gffs A \code{vector} of gff3 file paths as retrieved by prokka (Seeman,
#' 2014). These files must have a ".gff" extension.
#' @param hmm_pfam \code{character} with the path to Pfam-A.hmm file. If any of
#' \code{hmmPfam} or \code{datPfam} is left \code{NULL}, then the pfam step is
#' skipped and sequences will be preclustered by \code{phmmer} and \code{mcl}.
#' @param dat_pfam \code{character} with the path to Pfam-A.hmm.dat file. If any
#'  of \code{hmm_pfam} or \code{dat_pfam} is left \code{NULL}, then the pfam step
#'  is skipped and sequences will be preclustered by \code{phmmer} and \code{mcl}.
#' @param n_threads \code{integer}. The number of threads to use.
#' @param group_prefix \code{character}. The prefix you want to use to name the
#' cluster of orthologues. Default is "group". Cluster names will be formed by
#' \code{paste}ing this prefix to a number identifiying each cluster, in the
#' form [group_prefix]\code{0001}.
#' @param sep \code{character}. A separator to form unique gene names using the
#' name of the organism, in the form (if \code{sep='___'}) \code{organism___geneid}. This is
#' mostly used inside the function, although is used to form unique sequence
#' names in the final output. Default is "___" (3 underscores), and should be
#' changed if any of your input files or gene identifiers already contain this
#' separator string, so it not conflicts the process.
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
#' followed by all-vs-all distance calculation is performed. To align DNA
#' sequences, function \link[DECIPHER]{AlignTranslation} (DECIPHER package) is
#' used. A neighbour-joining tree is then generated and midpoint rooted.
#' If recent paralogues are detected (this is, nodes which all its descendants
#' belong to the same organism), then all but one tips are removed from the
#' gene-tree and saved.
#' Finally, the gene-tree is splited in as many sub-trees as necessary to have
#' the minimum set of sub-trees with just one gene of each organism. Paralogues
#' are then assigned to the cluster with the copheneticaly closer tip.
#'
#' In that way true-orthologue clusters are retrieved and saved in an object of
#' class "PgR6MS".
#' @return An object of class \link[pagoo]{PgR6MS} (pagoo package). For details
#' of this type of object, consult \code{help('PgR6MS')}, or go to
#' \url{https://github.com/iferres/pagoo} for a tutorial showing how to use this
#' class. A bunch of statistical and visualization methods are provided, as well
#' as examples of pangenome manipulation, and how to use it with popular third
#' party R packages to perform comparative genomics analyses, infer phylogenies,
#' or create highly customized figures.
#' @note External dependencies are HMMER3 and MCL, without them, running the
#' function will stop with an error.
#' @importFrom parallel mclapply
#' @importFrom S4Vectors mcols mcols<- DataFrame List
#' @importFrom Biostrings DNAStringSetList translate
#' @importFrom reshape2 melt
#' @import pagoo
#' @export
pangenome <- function(gffs,
                      hmm_pfam,
                      dat_pfam,
                      n_threads = 1,
                      group_prefix = 'group',
                      sep = '___',
                      verbose = TRUE){


  if (Sys.which("hmmsearch")==""){
    stop("\n\tHMMER (v.3) is not installed. (Couldn't find 'hmmsearch' in $PATH)
         \tPlease install it before re-running pangenome().\n\n")
  }

  if (Sys.which("mcl")==""){
    stop("\n\tMCL is not installed. (Couldn't find 'mcl' in $PATH)
         \tPlease install it before re-running pangenome().\n\n")
  }

  if(any(grepl('\\.gff$',gffs)==F)){
    stop("Files ('gffs') must have '.gff' extension.")
  }


  ## For debug:
  # tgz <- system.file('extdata', 'Hinfluenzae_example.tar.gz', package = 'pewit2')
  # untar(tarfile = tgz, exdir = tempdir())
  # gffs <- list.files(path = tempdir(), pattern = '[.]gff$', full.names = TRUE)

  if (verbose) message('Extracting sequences from gff files.')
  fastas <- mclapply(gffs,function(x){
    extractSeqsFromGff3(infile = x)
  }, mc.cores = n_threads, mc.preschedule = FALSE)

  mcls <- lapply(fastas, mcols)
  mcls <- do.call(rbind, mcls)
  fastas <- unlist(DNAStringSetList(fastas))
  mcols(fastas) <- mcls
  names(fastas) <- paste(mcls$organism, names(fastas), sep = sep)
  faas <- translate(fastas, if.fuzzy.codon = 'solve')
  aa_factor <- as.integer(factor(as.character(faas)))
  mcols(fastas)$aa_factor <- aa_factor
  faas <- unique(faas)
  tap <- tapply(names(fastas), aa_factor, c, simplify = FALSE)[unique(aa_factor)]
  attr(tap, 'dim') <- NULL
  mcols(faas) <- List(tap)


  lfs <- length(fastas)
  lfa <- length(faas)
  if(lfs!=lfa & verbose){
    redundant <- lfs - lfa
    mssg <- paste('   ', redundant, 'out of', lfs, 'sequences are redundant at amino acid level.')
    message(mssg)
  }

  if (!any(sapply(list(hmm_pfam, dat_pfam), is.null))){

    if (verbose) message('Searching Pfam domains, resolving overlap, and clustering.')
    faas <- domainSearch(faas = faas,
                         hmm_pfam = hmm_pfam,
                         dat_pfam = dat_pfam,
                         n_threads = n_threads)

  }else{

    warning('Either hmmPfam or/and datPfam arguments were left NULL.
         Pfam search will be skipped.\n', immediate. = TRUE)
    mcols(faas)$arch <- NA

  }


  if (verbose) message('Clustering orphan sequences (phmmer + mcl).')
  faas <- clusterOrphans(faas = faas,
                         n_threads = n_threads,
                         sep = sep,
                         verbose = verbose)


  if (verbose) message('Merging.')
  names(mcols(faas)$X) <- mcols(faas)$arch
  df <- DataFrame(melt(mcols(faas)$X)[, 2:3])
  colnames(df) <- c('Arch', 'Gid')
  mcols(fastas[df$Gid])$Arch <- df$Arch

  if (verbose) message('Resolving coarse clusters.')
  fastas <- splitPreClusters(fastas, n_threads, sep, verbose)
  mcls <- mcols(fastas)

  if (verbose) message('Preparing output.')
  DF <- list()
  # Gene name
  DF$gene <- sapply(strsplit(names(fastas), sep), '[', 2)
  # Group name
  groups <- mcls$Cluster
  wdth <- nchar(max(groups))
  DF$group <- paste0(group_prefix,
                     formatC(groups,
                             width = wdth,
                             format = 'd',
                             flag = '0'))
  # Organism
  DF$org <- mcls$organism
  # Gene name (as annotated)
  DF$geneName <- mcls$geneName
  # Domain architecture
  DF$Pfam_Arch <- mcls$Arch
  # Annotation
  DF$product <- mcls$product
  # Contig
  DF$contig <- mcls$Contig
  # Possition from
  DF$from <- mcls$From
  # Possition to
  DF$to <- mcls$To
  # Strand
  DF$strand <- mcls$Strand
  # DF
  DF <- DataFrame(DF)
  # Sequences
  seqs <- split(fastas, mcols(fastas)$organism)
  seqs <- DNAStringSetList(lapply(seqs, function(x){
    names(x)<- sapply(strsplit(names(x), sep), '[', 2)
    x
  }))

  pagoo_object <- PgR6MS$new(DF = DF,
                             sep = sep,
                             sequences = seqs)

  if (verbose) message('FINISH!
Returning an object of class "PgR6MS" (pagoo package, R6 class system)')
  return(pagoo_object)

}