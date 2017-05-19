#Heuristic phmmer comparisons in order to cluster proteins without Pfam-A
# domains asigned.
#' @name clusterOrphans
#' @title Cluster 'orphan' sequences (i.e. without Pfam-A domains of classes
#' 'Domain' or 'Family' found).
#' @description In order to avoid exponencial number of comparisons with the
#' number of sequences, an heuritic is implemented. It first takes all the
#' proteins of the 3 organisms which have more orphans and performs phmmer
#' comparisons against all orphans. The phmmer output is saved. Then the Markov
#' Clustering algorithm (MCL) clusters the sequences and those groups which
#' contain proteins of the first 3 organisms are provisionally removed from the
#' analysis. Then another round of phmmer comparisons is performed with the
#' proteins of the 3 organisms with more orphans left against all orphans left,
#' saves the phmmer output, and MCL clusters the sequences.
#'
#' The same is repeated until there are no more sequences left. Then it merge
#' all the saved phmmer output tables and performs the final MCL clustering
#' over it.
#' @param tout A \code{data.frame} with non overlapping domain hmmscan output
#' as passed by \link{processHmmsearch}.
#' @param fastas A \code{list} of protein sequeces as \code{SeqFastaAA} objects
#' @param n_threads \code{integer} Number of threads to use.
#' @return A \code{list} of protein clusters.
#' @author Ignacio Ferres
#' @importFrom parallel splitIndices
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach '%dopar%'
#' @importFrom seqinr write.fasta
#' @importFrom utils write.table
clusterOrphans <- function(tout,fastas,n_threads){

  #Retrieve all orphan sequences from previous steps.
  c(rownames(tout)[which(tout$Domain=='' & tout$Family=='')],
    names(fastas)[which(!names(fastas)%in%rownames(tout))])->orph.n

  orph.n -> orps

  res <- list()
  while (length(orps)>0){

    table(sapply(orps,function(x){strsplit(x,';')[[1]][1]})) -> tcont
    sort(tcont,decreasing = T) -> so

    if(length(so)>=3){
      tk <- 3
    }else{
      tk <- length(so)
    }

    names(so[1:tk]) -> ma
    unlist(sapply(ma,function(x){
      grep(paste0(x,';'),orps,fixed = T,value = T)
    })) -> fi

    if (length(fi)>=n_threads){
      parts <- n_threads
    }else{
      parts <- length(fi)
    }

    splitIndices(length(fi),parts) -> spin
    unlist(lapply(spin,function(x){
      tmp1<-tempfile()
      # write.fasta(lapply(fastas[fi[x]],memDecompress,'gzip',TRUE),
      #             names = fi[x],
      #             file.out = tmp1)
      write.fasta(fastas[fi[x]],
                  names = fi[x],
                  file.out = tmp1)
      tmp1
    })) -> temps

    tmp2 <- tempfile()
    # write.fasta(lapply(fastas[orps],memDecompress,'gzip',TRUE),
    #             names = orps,
    #             file.out = tmp2)
    write.fasta(fastas[orps],
                names = orps,
                file.out = tmp2)

    registerDoParallel(cores = n_threads)
    abc<-foreach(i=seq_along(spin),.combine = rbind)%dopar%{
      phmm.tmp<-tempfile()
      system2(command = "phmmer",
              args = c("-o /dev/null",
                       paste("--tblout",phmm.tmp),
                       "--cpu 0 --mx BLOSUM45",
                       temps[i],tmp2))

      #Return phmmer results
      outphmmer(pouti = phmm.tmp)
    }

    abc -> res[[length(res)+1]]

    #MCL
    tempfile() -> tmp
    write.table(abc[,c(1,2,4)],file = tmp,sep = "\t",
                quote = F,row.names = F,col.names = F)
    strsplit(runMCL(abc = tmp,neg.log10 = F,infl = 4),split = '\t') -> mclust

    which(sapply(mclust,length)==1) -> sin
    sapply(unlist(mclust[sin]),function(x){
      strsplit(x,';')[[1]][1]
    }) -> b
    which(!b%in%names(so)[1:3]) -> w
    if(length(w)==0){
      # mclust -> res[[length(res)+1]]
      orps[-which(orps%in%unlist(mclust))] -> orps
    }else{
      # mclust[-sin[w]] -> res[[length(res)+1]]
      orps[-which(orps%in%unlist(mclust[-sin[w]]))] -> orps
    }

  }

  do.call(rbind,res) -> res
  tempfile() -> restmp
  write.table(res[,c(1,2,4)],file = restmp,sep = '\t',
              quote = F,row.names=F,col.names = F)
  strsplit(runMCL(abc = restmp, neg.log10 = F,infl = 4),split = '\t') -> mclust
  mclust

}



#' @name outphmmer
#' @title Process phmmer output to make ir 'R'eadable.
#' @description Process phmmer output to make it readable by R. Taken and
#' adapted from \code{micropan} package (Lars Snipen and Kristian Hovde
#' Liland).
#' @param pouti \code{character}. phmmer output temporary file.
#' @return A \code{data.frame} with the phmmer output.
#' @author Ignacio Ferres
outphmmer<-function(pouti){
  readLines(pouti)->rl
  rl[which(!grepl("^\\#",rl))]->rl
  gsub("[ ]+"," ",rl)->rl
  strsplit(rl," ")->lst

  hit<-sapply(lst,function(x){x[1]})
  query<-sapply(lst,function(x){x[3]})
  eval<-as.numeric(sapply(lst,function(x){x[5]}))
  score<-as.numeric(sapply(lst,function(x){x[6]}))

  hmmer.table<-data.frame(Query=query,Hit=hit,
                          Evalue=eval,Score=score,
                          stringsAsFactors = F)
  return(hmmer.table)
}


#' @name runMCL
#' @title Run MCL
#' @description Cluster protein sequences by running the Markov Clustering
#' algorithm (MCL) over the phmmer comparison scores.
#' @param abc \code{character}. Temporary files of abc formatted phmmer
#' comparisons.
#' @param neg.log10 \code{logical}. Should the negative of the log10 value be
#' used as similarity measure?
#' @param infl \code{integer}. Inflacion value.
#' @return MCL output.
#' @author Ignacio Ferres
runMCL<-function(abc,neg.log10=TRUE,infl=6){

  if(neg.log10){
    arg <- paste("--abc --abc-neg-log10 -discard-loops y -abc-tf 'ceil(300)' -I",infl)
  }else{
    arg <- paste("--abc -discard-loops y -I",infl)
  }

  #deprecated
  #Run MCL
  # system2(command = "mcl",
  #         args = c(abc,arg,
  #                  "-q x -o -"),
  #         stdout = TRUE,
  #         stderr = FALSE)->rl
  #
  # rl

  #Run MCL
  tmpmcl <- tempfile()
  system(paste0('mcl ',
                abc,
                ' ',
                arg,
                ' -q x -o ',
                tmpmcl))
  readLines(tmpmcl) -> rl
  file.remove(tmpmcl)
  rl
}
