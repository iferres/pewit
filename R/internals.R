#Internal functions.

#' @name runOnExit
#' @title Removes everything if an error occur
#' @description Removes the output folder if an error occurs after finnishing.
#' @param outdir The normalized path of the output directory
#' @return Warning message if exits before finnishing all the process.
#' @author Ignacio Ferres
runOnExit <- function(outdir){
  paste0(outdir,'pangenome.rds') -> fi
  if (!file.exists(fi)){
    warning('Something gone wrong: removing output directory.')
    unlink(outdir,recursive = TRUE)
  }
}


#' @name extractSeqsFromGff3
#' @title Extract Sequences from Gff3 file
#' @description Extract Sequences from Gff3 file and outputs a list of element,
#' one per each CDS. Each element contain the DNA gene sequence and the AA
#' protein sequence.
#' @param infile \code{character}. The gff3 filename.
#' @param in.path Where the output (fasta files) will be written.
#' @param keep.aa \code{logical}. Should AA sequences be kept in memory?
#' @return A \code{list} with CDS gene and protein sequences.
#' @author Gregorio Iraola and Ignacio Ferres.
#' @importFrom seqinr as.SeqFastadna s2c
extractSeqsFromGff3 <- function(infile,
                                in.path=NULL,
                                keep.aa=T){

  if(!dir.exists(in.path)){
    stop('Directory does not exist.')
  }

  readLines(infile) -> rl
  if(rl[1]!="##gff-version 3"){
    stop('One or more gff files seems not to be gff3 (version 3).')
  }

  if(grepl('#',infile)){
    warning('Conflictive character "#", will be substituted by "_".')
    gsub('#','_',infile) -> infile
  }

  #Save sequence in SeqFastadna object
  if(!any(grepl("##FASTA",rl))){
    stop("One or more gff files do(es) not contain the fasta genome sequences.")
  }
  grep("##FASTA",rl) + 1 -> fasta.ini
  length(rl) -> fasta.end
  rl[fasta.ini:fasta.end] -> fasta
  grep('>',fasta,fixed = T) -> fheaders
  matrix(nrow = length(fheaders),ncol = 2) -> t
  fheaders + 1 -> t[,1]
  c((fheaders - 1)[-1],length(fasta)) -> t[,2]
  apply(t,1,function(x){
    paste(fasta[x[1]:x[2]],collapse = '')
  }) -> ap
  names(ap) <- gsub('>','',fasta[fheaders])
  lapply(ap,function(x){
    as.SeqFastadna(s2c(x),name = names(x))
  }) -> fnas

  extractGffTable(rl = rl) -> gfftable

  gfftable[which(gfftable$Type=='CDS'),] -> cds

  apply(cds,1,function(x){
    getFfnFaa(fnas=fnas,
              contig = x[1],
              strand = x[9],
              from = as.numeric(x[7]),
              to = as.numeric(x[8]),id = x[2],
              product = x[5])
  }) -> fin
  nam <- sapply(fin,function(x){attr(x[[1]],'name')})
  if(any(grepl('#',nam))){
    warning('Conflictive character "#", will be substituted by "_".')
    gsub('#','_',nam) -> nam
  }
  names(fin) <- nam

  if(is.null(in.path)){
    fin
  }else{
    sapply(fin,function(x){x[[1]]}) -> ffn
    names(ffn) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(ffn))
    write.fasta(ffn,
                names = names(ffn),
                file.out = paste0(in.path,sub('.gff','.ffn',infile),collapse = '/'))
    sapply(fin,function(x){x[[2]]}) -> faa
    names(faa) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(faa))
    write.fasta(faa,
                names = names(faa),
                file.out = paste0(in.path,sub('.gff','.faa',infile),collapse = '/'))
    if (keep.aa){
      lapply(faa,function(x){paste0(x,collapse = '')}) -> faa
      faa
    }
  }

}

#' @name getFfnFaa
#' @title Extract CDS gene (DNA) and protein (AA) sequences
#' @description Extract CDS gene (DNA) and protein (AA) sequences from both a
#' gff3 table and the fasta sequence.
#' @param fnas A \code{list} of \code{SeqFastadna} genome sequences.
#' @param contig \code{character}. The name of the contig.
#' @param strand \code{character}. Either "+" or "-".
#' @param from \code{numeric}. The coordinate "from".
#' @param to \code{numeric}. The coordinate "to".
#' @param id \code{character}. The ID of the CDS.
#' @param product \code{character}. The gene product.
#' @return A \code{list} with two entries. The first is the CDS's DNA sequence,
#' and the last is the protein sequence. Object classes are \code{SeqFastadna}
#' and \code{SeqFastaAA}, respectively.
#' @author Ignacio Ferres
#' @importFrom seqinr getFrag.SeqFrag comp as.SeqFastadna as.SeqFastaAA translate
getFfnFaa <- function(fnas,contig,strand,from,to,id,product){
  seqinr::getFrag.SeqFrag(fnas[contig][[1]],begin = from,end = to) -> ffn
  attr(ffn,'begin') <- NULL
  attr(ffn,'end') <- NULL
  if(strand=='-'){
    rev(seqinr::comp(ffn)) -> ffn
  }
  seqinr::as.SeqFastadna(toupper(ffn),
                         name = id,
                         Annot = product) -> ffn
  seqinr::as.SeqFastaAA(toupper(seqinr::translate(ffn,numcode = 11)),
                        name = id,
                        Annot = product) -> faa

  list(ffn,faa) -> out
  out
}

#' @name extractGffTable
#' @title Extract the gff3 table and make it 'R'eadable.
#' @description Read a gff3 file and transforms the table in a \code{data.frame}.
#' @param rl \code{character}. A vector of character strings as passed by
#' \code{readLines()}, reading the gff3 file.
#' @return A \code{data.frame}.
#' @author Ignacio Ferres
extractGffTable <- function(rl){
  which(grepl('^\\#\\#',rl)) -> w
  rev(w)[1] - 1 -> upto
  rev(w)[2] + 1 -> from
  rl[from:upto] -> o

  strsplit(o,'\t') -> lst

  contig <- sapply(lst,function(x){x[1]})
  type <- sapply(lst,function(x){x[3]})
  from <- sapply(lst,function(x){x[4]})
  to <- sapply(lst,function(x){x[5]})
  strand <- sapply(lst,function(x){x[7]})
  phase <- sapply(lst,function(x){x[8]})
  attrib <- sapply(lst,function(x){x[9]})

  metadata <- strsplit(attrib,';')

  id <- sapply(metadata,function(x){
    gp<-grep('ID=',x,value = T)
    if(length(gp)>0){sub('ID=','',gp)}else{''}
  })
  locustag <- sapply(metadata,function(x){
    gp<-grep('locus_tag=',x,value = T)
    if(length(gp)>0){sub('locus_tag=','',gp)}else{''}
  })
  gene <- sapply(metadata,function(x){
    gp<-grep('gene=',x,value = T)
    if(length(gp)>0){sub('gene=','',gp)}else{''}
  })
  product <- sapply(metadata,function(x){
    gp<-grep('product=',x,value = T)
    if(length(gp)>0){sub('product=','',gp)}else{''}
  })

  out <- data.frame(Contig=contig,
                    ID=id,
                    LocusTag=locustag,
                    Gene=gene,
                    Product=product,
                    Type=type,
                    From=from,
                    To=to,
                    Strand=strand,
                    Phase=phase,
                    stringsAsFactors = F)
  out
}


#' @name getSeqOfType
#' @title Get sequences (nucleotidic or aminoacidic)
#' @description Extract cds either nucleotidic (gene) or aminoacidic (protein).
#' @param seqs \code{list} of sequences as passed by \link{extractSeqsFromGff3}.
#' @param type \code{character}. Either 'AA' or 'DNA'.
#' @return A \code{list} of sequences.
#' @author Ignacio Ferres
getSeqOfType <- function(seqs,type='AA'){
  names(seqs) -> n

  unlist(seqs,recursive = F) -> un
  if (type=='AA'){
    lapply(un,function(x){x[[2]]}) -> lap
  }else{
    lapply(un,function(x){x[[1]]}) -> lap
  }
  names(lap) <- sub('.gff.',';',names(lap))
  lap
}

#Distributes all proteins in many files as n_threads set in order to optimize
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
splitAndWriteFastas<-function(fastas,n_threads){

  splitIndices(length(fastas),n_threads)->spi
  names(spi)<-seq_along(spi)
  tdir<-tempdir()
  registerDoParallel(cores = n_threads)
  temps<-foreach(i=seq_along(spi))%dopar%{
    tpo<-tempfile(pattern = "tmpo",tmpdir = tdir,fileext = ".fasta")
    # lapply(fastas[spi[[i]]],function(x){
    #   memDecompress(x,type = 'gzip',asChar = T)
    #   }) -> sq
    # names(sq) -> nsq
    # write.fasta(sq,
    #             names = nsq,
    #             file.out = tpo)
    write.fasta(fastas[spi[[i]]],
                names=names(fastas[spi[[i]]]),
                file.out = tpo)
    tpo
  }

  temps<-unlist(temps)
  return(temps)
}

#' @name runHmmsearch
#' @title Run Hmmsearch (HMMER 3)
#' @description Takes a fasta file and a Hidden Markov Model profile and
#' performs a search of the former over the latter.
#' @param fasta A protein fasta file.
#' @param pfam A hmm file. Must be pressed (see hmmpress from HMMER manual).
#' @param n_threads An \code{integer}. The number of cores to use.
#' @return The path to a temporary file where the hmmsearch output is placed.
#' @author Ignacio Ferres
runHmmsearch <- function(fasta,
                         pfam,
                         n_threads=1L){
  #run hmmsearch
  tempfile(pattern = 'tmpo',fileext = '.tab') -> domtblout
  paste0('hmmsearch -o /dev/null --domtblout ',
         domtblout,
         ' --noali --cut_ga --cpu ',
         as.character(n_threads),
         ' ',
         pfam,
         ' ',
         fasta) -> pfm
  system(pfm)
  file.remove(fasta)
  domtblout

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
outhmmsearch<-function(pouti,ref){
  readLines(pouti)->rl
  rl[which(!grepl("^\\#",rl))]->rl
  gsub("[ ]+"," ",rl)->rl
  strsplit(rl," ")->lst

  query<-sapply(lst,function(x){x[1]})
  hit<-sapply(lst,function(x){x[4]})
  pfmID<-sapply(lst,function(x){x[5]})
  eval<-as.numeric(sapply(lst,function(x){x[13]}))
  score<-as.numeric(sapply(lst,function(x){x[14]}))
  st<-as.numeric(sapply(lst, function(x){x[18]}))
  en<-as.numeric(sapply(lst,function(x){x[19]}))
  desc<-sapply(lst,function(x){paste(x[23:length(x)],collapse = " ")})

  hmmer.table<-data.frame(Query=query,Hit=hit,PfamID=pfmID,Evalue=eval,Score=score,
                          Start=st,End=en,Description=desc,stringsAsFactors = F)
  hmmer.table<-hmmer.table[-which(hmmer.table$Evalue>0.1),]

  sub("\\.\\d+","",ref$ID) -> codpf
  sub("\\.\\d+","",hmmer.table$PfamID) -> codhit
  lapply(codhit,function(x){which(codpf==x)})->ind
  unlist(ind) -> ind
  ref$TP[ind]->hmmer.table$Type
  ref$CL[ind]->hmmer.table$Clan
  return(hmmer.table)
}

#' @name processPfam_A_Dat
#' @title Process Pfam-A.dat file
#' @description Process Pfam-A.dat file.
#' @param datPfam \code{character}. The path to Pfam-A.dat file.
#' @param n_threads \code{integer}. The number of threads to use.
#' @return A \code{data.frame} with information of each Pfam entry.
#' @author Ignacio Ferres
#' @importFrom parallel mclapply
processPfam_A_Dat<-function(datPfam,n_threads){
  readLines(datPfam,skipNul = T)->rl
  grep("STOCKHOLM 1.0",rl)+1 -> pr
  grep("//",rl)-1 -> fn
  cbind(pr,fn)->pa
  apply(pa,1,function(x){rl[x[1]:x[2]]})->li
  do.call(rbind,mclapply(li,function(y){
    c(regmatches(y[grep('GF AC',y)],regexpr("PF\\d+[.]*\\d{0,1}",y[grep('GF AC',y)])),
      regmatches(y[grep('GF TP',y)],regexpr("Domain|Family|Motif|Repeat|Disordered|Coiled-coil",y[grep('GF TP',y)])),
      if(length(grep('GF CL',y))!=0){regmatches(y[grep('GF CL',y)],regexpr("CL\\d{4}",y[grep('GF CL',y)]))}else{'No-Clan'})
  },mc.cores = n_threads))->ref
  colnames(ref)<-c("ID","TP","CL")
  ref<-data.frame(ref,stringsAsFactors=F)
  return(ref)
}

#' @name outhmmscan
#' @title Process hmmscan output to make it 'R'eadable.
#' @description Process hmmscan output to make it readable by R.
#' @param pouti \code{character}. hmmscan output temporary file.
#' @param ref \code{data.frame}. Information about each Pfam-A entry.
#' @return A \code{data.frame} with the hmmscan output plus information about
#' Pfam-A hits.
#' @note Taken and adapted from \code{micropan} package (Lars Snipen and
#' Kristian Hovde Liland).
#' @author Ignacio Ferres
outhmmscan<-function(pouti,ref){
  readLines(pouti)->rl
  rl[which(!grepl("^\\#",rl))]->rl
  gsub("[ ]+"," ",rl)->rl
  strsplit(rl," ")->lst

  hit<-sapply(lst,function(x){x[1]})
  pfmID<-sapply(lst,function(x){x[2]})
  query<-sapply(lst,function(x){x[4]})
  eval<-as.numeric(sapply(lst,function(x){x[13]}))
  score<-as.numeric(sapply(lst,function(x){x[14]}))
  st<-as.numeric(sapply(lst, function(x){x[18]}))
  en<-as.numeric(sapply(lst,function(x){x[19]}))
  desc<-sapply(lst,function(x){paste(x[23:length(x)],collapse = " ")})

  hmmer.table<-data.frame(Query=query,Hit=hit,PfamID=pfmID,Evalue=eval,Score=score,
                          Start=st,End=en,Description=desc,stringsAsFactors = F)
  hmmer.table<-hmmer.table[-which(hmmer.table$Evalue>0.1),]

  sub("\\.\\d+","",ref$ID) -> codpf
  sub("\\.\\d+","",hmmer.table$PfamID) -> codhit
  sapply(codhit,function(x){which(codpf==x)})->ind
  ref$TP[ind]->hmmer.table$Type
  ref$CL[ind]->hmmer.table$Clan
  return(hmmer.table)
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


#' @name determineOverlap
#' @title Determine which Pfam domains overlaps.
#' @description Takes an overlappling-domain logical matrix \code{m} and
#' through a 'trace-forward' precedure from \code{m[1,1]} to \code{m[ ,ncol(m)]}
#' outputs a list with overlappling-domains row names.
#' @param m A logical \code{matrix} indicating if 'x' domain overlaps with 'y'.
#' @return A \code{list} of overlapping domains.
#' @author Ignacio Ferres
determineOverlap<-function(m){
  st<-c(1,1)
  li<-list()
  #start from m[1,1]
  li[[1]]<-colnames(m)[st[2]]
  #while 'st' doesn't reach the right edge of the matrix...
  while(st[2]!=ncol(m)){
    #If next (right) cell is TRUE..
    if(m[st[1],st[2]+1]==T){
      #..then move to next cell.
      st[2]<-st[2]+1
      li[[length(li)]][length(li[[length(li)]])+1]<-colnames(m)[st[2]]
      #If next cell is FALSE but bottom cell is TRUE...
    }else if(m[st[1]+1,st[2]]==T){
      #..then move to bottom cell.
      st[1]<-st[1]+1
      #If both next and bottom cells are FALSE...
    }else{
      #..then move to the bottom-right (diagonal) cell.
      st[1]<-st[1]+1
      st[2]<-st[2]+1
      li[length(li)+1]<-colnames(m)[st[2]]
    }
  }
  return(li)
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
rm_overlaping_clans<-function(tping,type){

  #Just consider one type
  # tpingo<-tping[which(tping$Type==type),,drop=F]
  tpingo<-as.data.frame(tping)[which(tping["Type"]==type),]

  if (dim(tpingo)[1]!=0){
    #Subset by clan
    lping<-split.data.frame(tpingo,as.factor(tpingo$Clan))

    for (i in 1:length(lping)){

      if (dim(lping[[i]])[1]>1){

        #Order by start (this is in order to allow determineOverlap() to work)
        lping[[1]][order(lping[[1]]$Start),]->lping[[1]]

        #Create list with sequence positions
        unlist(apply(lping[[i]],1,function(p){
          list(seq(as.numeric(p[6]),as.numeric(p[7]),1))
        }),recursive = F)->inc

        #Logical matrix of (non-)overlaping PF -determines shared positions-.
        matrix(unlist(lapply(inc,function(y){
          lapply(inc,function(p){any(p%in%y)})
        })),nrow = length(inc),
        dimnames = list(rownames(lping[[i]]),rownames(lping[[i]])))->m

        #Output list of overlapping domains (row names)
        determineOverlap(m = m)->ovlp


        #Process those with same-clan overlapping domains
        which(sapply(ovlp,length)>1)->multiclan
        if(length(multiclan)>=1){
          dej<-list()
          for (k in seq_along(multiclan)){
            which.min(lping[[i]][ovlp[[multiclan[k]]],4])->indx#columna de evalue
            #Remove the rest
            dej[[k]]<- lping[[i]][ovlp[[multiclan[k]]][indx],,drop=F]
          }
          do.call(rbind,dej)->dej
          do.call(rbind,list(lping[[i]][unlist(ovlp[-multiclan]),],dej))->lping[[i]]
        }


      }

    }
    do.call(rbind,lping)->lping

  }else{
    lping<-NULL
  }
  return(lping)
}

#' @name grouping
#' @title Split hmmscan output by query.
#' @description Split hmmscan output by query.
#' @param t \code{data.frame}. hmmscan output passed by \link{outhmmscan}.
#' @param singlehit \code{logical}. Has this query just one hit?
#' @param type \code{character}. One of 'Domain', 'Family', 'Motif', 'Repeat',
#' 'Disordered' or 'Coiled-coil'.
#' @return A \code{list} of \code{data.frame} with hits by query.
#' @author Ignacio Ferres
grouping<-function(t,singlehit=FALSE,type="Domain"){
  if (!singlehit){
    split.data.frame(t,t$Query)->x
    unlist(lapply(x,function(x){nrow(x)>1}))->y
    return(x[which(y)])
  }else{
    split.data.frame(t,t$Query)->x
    which(unlist(lapply(x,function(z){nrow(z)==1})))->y
    which(unlist(lapply(x[y],function(z){z$Type==type})))->j
    return(x[y[j]])
  }
}


#' @name processHmmsearch
#' @title Main hmmsearch output processing function
#' @description Outputs a data.frame with the domain structure of each sequence
#' (rows). It has 6 columns, each one with the semi-colon separated found
#' sequences of Pfam (PF) entries of the following 6 categories: 'Domain',
#' 'Family', 'Motif', 'Repeat','Disordered' or 'Coiled-coil'.
#' @param pout \code{character}. Temporary file name (hmmscan output).
#' @param ref \code{data.frame}. Information about each Pfam-A entry.
#' @return A \code{data.frame} with the domain structure of each sequence.
#' @author Ignacio Ferres
processHmmsearch <- function(pout,ref){
  writeLines('Loading hmmsearch output..')
  # outhmmscan(pouti = pout,ref = ref)->t
  outhmmsearch(pout,ref) -> t


  #### Keep single hit ####

  #Single hit (Domain)
  grouping(t = t,singlehit = TRUE,type = "Domain")->dom.sing
  #Single hit (Family)
  grouping(t = t,singlehit = TRUE,type = "Family")->fam.sing
  #Single hit (Motif)
  grouping(t = t,singlehit = TRUE,type = "Motif")->mot.sing
  #Single hit (Repeat)
  grouping(t = t,singlehit = TRUE,type = "Repeat")->rep.sing
  #Single hit (Disordered)
  grouping(t = t,singlehit = TRUE,type = "Disordered")->dis.sing
  #Single hit (Coiled-coil)
  grouping(t = t,singlehit = TRUE,type = "Coiled-coil")->cco.sing



  #### Group by query those which have more than one hit ####

  grouping(t = t,singlehit = FALSE)->gping


  #### Remove same-clan overlapping hits ####
  # cat('Removing same-clan overlapping domains..')

  #Domain
  lapply(gping,function(x){
    rm_overlaping_clans(tping=x,type = "Domain")
  })->dom.novlp
  #Family
  lapply(gping,function(x){
    rm_overlaping_clans(tping=x,type = "Family")
  })->fam.novlp
  #Motif
  lapply(gping,function(x){
    rm_overlaping_clans(tping=x,type = "Motif")
  })->mot.novlp
  #Repeat
  lapply(gping,function(x){
    rm_overlaping_clans(tping=x,type = "Repeat")
  })->rep.novlp
  #Disordered
  lapply(gping,function(x){
    rm_overlaping_clans(tping=x,type = "Disordered")
  })->dis.novlp
  #Coiled-coil
  lapply(gping,function(x){
    rm_overlaping_clans(tping=x,type = "Coiled-coil")
  })->cco.novlp
  # cat(' DONE!')

  #### Merge ####

  #List of lists
  list(dom.sing,fam.sing,mot.sing,
       rep.sing,dis.sing,cco.sing,
       dom.novlp,fam.novlp,mot.novlp,
       rep.novlp,dis.novlp,cco.novlp)->lpm
  #Index of non-empthy lists
  which(unlist(lapply(lpm,function(x){
    !is.null(unlist(x))
  })))->ind.lpm

  #Filtering done!
  do.call(rbind,unlist(lpm[ind.lpm],recursive = F))->t
  t[order(t$Query),]->t
  rownames(t)<-seq_len(nrow(t))


  #### Build table (data.frame) with domain structure of each secuence ####

  tout<-as.data.frame(matrix(nrow = length(unique(t$Query)),ncol = 6),
                      row.names = unique(as.character(t$Query)))
  colnames(tout)<-c('Domain',
                    'Family',
                    'Motif',
                    'Repeat',
                    'Disordered',
                    'Coiled-coil')

  tdom<-t[which(t$Type=="Domain"),]
  tfam<-t[which(t$Type=="Family"),]
  tmot<-t[which(t$Type=="Motif"),]
  trep<-t[which(t$Type=="Repeat"),]
  tdis<-t[which(t$Type=="Disordered"),]
  tcoi<-t[which(t$Type=="Coiled-coil"),]

  #Dom
  do.call(rbind,lapply(unique(t$Query),function(x){
    paste(tdom$Hit[grep(x,tdom$Query,fixed = T)],collapse =  ";")
  }))[,1]->tout$Domain
  #Fam
  do.call(rbind,lapply(unique(t$Query),function(x){
    paste(tfam$Hit[grep(x,tfam$Query,fixed = T)],collapse =  ";")
  }))[,1]->tout$Family
  #Mot
  do.call(rbind,lapply(unique(t$Query),function(x){
    paste(tmot$Hit[grep(x,tmot$Query,fixed = T)],collapse =  ";")
  }))[,1]->tout$Motif
  #Rep
  do.call(rbind,lapply(unique(t$Query),function(x){
    paste(trep$Hit[grep(x,trep$Query,fixed = T)],collapse =  ";")
  }))[,1]->tout$Repeat
  #Dis
  do.call(rbind,lapply(unique(t$Query),function(x){
    paste(tdis$Hit[grep(x,tdis$Query,fixed = T)],collapse =  ";")
  }))[,1]->tout$Disordered
  #Coi
  do.call(rbind,lapply(unique(t$Query),function(x){
    paste(tcoi$Hit[grep(x,tcoi$Query,fixed = T)],collapse =  ";")
  }))[,1]->tout$`Coiled-coil`



  return(tout)
}

#' @name align
#' @title Multiple alignment with mafft.
#' @description Takes a list of \code{SeqFastaAA} objects from \link{seqinr}
#' package and performs a multiple alignment with mafft aligner. Outputs a
#' \code{matrix} with columns as positions in the alignment and rows as
#' sequences.
#' @param rf A \code{list} of \code{SeqFastaAA} objects.
#' @param type \code{character}. One of 'AA' or 'DNA'.
#' @param n_threads \code{integer}. The number of threads to use.
#' @param accu \code{logical}. Accurate alignment?
#' @param mxit1000 \code{logical}. Set \code{mafft} parameter
#' \code{--maxiterate 1000} .
#' @return A \code{matrix} with columns as positions in the alignment and rows
#' as sequences.
#' @author Ignacio Ferres
#' @references Katoh, Misawa, Kuma, Miyata (2002). MAFFT: a novel method for
#' rapid multiple sequence alignment based on fast Fourier transform. \emph{Nucleic
#' Acids Res} \strong{30}:3059-3066.
#' @importFrom seqinr write.fasta
align<-function(rf,type='AA',n_threads=1L,accu=TRUE,mxit1000=TRUE){
  tempfile(fileext = ".faa")->tmp
  write.fasta(rf,names = names(rf),file.out = tmp)

  type <- match.arg(type,c('AA','DNA'))

  if(type=='AA'){

    if(accu){
      args <- c("--quiet --bl 30","--globalpair")
      if(mxit1000){
        args <- c(args,"--maxiterate 1000")
      }

    }else{
      if(length(rf)<2000){
        args <- c("--quiet --bl 30")
      }else{
        args <- c('--quiet --bl 30 --retree 1 --maxiterate 0')
      }

    }

  }else{
    if(accu){
      args <- c("--quiet","--globalpair")
      if(mxit1000){
        args <- c(args,"--maxiterate 1000")
      }

    }else{
      if(length(rf)<2000){
        args <- c("--quiet")
      }else{
        args <- c('--quiet --retree 1 --maxiterate 0')
      }
    }
  }


  if (n_threads>1){
    args <- c(args,"--thread",n_threads)
  }
  args <- c(args,tmp)


  #RUN MAFFT
  if (length(rf) <= 100){

    system2(command = "mafft",
            args = args,
            stdout = T)->rl
  }else{

    tmpo <-tempfile()
    paste0('mafft ',paste0(args,collapse = ' '),' > ',tmpo) -> run
    system(run)
    readLines(tmpo) -> rl
    file.remove(tmpo)

  }

  #Transform alignment in matrix
  grep(">",rl)->gp
  cbind(gp+1,c(gp[-1]-1,length(rl)))->coorsq
  strsplit(apply(coorsq,1,function(x){
    paste(toupper(rl[x[1]:x[2]]),collapse = '')
  }),split = '')->sqs
  names(sqs)<-sub(">","",rl[gp])
  do.call(rbind,sqs)->m
  file.remove(tmp)
  #apply(m,1:2,charToRaw)->m
  return(m)
}


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
#' as passed by \link{processHmmscan}.
#' @param fastas A \code{list} of protein sequeces as \code{SeqFastaAA} objects
#' @param n_threads \code{integer} Number of threads to use.
#' @return A \code{list} of protein clusters.
#' @author Ignacio Ferres
#' @importFrom parallel splitIndices
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach '%dopar%'
#' @importFrom seqinr write.fasta
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


#' @name checkIfParalogues
#' @title Check if cluster contains paralogues
#' @description Check if a cluster contains paralogues.
#' @param p A vector of proteins
#' @return \code{logical}, if contains paralogues (i.e. proteins of the same
#' organism).
#' @author Ignacio Ferres
checkIfParalogues<-function(p){
  sapply(p,function(y){strsplit(y,";")[[1]][1]})-> st
  if (any(duplicated(st))){
    TRUE
    if (length(unique(st))==1){
      FALSE
    }else{
      TRUE
    }
  }else{
    FALSE
  }
}


#' @name pdist.aa
#' @title Calculate p distance from a multiple alignment
#' @description Calculate p distance from a multiple alignment.
#' @param x \code{matrix}. A multiple alignment as passed by
#' \link{align}.
#' @param scaled \code{logical}. Return scaled by alignment length?
#' @return A \code{dist} object.
#' @note Taken and adapted from \code{\link[ape]{dist.aa}}, from ape package by
#' Paradis et al.
#' @author Ignacio Ferres
pdist.aa <- function (x, scaled = TRUE) {
  apply(x,1:2,function(i){charToRaw(toupper(i))}) -> x
  n <- nrow(x)
  d <- numeric(n * (n - 1)/2)
  X <- charToRaw("-")
  k <- 0L
  for (i in 1:(n - 1)) {
    a <- x[i, ]
    for (j in (i + 1):n) {
      b <- x[j, ]
      del <- a == X & b == X
      p <- length(b <- b[!del])
      tmp <- sum(a[!del] != b)
      k <- k + 1L
      d[k] <- if (!scaled){
        tmp
      }else{ tmp/p }
    }
  }
  attr(d, "Size") <- n
  attr(d, "Labels") <- rownames(x)
  attr(d, "Diag") <- attr(d, "Upper") <- FALSE
  attr(d, "call") <- match.call()
  class(d) <- "dist"
  d
}

#' @name splitClusters
#' @title Split clusters
#' @description Takes coarse protein groups clustered by ethier domain
#' structure or by phmmer comparison plus MCL, and split them in true
#' orthologues groups.
#' @param clstr A \code{list} of \code{SeqFastaAA} objects.
#' @details A multiple alignment is performed by mafft aligner. The alignment
#' is trimmed and a p-distance matrix is calculated. A neighbour joining (NJ)
#' tree is inferred, and midpoint rooted.
#'
#' Recent paralogues are detected by identifying those nodes which only contain
#' tips from the same organism, and just the one with the smaller distance to
#' its neighbours is kept.
#'
#' Maximum subtrees which contain at most one tip of each organism are
#' considered true orthologues.
#' @return A \code{list} of true orthologues protein clusters.
#' @author Ignacio Ferres
#' @importFrom phangorn midpoint
#' @importFrom ape njs subtrees
splitClusters <- function(clstr){
  #Build gene nj tree
  if (length(clstr)>600){
    accu <- FALSE
  }else{
    accu <- TRUE
  }
  align(rf = clstr,type = 'AA',accu = accu) -> al
  #Trim alignment
  round(nrow(al)*0.8) -> trtr
  which(apply(al,2,function(x){length(which(x=='-'))})>=trtr) -> trimal
  if (length(trimal)!=0){
    al[,-trimal] -> al
  }

  #Calculate p distances
  pdist.aa(al) -> d
  #Infer NJ tree
  try(phangorn::midpoint(ape::njs(d)),silent = T) -> tree
  if(class(tree)=='try-error'){
    ape::njs(d) -> tree
  }


  #Extract true orthologues
  ape::subtrees(tree) -> l
  whichNodesHaveRecentParalogues(subtrees = l) -> w
  if(!is.null(w)){
    #..and then remove them leaving just one (the closer one).
    rmvRecentParalogues(njtree = tree,w = w,d = d) -> tree.trim
    ape::subtrees(tree.trim$tree) -> l
  }
  whichSubtreesAreTrueOrthologues(subtrees = l) -> tru.or

  maxTrueOGs(tru.or = tru.or,l = l) -> tru.or2
  sapply(l[which(tru.or2)],function(x){x$tip.label},simplify = F) -> fin

  if(!is.null(w)){
    #Some singletones may have out of 'fin' because they doesn't form a
    # subtree, so..
    sapply(names(tree.trim$paralogs),function(x){any(grep(x,fin))}) -> sap
    if(length(which(!sap))>0){
      as.list(names(which(!sap))) -> sng
      append(fin,sng) -> fin
    }
    # attr(fin,'paralogues') <- tree.trim$paralogs
    sapply(names(tree.trim$paralogs),function(x){grep(x,fin)}) -> where
    for (i in 1:length(fin)){
      if (any(i==where)){
        which(names(tree.trim$paralogs)%in%fin[[i]]) -> prl
        attr(fin[[i]],'paralogues') <- as.vector(unlist(tree.trim$paralogs[prl]))
      }
    }

  }
  # cat('.')
  fin
}

#' @name whichNodeHaveRecentParalogues
#' @title Which nodes have recent paralogues
#' @description Detect which nodes contain only tips from the same organism.
#' @param subtrees A \code{list} of trees of class "phylo".
#' @return A \code{vector} of nodes that contain recen paralogues.
#' @author Ignacio Ferres
whichNodesHaveRecentParalogues <- function(subtrees){
  unlist(lapply(subtrees,function(x){
    sapply(x$tip.label,function(y){strsplit(y,';')[[1]][1]}) -> s
    if(length(unique(s))==1 & length(s)>1){
      x$name
    }
  }))
}


#' @name rmvRecentParalogues
#' @title Remove recent paralogues
#' @description If recent paralogues are detected, just the one with the closer
#' distance to the tips which are descendants of the immediate previous node,
#' is kept.
#' @param njtree A tree of class \code{phylo}.
#' @param w A \code{vector} of nodes which contain recent paralogues, as passed
#' by \link{whichNodesHaveRecentParalogues}.
#' @param d A \code{dist} object with the pairwise p-distance between proteins
#' as passed by \link{pdist.aa}.
#' @return A \code{list} with two elements: A tree of class \code{phylo} with
#' paralogues removed, and a \code{list} of paralogues removed.
#' @author Ignacio Ferres
#' @importFrom phangorn Descendants Ancestors
#' @importFrom ape drop.tip
rmvRecentParalogues <- function(njtree,w,d){
  #Removes subtrees and nodes contained in other subtrees
  sapply(w,function(x){
    phangorn::Descendants(njtree,x,type='tips')
  }) -> tps
  if(length(w)>1){
    rmv <- c()
    for (i in 1:length(tps)){
      any(sapply(tps[-i],function(x){all(tps[[i]]%in%x)})) -> rr
      rmv <- c(rmv,rr)
    }

    if(any(rmv==TRUE)){
      w[-which(rmv)] -> w
      tps[-which(rmv)] -> tps
    }
  }


  as.matrix(d) -> d2

  #Removes recent paralogues until just one is kept in each case
  keep <- c()
  paralogs <- list()
  rem <- c()
  for (i in 1:length(w)){
    njtree$tip.label[tps[[i]]] -> gns
    #Ancestral node
    phangorn::Ancestors(njtree,w[i],'parent') -> nd
    #Descendantes of the ancestral node
    phangorn::Descendants(njtree,nd,'tips') -> dsc
    #Genes of the closest node (descendants from the parent node)
    njtree$tip.label[dsc[[1]][-which(dsc[[1]]%in%tps[[i]])]] -> bros
    #Kept the one with the smaller distance between it and their brothers
    which.min(apply(d2[gns,bros,drop=F],1,sum)) -> k
    paralogs[[i]] <- gns[-k]
    keep <- c(keep,k)
    rem <- c(rem,gns[-k])
  }
  #list of paralogues. The chosen one is keept as element list names
  names(paralogs) <- names(keep)
  #Descard the rest of the tips (recent paralogues)
  ape::drop.tip(njtree,rem) -> njtree
  list(njtree,paralogs) -> out
  names(out) <- c('tree','paralogs')
  return(out)
}


#' @name whichSubtreesAreTrueOrthologues
#' @title Which subtrees are true orthologues
#' @description Which subtrees are true orthologues, i.e. which nodes contain
#' at most only one tip of each organism.
#' @param subtrees A \code{list} of trees of class \code{phylo}.
#' @return A \code{vector} of \code{logical} values indicating if subtrees
#' contain at most one tip of each organism.
#' @author Ignacio Ferres
whichSubtreesAreTrueOrthologues<-function(subtrees){
  res<-c()
  for (i in 1:length(subtrees)){
    sapply(subtrees[[i]]$tip.label,function(x){strsplit(x,';')[[1]][1]}) -> p
    table(p) -> t
    if(all(t==1)){
      res<-c(res,TRUE)
    }else{
      res<-c(res,FALSE)
    }
  }
  res
}


#' @name maxTrueOGs
#' @title Maximum subtrees which contains only true orthologues
#' @description Detect maximum subtrees which contains only true orthologues.
#' @param tru.or A \code{vector} of logical values.
#' @param l A \code{list} of trees of class \code{phylo}.
#' @return A logical \code{vector} indicating if true orthologue subtree is
#' maximum or not.
#' @author Ignacio Ferres
maxTrueOGs <- function(tru.or,l){
  res2<-tru.or
  for(i in 1:length(tru.or)){
    if(tru.or[i]==FALSE){
      next
    }else{
      l[[i]]$tip.label -> tps
      res2[i] <- FALSE
      sapply(l[which(res2)],function(x){x$tip.label},simplify = F)->s
      any(sapply(s,function(x){all(tps%in%x)})) -> r
      if(!r){
        res2[i] <- TRUE
      }
    }
  }
  res2
}

#' @name setClusterNames
#' @title Set cluster names
#' @description Creates othologue cluster names.
#' @param final.clusters A \code{list} of protein clusters.
#' @return A \code{vector} of names for the clusters
#' @author Ignacio Ferres
setClusterNames<-function(final.clusters){
  np<-nchar(as.character(length(final.clusters)))
  np<-paste("%0",np,"d",sep="")
  npnam<-paste("OG",sprintf(np,1:length(final.clusters)),sep="")
  npnam
}


#' @name writeClusters
#' @title Write Clusters
#' @description Write clusters in a fasta-like format.
#' @param outdir Output directory
#' @param final.clusters A list of clusters
#' @return A file 'clusters.txt' in \code{outdir}.
#' @author Ignacio Ferres
writeClusters <- function(outdir,final.clusters){
  sink(paste0(outdir,'clusters.txt'))
  for(i in 1:length(final.clusters)){
    attr(final.clusters[[i]],'pfamStr') -> pfm
    if(length(pfm)>0){
      paste('[',paste(pfm,collapse = ';'),']',sep = '') -> pfm
    }else{
      '' -> pfm
    }
    cat(paste0('>',names(final.clusters[i]),' ',pfm,'\n'))
    cat(final.clusters[[i]],sep = ' ',fill = F);cat('\n')
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
writeParalogues <- function(outdir,final.clusters){
  sink(paste0(outdir,'paralogues.txt'))
  for(i in 1:length(final.clusters)){
    attr(final.clusters[[i]],'pfamStr') -> pfm
    if(length(pfm)>0){
      paste('[',paste(pfm,collapse = ';'),']',sep = '') -> pfm
      if(!is.null(attributes(final.clusters[[i]])$paralogues)){
        cat(paste0('>',names(final.clusters[i]),' ',pfm,'\n'))
        cat(attr(final.clusters[[i]],'paralogues'),sep = ' ',fill = F);cat('\n')
      }
    }else{
      '' -> pfm
      if(!is.null(attributes(final.clusters[[i]])$paralogues)){
        cat(paste0('>',names(final.clusters[i]),' ',pfm,'\n'))
        cat(attr(final.clusters[[i]],'paralogues'),sep = ' ',fill = F);cat('\n')
      }
    }

  }
  sink()
}
