# Title     : motif_filtering
# Objective : Filter motif for each condition.
# Created by: viyt
# Created on: 08.10.2019

###############################################################################
# 1) Required packages
###############################################################################
suppressMessages(suppressWarnings(require(ShortRead)))
suppressMessages(suppressWarnings(require(stringr)))

###############################################################################
# 2) Functions
###############################################################################
motif.finding <- function(names, motif.seq, r1, r2, MM=0, A="rev", B="fwd"){
  # names = {charater} Name of the motif (+ wt or scramble)
  # motif.seq = {Biostrings} complet motif sequence from xlsx file
  # r1 = {ShortReadQ} fastq file for reads of interest
  # r2 = {ShortReadQ} same as r1 but for r2
  # MM = {integer} Number of MissMatch allowed [default:0]
  # A = {string} Reverse or forward indication [default:fwd]
  # B = {string} Reverse or forward indication [default:fwd]

  # Remove barcode "taatag" and restriction site of complete motif sequence column's
  seq.m <- subseq(motif.seq, start = 13, end = -13)
  rc.motif <- reverseComplement(seq.m)

  # motif search (side A)
  if(A == "fwd"){
    # If motif on reverse reads
    motif.irA <- vmatchPattern(gsub("G", "R", rc.motif), sread(r1), fixed=FALSE, max.mismatch=MM)
    rA.i <- (elementNROWS(motif.irA) == 1)
  }else{
    # If motif on foward reads
    motif.irA <- vmatchPattern(gsub("C", "Y", rc.motif), sread(r2), fixed=FALSE, max.mismatch=MM)
    rA.i <- which(elementNROWS(motif.irA) == 1)
  }


  # motif search (side B)
  if(B == "fwd"){
    # If motif on reverse reads
    motif.irB <- vmatchPattern(gsub("C", "Y", seq.m), sread(r2), fixed=FALSE, max.mismatch=MM)
    rB.i <- which(elementNROWS(motif.irB) == 1)
  }else{
    #  If motif on foward reads
    motif.irB <- vmatchPattern(gsub("G", "R", seq.m), sread(r1), fixed=FALSE, max.mismatch=MM)
    rB.i <- which(elementNROWS(motif.irB) == 1)
  }

  # Number of motif find
  out.vec <- c(names, length(rA.i), length(rB.i)); print(out.vec)

  # Write output file for each side
  if(length(rA.i) != 0 && length(rB.i) == 0){
    r1.A <- r1[rA.i,]
    r2.A <- r2[rA.i,]
    r1.B <- NULL
    r2.B <- NULL
  }else if(length(rA.i) == 0 && length(rB.i) != 0){
    r1.A <- NULL
    r2.A <- NULL
    r1.B <- r1[rB.i,]
    r2.B <- r2[rB.i,]
  }else if(length(rA.i) != 0 && length(rB.i) != 0){
    r1.A <- r1[rA.i,]
    r2.A <- r2[rA.i,]
    r1.B <- r1[rB.i,]
    r2.B <- r2[rB.i,]
  }else{
    r1.A <- NULL
    r2.A <- NULL
    r1.B <- NULL
    r2.B <- NULL
  }

  rm(seq.m, rc.motif, motif.irA, motif.irB, rA.i, rB.i); gc()

  return(list("r1.A"=r1.A, "r2.A"=r2.A, "r1.B"=r1.B, "r2.B"=r2.B, "info"=out.vec))
}


###############################################################################
# 3) Mains
###############################################################################
## Steps 1 : Variable attribution
r1 <- readFastq(snakemake@input[["r1"]])
r2 <- readFastq(snakemake@input[["r2"]])
sA <- snakemake@params[["sideA"]]
sB <- snakemake@params[["sideB"]]
missma <- snakemake@params[["missmatch"]]
TF <- read.table(file=snakemake@input[["motifs"]], sep="\t", header=TRUE)
bc.name <- snakemake@input[["condition"]]
outfilename <- snakemake@ouput[["output"]]
system(paste0("mkdir -pv ", DIRout), ignore.stdout=TRUE)


## Steps 2 : Motif filter
m.table<-c("Motifs", "SideA", "SideB")

for(m in 1:length(TF[,1])){

  m.name <- as.character(TF[m,1])
  m.seq <- DNAString(as.character(TF[m,4])) # Tag skipping in the function

  m.f <- motif.finding(names=m.name, motif.seq=m.seq, r1=r1, r2=r2, MM=missma, A=sA, B=sB)

  m.table <- rbind(m.table, m.f$info)
  somme.m <- as.numeric(m.f$info[2]) + as.numeric(m.f$info[3])

  # Save fastq
  if(somme.m !=0 ){

    if(!is.null(m.f$r1.A) && is.null(m.f$r1.B)){
      writeFastq(object=m.f$r1.A, file=paste0(m.o, "/",bc.name,"_",m.name,".R1.fq.gz"), mode='w', full=FALSE, compress=TRUE)
      writeFastq(object=m.f$r2.A, file=paste0(m.o, "/",bc.name,"_",m.name,".R2.fq.gz"), mode='w', full=FALSE, compress=TRUE)

    }else if(is.null(m.f$r1.A) && !is.null(m.f$r1.B)){
      writeFastq(object=m.f$r1.B, file=paste0(m.o, "/",bc.name,"_",m.name,".R1.fq.gz"), mode='w', full=FALSE, compress=TRUE)
      writeFastq(object=m.f$r2.B, file=paste0(m.o, "/",bc.name,"_",m.name,".R2.fq.gz"), mode='w', full=FALSE, compress=TRUE)

    }else{# if(!is.null(m.f$r1.A) && !is.null(m.f$r1.B))
      # Side A
      writeFastq(object=m.f$r1.A, file=paste0(m.o, "/",bc.name,"_",m.name,"_A.R1.fq.gz"), mode='w', full=FALSE, compress=TRUE)
      writeFastq(object=m.f$r2.A, file=paste0(m.o, "/",bc.name,"_",m.name,"_A.R2.fq.gz"), mode='w', full=FALSE, compress=TRUE)
      # Side B
      writeFastq(object=m.f$r1.B, file=paste0(m.o, "/",bc.name,"_",m.name,"_B.R1.fq.gz"), mode='w', full=FALSE, compress=TRUE)
      writeFastq(object=m.f$r2.B, file=paste0(m.o, "/",bc.name,"_",m.name,"_B.R2.fq.gz"), mode='w', full=FALSE, compress=TRUE)
      # Merge reads for side A & B
      system(paste0("cat ",paste0(m.o, "/",bc.name,"_",m.name,"_A.R1.fq.gz")," ",paste0(m.o, "/",bc.name,"_",m.name,"_B.R1.fq.gz")," > ", paste0(m.o, "/",bc.name,"_",m.name,".R1.fq.gz")))
      system(paste0("cat ",paste0(m.o, "/",bc.name,"_",m.name,"_A.R2.fq.gz")," ",paste0(m.o, "/",bc.name,"_",m.name,"_B.R2.fq.gz")," > ", paste0(m.o, "/",bc.name,"_",m.name,".R2.fq.gz")))

    }
  }else{
    # nothing
  }
}

# Motifs
#colnames(m.table) <-  m.table[1,]; m.table<-m.table[2:nrow(m.table),]
#write.table(x=m.table, file=paste0("Summary_number_of_motifs.tsv"), sep="\t", row.names=FALSE)