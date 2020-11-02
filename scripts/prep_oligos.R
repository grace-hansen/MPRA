#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(data.table)
library(biomaRt)
dir<-args[1]
bed<-args[2]
prefix<-args[3]

setwd(paste("~/projects/MPRA/",dir,sep=''))
rsids<-fread(bed,data.table=FALSE)
colnames(rsids)<-c("chr","start",'stop',"name","score","strand","source","alleles")
barcodes<-scan("../pseudorandom_barcodes.txt",what='character',sep='\n')
barcodes<-unique(barcodes)
barcodes<-barcodes[1:100000]

#Get sequence surrounding each variant
seq<-rsids
seq$chr<-paste("chr",seq$chr,sep='')
seq$start<-rsids$start-87
seq$stop<-rsids$start+88
write.table(seq,paste(prefix,"_getfasta.bed",sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

cmd=paste("bedtools getfasta -fi ~/midway/genos/hg38/GRCh38.primary_assembly.genome.fa -bed ",prefix,"_getfasta.bed -tab -fo temp",sep='')
system(cmd)            
cmd=paste("cat temp | tr ':' '\t' | tr '-' '\t' > ",prefix,"_seqs.bed",sep='')
system(cmd)
sequences<-fread(paste(prefix,"_seqs.bed",sep=''),data.table=FALSE)
colnames(sequences)<-c("chr",'start',"stop","seq")
sequences<-sequences[,c('chr','start','seq')]
seq<-inner_join(seq,sequences,by=c('chr','start'))
seq<-na.omit(seq,cols=alleles)


#Create sequence for each allele of each variant: exclude repeats >5bp, exclude indels
new_seqs<-data.table(matrix(nrow=0,ncol=10),stringsAsFactors = FALSE)
for (i in 1:nrow(seq)) {
  alleles<-strsplit(seq$alleles[i],'/')[[1]]
  if (sum(nchar(alleles)>1)==0 && sum(alleles=='-')==0) { #Only SNPs, no indels
    if (grepl("GGTACC",seq$seq[i])==FALSE && grepl("TCTAGA",seq$seq[i])==FALSE)
    for (j in 1:length(alleles)) {
      new_seq<-seq$seq[i]
      substr(new_seq,88,88)<-alleles[j]
      row<-c(as.character(seq[i,1:8]),alleles[j],new_seq)
      new_seqs<-rbind(new_seqs,t(row),make.row.names=FALSE)
    }
  }
}
colnames(new_seqs)<-c("chr","seq_start","seq_stop","rsid","score","strand","source","alleles","allele","enhancer_seq")
new_seqs<-inner_join(new_seqs,rsids[,c(4,2)],by=c("rsid"="name"))
colnames(new_seqs)[colnames(new_seqs) == "start"] <- "pos"
new_seqs<-new_seqs %>% dplyr::select(.,chr,seq_start,seq_stop,pos,everything())

#Remove rsids where RE binding site is created by one of the variants
new_seqs<-new_seqs[!(new_seqs$rsid %in% new_seqs$rsid[which(grepl("GGTACC",new_seqs$enhancer_seq))]),]
new_seqs<-new_seqs[!(new_seqs$rsid %in% new_seqs$rsid[which(grepl("TCTAGA",new_seqs$enhancer_seq))]),]

#Associate sequence to random barcodes
n_barcodes=floor(length(barcodes)/nrow(new_seqs))
write(paste(c("chr","seq_start","seq_stop","pos","rsid","source","alleles","allele","enhancer_seq","barcode","oligo"),collapse='\t'),paste(prefix,"_oligos.txt",sep=''))

row=1
for (i in 1:nrow(new_seqs)) {
  for (n in 1:n_barcodes) {
    barcode<-toupper(barcodes[row])
    oligo<-paste("ACTGGCCGCTTCACTG",new_seqs$enhancer_seq[i],"GGTACCTCTAGA",barcode,"AGATCGGAAGAGCGTCG",sep='')
    line<-paste(paste(new_seqs[i,c(1:5,8:11)],collapse='\t'),barcode,oligo,sep='\t')
    write(line,paste(prefix,"_oligos2.txt",sep=''),append=TRUE)
    row<-row+1
  }
}

cmd=paste("rm ",prefix,"_getfasta.bed ",prefix,"_seqs.bed",sep='')
system(cmd)

