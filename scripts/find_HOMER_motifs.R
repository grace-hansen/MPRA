#!/usr/bin/Rscript
if (length(args)==0) {
  stop("USAGE: MPRA_motifenrichment \n
       This script runs HOMER to perform motif enrichment on the results of the MPRA analysis. \n
       It runs HOMER in two ways: \n
       1) performing motif enrichment on significant enhancers, using nonsignificant enhancers as background, \n
       2) performing motif enrichment on significant enhancers, using the genome as background, \n
       3) performing motif enrichment on EMVars, using nonsignificant enhancers as background. \n"
       , call.=FALSE)
}
args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(data.table)
setwd("~/projects/MPRA/WHR/results/")

### STEPS ###
#1) Extract rsids for nonsignificant enhancers, significant enhancers, and significant EMVars
#2) Make bed files for '' ''
#3) Expand beds to 175bp, merge to collapse neighboring rsids into single locus
#4) Run HOMER to identify motifs:
#a) Significant enhancers vs. nonsignificant enhancers
#b) Significant enhancers vs. genome (unmasked if enhancers are in repeats)
#c) Significant EMVars vs nonsignificant enhancers

#Make beds of EMVars, significant enhancers, and nonsignificant enhancers.
seqs<-as.data.frame(fread("~/projects/MPRA/WHR/variants/WHR_adipose_subcutaneous_rsids_oligos.txt"),stringsAsFactors=FALSE)
barcodes<-as.data.frame(fread("~/projects/MPRA/WHR/results/barcode_activity.txt"),stringsAsFactors=FALSE)
seqs<-seqs[seqs$barcode %in% barcodes$barcode,]

enhancers<-as.data.frame(fread("~/projects/MPRA/WHR/results/sig_enhancers.bed"),stringsAsFactors=FALSE)
colnames(enhancers)<-c("chr","pos","pos1","rsid","dot","dot1","source")
enhancers<-enhancers %>% distinct(rsid,.keep_all=TRUE)

EMVars<-as.data.frame(fread("~/projects/MPRA/WHR/results/sig_EMVars.bed"),stringsAsFactors=FALSE)
colnames(EMVars)<-c("chr","pos","pos1","rsid","dot","dot1","source")
EMVars<- EMVars %>% distinct(rsid,.keep_all=TRUE)

nonsig_rsids<-barcodes[!(barcodes$rsid %in% enhancers$rsid),] %>% distinct(rsid,.keep_all=TRUE)

EMVar_seqs<-merge(seqs,EMVars,by=c("rsid","chr")) %>% distinct(enhancer_seq,.keep_all=TRUE)
enhancer_seqs<-merge(seqs,enhancers,by=c("rsid","chr")) %>% distinct(enhancer_seq,.keep_all=TRUE)
nonsig_seqs<-merge(seqs,nonsig_rsids,by=c("rsid","chr","seq_start","seq_stop")) %>% distinct(enhancer_seq,.keep_all=TRUE)

make_beds_fastas<-function(df,prefix) { #df must hace chr,pos,pos1,rsid as columns
  #Write beds: first write bed with enhancer seq location, then do bedtools merge to merge overlapping regions
  write.table(df[,c('chr','seq_start','seq_stop','rsid')],paste("motif/",prefix,"_unmerged.bed",sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')
  system(paste("sort -k1,1 -k2,2n motif/",prefix,"_unmerged.bed > motif/",prefix,"_sorted.bed",sep=''))
  system(paste("bedtools merge -i motif/",prefix,"_sorted.bed > motif/",prefix,".bed",sep=''))
  system(paste("bedtools getfasta -fi ~/midway/genos/hg38/GRCh38.primary_assembly.genome.fa -bed motif/",prefix,".bed -tab -fo motif/",prefix,".temp",sep=''))

  temp<-fread(paste("motif/",prefix,".temp",sep=''),header=FALSE)
  rev_temp<-chartr("ATGCN","TACGN",sapply(lapply(strsplit(temp, NULL), rev),paste, collapse="")) #Get reverse complement for reverse motif search
  fa<-rep(">",nrow(temp)*2)
  rev_fa<-rep(">",nrow(rev_temp)*2)
  fa[seq(2,length(fa),2)]<-temp$V2
  rev_fa[seq(2,length(rev_fa),2)]<-rev_temp$V2
  write(fa,paste("motif/",prefix,".fa",sep=''))
  write(rev_fa,paste("motif/",prefix,"_rev.fa",sep=''))
  system(paste("rm motif/",prefix,"_unmerged.bed motif/",prefix,"_sorted.bed motif/",prefix,".temp",sep=''))
}

make_beds_fastas(EMVar_seqs,"EMVars")
make_beds_fastas(enhancer_seqs,"enhancers")
make_beds_fastas(nonsig_seqs,"nonsig")

#Run HOMER

#sig enhancers vs genome:
cmd="findMotifsGenome.pl motif/enhancers.bed hg38 motif/enhancers_vs_hg38 -size given &> motif/enhancers_vs_hg38_HOMER.log"
system(cmd)

#sig enhancers vs nonsig enhancers:
cmd="findMotifs.pl motif/enhancers.fa fasta motif/enhancers_vs_nonsig -fasta motif/nonsig.fa &> motif/enhancers_vs_nonsig_HOMER.log"
system(cmd)
cmd="findMotifs.pl motif/enhancers_rev.fa fasta motif/enhancers_rev_vs_nonsig -fasta motif/nonsig_rev.fa &> motif/enhancers_rev_vs_nonsig_HOMER.log"
system(cmd)

#EMVars vs nonsig enhancers
cmd="findMotifs.pl motif/EMVars.fa fasta motif/EMVars_vs_nonsig -fasta motif/nonsig.fa &> motif/EMVars_vs_nonsig_HOMER.log"
system(cmd)
cmd="findMotifs.pl motif/EMVars_rev.fa fasta motif/EMVars_rev_vs_nonsig -fasta motif/nonsig_rev.fa &> motif/EMVars_rev_vs_nonsig_HOMER.log"
system(cmd)