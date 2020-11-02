#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: select_variants_MPRA.R <trait> <tissue> <type> <data_source>\n", call.=FALSE)
}
library(data.table)
library(tidyverse)
library(biomaRt)

trait<-"WHR"
tissue<-"adipose_subcutaneous"
type<-"expression"
data_source<-"GTEx_v8"


get_variants<-function(sex) {
  if (sex=="F") {
    trait=paste(trait,"_F",sep='')
  } else if (sex=="M") {
    trait=paste(trait,"_M",sep='')
  }
  setwd(paste("~/midway/",type,"/",trait,"/",tissue,"/results/",sep=''))
  dir.create("posthoc/MPRA")
  if (type=="expression") {
    genes<-fread(paste(data_source,".all.dat.top",sep=''),data.table=FALSE)  
  } else {
    genes<-fread(paste(data_source,".all.dat.top.gene",sep=''),data.table=FALSE)
    genes$ID<-genes$GENE
  }
  
  #Get core list of variants: best eQTL, colocalizing rsid if present, top variants contributing to model
  rsids<-data.table(matrix(nrow=0,ncol=5),stringsAsFactors = FALSE)
  coloc<-fread("posthoc/coloc/colocalizing_genes_rsids",data.table=FALSE,header=FALSE)
  for (i in 1:nrow(genes)) {
    rsids<-rbind(rsids,t(c(genes$CHR[i],genes$BEST.GWAS.ID[i],".","+",paste(genes$ID[i],"_GWAS_ID",sep=''))))
    rsids<-rbind(rsids,t(c(genes$CHR[i],genes$EQTL.ID[i],".","+",paste(genes$ID[i],"_EQTL_ID",sep=''))))
    if (genes$ID[i] %in% coloc$V1) {
      rsids<-rbind(rsids,t(c(genes$CHR[i],coloc$V2[coloc$V1==genes$ID[i]],".","+",paste(genes$ID[i],"_colocalizing_rsid",sep=''))))
    }
  }
  colnames(rsids)<-c("chr","rsid","score","strand","source")
  rsids<-rsids %>% group_by(rsid) %>% mutate(source=paste(source,collapse=',')) %>% distinct(rsid,.keep_all=TRUE)
  rsids$source<-paste(trait,':',rsids$source,sep='')
  
  #Get locations of rsids
  ensembl<-useMart("ENSEMBL_MART_SNP",dataset="hsapiens_snp")
  locs<-getBM(attributes=c('refsnp_id','chrom_start','chrom_end',"allele"),filters='snp_filter',values=rsids$rsid,mart=ensembl)
  rsids<-left_join(rsids,locs,by=c("rsid"="refsnp_id")) %>% distinct(rsid,.keep_all=TRUE)
  rsids<-rsids[c("chr","chrom_start","chrom_end","rsid","score","strand","source","allele")]
  write.table(rsids,paste("posthoc/MPRA/",trait,"_",tissue,"_rsids_core.bed",sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
  
  
  #Get variants in high LD with core variants
  ld_rsids<-data.table(matrix(nrow=0,ncol=7),stringsAsFactors = FALSE)
  for (i in 1:nrow(rsids)) {
    rsid<-rsids$rsid[i]
    chrom<-rsids$chr[i]
    pos<-rsids$chrom_start[i]
    cmd=paste("tabix -h ~/midway/genos/GTEx_v8/vcf/chr",chrom,"_rsid.vcf.gz chr",chrom,":",pos-200000,"-",pos+200000," > ",rsid,".vcf",sep='')
    system(cmd)
    
    cmd=paste("plink --r2 --ld-window-r2 0.95 --ld-snp ",rsid," --ld-window 10000 --maf 0.01 --vcf ",rsid,".vcf --out ",rsid,sep='')
    system(cmd)
    tryCatch(
      expr= {
        ld<-read.table(paste(rsid,".ld",sep=''),stringsAsFactors = FALSE,header=TRUE)
        if (nrow(ld)>500) {
          print(paste("More than 500 variants in ld with ",rsids$rsid[i],sep=''))
          rsids$source[i]<-"LD_REMOVE"
        } else {
          ld1<-ld[,c("CHR_A","BP_A","SNP_A","R2")]
          ld2<-ld[,c("CHR_B","BP_B","SNP_B","R2")]
          colnames(ld1)<-colnames(ld2)<-c("chr","start","rsid","ld")
          ld<-rbind(ld1,ld2)
          ld<-ld[which(ld$rsid!=rsids$rsid[i]),]
          variants<-ld
          variants$stop<-variants$start+1
          variants$score<-rep('.',nrow(variants))
          variants$strand<-rep('+',nrow(variants))
          variants$source<-paste(rsids$rsid[i],":ld=",variants$ld,sep='')
          variants<-variants[,c("chr","start","stop","rsid","score","strand","source")]
          ld_rsids<-rbind(ld_rsids,variants,use.names=FALSE)
        }
      }, error=function(e) {print(paste("no variants in ld with ",rsids$rsid[i],sep=''))}
    )
    system(paste("rm ",rsid,".*",sep=''))
  }
  colnames(ld_rsids)<-c("chr","chrom_start","chrom_end","rsid","score","strand","source")
  ld_rsids<-ld_rsids %>% group_by(rsid) %>% mutate(source=paste(source,collapse=',')) %>% distinct(rsid,.keep_all=TRUE)
  ld_rsids$source<-paste(trait,':',ld_rsids$source,sep='')
  alleles<-getBM(attributes=c('refsnp_id',"allele"),filters='snp_filter',values=ld_rsids$rsid,mart=ensembl)
  ld_rsids<-left_join(ld_rsids,alleles,by=c("rsid"="refsnp_id")) %>% distinct(rsid,.keep_all=TRUE)
  write.table(ld_rsids,paste("posthoc/MPRA/",trait,"_",tissue,"_rsids_ld.bed",sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
  
  #Combine core rsids and rsids in ld
  rsids$chr<-as.integer(rsids$chr)
  rsids_all<-rbind(rsids,ld_rsids)
  rsids_all<-rsids_all %>% group_by(rsid) %>% mutate(source=paste(source,collapse=',')) %>% distinct(rsid,.keep_all=TRUE)
  rsids_all<-rsids_all %>% arrange(chr,chrom_start)
  return(rsids_all)
}

############# Female ############
#Get significant genes list
setwd(paste("~/midway/",type,"/",trait,"_F/",tissue,"/results/",sep=''))
rsids_F<-get_variants("F")
setwd(paste("~/midway/",type,"/",trait,"_M/",tissue,"/results/",sep=''))
rsids_M<-get_variants("M")

rsids_all<-rbind(rsids_F,rsids_M)
rsids_all<-rsids_all %>% group_by(rsid) %>% mutate(source=paste(source,collapse=',')) %>% distinct(rsid,.keep_all=TRUE)
rsids_all<-rsids_all[which(rsids_all$source!="LD_REMOVE"),]
rsids_all<-rsids_all %>% arrange(chr,chrom_start)
write.table(rsids_all,paste("~/projects/MPRA/WHR/variants/",trait,"_",tissue,"_rsids.bed",sep=''),col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
