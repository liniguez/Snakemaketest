#!/usr/bin/env Rscript

library(tximport)

args <- commandArgs(trailingOnly=TRUE)

samples <- read.table(args[1], header = TRUE,sep="\t",stringsAsFactors = F)
fold <- args[2]
tx2gene <-read.table(args[3])

samples<-subset(samples, outFolder==fold)
samp<-unique(samples$SampleName)

files<-file.path(fold,"Salmon",samp,"quant.sf")
txi.salmon <- tximport(files, type = "salmon",
                  txOut = TRUE,ignoreTxVersion=TRUE,
                  dropInfReps=T)$counts
salmon_t<-round(txi.salmon)
mode(salmon_t)<-"integer"
colnames(salmon_t)<-samp
fileout<-file.path(fold,"Salmon_Transcript_table_counts.csv")
write.table(salmon_t,file=fileout, sep=",",quote=F)

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene,dropInfReps=T,ignoreTxVersion=T)
salmon_g<-round(txi.salmon$counts)
mode(salmon_g)<-"integer"
colnames(salmon_g)<-samp
fileout<-file.path(fold,"Salmon_Gene_table_counts.csv")
write.table(salmon_g,file=fileout, sep=",",quote=F)


salmon_g<-txi.salmon$abundance
colnames(salmon_g)<-samp
fileout<-file.path(fold,"Salmon_Gene_table_TPM.csv")
write.table(salmon_g,file=fileout, sep=",",quote=F)
