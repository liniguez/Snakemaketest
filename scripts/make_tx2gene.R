#!/usr/bin/env Rscript

library(GenomicFeatures)

args <- commandArgs(trailingOnly=TRUE)

txdb<-GenomicFeatures::makeTxDbFromGFF(file=args[1], format="gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

write.table(tx2gene,args[2],sep="\t",quote=F)
