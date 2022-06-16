#~/usr/bin/env Rscript
#Script for generating adjusted gene windows for read counts with enhancer annotation data
#Andrew Field - Last edited 3/30/2022
#To run, use the following command in the desired output directory:
#Rscript --vanilla dice.R -g gga.active.bed -i intra.txt -o out -s 250 -e 2250 -c 500
#Name Ideals:
#dice.R = Distal Intragenic Cropping of Elements
#mAUI = mRNA Annotation Useful Indexing
#PAWS = PRO-seq Annotation WindowS

#####LOAD LIBRARIES AND OPTIONS#####
message("Loading R libraries...")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr")

library("optparse")

opt <- data.frame(input="")
#Accept and parse arguments from command line
message("Loading arguments...")
option_list <- list(
  make_option(c("-g", "--GGA"), type="character", default=NULL, help="GGA Active BED annotation output (*_Dominant.TSS.TES.calls.bed)"),
  make_option(c("-i", "--intra"), type="character", default=NULL, help="Intragenic dREG peak output from CoGENT output (*_intragenic.txt)"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output prefix"),
  make_option(c("-s", "--start"), type="numeric", default=250, help="Gene body start coordinates relative to TSS (default: 250bp)."),
  make_option(c("-e", "--end"), type="numeric", default=2250, help="Gene body maximum end coordinates relative to TSS (default: 2250bp)."),
  make_option(c("-t", "--trim"), type="numeric", default=500, help="Trim distance for intragenic enhancers. Will trim gene body annotations by this distance upstream of overlapping dREG peaks (default: 500bp)."),
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

#Check for required arguments
message("Checking for required arguments...")
if(is.null(opt$GGA)){
  print_help(opt_parser)
  stop("Must provide options for active genes bed from GGA, intragenic peak calls from CoGENT, and desired output prefix.n", call.=FALSE)
} else if(is.null(opt$intra)){
  print_help(opt_parser)
  stop("Must provide options for active genes bed from GGA, intragenic peak calls from CoGENT, and desired output prefix.n", call.=FALSE)
} else if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Must provide options for active genes bed from GGA, intragenic peak calls from CoGENT, and desired output prefix.n", call.=FALSE)
}

library(dplyr)
library(tidyr)


#Load data tables
message("Loading annotations...")
gga <- read.table(opt$GGA, sep="\t", stringsAsFactors=F, header=F)
colnames(gga) <- c("gene_chr","gene_start","gene_end","geneID","gene_score","gene_strand")

intra <- read.table(opt$intra, sep="\t", stringsAsFactors=F, header=T)


#Initial trimming of gene annotations
message("Initial trimming of active transcripts...")
gga_trim <- data.frame("geneID"=gga$geneID,"chr"=gga$gene_chr)
gga_trim$start <- ifelse(gga$gene_strand == "+", 
                         gga$gene_start+opt$start, 
                         ifelse(gga$gene_end-gga$gene_start >= opt$end,
                                gga$gene_end-opt$end,
                                gga$gene_start))
gga_trim$end <- ifelse(gga$gene_strand == "+", 
                       ifelse(gga$gene_end-gga$gene_start >= opt$end, 
                              gga$gene_start+opt$end,
                              gga$gene_end),
                       gga$gene_end-opt$start)
gga_trim$strand <- gga$gene_strand

gga_trim <- gga_trim %>% filter(end > start)

message(paste0(nrow(gga)-nrow(gga_trim)," out of ",nrow(gga)," (", format(round(100*(nrow(gga)-nrow(gga_trim))/nrow(gga), 2), nsmall = 2),"%) active transcripts lost due to size restraints (annotation window less than ",opt$start,"bp long)."))


#Defining Bedtools function
message("Defining script functions...")
bedTools <- function(functionstring, bed1, bed2, opt.string)
{
  #create temp files (to be deleted later)
  a.file=tempfile()
  b.file=tempfile()
  out=tempfile()
  options(scipen = 99) #does not output as scientific notation
  
  #write bed formatted dataframes to tempfiles
  write.table(bed1,file=a.file, quote=F, sep="\t", col.names=F, row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  #create the command string and call the command using system()
  command=paste("bedtools",functionstring,opt.string,"-a",a.file,"-b",b.file,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  #read results table, save as res, delete temp files
  res=read.table(out,header=F)
  colnames(res)=c("gene_chr","gene_start","gene_end","geneID","gene_score","gene_strand","peak_chr","peak_start","peak_end","peakID","peak_score","peak_strand","overlap")
  unlink(a.file); unlink(b.file); unlink(out)
  return(res)
}

#Overlap gene windows with intragenic enhancers
message("Overlapping trimmed gene windows with intragenic enhancer annotations...")

overlaps <- bedTools("intersect",
                     gga_trim[,c("chr","start","end","geneID","start","strand")],
                     intra[,c("peak_chr","peak_start","peak_end","peakID","peak_score","strand")],
                     "-wao")

overlaps <- overlaps[overlaps$peakID!=".",] #filter for genes with overlap

message(paste0(length(unique(overlaps$geneID)), " active transcript windows overlap intragenic dREG peaks."))


#Truncating gene coordinates based on intragenic peak overlap
adj_coords <- t(apply(gga_trim,MAR=1,function(x){
  temp <- x[c("start","end")]
  dREG_starts <- overlaps$peak_start[overlaps$geneID==x["geneID"]]
  dREG_ends <- overlaps$peak_end[overlaps$geneID==x["geneID"]]
  if(x["strand"]=="+") {
    temp[2] <- min(c(x["end"],dREG_starts-opt$trim))
  }else{
    temp[1] <- max(c(x["start"],dREG_ends+opt$trim))
  }
  temp
}))

gga_trim[c("start","end")] <- adj_coords

gga_trim2 <- gga_trim %>% filter(end > start)

message(paste0(nrow(gga_trim)-nrow(gga_trim2)," out of ",nrow(gga_trim)," (", format(round(100*(nrow(gga_trim)-nrow(gga_trim2))/nrow(gga_trim), 2), nsmall = 2),"%) active transcripts lost due intragenic peak trimming."))

stat_table <- data.frame("Active BED Entries"=nrow(gga),
                         "Lost to Coordinate Cropping"=nrow(gga)-nrow(gga_trim),
                         "Lost to Intragenic Peak Overlap"=nrow(gga_trim)-nrow(gga_trim2),
                         "Remaining BED Entries"=nrow(gga_trim2),
                         "Percent Remaining BED Entries"=format(round(nrow(gga_trim2)*100/nrow(gga), 1), nsmall = 1))

#Make some output tables
Message("Making output tables...")

formakeheatmap
bed
saf

lost transcript list
- point at which it is lost?

stat table