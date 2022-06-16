#~/usr/bin/env Rscript
#Script for converting dREG output to standard BED file format and filtering peaks by score and read counts
#Andrew Field - Last edited 3/25/2022
#Now accepts bedGraphs for dREG peak read counts
#Added option for requiring a number of samples to pass the read threshold
#To run, use the following command in the desired output directory:
#Rscript --vanilla dRIP_filt.R -i input.dREG.peak.full.bed -o out -s 0.5 -p 0.025 -b bedGraph_manifest.csv

if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("matrixStats", quietly = TRUE))
  install.packages("matrixStats")
if (!requireNamespace("hexbin", quietly = TRUE))
  install.packages("hexbin")

library("optparse")
library(dplyr)
library("ggplot2")
library(matrixStats)
library(hexbin)

#Debugging options
#opt <- data.frame(
#  "input"="../../../Data Analysis/220316_EpiLC_PROseq/dREG/EpiLC.dreg.dREG.peak.full.bed", 
#  "output"="EpiLC", 
#  "score"=0.5, 
#  "pval"=0.025, 
#  "bgs"="bglist_test.csv",
#  "count"=10,
#  "nsamp"=2)

#Accept and parse arguments from command line
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Decompressed dREG output bed file with extension *.dREG.peak.full.bed", metavar="input.dREG.peak.full.bed"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file prefix", metavar="out"),
  make_option(c("-s", "--score"), type="numeric", default=0.50, help="Minimum dREG peak score [default=0.50]", metavar="0.50"),
  make_option(c("-p", "--pval"), type="numeric", default=0.025, help="Maximum dREG adjusted p-value [default=0.025]", metavar="0.025"),
  make_option(c("-b", "--bgs"), type="character", default=NULL, help="Comma-separated table of 3pr mapped bedGraphs and normalization factors (see GitHub for details)", metavar="bg_manifest.csv"),
  make_option(c("-c", "--count"), type="numeric", default=10, help="Minimum plus and minus strand read cutoff [default=10]", metavar="10"),
  make_option(c("-n", "--nsamp"), type="numeric", default=2, help="Minimum number of samples to meet read threshold [default=1]", metavar="1")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Check provided arguments
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Must provide input file, output prefix, and bedGraph manifest.n", call.=FALSE)
} else if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Must provide input file, output prefix, and bedGraph manifest.n", call.=FALSE)
} else if(is.null(opt$bgs)){
  print_help(opt_parser)
  stop("Must provide input file, output prefix, and bedGraph manifest.n", call.=FALSE)
}

#Load input bed file and variables
dREG.df = read.table(opt$input, header=FALSE)
colnames(dREG.df) <- c("Chr", "Start", "End", "Score", "adjPval", "CenterBin")
dREG.df["CenterStart"] <- dREG.df$CenterBin-10

#Adjust dREG peak coordinates so that centers fall within peak boundaries (Center bins sometimes fall on the edge of peaks)
dREG.df <- dREG.df %>%
  rowwise() %>%
  mutate(newStart = min(c(Start,CenterStart)))

dREG.df <- dREG.df %>%
  rowwise() %>%
  mutate(newEnd = max(c(End, CenterBin)))

#Generate QC plots
pdfname <- paste0(opt$output,"_scatter.pdf")
pdf(pdfname)
smoothScatter(dREG.df$Score, dREG.df$adjPval, colramp= colorRampPalette(c("white", "blue", "yellow", "red")), nrpoints=1000, xlab="dREG Score", ylab="Adjusted p-value") + abline(v=opt$score,lty="dashed") + abline(h=opt$pval,lty="dashed")
ggplot(dREG.df, aes(x=Score, y=adjPval), xlab="dREG Score", ylab="Adjusted p-value") + geom_hex(bins=200) + scale_fill_continuous(type="viridis", limits= c(1,50), na.value="yellow") + geom_vline(xintercept=opt$score) + geom_hline(yintercept=opt$pval) + theme_bw()
dev.off()

#Filter by score and p-value
dREG.df <- filter(dREG.df, Score >= opt$score)
dREG.df <- filter(dREG.df, adjPval <= opt$pval)

#Assign peakIDs
peakIDnum <- sprintf("%06d", seq(1,nrow(dREG.df), by=1))
peakIDs <- paste0("dREGpeak", peakIDnum)

#Make modified bed file
bed <- data.frame("Chr"=dREG.df$Chr)
bed$Start <- as.integer(format(dREG.df$newStart, scientific = FALSE))
bed$End <- as.integer(format(dREG.df$newEnd, scientific = FALSE))
bed$Name <- peakIDs
bed$Score <- format(round(dREG.df$Score*1000), scientific = FALSE)
bed$Strand <- "+"
bed$thickStart <- as.integer(format(dREG.df$CenterStart, scientific = FALSE))
bed$thickEnd <- as.integer(format(dREG.df$CenterBin, scientific = FALSE))

#Chunk up input bams into a list of strings
bgs <- read.table(opt$bgs, header=F, stringsAsFactors=F, sep=",")
colnames(bgs) <- c("forward","reverse","norm")

#Bedtools function for union bedgraph
unionbedg <- function(bg_list)
{
  #create temp files (to be deleted later)
  out=tempfile()
  options(scipen = 99) #does not output as scientific notation
  #list <- gsub(" ","\\ ",bg_list,fixed = TRUE) #allows for spaces in directory names
  #list <- gsub("(","\\(",bg_list,fixed = TRUE) #allows for "(" in directory names
  #list <- gsub(")","\\)",bg_list,fixed = TRUE) #allows for ")" in directory names
  
  #create the command string and call the command using system()
  command=paste("bedtools unionbedg -i",paste(bg_list, collapse = " "),">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  #read results table, save as res, delete temp files
  res=read.table(out,header=F)
  
  unlink(out)
  return(res)
}

#Bedtools function for map (read counting)
makeitcount <- function(ubg,b)
{
  #create temp files (to be deleted later)
  a.file=tempfile(fileext=".bed")
  b.file=tempfile(fileext=".bedgraph")
  out=tempfile()
  options(scipen = 99) #does not output as scientific notation
  
  #write bed formatted dataframes to tempfiles
  write.table(b,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(ubg,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  list.of.cols <- paste0(4:ncol(ubg),collapse=",")
  
  #create the command string and call the command using system()
  command=paste("bedtools map -a",a.file,"-b",b.file,"-c",list.of.cols,"-o sum >",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  #read results table, save as res, delete temp files
  res=read.table(out,header=F,stringsAsFactors=F)
  colnames(res)[1:8] <- c("Chr","Start","End","peakID","peakScore","peakStrand","thickStart","thickEnd")

  res[res=="."] <- 0 #replace blank chars in rows with no counts
  res[,9:ncol(res)] <- apply(res[,9:ncol(res)], 2, as.numeric) #coerce counts to integers
  
  unlink(out)
  unlink(a.file)
  unlink(b.file)
  unlink(list.of.cols)
  return(res)
}

#Make union bedgraphs
ubg_f <- unionbedg(bgs$forward)
ubg_r <- unionbedg(bgs$reverse)

#Generate dREG peak counts
counts_f <- makeitcount(ubg_f,bed)
counts_r <- makeitcount(ubg_r,bed)

#normalize sample counts
for (samp in nrow(bgs)){
  counts_f[,samp+8]<-counts_f[,samp+8]/bgs$norm[samp]
  counts_r[,samp+8]<-counts_r[,samp+8]/bgs$norm[samp]
}

#Count number of samples over read count threshold per peak
bool <- data.frame(row.names=bed$Name)
for (samp in 1:nrow(bgs)){
  boolcount<-c()
  for (peak in 1:nrow(bed)){
    if(counts_f[peak,samp+8] >= opt$count && counts_r[peak,samp+8] >= opt$count){
      boolcount[peak] <- 1
    }
    else{
      boolcount[peak] <- 0
    }
  }
  bool[,samp] <- boolcount
}

bool$sum <- rowSums(bool[,c(1:nrow(bgs))])
bool <- bool %>% filter(sum >= opt$nsamp)

bed_filt <- bed %>% filter(Name %in% rownames(bool))

#Write filtered bed file
filename <- paste0(opt$out,"_s",opt$score,"p",opt$pval,"c",opt$count,"n",opt$nsamp,"_filtered.bed")
write.table(bed_filt, file=filename, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
