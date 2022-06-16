#~/usr/bin/env Rscript
#Script for assigning dominant TSS's to dREG and ChIP-seq peaks
#Andrew Field - Last edited 3/31/2022
#Update: fixed a typo :)
#To run, use the following command in the desired output directory:
#Rscript --vanilla TSScenR.R -i input.filtered.bed -o out -t TSSclassify_detail_file.txt

#####LOAD LIBRARIES AND OPTIONS#####
message("Loading R libraries...")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse")
if (!requireNamespace("IRanges", quietly = TRUE))
  install.packages("IRanges")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("vroom", quietly = TRUE))
  install.packages("vroom")

library("optparse")
library("IRanges")
library(dplyr)
library("vroom")
opt <- data.frame(input="")

#Accept and parse arguments from command line
message("Loading arguments...")
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="dREG or ChIP-seq peaks in bed6 or bed8 format"),
  make_option(c("-t", "--TSS"), type="character", default=NULL, help="Detailed TSScall output run through TSSclassify.pl"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file prefix")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

#Check for required arguments
message("Checking for required arguments...")
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Must provide input peak file in bed6 or bed8 format.n", call.=FALSE)
} else if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Must provide output prefix.n", call.=FALSE)
} else if(is.null(opt$TSS)){
  print_help(opt_parser)
  stop("Must provide detailed ouput from TSScall run through TSSclassify.n", call.=FALSE)
}


#####LOAD INPUT FILES#####
message("Loading input files...")
#Read input files
bed <- read.table(opt$input, sep="\t")
if(ncol(bed)==6){
  colnames(bed) <- c("chr", "start", "end", "peakID", "score", "strand")
  bed$centerStart <- "NoCenter"
  bed$centerEnd <- "NoCenter"
} else if(ncol(bed)==8){
  colnames(bed) <- c("chr", "start", "end", "peakID", "score", "strand", "centerStart", "centerEnd")
} else{
  print_help(opt_parser)
  stop("Must provide input peak file in bed6 or bed8 format.n", call.=FALSE)
}
bed <- bed[order(bed$chr,bed$start),]

#Vroom is used to index large files prior to loading. as.data.frame() fnc converts from vroom output (tibble) to a dataframe.
tss <- as.data.frame(vroom(opt$TSS, delim="\t", col_select = c(1:8)))

tss_bed <- data.frame("chr"=tss$Chromosome, "start"=tss$Position, "end"=tss$Position, "TSS_ID"=tss$`TSS ID`, "reads"=tss$Reads, "strand"=tss$Strand)
tss_bed <- tss_bed[order(tss_bed$chr,tss_bed$start),]

#####FUNCTION DEFINITIONS#####
#Bedtools function, adapted from A. Akalin http://zvfak.blogspot.com/2011/02/calling-bedtools-from-r.html
bedTools <- function(functionstring, bed1, bed2, opt.string)
{
  #create temp files (to be deleted later)
  a.file=tempfile()
  b.file=tempfile()
  out=tempfile()
  options(scipen = 99) #does not output as scientific notation
  
  #writ bed formatted dataframes to tempfiles
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  #create the command string and call the command using system()
  if(functionstring=="merge")
    command=paste("bedtools",functionstring,opt.string,"-i",a.file,">",out,sep=" ")
  else
    command=paste("bedtools",functionstring,opt.string,"-a",a.file,"-b",b.file,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  #read results table, save as res, delete temp files
  res=read.table(out,header=F)
  unlink(a.file); unlink(b.file); unlink(out)
  return(res)
}


#####DETERMINE OVERLAPPING TSSs#####
overlap <- bedTools("intersect", bed, tss_bed, "-wao")
colnames(overlap) <- c("peak_chr","peak_start","peak_end","peak_ID","peak_score","peak_strand","peak_centerStart","peak_centerEnd","TSS_chr","TSS_start","TSS_end","TSS_ID","TSS_reads","TSS_strand","overlap")
overlap <- filter(overlap, TSS_ID!=".")

#####TSS WINNERS#####
message("Picking TSS winners within each peak...")
peaklist <- bed$peakID
peakout <- data.frame(peak_ID=character(),
                      TSS_count=double(),
                      TSS_IDs=character(),
                      TSS_winner_ID=character(),
                      TSS_winner_position=double(),
                      TSS_winner_strand=character(),
                      TSS_winner_reads=double(),
                      TSS_plus_ID=character(),
                      TSS_plus_position=double(),
                      TSS_plus_reads=double(),
                      TSS_minus_ID=character(),
                      TSS_minus_position=double(),
                      TSS_minus_reads=double(),
                      stringsAsFactors=FALSE)

i=1
while(i<=length(peaklist))
{
  if(is.integer(i/2000))
  {
    message(".")
  }
  peakfilt <- filter(overlap, peak_ID==peaklist[i])
  maxplus=NA
  maxplus_pos=NA
  maxplus_reads=NA
  maxminus=NA
  maxminus_pos=NA
  maxminus_reads=NA
  if(nrow(peakfilt)>=1)
  {
    peakline <- data.frame(peak_ID=peaklist[i])
    peakline$TSS_count <- nrow(peakfilt)
    peakline$TSS_IDs <- paste0(peakfilt$TSS_ID,collapse=",")
    peakline$TSS_winner_ID <- peakfilt[which.max(peakfilt$TSS_reads),]$TSS_ID
    peakline$TSS_winner_position <- peakfilt[which.max(peakfilt$TSS_reads),]$TSS_start
    peakline$TSS_winner_strand <- peakfilt[which.max(peakfilt$TSS_reads),]$TSS_strand
    peakline$TSS_winner_reads <- peakfilt[which.max(peakfilt$TSS_reads),]$TSS_reads
    
    
    peakfiltplus <- filter(peakfilt, TSS_strand=="+")
    if(nrow(peakfiltplus)>=1) {
      maxplus=peakfiltplus[which.max(peakfiltplus$TSS_reads),]$TSS_ID
      maxplus_pos=peakfiltplus[which.max(peakfiltplus$TSS_reads),]$TSS_start
      maxplus_reads=peakfiltplus[which.max(peakfiltplus$TSS_reads),]$TSS_reads
    }
    peakline$TSS_plus_ID <- maxplus
    peakline$TSS_plus_position <- maxplus_pos
    peakline$TSS_plus_reads <- maxplus_reads
    
    peakfiltminus <- filter(peakfilt, TSS_strand=="-") 
    if(nrow(peakfiltminus)>=1) {
      maxminus=peakfiltminus[which.max(peakfiltminus$TSS_reads),]$TSS_ID
      maxminus_pos=peakfiltminus[which.max(peakfiltminus$TSS_reads),]$TSS_start
      maxminus_reads=peakfiltminus[which.max(peakfiltminus$TSS_reads),]$TSS_reads
    }
    peakline$TSS_minus_ID <- maxminus
    peakline$TSS_minus_position <- maxminus_pos
    peakline$TSS_minus_reads <- maxminus_reads
    
  } else {
    peakline <- data.frame(peak_ID=peaklist[i])
    peakline$TSS_count <- 0
    peakline$TSS_IDs <- NA
    peakline$TSS_winner_ID <- NA
    peakline$TSS_winner_position <- NA
    peakline$TSS_winner_strand <- NA
    peakline$TSS_winner_reads <- NA
    peakline$TSS_plus_ID <- NA
    peakline$TSS_plus_position <- NA
    peakline$TSS_plus_reads <- NA
    peakline$TSS_minus_ID <- NA
    peakline$TSS_minus_position <-NA
    peakline$TSS_minus_reads <-NA
  }
  peakout <- rbind(peakout,peakline)
  i=i+1
}

message("Calculating gap distances...")

peakout$gap_distance <- peakout$TSS_plus_position - peakout$TSS_minus_position
peakout$divergent <- ifelse(peakout$gap_distance>0, TRUE, FALSE)

message("Preparing output...")
outtable <- merge(bed,peakout, by.x="peakID", by.y="peak_ID", all=T)
write.table(outtable, file=paste0(opt$output,"_TSSwinners.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
