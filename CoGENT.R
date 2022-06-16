#~/usr/bin/env Rscript
#Script for transcribed locus annotation from PRO-seq data
#Andrew Field - Last edited 3/30/2022
#Fixed some duplicated lines in outputs, fixed transcript filtering issues
#Updated annotation output to include strand information for intragenic and promoter proximal peaks
#To run, use the following command in the desired output directory:
#Rscript --vanilla CoGENT.R -i input.filtered.bed -o out -c basic.gtf -a annotations_Dominant.TSS.TES.calls.bed -T Annotated_Dominant_and_Nondominant_obsTSS_fordREG.txt -d [TSS/GB] -p 1000 -m 1000000 -g 1000 -r rmsk.txt

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
if (!requireNamespace("tidyr", quietly = TRUE))
  install.packages("tidyr")
if (!requireNamespace("vroom", quietly = TRUE))
  install.packages("vroom")

library("optparse")

opt <- data.frame(input="")
#Accept and parse arguments from command line
message("Loading arguments...")
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Filtered dREG bed file (filtered by both dREGfilter.R and by max read counts [suggested: >30 max reads])"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file prefix"),
  make_option(c("-c", "--comp_annotation"), type="character", default=NULL, help="Basic annotation gtf (unfiltered, UCSC coordinate compatible [eg. Gencode v33 https://www.gencodegenes.org/human/release_33.html])"),
  make_option(c("-a", "--active"), type="character", default=NULL, help="Filtered transcript annotation bed from geneAnnotation pipeline (*_Dominant.TSS.TES.calls.bed)"),
  make_option(c("-T", "--TSScall"), type="character", default=NULL, help="Observed transcription start sites in PRO-seq from GGA or detailed output file from TSScall/TSSclarify (*_Annotated_Dominant_and_Nondominant_obsTSS_fordREG.txt)"),
  make_option(c("-d", "--distance"), type="character", default="TSS", help="Coordinate to which distance is measured for associated genes. Options: TSS or GB (gene body) [default=TSS]"),
  make_option(c("-p", "--promoter"), type="numeric", default=1000, help="Distance in bp from TSS for promoter proximal regions [default=1000]"),
  make_option(c("-m", "--max"), type="numeric", default=1000000, help="Maximum distance for gene association [default=1000000]"),
  make_option(c("-g", "--gap"), type="numeric", default=1000, help="Minimum gap distance between dREG peaks for separate enhancer clusters. Peaks closer together than this value (in bp) will be merged into enhancer clusters [default=1000]"),
  make_option(c("-r", "--rmsk"), type="character", default=NULL, help="Unzipped RepeatMasker track from UCSC Genome Browser (eg http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz)"),
  make_option(c("-R", "--rfilt"), type="character", default="rRNA,snRNA,scRNA,srpRNA,tRNA", help="Filtered repeats as a comma separated list [default=rRNA,snRNA,scRNA,srpRNA,tRNA]"),
  make_option(c("-f", "--rnas"), type="character", default="snoRNA,Mt_rRNA,miRNA,scaRNA,scRNA,snRNA,Mt_tRNA,ribozyme", help="Filtered functional RNA biotypes as a comma separated list [default=rRNA,snRNA,scRNA,srpRNA,tRNA]"),
  make_option(c("-s", "--supp"), type="character", default="protein_coding,snRNA", help="Biotypes suppressed from transcript support level (TSL) filter entered as comma separated list [default=protein_coding,snRNA]"),
  make_option(c("-t", "--TSL"), type="numeric", default=3, help="Minimum transcript support level (TSL) for non-coding transcript filter [default=3]")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

#Check for required arguments
message("Checking for required arguments...")
if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Must provide filtered dREG bed, output prefix, comprehensive gtf, active genes bed, obsTSS calls, and RepeatMasker track.n", call.=FALSE)
} else if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Must provide filtered dREG bed, output prefix, comprehensive gtf, active genes bed, obsTSS calls, and RepeatMasker track.n", call.=FALSE)
} else if(is.null(opt$comp_annotation)){
  print_help(opt_parser)
  stop("Must provide filtered dREG bed, output prefix, comprehensive gtf, active genes bed, obsTSS calls, and RepeatMasker track.n", call.=FALSE)
} else if(is.null(opt$active)){
  print_help(opt_parser)
  stop("Must provide filtered dREG bed, output prefix, comprehensive gtf, active genes bed, obsTSS calls, and RepeatMasker track.n", call.=FALSE)
} else if(is.null(opt$TSScall)){
  print_help(opt_parser)
  stop("Must provide filtered dREG bed, output prefix, comprehensive gtf, active genes bed, obsTSS calls, and RepeatMasker track.n", call.=FALSE)  
} else if(is.null(opt$rmsk)){
  print_help(opt_parser)
  stop("Must provide filtered dREG bed, output prefix, comprehensive gtf, active genes bed, obsTSS calls, and RepeatMasker track.n", call.=FALSE) 
}

library("IRanges")
library(dplyr)
library(tidyr)
library("vroom")

#####LOAD ANNOTATION FILES#####
message("Loading annotation files...")
#Read annotation files
active_bed <- read.table(opt$active, sep="\t")
colnames(active_bed) <- c("chr", "start", "end", "geneID", "score", "strand")
conflict_bed <- active_bed[(active_bed$start>=active_bed$end),] #Add ambiguous category where TSS is after TES
write.table(conflict_bed, paste0(opt$output,"_conflicted_annotation_entries.bed"), row.names = FALSE)

rmsk <- read.table(opt$rmsk, sep="\t")
colnames(rmsk) <- c("bin", "swScore", "milliDiv", "milliDel", "milliIns", "genoName",	"genoStart", "genoEnd", "genoLeft",	"strand",	"repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft", "id")
rmsk_bed <- rmsk[,c(6:8)]
rmsk_bed$name <- paste(rmsk$repName,rmsk$repClass,rmsk$repFamily,sep=":")
rmsk_bed$score <- rmsk$swScore
rmsk_bed$strand <- rmsk$strand

#Remove conflict entries
active_bed <- active_bed[!(active_bed$geneID %in% conflict_bed$geneID),]
active_bed <- active_bed[order(active_bed$chr,active_bed$start),]

#Pull comprehensive annotation list from original gtf (UCSC Genome compatible chromosome names)
comp_gtf <- read.table(opt$comp_annotation, sep="\t") #comprehensive annotation list (unfiltered -- used to create "comp" dataframes below)
colnames(comp_gtf) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
comp_tx <- filter(comp_gtf, type == "transcript")
comp_gene <- filter(comp_gtf, type == "gene")

#Load TSScall results from GGA for comprehensive obsTSS list, switched to loading with vroom to accomodate larger TSScall output
#TSScall <- read.table(opt$TSScall, header=T, sep="\t")
TSScall <- as.data.frame(vroom(opt$TSScall, delim="\t", show_col_types=F))
TSScall <- unique(separate_rows(TSScall, `Gene.ID`, sep=";"))

#####FUNCTION DEFINITIONS#####
message("Defining R functions...")
#GTF attribute extraction function
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(as.character(gtf_attributes), "; ")
  att <- gsub("\"","",unlist(att))
  att <- gsub(";","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

#Bedtools function
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
  if(functionstring=="merge")
    command=paste("bedtools",functionstring,opt.string,"-i",a.file,">",out,sep=" ")
  else
    command=paste("bedtools",functionstring,opt.string,"-a",a.file,"-b",b.file,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  #read results table, save as res, delete temp files
  res=read.table(out,header=F)
  if(functionstring=="intersect") #assign output column names
    colnames(res)=c("peak_chr","peak_start","peak_end","peakID","peak_score","peak_strand","peak_center_start","peak_center_end","gene_chr","gene_start","gene_end","geneID","gene_score","gene_strand","overlap")
  if(functionstring=="closest")
    colnames(res)=c("peak_chr","peak_start","peak_end","peakID","peak_score","peak_strand","peak_center_start","peak_center_end","gene_chr","gene_start","gene_end","geneID","gene_score","gene_strand","distance")
  #  if(functionstring=="merge")
  #    colnames(res)=c("cluster_chr","cluster_start","cluster_end","peak_count","peak_IDs")
  unlink(a.file); unlink(b.file); unlink(out)
  return(res)
}


#####GENERATE GENE MODELS AND ANNOTATION LOOK-UP TABLES#####
message("Generating gene models and annotation look-up tables...")

#Extract desired attributes for transcript and gene look-up table
comp_tx$geneID <- unlist(lapply(comp_tx$attributes, extract_attributes, "gene_id"))
comp_tx$geneID <- sapply(strsplit(comp_tx$geneID, split="[.]"), head, 1) #remove gene version number if present
comp_tx$txID <- unlist(lapply(comp_tx$attributes, extract_attributes, "transcript_id"))
comp_tx$geneName <- unlist(lapply(comp_tx$attributes, extract_attributes, "gene_name"))
comp_tx$biotype <- unlist(lapply(comp_tx$attributes, extract_attributes, "gene_biotype"))
comp_tx$TSL <- unlist(lapply(comp_tx$attributes, extract_attributes, "transcript_support_level"))

#Removed suppressed biotypes and TSLs from comp_tx list
keep_biotypes <- unlist(strsplit(opt$supp,",")) #parse kept biotypes
yeet_tx <- comp_tx %>% filter((!biotype%in%keep_biotypes & TSL>opt$TSL) | (!biotype%in%keep_biotypes & TSL=="NA") | (!biotype%in%keep_biotypes & is.na(TSL))) #define suppressed transcripts
keep_tx <- comp_tx %>% filter(!txID%in%yeet_tx$txID) #filter out suppressed transcripts

#comp_gene coords defined by kept transcripts (TSS to TES)
#supp_gene coords defined by suppressed transcripts (-1kb to TES)
keep_coords <- keep_tx %>%
  group_by(geneID) %>%
  summarize(start = min(start), end = max(end))

yeet_coords <- yeet_tx %>%
  group_by(geneID) %>%
  summarize(start = min(start), end = max(end))

comp_gene$geneID <- unlist(lapply(comp_gene$attributes, extract_attributes, "gene_id"))
comp_gene$geneID <- sapply(strsplit(comp_gene$geneID, split="[.]"), head,1) #remove gene version number if present
comp_gene$geneName <- unlist(lapply(comp_gene$attributes, extract_attributes, "gene_name"))
comp_gene$biotype <- unlist(lapply(comp_gene$attributes, extract_attributes, "gene_biotype"))

#Replace start and end coordinates with min and max coordinates for kept and suppressed transcripts
keep_gene <- comp_gene %>% 
  filter(geneID %in% keep_coords$geneID)

keep_gene[,c('start','end')] <- keep_coords[match(keep_gene$geneID,keep_coords$geneID),c('start','end')]

yeet_gene <- comp_gene %>% 
  filter(geneID %in% yeet_coords$geneID)

yeet_gene[,c('start','end')] <- yeet_coords[match(yeet_gene$geneID,yeet_coords$geneID),c('start','end')]

#Set up inactive gene bed catch (from keep_gene)
keep_gb <- keep_gene[,c("chr","start","end","geneID","score","strand")]

#Grab observed TSS's
obs_tss <- data.frame("chr"=TSScall$Chromosome)
obs_tss$chromStart <- ifelse(TSScall$Strand=="+",TSScall$Position-opt$promoter-1,TSScall$Position-opt$promoter) #Adjust TSScall coordinates to 0-based from 1-based
obs_tss$chromEnd <- ifelse(TSScall$Strand=="+",TSScall$Position+opt$promoter-1,TSScall$Position+opt$promoter)
obs_tss$geneID <- TSScall$`Gene.ID`
obs_tss$score <- TSScall$Reads
obs_tss$strand <- TSScall$Strand
obs_tss <- obs_tss[order(obs_tss$chr,obs_tss$chromStart),]

#Define TSS regions (+/- promoter distance)
active_tss <- data.frame("chr"=active_bed$chr)
active_tss$chromStart <- ifelse(active_bed$strand=="+", active_bed$start-opt$promoter, active_bed$end-opt$promoter)
#active_tss$chromStart <- ifelse(active_bed$chromStart>=0, active_bed$chromStart, 0)
active_tss$chromEnd <- ifelse(active_bed$strand=="+", active_bed$start+opt$promoter, active_bed$end+opt$promoter)
active_tss$geneID <- active_bed$geneID
active_tss <- merge(active_tss[,c(1:4)], obs_tss, all.x=T) #grab read scores from obs_tss object
active_tss$strand <- active_bed$strand
active_tss <- active_tss[order(active_tss$chr,active_tss$chromStart),]

unstable_tss <- unique(obs_tss[!(obs_tss$geneID %in% active_tss$geneID),])

nondom_tss <- obs_tss[(obs_tss$geneID %in% active_tss$geneID),]
nondom_tss <- anti_join(nondom_tss, active_tss)

#Add in kept and yeeted TSS's
keep_tss <- data.frame("chr"=keep_tx$chr)
keep_tss$chromStart <- ifelse(keep_tx$strand=="+", keep_tx$start-opt$promoter, keep_tx$end-opt$promoter)
keep_tss$chromStart <- ifelse(keep_tss$chromStart>=0, keep_tss$chromStart, 0)
keep_tss$chromEnd <- ifelse(keep_tx$strand=="+", keep_tx$start+opt$promoter, keep_tx$end+opt$promoter)
keep_tss$geneID <- keep_tx$geneID
keep_tss$score <- 0
keep_tss$strand <- keep_tx$strand
keep_tss <- keep_tss[order(keep_tss$chr,keep_tss$chromStart),]

yeet_tss <- data.frame("chr"=yeet_tx$chr)
yeet_tss$chromStart <- ifelse(yeet_tx$strand=="+", yeet_tx$start-opt$promoter, yeet_tx$end-opt$promoter)
yeet_tss$chromStart <- ifelse(yeet_tss$chromStart>=0, yeet_tss$chromStart, 0)
yeet_tss$chromEnd <- ifelse(yeet_tx$strand=="+", yeet_tx$start+opt$promoter, yeet_tx$end+opt$promoter)
yeet_tss$geneID <- yeet_tx$geneID
yeet_tss$score <- NA
yeet_tss$strand <- yeet_tx$strand
yeet_tss <- yeet_tss[order(yeet_tss$chr,yeet_tss$chromStart),]

#Define just TSS positions
active_tssP <- data.frame("chr"=active_bed$chr)
active_tssP$chromStart <- ifelse(active_bed$strand=="+", active_bed$start, active_bed$end)
active_tssP$chromEnd <- active_tssP$chromStart
active_tssP$geneID <- active_bed$geneID
active_tssP$score <- active_bed$score
active_tssP$strand <- active_bed$strand
active_tssP <- active_tssP[order(active_tssP$chr,active_tssP$chromStart),]

#Create lookup table for obsTSS positions. Genes without an obsTSS call will be listed with a position value of NA.
obs_tssP <- data.frame("geneID"=TSScall$`Gene.ID`, 
                       "obsTSS"=ifelse(TSScall$Strand=="+",TSScall$Position-1,TSScall$Position), #Adjust TSScall coordinates to 0-based from 1-based
                       "reads.TSS"=TSScall$Reads)
obs_tssP <- unique(obs_tssP)
no_obsTSS <- data.frame("geneID"=comp_gene[!(comp_gene$geneID %in% obs_tssP$geneID),]$geneID,
                        "obsTSS"=NA,
                        "reads.TSS"=0)
no_obsTSS <- unique(no_obsTSS)
all_tssP <- rbind(obs_tssP, no_obsTSS)

#Pull in filtered dREG results
dREG <- read.table(opt$input, sep="\t")
colnames(dREG) <- c("chr", "start", "end", "peakID", "score", "strand", "centerStart", "centerEnd")


#####REMOVING PROBLEMATIC LOCI#####
#Filter dREG peaks overlapping suppressed repeat types
message("Removing suppressed RepeatMasker elements...")
banned_rmsk <- rmsk %>% filter(repClass %in% unlist(strsplit(opt$rfilt,",")))
banned_rmsk_bed <- banned_rmsk[,c(6:8)]
banned_rmsk_bed$name <- paste(banned_rmsk$repName,banned_rmsk$repClass,banned_rmsk$repFamily,sep=":")
banned_rmsk_bed$score <- banned_rmsk$swScore
banned_rmsk_bed$strand <- banned_rmsk$strand

banned_overlap <- bedTools("intersect", dREG, banned_rmsk_bed, "-wao")
catch <- banned_overlap[banned_overlap$geneID==".",]
banned_overlap <- banned_overlap[banned_overlap$geneID!=".",]

catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)] #Drop null columns

#Removing remaining function RNA elements from consideration as putative enhancers
message("Removing functional RNA elements...")
banned_rnas <- comp_tx %>% filter(biotype %in% unlist(strsplit(opt$rnas,",")))
banned_rnas_bed <- banned_rnas[,c(1,4,5)]
banned_rnas_bed$name <- paste(banned_rnas$geneID,banned_rnas$txID,banned_rnas$geneName,banned_rnas$biotype,sep=":")
banned_rnas_bed$score <- banned_rnas$score
banned_rnas_bed$strand <- banned_rnas$strand

banned_overlap2 <- bedTools("intersect", catch, banned_rnas_bed, "-wao")
catch <- banned_overlap2[banned_overlap2$geneID==".",]
banned_overlap_merged <- unique(rbind(banned_overlap, banned_overlap2[banned_overlap2$geneID!=".",]))

#Report filtered peaks
colnames(banned_overlap_merged)[9:14] <- c("element_chr","element_start","element_end","element","element_score","element_strand")
write.table(banned_overlap_merged, file=paste0(opt$output,"_RNA_rmsk_filtered_peaks.txt"), sep="\t", quote=F, row.names=F, col.names=T)

catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)] #Drop null columns


#####PULL dREG PEAKS OVERLAPPING PROMOTERS#####
message("Pulling promoter proximal dREG peaks...")
#Find promoter proximal and distal peaks (category = active gene promoter)
active_promprox <- bedTools("intersect", catch, active_tss, "-wao")

#Catch unstable gene promoters (category = non-dominant gene promoter)
catch <- active_promprox[active_promprox$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)] #Drop null columns

nondom_promprox <- bedTools("intersect", catch, nondom_tss, "-wao")
active_promprox <- active_promprox[active_promprox$geneID!=".",]

#Catch unstable gene promoters (category = unstable gene promoter)
catch <- nondom_promprox[nondom_promprox$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)] #Drop null columns

unstable_promprox <- bedTools("intersect", catch, unstable_tss, "-wao")
nondom_promprox <- nondom_promprox[nondom_promprox$geneID!=".",]

#Catch for the comp_tss list (category = inactive gene promoter)
catch <- unstable_promprox[unstable_promprox$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)] #Drop null columns

#Inactive TSS's (not suppressed)
inactive_promprox <- bedTools("intersect", catch, keep_tss, "-wao")
unstable_promprox <- unstable_promprox[unstable_promprox$geneID!=".",]

catch <- inactive_promprox[inactive_promprox$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)] #Drop null columns

#Suppressed transcript TSS's
suppress_promprox <- bedTools("intersect", catch, yeet_tss, "-wao")
suppress_promprox <- suppress_promprox[suppress_promprox$geneID!=".",]

#Make bed for distal peaks (inter-/intragenic) -- Note: INCLUDES peaks overlapping suppressed transcript promoters
dist.bed <- unique(inactive_promprox[inactive_promprox$geneID==".",])
dist.bed <- dist.bed[,!(names(dist.bed) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
dist.bed <- dist.bed[-c(9:15)] #Drop null columns

#Finalize inactive TSS list
inactive_promprox <- inactive_promprox[inactive_promprox$geneID!=".",]

#Compile active promoter proximal peak list and add to annotation lookup table
annot_list <- unique(active_promprox[c(1:5,14,7:8)]) #Start master dREG peak annotation list
annot_list$localGeneStatus <- "Active" 
annot_list$TSSproximity <- "Dominant Promoter"
annot_list$genomicCategory <- "Promoter Proximal Element"

#Compile non-dominant promoter proximal peak list and add to annotation lookup table
build_list <- unique(nondom_promprox[c(1:5,14,7:8)])
build_list$localGeneStatus <- "Active" 
build_list$TSSproximity <- "Non-Dominant Promoter"
build_list$genomicCategory <- "Promoter Proximal Element"
annot_list <- rbind(build_list, annot_list)

#Compile unstable gene promoters and add to annotation lookup table
build_list <- unique(unstable_promprox[c(1:5,14,7:8)])
build_list$localGeneStatus <- "Unstable"
build_list$TSSproximity <- "Unstable Gene Promoter"
build_list$genomicCategory <- "Promoter Proximal Element"
annot_list <- rbind(build_list, annot_list)

#Compile inactive gene promoters and add to annotation lookup table
build_list <- unique(inactive_promprox[c(1:5,14,7:8)])
build_list$localGeneStatus <- "Inactive"
build_list$TSSproximity <- "Inactive Gene Promoter"
build_list$genomicCategory <- "Promoter Proximal Element"
annot_list <- rbind(build_list, annot_list)

#Combine promprox dREG peaks into one master list
active_promprox$localGeneStatus <- "active"
nondom_promprox$localGeneStatus <- "active"
unstable_promprox$localGeneStatus <- "unstable"
inactive_promprox$localGeneStatus <- "inactive"

active_promprox$TSSdominance <- "dominant"
nondom_promprox$TSSdominance <- "non-dominant"
unstable_promprox$TSSdominance <- NA
inactive_promprox$TSSdominance <- NA
#Note: Final compilation of promoter proximal list to come after intragenic peak assignment


#####PULL INTERGENIC AND INTRAGENIC PEAKS#####
message("Categorizing distal dREG peaks...")
#Inter/intragenic (distal) peak to gene assignments and completing annotation list
active_intra <- bedTools("intersect", dist.bed, active_bed, "-wao")
catch <- active_intra[active_intra$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)] #Drop null columns
inactive_intra <- bedTools("intersect", catch, keep_gb, "-wao")

active_intra <- active_intra[active_intra$geneID!=".",]

interg <- inactive_intra[inactive_intra$geneID==".",]
interg <- interg[,!(names(interg) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
interg <- interg[-c(9:15)] #Drop null columns

colnames(interg)[6] <- "gene_strand" #Replace with overlapping gene strand place holder
interg$gene_strand <- "." #No strand char assigned

inactive_intra <- inactive_intra[inactive_intra$geneID!=".",]

#Compile suppressed gene promoters and generate annotation table entries
build_list <- unique(suppress_promprox[c(1:5,14,7:8)])
build_list$localGeneStatus <- NA

build_list <- build_list %>% mutate(localGeneStatus = case_when(peakID %in% active_intra$peakID ~ "Active",
                                                                peakID %in% inactive_intra$peakID ~ "Inactive",
                                                                peakID %in% interg$peakID ~ "Suppressed"))

build_list$TSSproximity <- "Suppressed Gene Promoter"

build_list$genomicCategory <- NA
build_list <- build_list %>% mutate(genomicCategory = case_when(peakID %in% active_intra$peakID ~ "Intragenic Transcribed Locus",
                                                                peakID %in% inactive_intra$peakID ~ "Intragenic Transcribed Locus",
                                                                peakID %in% interg$peakID ~ "Intergenic Transcribed Locus"))

annot_list <- unique(rbind(build_list, annot_list))

suppress_promprox$localGeneStatus <- NA
suppress_promprox <- suppress_promprox %>% mutate(localGeneStatus = case_when(peakID %in% active_intra$peakID ~ "Active",
                                                                peakID %in% inactive_intra$peakID ~ "Inactive",
                                                                peakID %in% interg$peakID ~ "Suppressed"))

suppress_promprox$TSSdominance <- NA

#Create master promoter proximal peak table
promprox <- unique(rbind(active_promprox, nondom_promprox, unstable_promprox, inactive_promprox, suppress_promprox))
names(promprox)[names(promprox) == 'gene_score'] <- 'reads@obsTSS' #this gene score column was derived from the reads featuring this 5' in the obsTSS script
names(promprox)[names(promprox) == 'overlap'] <- 'dREGpeak_overlap' #basepair overlap of promoter proximal region with associated dREG peak
names(promprox)[names(promprox) == 'attributes'] <- 'geneAttributes'

#Prepare detailed promoter proximal sheet
promprox <- merge(promprox, comp_gene, by.x="geneID", by.y="geneID") #Add gene info
promprox <- merge(promprox, all_tssP, by.x="geneID", by.y="geneID") #Add obsTSS's
promprox$TSS_in_dREGpeak <- ifelse(promprox$obsTSS>=promprox$peak_start & promprox$obsTSS<=promprox$peak_end, TRUE, FALSE)
promprox <- promprox[,!(names(promprox) %in% c("peak_strand", "gene_chr", "gene_start", "gene_end", "gene_score", "gene_strand", "type", "score", "phase"))] #remove excess/redundant columns
promprox$geneName <- paste0("=\"",promprox$geneName,"\"") #make gene names Excel safe
promprox <- promprox[c("peakID", "peak_chr", "peak_start", "peak_end", "peak_center_start", "peak_center_end", "peak_score", "localGeneStatus", "TSSdominance", "geneID", "geneName", "source", "biotype", "chr", "start", "end", "strand", "obsTSS", "reads@obsTSS", "dREGpeak_overlap", "TSS_in_dREGpeak", "attributes")] #reorder columns
promprox <- promprox[order(promprox$peakID),]
row.names(promprox) <- 1:nrow(promprox)
promprox <- unique(promprox)
write.table(promprox, file=paste0(opt$output,"_p",opt$promoter,"_promprox.txt"), sep="\t", quote=F, row.names=F, col.names=T)

#Add intragenic and intergenic categories to annotation list
build_list <- unique(active_intra[c(1:5,14,7:8)])
build_list <- build_list %>% filter(!peakID %in% suppress_promprox$peakID) #exclude suppressed from this list
build_list$localGeneStatus <- "Active"
build_list$TSSproximity <- "Distal"
build_list$genomicCategory <- "Intragenic Transcribed Locus"
annot_list <- rbind(build_list, annot_list)

build_list <- unique(inactive_intra[c(1:5,14,7:8)]) #exclude suppressed from this list
build_list <- build_list %>% filter(!peakID %in% suppress_promprox$peakID) #exclude suppressed from this list
build_list$localGeneStatus <- "Inactive"
build_list$TSSproximity <- "Distal"
build_list$genomicCategory <- "Intragenic Transcribed Locus"
annot_list <- rbind(build_list, annot_list)

build_list <- interg
build_list <- build_list %>% filter(!peakID %in% suppress_promprox$peakID) #exclude suppressed from this list
build_list$localGeneStatus <- NA
build_list$TSSproximity <- "Distal"
build_list$genomicCategory <- "Intergenic Transcribed Locus"
annot_list <- rbind(build_list, annot_list)

annot_list <- annot_list[order(annot_list$peak_chr,annot_list$peak_start),]
annot_list <- unique(annot_list)
row.names(annot_list) <- 1:nrow(annot_list)

#Overlaps with RepeatMasker unsuppressed elements
message("Overlapping with RepeatMasker elements...")
rmsk_overlap <- bedTools("intersect", annot_list[1:8], rmsk_bed, "-wao")
rmsk_match <- unique(rmsk_overlap[,c(4,12)])
rmsk_list <- t(sapply(sort(unique(rmsk_match$peakID)),function(x){
  subset_table <- rmsk_match[rmsk_match$peakID==x,]
  rname <- paste(subset_table$geneID, collapse=";")
  c(peakID=as.character(x),rmskOverlap=rname)
}))

annot_list <- merge(annot_list, rmsk_list, by.x="peakID", by.y="peakID")
annot_list <- unique(annot_list) #remove duplicate entries (from multiple overlaps with annotations)

write.table(annot_list, file=paste0(opt$output,"_p",opt$promoter,"_dREGpeak_annotation.txt"), sep="\t", quote=F, row.names=F, col.names=T)


#Prepare detailed intragenic sheet
active_intra$localGeneStatus <- "active"
inactive_intra$localGeneStatus <- "inactive"

active_intra$TSSdominance <- "dominant"
inactive_intra$TSSdominance <- "below threshold"

inactive_intra <- inactive_intra %>% mutate(TSSdominance = case_when(peakID %in% suppress_promprox$peakID ~ "suppressed")) #if it overlaps with a suppressed biotype/TSL transcript promoter

intra <- rbind(active_intra,inactive_intra)
intra <- merge(intra, keep_gene, by.x="geneID", by.y="geneID") #Add gene info (with adjusted coords)
intra <- merge(intra, all_tssP, by.x="geneID", by.y="geneID") #Add obsTSS's
intra <- intra[,!(names(intra) %in% c("peak_strand", "gene_chr", "gene_start", "gene_end", "gene_score", "gene_strand", "type", "score", "phase"))] #remove excess/redundant columns
intra$geneName <- paste0("=\"",intra$geneName,"\"") #make gene names Excel safe
intra <- intra[c("peakID", "peak_chr", "peak_start", "peak_end", "peak_center_start", "peak_center_end", "peak_score", "localGeneStatus", "TSSdominance", "geneID", "geneName", "source", "biotype", "chr", "start", "end", "strand", "obsTSS", "attributes")] #reorder columns
intra <- intra[order(intra$peakID),]
intra <- unique(intra) #remove duplicate entries from multiple overlaps
row.names(intra) <- 1:nrow(intra)
write.table(intra, file=paste0(opt$output,"_p",opt$promoter,"_intragenic.txt"), sep="\t", quote=F, row.names=F, col.names=T)

message("dREG peak categorization complete.")


#####DISTANCE METRICS AND QC PLOTS#####
message("Generating distance metrics and QC plots...")
#Distance to closest features for QA/QC plots
dist_tss <- bedTools("closest", dREG, active_tssP, "-k 1 -d -t first")
dist_gb <- bedTools("closest", dREG, active_bed, "-k 1 -d -t first")

#Plot distance to closest peak
pdf(paste0(opt$output,"_p",opt$promoter,"_m", opt$max,"_gene_histograms.pdf"))
hist(dist_tss[which(dist_tss$distance>=0),]$distance, main="Promoter Proximal Distance to TSS", xlab="Distance to TSS (bp)", ylab="Number of dREG Peaks (50bp bins)", breaks=seq(0,max(dist_tss$distance)+50, by=50), xlim=c(0,5000), ylim=c(0,5000)); abline(v=opt$promoter, col="red", lwd=2, lty=2)
hist(dist_tss[which(dist_tss$distance>=0),]$distance, main="Distal Peak Distance to TSS", xlab="Distance to TSS (bp)", ylab="Number of dREG Peaks (5000bp bins)", breaks=seq(0,max(dist_tss[which(dist_tss$distance>=0),]$distance)+5000, by=5000), xlim=c(0,opt$max*2), ylim=c(0,5000)); abline(v=opt$max, col="red", lwd=2, lty=2)
hist(dist_gb[which(dist_gb$distance>=0),]$distance, main="Proximal Peaks to Gene Bodies", xlab="Distance to GB (bp)", ylab="Number of dREG Peaks (50bp bins)", breaks=seq(0,max(dist_tss[which(dist_tss$distance>=0),]$distance)+50, by=50), xlim=c(0,5000), ylim=c(0,1000))
hist(dist_gb[which(dist_gb$distance>=0),]$distance, main="Distal Peaks to Gene Bodies", xlab="Distance to GB (bp)", ylab="Number of dREG Peaks (5000bp bins)", breaks=seq(0,max(dist_tss[which(dist_tss$distance>=0),]$distance)+5000, by=5000), xlim=c(0,opt$max*2), ylim=c(0,1000)); abline(v=opt$max, col="red", lwd=2, lty=2)
dev.off()


#####TWO CLOSEST ACTIVE GENES TO EACH dREG PEAK#####
message("Defining closest genes to distal dREG peaks...")
#Distal peak to gene assignments (closest 2)
if(opt$distance == "TSS")
{
  dist1 <- bedTools("closest", dist.bed, active_tssP, "-k 1 -D b -t first")
  dist2 <- bedTools("closest", dist.bed, active_tssP, "-k 2 -D b -t all")
}
if(opt$distance == "GB")
{
  dist1 <- bedTools("closest", dist.bed, active_bed, "-k 1 -D b -t first")
  dist2 <- bedTools("closest", dist.bed, active_bed, "-k 2 -D b -t all")
}

#Filter distances greater than max
dist1 <- filter(dist1, distance <= opt$max)
dist2 <- filter(dist2, distance <= opt$max)

#Filter out entries in dist2 matching dist1 (dist2nd) then filter by first entry per peakID
dist2nd <- dplyr::anti_join(dist2, dist1)
dist2nd <- dist2nd[!(duplicated(dist2nd$peakID)),]

#Produce enhancer-centric list of gene associations (2 closest genes or NA if 2nd closest is further away than max)
gene_ann <- subset(keep_gene, select = c("geneID","geneName","chr","start","end","strand","biotype","attributes")) #gene annotation features

#Closest gene table converted to output format
dist_out <- dist1[,(names(dist1) %in% c("peakID","peak_chr","peak_start","peak_end","peak_center_start","peak_center_end","peak_score","geneID","distance"))]
names(dist_out)[names(dist_out) == "peakID"] <- "peak_ID"
names(dist_out)[names(dist_out) == "geneID"] <- "gene1_geneID"
names(dist_out)[names(dist_out) == "distance"] <- paste0("gene1_peak_distance_to_",opt$distance)

dist_out <- merge(dist_out, gene_ann, all.x = TRUE, by.x="gene1_geneID", by.y="geneID", sort=FALSE) #add gene features to each row
dist_out <- merge(dist_out, obs_tssP, all.x = TRUE, by.x="gene1_geneID", by.y="geneID", sort=FALSE)

names(dist_out)[names(dist_out) == "geneName"] <- "gene1_name"
names(dist_out)[names(dist_out) == "biotype"] <- "gene1_biotype"
names(dist_out)[names(dist_out) == "attributes"] <- "gene1_attributes"
names(dist_out)[names(dist_out) == "obsTSS"] <- "gene1_obsTSS"
names(dist_out)[names(dist_out) == "reads.TSS"] <- "gene1_reads@obsTSS"
names(dist_out)[names(dist_out) == "chr"] <- "gene1_chr"
names(dist_out)[names(dist_out) == "start"] <- "gene1_start"
names(dist_out)[names(dist_out) == "end"] <- "gene1_end"
names(dist_out)[names(dist_out) == "strand"] <- "gene1_strand"

#Second closest gene table
dist2nd <- dist2nd[,(names(dist2nd) %in% c("peakID","peak_chr","peak_start","peak_end","peak_center_start","peak_center_end","peak_score","geneID","distance"))]
names(dist2nd)[names(dist2nd) == "peakID"] <- "peak_ID"
names(dist2nd)[names(dist2nd) == "geneID"] <- "gene2_geneID"
names(dist2nd)[names(dist2nd) == "distance"] <- paste0("gene2_peak_distance_to_",opt$distance)

dist2nd <- merge(dist2nd, gene_ann, all.x = TRUE, by.x="gene2_geneID", by.y="geneID", sort=FALSE) #add gene features to each row
dist2nd <- merge(dist2nd, obs_tssP, all.x = TRUE, by.x="gene2_geneID", by.y="geneID", sort=FALSE)

names(dist2nd)[names(dist2nd) == "geneName"] <- "gene2_name"
names(dist2nd)[names(dist2nd) == "biotype"] <- "gene2_biotype"
names(dist2nd)[names(dist2nd) == "attributes"] <- "gene2_attributes"
names(dist2nd)[names(dist2nd) == "obsTSS"] <- "gene2_obsTSS"
names(dist2nd)[names(dist2nd) == "reads.TSS"] <- "gene2_reads@obsTSS"
names(dist2nd)[names(dist2nd) == "chr"] <- "gene2_chr"
names(dist2nd)[names(dist2nd) == "start"] <- "gene2_start"
names(dist2nd)[names(dist2nd) == "end"] <- "gene2_end"
names(dist2nd)[names(dist2nd) == "strand"] <- "gene2_strand"

dist_out <- merge(dist_out, dist2nd, all.x = TRUE, sort=FALSE)

dist_out <- dist_out[order(dist_out$peak_chr,dist_out$peak_start),]
row.names(dist_out) <- 1:nrow(dist_out)


dist_out$peak_genomicCategory <- NA #Should not happen
dist_out <- dist_out %>% mutate(peak_genomicCategory = case_when(peak_ID %in% interg$peakID ~ "Intergenic Transcribed Locus",
                                                                peak_ID %in% intra$peakID ~ "Intragenic Transcribed Locus",
                                                                peak_ID %in% promprox$peakID ~ "Suppressed Promoter Proximal"))

dist_out$peak_overlappingGeneStatus <- NA #Can happen
dist_out <- dist_out %>% mutate(peak_overlappingGeneStatus = case_when(peak_ID %in% active_promprox$peakID ~ "Active",
                                                                                                  peak_ID %in% nondom_promprox$peakID ~ "Non-Dominant",
                                                                                                  peak_ID %in% inactive_promprox$peakID ~ "Inactive",
                                                                                                  peak_ID %in% active_intra$peakID ~ "Active",
                                                                                                  peak_ID %in% inactive_intra$peakID ~ "Inactive"))

#Make gene names Excel safe
dist_out$gene1_name <- paste0("=\"",dist_out$gene1_name,"\"")
dist_out$gene2_name <- paste0("=\"",dist_out$gene2_name,"\"")

#Reorder data table
dist_out <- dist_out[c("peak_ID","peak_chr","peak_start","peak_end","peak_center_start","peak_center_end","peak_score","peak_genomicCategory","peak_overlappingGeneStatus",
                       "gene1_geneID","gene1_name","gene1_chr","gene1_start","gene1_end","gene1_strand","gene1_obsTSS","gene1_reads@obsTSS",paste0("gene1_peak_distance_to_",opt$distance),"gene1_biotype","gene1_attributes",
                       "gene2_geneID","gene2_name","gene2_chr","gene2_start","gene2_end","gene2_strand","gene2_obsTSS","gene2_reads@obsTSS",paste0("gene2_peak_distance_to_",opt$distance),"gene2_biotype","gene2_attributes")] #reorder columns

dist_out <- unique(dist_out) #remove duplicate entries introduced by multiple merges

write.table(dist_out, file=paste0(opt$output,"_p",opt$promoter,"_m",opt$max,"_distal_dREGpeak_closestGenes.txt"), sep="\t", quote=F, row.names=F, col.names=T) #Save dREG peak-centric data table


#####GENE CENTRIC DISTANCE ANALYSIS#####
message("Generating gene-centric analysis...")
#Produce gene-centric list of first and second closest genes to each dREG peak
genelist <- gene_ann$geneID
grepout <- data.frame(peakIDs=character(),stringsAsFactors=FALSE)
i=1
while(i<=length(genelist))
{
  genegrep <- dist2[grep(genelist[i], dist2$geneID),]
  if(nrow(genegrep)>=1)
  {
    geneline <- data.frame(geneID=genelist[i])
    geneline$peakIDs <- paste0(genegrep$peakID,collapse=",")
    geneline$peak_count <- nrow(genegrep)
    geneline$closest_peak <- genegrep[which.min(abs(genegrep$distance)),]$peakID
    geneline$closest_distance <- min(abs(genegrep$distance))
    grepout <- rbind(grepout,geneline)
  }
  i=i+1
}
grepout <- merge(grepout,gene_ann)
grepout <- grepout[c("geneID","geneName","biotype","attributes","chr","start","end","strand","closest_peak","closest_distance","peak_count","peakIDs")]
grepout$geneName <- paste0("=\"",grepout$geneName,"\"")
grepout <- unique(grepout) #remove duplicated entries from grep output
write.table(grepout, file=paste0(opt$output,"_p",opt$promoter,"_m",opt$max,"_gene-centric_closest-dREGpeaks.txt"), sep="\t", quote=F, row.names=F, col.names=T)


#####dREG PEAK CLUSTER ANALYSIS#####
message("Beginning dREG peak cluster analysis...")
#Calculate distance to closest neighbor dREG peak and produce histogram
clust_dist <- bedTools("closest", dREG, dREG, "-k 1 -d -io -t first") #Note: columns are mislabeled due to comparing dREG peak list to itself. Column 17 represents distances to closest neighbor.
clust_dist <- clust_dist[which(clust_dist[,17]>=0),][,17] #-1 is assigned to those with no closest neighbor (coordinates on chromosome fragments)
pdf(paste0(opt$output,"_g",opt$gap,"_neighbor_histograms.pdf"))
hist(clust_dist, main="Distance to Closest dREG Peak Neighbor", xlab="Distance to Neighbor (bp)", ylab="Number of dREG Peaks (50bp bins)", breaks=seq(0,max(clust_dist)+50,by=50), xlim=c(0,2000), ylim=c(0,1000)); abline(v=1000, col="red", lwd=2, lty=2)
hist(clust_dist, main="Distance to Closest dREG Peak Neighbor", xlab="Distance to Neighbor (bp)", ylab="Number of dREG Peaks (100bp bins)", breaks=seq(0,max(clust_dist)+100,by=100), xlim=c(0,5000), ylim=c(0,5000)); abline(v=1000, col="red", lwd=2, lty=2)
hist(clust_dist, main="Distance to Closest dREG Peak Neighbor", xlab="Distance to Neighbor (bp)", ylab="Number of dREG Peaks (500bp bins)", breaks=seq(0,max(clust_dist)+500,by=500), xlim=c(0,20000), ylim=c(0,10000)); abline(v=1000, col="red", lwd=2, lty=2)
dev.off()

#Generate peak clusters
message("Defining dREG peak clusters...")
clusters <- bedTools("merge",dREG,dREG,paste0("-d ",opt$gap," -c 4,4 -o count,distinct"))
colnames(clusters)=c("cluster_chr","cluster_start","cluster_end","peak_count","peak_IDs")

pdf(paste0(opt$output,"_g", opt$gap,"_peak_clusters.pdf"))
hist(clusters$peak_count, main="Number of Peaks per Cluster", xlab="Number of Peaks in Cluster", ylab="Number of Clusters", breaks=seq(0,max(clusters$peak_count),by=1),labels=TRUE) #peaks/cluster count histogram
dev.off()

clusters <- clusters[which(clusters$peak_count>1),] #filter out clusters with only 1 peak

#Name clusters
clusterIDnum <- sprintf("%05d", seq(1,nrow(clusters), by=1))
clusterIDs <- paste0("dREGcluster", clusterIDnum)
clusters$clusterID <- clusterIDs
clusters <- clusters[c("clusterID","cluster_chr","cluster_start","cluster_end","peak_count","peak_IDs")]
clusters.bed <- data.frame("chr"=clusters$cluster_chr)
clusters.bed$start <- clusters$cluster_start
clusters.bed$end <- clusters$cluster_end
clusters.bed$id <- clusters$clusterID
clusters.bed$score <- clusters$peak_count
clusters.bed$stand <- "+"
clusters.bed$thickstart <- clusters$cluster_start
clusters.bed$thickend <- clusters$cluster_end
clusters.peakIDs <- data.frame("clusterID"=clusters$clusterID,"peak_IDs"=clusters$peak_IDs) #peakID lookup table

#Annotate clusters by overlap with annotations
#active_tss > nondom_tss > unstable_tss > comp_tss > active_bed > comp_gb
active_promClust <- bedTools("intersect", clusters.bed, active_tss, "-wao")
catch <- active_promClust[active_promClust$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)]
active_promClust <- active_promClust[active_promClust$geneID!=".",]

nondom_promClust <- bedTools("intersect", catch, nondom_tss, "-wao")
catch <- nondom_promClust[nondom_promClust$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)]
nondom_promClust <- nondom_promClust[nondom_promClust$geneID!=".",]

unstable_promClust <- bedTools("intersect", catch, unstable_tss, "-wao")
catch <- unstable_promClust[unstable_promClust$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)]
unstable_promClust <- unstable_promClust[unstable_promClust$geneID!=".",]

inactive_promClust <- bedTools("intersect", catch, keep_tss, "-wao")
catch <- inactive_promClust[inactive_promClust$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)]
inactive_promClust <- inactive_promClust[inactive_promClust$geneID!=".",]

active_intraClust <- bedTools("intersect", catch, active_bed, "-wao")
catch <- active_intraClust[active_intraClust$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)]
active_intraClust <- active_intraClust[active_intraClust$geneID!=".",]
catch2 <- active_intraClust[,!(names(active_intraClust) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch2 <- catch2[-c(9:15)]

inactive_intraClust <- bedTools("intersect", catch, keep_gb, "-wao")
catch <- inactive_intraClust[inactive_intraClust$geneID==".",]
catch <- catch[,!(names(catch) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch <- catch[-c(9:15)]
inactive_intraClust <- inactive_intraClust[inactive_intraClust$geneID!=".",]
catch3 <- inactive_intraClust[,!(names(inactive_intraClust) %in% c("gene_chr", "gene_start", "gene_end", "geneID", "gene_score", "gene_strand", "overlap"))]
catch3 <- catch3[-c(9:15)]

if(opt$distance=="TSS")
{
  intergClust <- unique(bedTools("closest", catch, active_tssP, "-k 1 -D b -t first"))
  active_intraClosestActive <- unique(bedTools("closest", catch2, active_tssP, "-k 1 -D b -t first"))
  inactive_intraClosestActive <- unique(bedTools("closest", catch3, active_tssP, "-k 1 -D b -t first"))
}
if(opt$distance == "GB")
{
  intergClust <- unique(bedTools("closest", catch, active_bed, "-k 1 -D b -t first"))
  active_intraClosestActive <- unique(bedTools("closest", catch2, active_bed, "-k 1 -D b -t first"))
  inactive_intraClosestActive <- unique(bedTools("closest", catch3, active_bed, "-k 1 -D b -t first"))  
}

#Make concise gene annotation look-up table
gene_ann_tab <- unique(gene_ann[,c(1:2,7)])
colnames(gene_ann_tab) <- c("overlapping_geneID","overlapping_geneName","overlapping_biotype")
gene_ann_tab$overlapping_geneName <- paste0("=\"",gene_ann_tab$overlapping_geneName,"\"")

#Create cluster annotation and data file
annotClust <- active_promClust[1:5]
annotClust$TSSproximity <- "Dominant"
annotClust$overlapping_geneStatus <- "Active"
annotClust$genomicCategory <- "Promoter Proximal Cluster"
annotClust$overlapping_geneID <- active_promClust$geneID
annotClust <- merge(annotClust, gene_ann_tab, sort=FALSE)
annotClust$closestActive_geneID <- annotClust$overlapping_geneID
annotClust$closestActive_geneName <- annotClust$overlapping_geneName
annotClust$closestActive_biotype <- annotClust$overlapping_biotype
annotClust$closestActive_distance <- 0

listClust <- nondom_promClust[1:5]
listClust$TSSproximity <- "Non-Dominant"
listClust$overlapping_geneStatus <- "Active"
listClust$genomicCategory <- "Promoter Proximal Cluster"
listClust$overlapping_geneID <- nondom_promClust$geneID
listClust <- merge(listClust, gene_ann_tab, sort=FALSE)
listClust$closestActive_geneID <- listClust$overlapping_geneID
listClust$closestActive_geneName <- listClust$overlapping_geneName
listClust$closestActive_biotype <- listClust$overlapping_biotype
listClust$closestActive_distance <- 0
annotClust <- rbind(listClust, annotClust)

listClust <- unstable_promClust[1:5]
listClust$TSSproximity <- "Non-Dominant"
listClust$overlapping_geneStatus <- "Unstable"
listClust$genomicCategory <- "Promoter Proximal Cluster"
listClust$overlapping_geneID <- unstable_promClust$geneID
listClust <- merge(listClust, gene_ann_tab, sort=FALSE)
listClust$closestActive_geneID <- NA
listClust$closestActive_geneName <- NA
listClust$closestActive_biotype <- NA
listClust$closestActive_distance <- NA
annotClust <- rbind(listClust, annotClust)

listClust <- inactive_promClust[1:5]
listClust$TSSproximity <- "Non-Dominant"
listClust$overlapping_geneStatus <- "Inactive"
listClust$genomicCategory <- "Promoter Proximal Cluster"
listClust$overlapping_geneID <- inactive_promClust$geneID
listClust <- merge(listClust, gene_ann_tab, sort=FALSE)
listClust$closestActive_geneID <- NA
listClust$closestActive_geneName <- NA
listClust$closestActive_biotype <- NA
listClust$closestActive_distance <- NA
annotClust <- rbind(listClust, annotClust)

#Updated 11/11/21 AF
listClust <- active_intraClust[1:5]
listClust$TSSproximity <- "Distal"
listClust$overlapping_geneStatus <- "Active"
listClust$genomicCategory <- "Intragenic Cluster"
listClust$overlapping_geneID <- active_intraClust$geneID
listClust <- merge(listClust, gene_ann_tab, sort=FALSE)
listClust <- unique(listClust) #remove redundant entries
listClust$closestActive_geneID <- merge(listClust, active_intraClosestActive, by.x="peakID", by.y="peakID", all.x=T, sort=F)$geneID
listClust$closestActive_geneName <- merge(listClust, gene_ann_tab, by.x="closestActive_geneID", by.y="overlapping_geneID", all.x=T, sort=F)$overlapping_geneName.y
listClust$closestActive_biotype <- merge(listClust, gene_ann_tab, by.x="closestActive_geneID", by.y="overlapping_geneID", all.x=T, sort=F)$overlapping_biotype.y
listClust$closestActive_distance <- merge(listClust, active_intraClosestActive, by.x="peakID", by.y="peakID", all.x=T, sort=F)$distance
listClust <- unique(listClust) #remove redundant entries
annotClust <- rbind(listClust, annotClust)

#Updated 11/11/21 AF
listClust <- inactive_intraClust[1:5]
listClust$TSSproximity <- "Distal"
listClust$overlapping_geneStatus <- "Inactive"
listClust$genomicCategory <- "Intragenic Cluster"
listClust$overlapping_geneID <- inactive_intraClust$geneID
listClust <- merge(listClust, gene_ann_tab, sort=FALSE)
listClust$closestActive_geneID <- merge(listClust, inactive_intraClosestActive, by.x="peakID", by.y="peakID", all.x=T, sort=F)$geneID
listClust$closestActive_geneName <- merge(listClust, gene_ann_tab, by.x="closestActive_geneID", by.y="overlapping_geneID", all.x=T, sort=F)$overlapping_geneName.y
listClust$closestActive_biotype <- merge(listClust, gene_ann_tab, by.x="closestActive_geneID", by.y="overlapping_geneID", all.x=T, sort=F)$overlapping_biotype.y
listClust$closestActive_distance <- merge(listClust, inactive_intraClosestActive, by.x="peakID", by.y="peakID", all.x=T, sort=F)$distance
listClust <- unique(listClust) #remove redundant entries
annotClust <- rbind(listClust, annotClust)

#Remove intergenic clusters on chromosome fragments (unable to find closest active gene)
intergClust <- intergClust[intergClust$geneID!=".",]

#Updated 11/11/21 AF
listClust <- intergClust[1:5]
listClust$TSSproximity <- "Distal"
listClust$overlapping_geneStatus <- NA
listClust$genomicCategory <- "Intergenic Cluster"
listClust$overlapping_geneID <- NA
listClust$overlapping_geneName <- NA
listClust$overlapping_geneName <- NA
listClust$overlapping_biotype <- NA
listClust$closestActive_geneID <- merge(listClust, intergClust, by.x="peakID", by.y="peakID", all.x=T, sort=F)$geneID
listClust$closestActive_geneName <- merge(listClust, gene_ann_tab, by.x="closestActive_geneID", by.y="overlapping_geneID", all.x=T, sort=F)$overlapping_geneName.y
listClust$closestActive_biotype <- merge(listClust, gene_ann_tab, by.x="closestActive_geneID", by.y="overlapping_geneID", all.x=T, sort=F)$overlapping_biotype.y
listClust$closestActive_distance <- merge(listClust, intergClust, by.x="peakID", by.y="peakID", all.x=T, sort=F)$distance
listClust <- unique(listClust) #remove redundant entries
annotClust <- rbind(listClust, annotClust)

#Rename columns peaks->clusters
names(annotClust)[names(annotClust) == "peak_chr"] <- "cluster_chr"
names(annotClust)[names(annotClust) == "peak_start"] <- "cluster_start"
names(annotClust)[names(annotClust) == "peak_end"] <- "cluster_end"
names(annotClust)[names(annotClust) == "peakID"] <- "clusterID"
names(annotClust)[names(annotClust) == "peak_score"] <- "peak_count"
annotClust <- unique(annotClust)
annotClust <- merge(annotClust, clusters.peakIDs, by.x = "clusterID", by.y = "clusterID", sort=FALSE, all=TRUE)

annotClust <- annotClust[c("clusterID", "cluster_chr", "cluster_start", "cluster_end", "peak_count", "peak_IDs", "TSSproximity", "genomicCategory", "overlapping_geneStatus", "overlapping_geneID", "overlapping_geneName", "overlapping_biotype", "closestActive_geneID", "closestActive_geneName", "closestActive_biotype", "closestActive_distance")]
annotClust <- annotClust[order(annotClust$cluster_chr,annotClust$cluster_start),]
annotClust <- unique(annotClust) #Shouldn't be duplicates here, but just in case they are...
write.table(annotClust, file=paste0(opt$output,"_p",opt$promoter,"_m",opt$max,"_g",opt$gap,"_dREGpeakClusters.txt"), sep="\t", quote=F, row.names=F, col.names=T)

save.image(paste0(opt$output,"_dREGannotate.RData"))
message("dREG annotation script completed successfully!")
