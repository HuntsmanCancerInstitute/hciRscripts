#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

opts <-  list(
   make_option(c("-d", "--dir"),
         help="Directory name in /tomato/dev/data"),
    make_option(c("-r", "--release"),
       help="Release number")
)

parser <- OptionParser(option_list=opts, description = "
Adds new assembly directory, downloads genome, creates chrom.size and sequence dictionary file.
Be sure to add the new assembly name to STAR_ref_dbs.txt before running.
")
opt <- parse_args(parser)

if( is.null(opt$dir) ){
   print_help(parser)
   quit(status=1)
}

x <- read.delim("/home/BioApps/hciR/STAR_ref_dbs.txt", stringsAsFactors=FALSE)
i <- which(x$dir == opt$dir)
if(length(i)==0) stop("No match to ", opt$dir, " in STAR_ref_dbs.txt")

release <- opt$release   # 90
dir <- opt$dir           # "Pig"
species <- x$species[i]   # "Sus_scrofa"
assembly <- x$assembly[i] # "Sscrofa11.1"

#1. Create directory (as tomatosrvs??)
datadir <- paste0("/tomato/dev/data/", dir, "/", assembly)
if(!file.exists(datadir)){
   message("Creating new directory ", datadir)
   system(paste("mkdir", datadir))
}
setwd(datadir)
# add options to STAR_ref_dbs.txt file in  /home/BioApps/hciR?

#2.   Download genome
#ftp://ftp.ensembl.org/pub/release-90/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
ftp <- paste0("ftp://ftp.ensembl.org/pub/release-", release, "/fasta")
# toplevel or primary ?
dna <- "dna.toplevel"
if(dir %in% c("Human", "Mouse", "Zebrafish")) dna<- "dna.primary_assembly"
 dnafile <- paste0(x$species[i], ".", x$assembly[i], ".", dna, ".fa.gz")
 ftpdir <- paste0(ftp, "/", tolower(x$species[i]), "/dna/",  dnafile)
 message( "Downloading ", dnafile)
 download.file(ftpdir, dnafile,  quiet = TRUE)
 system(paste("gunzip",  dnafile))

##3. Create chrom.size  for bedgraph to bigwig
message("Creating chrom.sizes file")
dnafile <- gsub(".gz$", "", dnafile)

system(paste("/tomato/dev/app/samtools/1.5/samtools faidx", dnafile))  
system(paste0("cut -f1,2 ", dnafile, ".fai |  sort -V  > chrom.sizes"))

##4. Create Sequence Dictionary file

message("Creating sequence dictionary file")
system(paste0("java -jar /tomato/dev/app/picard/1.137/picard.jar CreateSequenceDictionary R=", dnafile,
   " O=", gsub("fa$", "dict", dnafile), " GENOME_ASSEMBLY=", assembly, " SPECIES=", species))
