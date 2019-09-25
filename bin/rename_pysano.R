#!/usr/bin/env Rscript



suppressPackageStartupMessages(library("optparse"))

opts <-  list(
   make_option(c("-d", "--directory"), default="NA",
         help="New directory name, required"),
   make_option(c("-s", "--separator"), default="._",
         help="Parse ID to separator, default . or _"),
   make_option(c("-p", "--pattern"), default="^20",
         help="Only move analysis directories matching this pattern, default ^20")
)

parser <- OptionParser(option_list=opts, description = "
Running pysano with the #a option will output to separate directories like
2017_04_21_16_15_2.892994.  This script will rename that directory using the
first part of the bam file name to <directory>/<bam> like Alignments/14180X1"  )
  opt <- parse_args(parser)

  if( opt$directory == "NA" ){
     print_help(parser)
     quit(status=1)
  }

x <- list.files(".", "*.bam$", recursive=TRUE)
x <- grep("Aligned.*bam$", x, value=TRUE, invert=TRUE)
if(length(x)==0 ) stop("No unique bam files found")

x <- grep(opt$pattern, x, value=TRUE)
if(length(x)==0 ) stop("No directories matching pattern",  opt$pattern)

y <- strsplit(x, "/")

old_dir <- sapply(y, "[[", 1)
if( any(duplicated(old_dir)) ) stop("No unique ID to rename directory since more than 1 bam file present")

old_file <- sapply(y, "[[", 2)
## use id before _ or .
id  <- gsub( paste0("[", opt$separator, "].*"), "", old_file)
if(length(id)==0) stop("Problem parsing IDs")
if(any(nchar(id)==0)) stop("Some IDs are empty strings")
new_dir <- paste0(opt$directory,  "/", id)

if(!dir.exists(opt$directory))  dir.create( opt$directory)
for(i in seq_along(x)){
message("Moving ", old_dir[i], " to ", new_dir[i] )
   file.rename( old_dir[i], new_dir[i])
}
