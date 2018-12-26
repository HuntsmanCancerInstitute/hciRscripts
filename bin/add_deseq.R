#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))

opts <-  list(
   make_option(c("-r", "--run"), default=NULL,
        help="Run ID for report title"),
   make_option(c("-d", "--database"), default="human",
        help="Annotation database, either human, mouse, fly, worm, pig, rat, rabbit,
                sheep, yeast or zebrafish, default human"),
   make_option(c("-v", "--version"), default="92",
        help="Ensembl release, default 92"),
   make_option(c("-s", "--samples"), default="samples.txt",
        help="Tab-delimited file with ids in column 1 matching count column names
                and a treatment column for contrasts, default samples.txt"),
   make_option(c("-t", "--trt"), default="trt",
        help="Name of treatment column in sample table, default trt."),
   make_option(c("-c", "--counts"), default="counts.txt",
        help="Tab-delimited count matrix file, default counts.txt"),
    make_option(c("-f", "--filter"), default=5,
        help="Low count cutoff to filter counts, default 5"),
    make_option(c("-p", "--padj"), default=0.05,
        help="Adjusted p-value cutoff, default 0.05"),
    make_option(c("-x", "--vs"), default="all",
        help="Compare groups to a specific treatment, default is all vs. all"),
    make_option(c("-l", "--levels"), default=NULL,
        help="A comma-separated list to reorder treatments.  By default treatments are sorted
                alphabetically, so use 'C,B,A' to compare C vs A, C vs B and B vs A"),
    make_option(c("-m", "--mouseover"), default="id",
        help="A comma-separated list of sample column names for tooltips in PCA plot, default is id column")
)

parser <- OptionParser(option_list=opts, description = "
Create a DESeq markdown file with command to run DESeq2
")
opt <- parse_args(parser)

if(is.null(opt$run)){
   print_help(parser)
   quit(status=1)
}
## use Rmd file in hciR package if missing
if(file.exists("DESeq.Rmd")) stop("DESeq.Rmd file already exists")

# Path to hciR script
## https://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script

cmdArgs <- commandArgs(trailingOnly = FALSE)
n <- grep("--file=", cmdArgs)
install_dir <-  normalizePath(gsub("--file=|/bin/.*", "", cmdArgs[n]))

## TO DO - allow 2 or more columns for trt (and paste together in design formula ~trt + cell)

# check number of samples for rlog or vst?
useVST <- FALSE
samples <- suppressMessages(read_tsv(opt$samples))
if(nrow(samples)>40){
   message("More than 40 samples, using variance stabilizing transform and not regularized logs")
   useVST <- TRUE
}
## check number of contrasts and add more result sections?


## set fig.height and width in sample dist plot?
figwidth <- round(nrow(samples)/3) + 2
figwidth[figwidth > 8] <- 8
figheight <- figwidth * 0.75

## convert "A,B,C" to c("A", "B", "C") for markdown
if(!is.null(opt$levels)) opt$levels <- capture.output(dput(strsplit(opt$levels, ", *")[[1]]))

if(!opt$trt %in% names(samples)) stop("Trt column name is missing from sample table")
x <- strsplit(opt$mouseover, ", *")[[1]]
if(!all(x %in% names(samples))) stop("Mouseover names are missing from sample table")
opt$mouseover <- capture.output(dput(x))

db <- c("human", "mouse", "fly", "sheep", "pig", "rat", "rabbit", "worm", "yeast", "zebrafish")
if(!opt$database %in% db){
   stop("Database should match human, mouse, fly, pig, rat, rabbit, sheep,  worm, yeast or zebrafish.")
}
db1 <- paste0(opt$database, opt$version)

if(!opt$version %in% c("90", "92")) message("WARNING: only Ensembl version 90 and 92 are available in hciR_data.
You may need to update the gene annotation section and use read_biomart instead")
if(is.null(opt$design)) opt$design <- opt$trt

rmd <- readr::read_lines(paste0(install_dir, "/templates/DESeq_template.Rmd"))

## if null, delete line 23  with sample$@trt <- factor(samples$@trt, levels = @levels)
if(is.null(opt$levels))  rmd <- rmd[!grepl("^samples\\$@trt", rmd)]
## replace rlog with vst
if(useVST){
   rmd <- stringr::str_replace_all(rmd, c(
     `regularized log` = "variance stabilizing",
     rlog  = "vst",
     rld   = "vsd"))
}
rmd <- stringr::str_replace_all(rmd, c(
`@run`       = opt$run,
`@name`      = opt$database,
`@db`        = db1,
`@samples`   = opt$samples,
`@counts`    = opt$counts,
`@version`   = opt$version,
`@trt`       = opt$trt,
`@filter`    = opt$filter,
`@padj`      = opt$padj,
`@fdr`       = opt$padj*100,
`@vs`        = opt$vs,
`@relevel`   = opt$levels,
`@mouseover` = opt$mouseover,
`@figwidth`  = figwidth,
`@figheight` = figheight))

readr::write_lines(rmd, "DESeq.Rmd")
message("Created DESeq.Rmd markdown file with ", opt$database, " annotations from Ensembl version ", opt$version)
message("Run 'render.R -f DESeq.Rmd' to create an html report and Excel file")
