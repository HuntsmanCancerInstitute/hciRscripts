#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("optparse"))

opts <-  list( make_option(c("-f", "--file"), default=NULL,
           help="R markdown file to render"))
parser <- OptionParser(option_list=opts, description = "
Render an R markdown file
")
opt <- parse_args(parser)

if( is.null(opt$file )){
   print_help(parser)
   quit(status=1)
}
if( !file.exists( opt$file)) stop("No file matching ", opt$file, " found")
render(opt$file)
