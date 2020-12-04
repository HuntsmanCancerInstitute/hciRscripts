#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

## using -g returns WARNING: unknown gui ...  so no --genome option

opts <-  list(
    make_option(c("-s", "--sequencing"), default="single",
          help="single, paired or miRNA sequencing"),
    make_option(c("-r", "--run"), default=NULL,
          help="Run ID for report title"),
    make_option(c("-i", "--id"), default=NULL,
          help="Sample ID prefix, will use run ID with X1 ending if missing"),
    make_option(c("-f", "--fastq"), default=NULL,
          help="Sample fastq file name, will check /Repository/MicroarrayData with run ID if missing"),
    make_option(c("-d", "--database"), default="human",
          help="Reference database, default human or mouse, elephant, fly, worm, pig,
     rat, rabbit, sheep, vervet, yeast, zebrafish"),
    make_option(c("-v", "--version"), default="102",
          help="Ensembl release, default 102"),
   make_option(c("-p", "--adapters"), default="TruSeq",
         help="Adapter sequences, default TruSeq or Nextera for new RiboZero kits"),
    make_option(c("-n", "--ncpu"), default="24",
          help="Number CPUs, pysano will use the maximum number, default 24"),
    make_option(c("-a", "--align"), default="Alignments",
          help="Directory name with cmd.txt output, default Alignments"),
     make_option(c("-l", "--length"), default="50",
        help="Read length, default 50")
)

parser <- OptionParser(option_list=opts, description = "
Creates a README.Rmd file for single, paired or miRNA sequencing workflows at HCI using defaults
in the cmd.txt files from setup_jobs.R.  Run 'render.R -f README.Rmd' to create the html report.
")
opt <- parse_args(parser)

if( is.null(opt$run) ){
   print_help(parser)
   quit(status=1)
}

# Path to hciR script
install_dir <- "/home/BioApps/hciR"

## checks
if(file.exists("README.Rmd")) stop("README.Rmd file already exists")
if(opt$sequencing == "novaseq") opt$sequencing <- "clumpify"
if(!opt$sequencing %in% c("single", "paired", "miRNA", "qiagen", "metagenome", "metatranscriptome", "clumpify", "microbe")) stop("No README template found")
## if sample ID is missing, convert run 14980R to 14980X1
if(is.null(opt$id)) opt$id <- gsub("R.*", "X1", opt$run)
## if fastq is missing, search /Repository/MicroarrayData
if(is.null(opt$fastq)){
   opt$fastq <- Sys.glob( paste0("/Repository/MicroarrayData/*/", opt$run, "/Fastq/", opt$id, "_*gz"))
   if(length(opt$fastq) == 0) stop("No match to fastq files in /Repository/MicroarrayData, please add --fastq option")
   opt$fastq <- gsub(".*/", "", opt$fastq)  # remove path
   ## combine paried end into "fastq1.fq fastq2.fq"
   if(length(opt$fastq) == 2 ) opt$fastq <- paste(opt$fastq, collapse= " @NEWLINE ")
}
if(!opt$length %in% c("50", "125")) message("Length should be 50 or 125.  Please check if star", opt$length, " exists")
sjdb <- as.numeric(opt$length) -1

## Load STAR reference table
x <- read.delim( paste0( install_dir, "/STAR_ref_dbs.txt"), stringsAsFactors=FALSE)
n <- grep(opt$database, x$name, ignore=TRUE)
if(length(n)!=1) stop("Database should match human, mouse, fly, worm, pig, rat, rabbit, sheep, yeast or zebrafish.")
name <- x$name[n]
species <- x$species[n]
assembly <- x$assembly[n]
release <- as.numeric(opt$version)
STAR_version <- "2.7.6a"
if(release == 100) STAR_version <- "2.7.3a"
if(release == 98) STAR_version <- "2.7.2c"
if(release == 96) STAR_version <- "2.7.0f"
if(release == 94) STAR_version <- "2.6.1b"
if(release == 92) STAR_version <- "2.5.4a"
if(release == 90) STAR_version <- "2.5.2b"

# adapters
if(tolower(opt$adapters) == "nextera"){
   adapt1 <- "CTGTCTCTTATACACATCT"
   adapt2 <- "CTGTCTCTTATACACATCT"
}else{
   # Illumina TruSeq adapters
   adapt1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
   adapt2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
}

## fasta file name on FTP
dna <- "toplevel"
if(name %in% c("mouse", "human", "zebrafish")) dna <- "primary_assembly"
## old assemblies
if(assembly == "GRCz11" & release == 90) assembly <- "GRCz10"
if(assembly == "BDGP6.22" & release < 96) assembly <- "BDGP6"


## Load R markdown template and replace @Variables
rmd_txt <- paste0(install_dir, "/templates/README_", opt$sequencing, ".Rmd")
rmd <- readr::read_lines(rmd_txt, skip=2)
rmd <- stringr::str_replace_all(rmd, c(
   `@NAME`     = name,
   `@VERSION`  = opt$version,
   `@LEVEL`    = dna,
   `@SAMPLE`   = opt$id,
   `@FASTQ`    = opt$fastq,
   `@NEWLINE`  = "\\\\
 ",
   `@RUN`      = opt$run,
   `@LENGTH`   = opt$length,
   `@SJDB`     = sjdb,
   `@SPECIES`  = species ,
   `@LOWER_SPECIES`= tolower(species),
   `@ASSEMBLY` = assembly ,
   `@ADAPT1`   = adapt1,
   `@ADAPT2`   = adapt2,
   `@STAR`     = STAR_version,
   `@NCPU`     = opt$ncpu,
   `@ALIGNDIR` = opt$align  ))

readr::write_lines(rmd, "README.Rmd")
message("Created README.Rmd markdown file for ", name, " star", opt$length)
message("Run 'render.R -f README.Rmd' to create an html report")
