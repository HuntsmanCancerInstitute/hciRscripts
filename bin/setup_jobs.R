#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

opts <-  list(
   make_option(c("-e", "--email"), default="NA",
         help="#e email or hci user name for pysano directive, required"),
   make_option(c("-c", "--cluster"), default="kingspeak",
         help="#c cluster name for pysano with at least 64 GB of RAM, default kingspeak"),
   make_option(c("-s", "--sequencing"), default="single",
         help="single, paired or miRNA sequencing"),
   make_option(c("-i", "--input"), default="NA",
         help="Input directory with Fastq files (either absolute path or relative path in current directory)"),
   make_option(c("-f", "--fastq"), default=".gz",
         help="Match fastq files ending in default .gz"),
   make_option(c("-m", "--modified"), default="NA",
         help="Only link fastq files with modified time >= YYYY-MM-DD"),
   make_option(c("-r", "--run"), default="NA",
         help="Run ID in /Repository/MicroarryData, optional"),
   make_option(c("-a", "--analysis"), default="",
         help="Save files to /Repository/AnalysisData, optional"),
   make_option(c("-v", "--version"), default="94",
         help="Ensembl release, default 94, only 90, 92, and 94 available"),
   make_option(c("-d", "--database"), default="human",
         help="Reference database, either human, mouse, fly, worm, pig, rat, rabbit, sheep, yeast or zebrafish, default human"),
    make_option(c("-l", "--length"), default="50",
       help="Read length for STAR reference, default 50 or 125")
)

parser <- OptionParser(option_list=opts, description = "
Creates a cmd.txt file and sample directories with Fastq file links in order to run STAR,
featureCounts and quality metrics on the CHPC clusters using pysano.  The default is to align
to the human reference configured for 50 bp reads in /tomato/dev/data/Human/GRCh38/release92/star50")
  opt <- parse_args(parser)
  if( opt$email == "NA" ){
     print_help(parser)
     quit(status=1)
  }

if(file.exists( "cmd.txt")){
   message("Note: cmd.txt file already exists")
}else{
   if( !grepl("@", opt$email )) opt$email <- paste0( opt$email, "@hci.utah.edu")
   ## STAR version should match version used to create index
   release <- as.numeric(opt$version)
   if(!release %in% c(90, 92, 94)) stop("Version should be 90, 92 or 94")
   STAR_version <- "2.6.1b"
   if(release == 92) STAR_version <- "2.5.4a"
   if(release == 90) STAR_version <- "2.5.2b"

   if( opt$analysis != "" ) opt$analysis <- paste0("#a ", opt$analysis)
   if(!opt$length %in% c("50", "125")) message("Length should be 50 or 125.  Please check if star", opt$length, " exists")
   if(!opt$sequencing %in% c("single", "paired", "miRNA")) stop("Sequencing should be single, paired or miRNA")

   x <- read.delim("/home/BioApps/hciR/STAR_ref_dbs.txt", stringsAsFactors=FALSE)
   refdb <- stringr::str_replace_all(opt$database, c(fly = "Drosophila", fruitfly = "Drosophila",
                  worm = "C_elegans", yeast = "S_cerevisiae", Homo ="Human", Mus = "Mouse"))

   n <- grep(refdb, x$dir, ignore=TRUE)
   if(length(n)!=1) stop("Database should match human, mouse, fly, worm, pig, rat, rabbit, sheep, yeast or zebrafish.")
   species <- x$species[n]
   assembly <- x$assembly[n]
   tomato_dir <- x$dir[n]
   ## old GRCz10 assembly for version 90
   if( release ==90){
       if( assembly == "GRCz11") assembly <- "GRCz10"
    }
   cmd_txt <- paste0("/home/BioApps/hciR/templates/cmd_", opt$sequencing, ".txt")
   cmd <- readr::read_lines(cmd_txt, skip=2)
   cmd <- stringr::str_replace_all(cmd, c(`@EMAIL` = opt$email, `@CLUSTER` = opt$cluster,
    `@ANALYSIS` = opt$analysis, `@LENGTH` =  opt$length, `@DIR` = tomato_dir, `@SPECIES`= species ,
    `@ASSEMBLY`= assembly ,`@VERSION` = release, `@STAR` = STAR_version, `@FASTQ` = opt$fastq  ))
   readr::write_lines(cmd, "cmd.txt")
    message("Created cmd.txt file for ", refdb, " star", opt$length)
}

fastq_end <- paste0(opt$fastq, "$")

if( opt$run == "NA" & opt$input == "NA"){
   message("Note: add a run ID using -r or fastq file directory using -i to link fastq files")
}else{
   if(opt$run != "NA"){
      repo <- Sys.glob( paste0("/Repository/MicroarrayData/*/", opt$run))
      if(length(repo) == 0) stop("No match to ", opt$run, " in /Repository/MicroarrayData")
      fastq_full <- list.files(paste0(repo, "/Fastq"), pattern = fastq_end, full.names =TRUE)
      if(length(fastq_full) == 0) stop("No Fastq files found in Repository" )
   }else{
      fastq_full <- list.files(opt$input, pattern = fastq_end, full.names=TRUE, recursive=TRUE)
      if(length(fastq_full) == 0) stop("No Fastq files found in ", opt$input)
      ## get absolute path for links
      if( !grepl("^/", opt$input)) fastq_full <- paste0(getwd(), "/", fastq_full)
   }
   ### check Modified time
   if(opt$modified != "NA"){
      mtime <- lapply(fastq_full, file.mtime)
      n <- sapply(mtime, substr, 1, 10) >= opt$modified
      if(sum(n)==0) stop("No fastq files modified after ", opt$modified)
      fastq_full <- fastq_full[n]
   }
   
   fastq <- basename(fastq_full)

   ## add parsing option?  this gets string before _ or - or . for directory name
  ids <- gsub("([^_.-]+).*", "\\1", fastq)

   for(i in seq_along(fastq)){
     if( !dir.exists(ids[i])) dir.create(ids[i])
     file.symlink( fastq_full[i], ids[i])
     cmd_txt <- paste0(ids[i], "/cmd.txt")
     ## avoid error for paired end
     if(!file.exists(cmd_txt))  file.link( "cmd.txt", cmd_txt)
   }
   message("Linked ",  length(fastq), " Fastq files")
}