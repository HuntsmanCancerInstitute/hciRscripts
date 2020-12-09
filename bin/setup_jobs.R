#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

opts <-  list(
   make_option(c("-e", "--email"), default="NA",
         help="#e email or hci user name for pysano directive, required"),
   make_option(c("-c", "--cluster"), default="kingspeak",
         help="#c cluster name for pysano with at least 64 GB of RAM, default kingspeak"),
   make_option(c("-s", "--sequencing"), default="novaseq",
         help="single, paired, novaseq, qiagen, NEB or metagenome.  Single and paired options are for HiSeq runs"),
   make_option(c("-i", "--input"), default="NA",
         help="Input directory with Fastq files (either absolute path or relative path in current directory)"),
   make_option(c("-f", "--fastq"), default=".gz",
         help="Match fastq files ending in default .gz"),
   make_option(c("-m", "--modified"), default="NA",
         help="Only link fastq files with modified time >= YYYY-MM-DD"),
   make_option(c("-p", "--adapters"), default="Illumina",
         help="Adapter sequences, default TruSeq or Nextera for new RiboZero kits"),
   make_option(c("-r", "--run"), default="NA",
         help="Run ID in /Repository/MicroarryData, optional"),
   make_option(c("-a", "--analysis"), default="",
         help="Save files to /Repository/AnalysisData, optional"),
   make_option(c("-v", "--version"), default="102",
         help="Ensembl release, default 102, version 92, 94, 96, 98 and 100 are also available for human and mouse"),
   make_option(c("-d", "--database"), default="human",
         help="Reference database, default human or mouse, elephant, fly, worm, pig,
     rat, rabbit, sheep, vervet, yeast, zebrafish"),
    make_option(c("-l", "--length"), default="50",
       help="Read length for STAR reference, default 50 or 125")
)

parser <- OptionParser(option_list=opts, description = "
Creates a cmd.txt file and sample directories with Fastq file links in order to run STAR,
featureCounts and quality metrics on the CHPC clusters using pysano.  The default is to align
paied-end 50 bp reads from NovaSeq runs to the latest human reference. Check /tomato/dev/data to see if older
references are available")
  opt <- parse_args(parser)
  if( opt$email == "NA" ){
     print_help(parser)
     quit(status=1)
  }

# logging
write(as.character(Sys.time()), file = "./log.txt", append = TRUE)
write(paste("email: ", opt$email), file = "./log.txt", append = TRUE)
write(paste("cluster: ", opt$cluster), file = "./log.txt", append = TRUE)
write(paste("sequencing: ", opt$sequencing), file = "./log.txt", append = TRUE)
write(paste("input: ", opt$input), file = "./log.txt", append = TRUE)
write(paste("fastq: ", opt$fastq), file = "./log.txt", append = TRUE)
write(paste("modified: ", opt$modified), file = "./log.txt", append = TRUE)
write(paste("run: ", opt$run), file = "./log.txt", append = TRUE)
write(paste("analysis: ", opt$analysis), file = "./log.txt", append = TRUE)
write(paste("version: ", opt$version), file = "./log.txt", append = TRUE)
write(paste("database: ", opt$database), file = "./log.txt", append = TRUE)
write(paste("length: ", opt$length), file = "./log.txt", append = TRUE)

if(file.exists( "cmd.txt")){
   message("Note: cmd.txt file already exists")
}else{
   if( !grepl("@", opt$email )) opt$email <- paste0( opt$email, "@hci.utah.edu")
   ## STAR version should match version used to create index
   release <- as.numeric(opt$version)
   if(!release %in% c(100, 102)) message("Note: Version may not have a reference, please check /tomato/dev/data")
   STAR_version <- "2.7.6a"
   if(release == 100) STAR_version <- "2.7.3a"
   if(release == 98)  STAR_version <- "2.7.2c"
   if(release == 96)  STAR_version <- "2.7.0f"
   if(release == 94)  STAR_version <- "2.6.1b"
   if(release == 92)  STAR_version <- "2.5.4a"
   if(release == 90)  STAR_version <- "2.5.2b"

   if( opt$analysis != "" ) opt$analysis <- paste0("#a ", opt$analysis)
   if(!opt$length %in% c("50", "125")) message("Length should be 50 or 125.  Please check if star", opt$length, " exists")
   if(!opt$sequencing %in% c("single", "paired", "qiagen", "clumpify", "novaseq",
              "NEB", "metagenome", "metatranscriptome", "microbe", "screen")){
       stop("Sequencing should be single, paired, novaseq, qiagen, NEB, microbe, metagenome or metatranscriptome")
   }
   if(opt$sequencing == "novaseq") opt$sequencing <- "clumpify"
   x <- read.delim("/home/BioApps/hciR/STAR_ref_dbs.txt", stringsAsFactors=FALSE)
   refdb <- stringr::str_replace_all(opt$database, c(fly = "Drosophila", fruitfly = "Drosophila",
                  worm = "C_elegans", yeast = "S_cerevisiae", Homo ="Human", Mus = "Mouse"))

   n <- grep(refdb, x$dir, ignore=TRUE)
   if(length(n)!=1) stop("Database should match human, mouse, elephant, fly, worm, pig, rat, rabbit, sheep, vervet, yeast or zebrafish.")
   species <- x$species[n]
   assembly <- x$assembly[n]
   tomato_dir <- x$dir[n]
   ## old assemblies
   if(assembly == "GRCz11" & release == 90) assembly <- "GRCz10"
   if(assembly == "BDGP6.28" & release == 98) assembly <- "BDGP6.22"
   if(tolower(opt$adapters) == "nextera"){
	   adapt1 <- "CTGTCTCTTATACACATCT"
	   adapt2 <- "CTGTCTCTTATACACATCT"
   }else{
	   # Illumina TruSeq adapters
	   adapt1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
	   adapt2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
   }

   # set fastq screen conf file based on host cluster
   if(opt$cluster == "redwood"){
       fastq_screen_conf = "/tomato/dev/data/FastQ_Screen_Genomes/redwood_fastq_screen.conf"
   } else {
       fastq_screen_conf = "/tomato/dev/data/FastQ_Screen_Genomes/eukrRNA/fastq_screen_rRNA.conf"
   }

   cmd_txt <- paste0("/home/BioApps/hciR/templates/cmd_", opt$sequencing, ".txt")
   cmd <- readr::read_lines(cmd_txt, skip=2)
   cmd <- stringr::str_replace_all(cmd, c(`@EMAIL` = opt$email, `@CLUSTER` = opt$cluster,
    `@ANALYSIS` = opt$analysis, `@LENGTH` =  opt$length, `@DIR` = tomato_dir, `@SPECIES`= species,
    `@ASSEMBLY`= assembly ,`@VERSION` = release, `@STAR` = STAR_version, `@FASTQ` = opt$fastq, `@SCREEN_CONF` = fastq_screen_conf,
     `@ADAPT1`= adapt1, `@ADAPT2`= adapt2 ))
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
   # ids <- gsub("([^_.-]+).*", "\\1", fastq)
   ids <- gsub("([^_.]+).*", "\\1", fastq)

   for(i in seq_along(fastq)){
     if( !dir.exists(ids[i])) dir.create(ids[i])
     file.symlink( fastq_full[i], ids[i])
     cmd_txt <- paste0(ids[i], "/cmd.txt")
     ## avoid error for paired end
     if(!file.exists(cmd_txt))  file.link( "cmd.txt", cmd_txt)
   }
   message("Linked ",  length(fastq), " Fastq files")
}
