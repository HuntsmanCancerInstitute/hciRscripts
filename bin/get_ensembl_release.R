#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

opts <-  list(
   make_option(c("-d", "--dir"), default="all",
         help="Directory name in /tomato/dev/data.  Default is all genomes"),
   make_option(c("-r", "--release"),
          help="Ensembl release number"),
    make_option(c("-v", "--version"), default="2.6.1b",
          help="STAR version number")
)


parser <- OptionParser(option_list=opts, description = "
Downloads Ensembl GTF files for genomes listed in STAR_ref_dbs.txt,
creates refflat and ribosome interval files for collectRNAseq metrics
and cmd.txt files to build new STAR reference databases.
")
opt <- parse_args(parser)

if( is.null(opt$release) ){
   print_help(parser)
   quit(status=1)
}

release <- opt$release
ftp <- paste0("ftp://ftp.ensembl.org/pub/release-", release, "/gtf")
x <- read.delim("/home/BioApps/hciR/STAR_ref_dbs.txt", stringsAsFactors=FALSE)

if(opt$dir == "all"){
   genome_rows <- 1:nrow(x)
}else{
   i <- which(x$dir == opt$dir)
   if(length(i)==0) stop("No match to ", opt$dir, " in STAR_ref_dbs.txt")
   genome_rows <- i
}
for(i in genome_rows){
   message(x$dir[i], " release ", release)
   # /tomato/dev/data/Mouse/GRCm38/release90
   datadir <- paste0("/tomato/dev/data/", x$dir[i], "/", x$assembly[i], "/release", release)
   if(file.exists( datadir)){
      message(x$dir[i], " release", release, " directory already exists")
   }else{
      system(paste("mkdir", datadir))
      system(paste("chmod 777", datadir))

      # DOWNLOAD gtf
      gtffile <- paste0(x$species[i], ".", x$assembly[i], ".", release, ".gtf.gz")
      # ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz
      ftpdir <- paste0(ftp, "/", tolower(x$species[i]), "/",  gtffile)
      message( " Downloading GTF file")
      setwd(datadir)
      download.file(ftpdir, gtffile,  quiet = TRUE)
      system(paste("gunzip",  gtffile))
      ## CREATE refflat
      gtf <- gsub(".gz$", "", gtffile)
      org <- gsub(".gtf$", "", gtf)
      message(" Creating RefFlat file")
      system(paste0("/home/BioApps/UCSCExe/gtfToGenePred -genePredExt ", gtf, " tmp.out" ))
      ## paste column 12 and 1-10,
      system("cut -f1-10 tmp.out > tmp10.out")
      system("cut -f12 tmp.out > tmp12.out")
      system(paste0("paste tmp12.out tmp10.out > ", org, ".refflat"))
      file.remove("tmp.out", "tmp10.out", "tmp12.out")
      message(" Getting Ribosome intervals")
      system(paste("grep 'gene_biotype \"[Mt_]*rRNA\"'", gtf, "| awk '$3 == \"transcript\"' | cut -f1,4,5,7,9 > tmp.out"))
      system("perl -lane '/transcript_id \"([^\"]+)\"/; print join \"\t\", (@F[0,1,2], $1, 1, @F[3])' tmp.out | sort -k1V -k2n -k3n > rRNA.bed")
      dict <- list.files("..", pattern = "dict$")
      if(length(dict) != 1) stop("Sequence dictionary file not found")
      system(paste0("java -jar /tomato/dev/app/picard/2.9.0/picard.jar BedToIntervalList quiet=TRUE I=rRNA.bed SD=../", dict, " O=", org, ".rRNA.interval"))
      file.remove("tmp.out")
      ## CMD.txt files
      message(" Creating cmd.txt file")
      # link fasta for pysano
      fa <- list.files("..", pattern = "fa$")
      file.symlink(paste0("../",fa), fa)
      system("mkdir star50")
      system("mkdir star125")


      cat(paste0("#e chris.stubben@hci.utah.edu
#c kingspeak_24

STAR=/tomato/dev/app/STAR/", opt$version, "/STAR

$STAR   --runMode genomeGenerate \\
        --genomeDir star50 \\
        --genomeFastaFiles ", fa, " \\
        --runThreadN $NCPU \\
        --sjdbGTFfile ", gtf, " \\
        --sjdbOverhang 49

$STAR   --runMode genomeGenerate \\
       --genomeDir star125 \\
       --genomeFastaFiles ", fa, " \\
       --runThreadN $NCPU \\
       --sjdbGTFfile ", gtf, " \\
       --sjdbOverhang 124
"), file="cmd.txt")

   ## pstart?
   }
}
