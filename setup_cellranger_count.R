#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

opts <- list(
  make_option(c("-e", "--email"),
    default = "brian.lohman@hci.utah.edu",
    help = "#e email or hci user name for pysano directives, required"
  ),
  make_option(c("-c", "--cluster"),
    default = "kingspeak_24",
    help = "#c cluster name for pysano with at least 128 GB of RAM, default kingspeak_24"
  ),
  make_option(c("-n", "--ncells"),
    default = "3000",
    help = "Argument to pass to 'cellranger count --expect-cells' . Depends on how many cells are expected to be recovered."
  ),
  make_option(c("-i", "--input"),
    default = "NA",
    help = "Input directory with Fastq files (either absolute path or relative path in current directory)"
  ),
  make_option(c("-f", "--fastq"),
    default = ".gz",
    help = "Match fastq files ending in default .gz"
  ),
  make_option(c("-r", "--run"),
    default = "NA",
    help = "Run ID in /Repository/MicroarryData, optional"
  ),
  make_option(c("-d", "--database"),
    default = "human",
    help = "Reference database. human (default) and mouse have preset paths but a full path to a reference in /tomato can be provided"
  ),
  make_option(c("-v", "--vdj"),
    action = "store_true", default = FALSE,
    help = "Use this flag for 'cellranger vdj' call."
  ),
  make_option(c("-a", "--atac"),
    action = "store_true", default = FALSE,
    help = "Use this flag for 'cellranger-atac count' call."
  )
)

parser <- OptionParser(option_list = opts, description = "
Soft links Cell Ranger 'mkfastq' output to current working directory in order to run Cellranger pipelines with Pysano.")
opt <- parse_args(parser)
if (opt$email == "NA") {
  print_help(parser)
  quit(status = 1)
}

# Logging
write(as.character(Sys.time()), file = "./log.txt", append = TRUE)
write(paste("email: ", opt$email), file = "./log.txt", append = TRUE)
write(paste("cluster: ", opt$cluster), file = "./log.txt", append = TRUE)
write(paste("ncells: ", opt$ncells), file = "./log.txt", append = TRUE)
write(paste("input: ", opt$input), file = "./log.txt", append = TRUE)
write(paste("fastq: ", opt$fastq), file = "./log.txt", append = TRUE)
write(paste("run: ", opt$run), file = "./log.txt", append = TRUE)
write(paste("database: ", opt$database), file = "./log.txt", append = TRUE)
if (opt$vdj) {
  write("run type = vdj", file = "./log.txt", append = TRUE)
} else if (opt$atac) {
  write("run type = atac", file = "./log.txt", append = TRUE)
} else {
  write("run type = count", file = "./log.txt", append = TRUE)
}

# Check if output already exists and prompt user will be overwritten
if (file.exists("cmd.txt")) {
  message("Note: cmd.txt file already exists")
} else {
  # Convert an HCI user name to email address
  if (!grepl("@", opt$email)) opt$email <- paste0(opt$email, "@hci.utah.edu")

  # Cast ncells as integer
  opt$ncells <- as.integer(opt$ncells)

  # Check that ncells is present
  if (is.na(opt$ncells)) {
    stop("Option -n (--ncells) must be an integer, e.g. 3000.")
  }

  # Check that cluster has memory needed
  if (opt$cluster != "kingspeak_24") {
    stop("Cellranger requires at least 128 Gb of RAM. Jobs can be run on Redwood with -c 64 slurm option")
  }

  # Select the correct cmd.txt template based on options
  # VDJ
  if (opt$vdj) {
    cmd_txt <- "/home/BioApps/hciSingleCellScripts/templates/cmd_cellranger_vdj.txt"
    if (opt$database == "human") {
      refversion <- "/tomato/dev/data/Human/GRCh38/10x_atac/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"
    } else if (opt$database == "mouse") {
      refversion <- "/tomato/dev/data/Mouse/GRCm38/10x_star/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0"
    } else {
      message(paste("Database", opt$database, "may not be supported. Check that valid database is available on desired cluster."))
      refversion <- opt$database
    }
  }
  # ATAC
  else if (opt$atac) {
    cmd_txt <- "/home/BioApps/hciSingleCellScripts/templates/cmd_cellranger-atac_count.txt"
    if (opt$database == "human") {
      refversion <- "/tomato/dev/data/Human/GRCh38/10x_atac/refdata-cellranger-atac-GRCh38-1.0.1"
    } else if (opt$database == "mouse") {
      refversion <- "/tomato/dev/data/Mouse/GRCm38/10x_star/refdata-cellranger-atac-mm10-1.2.0"
    } else {
      message(paste("Database", opt$database, "may not be supported. Check that valid database is available on desired cluster."))
      refversion <- opt$database
    }
  }
  # Standard gene expression
  else {
    cmd_txt <- "/home/BioApps/hciSingleCellScripts/templates/cmd_cellranger_count.txt"
    if (opt$database == "human") {
      refversion <- "/tomato/dev/data/Human/GRCh38/10x_star/refdata-gex-GRCh38-2020-A"
    } else if (opt$database == "mouse") {
      refversion <- "/tomato/dev/data/Mouse/GRCm38/10x_star/refdata-gex-mm10-2020-A"
    } else if (opt$database == "zebrafish") {
      refversion <- "/tomato/dev/data/Zebrafish/GRCz11/10x_star/GRCz11_v104"
    } else {
      message(paste("Database", opt$database, "may not be supported. Check that valid database is available on desired cluster."))
      refversion <- opt$database
    }
  }

  # Read in the template and change variables
  cmd <- readr::read_lines(cmd_txt, skip = 2)
  cmd <- stringr::str_replace_all(cmd, c(
    `@EMAIL` = opt$email, `@CLUSTER` = opt$cluster,
    `@REFVERSION` = refversion, `@DIR` = refversion,
    `@EXPECTED` = opt$ncells
  ))
  # Write cmd.txt to file
  readr::write_lines(cmd, "cmd.txt")
  message("Created cmd.txt file.")
}

# Setup files for Pysano
if (opt$run == "NA" & opt$input == "NA") {
  message("Note: add a run ID using -r or fastq file directory using -i to link fastq files")
} else {
  if (opt$run != "NA") {
    repo <- Sys.glob(paste0("/Repository/MicroarrayData/*/", opt$run))
    if (length(repo) == 0) stop("No match to ", opt$run, " in /Repository/MicroarrayData")
    sample_full <- list.files(paste0(repo, "/Fastq"), pattern = "*X*", full.names = TRUE)
    if (length(sample_full) == 0) stop("No sample IDs found in Repository/Fastq.")
  } else {
    sample_dirs <- list.dirs(opt$input, full.names = TRUE)[-1]
    if (length(sample_dirs) == 0) {
      stop("No sample IDs files found in ", opt$input)
    }
    ## get absolute path for links
    if (!grepl("^/", opt$input)) {
      # drop any leading "."
      if (grepl("^.", opt$input)) {
        sample_dirs <- sub("^.", "", sample_dirs)
      }
      sample_dirs <- paste0(getwd(), sample_dirs)
    }
  }

  # Check 4 fastq files per library
  nfastqs <- sapply(sample_dirs, function(d) length(list.files(path = d, pattern = "fastq.gz$")))
  if (any(nfastqs < 4)) {
    stop("Each 10x library must have 4 fastq files associated. This is not the case for at least one.")
  }
  names(sample_dirs) <- basename(sample_dirs)

  # Get the fastq files by library
  fastq_path <- lapply(sample_dirs, list.files, pattern = "fastq.gz$", full.names = T)

  # Verify I1, I2, R1, R2 files are present for each sample.
  check_types <- function(f, type) {
    any(grepl(pattern = type, x = f))
  }
  all_types <- sapply(c("I1", "I2", "R1", "R2"), function(t_) {
    all(sapply(fastq_path, check_types, type = t_))
  })
  if (!all(all_types)) {
    stop("The input FASTQ files are missing one of the required types: I1, I2, R2, R3")
  }

  # For each 10x library link the FASTQ files.
  # Produce the following directory tree
  # lib/
  #    cmd.txt
  #    MKFASTQ/I1.fastq.gz
  #    MKFASTQ/I2.fastq.gz
  #    MKFASTQ/R1.fastq.gz
  #    MKFASTQ/R2.fastq.gz

  for (i in names(fastq_path)) {
    mkfastqdir <- file.path(i, paste0("MKFASTQ_", i))
    if (!dir.exists(i)) {
      dir.create(i)
      dir.create(mkfastqdir)
    }
    for (fastq in fastq_path[[i]]) {
      linked <- file.symlink(fastq, mkfastqdir)
      if (!linked) {
        message("Warning: file ", fastq, " was not linked to dir.", i)
      } else {
        message(i, ": soft linked file ", fastq)
      }
    }
    cmd_txt <- file.path(i, "/cmd.txt")

    # link the master cmd.txt file to the library directory
    if (!file.exists(cmd_txt)) file.link("cmd.txt", cmd_txt)
  }
}
