#!/usr/local/bin/Rscript

#Script for filtering of FASTQ files using dada2
  #According to user-set preferences this script
  #quality-filters selected FASTQ files and saves
  #the output in a new folder called "filtered".

#written by Christoph Schmid, February 2017

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ------------------------------------
  

    library(optparse)

  # evaluate supplied arguments
    option_list = list(
      make_option(c("-f", "--forward"), type = "character", default = NULL,
                  help = "path file to forward read FASTQs"),
      make_option(c("-r", "--reverse"), type = "character", default = NULL, 
                  help = "path file to reverse read FASTQs"),
      make_option(c("-x", "--truncRfwd"), type = "integer", default = NULL,
                  help = "cut first [x] bases of fwd. reads"),
      make_option(c("-y", "--truncRrev"), type = "integer", default = NULL,
                  help = "cut first [y] bases of rev. reads"),
      make_option(c("--truncLfwd"), type = "integer", default = 10,
                  help = "cut forward reads at length [truncLfwd]"),
      make_option(c("--truncLrev"), type = "integer", default = 10,
                  help = "cut reverse reads at length [truncLrev]"),
      make_option(c("-o", "--output"), type = "character", default = NULL,
                  help = "output directory for dereplication"),
      make_option(c("--minLenF"), type = "integer", default = NULL,
                  help = "Set minimum length of forward reads"),
      make_option(c("--minLenR"), type = "integer", default = NULL,
                  help = "Set minimum length of reverse reads"), 
      make_option(c("--maxLenF"), type = "integer", default = NULL,
                  help = "Set maximum length of forward reads"),
      make_option(c("--maxLenR"), type = "integer", default = NULL,
                  help = "Set maximum length of reverse reads"),
      make_option(c("-e", "--maxError"), type = "integer", default = NULL,
                  help = "Maximum expected errors per sequence"),
      make_option(c("-q", "--quality"), type = "integer", default = 2,
                  help = "Minimum quality of read. [default %default]"),
      make_option(c("-c", "--compress"), action = "store_false", default= TRUE,
                  help = "If set, compression of filter output is omitted"),
      make_option(c("-v", "--verbose"), action = "store_false", default = TRUE,
                  help = "If set, verbose output is turned off"),
      make_option(c("-V", "--version"), type = "character", default = NULL,
                  help = "DADA2 version to be used. Unknown versions will be replaced by latest stable."),
      make_option(c("--path"), type = "character", default = NULL,
                  help = "The installation path of the pipeline.")
    )
    
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)
    
  # check if a valid installation path was provided
    if(is.null(opt$path)) {
      print_help(opt_parser)
      stop("No installation path was provided to the --path option")
    } else if(!file.exists(file.path(opt$path, "versionsDADA2.txt"))) {
      stop("The versions file for dada2 installations was not found.")
    }
  
  # check and read in arguments supplied
    if(is.null(opt$truncRfwd) | is.null(opt$truncRrev)) {
      print_help(opt_parser)
      stop("Truncation length missing", call. = TRUE)
    }
    
    if(is.null(opt$forward)) {
      print_help(opt_parser)
      stop("Path file for forward reads missing.", call. = TRUE)
    }else {
      try(inputFile <- file(opt$forward, open = "r"))
      fnFs <- readLines(inputFile)
      close(inputFile)
    }
    
    if(is.null(opt$reverse)) {
      message("No reverse reads supplied, using forward reads only.")
    }else {
      try(inputFile <- file(opt$reverse, open = "r"))
      fnRs <- readLines(inputFile)
      close(inputFile)
    }
    
  # check output directory
    if(is.null(opt$output)) {
      stop("Output directory missing", call. = TRUE)
    }else {
      if(!dir.exists(opt$output)) {
        message(paste0("Output directory created: ", opt$output))
        try(dir.create(opt$output))
      }
    }
    
  # check dada2 version requested
    if(file.exists(file.path(opt$path, "versionsDADA2.txt"))) {
      versAvlb <- read.delim(file.path(opt$path, "versionsDADA2.txt"), 
                             header = T, stringsAsFactors = F)
    } else{
      stop("Did not find file: versionsDADA2.txt")
    }
    
    if(is.null(opt$version)) {
      opt$version <- max(numeric_version(versAvlb[versAvlb$status == "stable",]$version))
      message("No DADA2 version requested, using latest stable.: ", opt$version)
    }else if(!opt$version %in% versAvlb$version) {
      opt$version <- max(numeric_version(versAvlb[versAvlb$status == "stable",]$version))
      message("DADA2 version requested not available, using latest stable: ", opt$version)
    }

# FILTERING AND TRIMMING --------------------------------------------------------

  # load necessary libraries
    if(versAvlb[versAvlb$version == opt$version,]$path == "[default]") {
      suppressPackageStartupMessages(library(dada2))
    } else {
      suppressPackageStartupMessages(library(dada2, lib.loc = versAvlb[versAvlb$version == opt$version,]$path))
    }
    suppressPackageStartupMessages(library(ShortRead))
    
  # check which DADA2 version has been loaded
    message(paste0("You are using dada2 version: ", getNamespaceVersion("dada2")))
    
    message("Quality filtering of FASTQs ...")
    message("-------------------------------")
    
  # Get sample names from the first part of the forward read filenames
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  # Set file names and directory for filtered FASTQs
    filt_path <- file.path(opt$output, "filtered")
    if(!file_test("-d", filt_path)) dir.create(filt_path)
    
    filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
    if(!is.null(opt$reverse)) {
      filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
    }
    
    
  # Filter according to specified settings
    filtArgs <- list(fwd = fnFs, 
                     filt = filtFs, 
                     rm.phix = TRUE,
                     maxN = 0, 
                     truncQ = opt$quality, 
                     compress = opt$compress, 
                     verbose = opt$verbose,
                     multithread = TRUE)
  # add optionally supplied filter settings
    if(!is.null(opt$reverse)) {
      filtArgs$rev = fnRs
      filtArgs$filt.rev = filtRs
      filtArgs$trimLeft = c(opt$truncLfwd, opt$truncLrev)
      filtArgs$truncLen = c(opt$truncRfwd, opt$truncRrev)
      if(!is.null(opt$minLen)) filtArgs$minLen <- c(opt$minLenF, opt$minLenR)
      if(!is.null(opt$maxLen)) filtArgs$maxLen <- c(opt$maxLenF, opt$maxLenR)
    } else {
      filtArgs$trimLeft = opt$truncLfwd
      filtArgs$truncLen = opt$truncRfwd
      if(!is.null(opt$minLen)) filtArgs$minLen <- opt$minLenF
      if(!is.null(opt$maxLen)) filtArgs$maxLen <- opt$maxLenF
    }
    if(!is.null(opt$maxError)) filtArgs$maxEE <- opt$maxError
    
  # perform filtering and save results to file
    out <- do.call("filterAndTrim", filtArgs)
    (cbind(reads.in = out[,1], reads.out =  out[,2], proportion = out[,2] / out[,1]))
    
  # write report file for filtering
    write.table(out, file = file.path(opt$output, "filterReport.txt"), sep = "\t", quote = F)
