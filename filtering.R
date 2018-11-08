#$ -S /usr/bin/env Rscript
#$ -cwd 
#$ -e log
#$ -o log
# name for job
#$ -N filtering.R
#$ -pe hmp 10
#
# V1.00 written by Christoph Schmid, February 2017
# V1.1, September 2018
# ---------------------------

#Script for filtering of FASTQ files using dada2
  #According to user-set preferences this script
  #quality-filters selected FASTQ files and saves
  #the output in a new folder called "filtered".
  #Furthermore, dereplication is done to prepare
  #denoising.

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ---------------------------------
  
    library(optparse)

  #parse supplied arguments using optparse
    option_list = list(
      make_option(c("-f", "--forward"), type = "character", default = NULL,
                  help = "path file to forward read FASTQs"),
      make_option(c("-r", "--reverse"), type = "character", default = NULL, 
                  help = "path file to reverse read FASTQs"),
      make_option(c("-x", "--truncRfwd"), type = "integer", default = NULL,
                  help = "cut first ... bases of fwd. reads"),
      make_option(c("-y", "--truncRrev"), type = "integer", default = NULL,
                  help = "cut first ... bases of rev. reads"),
      make_option(c("--truncLfwd"), type = "integer", default = 10,
                  help = "cut forward reads at length ..."),
      make_option(c("--truncLrev"), type = "integer", default = 10,
                  help = "cut reverse reads at length ..."),
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
      make_option(c("-d", "--derep"), action = "store_true", default = FALSE,
                  help= "If set, dereplication of sequences is omitted")
#      make_option(c("-V", "--version"), type = "character", default = NULL,
#                  help = "DADA2 version to be used. Unknown versions will be replaced by latest stable.")
    )
    
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)
  
  #check and read in arguments supplied
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
      print_help(opt_parser)
      stop("Path file for reverse reads missing.", call. = TRUE)
    }else {
      try(inputFile <- file(opt$reverse, open = "r"))
      fnRs <- readLines(inputFile)
      close(inputFile)
    }
    
  #check output directory
    if(is.null(opt$output)) {
      stop("Output directory missing", call. = TRUE)
    }else {
      if(!dir.exists(opt$output)) {
        message(paste0("Output directory created: ", opt$output))
        try(dir.create(opt$output))
      }
    }
    
  #check amount of bp cut from beginning of reads
    if(opt$truncLfwd < 10) {
      message("Truncation of 5'-end of forward reads changed to 10. (recommended minimum)")
      opt$truncLfwd = 10
    }
    if(opt$truncLrev < 10) {
      message("Truncation of 5'-end of reverse reads changed to 10. (recommended minimum)")
      opt$truncLrev = 10
    }
    
  #check dada2 version requested
# VERSION MANAGEMENT CURRENTLY NOT IMPLEMENTED
#    if(file.exists("/project/genomics/Christoph/DADA2/package/versionsDADA2.txt")) {
#      versAvlb <- read.delim("/project/genomics/Christoph/DADA2/package/versionsDADA2.txt", 
#                             header = T, stringsAsFactors = F)
#    } else{
#      stop("Did not find DADA2 installation.")
#    }
#    
#    if(is.null(opt$version)) {
#      opt$version <- max(versAvlb[versAvlb$status == "stable",]$version)
#      message("No DADA2 version requested, using latest stable.: ", opt$version)
#    }else if(!opt$version %in% versAvlb$version) {
#      opt$version <- max(versAvlb[versAvlb$status == "stable",]$version)
#      message("DADA2 version requested not available, using latest stable: ", opt$version)
#    }
    
    setwd(opt$output)

# FILTERING AND TRIMMING --------------------------------------------------------

    
  #load necessary libraries
    capture.output(
      {library(grDevices)
#        library("dada2", lib.loc = versAvlb[versAvlb$version == opt$version,]$path)
        library(dada2)
        library(ShortRead)
        library(graphics)},
      type = c("message"),
      file = "/dev/null")
    
    message(paste0("This is dada2 version: ", getNamespaceVersion("dada2")))
    
    message("Quality filtering of FASTQs ...")
    message("-------------------------------")
  # Get sample names from the first part of the forward read filenames
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
  # Set file names and directory for filtered FASTQs
    filt_path <- file.path(getwd(), "filtered")
    if(!file_test("-d", filt_path)) dir.create(filt_path)
    filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
    
    # Filter according to specified settings
    filtArgs <- list(fwd = fnFs, 
                     rev = fnRs, 
                     filt = filtFs, 
                     filt.rev = filtRs,
                     trimLeft = c(opt$truncLfwd, opt$truncLrev), 
                     truncLen = c(opt$truncRfwd, opt$truncRrev), 
                     rm.phix = TRUE,
                     maxN = 0, 
                     truncQ = opt$quality, 
                     compress = opt$compress, 
                     verbose = opt$verbose,
                     multithread = TRUE)
    #add optionally supplied filter settings
    if(!is.null(opt$minLen)) filtArgs$minLen <- c(opt$minLenF, opt$minLenR)
    if(!is.null(opt$maxLen)) filtArgs$maxLen <- c(opt$maxLenF, opt$maxLenR)
    if(!is.null(opt$maxError)) filtArgs$maxEE <- opt$maxError
    
    #perform filtering and save results to file
    out <- do.call("filterAndTrim", filtArgs)
    (cbind(reads.in = out[,1], reads.out =  out[,2], proportion = out[,2] / out[,1]))
    write.table(out, file = file.path(opt$output, "filterReport.txt"), sep = "\t", quote = F)

# DEREPLICATION -----------------------------------------------------------------
    
  # quit session if dereplication is turned off by user
    if(opt$derep) quit(save = "no")
    
    message("Dereplication of FASTQs in process ...")
    message("--------------------------------------")
  # dereplicate filtered FASTQs
    derepFs <- derepFastq(filtFs, verbose=opt$verbose)
    derepRs <- derepFastq(filtRs, verbose=opt$verbose)
  # name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    names(derepRs) <- sample.names
    
    saveName <- file.path(opt$output, "derep.RData")
    save(derepFs, derepRs, sample.names, file = saveName)
