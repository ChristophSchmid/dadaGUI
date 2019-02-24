#!/usr/local/bin/Rscript

#Script for processing of selected input FASTQ files
  #The script will read FASTQ files and produce quality
  #plots for forward and reverse reads separately. The
  #number of produced plots can be changed using the -p 
  #argument. Produced plots will be written to a subfolder 
  #within the input path.

#written by Christoph Schmid, February 2017

# CHECK PASSED ARGUMENTS AND READ INPUT FILE ------------------------------------
    library(optparse)
  #evaluate supplied arguments
    option_list = list(
      make_option(c("-i", "--input"), type = "character", default = NULL, 
                  help = "input file with paths to FASTQs"),
      make_option(c("-o", "--output", type = "character"), default = NULL,
                  help = "output directory"),
      make_option(c("-p", "--plot"), type = "integer", default = 5,
                  help = "number of produced quality profile plots [default %default]"),
      make_option(c("-V", "--version"), type = "character", default = NULL,
                  help = "DADA2 version to be used. Unknown versions will be replaced by latest stable."),
      make_option(c("--path"), type = "character", default = NULL,
                  help = "The installation path of the pipeline.")
      )
    
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)
    
  #check if a valid installation path was provided
    if(is.null(opt$path)) {
      print_help(opt_parser)
      stop("No installation path was provided to the --path option")
    } else if(!file.exists(file.path(opt$path, "versionsDADA2.txt"))) {
      stop("The installation path was not found.")
    }
    
  #open file connection to input file and read in paths to FASTQ files
    #check supplied number of plots
    if(is.null(opt$plot)) {
      print_help(opt_parser)
      stop("Invalid number of plots", call. = TRUE)
    }
  #check input file (existing, paths readable?)
    if(is.null(opt$input)) {
      print_help(opt_parser)
      stop("Input file missing", call. = TRUE)
    }else {
      try(inputFile <- file(opt$input, open = "r"))
      fns <- readLines(inputFile)
      close(inputFile)
      fileStatus <- file.exists(fns)
      if (!all(fileStatus)) {
        stop(paste0("Input file contains unreadable file paths:\n", fns[which(fileStatus)]))
      }
    }
  #check output directory
    if(is.null(opt$output)) {
      stop("Output directory missing", call. = TRUE)
    }else {
      if(!dir.exists(opt$output)) {
        print(paste0("Output directory created: ", opt$output))
        try(dir.create(opt$output))
      }
    }
  #check dada2 version requested
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

# READ IN FASTQ-FILES (after adapters were removed) ----------------------

  # load necessary libraries
    suppressPackageStartupMessages(library(ShortRead))
    
  # check which DADA2 version has been loaded
    message(paste0("You are using dada2 version: ", getNamespaceVersion("dada2")))
    
  #forward and reverse reads per sample in separate files
  #adapters should be removed before starting with the pipeline
    
  #grep FASTQs from input file paths
    fastqs <- fns[grepl("pair[[:digit:]]", fns)]
    fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
    fnFs <- fastqs[grepl("pair1", fastqs)] # Just the forward read files
    fnRs <- fastqs[grepl("pair2", fastqs)] # Just the reverse read files
    
  # Get sample names from the first part of the forward read filenames
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    
  #write paths to FASTQs in text files for further use
    selectedFqsF <- file.path(opt$output, "selectedFilesF.txt")
    selectedFqsR <- file.path(opt$output, "selectedFilesR.txt")
    exportSel <- file(selectedFqsF, open = "w")
    writeLines(fnFs, exportSel)
    close(exportSel)
    exportSel <- file(selectedFqsR, open = "w")
    writeLines(fnRs, exportSel)
    close(exportSel)
    
# PRODUCE QUALITY PLOTS FOR CHOSEN AMOUNT OF FASTQ FILES ------------------------
    
    if(opt$plot > 0) {
      #create path for writing quality plots
      plotPath <- file.path(opt$output, "/qualityPlots/")
      if(!dir.exists(plotPath)) {
        dir.create(plotPath)       
      }
      
      for(i in ceiling(seq(from = 1, to = length(fnFs), length.out = min(length(fnFs), opt$plot)))) {
        message(paste0("Processing sample: ", sample.names[i]))
        
        #create quality profile plots
        plotF <- dada2::plotQualityProfile(fnFs[[i]])
        plotR <- dada2::plotQualityProfile(fnRs[[i]])
        
        #save plots as png files
        ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_F.png")), plot = plotF, 
                        device = "png", width = 15, height = 12, units = "cm")
        ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_R.png")), plot = plotR,
                        device = "png", width = 15, height = 12, units = "cm")
      }
    }