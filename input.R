#$ -S /usr/bin/env Rscript
#$ -cwd 
# set log file name(s)
#$ -e log
#$ -o log
# name for job
#$ -N input.R
#$ -pe hmp 10
#
# V1.00 written by Christoph Schmid, February 2017
# V1.1, September 2018
# ---------------------------

#Script for processing of selected input FASTQ files
  #The script will read FASTQ files and produce quality
  #plots for forward and reverse reads separately. The
  #number of produced plots can be changed using the -p 
  #argument. Produced plots will be written to a subfolder 
  #within the input path.

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ---------------------------------
    library(optparse)
  #evaluate supplied arguments
    option_list = list(
      make_option(c("-i", "--input"), type = "character", default = NULL, 
                  help = "input file with paths to FASTQs"),
      make_option(c("-o", "--output", type = "character"), default = NULL,
                  help = "output directory"),
      make_option(c("-p", "--plot"), type = "integer", default = 5,
                  help = "number of produced quality profile plots [default %default]")
#      make_option(c("-V", "--version"), type = "character", default = NULL,
#                  help = "DADA2 version to be used. Unknown versions will be replaced by latest stable.")
      )
    
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)
    
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
#VERSION MANAGEMENT CURRENTLY NOT IMPLEMENTED
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

# READ IN FASTQ-FILES (after adapters are removed) -----------------------

  #load necessary libraries
    capture.output(
      {library(grDevices)
#        library("dada2", lib.loc = versAvlb[versAvlb$version == opt$version,]$path)
        library(dada2)
        library(ShortRead)},
      type = c("message"),
      file = "/dev/null")
    
    message(paste0("This is dada2 version: ", getNamespaceVersion("dada2")))
    
  #forward and reverse reads per sample in separate files
  #adapters should be removed before starting with the pipeline

  #set WD to output path
    pathWD <- opt$output
    setwd(pathWD)
    
  #grab FASTQs from input file paths
    fastqs <- fns[grepl("pair[[:digit:]].truncated", fns)]
    fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
    fnFs <- fastqs[grepl("pair1", fastqs)] # Just the forward read files
    fnRs <- fastqs[grepl("pair2", fastqs)] # Just the reverse read files
    
  # Get sample names from the first part of the forward read filenames
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    
  #write paths to FASTQs in text files for further use
    selectedFqsF <- file.path(pathWD, "selectedFilesF.txt")
    selectedFqsR <- file.path(pathWD, "selectedFilesR.txt")
    exportSel <- file(selectedFqsF, open = "w")
    writeLines(fnFs, exportSel)
    close(exportSel)
    exportSel <- file(selectedFqsR, open = "w")
    writeLines(fnRs, exportSel)
    close(exportSel)
    
# PRODUCE QUALITY PLOTS FOR SELECTED FASTQ FILES --------------------------------------------
    
    if(opt$plot > 0) {
      #create path for writing quality plots
      plotPath <- file.path(pathWD, "/qualityPlots/")
      if(!dir.exists(plotPath)) {
        dir.create(plotPath)       
      }
      
      for(i in ceiling(seq(from = 1, to = length(fnFs), length.out = min(length(fnFs), opt$plot)))) {
        message(paste0("Processing sample: ", sample.names[i]))
        
        #create quality profile plots
        plotF <- plotQualityProfile(fnFs[[i]])
        plotR <- plotQualityProfile(fnRs[[i]])
        
        #save plots as png files
        ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_F.png")), plot = plotF, 
                        device = "png", width = 15, height = 12, units = "cm")
        ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_R.png")), plot = plotR,
                        device = "png", width = 15, height = 12, units = "cm")
      }
    }
