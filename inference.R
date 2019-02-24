#!/usr/local/bin/Rscript

# Script for denoising filtered and dereplicated FASTQ sequences
# The script will read in previously filtered FASTQs, before it
# does the dereplication and denoising and produces plots of the error model.
# As no further user interaction is required, merging of paired reads,
# construction a raw sequence table and removal of chimeric sequences
# is done automatically.

#written by Christoph Schmid, July 2018

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ---------------------------------

  #load optparse library
    library(optparse)

  #evaluate supplied arguments
    option_list = list(
      make_option(c("-f", "--filterpath"), type = "character", default = NULL,
                  help = "path to folder with filtered fastq files"),
      make_option(c("-p", "--plot"), type = "integer", default = 5,
                  help = "number of produced error model plots [default %default]"),
      make_option(c("-o", "--output"), type = "character", default = NULL,
                  help = "output path"),
      make_option(c("--seqtab"), action = "store_true", default = FALSE,
                  help = "If set, construction of sequence table is turned off"),
      make_option(c("--chimera"), action = "store_true", default = FALSE,
                  help = "If set, chimera removal is turned off"),
      make_option(c("--pool"), type = "integer", default = 0,
                  help = "If 0, pooling is turned off.\nPositive values: min. prevalence for priors."),
      make_option(c("--concat"), action = "store_true", default = FALSE,
                  help = "If set, forward and reverse reads will be concatenated instead of merged."),
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
      stop("The installation path was not found.")
    }

  # open file connection to input file and read inputs
  # check path folder (existing, paths readable?)
    if(is.null(opt$filterpath)) {
      print_help(opt_parser)
      stop("Path to /filtered folder missing", call. = TRUE)
    } else {
      try(inputPathFiles <- list.files(opt$filterpath, full.names = TRUE))
    }
  # check if output path was given and create directory
    if(is.null(opt$output)) {
      print_help(opt_parser)
      stop("Output path missing", call. = TRUE)
    } else {
      if(!dir.exists(opt$output)) dir.create(opt$output)
    }
  # check if plot argument was given and create directory for plots
    if(is.null(opt$plot)) {
      print_help(opt_parser)
      stop("Invalid number of plots", call. = TRUE)
    } else {
      plotPath <- file.path(opt$output, "/errorPlots/")
      if(!dir.exists(plotPath)) {
        message("Creating path for error model plots...")
        dir.create(plotPath)
      }
    }
  # check if pool argument was given and check value
    if(is.null(opt$pool)) {
      message("No pooling argument given. Assuming 0 ...")
      opt$pool <- 0
    } else {
      if(opt$pool < 0) {
        message("Pool argument < 0, changing to 0 ...")
        opt$pool <- 0
      }
    }
  # check dada2 version requested
    if(file.exists(file.path(opt$path, "versionsDADA2.txt"))) {
      versAvlb <- read.delim(file.path(opt$path, "versionsDADA2.txt"), 
                             header = T, stringsAsFactors = F)
    } else {
      stop("Did not find file: versionsDADA2.txt")
    }

    if(is.null(opt$version)) {
      opt$version <- max(numeric_version(versAvlb[versAvlb$status == "stable",]$version))
      message("No DADA2 version requested, using latest stable: ", opt$version)
    } else if(!opt$version %in% versAvlb$version) {
      opt$version <- max(numeric_version(versAvlb[versAvlb$status == "stable",]$version))
      message("DADA2 version requested not available, using latest stable: ", opt$version)
    }

# SAMPLE DENOISING FIRST ROUND --------------------------------------------------

  # load necessary libraries
    if(versAvlb[versAvlb$version == opt$version,]$path == "[default]") {
      suppressPackageStartupMessages(library(dada2))
    } else {
      suppressPackageStartupMessages(library(dada2, lib.loc = versAvlb[versAvlb$version == opt$version,]$path))
    }
    suppressPackageStartupMessages(library(ShortRead))

  # check which DADA2 version has been loaded
    message(paste0("You are using dada2 version: ", getNamespaceVersion("dada2")))

  # create filtF and filtR
    fastqs <- inputPathFiles[grepl(".fastq", inputPathFiles)]
    fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
    filtFs <- fastqs[grepl("_F_", fastqs)] # Just the forward read files
    filtRs <- fastqs[grepl("_R_", fastqs)] # Just the reverse read files

  # create indicator boolean for "forward only"
    if(length(filtRs) == 0)  {
      fwdOnly = TRUE
    } else if (length(filtRs) != length(filtFs)) {
      stop("Different number of samples for forward and reverse reads")
    } else {
      fwdOnly = FALSE
    }

  # Get sample names from the first part of the forward read filenames
  # And assign to file paths
    sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
    names(filtFs) <- sample.names
    if(!fwdOnly) names(filtRs) <- sample.names

    message("Calculating error models for sequence reads ...")

  #learn read errors from 1e8 bp / 1e6 reads
    if(numeric_version(getNamespaceVersion("dada2")) >= numeric_version("1.8.0")) {
      errF <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE)
      if(!fwdOnly) errR <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)
    } else {
      errF <- learnErrors(filtFs, nread = 1e6, multithread = TRUE, randomize = TRUE)
      if(!fwdOnly) errR <- learnErrors(filtRs, nread = 1e6, multithread = TRUE, randomize = TRUE)
    }

    if(!fwdOnly) {
      save(errF, errR, file = file.path(opt$output, "errorRates.RData"))
    } else {
      save(errF, file = file.path(opt$output, "errorRates.RData"))
    }

  # create function to get amount of sequences per sample
    getN <- function(x) sum(getUniques(x))
  # set random seed
    set.seed(42)
    
    message("Performing denoising of sequence reads ...")

  # forward reads
    derepF <- derepFastq(filtFs, verbose = TRUE)
    ddFs <- dada(derepF, err=errF, multithread = TRUE)

  # reverse reads
    if(!fwdOnly) {
      derepR <- derepFastq(filtRs, verbose = TRUE)
      ddRs <- dada(derepR, err=errR, multithread = TRUE)
    }
    
  # pack in list in case a single sample is used
    if(length(sample.names) == 1) {
      derepF <- list(derepF)
      ddFs <- list(ddFs)
      
      if(!fwdOnly) {
        derepR <- list(derepR)
        ddRs <- list(ddRs)
      }
    }

  # plot estimated error rates for a sub-sample of samples
    
    message("Plotting error models ...")
    
    for(i in seq(from = 1, to = length(ddFs), length.out = min(length(ddFs), opt$plot))) {
      
      #create selected error plots
      errPlotF <- plotErrors(ddFs[[i]], nominalQ=TRUE)
      #save plot as png files
      ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_F.png")), plot = errPlotF, 
                      device = "png", width = 15, height = 12, units = "cm")
      if (!fwdOnly) {
        errPlotR <- plotErrors(ddRs[[i]], nominalQ=TRUE)
        #save plot
        ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_R.png")), plot = errPlotR,
                        device = "png", width = 15, height = 12, units = "cm")
      }
    }

  # if pseudo-pooling is to be performed, extract sequences for pooling
    if(opt$pool != 0 & numeric_version(getNamespaceVersion("dada2")) >= numeric_version("1.8.0")) {
      #if pseudo-pooling is to be performed, seqtables are stored temporarily
      seqtabF <- makeSequenceTable(ddFs, derepF)
      rownames(seqtabF) <- sample.names
      if(!fwdOnly) {
        seqtabR <- makeSequenceTable(ddRs, derepR)
        rownames(seqtabR) <- sample.names
      }
      
  # otherwise merge samples (paired reads) or pass ddFs (fwd reads only) and create report
    } else {
      
      message("Merging reads ...")
      
      #if no pseudo-pooling performed: sequences are merged
      if(!fwdOnly) {
        mergers <- mergePairs(ddFs, derepF, ddRs, derepR, justConcatenate = opt$concat)
      } else {
        mergers <- ddFs
      }

      #report amount of sequences left
      if(!fwdOnly) {
        report <- cbind(sapply(ddFs, getN), sapply(ddRs, getN), sapply(mergers, getN))
        
      } else {
        report <- cbind(sapply(ddFs, getN), sapply(mergers, getN))
      }
      
      rownames(report) <- sample.names
    }

	# before continueing return memory to OS
  	rm(ddFs)
  	if(!fwdOnly) {
  	  rm(ddRs)
  	}
    
# PERFORM PSEUDO-POOLING --------------------------------------------------------
    
  #if --pool option is set, pseudo-pooling of samples is performed
  #This enhances the resolution of DADA2 (only available from V. 1.8.0)
    if(opt$pool != 0 & numeric_version(getNamespaceVersion("dada2")) >= numeric_version("1.8.0")) {
      
    #extract prior sequences to be used for pooling
      priorsF <- getSequences(seqtabF)[colSums(seqtabF > 0) >= opt$pool]
      if(!fwdOnly) priorsR <- getSequences(seqtabR)[colSums(seqtabR > 0) >= opt$pool]
      
    #Infer sequence variants NOW WITH PRIORS
      message("Repeating denoising with pseudo-pooled sequences ...")
          
    #repeat denoising with chosen prior sequences
      #forward reads
      ddFs <- dada(derepF, err=errF, priors = priorsF)
      
      #reverse reads
      if(!fwdOnly) {
        ddRs <- dada(derepR, err=errR, priors = priorsR)
      }
      
    #merge sequences after pseudo-pooling (paired reads) or pass ddFs (fwd reads only)
      if(!fwdOnly) {
        mergers <- mergePairs(ddFs, derepF, ddRs, derepR, justConcatenate = opt$concat)
      } else {
        mergers <- ddFs
      }

    # report amount of sequences left
      if(!fwdOnly) {
        report <- cbind(sapply(ddFs, getN), sapply(ddRs, getN), sapply(mergers, getN))
      } else {
        report <- cbind(sapply(ddFs, getN), sapply(mergers, getN))
      }

      rownames(report) <- sample.names
      

	  #free memory
  	  rm(derepF)
  	  rm(ddFs)
  	  if(!fwdOnly) {
  	    rm(derepR)
  	    rm(ddRs)
  	  }
    }
    
  #save denoising results to file in output directory
    save(mergers, file = file.path(opt$output, "mergedReads.RData"))
    
  
# CONSTRUCT SEQUENCE TABLE ------------------------------------------------------

    message("Constructing raw sequence table ...")
    
  #if -s option is set, sequence table generation is omitted
    if(!opt$seqtab) {
    #construct table from merged pairs
    seqtab <- makeSequenceTable(mergers)
    #save raw sequence table to file
    outSeq <- seqtab
    concat <- opt$concat
    save(outSeq, concat, file = file.path(opt$output, "seqTabRaw.RData"))
    }

        
# REMOVE CHIMERAS ---------------------------------------------------------------

    message("Identifying chimeric sequences ...")
    
  #if -c option is set, program is terminated here
    if(!opt$chimera | !opt$seqtab) {
      #remove chimeric sequences from sequence table
      seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose=TRUE, multithread = TRUE)
      #read fraction of non-chimeric sequences
      message(paste0("Fraction of non-chimeras is: ", sum(seqtab.nochim)/sum(seqtab)))
      #save cleaned sequence table to file
      outSeq <- seqtab.nochim
      concat <- opt$concat
      save(outSeq, concat, file = file.path(opt$output, "seqTabClean.RData"))
      write.table(t(outSeq), file = file.path(opt$output, "seqTabClean_wo_taxonomy.csv"), 
                sep = "\t", quote = F)
    }

        
# REPORT READ NUMBERS -----------------------------------------------------------
  #report number of sequences left after removing chimeras
    
    message("Writing read number summary ...")
    
    if(!fwdOnly) {
      if(exists("seqtab.nochim")) {
        report <- cbind(report, rowSums(seqtab.nochim))
        colnames(report) <- c("denoisedF", "denoisedR", "merged", "non-chimeras")
      } else {
        colnames(report) <- c("denoisedF", "denoisedR", "merged")
      }
    } else {
      if(exists("seqtab.nochim")) {
        report <- cbind(report, rowSums(seqtab.nochim))
        colnames(report) <- c("denoisedF", "merged", "non-chimeras")
      } else {
        colnames(report) <- c("denoisedF", "merged")
      }
    }

    write.table(report, file = file.path(opt$output, "readReport.txt"), sep = "\t", quote = F)
    
