#Script for denoising filtered and dereplicated FASTQ sequences
  #The script will read in previously filtered and dereplicated
  #FASTQs in the form of a .RData file. The .RData file contains three objects:
  #derepF, derepR, sample.names
  #Then it does the denoising and produces plots of the error model.
  #As no further user interaction is required, merging of paired reads,
  #construction a raw sequence table and removal of chimeric sequences
  #is done automatically.

#written by Christoph Schmid, February 2017

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ---------------------------------

  #load optparse library
    library(optparse)
  
  #evaluate supplied arguments
    option_list = list(
      make_option(c("-i", "--input"), type = "character", default = NULL, 
                  help = "path to input .RData file"),
      make_option(c("-p", "--plot"), type = "integer", default = 5,
                  help = "number of produced error model plots [default %default]"),
      make_option(c("-o", "--output"), type = "character", default = NULL,
                  help = "output path"),
      make_option(c("-m", "--merge"), action = "store_true", default = FALSE,
                  help = "If set, read merging is turned off"),
      make_option(c("-s", "--seqtab"), action = "store_true", default = FALSE,
                  help = "If set, construction of sequence table is turned off"),
      make_option(c("-c", "--chimera"), action = "store_true", default = FALSE,
                  help = "If set, chimera removal is turned off")
    )
    
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)
  
  #open file connection to input file and read inputs
    if(is.null(opt$plot)) {
      print_help(opt_parser)
      stop("Invalid number of plots", call. = TRUE)
    }
    if(is.null(opt$output)) {
      print_help(opt_parser)
      stop("Output path missing", call. = TRUE)
    }else {
      if(!dir.exists(opt$output)) dir.create(opt$output)
    }
    if(is.null(opt$input)) {
      print_help(opt_parser)
      stop("Input file missing", call. = TRUE)
    }else {
      #Loads input file containing derepF, derepR and sample.names
      try(load(opt$input))
    }
    
# SAMPLE INFERENCE ---------------------------------------------------------------
    
  #load necessary libraries
    capture.output(
      {library(grDevices)
        library("dada2", lib.loc = "/project/genomics/Christoph/DADA2/package/")
        library(ShortRead)
        library(graphics)},
      type = c("message"),
      file = "/dev/null")
    
    message(paste0("This is dada2 version: ", packageVersion("dada2")))
    
  #start error model computation and denoising
  #objects derepFs, derepRs and sample.names stem from input file derep.Rdata!
    dadaFs <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread = TRUE)
    dadaRs <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread = TRUE)
  
  #save denoising results to file in output directory
    save(dadaFs, dadaRs, file = file.path(opt$output, "dada.RData"))
    
  #plot estimated error rates for a sub-sample of samples
    for(i in seq(from = 1, to = length(dadaFs), length.out = min(length(dadaFs), opt$plot))) {
      
      #create plot directory
      if(i == 1) {
        plotPath <- file.path(opt$output, "/errorPlots/")
        if(!dir.exists(plotPath)) dir.create(plotPath)
      }
      
      #create selected error plots
      errPlotF <- plotErrors(dadaFs[[i]], nominalQ=TRUE)
      errPlotR <- plotErrors(dadaRs[[i]], nominalQ=TRUE)
      
      #save plots as png files
      ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_F.png")), plot = errPlotF, 
             device = "png", width = 15, height = 12, units = "cm")
      ggplot2::ggsave(filename = file.path(paste0(plotPath, sample.names[i], "_R.png")), plot = errPlotR,
             device = "png", width = 15, height = 12, units = "cm")
    }

  
  #optional, if convergence not reached after 10 rounds:
  #check convergence manually
    dada2:::checkConvergence(dadaFs[[1]])
    dada2:::checkConvergence(dadaRs[[1]])

# MERGE PAIRED READS ------------------------------------------------------------
    
  #if -m option is set, program is terminated here
    if(opt$merge) quit(save = "no")
  #merge reads
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
    head(mergers[[1]])
  #save merged reads to file
    save(mergers, file = file.path(opt$output, "mergedReads.RData"))

# CONSTRUCT SEQUENCE TABLE ------------------------------------------------------

  #if -s option is set, program is terminated here
    if(opt$seqtab) quit(save = "no")
  #construct table from merged pairs
    seqtab <- makeSequenceTable(mergers)
  #save raw sequence table to file
    outSeq <- seqtab
    save(outSeq, file = file.path(opt$output, "seqTabRaw.RData"))
    # write.table(seqtab, file = file.path(opt$output, "seqTabRaw_wo_taxonomy.csv"), sep = "\t", quote = F)

# REMOVE CHIMERAS ---------------------------------------------------------------

  #if -c option is set, program is terminated here
    if(opt$chimera) quit(save = "no")
  #remove chimeric sequences from sequence table
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", verbose=TRUE, multithread = TRUE)
  
  #read fraction of non-chimeric sequences
    print("Fraction of non-chimeras is:")
    print(sum(seqtab.nochim)/sum(seqtab))
  #save cleaned sequence table to file
    outSeq <- seqtab.nochim
    save(outSeq, file = file.path(opt$output, "seqTabClean.RData"))
    write.table(t(outSeq), file = file.path(opt$output, "seqTabClean_wo_taxonomy.csv"), sep = "\t", quote = F)
