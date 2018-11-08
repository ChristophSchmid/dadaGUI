#Script for calculating a taxonomic trees from a set of input sequences
#The script calculates a phylogenetic tree ( i.e. a generalized time-reversible 
#with Gamma rate variation maximum likelihood tree) from a given input
#of sequences. This is meant for the creation of a tree after filtering the ASVs
#by the user.

#written by Christoph Schmid, July 2018

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ---------------------------------

  #load optparse library
    library(optparse)
  
  #evaluate supplied arguments
    option_list = list(
      make_option(c("-i", "--input"), type = "character", default = NULL, 
                  help = "FASTA file with input sequences"),
      make_option(c("-o", "--output"), type = "character", default = NULL,
                  help = "output path")
    )
    
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)
  
  #check output directory and create if necessary
    if(is.null(opt$output)) {
      print_help(opt_parser)
      stop("Output path missing", call. = TRUE)
    }else {
      if(!dir.exists(opt$output)) dir.create(opt$output)
    }
  #open file connection to input file and read inputs
    library(Biostrings, quietly = T)
    
    if(is.null(opt$input)) {
      print_help(opt_parser)
      stop("Input file missing", call. = TRUE)
    }else {
      #Loads input FASTA file
      try(seqs <- readDNAStringSet(opt$input))
    }

# CONSTRUCT PHYLOGENETIC TREE ---------------------------------------------------

  #using package DECIPHER to align sequences
    library(DECIPHER, quietly = T, lib.loc = "/project/genomics/Christoph/DADA2/package/")
    
    alignment <- AlignSeqs(seqs, anchor = NA)
    writeXStringSet(alignment, file = file.path(opt$output, "sequenceAlignment.fasta"))
  
  #using package phangorn to construct phylogenetic starting tree
    library(phangorn, quietly = T)
    
    phang.align <- phyDat(as.matrix(alignment, use.names = TRUE), type = "DNA")
    dm <- dist.ml(phang.align)
  
  #neighbour joining tree as starting point
    treeNJ <- NJ(dm) #tip order != sequence order !!
    # fit <- pml(treeNJ, data = phang.align)
    
  #remove negative edges
    if (any(treeNJ$edge.length < 0)) {
      treeNJ$edge.length[treeNJ$edge.length < 0] <- 1e-08
      message("negative edges length changed to 0!")
    }
  
  # #fit a Generalized time-reversible with Gamma rate variation maximum likelihood tree
  #   fitGTR <- update(fit, k = 4, inv = 0.2)
  #   fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
  #                       rearrangement = "stochastic", control = pml.control(trace = 0))
    
    write.tree(treeNJ, file = file.path(opt$output, "startingTreeNJ.tre"))