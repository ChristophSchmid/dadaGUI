#$ -S /usr/bin/env Rscript
#$ -cwd 
# set log file(s) name
#$ -e log
#$ -o log
# name for job
#$ -N taxonomy.R
#$ -pe hmp 10
#
# V1.00 written by Christoph Schmid, February 2017
# V1.1, September 2018
# ---------------------------

#Script for assigning taxonomy to a sequence table
  #The script assigns taxonomic units to a previously computed
  #sequence table. The database used for that can be specified
  #by the user.

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ---------------------------------

#load optparse library
    library(optparse)
  
  #evaluate supplied arguments
    option_list = list(
      make_option(c("-i", "--input"), type = "character", default = NULL, 
                  help = "path to input .RData file"),
      make_option(c("-o", "--output"), type = "character", default = NULL,
                  help = "output path"),
      make_option(c("-d", "--database"), type = "character", default = "silva",
                  help = "database used for assignments. Must be either of silva, rdp or gg [default: %default]"),
      make_option(c("--tree"), action = "store_true", default = FALSE,
                  help = "If set, creation of phylogenetic tree is turned on"),
      make_option(c("--noPS"), action = "store_true", default = FALSE,
                  help = "If set, creation of phyloseq object is turned off")
#      make_option(c("-V", "--version"), type = "character", default = NULL,
#                  help = "DADA2 version to be used. Unknown versions will be replaced by latest stable.")
    )
    
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)

  #open file connection to input file and read inputs
    if(is.null(opt$database)) {
      print_help(opt_parser)
      stop("No database specified.", call. = TRUE)
    } else if(is.element(tolower(opt$database), c("silva", "rdp", "gg", "unite"))) {
      if(tolower(opt$database) == "silva") databasePath = file.path("/project/genomics/Christoph/DADA2/taxonomy/silva")
      if(tolower(opt$database) == "rdp") databasePath = file.path("/project/genomics/Christoph/DADA2/taxonomy/rdp")
      if(tolower(opt$database) == "gg") databasePath = file.path("/project/genomics/Christoph/DADA2/taxonomy/gg")
      if(tolower(opt$database) == "unite") databasePath = file.path("/project/genomics/Christoph/DADA2/taxonomy/unite")
    } else {
      print_help(opt_parser)
      stop("Unknown database.", call. = TRUE)
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
      #Loads input file containing outSeq and concat objects
      try(load(opt$input))
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

# ASSIGN TAXONOMY -------------------------------------------------------------
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
    
  #to genus level (== to species for GG / unite databases)
  #outSeq object stems from previous script
    taxa <- assignTaxonomy(outSeq, file.path(databasePath, "toGenus.fa.gz"), 
                           multithread = TRUE, tryRC = TRUE)
  #add species (skipped for GG and unite database as well as if sequences were concatenated)
    if((tolower(opt$database) != "gg") & (tolower(opt$database) != "unite") & !concat) {
      taxa.plus <- addSpecies(taxa, file.path(databasePath, "toSpecies.fa.gz"), verbose=TRUE)
    }
    
  #add taxonomic units as column names
    if((tolower(opt$database) == "gg") | (tolower(opt$database) == "unite")) {
      colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      taxaOut <- taxa
    } else if(concat) {
      colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      taxaOut <- taxa
    } else {
      colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      taxaOut <- taxa.plus
    }
    
  #save taxonomy table to files
    save(taxaOut, file = file.path(opt$output, "taxonomyTable.RData"))
    seqTaxTable <- cbind(t(outSeq), taxaOut)
    write.table(seqTaxTable, file = file.path(opt$output, "seqTabClean_taxonomy.csv"), sep = "\t", quote = F)
    
# CONSTRUCT PHYLOGENETIC TREE ---------------------------------------------------
  
  #if --tree option is set, phylogenetic tree is not created
    if(opt$tree) 
    {
    #using package DECIPHER to align sequences
      library(DECIPHER, quietly = T)
      
      seqs <- getSequences(outSeq)
      names(seqs) <- seqs #necessary for tip labels of tree
      alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)
    
    #using package phangorn to construct phylogenetic tree
      library(phangorn, quietly = T)
      
      phang.align <- phyDat(as.matrix(alignment, use.names = TRUE), type = "DNA")
      dm <- dist.ml(phang.align)
      
    #neighbour joining tree as starting point
      treeNJ <- NJ(dm) #tip order != sequence order !!
      fit <- pml(treeNJ, data = phang.align)
    
    #fit a Generalized time-reversible with Gamma rate variation maximum likelihood tree
      fitGTR <- update(fit, k = 4, inv = 0.2)
      fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                          rearrangement = "stochastic", control = pml.control(trace = 0))
      
      write.tree(fitGTR$tree, file = file.path(opt$output, "GTR_phylotree.tre"))
    }
    
# COMBINE DATA INTO PHYLOSEQ OBJECT FOR FURTHER USE -----------------------------
    
    #if --noPS option is set, phyloseq object is not created
    if(!opt$noPS)
    {
    #load phyloseq
      library(phyloseq)
      
    #create phyloseq object with or without tree
      if(!opt$tree) 
        {
        RSVs <- phyloseq(tax_table(taxaOut), otu_table(outSeq, taxa_are_rows = FALSE))
      }else(
        RSVs <- phyloseq(tax_table(taxaOut), otu_table(outSeq, taxa_are_rows = FALSE),
                       phy_tree(fitGTR$tree))
      )
    
    #save taxonomy table to files
      save(RSVs, file = file.path(opt$output, "forPhyloseq.RData"))
    }
    
