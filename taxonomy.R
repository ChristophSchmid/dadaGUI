#!/usr/local/bin/Rscript

#Script for assigning taxonomy to sequence table
#The script assigns taxonomic units to a previously computed
#sequence table. The database used for that can be specified
#by the user.

#written by Christoph Schmid, February 2017

# CHECK ARGUMENTS PASSED AND READ INPUT FILE ---------------------------------

  # load optparse library
    library(optparse)

  # evaluate supplied arguments
    option_list = list(
      make_option(c("-i", "--input"), type = "character", default = NULL, 
                  help = "path to input .RData file"),
      make_option(c("-o", "--output"), type = "character", default = NULL,
                  help = "output path"),
      make_option(c("-d", "--database"), type = "character", default = "silva",
                  help = "database used for assignments. Must be either of silva, rdp or gg [default: %default]"),
      make_option(c("--noPS"), action = "store_true", default = FALSE,
                  help = "If set, creation of phyloseq object is turned off"),
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
    if(is.null(opt$database)) {
      print_help(opt_parser)
      stop("No database specified.", call. = TRUE)
    } else if(tolower(opt$database) %in% c("silva", "rdp", "gg", "unite")) {
      if(tolower(opt$database) == "silva") databasePath = file.path(opt$path, "taxonomy/silva")
      if(tolower(opt$database) == "rdp") databasePath = file.path(opt$path, "taxonomy/rdp")
      if(tolower(opt$database) == "gg") databasePath = file.path(opt$path, "taxonomy/gg")
      if(tolower(opt$database) == "unite") databasePath = file.path(opt$path, "taxonomy/unite")
      
      if(tolower(opt$database) %in% c("silva", "rdp", "gg")) {
        tmpFiles <- list.files(databasePath, full.names = TRUE)
        toGenus <- tmpFiles[grepl(pattern = "train_set", x = tmpFiles)]
        if(tolower(opt$database) != "gg") toSpecies <- tmpFiles[grepl(pattern = "species", x = tmpFiles)]
      } else {
        toGenus <- tmpFiles[grepl(pattern = "general_release", x = tmpFiles)]
      }
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
      # Loads input file containing outSeq and concat objects
      try(load(opt$input))
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

# ASSIGN TAXONOMY -------------------------------------------------------------

    message("Loading necessary R packages ...")
  # load necessary libraries
    if(versAvlb[versAvlb$version == opt$version,]$path == "[default]") {
      suppressPackageStartupMessages(library(dada2))
    } else {
      suppressPackageStartupMessages(library(dada2, lib.loc = versAvlb[versAvlb$version == opt$version,]$path))
    }
    suppressPackageStartupMessages(library(ShortRead))
    
  # to genus level (== to species for GG / unite databases)
  # outSeq object stems from previous script
    message(paste0("Assigning taxonomy using database file: ", toGenus))
    taxa <- assignTaxonomy(outSeq, toGenus, multithread = TRUE, tryRC = TRUE)
  # add species (skipped for GG and unite database as well as if sequences were concatenated)
    if(tolower(opt$database) %in% c("silva", "rdp") & !concat) {
      message(paste0("Adding species assignments using database file: ", toSpecies))
      taxa.plus <- addSpecies(taxa, toSpecies, verbose=TRUE)
    }
    
  # add taxonomic units as column names
    if(tolower(opt$database) %in% c("gg", "unite")) {
      colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      taxaOut <- taxa
    } else if(concat) {
      colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      taxaOut <- taxa
    } else {
      colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      taxaOut <- taxa.plus
    }
    
  # save taxonomy table to files
    save(taxaOut, file = file.path(opt$output, "taxonomyTable.RData"))
    seqTaxTable <- cbind(t(outSeq), taxaOut)
    write.table(seqTaxTable, file = file.path(opt$output, "seqTabClean_taxonomy.csv"), sep = "\t", quote = F)
    
    
# COMBINE DATA INTO PHYLOSEQ OBJECT FOR FURTHER USE -----------------------------
    
  # if --noPS option is set, phyloseq object is not created
    if(!opt$noPS)
    {
    #load phyloseq
      suppressPackageStartupMessages(library(phyloseq))
      
  # create phyloseq object with or without tree
      message("Creating phyloseq object ...")
      RSVs <- phyloseq(tax_table(taxaOut), otu_table(outSeq, taxa_are_rows = FALSE))
    
  # save taxonomy table to files
      save(RSVs, file = file.path(opt$output, "forPhyloseq.RData"))
    }
    