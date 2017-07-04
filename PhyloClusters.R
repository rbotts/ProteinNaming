# Code finds shortest representative sequence from each cluster of sequences in a family
# Extracts the fasta AA sequence from the file and performs an alignment all of those in each family and creates a phylo tree
#
# Author:  Ryan Botts
# Date:  Sep 2016
#
# Input:
#   .fa file - all fasta AA sequences
#   ClusterGroups.csv - file with lists of all sequences in each group
#
# Output:
#   .fa - file for each protein family with centroids

library("seqinr") # used for sequence manipulations
library("msa")  # used for sequence alignments
library("ape") # used for phylogenetic tree
#library("data.table") # used for fast table subsetting

gatherSeqsByCluster <- function (projname)
{
  allSeqs <- read.fasta(file = paste0(projname,".fa"),seqtype = "AA")
  CG <<-
    read.csv(file = "ClusterGroups_6-27-2017_1.csv", head = FALSE, sep = ",")
  CG <<- na.omit(CG)
  
  CG[which(CG$V2 %in% CG[which(CG$V9==1),"V2"]),"V9"] <<- ( 1 + CG[which(CG$V2 %in% CG[which(CG$V9==1),"V2"]),"V9"] )
  
  for (name in unique(CG$V1))
  {
    rows <- CG[which(CG$V1==name),]
    
    # Criteria for which ones are to be written to the file
    # rows <- rows[which(rows$V8 == '*'),]
    rows <- rows[which(rows$V9 >= 1),]
    
    
    seqNames <- paste(sep="_",rows$V3,rows$V4,rows$V5,rows$V6)
    seqs <- allSeqs[names(allSeqs) %in% seqNames]
    write.fasta(seqs,names(seqs),file.out = paste0("Seqs_ex/",name,".faa"))
  }
}

  alignSeqs <- function() {
    #
    # function reads all sequences in a folder and computes a muscle alignment for each.
    # seqs are stored in a new folder in fasta alignment format
    
    fastanames <- list.files(path = "Seqs/", full.names = FALSE, recursive = FALSE)
    for (name in fastanames){
      print(paste("Aligning seqs for group:",name))
      seqs <- readAAStringSet(paste0("Seqs/",name),format = "fasta")
      # perform the multiple sequence alignment
      # syntax here specifies using muscle from the muscle package which allows writing a log file
      if (length(seqs)<251){
        algn <- muscle::muscle(seqs, log = paste0("SeqAlignments/",gsub(".fa","log.txt",name)), verbose = TRUE)
        # write the output to a fast file
        writeXStringSet(as(algn, "AAStringSet"), file = paste0("SeqAlignments/",name))
      } else {
        print(paste("Number of seqs=",length(seqs),"Too many seqs to align."))}
    }
    
  }

phyloTree <- function() 
{
  #
  # function reads all sequences in a folder and computes a muscle alignment for each.
  # seqs are stored in a new folder in fasta alignment format
  fastanames <- list.files(path = "./SeqAlignments/",pattern = ".faa",full.names = FALSE, recursive = FALSE)
  for (name in fastanames)
  {
    print(paste("Building tree for:",name))
    align <- read.alignment(file = paste0("./","SeqAlignments/",name),format = "fasta")
    matrix <- dist.alignment(align)
    if(length(align$seq)>2){
      write.tree(njs(matrix), file = file.path("./","SeqTreesTest/",gsub(".faa",".nwk",name)), tree.names = TRUE)
    } else print(paste(name, "has too few of clusters"))
  }
}

"Plasmids20-200kb-6-9-2016AA" -> proj

gatherSeqsByCluster(proj)
#alignSeqs()
#phyloTree()
