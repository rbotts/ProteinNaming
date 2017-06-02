library(shiny)
library(datasets)
library(ggplot2)
library(seqinr)
library(ggtree)
library(msa)  # used for sequence alignments
# msa will only output a nice alignment as a pdf
# png and animation packages allow us to convert to a png and then open as a plot
library(png)
library(animation)
library(ape)  # used for phylogenetic trees

shinyServer(function(input, output) {
  output$description <- renderUI({
    altNames <-
      as.character(Descriptions$V3[which(Descriptions$V1 == input$variable)])
    
    if (grepl("\\^", altNames))
      altNames <- gsub("\\^", ", ", altNames)
    else
      altNames <- altNames
    
    tagList(h1(input$variable, color = "blue"),
            h4(Descriptions$V2[which(Descriptions$V1 == input$variable)]),
            HTML(
              paste0("<b>Alternative Names: </b>", "<i>", altNames, "</i>")
            ))
  })
  output$proteinPlot <- renderPlotly({
    backbone <- input$variable
    
    # Get the cluster items relevant to the choice of backbone
    CG2 <- ClusterGroups[which(ClusterGroups$V1 == backbone),]
    
    # Generate a unique color for each cluster
    cols <- rainbow(length(unique(CG2$V2)))
    
    # V6 used to assign colors to each cluster
    CG2$V6 <- c(1:nrow(CG2))
    
    # Assign colors to each clusters
    index = 1
    for (i in unique(CG2$V2))
    {
      # Note that, as we would expect, all of the ones in the same cluster
      # will have the same color
      # CG2$V6[which(CG2$V2 == i)] <- cols[index] # old: give unique color to each
      CG2$V6[which(CG2$V2 == i)] <- index # New: give a unique ID
      index = index + 1
    }
    
    # Chunk to ensure correct data formatting
    # These few lines may be unnecessary
    CG2$V5[which(CG2$V5 == "*")] <- "100.0" #"100.0"
    CG2$V5 <- as.double(as.character(CG2$V5))
    CG2$V5[which(is.na(CG2$V5))] <- 0.0
    CG2$V4 <- as.double(CG2$V4)
    CG2$V5 <- CG2$V5 / 100.00
    
    # Order by either sequence length or cluster
    if (input$arrange == "Sequence Length")
    {
      CG3 <-
        CG2[order(CG2$V4),]
    }
    else
    {
      CG3 <- CG2 # This case will "order by cluster."
    }
    
    # V7 is used for labeling
    CG3$V7 <-
      c(1:nrow(CG3))
    
    # V8 says whether there is a reference plasmid in the *cluster*
    # V9 says whether the *specific protein* comes from a reference plasmid
    CG3$V8 <- FALSE
    CG3$V9 <- FALSE
    references <- tolower(as.character(References$V1))
    cids <- c() # Contains cluster IDs
    for (cRow in c(1:nrow(CG3)))
    {
      plasmid <- strsplit(as.character(CG3$V3[cRow]), "_")[[1]]
      if (length(plasmid) > 2)
      {
        plasmid <- plasmid[3:length(plasmid)]
        plasmid <- paste(plasmid, collapse = "_")
      }
      else
      {
        plasmid <- plasmid[2]
      }
      if (tolower(plasmid) %in% references)
      {
        cids <- c(cids, CG3$V2[cRow])
        CG3$V9[cRow] <-
          TRUE #This *specific protein* comes from a reference
      }
    }
    # Note the ones that are in a *cluster* containing a reference plasmid
    CG3[which(CG3$V2 %in% cids), "V8"] <- TRUE
    
    # Take a subset (method from global.R)
    if (input$Subset)
    {
      CG3 <- subsetData(CG3, input$ranges)
    }
    
    # seqData is only used if user tries to download
    seqData <<-
      as.character(CG3$V3)
    
    pp <- plot_ly(
      CG3,
      y = ~ c(1:nrow(CG3)),
      x = ~ CG3$V4,
      type = "bar",
      orientation = "h",
      color = ~ CG3$V6,
      hoverinfo = "text",
      text = ifelse(CG3$V8, paste0(CG3$V3, ""), paste0("<b>", CG3$V3, "</b>"))
    ) %>%
      layout(xaxis = list(title = "Sequence Length"),
             yaxis = list(title = "Protein"))
    hide_colorbar(pp)
    
  })
  
  output$ptrees <- renderPlot({
    # read phylogenetic trees from previously computed newick tree
    # this code may fail if the family only has one cluster
    
    fpath <-
      file.path("SeqTreesTest", paste0(input$variable, ".nwk"))
    validate(need(
      file.exists(fpath),
      paste(
        input$variable,
        " tree unavailable due to memory constraints."
      )
    ))
    tree <- read.tree(fpath)
    ggtree(tree) + geom_tiplab(color = "purple",
                               hjust = 1.0,
                               vjust = -0.75) +
      geom_nodepoint(color = "#b5e521",
                     alpha = 3 / 8,
                     size = 8) +
      geom_tippoint(color = "purple",
                    shape = 20,
                    size = 4)
    
  })
  
  # This is supposed to be the download feature.
  # Offline, it gives the correct file output, but also throws an error.
  # Online, it just crashes the shinyapp server. (5/23/2017) -ZL
  output$downloadData <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste(input$variable, "faa", sep = ".")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      print("READING FASTA... (MAY TAKE UP TO 5 MINUTES)")
      plasmidAA <-
        read.fasta("Plasmids20-200kb-6-9-2016AA.fa", "AA")
      print("FINISHED FASTA READ")
      
      writeSeqs <- function(seq)
      {
        seqName <- attr(seq, "name")
        if (grepl("<unknown description>", seqName))
        {
          seqName <- strsplit(seqName, " ")[[1]][1]
        }
        if (seqName %in% seqData)
        {
          print("RETURNING...")
          return(seq)
          print("DONE WITH A RETURN")
        }
      }
      Seqs <<- lapply(plasmidAA, writeSeqs)
      Seqs <- Seqs[!sapply(Seqs, is.null)] # REMOVE NULL seqs
      # (since lapply return same length)
      print("DONE")
      print(input$variable)
      fourName <- input$variable
      for (i in c(1:length(Seqs)))
      {
        seqName <- attr(Seqs[i], "name")
        print("1")
        if (length(seqName) > 0)
        {
          if (grepl("<unknown description>", seqName))
          {
            print("2")
            seqName <- strsplit(seqName, " ")[[1]][1]
            print("3")
          }
        }
        if (i != 1)
        {
          sNames <- c(sNames, paste(fourName, seqName, sep = "_"))
        }
        else
        {
          sNames <- c(paste(fourName, seqName, sep = "_"))
        }
      }
      print(sNames)
      print(Seqs)
      write.fasta(Seqs, sNames,
                  file.out = paste(fourName, "faa", sep = "."))
    }
  )
  
  output$Align <- downloadHandler(filename <- function() {
    return(paste0(input$variable, ".faa"))
  },
  
  content <- function(file) {
    seqAlign <-
      read.fasta(paste0("SeqAlignments/", input$variable, ".fa"), "AA")
    write.fasta(seqAlign, names = attr(seqAlign, "name"), file)
  })
})