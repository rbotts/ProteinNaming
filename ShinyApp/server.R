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
  
  processData <- function(bigPlot)
  {
    backbone <- input$variable
    
    # Get the cluster items relevant to the choice of backbone
    ProtData <- ClusterGroups[which(ClusterGroups$V1 == backbone), ]
    
    names(ProtData) <- c("4LetterName","cid","protName","IncGroup","Plasmid","SeqLength","Identity")
    # Generate a unique color for each cluster
    
    cols <- rainbow(length(unique(ProtData$cid)))
    
    # V6 used to assign colors to each cluster
    ProtData$V6 <- c(1:nrow(ProtData))
    
    # Assign colors to each clusters
    index = 1
    for (i in unique(ProtData$cid))
    {
      # Note that, as we would expect, all of the ones in the same cluster
      # will have the same color
      if (bigPlot)
      {
      ProtData$V6[which(ProtData$cid == i)] <- cols[index] # old: give unique color to each
      }
      else
      {
      ProtData$V6[which(ProtData$cid == i)] <-
        index # New: give a unique ID
      }
      index = index + 1
    }
    
    if (input$highlightCentroids)
    {
      ProtData$V6[which(ProtData$Identity == "*")] <- 2 * nrow(ProtData)
    }
    # Chunk to ensure correct data formatting
    # These few lines may be unnecessary
    ProtData$Identity <- (as.character(ProtData$Identity))
    ProtData$Identity[which(ProtData$Identity == "*")] <- "This is a centroid" #"100.0"
    
    ProtData$Identity[which(is.na(ProtData$Identity))] <- 0.0
    ProtData$SeqLength <- as.double(ProtData$SeqLength)
    # ProtData$Identity <- ProtData$Identity / 100.00
    
    # V8 says whether there is a reference plasmid in the *cluster*
    # V9 says whether the *specific protein* comes from a reference plasmid
    ProtData$V8 <- FALSE
    ProtData$V9 <- FALSE
    ProtData$Ex <- FALSE
    ProtData$Ex2 <- FALSE
    references <- tolower(as.character(References$V1))
    cids <- c() # Contains cluster IDs
    cids_ex <- c()
    for (cRow in c(1:nrow(ProtData)))
    {
      plasmid <- ProtData$Plasmid[cRow]
      if (tolower(plasmid) %in% references)
      {
        cids <- c(cids, ProtData$cid[cRow])
        ProtData$V9[cRow] <-
          TRUE #This *specific protein* comes from a reference
      }
      if (plasmid %in% Exemplars)
      {
        cids_ex <- c(cids_ex, ProtData$cid[cRow])
        ProtData$Ex2[cRow] <- TRUE
      }
    }
    # Note the ones that are in a *cluster* containing a reference plasmid
    ProtData[which(ProtData$cid %in% cids), "V8"] <- TRUE
    ProtData[which(ProtData$cid %in% cids_ex), "Ex"] <- TRUE
    
    if (input$refOnly == "proteins grouped with reference plasmids")
    {
      ProtData <- ProtData[which(ProtData$V8 == TRUE), ]
    }
    else if (input$refOnly == "proteins directly from reference plasmids")
    {
      ProtData <- ProtData[which(ProtData$V9 == TRUE), ]
    }
    else if (input$refOnly == "proteins grouped with well-studied plasmids")
    {
      ProtData <- ProtData[which(ProtData$Ex == TRUE), ]
    }
    else if (input$refOnly == "proteins directly from well-studied plasmids")
    {
      ProtData <- ProtData[which(ProtData$Ex2 == TRUE), ]
    }
    else if (input$refOnly == "all EXCEPT proteins grouped with reference plasmids")
    {
      ProtData <- ProtData[which(ProtData$V8 != TRUE), ]
    }
    else if (input$refOnly == "all EXCEPT proteins grouped with well-studied plasmids")
    {
      ProtData <- ProtData[which(ProtData$Ex != TRUE), ]
    }
    else if (input$refOnly == "All") {}
    
    if (input$incgroupOnly)
    {
      ProtData <- ProtData[which(ProtData$IncGroup == tolower(input$incgroup)),]
    }
    
    # Re-order from the default (AA Cluster)
    if (input$arrange == "Sequence Length")
    {
      ProtData <-
        ProtData[order(ProtData$SeqLength), ]
    }
    else if (input$arrange == "Plasmid Name")
    {
      ProtData <- ProtData[order(ProtData$Plasmid), ]
    }
    else if (input$arrange == "Incompatibility Group")
    {
      ProtData <- ProtData[order(ProtData$IncGroup), ]
    }
    
    # seqData is only used if user tries to download
    seqData <<-
      as.character(ProtData$protName)
    
    return(ProtData)
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
  
  
  
  output$refTable <- renderTable({
    names(References) <- c("Plasmid","IncGroup","Class","Order","Genus","Species","Locus")
    (References)
  })
  
  
  output$proteinPlot <- renderPlotly({
    
    ProtData <- processData(FALSE)
    
    pp <- plot_ly(
      ProtData,
      y = ~ c(1:nrow(ProtData)),
      x = ~ ProtData$SeqLength,
      type = "bar",
      orientation = "h",
      color =  ~ProtData$V6 ,
      hoverinfo = "text",
      text = paste0(
        "NCBI Label: ",
        ProtData$protName,
        "<br>",
        "IncGroup: ",
        ProtData$IncGroup,
        "<br>",
        "Plasmid: ",
        ProtData$Plasmid,
        "<br>",
        "% Identity w/ Centroid: ",
        ProtData$Identity
      )
    ) %>%
      layout(xaxis = list(title = "Sequence Length"),
             yaxis = list(title = "Protein"))
    hide_colorbar(pp)
    
  })
  
  
  
  
  
  
  
  
  output$heatMap <- renderPlotly(
    {
      ProtData <- processData(FALSE)
      
      Data <- ProtData[which(ProtData$`4LetterName`==input$variable),]
      
      Data2 <- intersect(Data$Plasmid,References$V1)
      
      subGenome <- Genome[Data2,]
      
      p <- plot_ly(type="heatmap", z = subGenome,
                   y = Data2, x = Descriptions$V1)
      
      hide_colorbar(p)
    }
  )
  
  
  
  
  
  
  
  
  output$bigPlot <- renderImage({
    
    ProtData <- processData(TRUE)
    
    # Generate the plot.
    ggplot(ProtData, aes(x = c(1:nrow(ProtData)), y = ProtData$SeqLength)) + # Notice X position is determined by V7
      geom_bar(stat = "identity",
               fill = ProtData$V6,
               alpha = ProtData$Identity) +
      scale_y_continuous() +
      scale_x_discrete() + # Remove extra space
      
      # Protein name labels
      geom_label(
        size = .90,
        y = max(ProtData$SeqLength) / 2.0,
        fontface = ifelse(ProtData$Ex, "bold", "plain"),
        label = paste(ProtData$protName,ProtData$IncGroup,ProtData$Plasmid,sep="_"),
        color = ifelse(ProtData$Ex2, "red", "black"),
        label.padding = unit(0.05, "lines"),
        label.r = unit(0.05, "lines"),
        label.size = .01,
        alpha = .5
      ) +
      
      # Bar number labels (somewhat arbitary)
      # geom_text(size = 1.5,
      #           y = -5,
      #           aes(label = rev(ProtData$V7)),
      #           color = "black") +
      coord_flip() +
      labs(y = "Sequence Length", x = "") +
      theme(axis.text.y = element_text(size = 1),
            axis.title.y = element_text(size = 2))
    
    # The weird height values are mostly calculated from trial and error
    ggsave(
      filename = "temp.png",
      height = 25 + as.integer(25 * (10 * nrow(ProtData)) / 150.0),
      width = 3.0 * 25,
      units = "mm",
      limitsize = FALSE
    )
    
    # Return the graph
    list(src = "temp.png")
  }, deleteFile = TRUE)
  
  
  
  
  
  
  
  
  
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