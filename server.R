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
  
  processData <- function(bigPlot,forceSub=FALSE)
  {
    backbone <- input$variable
    
    # Get the cluster items relevant to the choice of backbone
    ProtData <- CG[which(CG$V1 == backbone), ]
    
    names(ProtData) <- c("4LetterName","cid","protName","IncGroup","Plasmid","uid","SeqLength","Identity","Ex")
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
    
    if (input$refOnly == "proteins grouped with well-studied plasmids")
    {
      ProtData <- ProtData[which(ProtData$Ex != "NONE"), ]
    }
    else if (input$refOnly == "proteins directly from well-studied plasmids")
    {
      ProtData <- ProtData[which(ProtData$Ex != "NONE"), ]
      ProtData <- ProtData[which( !grepl("~",ProtData$Ex) ), ]
    }
    else if (input$refOnly == "all EXCEPT proteins grouped with well-studied plasmids")
    {
      ProtData <- ProtData[which(ProtData$Ex == "NONE"), ]
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
    
    tagList(
            strong(style="color:blue;font-size:30px",input$variable),
            HTML("&nbsp&nbsp&nbsp&nbsp"),
            em(Descriptions$V2[which(Descriptions$V1 == input$variable)]),
            br(),
            HTML(
              paste0("<b>Alternative Labels: </b>", "<i>", altNames, "</i>")
            ))
  })
  
  output$description2 <- renderUI({
    altNames <-
      as.character(Descriptions$V3[which(Descriptions$V1 == input$variable)])
    
    if (grepl("\\^", altNames))
      altNames <- gsub("\\^", ", ", altNames)
    else
      altNames <- altNames
    
    tagList(
      strong(style="color:blue;font-size:30px",input$variable),
      HTML("&nbsp&nbsp&nbsp&nbsp"),
      em(Descriptions$V2[which(Descriptions$V1 == input$variable)]),
      br(),
      HTML(
        paste0("<b>Alternative Labels: </b>", "<i>", altNames, "</i>")
      ))
  })
  
  output$refTable <- renderTable({
    read.csv("References.csv", sep="\t", head = T)
    
  })
  
  
  output$proteinPlot <- renderPlotly({
    
    event.data <- event_data("plotly_click", source = "select")
    
    if(!is.null(event.data))
    {
      
    }
      
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
        "Product: ",
        ProtData$protName,
        "<br>",
        "IncGroup: ",
        ProtData$IncGroup,
        "<br>",
        "Plasmid: ",
        ProtData$Plasmid,
        "<br>",
        "Cluster ID: ",
        ProtData$cid,
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
      
      ProtData <- CG[which(CG$V1 == input$variable), ]
      
      names(ProtData) <- c("4LetterName","cid","protName","IncGroup","Plasmid","uid","SeqLength","Identity","Ex")
      
      if (input$viewProteins != "All")
      {
        ProtData <- ProtData[which(grepl(input$viewProteins,ProtData$Ex)),"Plasmid"]
      }
      else
      {
        ProtData <- ProtData[which(ProtData$Ex != "NONE"), "Plasmid"]
      }
      #rownames(Genome)[rownames(Genome) %in% unique(ProtData)]
      #Genome [rownames(Genome) %in% unique(ProtData),]
      #Genome[(as.character(rownames(Genome)) %in% ProtData),]
      
      p <- plot_ly( 
                    type="heatmap", z = Genome[(as.character(rownames(Genome)) %in% ProtData),],
                    y = ProtData,   x = colnames(Genome)   )
      
      hide_colorbar(p)
    }
  )
  
  
  
  
  
  
  
  
  output$bigPlot <- renderImage({
    
    ProtData <- processData(TRUE)
    
    # Generate the plot.
    ggplot(ProtData, aes(x = c(1:nrow(ProtData)), y = ProtData$SeqLength)) + # Notice X position is determined by V7
      geom_bar(stat = "identity",
               fill = ProtData$V6,
               alpha = as.double(ProtData$Identity) / 100.0) +
      scale_y_continuous() +
      scale_x_discrete() + # Remove extra space
      
      # Protein name labels
      geom_label(
        size = .90,
        y = max(ProtData$SeqLength) / 2.0,
        fontface = ifelse(ProtData$Ex >= 1, "bold", "plain"),
        label = paste(ProtData$protName,ProtData$IncGroup,ProtData$Plasmid,sep="_"),
        color = ifelse(ProtData$Ex == 2, "red", "black"),
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
  
  
  
  
  
  
  
  
  getHeight <- function(){
   
    return(
      exprToFunction(
        16*read.tree(file.path("SeqTreesTest", paste0(input$variable,"_",input$clusterNum, ".nwk")))$Nnode
        )
      )
    }
  
  output$ptrees <- renderPlot({
    # read phylogenetic trees from previously computed newick tree
    # this code may fail if the family only has one cluster
    
    fpath <-
    file.path("SeqTreesTest", paste0(input$variable,"_",input$clusterNum, ".nwk"))
      #file.path("SeqTreesTest", "Dtr_16210.nwk")
    validate(need(
      file.exists(fpath),
      paste(
        input$variable,
        " tree unavailable due to memory constraints."
      )
    ))
    tree <- read.tree(fpath)
    ggtree(tree) + geom_tiplab(color = "purple",
                               hjust = 0.0,
                               vjust = -0.75) +
      geom_nodepoint(color = "#b5e521",
                     alpha = 3 / 8,
                     size = 8) +
      geom_tippoint(color = "purple",
                    shape = 20,
                    size = 4) + xlim_expand(0.4, panel = "default")
    
  },width = 'auto', height = getHeight())
  
  
  
  
  
  output$Align <- downloadHandler(filename = paste0(input$variable, ".faa"),
  
  content = function(file) {
    seqAlign <-
      read.fasta(paste0("SeqAlignments_bycluster/", input$variable, "_", input$clusterNumMSA,".faa"), "AA")
    write.fasta(seqAlign, names = attr(seqAlign, "name"), file)
  }, contentType = "text/faa")


  output$Seqs <- downloadHandler(filename = function() {paste0(input$variable,".faa")},
                                
                                content = function(file) {
                                  if (input$refOnly == "All")
                                  {
                                  seqAlign <-
                                    read.fasta(paste0("Seqs_all/", input$variable, ".faa"), "AA")
                                  }
                                  else if (input$refOnly == "proteins grouped with well-studied plasmids")
                                  {
                                    seqAlign <-
                                      read.fasta(paste0("Seqs_groupedwith_ex/", input$variable, ".faa"), "AA")
                                  }
                                  else if (input$refOnly == "proteins directly from well-studied plasmids")
                                  {
                                    seqAlign <-
                                      read.fasta(paste0("Seqs_only_ex/", input$variable, ".faa"), "AA")
                                  }
                                  else if (input$refOnly == "all EXCEPT proteins grouped with well-studied plasmids")
                                  {
                                    seqAlign <-
                                      read.fasta(paste0("Seqs_allbut_ex/", input$variable, ".faa"), "AA")
                                  }
                                  write.fasta(seqAlign, names = attr(seqAlign, "name"), file)
                                }, contentType = "text/faa")
  
  output$clustNum <- renderUI(
    {
      filenames <- list.files("./SeqTreesTest/")
      filenames <- filenames[which(grepl(input$variable,filenames))]
      filenames <- gsub(".nwk","",filenames)
      numOptions <- gsub("^.*_","",filenames)
      selectInput("clusterNum", "Cluster ID:",choices=numOptions)
    }
  )
  
  output$clustNumMSA <- renderUI(
    {
      filenames <- list.files("./SeqAlignments_bycluster/")
      filenames <- filenames[which(grepl(input$variable,filenames))]
      filenames <- gsub(".faa","",filenames)
      numOptions <- gsub("^.*_","",filenames)
      selectInput("clusterNumMSA", "Cluster ID:",choices=numOptions)
    }
  )
})