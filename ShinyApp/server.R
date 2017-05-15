library(shiny)
library(datasets)
library(ggplot2)
library(seqinr)
#library(msa)  # used for sequence alignments
# msa will only output a nice alignment as a pdf
# png and animation packages allow us to convert to a png and then open as a plot
library(png) 
library(animation)
library(ape)  # used for phylogenetic trees

shinyServer(function(input, output) {
  
  getHeight<-function() 
  {
      return(
        exprToFunction(
          
          if(
          (50+13*nrow(
            if(input$Subset)
            {
              subsetData(
                ClusterGroups[which(ClusterGroups$V1==input$variable),],input$ranges)
            }
            else
            {
              ClusterGroups[which(ClusterGroups$V1==input$variable),]
            }
            )
           )<32000
          )
          {
            (50+13*nrow(
              if(input$Subset)
              {
                subsetData(
                  ClusterGroups[which(ClusterGroups$V1==input$variable),],input$ranges)
              }
              else
              {
                ClusterGroups[which(ClusterGroups$V1==input$variable),]
              }
            ))
          }
          else
          {
            32000
          }
          )
        )
  }
  
  output$proteinPlot <- renderImage({
    backbone <- input$variable
    if (!input$Subset)
    {
      if (input$arrange == "Sequence Length")
      {
      filename <- normalizePath(file.path('./BarGraphs',
                                          paste(backbone, '.png', sep='')))
      list(src = filename)
      }
      else
      {
        filename <- normalizePath(file.path('./BarGraphs',
                                            paste(backbone, '_AA.png', sep='')))
        list(src = filename)
      }
      
    }
    else if (input$Subset)
    {
      
    # x=c()
    # 
    # group <- unique(ClusterGroups$V1)
    # 
    # for (id in group)
    # {
    #   x<-c(x,(nrow(ClusterGroups[which(ClusterGroups$V1==id),])))
    # }
    # maxx <- max(x)
    # 
      CG2 <- ClusterGroups[which(ClusterGroups$V1==backbone),]
      cols <- rainbow(length(unique(CG2$V2)))
      CG2$V6 <- c(1:nrow(CG2))
      rC<<-nrow(ClusterGroups)
      
      index = 1
      for (i in unique(CG2$V2))
      {
        CG2$V6[which(CG2$V2==i)] <- cols[index]
        index = index + 1
      }

      # Chunk to ensure correct data formatting
      # These few lines may be unnecessary
      CG2$V5[which(CG2$V5=="*")] <- "100.0" #"100.0" 
      CG2$V5<-as.double(as.character(CG2$V5))
      CG2$V5[which(is.na(CG2$V5))] <- 0.0
      CG2$V4 <- as.double(CG2$V4)
      CG2$V5 <- CG2$V5/100.00

      if (input$arrange == "Sequence Length")
      {
        CG3 <- CG2[order(CG2$V4),] # CRITICAL: Allow the reordering of colors!
      }
      else
      {
        CG3 <- CG2
      }
      
      CG3$V7<-c(1:nrow(CG3)) # CRITICAL: Used for labeling and position!
      # Will be the same as V7 if no subset
      
      references <- tolower(c("R7K","R46","pRA3","pKJK5","RK2","RP4","pNDM-1_Dok01"))
      cids <- c()
      for (cRow in c(1:nrow(CG3)))
      {
        plasmid <- strsplit(as.character(CG3$V3[cRow]),"_")[[1]]
        if(length(plasmid)>2)
        {
          plasmid <- plasmid[2:length(plasmid)]
          plasmid <- paste(plasmid,collapse="_")
        }
        else
        {
          plasmid <- plasmid[2]
        }
        if (tolower(plasmid) %in% references)
        {
          cids <- c(cids, CG3$V2[cRow])
        }
      }
      cids2 <- unique(cids)
      CG3$V8 <- FALSE
      CG3[which(CG3$V2 %in% cids),"V8"] <- TRUE
      
      CG3<-subsetData(CG3,input$ranges)
      
      #CG3$V8<-c(1:nrow(CG3)) # CRITICAL: Used for labeling and position!
      
      seqData <<- as.character(CG3$V3) # CRITICAL: Used to write file
      
      ggplot(CG3,aes(x=c(1:nrow(CG3)),y=CG3$V4))+ # Notice X position is determined by V7
               geom_bar(stat="identity",
                        fill=CG3$V6,
                        alpha=CG3$V5
               )+
               scale_y_continuous()+
               scale_x_discrete()+ # Remove extra space
               geom_text(size=1.75,y=max(CG3$V4)/2.0,
                         fontface=ifelse(CG3$V8,"bold","plain"),
                         label=CG3$V3,
                         color="black") +
               geom_label(size=1.5,y=-5,
                         aes(label=rev(CG3$V7)),
                         color="black") +
              coord_flip() +
               ggtitle(backbone)+
              labs(y="Sequence Length",x="")+
               theme(
                     axis.text.y=element_text(size=1),
                     axis.title.y=element_text(size=2))

      ggsave(filename = paste("BarGraphs/",backbone,"_temp.png",sep=''), 
             height = 25+as.integer(25*(13*nrow(CG3))/150.0), 
             width=2.2*25, 
             units="mm", 
             limitsize = FALSE)
      
      list(src = paste("BarGraphs/",backbone,"_temp.png",sep=''))
    }
  },
  deleteFile = FALSE
  #height=getHeight() # WEIRD!!!
  )
  
  output$alignment <- renderPrint({
    # # read sequence alignments
    # # this code may fail if the family only has one cluster
    # fpath <- file.path(".","SeqAlignments",paste0(input$variable,".fa"),fsep = .Platform$file.sep)
    # validate(
    #   need(file.exists(fpath), paste(input$variable, "failed to align, may be too many clusters, or single clusters."))
    # )
    # #align <- read.alignment(file = fpath,format = "fasta")
    # #print(align)
    # align <- readAAStringSet(fpath, format = "fasta")
    # #print(msa(align))#, show = "complete")
    # msaPrettyPrint(msa(align),output = "pdf", file = "TempAlign.pdf", showNames="left", 
    #                showLogo="none", shadingMode = "similar", showConsensus = "top", 
    #                askForOverwrite=FALSE)
    # file.rename("TempAlign.pdf", "./www/TempAlign.pdf")
  })
  
  output$ptrees <- renderPlot({
    # # read phylogenetic trees from previously computed newick tree
    # # this code may fail if the family only has one cluster
    # fpath <- file.path(".","SeqTrees",paste0(input$variable,".nwk"),fsep = .Platform$file.sep)
    # validate(
    #   need(file.exists(fpath), paste(input$variable, "has 2 or fewer clusters, tree is meaningless"))
    # )  
    # tree <- read.tree(file = fpath)
    # plot.phylo(tree)
   })
  
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
      plasmidAA <- read.fasta("Plasmids20-200kb-6-9-2016AA.fa","AA")
      print("FINISHED FASTA READ")
      
      writeSeqs <- function(seq)
      {
        seqName <- attr(seq,"name")
        if (grepl("<unknown description>",seqName))
        {
          seqName <- strsplit(seqName," ")[[1]][1]
        }
        if (seqName %in% seqData)
        {
          print("RETURNING...")
          return(seq)
          print("DONE WITH A RETURN")
        }
      }
      Seqs <<- lapply(plasmidAA,writeSeqs)
      Seqs<-Seqs[!sapply(Seqs, is.null)] # REMOVE NULL seqs 
                                         # (since lapply return same length)
      print("DONE")
      print(input$variable)
      fourName <- input$variable
      for (i in c(1:length(Seqs)))
      {
        seqName <- attr(Seqs[i],"name")
        print("1")
        if (length(seqName)>0)
        {
          if (grepl("<unknown description>",seqName))
          {
            print("2")
            seqName <- strsplit(seqName," ")[[1]][1]
            print("3")
          }
        }
        if (i!=1)
        {
        sNames <- c(sNames,paste(fourName,seqName,sep="_"))
        }
        else
        {
        sNames <- c(paste(fourName,seqName,sep="_"))
        }
      }
      print(sNames)
      print(Seqs)
      write.fasta(Seqs,sNames,
                  file.out=paste(fourName, "faa", sep = "."))
    }
  )
})