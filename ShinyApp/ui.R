library(shiny)
#setwd("~/Documents/Projects/BioPythonTools/ProteinNaming-master")

ClusterGroups <<- read.csv(file="ClusterGroups.csv",head=FALSE,sep=",")
ClusterGroups <<- na.omit(ClusterGroups)
ClusterGroups <<- ClusterGroups[which(ClusterGroups$V1!="Fip"),]
ClusterGroups <<- ClusterGroups[which(ClusterGroups$V1!="mpR"),]

shinyUI(fluidPage(
  titlePanel("Plasmid Backbone Families"),
  sidebarLayout(
    
    sidebarPanel(
      selectInput("variable", "Backbone Protein Name:",
                choices=names(
                  setNames(
                    unique(ClusterGroups$V1),
                    unique(ClusterGroups$V1)
                  )
                )
      ),
    
      radioButtons("arrange","Order By:",c("AA Cluster","Sequence Length")),
      textInput("ranges",label="Select Ranges (1-5, 10-20, ...):"),
      checkboxInput("Subset", "Only Select These Ranges", FALSE),
    
      # To use the local download feature, uncomment the following line
      # AND add a comma after checkboxInput
      downloadButton('downloadData', 'Download (may take ~5 min)')
    ),
  
    mainPanel(
      tabsetPanel(
        tabPanel("Plot of protein families",
          h3("Plot of this protein (
             may take a few min to render)"),
          h5("Bold labeled font = matched a reference plasmid"),
          plotOutput("proteinPlot")
        ),
        tabPanel("Alignments", 
                 h5("Download and reopen the linked alignment after updating the protein family."),
                 h5("Proteins with more than 100 protein families, may not have alignments due to memory constraints."),
                 shiny::a("Open alignment.",target="_blank",href="TempAlign.pdf"),
                 textOutput("alignment")
        ),
        tabPanel("Phylogenetic trees",
                 h5("Computing phylogenetic tree may take a few minutes.  "),
                 plotOutput("ptrees"))
      )
    )
  )
)
)
