library(shiny)
library(plotly)
#setwd("~/Documents/Projects/BioPythonTools/ProteinNaming-master")

# Load items to be used on both the server and UI, save them to global space
References <<- read.csv(file="90PlasmidsAndIncGroups.csv",head=FALSE,sep=",")
ClusterGroups <<- read.csv(file="ClusterGroups_6-1-2017.csv",head=FALSE,sep=",")
ClusterGroups <<- na.omit(ClusterGroups)
Descriptions <<- read.csv(file="BackboneDescriptions_6-2-2017.csv",head=FALSE,sep=",")

shinyUI(fluidPage(
  titlePanel("Plasmid Backbone Families"),
  sidebarLayout(
    
    sidebarPanel(
      selectInput("variable", "Backbone Protein Name:",
                choices=names(
                  setNames(
                    unique(Descriptions$V1),
                    unique(Descriptions$V1)
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
                uiOutput("description"),
                 plotlyOutput("proteinPlot")
        ),
        tabPanel("Alignments", 
                 h5("Download and reopen the linked alignment after updating the protein family."),
                 h5("Proteins with more than 100 protein families, may not have alignments due to memory constraints."),
                 downloadButton("Align","Download alignment")
        ),
        tabPanel("Phylogenetic trees",
                 h5("Computing phylogenetic tree may take a few minutes.  "),
                 plotOutput("ptrees"))
      )
    )
  )
)
)
