library(shiny)
library(plotly)
#setwd("~/Documents/Projects/BioPythonTools/ProteinNaming-master")

# Load items to be used on both the server and UI, save them to global space
References <<- read.csv(file="90PlasmidsAndIncGroups.csv",head=FALSE,sep=",")
ClusterGroups <<- read.csv(file="ClusterGroups_5-30-17.csv",head=FALSE,sep=",")
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
          h4("Legend"), 
          strong("Red colored and bold font:",style="color:red"), br(),
          em("This sequence comes from a reference plasmid."), br(), br(),
          strong("Black colored and bold font:"), br(),
          em("This sequence is in a cluster that contains at least one sequence from a reference plasmid."), br(), br(),
          div("Black colored and plain font:"),
          em("This sequence neither came from a reference plasmid nor is in a cluster that contains one."), br(),
         
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
