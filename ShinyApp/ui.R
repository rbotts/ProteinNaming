library(shiny)
library(plotly)
#setwd("~/Documents/Projects/BioPythonTools/ProteinNaming-master")

# Load items to be used on both the server and UI, save them to global space
References <<-
  read.csv(file = "90PlasmidsAndIncGroups.csv", head = FALSE, sep = ",")
ClusterGroups <<-
  read.csv(file = "ClusterGroups_6-13-2017.csv", head = FALSE, sep = ",")
ClusterGroups <<- na.omit(ClusterGroups)
Descriptions <<-
  read.csv(file = "BackboneDescriptions_6-13-2017.csv", head = FALSE, sep =
             ",")
Exemplars <<- c("pNDM-1_Dok01","F","RP4","R751","pRA3","R7K","pSK41")


Genome <<- matrix(00,nrow=length(unique(References$V1)),
                               ncol = length(unique(Descriptions$V1)))
dimnames(Genome) <<- list(unique(References$V1),unique(Descriptions$V1))

# Go through each plasmid.
for (psmid in unique(References$V1))
{
  # Get ALL the 4-letter names associated with each plasmid.
  protNamesByPlasmid <- ClusterGroups[which(ClusterGroups$V5==psmid), "V1"]

  # Each 4-letter name that is associated with that plasmid is noted.
  # The way it is noted is by incrementing the value in a matrix.
  for (bckbn in unique(protNamesByPlasmid))
  {
    Genome[psmid,bckbn] <<- Genome[psmid,bckbn] + 1
  }
}

shinyUI(fluidPage(
  tags$head(
    tags$style(
      type = "text/css",
      "#loadmessage {
      position: fixed;
      top: 50%;
      left: 50%;
      ocacity: 0.50;
      text-align: center;
      font-weight: bold;
      font-size: 300%;
      color: #000000;
      z-index: 105;
      animation: blinker 1s linear infinite;
      }"
)
    ),

titlePanel("Plasmid Backbone Families"),
sidebarLayout(
  sidebarPanel(
    selectInput("variable", "Backbone Protein Name:",
                choices = names(setNames(
                  unique(Descriptions$V1),
                  unique(Descriptions$V1)
                ))),
    
    radioButtons("arrange", "Order By:", c("AA Cluster", "Sequence Length",
                                           "Plasmid Name",
                                           "Incompatibility Group")),
    radioButtons("refOnly","Only view:",
                 c("All",
                 "proteins grouped with well-studied plasmids",
                 "proteins directly from well-studied plasmids",
                 "all EXCEPT proteins grouped with well-studied plasmids"),
                 selected = "proteins grouped with well-studied plasmids"
    ),
    selectInput("incgroup", "Filter by Inc Group:",
                choices = names(setNames(unique(References$V2),unique(References$V2)))),
    checkboxInput("incgroupOnly","View this Inc Group only",FALSE),
    checkboxInput("highlightCentroids","Highlight centroids",FALSE)
    # textInput("ranges", label = "Select Ranges (1-5, 10-20, ...):"),
    # checkboxInput("Subset", "Only Select These Ranges", FALSE),
    # # To use the local download feature, uncomment the following line
    # # AND add a comma after checkboxInput
    # downloadButton('downloadData', 'Download (may take ~5 min)')
  ),
  
  mainPanel(tabsetPanel(
    tabPanel(
      "Plot of protein families",
      
      conditionalPanel(
        condition = "$('html').hasClass('shiny-busy')",
        tags$div("Loading...", id = "loadmessage"),
        tags$script(
          HTML(
            "
            (function blink() {
            $('#loadmessage').fadeOut(500).fadeIn(500, blink);
            })();
            "
          )
        )
        ),
      
      uiOutput("description"),
      plotlyOutput("proteinPlot")
          ),
    tabPanel("Heat Map", plotlyOutput("heatMap")),
    tabPanel("User Guide",
             h4("1. To see the proteins, click the tab above, 'Plot of protein families.'"),
              h4("2. On the left sidebar, use the dropdown box to select a backbone protein of interest."),
             h4("3. Hover the mouse over the graph to see information about each protein."),
             h5("3a. Each bar represents a plasmid backbone protein from NCBI."),
             h5("3b. Each color represents a group of proteins which had similar amino acid sequence identity."),
             h4("4. Drag and drop a selection box to zoom in on a group of proteins."),
             h5("4a. Use the 'Autoscale' button to zoom back out.")
             
             ),
    tabPanel(
      "Large plot with labels",
      conditionalPanel(
        condition = "$('html').hasClass('shiny-busy')",
        tags$div("Loading...", id = "loadmessage"),
        tags$script(
          HTML(
            "
            (function blink() {
            $('#loadmessage').fadeOut(500).fadeIn(500, blink);
            })();
            "
          )
        )
        ),
      #uiOutput("description"),
      plotOutput("bigPlot")
    ),
    tabPanel(
      "Alignments",
      h5(
        "Proteins with more than 100 families may not have alignments due to memory constraints."
      ),
      downloadButton("Align", "Download alignment")
    ),
    tabPanel(
      "Phylogenetic trees",
      conditionalPanel(
        condition = "$('html').hasClass('shiny-busy')",
        tags$div("Loading...", id = "loadmessage"),
        tags$script(
          HTML(
            "
            (function blink() {
            $('#loadmessage').fadeOut(500).fadeIn(500, blink);
            })();
            "
          )
        )
        ),
      plotOutput("ptrees")
          ),
    tabPanel(
      "Table of reference plasmids",
      h5("List of well-studied plasmids: pNDM-1_Dok01, F, RP4, R751, pRA3, R7K, pSK41"),
      tableOutput("refTable")
    )
        ))
    )
    ))
