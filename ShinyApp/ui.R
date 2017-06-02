library(shiny)
library(plotly)
#setwd("~/Documents/Projects/BioPythonTools/ProteinNaming-master")

# Load items to be used on both the server and UI, save them to global space
References <<-
  read.csv(file = "90PlasmidsAndIncGroups.csv", head = FALSE, sep = ",")
ClusterGroups <<-
  read.csv(file = "ClusterGroups_6-1-2017.csv", head = FALSE, sep = ",")
ClusterGroups <<- na.omit(ClusterGroups)
Descriptions <<-
  read.csv(file = "BackboneDescriptions_6-2-2017.csv", head = FALSE, sep =
             ",")

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
    
    radioButtons("arrange", "Order By:", c("AA Cluster", "Sequence Length")),
    textInput("ranges", label = "Select Ranges (1-5, 10-20, ...):"),
    checkboxInput("Subset", "Only Select These Ranges", FALSE),
    
    # To use the local download feature, uncomment the following line
    # AND add a comma after checkboxInput
    downloadButton('downloadData', 'Download (may take ~5 min)')
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
          )
        ))
    )
    ))
