library(shiny)
library(plotly)
#setwd("~/Documents/Projects/BioPythonTools/ProteinNaming-master")

# Load items to be used on both the server and UI, save them to global space
Exemplars <<-
  c("pNDM-1_Dok01", "F", "RP4", "R751", "pRA3", "R7K", "pSK41")
References <<-
  read.csv(file = "90PlasmidsAndIncGroups.csv", head = F, sep = ",")
CG <<-
  read.csv(file = "ClusterGroups_2-5-2018_1.csv", head = F, sep = ",", stringsAsFactors = F)
CG <<- na.omit(CG)

# # Sort out which ones are grouped with well-studied plasmids and label
# findRefs <- function(psmid)
# {
#   paste0(CG[which(CG$V2 %in% CG[which(CG$V9 == psmid),"V2"]),"V9"],"~") ->> CG[which(CG$V2 %in% CG[which(CG$V9 == psmid),"V2"]),"V9"]
#   gsub(paste0(psmid,"~"),psmid,CG$V9) ->> CG$V9
#   gsub("NONE~",paste0(psmid,"~"),CG$V9) ->> CG$V9
# }
# lapply(Exemplars,findRefs)

Descriptions <<-
  read.csv(file = "BackboneDescriptions_11-9-2017.csv", head = FALSE, sep =",")

Genome <<- as.matrix( read.table(file = "matrixData_2-5-2018_1.csv",sep=",",header = TRUE))

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

h4("Plasmid Backbone Groups Visualization Tool (last updated 3/7/18)"),
selectInput("variable", label = NULL,
            choices = names(setNames(
              unique(Descriptions$V1),
              unique(Descriptions$V1)
            ))),
tabsetPanel(
  tabPanel(
    HTML("<b>1. Plot of protein families</b>"),
    
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
    tabsetPanel(
      tabPanel(
        HTML("<i>1a. Plot</i>"),
        sidebarLayout(
          sidebarPanel(
            
            radioButtons(
              "refOnly",
              "Filters",
              c(
                "All",
                "proteins grouped with well-studied plasmids",
                "proteins directly from well-studied plasmids",
                "all EXCEPT proteins grouped with well-studied plasmids"
              ),
              selected = "proteins grouped with well-studied plasmids"
            ),
            
            downloadButton("Seqs", HTML("<b>Download selected sequences</b>")), HTML("<br><br>"),
            
             selectInput(
              "arrange",
              "Order By:",
              choices = c(
                "AA Cluster",
                "Sequence Length",
                "Plasmid Name",
                "Incompatibility Group"
              )
            ),
            selectInput("incgroup", "Filter by Inc Group:",
                        choices = names(setNames(
                          unique(References$V2), unique(References$V2)
                        ))),
            checkboxInput("incgroupOnly", "View this Inc Group only", FALSE),
            checkboxInput("highlightCentroids", "Highlight centroids", FALSE)
            # textInput("ranges", label = "Select Ranges (1-5, 10-20, ...):"),
            # checkboxInput("Subset", "Only Select These Ranges", FALSE),
            # # To use the local download feature, uncomment the following line
            # # AND add a comma after checkboxInput
            # downloadButton('downloadData', 'Download (may take ~5 min)')
          ),
          
          mainPanel(
            uiOutput("description"),
            plotlyOutput("proteinPlot")
          )
        )
      ),
      tabPanel(
        HTML("<i>1b. Description and User Guide</i>"),
        fluidRow(
          column(4,
        h3("Summary"),
        HTML(
          "<p>This plot shows multiple distinct protein families in which at least one member had the selected product name. </p>"
        ),
        h3("Colors"),
        
        p(
          "Each bar on the graph represents a different protein.
          The length of the bar represents the protein's sequence length.
          Each different color represents a distinct family.
          "
        ),
        
        h3("Hover Info"),
        p("Hovering the mouse over a protein will display the following information:"),
        tags$ul(tags$li(HTML("<b>Product:</b> The product name assigned to this protein in the GenBank record.")),
                tags$li(HTML("<b>IncGroup:</b> The incompatibility group, as found in the table of 90 reference plasmids. 
                         Often unavailable and displayed as a question mark.")),
                tags$li(HTML("<b>Plasmid:</b> The plasmid on which this protein is encoded.")),
                tags$li(HTML("<b>% Identity w/ Centroid:</b> The Amino Acid sequence identity of this protein with the centroid of this cluster."))
                ),
        h3("Filters"),
        p("On the lefthand sidebar, the user can select different groups of proteins to view:"),
        tags$ul(tags$li(HTML("<b>All:</b> Selecting this option will cause every protein associated with the selected four-letter label to be displayed.")),
                tags$li(HTML("<b>Proteins grouped with well studied plasmids (default):</b> Only displays proteins which either come from a well-studied plasmid OR are in a cluster (protein family) with a protein from a well-studied plasmid.")),
                tags$li(HTML("<b>Proteins directly from well-studied plasmids:</b> Only displays the proteins which come directly from a well-studied plasimd.")),
                tags$li(HTML("<b>All EXCEPT proteins grouped with well-studied plasmids:</b> The first option minus the second option."))
                )
          ),
        column(8,img(src = 'desc.png'))
          )
        ),
      tabPanel(
        HTML("<i>1c. Large plot with labels</i>"),
        h3("Description"),
        p("The following plot shows the same information as the plot in 1a, but larger so that all of the product names can be seen at a glance."),
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
            )
      # tabPanel(
      #   HTML("<i>1d. Reference Table</i>"),
      #   HTML("<b>Table from:</b> <a href=https://doi.org/10.1016/j.plasmid.2017.03.006>Annotation of plasmid genes</a>"),
      #   h5(
      #     "Plasmid core functions: generic names plus names of paralogs in examples of different well-studied plasmids."
      #   ),
      #   HTML(readChar("table.txt", file.info("table.txt")$size))
      #   #tableOutput("refTable")
      # )
        
    )
      ),
  
  
  tabPanel(HTML("<b>2. Heat Map</b>"),
           tabsetPanel(
             tabPanel(HTML("<i>2a. Heat map</i>"),
                      sidebarLayout(
                          sidebarPanel(
                            selectInput("viewProteins",label="View proteins grouped with:",
                                        choices = c(Exemplars,"All"),selected = "All")
                          ),
                          mainPanel(
                            plotlyOutput("heatMap")
                          )
                        )
                      ),
             tabPanel(HTML("<i>2b. Description</i>"),
                      h3("Summary"),
                      p("This heat map displays the potential presence of different backbone genes
                        within each well-studied plasmid. If at least one four-letter label is associated
                        with the plasmid, the map will highlight the protein-plasmid pair in yellow. Each
                        plamsid is displayed on the y-axis. Each four-letter label is displayed on the x-axis."),
                      br(),
                      p("The goal of the heat map is to give the user a snapshot of the overall genome for each plasmid.
                        This also allows users to identify which portions of the plasmid backbone are shared by the other plasmids.")
                      )
             )
           ),
  
  tabPanel(
    HTML("<b>3. Alignments</b>"),
    h5(
      "Proteins with more than 100 families may not have alignments due to memory constraints."
    ),
    downloadButton("Align", "Download alignment"),
    br(),
    h3("Summary"),
    HTML("<p>Each alignment file contains a <a href = 'https://doi.org/10.1093/nar/gkh340'>MUSCLE</a> alignment of the representative sequences from each protein family.<p>")
  ),
  tabPanel(
    HTML("<b>4. Phylogenetic trees</b>"),
    tabsetPanel(
      tabPanel( HTML ("<i> 4a. Tree </i>"),
                uiOutput("description2"),
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
    uiOutput("clustNum"),
    plotOutput("ptrees")
      ),
    tabPanel ( HTML ("<i> 4b. Description </i>"),
               p(HTML(readChar("phylo_exp.txt", file.info("phylo_exp.txt")$size)))
                 
               )
    )
        ),
  tabPanel(
    HTML("<b>5. Explanation of Methodology</b>"),
      tabsetPanel(
        tabPanel(HTML("<i>5a. Summary</i>"),
      p(HTML(readChar("explanation.txt", file.info("explanation.txt")$size)))
        ),
      tabPanel(HTML("<i>5b. Citations</i>"),
               HTML(readChar("citations.txt", file.info("citations.txt")$size))
               )
      )
    )
  # tabPanel(HTML("<b>6. Backbone gene reference file</b>"),
  #          )
    
      )))
