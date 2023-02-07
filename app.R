#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("bespokeVolcano_v1.R")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(shinyWidgets)
library(colourpicker)
library(tools)

#Adding + removing gene labels breaks it for some reason, idk why
datasets <- list(
    Day_20 = readRDS("data/day20deseq_verbose.rds"),
    Day_80 = readRDS("data/day80deseq_verbose.rds"),
    Day_109 = readRDS("data/day109deseq_verbose.rds"),
    userFile = NULL
)

# Named Character vector ("Name" = "Character"), where the characters
# represent the names in the datasets list, and the names are displayed in the GUI
datasetDisplayNames <- c("Day 20 Prosensory" = "Day_20",
                          "Day 80 Hair Cells" = "Day_80",
                          "Day 109 Hair Cells" = "Day_109",
                          "Custom Dataset" = "userFile")

# gene.set is a named list of named character vectors, containing the genes to highlight
# when a dataset is loaded in. The names in the list should correspond to the names
# in the dataset list
gene.sets <- list(
    Day_20 = c("NR2F2", "RSPO3", "EDN3", "LRP2", "SULF1", "GAS1", "PTCH1", "OTX1", "KDM7A","TAGLN",
               "CXCR4", "JAG1", "DLX5", "GPR155", "ACSL4", "MSX1"),
    Day_80 = c("NR2F1", "GATA3", "INSM1", "ZNF503", "FGF8", "GNG8", "LFNG", "FGFR3", "LGR5", "RPRM",
               "CD164L2", "ZBBX", "TEKT1", "SKOR1", "AMPD3", "VEPH1"),
    Day_109 = c("GATA3", "NR2F1", "INSM1", "HES6", "TMPRSS3", "GNG8", "ZNF503",
                 "TEKT1", "NEUROD6", "ZBBX", "CD164L2", "PCDH20", "SKOR1", "VEPH1", "TEKT2", "TCTEX1D1"),
    userFile = c()
)

# Server Side Functions ----

# Note: Toptable is a reactive (updates in response to changing inputs) data.frame
#       that is instantiated in the server. It contains gene names, logFC, and -log10 p.vals.
#       Toptable is used as the input to bespokeVolcano, a bastardization of EnhancedVolcano 
#       from the EnhancedVolcano package

# Ensures that selected genes (and only selected genes) are duplicated in toptable
# This allows highlighted genes to be plotted twice, ensuring they aren't buried
update.toptable <- function(geneList, toptable) {
  gene.reps <- toptable[grep("\\.r$", toptable$Gene), ]
  gene.reps$Gene <- sapply(gene.reps$Gene, FUN = function(x) substr(x, 1, nchar(x)-2))
  gene.reps <- gene.reps[gene.reps$Gene %in% geneList, ]
  gene.reps <- rbind(gene.reps, toptable[toptable$Gene %in% geneList[!geneList %in% gene.reps$Gene], ])
  gene.reps$Gene <- paste0(gene.reps$Gene, ".r")
  
  # This binds all rows in the toptable that aren't replicates ( -1*grep("\\.r$", toptable$Gene) )
  # with the new list of replicated genes. BUT !!! if there aren't any replicated genes,
  # it returns an empty dataframe. Hence, the if clause
  if ( grep("\\.r$", toptable$Gene) %>% length() == 0 ) { 
    toptable <- rbind(toptable, gene.reps)
  } else {
    toptable <- rbind(toptable[-1*grep("\\.r$", toptable$Gene),], gene.reps)
    }
  return(toptable)
  
}


# This classifies toptable entries as being above/below p-value and fold change 
# cutoffs
# Note: Nested ifelse statements are, hard to debug, and should be refactored eventually
categorize.toptable <- function(toptable, geneList, log2FC, log10P) {
  
  toptable$Category <- ifelse(
    toptable$Log2FC < log2FC*-1 & toptable$Log10P > log10P,
    "IWP2",
    ifelse(
      toptable$Log2FC > log2FC & toptable$Log10P > log10P,
      "CHIR",
      "Below Threshold"  
    )
  ) 
  
  toptable$Category <- ifelse(
    toptable$Category == "IWP2" & toptable$Gene %in% c(geneList,
                                                       paste0(geneList, ".r")),
    "IWP2highlight",
    ifelse(
      toptable$Category == "CHIR" & toptable$Gene %in% c(geneList,
                                                         paste0(geneList, ".r")),
      "CHIRhighlight",
      toptable$Category
    )
  )
  
  toptable <- toptable[ , c(1,4,2,3)]
  toptable <- toptable %>% arrange(Category)

  return(toptable)
}

makeGeneTable <- function(toptable, geneList) {
  subset(toptable, Gene %in% geneList) %>%
    mutate(Category = sanitizeCategory(Category))
}

updateGeneList <- function(gene, geneList) {
  if(length(gene) == 0) { 
    return(geneList) 
  } else if(gene %in% geneList) {
    geneList <- geneList[!geneList == gene]
  } else {
    geneList <- c(gene, geneList)
  }
  return(geneList)
}

sanitizeCategory <- function(category) {
  if(length(grep("highlight", category)) > 0) {
    category[grep("highlight", category)] <- substring(category[grep("highlight", category)], 1, 4)
    return( category )
  } else {
    return(category)
  }
}

textInputRow<-function (inputId, label, value = ""){
  div(style="display:inline-block; width:150px; padding-bottom:10px; padding-right:50px",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value, class="input-small"))
} 

validateTable <- function(table){
  if ("gene" %in% names(table) &
      "log2FoldChange" %in% names(table) & 
      "padj" %in% names(table)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# ----Define UI ----
ui <- fluidPage(
    
    #Suppress warnings
    tags$style(
        type="text/css"
    ),
    tags$head(tags$style(
      "#hover_info{color: black;
              font-size: 20px;
              font-style: bold;
              }"
    )),

    # Application title
    titlePanel("Volcano Plot Customization"),
    
    actionButton("geneCustomize", "Gene Table"), #actionButtons iterate by 1 whenever they're clicked
    
    actionButton("plotCustomize", "Plot Options"), #actionButtons iterate by 1 whenever they're clicked
    
    actionButton("colorCustomize", "Color Options"), #actionButtons iterate by 1 whenever they're clicked
    
    actionButton("downloadCustomize", "Download Options"), #actionButtons iterate by 1 whenever they're clicked

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          
          selectInput("dataset", "Dataset",
                      datasetDisplayNames,
                      width = '200px'
          ),
          
          fileInput("userFile", "Upload new data as .csv",
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")
                    ),
          
          
          conditionalPanel( # Gene Table Controls
            #-----------
            condition = "output.showGeneFinder",
            dataTableOutput("geneTable"),
            br(),
            selectizeInput("updateGene",
                           "Search Genes",
                           choices = NULL),
            actionButton("addGene",
                         "Add/Remove Gene Label")
          ), #----
          
          conditionalPanel( # Plot Controls
            #-----------
            condition = "output.showPlotControls",
            
            sliderInput("logpval",
                        "-Log10 p-value",
                        min = 0,
                        max = 300,
                        value = 10),
            
            sliderInput("log2FC",
                        "Fold Change Cutoff",
                        min = 0,
                        max = 3,
                        value = 1,
                        step = 0.1),
            
            sliderInput("pointSize",
                        "Point Size",
                        min = 0.5,
                        max = 12,
                        value = 5,
                        step = 0.25),
            
            sliderInput("pointAlpha",
                        "Point Transparency",
                        min = 0.05,
                        max = 1,
                        value = 0.5,
                        step = 0.05),
            h5("Note:"),
            helpText("Point Transparency must be set to 1 to be compatible with the .eps file format"),
            
            sliderInput("xlim",
                        "X-Axis Bounds",
                        min = 3,
                        max = 12,
                        value = 6,
                        step = 0.25),
            
            sliderInput("labelSize",
                        "Label Size",
                        min = 4,
                        max = 24,
                        value = 8,
                        step = 1),
            
            prettySwitch("drawLabels",
                         "Labels",
                          value = TRUE,
                          width = "125px")
              
          ), #----
          
          conditionalPanel( # Color Controls
            #----
            condition = "output.showColorControls",
            fluidRow(
              column(3,
                     helpText("CHIR"),
                     colourInput("CHIRcol", NULL, "#55A0FB"),
                     br(),
                     helpText("CHIR Highlight"),
                     colourInput("CHIRhighlight", NULL, "#333194")
              ),
              column(3,
                     helpText("IWP2"),
                     colourInput("IWP2col", NULL, "#B856D7"),
                     br(),
                     helpText("IWP2 Highlight"),
                     colourInput("IWP2highlight", NULL, "#8b0000")
              ),
            ) 
          ), #----
          
          conditionalPanel( # Download Controls
            #----
            condition = "output.showDownloadControls",
            textInputRow(inputId="downloadWidth", label="Download Width (in)", value = 8.0),
            textInputRow(inputId="downloadHeight", label="Download Height (in)", value = 6.0),
            br(),
            downloadButton("VolcanoPlotImage", "Download Plot as .png"),
            downloadButton("VolcanoPlotEPS", "Download Plot as .eps file"),
            br(),
          ), #----
        ),
      
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("volcanoPlot", height = "600px",
                      click = "plot_click",
                      hover = "plot_hover"),
           br(),
           textOutput("hover_info"),
           br(),
           helpText("Hover over a data point to show the gene name (above)"),
           helpText("Click a data point to highlight it on the plot"),
           helpText("Search for a gene (under the Gene Table tab) to highlight its position on the plot"),
           #verbatimTextOutput("debug")
        )
    )
)


# Define server logic required to create/update the Volcano Plot
server <- function(input, output) {
  
  #---- Define Reactive Values ----
  values <- reactiveValues()
  values$toptable <- data.frame()
  values$geneTable <- data.frame()
  values$geneList <- character()
  values$conditional <- "gene"
  values$datasets <- datasets
  values$validUserData <- FALSE
  values$errorMessage <- "Unknown Error Encountered"
  values$plotTitle <- " "
  
  observeEvent(input$userFile, {
    
    inFile <- input$userFile
    if(!file_ext(inFile$datapath) == "csv") {
      values$errorMessage <- "Selected file must be a .csv"
      values$validUserData <- FALSE
      return()
    } 
    
    inputData <- read.csv(inFile$datapath)
    if(!validateTable(inputData)){
      values$errorMessage <- "Input data must contain the following columns: \"gene\", \"log2FoldChange\", and \"padj\"."
      values$validUserData <- FALSE
      return()
    }
    
    values$errorMessage <- "Unknown Error Encountered"
    values$validUserData <- TRUE
    values$datasets$userFile <- inputData
  })
  
  #Update all fields when dataset is changed
  observeEvent(input$dataset, {
    
    values$geneList <- gene.sets[[input$dataset]]
    
    # For user input data, need to populate the geneList with at least one entry
    if (input$dataset == "userFile"){
      if (values$validUserData){
        defaultGene <- c(
          values$datasets[[input$dataset]][order(values$datasets[[input$dataset]]$padj, decreasing = FALSE)[1], "gene"])
        values$geneList <- updateGeneList(defaultGene, values$geneList)
        values$plotTile <- "Custom Volcano Plot"
      } else {
        return()
      }
    }
    
    values$toptable <- data.frame(
      Gene = values$datasets[[input$dataset]]$gene,
      Log2FC = values$datasets[[input$dataset]]$log2FoldChange,
      Log10P = ifelse(values$datasets[[input$dataset]]$padj == 0, .Machine$double.xmin, values$datasets[[input$dataset]]$padj) %>%
        log10() %>%
        "*"(-1)
    )
    
    if (input$dataset != "userFile"){
      values$plotTile <- paste0(
        "Day ",
        substr(input$dataset, 5, nchar(input$dataset)),
        if( input$dataset == "Day_20" ) { " Prosensory" } else {" Hair" },
        " Cells: IWP2 vs CHIR"
      )
    }
    
  })

  colorMap <- reactive({
    list(
      IWP2 = input$IWP2col,
      CHIR = input$CHIRcol,
      IWP2highlight = input$IWP2highlight,
      CHIRhighlight = input$CHIRhighlight,
      `Below Threshold` = 'lightgrey',
      textCol = "black"
    )
  })

  #---- Conditional Panel Control ----

  output$showGeneFinder <- reactive({
    values$conditional == "gene"
  })

  output$showPlotControls <- reactive({
    values$conditional == "plot" # Add whatever condition you want here. Must return TRUE or FALSE
  })

  output$showColorControls <- reactive({
    values$conditional == "color" # Add whatever condition you want here. Must return TRUE or FALSE
  })

  output$showDownloadControls <- reactive({
    values$conditional == "download" # Add whatever condition you want here. Must return TRUE or FALSE
  })
  
  observeEvent(
    input$geneCustomize, {
    values$conditional <- "gene" }
  )
  
  observeEvent(
    input$plotCustomize, {
      values$conditional <- "plot" }
  )
  
  observeEvent(
    input$colorCustomize, {
      values$conditional <- "color" }
  )  
  
  observeEvent(
    input$downloadCustomize, {
      values$conditional <- "download" }
  )

  outputOptions(output, "showGeneFinder", suspendWhenHidden = FALSE)

  outputOptions(output, 'showPlotControls', suspendWhenHidden = FALSE)

  outputOptions(output, "showColorControls", suspendWhenHidden = FALSE)

  outputOptions(output, "showDownloadControls", suspendWhenHidden = FALSE) # lmao idk what suspend when hidden even doessssss

  #---- Gene Updating ----
  #Updates DataTable - GT based on plot click input
  observeEvent(input$plot_click, {
    click_gene <- nearPoints(values$toptable,
                             input$plot_click,
                             xvar = "Log2FC",
                             yvar = "Log10P",
                             maxpoints = 1)$Gene
    
    values$geneList <- updateGeneList(click_gene, values$geneList)
    values$toptable <- update.toptable(values$geneList, values$toptable)
  })

  # Updates DataTable - GT based on selectizeInput + addGene combination
  observeEvent(input$addGene, {
    validate(
      need(input$dataset != "userFile" || values$validUserData,
           "User file is invalid")
    )
    values$geneList <- updateGeneList(input$updateGene, values$geneList)
    values$toptable <- update.toptable(values$geneList, values$toptable)
  })

  # Populates selectize input with all genes in toptable, except those that were duplicated
  # For aesthetic purposes (".r$")
  observe({
    updateSelectizeInput(
      inputId = "updateGene",
      choices = values$toptable$Gene[!str_detect(values$toptable$Gene, ".r$")],
      selected = "",
      server = T
    )
  })

  #---- Output ----
  
  output$geneTable <- renderDataTable({
    validate(
      need(input$dataset != "userFile" || values$validUserData,
           "User file is invalid")
    )
    
    geneTable <- makeGeneTable(
      categorize.toptable(values$toptable,
                          values$geneList,
                          input$log2FC,
                          input$logpval),
      values$geneList
    )
    geneTable
    data.frame(
      Gene = geneTable$Gene,
      Condition = geneTable$Category,
      Log2FC = geneTable$Log2FC %>% round(2),
      Log10P = geneTable$Log10P %>% round(0)
    )
  },
    options = list(# Javascript Options
      pageLength = 10,
      aoColumnDefs=list(
        list(
          sClass = "alignRight",
          aTargets = c(list(2), list(3))
        )
      )
    )
  )
  
  # Display the gene of plot point underneath mouse
  output$hover_info <- renderText({
    validate(
      need(input$dataset != "userFile" || values$validUserData,
           "User file is invalid")
    )
    
    gene.name <- nearPoints(values$toptable,
                            input$plot_hover,
                            xvar = "Log2FC",
                            yvar = "Log10P",
                            maxpoints = 1)[[1]]
    if(length(gene.name) > 0) {
      return(paste("GENE:", as.character(gene.name)))
    } else {
      return("GENE: ")
    }
  })
  
  reactiveVolcano <- reactive ({
    bespokeVolcano(
      toptable = categorize.toptable(
        values$toptable,
        values$geneList,
        input$log2FC,
        input$logpval
        ),
      geneTable = makeGeneTable(
        categorize.toptable(
          values$toptable,
          values$geneList,
          input$log2FC,
          input$logpval),
        values$geneList
      ),
      title = values$plotTitle,
      pCutoff = input$logpval,
      FCcutoff = input$log2FC,
      pointSize = input$pointSize,
      colAlpha = input$pointAlpha,
      xlim = input$xlim,
      labSize = input$labelSize,
      drawLabel = input$drawLabels,
      colorMap = colorMap()
    )
  })
  
  output$volcanoPlot <- renderPlot({
    validate(
      need(input$dataset != "userFile" || values$validUserData,
           values$errorMessage)
    )
    reactiveVolcano()
  })

  output$VolcanoPlotImage <- downloadHandler(
      filename = function() { paste(input$dataset, "_VolcanoPlot", '.png', sep='') },
      content = function(file) {
          ggsave(file,
                 plot = reactiveVolcano(),
                 device = "png",
                 units = "in",
                 width = as.numeric(input$downloadWidth),
                 height = as.numeric(input$downloadHeight),
                 dpi = 300
          )
      }
  )

  output$VolcanoPlotEPS <- downloadHandler(
      filename = function() { paste(input$dataset, "_VolcanoPlot", '.eps', sep='') },
      content = function(file) {
          ggsave(file,
                 plot = reactiveVolcano(),
                 device = "eps",
                 units = "in",
                 width = as.numeric(input$downloadWidth),
                 height = as.numeric(input$downloadHeight),
                 dpi = 300
          )
      }
  )
  
  output$debug <- renderText({
    values$conditional
    values$geneList
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
