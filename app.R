#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("data/EnhancedVolcano.R")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggalt)
library(ggrastr)
library(shinyWidgets)
library(colourpicker)
library(plotly)

IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

datasets <- list(
    Day_20 = readRDS("data/day20deseq_verbose.rds"),
    Day_80 = readRDS("data/day80deseq_verbose.rds"),
    Day_109 = readRDS("data/day109deseq_verbose.rds")
)

gene.sets <- list(
    Day_20 = c("NR2F2", "RSPO3", "EDN3", "LRP2", "SULF1", "GAS1", "PTCH1", "OTX1", "KDM7A",
      "TAGLN", "CXCR4", "JAG1", "DLX5", "GPR155", "ACSL4", "MSX1"),
    Day_80 = c("NR2F1", "GATA3", "INSM1", "ZNF503", "FGF8", "GNG8", "LFNG", "FGFR3", "LGR5", "RPRM",
               "CD164L2", "ZBBX", "TEKT1", "SKOR1", "AMPD3", "VEPH1"),
    Day_109 = c("GATA3", "NR2F1", "INSM1", "HES6", "TMPRSS3", "GNG8", "ZNF503",
                 "TEKT1", "NEUROD6", "ZBBX", "CD164L2", "PCDH20", "SKOR1", "VEPH1", "TEKT2", "TCTEX1D1")
)

textInputRow<-function (inputId, label, value = ""){
  div(style="display:inline-block; width:150px; padding-bottom:10px; padding-right:50px",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value, class="input-small"))
}             

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    #Suppress warnings
    tags$style(
        type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),

    # Application title
    titlePanel("Volcano Plot Customization"),
    
    actionButton("plotCustomize", "Show/Hide Plot Options"),
    
    actionButton("colorCustomize", "Show/Hide Color Options"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          
          selectInput("dataset", "Dataset",
                      c("Day 20 Prosensory" = "Day_20",
                        "Day 80 Hair Cells" = "Day_80",
                        "Day 109 Hair Cells" = "Day_109"),
                      width = '200px'
          ),
          
          conditionalPanel(
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
                        value = 4,
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
            
            
          ),
          #----
          
          conditionalPanel(
            #----
            condition = "output.showColorControls",
            fluidRow(
              column(3,
                     helpText("CHIR"),
                     colourInput("chirCol", NULL, "#55A0FB"),
                     br(),
                     helpText("CHIR Highlight"),
                     colourInput("chirHighlight", NULL, "#333194")
              ),
              column(3,
                     helpText("IWP2"),
                     colourInput("iwp2Col", NULL, "#B856D7"),
                     br(),
                     helpText("IWP2 Highlight"),
                     colourInput("iwp2Highlight", NULL, "#8b0000")
              ),
            ) #----
          ),
          
          column(width = 6,
                 verbatimTextOutput("clickTable")
          ),
          br(),
            
          #Managing plot downloads
          textInputRow(inputId="downloadWidth", label="Download Width (in)", value = 8.0),
          textInputRow(inputId="downloadHeight", label="Download Height (in)", value = 6.0),
          br(),
          downloadButton("VolcanoPlotImage", "Download Plot as .png"),
          downloadButton("VolcanoPlotEPS", "Download Plot as .eps file")
        ),
      
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("volcanoPlot", height = "600px",
                      click = "plot_click")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$showPlotControls <- reactive({
    input$plotCustomize%%2 # Add whatever condition you want here. Must return TRUE or FALSE
  })
  
  output$showColorControls <- reactive({
    input$colorCustomize%%2 # Add whatever condition you want here. Must return TRUE or FALSE
  })
  
  outputOptions(output, 'showPlotControls', suspendWhenHidden = FALSE)
  
  outputOptions(output, "showColorControls", suspendWhenHidden = FALSE)
    
  reactiveVolcano <- reactive ({
      vp.data <- datasets[[input$dataset]]
      genes.of.interest <- gene.sets[[input$dataset]]
      
      color.key <- ifelse(
          vp.data$log2FoldChange < input$log2FC*-1 & vp.data$padj < 10^(-input$logpval),
          input$iwp2Col, 
          ifelse(
              vp.data$log2FoldChange > input$log2FC & vp.data$padj < 10^(-input$logpval),
              input$chirCol,
              "lightgrey"
          )
      )
      
      names(color.key) <- ifelse(vp.data$log2FoldChange < input$log2FC*-1 & vp.data$padj < 10^(-input$logpval),
                                 "IWP2",
                                 ifelse(
                                     vp.data$log2FoldChange > input$log2FC & vp.data$padj < 10^(-input$logpval),
                                     "CHIR",
                                     "Below Threshold"  
                                 )
      )
      
      #Define color for highlighted genes
      color.key[which(vp.data$gene %in% c(genes.of.interest, paste0(genes.of.interest, ".r"))  &
                        vp.data$log2FoldChange > 0)] <- input$chirHighlight
      
      color.key[which(vp.data$gene %in% c(genes.of.interest, paste0(genes.of.interest, ".r"))  &
                        vp.data$log2FoldChange < 0)] <- input$iwp2Highlight
      
      names(color.key)[which(vp.data$gene %in% c(genes.of.interest, paste0(genes.of.interest, ".r")) &
                             vp.data$log2FoldChange > 0)] <- "Dorsal Marker Genes"
      
      names(color.key)[which(vp.data$gene %in% c(genes.of.interest, paste0(genes.of.interest, ".r")) &
                               vp.data$log2FoldChange < 0)] <- "Ventral Marker Genes"
      
      
      # construct the volcano plot
      EnhancedVolcano(
          toptable = vp.data,
          lab = vp.data$gene,
          x = "log2FoldChange",
          y = "padj",
          title = paste0("Day ", 
                         substr(input$dataset, 5, nchar(input$dataset)),
                         if(input$dataset == "Day_20") { " Prosensory" } else { " Hair" },
                         " Cells: IWP2 vs CHIR"),
          pCutoff = 10^(-input$logpval),
          FCcutoff = input$log2FC,
          pointSize = input$pointSize,
          colCustom = color.key,
          selectLab = if(input$drawLabels){genes.of.interest}else{"XYZ"},
          legendPosition = "right",
          labSize = input$labelSize,
          colAlpha = input$pointAlpha,
          drawConnectors = T,
          arrowheads = F
      )  + NoLegend() +
          xlim(c(input$xlim*-1, input$xlim))  
  }) %>% 
      bindCache(input$log2FC, 
                input$logpval,
                input$pointSize,
                input$drawLabels,
                input$labelSize,
                input$pointAlpha,
                input$xlim,
                input$dataset,
                input$iwp2Highlight,
                input$chirHighlight,
                input$iwp2Col,
                input$chirCol)

  output$volcanoPlot <- renderPlot({
      reactiveVolcano()
  }) %>% 
      bindCache(input$log2FC, 
                input$logpval,
                input$pointSize,
                input$drawLabels,
                input$labelSize,
                input$pointAlpha,
                input$xlim,
                input$dataset,
                input$iwp2Highlight,
                input$chirHighlight,
                input$iwp2Col,
                input$chirCol)
  
  output$clickTable <- renderDataTable(

    nearPoints(
      mutate(
        datasets[[input$dataset]],
        log10FC = -1*log10(padj)
        ),
      input$plot_click,
      maxpoints = 1)
  )
  
  output$click_info <- renderPrint({
    cat("input$plot_click:\n")
    str(input$plot_click)
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
}

# Run the application 
shinyApp(ui = ui, server = server)
