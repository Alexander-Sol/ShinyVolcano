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

IWP2 <- "#B856D7"
CHIR <- "#55A0FB"

datasets <- list(
    Day_20 = readRDS("data/day20deseq.rds"),
    Day_109 = readRDS("data/day109deseq.rds")
)

gene.sets <- list(
    Day_20 = c("NR2F2", "RSPO3", "EDN3", "LRP2", "SULF1", "GAS1", "PTCH1", "OTX1",
      "TAGLN", "CXCR4", "JAG1", "DLX5", "GPR155", "ACSL4", "MSX1"),
    Day_109 = c("GATA3", "NR2F1", "INSM1", "HES6", "TMPRSS3", "GNG8", "ZNF503",
                 "TEKT1", "NEUROD6", "ZBBX", "CD164L2", "PCDH20", "SKOR1", "VEPH1", "TEKT2", "TCTEX1D1")
)
                        

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    #Suppress warnings
    tags$style(
        type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),

    # Application title
    titlePanel("Day 20 Volcano Plot Customization"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            
            selectInput(
              "dataset",
              "Dataset",
              c("Day 20 Prosensory" = "Day_20",
                "Day 109 Hair Cells" = "Day_109")
            ),
            
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
            
            sliderInput("labelSize",
                        "Label Size",
                        min = 4,
                        max = 24,
                        value = 8,
                        step = 1),
            
            sliderInput("xlim",
                        "X-Axis Bounds",
                        min = 3,
                        max = 12,
                        value = 6,
                        step = 0.25),
            
            checkboxInput("drawLabels", "Labels", value = T),
            
            #Managing plot downloads
            textInput("downloadWidth", "Download Width (in)", 8),
            textInput("downloadHeight", "Download Height (in)", 6),
            downloadButton("VolcanoPlotImage", "Download Plot as .png"),
            downloadButton("VolcanoPlotEPS", "Download Plot as .eps file")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("volcanoPlot", height = "600px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    reactiveVolcano <- reactive ({
        vp.data <- datasets[[input$dataset]]
        genes.of.interest <- gene.sets[[input$dataset]]
        
        color.key <- ifelse(
            vp.data$log2FoldChange < input$log2FC*-1 & vp.data$padj < 10^(-input$logpval),
            IWP2, 
            ifelse(
                vp.data$log2FoldChange > input$log2FC & vp.data$padj < 10^(-input$logpval),
                CHIR,
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
        color.key[which(vp.data$gene %in% c(genes.of.interest, paste0(genes.of.interest, ".r")))] <- "darkred"
        names(color.key)[which(vp.data$gene %in% c(genes.of.interest, paste0(genes.of.interest, ".r")))] <- "Selected"
        
        
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
                  input$dataset)

    output$volcanoPlot <- renderPlot({
        print(reactiveVolcano())
    }) %>% 
        bindCache(input$log2FC, 
                  input$logpval,
                  input$pointSize,
                  input$drawLabels,
                  input$labelSize,
                  input$pointAlpha,
                  input$xlim,
                  input$dataset)
    
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
