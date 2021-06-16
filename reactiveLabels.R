#Bespoke Volcano
library(ggplot2)
library(dplyr)
library(purrr)



#--- Define temputs ----
temput <- list()
temput$dataset <- readRDS("data/day80deseq_verbose.rds")
temput$logpval <- 25
temput$log2FC <- 1
temput$IWP2col <- "#B856D7"
temput$CHIRcol <- "#55A0FB"

defaultGS <- c("NR2F1", "GATA3", "INSM1", "ZNF503", 
               "FGF8", "GNG8", "LFNG", "FGFR3", "LGR5", 
               "RPRM", "CD164L2", "ZBBX", "TEKT1",
               "SKOR1", "AMPD3", "VEPH1")

# Appending rows to toptable with ".r" appended to the gene name provides better 
# Visualization of highlighted points (ie they don't get buried in the mix)
update.toptable <- function(gene, toptable, geneTable) {
  if(length(gene) == 0) { 
    return(toptable) 
  } else if(gene %in% geneTable$Gene) {
    toptable <- toptable[!toptable$Gene == paste0(gene, ".r"), ]
    toptable$Category[toptable$Gene == gene] <- toptable$Category[toptable$Gene == gene] %>%
      sanitizeCategory() # LMAO good thing IWP2 and CHIR are the same number of characters otherwise this wouldnt work
  } else {
    new_row <- toptable[toptable$Gene %in% gene, ]
    new_row$Gene <- paste0(new_row$Gene, ".r")
    toptable <- rbind(toptable, new_row)
    if(toptable$Category[which(toptable$Gene == gene)] != "Below Threshold") {
      toptable$Category[which(toptable$Gene %in% c(gene, paste0(gene, ".r") ))] <- paste0(toptable$Category[toptable$Gene %in% gene], "highlight")
    }
  }
  return(toptable)
}

update.genetable <- function(gene, toptable, geneTable) {
  if(length(gene) == 0) { 
    return(geneTable) 
  } else if(gene %in% geneTable$Gene) {
    geneTable <- geneTable[!geneTable$Gene == gene, ]
  } else {
    new_row <- toptable[toptable$Gene %in% gene, ]
    new_row$Category <- sanitizeCategory(new_row$Category)
    geneTable <- rbind(geneTable, new_row) 
  }
}

sanitizeCategory <- function(category) {
  if(length(grep("highlight", category)) > 0) {
    category[grep("highlight", category)] <- substring(category[grep("highlight", category)], 1, 4)
    return( category )
  } else {
    return(category)
  }
}

toptable <- data.frame(
  Gene = temput$dataset$gene,
  Log2FC = temput$dataset$log2FoldChange,
  Log10P = ifelse(temput$dataset$padj == 0, .Machine$double.xmin, temput$dataset$padj) %>% 
    log10() %>% 
    "*"(-1)
)

toptable$Category <- ifelse(
  toptable$Log2FC < temput$log2FC*-1 & toptable$Log10P > temput$logpval,
  "IWP2",
  ifelse(
    toptable$Log2FC > temput$log2FC & toptable$Log10P > temput$logpval,
    "CHIR",
    "Below Threshold"  
  )
) 

toptable$Category <- ifelse(
  toptable$Category == "IWP2" & toptable$Gene %in% c(defaultGS, paste0(defaultGS, ".r")),
  "IWP2highlight",
  ifelse(
    toptable$Category == "CHIR" & toptable$Gene %in% c(defaultGS, paste0(defaultGS, ".r")),
    "CHIRhighlight",
    toptable$Category
  )
)

toptable <- toptable[ , c(1,4,2,3)]



geneTable <- data.frame(
  subset(toptable, Gene %in% defaultGS) %>%
    mutate(Category = sanitizeCategory(Category))
)


#---- create table ----

toptable <- data.frame(
  Gene = temput$dataset$gene,
  Log2FC = temput$dataset$log2FoldChange,
  Log10P = ifelse(temput$dataset$padj == 0, .Machine$double.xmin, temput$dataset$padj) %>% 
    log10() %>% 
    "*"(-1)
)

#Names the categories/colors
toptable$Category <- ifelse(
  toptable$Log2FC < temput$log2FC*-1 & toptable$Log10P > temput$logpval,
  "IWP2",
  ifelse(
    toptable$Log2FC > temput$log2FC & toptable$Log10P > temput$logpval,
    "CHIR",
    "Below Threshold"  
  )
) 

toptable$Category <- ifelse(
  toptable$Category == "IWP2" & toptable$Gene %in% c(defaultGS, paste0(defaultGS, ".r")),
  "IWP2highlight",
  ifelse(
    toptable$Category == "CHIR" & toptable$Gene %in% c(defaultGS, paste0(defaultGS, ".r")),
    "CHIRhighlight",
    toptable$Category
  )
)

toptable <- toptable[ , c(1,4,2,3)]

#---- Stolen from EnhancedVolcano -----

titleLabSize = 18
subtitleLabSize = 14
captionLabSize = 14
legendLabSize = 14
axisLabSize = 18
legendPosition = "right"

#---- Theme ----

th <- theme_bw(base_size = 24) +
  
  theme(
    legend.background = element_rect(),
    
    # title, subtitle, and caption
    plot.title = element_text(
      angle = 0,
      size = titleLabSize,
      face = 'bold',
      vjust = 1),
    plot.subtitle = element_text(
      angle = 0,
      size = subtitleLabSize,
      face = 'plain',
      vjust = 1),
    plot.caption = element_text(
      angle = 0,
      size = captionLabSize,
      face = 'plain',
      vjust = 1),
    
    # axis text
    axis.text.x = element_text(
      angle = 0,
      size = axisLabSize,
      vjust = 1),
    axis.text.y = element_text(
      angle = 0,
      size = axisLabSize,
      vjust = 0.5),
    axis.title = element_text(
      size = axisLabSize),
    
    # legend
    legend.position = legendPosition,
    legend.key = element_blank(),
    legend.key.size = unit(0.5, 'cm'),
    legend.text = element_text(
      size = legendLabSize),
    title = element_text(
      size = legendLabSize),
    legend.title = element_blank())



#---- App Starts Here ----

ui <- fluidPage(
  
  tags$head(tags$style(".table .alignRight {color:black; text-align:right;}")),
  
  sidebarLayout(
    
    sidebarPanel(
        
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
                     width = "125px"),
        
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
        ), #----


     dataTableOutput("geneTable"),
     br(),
     selectizeInput("updateGene",
                    "Search Genes",
                    choices = NULL
                    ),
     actionButton("addGene",
                  "Add/Remove Gene Label")
    ),
    
    mainPanel(
      
      plotOutput("volcanoPlot", height = "600px",
                 click = "plot_click",
                 hover = "plot_hover"),
      br(),
      textOutput("hover_info"),
      
      verbatimTextOutput("debug")
    )
  )
)

server <- function(input, output) {
  
  #Define a reactive dataframe that will be updated on click/search input
  values <- reactiveValues()  
  values$toptable <- data.frame()
  values$geneTable <- data.frame()
  
  values$toptable <- toptable
  values$geneTable <- geneTable
  
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
  
  #Updates DataTable - GT based on plot click input
  observeEvent(input$plot_click, {
    click_gene <- nearPoints(toptable, 
                          input$plot_click,
                          xvar = "Log2FC",
                          yvar = "Log10P",
                          maxpoints = 1)$Gene
    values$toptable <- update.toptable(click_gene, values$toptable, values$geneTable) #Be sure to update toptable values first
    values$geneTable <- update.genetable(click_gene, values$toptable, values$geneTable)
    })
  
  # Display the gene of plot point underneath mouse
  output$hover_info <- renderText({
    gene.name <- nearPoints(toptable, input$plot_hover,
                            xvar = "Log2FC",
                            yvar = "Log10P",
                            maxpoints = 1)[[1]]
    if(length(gene.name) > 0) {
      return(paste("GENE:",as.character(gene.name)))
    } else {
      return("GENE: ")
    }
  })
  
  #Updates DataTable - GT based on selectizeInput + addGene combination
  observeEvent(input$addGene, {
    values$toptable <- update.toptable(input$updateGene, values$toptable, values$geneTable)
    values$geneTable <- update.genetable(input$updateGene, values$toptable, values$geneTable)
    })
  
  # Populates selectize input with all genes in toptable, except those that were duplicated
  # For aesthetic purposes (".r$")
  updateSelectizeInput(
    inputId = "updateGene",
    choices = toptable$Gene[!str_detect(toptable$Gene, ".r$")],
    selected = "",
    server = T
    )

  # Data table that displays all labelled genes
  output$geneTable <- renderDataTable({
    data.frame(Gene = geneTable$Gene,
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

  
  # The meat and potatoes babyyyyy
  output$volcanoPlot <- renderPlot({
    bespokeVolcano(
      toptable = values$toptable,
      geneTable = values$geneTable,
      title = paste0(
        "Day ", "80 Hair Cells"
        # substr(input$dataset, 5, nchar(input$dataset)),
        # if( input$dataset == "Day_20" ) { " Prosensory" } else {" Hair" },
        # " Cells: IWP2 vs CHIR"
      ),
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
  
  output$debug <- renderText({
    paste(is.null(values$toptable) %>% as.character(), head(values$toptable) %>% as.character())
  })
}

shinyApp(ui = ui, server = server)


