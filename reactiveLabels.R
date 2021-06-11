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

defaultGS <- c("NR2F1", "GATA3", "INSM1", "ZNF503", "FGF8", "GNG8", "LFNG", "FGFR3", "LGR5", "RPRM",
                         "CD164L2", "ZBBX", "TEKT1", "SKOR1", "AMPD3", "VEPH1")

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


#Map colors to category, needed for scale_color_manual
#Probably colors should be reactive. something like
colorMap <- reactive({
  list(
    IWP2 = temput$IWP2col,
    CHIR = temput$CHIRcol,
    `Below Threshold` = 'lightgrey'
  )
})
colorMap <- list(
  IWP2 = temput$IWP2col,
  CHIR = temput$CHIRcol,
  `Below Threshold` = 'lightgrey',
  textCol = "black"
)


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
  
  sidebarLayout(
    
    sidebarPanel(
      dataTableOutput("geneTable")
      
    ),
    
    mainPanel(
      
      plotOutput("volcanoPlot", height = "600px",
                 click = "plot_click",
                 hover = "plot_hover"),
      
      textOutput("hover_info")
    )
  )
)

server <- function(input, output) {
  values <- reactiveValues()
  
  values$GT <- data.frame(
    subset(toptable, Gene %in% reactiveGS)
    )
  
  observeEvent(input$plot_click, {
    new_row <- nearPoints(toptable, 
                          input$plot_click,
                          xvar = "Log2FC",
                          yvar = "Log10P",
                          maxpoints = 1)
    
    values$GT <- rbind(values$GT, new_row)
    values$GT <- values$GT[!( duplicated(values$GT) | duplicated(values$GT, fromLast = T) ), ] #Removing all instances of duplicate rows disables labels on second click
  })
  
  # output$click_info <- renderPrint({
  
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
    # cat("input$plot_click:\n")
    # str(input$plot_click)
  })
  
  output$volcanoPlot <- renderPlot({
    ggplot(toptable, aes(x=Log2FC, y=Log10P, label = Gene, color = Category)) + th + 
      geom_point(
        size = 6, #Change to temput var
        alpha = 0.5,
        shape = 19
      ) +
      scale_color_manual(values = colorMap) +
      geom_text_repel(
                data = subset(toptable, Gene %in% values$GT$Gene),
                aes(Log2FC, Log10P, label = Gene, color = "textCol"),
                max.overlaps = 250,
                min.segment.length = 0.2,
                force_pull = 0.01,
                segment.colour = "black",
                segment.size = 1,
                size = 8
      )
  })
  
  output$geneTable <- renderDataTable({
    data.frame(Gene = values$GT$Gene,
               Log2FC = values$GT$Log2FC %>% round(2),
               Log10P = values$GT$Log10P %>% round(0),
               Condition = values$GT$Category)
  })
  
}

shinyApp(ui = ui, server = server)


