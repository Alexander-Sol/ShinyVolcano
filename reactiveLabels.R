#Bespoke Volcano
library(ggplot2)
library(dplyr)
library(purrr)

#--- Define inputs ----
input <- list()
input$dataset <- datasets$Day_80
input$logpval <- 25
input$log2FC <- 1
input$IWP2col <- "#B856D7"
input$CHIRcol <- "#55A0FB"

reactiveGS <- c("NR2F1", "GATA3", "INSM1", "ZNF503", "FGF8", "GNG8", "LFNG", "FGFR3", "LGR5", "RPRM",
                         "CD164L2", "ZBBX", "TEKT1", "SKOR1", "AMPD3", "VEPH1")

#---- create table ----

toptable <- data.frame(
  Log2FC = input$dataset$log2FoldChange,
  Log10P = ifelse(input$dataset$padj == 0, .Machine$double.xmin, input$dataset$padj) %>% 
    log10() %>% 
    "*"(-1),
  Gene = input$dataset$gene
)

#Names the categories/colors
toptable$Category <- ifelse(
  toptable$Log2FC < input$log2FC*-1 & toptable$Log10P > input$logpval,
  "IWP2",
  ifelse(
    toptable$Log2FC > input$log2FC & toptable$Log10P > input$logpval,
    "CHIR",
    "Below Threshold"  
  )
) 


#Map colors to category, needed for scale_color_manual
#Probably colors should be reactive. something like
colorMap <- reactive({
  list(
    IWP2 = input$IWP2col,
    CHIR = input$CHIRcol,
    `Below Threshold` = 'lightgrey'
  )
})
colorMap <- list(
                 IWP2 = input$IWP2col,
                 CHIR = input$CHIRcol,
                 `Below Threshold` = 'lightgrey',
                 textCol = "black"
                )

# TODO: make genes reactive

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



#---- Plot ----  
plot <- ggplot(toptable, aes(x=Log2FC, y=Log10P, label = Gene, color = Category)) + th + 
  geom_point(
    size = 2, #Change to input var
    alpha = 0.5,
    shape = 19
   ) +
  scale_color_manual(values = colorMap) +
  geom_text(data = subset(toptable, Gene %in% reactiveGS),
            aes(Log2FC, Log10P, label = Gene, color = "textCol"))



# over-ride legend icon sizes for colour and shape.
# guide_legends are separate for colour and shape;
# so, legends will be drawn separate IF shape is also
# included as aes to geom_point (it is not, here)
guides(
  colour = guide_legend(
    order = 1,
    override.aes = list(
      size = legendIconSize)),
  shape = guide_legend(
    order = 2,
    override.aes = list(
      size = legendIconSize))) +
