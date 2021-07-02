
bespokeVolcano <- function(
  toptable,
  geneTable,
  title,
  pCutoff,
  FCcutoff,
  pointSize,
  colAlpha,
  xlim,
  labSize,
  drawLabel,
  colorMap) {
  
  # ---- Static Arguments ----
  

  legendPosition <- "none"
  ylim <- c(0, max(toptable[["Log10P"]], na.rm<-TRUE) + 5)
  xlab <- bquote(~Log[2]~ "fold change") #Hardcode
  ylab <- bquote(~-Log[10]~italic(P)) #Hardcode
  axisLabSize <- 18
  caption <- paste0('total <- ', nrow(toptable), ' variables')
  titleLabSize <- 18
  captionLabSize <- 14
  cutoffLineType <- 'longdash'
  cutoffLineCol <- 'black'
  cutoffLineWidth <- 0.4
  labCol <- 'black'
  labFace <- 'plain'
  shape <- 19
  legendLabSize <- 8
  subtitleLabSize <- 8
  legendLabels <- c('NS', expression(Log[2]~FC),
                   'p-value', expression(p-value~and~log[2]~FC))      #Don't Need
  legendIconSize <- 5.0     #Don't Need
  legendDropLevels <- TRUE    #Don't Need
  encircle <- NULL
  encircleCol <- 'black'
  encircleFill <- 'pink'
  encircleAlpha <- 3/4
  encircleSize <- 2.5
  shade <- NULL
  shadeFill <- 'grey'
  shadeAlpha <- 1/2
  shadeSize <- 0.01
  shadeBins <- 2
  drawConnectors <- TRUE
  widthConnectors <- 0.5
  typeConnectors <- 'closed'
  endsConnectors <- 'first'
  lengthConnectors <- unit(0.02, 'npc')
  colConnectors <- 'grey10'
  maxoverlapsConnectors <- 150
  directionConnectors <- 'both'
  arrowheads <- FALSE
  hline <- NULL
  hlineType <- 'longdash'
  hlineCol <- 'black'
  hlineWidth <- 0.4
  vline <- NULL
  vlineType <- 'longdash'
  vlineCol <- 'black'
  vlineWidth <- 0.4
  gridlines.major <- TRUE
  gridlines.minor <- TRUE
  border <- 'partial'
  borderWidth <- 0.8
  borderColour <- 'black'
  raster <- FALSE
  
  
  
  # ---- Theme ----
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
      
      legend.position = legendPosition
      )
  
  # ---- Plot ----
  plot <- ggplot(toptable, aes(x=Log2FC, y=Log10P, label = Gene, color = Category)) + th +
    
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
    
    # include new colour encodings as aes.
    # 'shape' is included, but outside aes
    geom_point(
      alpha = colAlpha,
      shape = shape,
      size = pointSize,
      na.rm = TRUE) +
    
    #This allows simple mapping via a (reactive) list
    scale_colour_manual(
      values = colorMap
    ) +
    
    # 'shape' is not included as aes. Specifying guide = TRUE
    # here will result in legends merging
    # This was cribbed from Enhanced Volcano I have no idea what this is doing
    scale_shape_manual(guide = TRUE) 
  
  # Adds labels as text (if you want boxes behind the labels, use geom_label_repel)
  # Labels are constructed by subsetting toptable based on genes present in the reactive table 
  if(drawLabel) {
    plot <- plot +  geom_text_repel(
      data = subset(toptable, Gene %in% geneTable$Gene),
      aes(Log2FC, Log10P, label = Gene, color = "textCol"),
      max.overlaps = 250,
      min.segment.length = 0.2,
      force_pull = 0.01,
      segment.colour = "black",
      segment.size = 1,
      size = labSize
    )
  }
  
  
  # ---- Labels, Borders, and Captions ----
  plot <- plot +
    
    xlab(xlab) +
    ylab(ylab) +
    
    xlim(-xlim, xlim) +
    ylim(ylim[1], ylim[2]) +
    
    geom_vline(xintercept = c(-FCcutoff, FCcutoff),
               linetype = cutoffLineType,
               colour = cutoffLineCol,
               size = cutoffLineWidth) +
    
    geom_hline(yintercept = pCutoff,
               linetype = cutoffLineType,
               colour = cutoffLineCol,
               size = cutoffLineWidth) +
    
    labs(title = title, 
         caption = caption) +
    
    theme(axis.line = element_line(size = borderWidth, colour = borderColour),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  return(plot)
}




