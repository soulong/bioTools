
function (plot, object = NULL, ident = "SelectedCells", ...) {
  ui <- miniPage(gadgetTitleBar(title = "Cell Selector", 
                                left = miniTitleBarButton(inputId = "reset", label = "Reset")), 
                 miniContentPanel(
                   fillRow(
                     plotOutput(outputId = "plot", nheight = "100%", 
                                brush = brushOpts(id = "brush", delay = 100, 
                                                  delayType = "debounce", clip = TRUE, 
                                                  resetOnNew = FALSE)
                                )
                     ),
                   )
                 )
  
  if (inherits(x = plot, what = "patchwork")) {
    if (length(x = plot$patches$plots)) {
      warning("Multiple plots passed, using last plot", 
              call. = FALSE, immediate. = TRUE)
    }
    class(x = plot) <- grep(pattern = "patchwork", x = class(x = plot), 
                            value = TRUE, invert = TRUE)
  }
  
  xy.aes <- GetXYAesthetics(plot = plot)
  dark.theme <- !is.null(x = plot$theme$plot.background$fill) && 
    plot$theme$plot.background$fill == "black"
  plot.data <- GGpointToBase(plot = plot, do.plot = FALSE)
  plot.data$selected_ <- FALSE
  rownames(x = plot.data) <- rownames(x = plot$data)
  colnames(x = plot.data) <- gsub(pattern = "-", replacement = ".", 
                                  x = colnames(x = plot.data))
  
  server <- function(input, output, session) {
    plot.env <- reactiveValues(data = plot.data)
    observeEvent(eventExpr = input$done, handlerExpr = {
      PlotBuild(data = plot.env$data, dark.theme = dark.theme)
      selected <- rownames(x = plot.data)[plot.env$data$selected_]
      if (inherits(x = object, what = "Seurat")) {
        if (!all(selected %in% Cells(x = object))) {
          stop("Cannot find the selected cells in the Seurat object, please be sure you pass the same object used to generate the plot")
        }
        Idents(object = object, cells = selected) <- ident
        selected <- object
      }
      stopApp(returnValue = selected)
    })
    observeEvent(eventExpr = input$reset, handlerExpr = {
      plot.env$data <- plot.data
      session$resetBrush(brushId = "brush")
    })
    observeEvent(eventExpr = input$brush, handlerExpr = {
      plot.env$data <- brushedPoints(df = plot.data, brush = input$brush, 
                                     xvar = xy.aes$x, yvar = xy.aes$y, allRows = TRUE)
      plot.env$data$color <- ifelse(test = plot.env$data$selected_, 
                                    yes = "#DE2D26", no = "#C3C3C3")
    })
    output$plot <- renderPlot(expr = PlotBuild(data = plot.env$data, 
                                               dark.theme = dark.theme))
  }
  
  return(runGadget(app = ui, server = server))
}
