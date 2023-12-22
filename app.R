# source("J:/Mashayekhi/flow/MM137/scripts/scenith_visuals_ui.R") # source ui
# source("J:/Mashayekhi/flow/MM137/scripts/scenith_visuals_server.R") # source server
# shinyApp(ui = ui, server = server)

options(shiny.maxRequestSize = 5*1024^3)

library(shiny)
library(shinythemes)
library(scales)
library(ggplot2)
library(viridis)

ui <- navbarPage("SCENITH MM137", # Add navbarPage here
  tabPanel("Tab 1", # Add tabPanel for each page
    fluidPage(
      titlePanel(""),
      sidebarLayout(
        sidebarPanel(
          fileInput("file1", "Choose .RDS File",
                    accept = c(
                      "application/x-gzip",
                      ".rds")
          ),
          numericInput("sampleSize", "Sample Size", value = 10000, min = 1),
          selectInput("type", "Color by", choices = c("cluster", "cell"), selected = "cell"),
          uiOutput("coltype"), # Use uiOutput here
          selectInput("select_trim", "Trim variable by percentile", choices = c("no trim", as.character(1:10)), selected = "no trim"),
          uiOutput("clustype"),
          # selectInput("clustering_algo", "Clustering algorithm", choices = c("leiden","louvain"), selected = "leiden"),
          selectInput("method", "Center method", choices = c("Mean", "Median", "Hodges-Lehmann"), selected = "Mean"),
          actionButton("plotButton", "Update Plot")
          # Add more UI elements here...
        ),
        mainPanel(
          plotOutput("plot1"),
          verbatimTextOutput("matrixHead")#,
          #verbatimTextOutput("mapHead")#,
          #verbatimTextOutput("meanValues")
          # Add more output elements here...
        )
      )
    )
    # Add more tabPanels here...
  ), 
  tabPanel("Tab 2", 
    fluidPage(
      titlePanel("Tab 2"),
      sidebarLayout(
        sidebarPanel(
          # Add UI elements for Tab 2 here...
        ),
        mainPanel(
          # Add output elements for Tab 2 here...
        )
      )
    )
  )
)

server <- function(input, output) {
  data <- reactive({
    req(input$file1)
    readRDS(input$file1$datapath)
  })
  #recorded_data <- reactiveValues(data = NULL)
  #observe({
  #  recorded_data$data <- data()$data
  #  recorded_data$mean <- colMeans(recorded_data$data, na.rm = TRUE)
  #})
  #output$meanValues <- renderPrint({
  #  req(recorded_data$mean)
  #   recorded_data$mean
  #})
  output$coltype <- renderUI({ # Use renderUI to create the selectInput
    req(data())
    cnam <- colnames(data()$data)
    selectInput("coltype", "Variable", choices = c(cnam,"cluster"))
  })
  output$clustype <- renderUI({ # Use renderUI to create the selectInput
    req(data())
    index_algos <- c("leiden","louvain","phenograph","flowsom")
    clus_choices <- names(data())
    selectInput("clustype", "Clustering algorithm", choices = clus_choices[which(clus_choices %in% index_algos)])
  })
  output$matrixHead <- renderPrint({
    req(data())
    head(data()$data)
  })
  ds_coord <- reactive({
    req(data())
    map_coord <- data()$umap$coordinates
    map_coord$cluster <- as.character(data()[[input$clustype]]$clusters)
    if(input$coltype != "cluster"){
      map_coord$color_by <- data()$data[,input$coltype]
    } else {
      map_coord$color_by <- factor(map_coord$cluster)
    }
    set.seed(123)
    sample_iloc <- sample(seq_len(nrow(map_coord)), size = min(nrow(map_coord), input$sampleSize), replace = FALSE)
    map_coord[sample_iloc,]
  })
  output$mapHead <- renderPrint({
    #req(data())
    head(ds_coord())
  })
  calculated_value <- reactive({
    req(data())
    if (input$method == "Hodges-Lehmann") {
      # Calculate Hodges-Lehmann estimator...
      # Replace this with your own code...
      wilcox.test(ds_coord())$estimate
    } else if (input$method == "Median") {
      median(ds_coord())
    } else if (input$method == "Mean") {
      mean(ds_coord())
    }
  })
  #output$plot1 <- renderPlot({
    # initial plot
    # plot(ds_coord()[,1:2], pch = 19, col = scales::alpha("black", 0.1))
  #  plot(1:10,1:10)
  #})
  observeEvent(input$plotButton, {
    req(data())
    output$plot1 <- renderPlot({
      df <- ds_coord()
      if (input$select_trim != "no trim" && input$type == "cell" && input$coltype != "cluster") {
        st <- as.numeric(input$select_trim)
        trim_val <- quantile(x = df[,"color_by"], probs = c(st/100, 1-(st/100)))
        df <- df[df$color_by >= trim_val[1] & df$color_by <= trim_val[2], ]
      }
      if (input$type == "cluster" && input$coltype != "cluster") {
        unique_clusters <- unique(df$cluster); it_clus <- rep()
        df$color_by <- kmeans(df[, c("UMAP1", "UMAP2")], centers = 3)$cluster
      }
      plt <- ggplot(data = df, mapping = aes(x = UMAP1, y = UMAP2, color = color_by)) + 
      geom_point(size = 0.5, pch = 19, alpha = 0.2)
      if (input$coltype == "cluster") {
        plt <- plt + guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
        theme(legend.position = "right", 
              legend.text = element_text(size = 14), 
              legend.title = element_text(size = 16))
      } else {
        plt <- plt + scale_color_viridis(option = "D") + 
        theme(legend.position = "right")
      }
      plt <- plt + theme_minimal() + 
      ggtitle(input$coltype) + 
      theme(axis.title = element_text(size = 16), 
            axis.text = element_blank(), 
            plot.title = element_text(hjust = 0.5, size = 20))
      plt
    })
  })
}

shinyApp(ui = ui, server = server)
