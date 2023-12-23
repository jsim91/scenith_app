options(shiny.maxRequestSize = 5*1024^3)

library(shiny)
library(shinythemes)
library(scales)
library(ggplot2)
library(viridis)
library(FCSimple)
library(RColorBrewer)
library(circlize)
library(shinyBS)

ui <- navbarPage("SCENITH MM137", # Add navbarPage here
  tabPanel("UMAP", # Add tabPanel for each page
    fluidPage(
      titlePanel(""),
      sidebarLayout(
        sidebarPanel(
          fileInput("file1", "Choose .RDS File",
                    accept = c(
                      "application/x-gzip",
                      ".rds")
          ),
          numericInput("sampleSize", "Sample Size", value = 100000, min = 1),
          selectInput("type", "Color by", choices = c("cell", "cluster"), selected = "cell"),
          uiOutput("coltype"), 
          selectInput("select_trim", "Trim variable by percentile", choices = c("no trim", as.character(1:10)), selected = "no trim"),
          uiOutput("clustype"),
          selectInput("method", "Center method", choices = c("Mean", "Median", "Hodges-Lehmann"), selected = "Mean"),
          actionButton("plotButton1", "Update UMAP")
        ),
        mainPanel(
          plotOutput(outputId = "plot1")
        )
      )
    )
  ), 
  tabPanel("heatmap",
    fluidPage(
      titlePanel(""),
      sidebarLayout(
        sidebarPanel(
          numericInput("numInputs", "Number of Groups:", value = 1, min = 1),
          br(), br(), 
          bsTooltip(id = "format_button", 
            title = "groupID:1,2,3,5,8,13", 
            placement = "right", 
            trigger = "click"),
          actionButton("format_button", "Click for format guide"),
          uiOutput("stringInputs"),
          br(), 
          actionButton("plotButton2", "Update heatmap"), 
          verbatimTextOutput("result")
        ),
        mainPanel(
          plotOutput(outputId = "plot2")
        )
      )
    )
  ),
  tabPanel("plot settings", 
    fluidPage(
      titlePanel(""),
      sidebarLayout(
        sidebarPanel(
          numericInput("plot_width", "Plot width in pixels:", 900, min = 300, max = 2000),
          numericInput("plot_height", "Plot height in pixels:", 900, min = 300, max = 2000),
          numericInput("cluster_txt", "Plot annotation text size:", 8, min = 1, max = 20), 
          numericInput("axis_ti", "axis title size:", 20, min = 1, max = 40), 
          numericInput("plot_ti", "plot title size:", 20, min = 1, max = 40), 
          numericInput("legend_te", "legend text size:", 18, min = 1, max = 40),
          numericInput("legend_ti", "legend title size:", 20, min = 1, max = 40), 
          selectInput("clus_overlay", "annotate cluster numbers", choices = c("Yes","No"), selected = "Yes"), 
          selectInput("show_legend", "show legend", choices = c("Yes","No"), selected = "No"),
          selectInput("annotate_text_color", "cluster annotation text color", choices = c("black","red","white","blue"), selected = "black"), 
          numericInput("legend_key", "legend size:", 2, min = 1, max = 20)
        ),
        mainPanel(
          # Add output elements for Tab 2 here...
        )
      )
    )
  )
)

server <- function(input, output) {
  rv <- reactiveValues(cluster_centers = NULL)

  output$stringInputs <- renderUI({
    numInputs <- input$numInputs
    lapply(1:numInputs, function(i) {
      # textInput(paste0("input", i), paste0("Input ", i))
      textInput(paste0("input", i), "")
    })
  })
  observeEvent(input$plotButton2, {
    result <- lapply(1:input$numInputs, function(i) {
      pr_out <- strsplit(x = input[[paste0("input", i)]], split = "(,|\\:)")[[1]]
      grp_label <- pr_out[1]
      grp_members <- pr_out[-1]; grp_members <- gsub(" ","",grp_members)
      return(list(group_id = grp_label, group_members = grp_members))
    })
    output$result <- renderPrint(result)
  })
  #observeEvent(input$plotButton2, {
  #  stringInputsValues <- reactive({
  #    result <- lapply(1:input$numInputs, function(i) {
  #      pr_out <- strsplit(x = input[[paste0("input", i)]], split = "(,|\\:)")[[1]]
  #      grp_label <- pr_out[1]
  #      grp_members <- pr_out[-1]; grp_members <- gsub(" ","",grp_members)
  #      return(list(group_id = grp_label, group_members = grp_members))
  #    })
  #  })
  #})
  #observeEvent(input$plotButton2, {
  #  output$result <- renderPrint({
  #    stringInputsValues()
  #  })
  #})

  data <- reactive({
    req(input$file1)
    readRDS(input$file1$datapath)
    app_obj <- readRDS(input$file1$datapath)
    index_algos <- c("leiden","louvain","phenograph","flowsom")
    calc_algo <- index_algos[index_algos %in% names(app_obj)]
    for(i in seq_along(calc_algo)) {
      app_obj <- FCSimple::fcs_cluster_heatmap(fcs_join_obj = app_obj, algorithm = calc_algo[i])
    }
    return(app_obj)
  })

  observeEvent(c(input$sampleSize, input$plotButton1), {
    df <- ds_coord()
    split_df <- split(df, df$cluster)
    rv$cluster_centers <- lapply(split_df, function(data) {
      c(UMAP1 = median(data$UMAP1), UMAP2 = median(data$UMAP2))
    })
    rv$cluster_centers <- as.data.frame(do.call(rbind, rv$cluster_centers))
    rv$cluster_centers$cluster <- rownames(rv$cluster_centers)
    rownames(rv$cluster_centers) <- NULL
  })

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
  observeEvent(input$plotButton1, {
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

      if (input$coltype == "cluster") {
        plt <- ggplot(data = df, mapping = aes(x = UMAP1, y = UMAP2, color = cluster)) + 
        geom_point(size = 0.5, pch = 19, alpha = 0.2) + 
        guides(color = guide_legend(override.aes = list(size = 5, alpha = 1), title = "cluster"))
      } else {
        plt <- ggplot(data = df, mapping = aes(x = UMAP1, y = UMAP2, color = color_by)) + 
        geom_point(size = 0.5, pch = 19, alpha = 0.2) + scale_color_viridis(option = "D")
      }
      if (input$clus_overlay == "Yes") {
        plt <- plt + geom_text(data = rv$cluster_centers, mapping = aes(label = cluster, x = UMAP1, y = UMAP2), 
        size = input$cluster_txt, check_overlap = FALSE, color = input$annotate_text_color)
      }
      plt <- plt + theme_minimal() + 
      ggtitle(input$coltype) + 
      theme(axis.title = element_text(size = input$axis_ti),
            legend.text = element_text(size = input$legend_te),
            legend.title = element_text(size = input$legend_ti),
            legend.key.size = unit(input$legend_key, "line"),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5, size = input$plot_ti), 
            legend.position = ifelse(input$show_legend == "Yes", "right", "none"))
      plt
    }, width = function() { input$plot_width }, height = function() { input$plot_height })
  })
  observeEvent(input$plotButton2, {
    req(data())
    output$plot2 <- renderPlot({
      hm_pal <- viridis::viridis(n = 100,begin = 0, end = 1, alpha = 1)
      rel_block_height = 0.9
      tiles <- data()[[paste0(input$clustype,"_heatmap")]][["heatmap_tile_data"]]
      all_clus <- 0:max(data()[[input$clustype]]$clusters)
      heatmap_color_palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
      color.map.fun = circlize::colorRamp2(seq(min(tiles), max(tiles), l = n <- 100), colorRampPalette(heatmap_color_palette)(n))

      hm <- ComplexHeatmap::Heatmap(matrix = t(tiles), name = "MSE", 
                              na_col = "black", col = color.map.fun,
                              heatmap_legend_param=list(at=round(c(min(range(tiles)),max(range(tiles))),2),
                                                        labels = c("low", "high"),
                                                        legend_height=unit(10,"cm"),
                                                        grid_width=unit(0.6,"cm"),title_position="topleft",
                                                        labels_gp=gpar(fontsize=14),
                                                        title_gp=gpar(fontsize=0,fontface="bold")),
                              show_row_dend = FALSE, show_column_dend = FALSE, cluster_rows = TRUE,
                              cluster_columns = TRUE, show_heatmap_legend = FALSE,
                              column_title_gp = gpar(cex = 0, alpha = 0), column_title_side = "top",
                              column_names_side = "top",
                              column_names_rot = 0, column_names_gp = gpar(size = 12, fontface = "bold", just = "left"),
                              row_names_gp = gpar(size = 12, fontface = "bold"),
                              show_column_names = TRUE, row_names_side = "left", cluster_column_slices = FALSE)
      hm
    }, width = function() { input$plot_width }, height = function() { input$plot_height })
  })
}

shinyApp(ui = ui, server = server)
