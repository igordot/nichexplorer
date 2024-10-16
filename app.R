suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinycssloaders)
  library(dplyr)
  library(glue)
  library(Matrix)
  library(ggplot2)
  library(cowplot)
  library(markdown)
})

options(spinner.type = 6)
options(spinner.color = "#ababab")
theme_set(theme_cowplot())

# source data ----
# cluster/sample info per cell - tibble
meta_tbl <- readRDS("data.meta.rds")
# log-transformed expression values - sparse matrix
exp_mat <- readRDS("data.exp.rds")

# clusters color palette ----
# tableau_color_pal("tableau10") in ggthemes 3, but tableau_color_pal("Classic 10") in ggthemes 4
colors_v <- c("#ff7f0e", "#1f77b4")
colors_v <- setNames(colors_v, c("V1", "V2"))
colors_p <- c("#d62728", "#2ca02c", "#8c564b", "#e377c2", "#d95f02")
colors_p <- setNames(colors_p, c("P1", "P2", "P3", "P4", "P5"))
colors_o <- c("#7f7f7f", "#1b9e77", "#bcbd22")
colors_o <- setNames(colors_o, c("O1", "O2", "O3"))
colors_c <- c("#e7298a")
colors_c <- setNames(colors_c, c("C"))
colors_clusters <- c(colors_v, colors_p, colors_o, colors_c)
colors_clusters_long <- colors_clusters[levels(meta_tbl$cluster)]
names(colors_clusters_long) <- levels(meta_tbl$cluster_long)

# check that the cluster order is the same in the meta data table and the colors vector
clusters_ordered <- meta_tbl %>% pull(cluster) %>% levels()
if (!identical(clusters_ordered, names(colors_clusters))) stop("unexpected color order")

# ui: define UI for dataset viewer app ----
ui <- fluidPage(

  # layout: header ----
  tags$head(includeHTML("gtag.html")),
  theme = shinytheme("paper"),

  # layout: title ----
  titlePanel("nichExplorer"),

  # layout: line break ----
  hr(),

  # layout: main (data) row of content ----
  fluidRow(

    # layout: left (gene selector) panel ----
    column(
      width = 3,
      # input: cell treatment group selector
      radioButtons(
        inputId = "in_treatment",
        label = "cells:",
        choices = c(
          "steady state only" = "ctrl",
          "steady state and treated" = "all",
          "split by treatment" = "split"
        )
      ),
      # input: gene selector
      selectInput(
        inputId = "in_gene",
        label = "gene:",
        choices = NULL
      )
    ),

    # layout: right (plots) panel ----
    column(
      width = 9,
      # layout: tabs for different plot types
      tabsetPanel(
        tabPanel(
          title = "tSNE Plot (Per Cell)",
          fluidRow(
            column(
              width = 6,
              withSpinner(plotOutput("tsne_gene_plot"))
            ),
            column(
              width = 6,
              withSpinner(plotOutput("tsne_cluster_plot"))
            )
          )
        ),
        tabPanel(
          title = "Bar Plot (Per Cluster)",
          withSpinner(plotOutput("bar_plot", height = "300px"))
        ),
        tabPanel(
          title = "Violin Plot (Per Cluster)",
          withSpinner(plotOutput("vln_plot", height = "300px"))
        )
      )
    )

  ),

  # layout: line break ----
  hr(),

  # layout: info (bottom) row of content ----
  fluidRow(

    column(
      width = 8,
      includeMarkdown("text.about.md"),
      includeMarkdown("text.abstract.md")
    ),
    column(
      width = 4,
      includeMarkdown("text.data.md")
    )

  ),

  # layout: line break ----
  hr()

)

# server: define server logic to summarize and view selected dataset ----
server <- function(input, output) {

  updateSelectizeInput(inputId = "in_gene", choices = rownames(exp_mat), selected = "Lepr", server = TRUE)

  # generate single gene expression values table ----
  exp_tbl <- reactive({
    req(input$in_gene, input$in_treatment)
    set.seed(99)
    # adjust the meta data table based on the requested cells subset
    if (input$in_treatment == "ctrl") {
      # keep only steady state cells
      meta_tbl <- meta_tbl %>% filter(treatment == "CTRL")
    } else if (input$in_treatment == "all") {
      # rename treatment since some plots will summarize by treatment
      meta_tbl <- meta_tbl %>% mutate(treatment = ".")
    }
    tibble(cell = colnames(exp_mat), exp_log = exp_mat[input$in_gene, ]) %>%
      inner_join(meta_tbl, by = "cell") %>%
      sample_frac()
  })

  # generate tSNE plots ----
  # expression and clusters generated together so the contents can be aligned
  tsne_plot_reactive <- reactive({
    req(input$in_gene)
    tsne_gene_plot <-
      ggplot(exp_tbl(), aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(aes(color = exp_log), size = 1) +
      labs(title = paste("Gene Expression:", input$in_gene)) +
      guides(color = guide_colorbar(title = "Expr.\nLevel\n(Log)")) +
      theme(
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank()
      ) +
      scale_color_gradientn(colors = c("gray85", "red2"))
    tsne_cluster_plot <-
      ggplot(exp_tbl(), aes(x = tSNE_1, y = tSNE_2)) +
      geom_point(aes(color = cluster_long), size = 1, show.legend = TRUE) +
      labs(title = "Clusters") +
      theme(
        plot.title = element_text(hjust = 0.5),
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank()
      ) +
      guides(color = guide_legend(title = "Cluster", override.aes = list(size = 5))) +
      scale_color_manual(values = colors_clusters_long)
    align_plots(tsne_gene_plot, tsne_cluster_plot, align = "hv", axis = "tblr")
  })

  # generate bar plot ----
  bar_plot_reactive <- reactive({
    req(input$in_gene, input$in_treatment)
    # summarize per cluster
    exp_avg_tbl <-
      exp_tbl() %>%
      mutate(exp_norm = expm1(exp_log)) %>%
      group_by(cluster, treatment) %>%
      summarize(
        num_cells = n(),
        avg_exp_norm = mean(exp_norm),
        avg_exp_log = mean(exp_log),
        std_dev_norm = sd(exp_norm),
        std_dev_log = sd(exp_log)
      ) %>%
      mutate(
        std_err_norm = std_dev_norm / sqrt(num_cells),
        std_err_log = std_dev_log / sqrt(num_cells)
      )

    # manually set the y-axis limit to prevent cutting off top bar and have x-axis cross at 0
    y_limit <- exp_avg_tbl %>% mutate(max_val = avg_exp_norm + std_err_norm) %>% pull(max_val) %>% max()
    y_limit <- y_limit * 1.05

    # generate the plot
    bar_plot <-
      ggplot(exp_avg_tbl, aes(x = cluster, y = avg_exp_norm)) +
      geom_col(aes(fill = cluster), color = "black", size = 1) +
      geom_errorbar(
        aes(ymin = avg_exp_norm - std_err_norm, ymax = avg_exp_norm + std_err_norm),
        width = 0.3, size = 1
      ) +
      labs(title = input$in_gene, x = "Cluster", y = "Norm. Expr. Level") +
      scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
      scale_fill_manual(values = colors_clusters) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing.x = unit(0.1, "lines")
      )

    # split by treatment if requested
    if (input$in_treatment == "split") {
      bar_plot <-
        bar_plot +
        facet_wrap(vars(cluster, treatment), scales = "free_x", nrow = 1, strip.position = "bottom") +
        theme(axis.text.x = element_blank())
    }

    # return the plot
    bar_plot
  })

  # generate violin plot ----
  vln_plot_reactive <- reactive({
    req(input$in_gene, input$in_treatment)

    vln_plot <-
      ggplot(exp_tbl(), aes(x = cluster, y = exp_log)) +
      geom_violin(aes(fill = cluster, color = cluster), scale = "width") +
      labs(title = input$in_gene, x = "Cluster", y = "Norm. Expr. Level (Log)") +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing.x = unit(0.1, "lines")
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = colors_clusters) +
      scale_color_manual(values = colors_clusters)

    # split by treatment if requested
    if (input$in_treatment == "split") {
      vln_plot <-
        vln_plot +
        facet_wrap(vars(cluster, treatment), scales = "free_x", nrow = 1, strip.position = "bottom") +
        theme(axis.text.x = element_blank())
    }

    # return the plot
    vln_plot
  })

  # output: tSNE plot of expression levels ----
  output$tsne_gene_plot <- renderPlot({
    # adding ggdraw due to cowplot::align_plots() output format
    ggdraw(tsne_plot_reactive()[[1]])
  })

  # output: tSNE plot of clusters ----
  output$tsne_cluster_plot <- renderPlot({
    # adding ggdraw due to cowplot::align_plots() output format
    ggdraw(tsne_plot_reactive()[[2]])
  })

  # output: bar plot ----
  output$bar_plot <- renderPlot({
    bar_plot_reactive()
  })

  # output: violin plot ----
  output$vln_plot <- renderPlot({
    vln_plot_reactive()
  })

}

# create Shiny app ----
shinyApp(ui = ui, server = server)
