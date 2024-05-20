# ui.R
library(shiny)
library(shinyBS)
library(shinyFiles)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(graphics)
library(reticulate)
## used for EnSDD
library(EnSDD)
library(png)
options(shiny.maxRequestSize=100*1024^2)
shinyUI(navbarPage("EnSDD", id="navbar",
                   theme = shinytheme("flatly"),
                   fluidPage(
                     sidebarLayout(
                       sidebarPanel(

                         ###Parameters of twelve individual methods
                         # fileInput("count",
                         #           label = "Count matrix",
                         #           multiple = FALSE,
                         #           accept = c("text/csv",
                         #                      "text/comma-separated-values,text/plain",
                         #                      ".csv")),
                         shinyFilesButton("counts", "Spatial gene expression" ,
                                          title = "Please select a file:", multiple = FALSE,
                                          buttonType = "default", class = NULL),
                         shinyFilesButton("metadata", "Spatial location" ,
                                          title = "Please select a file:", multiple = FALSE,
                                          buttonType = "default", class = NULL),
                         shinyFilesButton("image", "H&E image" ,
                                          title = "Please select a file:", multiple = FALSE,
                                          buttonType = "default", class = NULL),

                         # fileInput("metadata",
                         #           label = "Spatial location",
                         #           multiple = FALSE,
                         #           accept = c("text/csv",
                         #                      "text/comma-separated-values,text/plain",
                         #                      ".csv")),
                         # fileInput("image",
                         #           label = "H&E image",
                         #           multiple = FALSE,
                         #           accept = c(".tif")),
                         textInput("python_env",
                                   label = "python environment",
                                   value = "~/.conda/envs/EnSDD/bin/python"),
                         numericInput("n_HVG", "Number of HVGs", value = 2000, step = 20, min = 1000),
                         numericInput("n_PCA", "Number of PCAs", value = 30, step = 5, min = 10),
                         radioButtons("saving_results", "Saving results?", c("TRUE", "FALSE"), inline = FALSE),

                         ### ensemble parameters
                         h4("Ensemble learning parameters"),
                         numericInput("en_pro", "Pro quantile", value = 0.5, min = 0, max = 1, step = 0.1),
                         numericInput("en_niter", "Iteration", value = 100, min = 0, max = 200, step = 1),
                         numericInput("en_epi", "Epsilon", value = 1e-5, min = 0, max = 1e-3, step = 1e-5),
                         radioButtons("en_parallel", "Parallel", c("TRUE", "FALSE"), inline = TRUE),
                         conditionalPanel("input.en_parallel === 'TRUE'",
                                          numericInput("en_cores", label = "Cores", value = 32, min =2, max = 32, step = 1)),
                         numericInput("res_increase", "Resolution interval", value = 0.05, min = 0.01, max = 0.1, step = 0.01),

                         #### base clustering parameters
                         h4("base clustering parameters"),
                         numericInput("clusters", label = "clusters",
                                      value = 7, min = 1, max = 50, step = 1),
                         h4("SDD parameters"),
                         helpText("The following parameters are default setting for individual SDD methods."),

                         # checkboxInput("BayesSpace", "BayesSpace", value = TRUE),

                         #BayesSpace
                         radioButtons("BayesSpace", "BayesSpace", c("TRUE", "FALSE"), inline = TRUE),
                         conditionalPanel("input.BayesSpace === 'FALSE'",
                                          helpText("The BayesSpace will not running for the ensemble learning")),
                         conditionalPanel("input.BayesSpace === 'TRUE'",
                                          numericInput("BayesSpace_nrep", label = "BayesSpace MCMC", value = 50000, step = 100)),

                         #DR.SC

                         radioButtons("DR_SC", "DR.SC", c("TRUE", "FALSE"), inline = TRUE),
                         conditionalPanel("input.DR_SC === 'FALSE'",
                                          helpText("The DR.SC will not running for the ensemble learning")),
                         conditionalPanel("input.DR_SC === 'TRUE'",
                                          numericInput("DR.SC_n_SVGs", "DR.SC n SVGs", value = 600, step = 20),
                                          numericInput("DR.SC_latent_q", "DR.SC latent q", value = 15, min = 5,step = 5)),


                         #SpaGCN
                         radioButtons("SpaGCN", "SpaGCN", c("TRUE", "FALSE"), inline = TRUE),
                         conditionalPanel("input.SpaGCN === 'FALSE'",
                                          helpText("The SpaGCN will not running for the ensemble learning")),
                         conditionalPanel("input.SpaGCN === 'TRUE'",
                                          numericInput("SpaGCN_beta", "SpaGCN beta", value = 55, min = 0, step = 5),
                                          numericInput("SpaGCN_alpha", "SpaGCN alpha", value = 1, min = 0, max = 1, step = 0.1),
                                          numericInput("SpaGCN_p", "SpaGCN p", value = 0.5, min = 0, max = 1, step = 1),
                                          numericInput("SpaGCN_l", "SpaGCN l", value = 0.5, min = 0, max = 1, step = 1),
                                          numericInput("SpaGCN_lr", "SpaGCN lr", value = 0.05, min = 0.01, max = 0.1, step = 0.01),
                                          numericInput("SpaGCN_epoches", "SpaGCN epoches", value = 200, min = 100, max = 500, step = 10)),


                         #stLearn
                         radioButtons("stLearn", "stLearn", c("TRUE", "FALSE"), inline = TRUE),
                         conditionalPanel("input.stLearn === 'FALSE'",
                                          helpText("The stLearn will not running for the ensemble learning")),
                         conditionalPanel("input.stLearn === 'TRUE'",
                                          numericInput("stLearn_n.PC", "stLearn PC", value = 30, min = 10, max = 50, step = 5),
                                          numericInput("stLearn_res", label = "stLearn resolution",
                                                       value = 0.7, min = 0.1, max = 5, step = 0.1),
                                          numericInput("stLearn_n_neig", "stLearn neighboors", value = 15, min = 10, max = 30, step = 1),
                                          selectInput('stLearn_normalize_type', 'stLearn normalization', c("weights_matrix_all","weights_matrix_pd_gd",
                                                                                                           "weights_matrix_pd_md",
                                                                                                           "weights_matrix_gd_md",
                                                                                                           "gene_expression_correlation",
                                                                                                           "physical_distance",
                                                                                                           "morphological_distance"), selected = "physical_distance"),
                                          helpText("spatial location (S), tissue morphological feature (M) and gene expression (E)")),

                         #GraphST
                         radioButtons("GraphST", "GraphST", c("TRUE", "FALSE"), inline = TRUE),
                         conditionalPanel("input.GraphST === 'FALSE'",
                                          helpText("The GraphST will not running for the ensemble learning")),
                         conditionalPanel("input.GraphST === 'TRUE'",
                                          numericInput("GraphST_lambda1", "GraphST lambda1", value = 10, min = 0, max = 20, step = 0.5),
                                          numericInput("GraphST_lambda2", "GraphST lambda2", value = 1, min = 0, max = 5, step = 0.1),
                                          selectInput("GraphST_tool", "GraphST tool", c('leiden', 'louvain', 'mclust'), selected = "mclust"),
                                          numericInput("GraphST_radius", "GraphST radius", value = 50, min = 10, max = 100, step = 5),
                                          checkboxInput("GraphST_refinement", "GraphST refinement", value = FALSE),
                                          numericInput("GraphST_n_PCs", "GraphST PC", value = 20, min = 10, max = 50, step = 5)),

                         #STAGATE
                         radioButtons("STAGATE", "STAGATE", c("TRUE", "FALSE"), inline = TRUE),
                         conditionalPanel("input.STAGATE === 'FALSE'",
                                          helpText("The STAGATE will not running for the ensemble learning")),
                         conditionalPanel("input.STAGATE === 'TRUE'",
                                          numericInput("STAGATE_alpha", "STAGATE alpha", value = 0, min = 0, max = 1, step = 0.1)),




                         actionButton("Run_EnSDD",
                                      label = "Run EnSDD",
                                      style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),

                         width = 4
                       ),
                       mainPanel(
                         fluidRow(
                           column(9,
                                  plotOutput("clustering_res_EnSDD"),
                                  br(),
                                  plotOutput("clustering_res_base_methods"),
                                  br(),
                                  textOutput("summary1"),
                                  br(),
                                  textOutput("summary2"),
                                  br(),
                                  textOutput("summary3"),
                                  br(),
                                  textOutput("summary4"),
                                  br(),
                                  textOutput("summary5"),
                                  br(),
                                  textOutput("summary6")
                           ),
                           column(3,
                                  radioButtons("fileformat",
                                               label = "Output File Format",
                                               choices = c(".txt (tab-delimited text)" = "txt",
                                                           ".csv (comma-separated values)" = "csv"),
                                               selected = "csv"),
                                  textInput("dlname",
                                            label = "Output File Name (Do not include file extension)"),
                                  downloadButton("download",
                                                 label = "Download")
                           )
                         ),
                         fluidRow(
                           DT::dataTableOutput("result_EnSDD")),

                       )
                     )
                   )
)
)
