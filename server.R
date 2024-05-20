options(shiny.maxRequestSize=100*1024^2)
# options(repos = BiocManager::repositories())
# getOption("repos")

function(input, output, session) {
  volumes = getVolumes()
  shinyFileChoose(input, "counts", roots = volumes, session = session)
  shinyFileChoose(input, "metadata", roots = volumes, session = session)
  shinyFileChoose(input, "image", roots = volumes, session = session)
  observeEvent(input$Run_EnSDD, {
    if(!is.null(input$counts) & !is.null(input$metadata) & !is.null(input$image)){
      # browser()
      file_selected_counts <-parseFilePaths(volumes, input$counts)
      file_selected_metadata <-parseFilePaths(volumes, input$metadata)
      file_selected_image <-parseFilePaths(volumes, input$image)
      # output$txt_file <- renderText(as.character(file_selected$datapath))
    }

    reticulate::use_python(input$python_env, require = T)
    reticulate::py_config()

    Seurat.data <- data_process(counts_path = file_selected_counts$datapath, loc_path = file_selected_metadata$datapath,
                                img_path = file_selected_image$datapath, n_HVG = input$n_HVG, n_PCA = input$n_PCA)
    load("~/Documents/R package/Shiny/EnSDD-shiny/res.clustering.RData")
    time.used = system.time({
      res = try(EnSDD::run_individual_cluster(counts_path = file_selected_counts$datapath,
                                              loc_path = file_selected_metadata$datapath, img_path = file_selected_image$datapath,
                                              python_env = input$python_env,
                                              n_HVG = input$n_HVG, n_PCA = input$n_PCA,
                                              number.setting.by.user = as.integer(input$clusters),
                                              saving_results = input$saving_results,
                                              BayesSpace = input$BayesSpace,
                                              DR.SC = input$DR_SC, SpaGCN = input$SpaGCN, stLearn = input$stLearn,
                                              GraphST = input$GraphST, STAGATE = input$STAGATE,
                                              BayesSpace_nrep = input$BayesSpace_nrep,
                                              DR.SC_n_SVGs = input$DR.SC_n_SVGs, DR.SC_latent_q = input$DR.SC_latent_q,
                                              SpaGCN_beta = input$SpaGCN_beta, SpaGCN_alpha = input$SpaGCN_alpha,
                                              SpaGCN_p = input$SpaGCN_p, SpaGCN_l = input$SpaGCN_l,
                                              SpaGCN_lr = input$SpaGCN_lr, SpaGCN_epoches = input$SpaGCN_epoches,
                                              stLearn_n.PCs = input$stLearn_n.PC,
                                              stLearn_normalize_type = input$stLearn_normalize_type,
                                              GraphST_lambda1 = input$GraphST_lambda1, GraphST_lambda2 = input$GraphST_lambda2,
                                              GraphST_tool = input$GraphST_tool, GraphST_radius = input$GraphST_radius,
                                              GraphST_refinement = input$GraphST_refinement, GraphST_n_PCs = input$GraphST_n_PCs,
                                              STAGATE_alpha = input$STAGATE_alpha))})

    output$summary1 = renderText("Run base clustering methods is finished")
    output$summary2 = renderText({
      paste("Time used (seconds):", round(as.numeric(time.used[3]),2))
    })


    time.used.ensemble = system.time({
      consensus_matrix = try(EnSDD::solve_ensemble(Results.clustering = res[[2]], lambda = NULL, prob.quantile = input$en_pro,
                                                   niter = input$en_niter, epsilon = input$en_epi,
                                                   parallel = input$en_parallel, ncors.used = input$en_cores))})
    output$summary3 = renderText("Genrating a consensus matrix based on multiple binary similarity matrix")
    output$summary4 = renderText({
      paste("Time used (seconds):", round(as.numeric(time.used.ensemble[3]),2))
    })

    time.used.clutering = system.time({
      EnSDD_clustering = try(EnSDD::Louvain_clustering(Seurat.data, similarity_mat = consensus_matrix[[1]],
                                                       setting_k = input$clusters, max.res = 5,
                                                       resolution_louvain = NULL, increase = input$res_increase))})

    output$summary5 = renderText("Spatial domain assignment")
    output$summary6 = renderText({
      paste("Time used (seconds):", round(as.numeric(time.used.clutering[3]),2))
    })


    meta_data <- Seurat.data@meta.data
    res.vec <- res[[1]]
    res.vec.df <- as.data.frame(res.vec)
    res.vec.df['EnSDD'] <- EnSDD_clustering[rownames(Seurat.data@meta.data)]
    meta_data_all <- cbind(meta_data, res.vec.df)

    methods_names <- c(names(res[[1]]), "EnSDD")

    write.table(meta_data_all, file = "res.clustering.txt", sep = "\t")
    ### running the code for the plotting
    # Seurat.data.plot <- Seurat.data
    # meta.data.all.plot <- meta_data_all
    reticulate::source_python(system.file("python", "Visualization_HBC_clustering_shiny.py", package = "EnSDD"))
    # cat("point1")
    visualization_clustering_shiny(counts_path = file_selected_counts$datapath,
                                   meta_path = "res.clustering.txt",
                                   img_path = file_selected_image$datapath,
                                   cluster_num = as.integer(input$clusters),
                                   spot_diameter_fullres = 150, spot_size = 1.5)
    unlink("res.clustering.txt")
    # cat("point2")
    output$clustering_res_EnSDD <- renderPlot({
      img_EnSDD <- readPNG("EnSDD.png")
      print(grid::grid.raster(img_EnSDD))
    })

    output$clustering_res_base_methods <- renderPlot({
      img_base_methods <- readPNG("base_methods.png")
      print(grid::grid.raster(img_base_methods))
    })


    # unlink("res.clustering.txt")

    output$result_EnSDD <- DT::renderDataTable({

      DT::datatable(meta_data_all,  extensions = 'Scroller', options = list(
        deferRender = TRUE,
        scrollY = 296,
        scroller = TRUE
      ))

    })

    output$download = downloadHandler(

      filename = function(){ paste(input$dlname, ".", input$fileformat, sep="") },

      content = function(file) {
        if (input$fileformat == "txt"){
          write.table(meta_data_all, file, sep = "\t", quote = FALSE)
        }
        else{
          write.csv(meta_data_all, file,  quote = FALSE)
        }
      }
    )

  })

}
