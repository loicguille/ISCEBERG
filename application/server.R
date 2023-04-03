library(shiny)
library(Seurat)
library(waiter)
library(shinyjs)
library(shinyalert)
library(shinycssloaders)
library(Matrix)
library(bslib)
library(ggplot2)
library(clustree)
library(DT)
library(cowplot)
library(RColorBrewer)
library(plotly)
library(grid)
library(dplyr)
library(shinymanager)
library(gridExtra)
library(colourpicker)
library(scales)
library(rlang)
library(rmarkdown)


source("Utils.R")
source("Pooled_Genes_mining.R")
source("Extraction_information.R")
source("Annotation.R")
source("Gene_mining.R")
source("QC.R")
source("ReadFile.R")
source("Automatic_annot_SCINA.R")
source("New_SCINA.R")
source("write_output_download.R")
options(shiny.maxRequestSize =100000*1024^2)


server <- function(input, output,session) {
    hideTab(inputId = "tabs", target = "filtering")
    hideTab(inputId = "tabs", target = "Cluster tree")
    hideTab(inputId = "tabs", target = "DE between clusters")
    hideTab(inputId = "tabs", target = "Data mining for one gene")
    hideTab(inputId = "tabs", target = "Data mining for a combination of gene")
    hideTab(inputId = "tabs", target = "Extract Information")
    hideTab(inputId = "tabs", target = "Subclustering")
    hideTab(inputId = "tabs", target = "Add annotation")
    hideTab(inputId = "tabs", target = "QC")
    hideTab(inputId = "tabs", target = "Add automatic annotation")
    shinyjs::hide("Refresh")
    ##### read survey to give the rds available on the application ######
    output$Isceberg <- renderUI({
      img(src = "iceberg2.png", height = "90%", width ="90%")
    })
    output$memoryCons <- renderUI({
        img(src = "nb_giga.jpeg", height = "90%", width ="90%")
        })
    output$timeCons <- renderUI({
        img(src = "time.jpeg", height = "90%", width = "90%")
        })
    observe({
        shinyjs::toggleState("createSeurat", !is.null(input$InputFile) & input$InputFile != "" |  !is.null(input$MatrixFile) & input$MatrixFile != "" |  !is.null(input$listFiles) & input$listFiles!= "")#allow to grey a button
    })
    observeEvent(input$file,{
        if(input$file =="RDS"){
            shinyalert("Information", "If you choose to give your proper RDS dataset, this has to be an output from seurat analyzed by yourself", type = "info")
        }
    })
    observeEvent(input$Refresh,{
        shinyalert("warning", "Do this operation will reset all what you have done before", type = "warning", showCancelButton = TRUE, callbackR = function(x){if(x == TRUE){shinyjs::refresh()}})
        
    })

    
    #####################################################################################################################################################
    ############################### Create or Upload Seurat object ######################################################################################
    #####################################################################################################################################################
    observeEvent(input$createSeurat ,{
        shinyjs::show("Refresh")
        shinyjs::disable("Refresh")
        #Hide all buttons from first page so the user will need to reload the page if he want to see another file
        for(i in c("createSeurat","choiceFile","file","listFiles","header","separator","GeneFile","CellsFile","MatrixFile","InputFile","Metadata","MinGeneByCells","MinCells","download_barplot_before_after","InfoToPlotFilter","BooleanColorsFilter","downloadColorFileFilter","download_cell_cycle","download_choice_filter","download_clustering_plot","showlabelcc", "showlabelfilter","showlabelcluster","downloadAfterFiltering","download_nb_cells_dt","download_UMAP_nb_genes_cells","downloadbatchVizu","showlabelBatch", "FeatureName","Scina_annot_DL","Probabilities_Dl","Proba_selection","pointSizeNbCell","pointSizeBatch")){
          shinyjs::hide(i)
        }
        ##Loading data and create seurat object
        rdsbool <- FALSE
        withProgress(message = "Reading file",value=0,{
            if(input$file == "RDS"){
              rdsbool <- TRUE
              seuratObj <- readRDS(input$InputFile$datapath)
            }
            if(input$file == "CSV" ){
              seuratObj <- ReadCSV_files(csvPath = input$InputFile, metadataFiles = input$Metadata, header = input$header, separator = input$separator, minGenebyCells = input$MinGeneByCells , minCells = input$MinCells)
            }
            if(input$file == "Mtx" ){
              seuratObj <- readMatrix_files(matrix_files = input$MatrixFile, metadataFiles = input$Metadata, feature_files = input$GeneFile, cells_files = input$CellsFile,minGenebyCells = input$MinGeneByCells, minCells = input$MinCells, FeatureName = input$FeatureName)
            }
            if(input$file == "H5"){
              seuratObj <- readH5_files(H5_files = input$InputFile, metadataFiles = input$Metadata, minGenebyCells = input$MinGeneByCells , minCells = input$MinCells)
            }
            if(input$file == "Txt"){
              seuratObj <- Readtxt_files(txt_files = input$InputFile, metadataFiles = input$Metadata, header = input$header, separator = input$separator, minGenebyCells = input$MinGeneByCells , minCells = input$MinCells)
            }
            seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^mt-")
             output$PreProcessPlot <- renderPlot({
                makeVlnGreatAgain(seuratObj, grouping = "orig.ident",var = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
            })
            output$downloadVln<- downloadHandler(
              filename="VlnPlot_QC.tiff",
              content=function(file){
                tiff(file, width = 900 , height = 600,res = 100)
                print(makeVlnGreatAgain(seuratObj,  grouping = "orig.ident",var = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
                dev.off()
              })
            output$BeforeFiltering <- renderPlot({
              MakeComplicatedVln(seuratObj,var = "orig.ident", "QC before subsetting datasets",  minGenePerCell = input$minGenePerCells, maxGenePerCell = input$maxGenePerCells, maxCountPerCells = input$maxCountPerCells,percentMito = input$percentMito)
            })
            output$downloadBeforeFiltering<- downloadHandler(
              filename="VlnPlot_before_filt.tiff",
              content=function(file){
                tiff(file,width = 900 , height = 600,res = 100)
                print(MakeComplicatedVln(seuratObj,var = "orig.ident", "QC before subsetting datasets", minGenePerCell = input$minGenePerCells, maxGenePerCell = input$maxGenePerCells, maxCountPerCells = input$maxCountPerCells,percentMito = input$percentMito))
                dev.off()
              })
            if(rdsbool == FALSE){
                showTab(inputId = "tabs", target = "filtering")
                shinyjs::hide("downloadSeurat")
                shinyjs::hide("downloadLogFile")
                shinyjs::hide("NameSeuratLog")
                shinyjs::enable("Refresh")
            }else{
                
                showTab(inputId = "tabs", target = "QC")
                showTab(inputId = "tabs", target = "Cluster tree")
                showTab(inputId = "tabs", target = "DE between clusters")
                showTab(inputId = "tabs", target = "Data mining for one gene")
                showTab(inputId = "tabs", target = "Data mining for a combination of gene")
                showTab(inputId = "tabs", target = "Extract Information")
                
                showTab(inputId = "tabs", target ="Subclustering")
                for(i in c("downloadUMAPaftersbt", "labelsAfterSubset")){
                  shinyjs::hide(i)
                }
                for(i in c("downloadUMAPcellKept","cellsHighlightSize","cellsSize")){
                  shinyjs::hide(i)
                }
                
                showTab(inputId = "tabs", target = "Add annotation")
                for(i in c("downloadUMAP_post_annotation","labelsPostAnnot")){
                  shinyjs::hide(i)
                }
                
                showTab(inputId = "tabs", target = "Add automatic annotation")
                shinyjs::enable("Refresh")
            }
            #shinyjs::show("Refresh")
        })
        if(rdsbool == FALSE){
            observeEvent(input$runFiltNorm,{
                shinyjs::disable("runFiltNorm")
                SeuratObjsubset <- subset(seuratObj, subset = nFeature_RNA > input$minGenePerCells & nFeature_RNA < input$maxGenePerCells & percent.mt < input$percentMito & nCount_RNA < input$maxCountPerCells)
                withProgress(message = "Normalizing data", value = 0,{
                    if(input$normalization == "SCTransform"){
                        
                        SeuratObjsubset <- SCTransform(SeuratObjsubset)
                        incProgress(1)
                        
                    }else{
                        SeuratObjsubset <- NormalizeData(SeuratObjsubset)
                        incProgress(1/3)
                        SeuratObjsubset <- FindVariableFeatures(SeuratObjsubset)
                        incProgress(1/3)
                        SeuratObjsubset <- ScaleData(SeuratObjsubset)
                        incProgress(1/3)
                    }
                })
                withProgress(message = "Running PCA", value = 0,{
                    SeuratObjsubset <- RunPCA(SeuratObjsubset)
                    incProgress(1/3,message ="Finish PCA")
                    SeuratObjsubset <- FindNeighbors(SeuratObjsubset)
                    incProgress(1/3,message ="Running Umap")
                    SeuratObjsubset <- RunUMAP(SeuratObjsubset, dims = 1:30)
                    incProgress(1/3,message ="Finish Umap")
                })
                withProgress(message = "Clustering data", value = 0, {
                    inc <- (input$Resolution[2]-input$Resolution[1])/input$ResolutionStep
                    for(i in seq(input$Resolution[1],input$Resolution[2], input$ResolutionStep)){
                        SeuratObjsubset <- FindClusters(SeuratObjsubset, res = i)
                        incProgress(1/inc)
                    } 
                })
                for(i in c("downloadSeurat","downloadLogFile","NameSeuratLog","download_barplot_before_after","InfoToPlotFilter","BooleanColorsFilter","downloadColorFileFilter","download_cell_cycle","download_choice_filter","download_clustering_plot","showlabelcc", "showlabelfilter","showlabelcluster","downloadAfterFiltering","download_nb_cells_dt","download_UMAP_nb_genes_cells","downloadbatchVizu","showlabelBatch","pointSizeNbCell","pointSizeBatch")){
                  shinyjs::show(i)
                }
                observe({
                  shinyjs::toggleState("downloadColorFileFilter", input$BooleanColorsFilter != FALSE)
                })
                observe({
                    shinyjs::toggleState("downloadSeurat",input$NameSeuratLog != "")
                })
                output$downloadSeurat <- downloadHandler(
                    filename = function(){
                        paste0(input$NameSeuratLog,".rds")
                    },
                    content = function(file){
                        saveRDS(SeuratObjsubset,file)
                    }
                    
                )
                observe({
                    shinyjs::toggleState("downloadLogFile",input$NameSeuratLog != "")
                })
                output$downloadLogFile <- downloadHandler(
                    filename = function(){
                        paste0(input$NameSeuratLog,"_log.html")
                    },
                    content = function(file){
                      tempReport <- file.path(tempdir(),"logFile.Rmd")
                      file.copy("logFile.Rmd",tempReport ,overwrite = TRUE)
                      command <- list()
                      line <- paste0("SeuratObjsubset <- subset(seuratObj, subset = nFeature_RNA > ",input$minGenePerCells," & nFeature_RNA < ",input$maxGenePerCells," & nCount_RNA < ",input$maxCountPerCells," & percent.mt < ",input$percentMito,")")
                      command <- append(command,line)
                      
                      if(input$normalization == "SCTransform"){
                          line <- paste0("SeuratObjsubset <- SCTransform(SeuratObjsubset)")
                          command <- append(command,line)
                      }else{
                          line <- paste0("SeuratObjsubset <- NormalizeData(SeuratObjsubset)")
                          command <- append(command,line)
                          line <- paste0("SeuratObjsubset <- FindVariableFeatures(SeuratObjsubset)")
                          command <- append(command,line)
                          line <- paste0("SeuratObjsubset <- ScaleData(SeuratObjsubset)")
                          command <- append(command,line)
                      }
                      command <- append(command,"SeuratObjsubset <- RunPCA(SeuratObjsubset)")
                      command <- append(command,"SeuratObjsubset <- FindNeighbors(SeuratObjsubset)")
                      command <- append(command,"SeuratObjsubset <- RunUMAP(SeuratObjsubset, dims = 1:30)")
                      for(i in seq(input$Resolution[1],input$Resolution[2], input$ResolutionStep)){
                        command <- append(command,paste0("SeuratObjsubset <- FindClusters(SeuratObjsubset, res =", i,")"))
                      }
                      params <- list(use = command)
                      rmarkdown::render(tempReport,output_file = file,params = params, envir = new.env(parent = globalenv()))
                    }
                )
                
                shinyjs::enable("runFiltNorm")
                showTab(inputId = "tabs", target = "Cluster tree")
                showTab(inputId = "tabs", target = "DE between clusters")
                showTab(inputId = "tabs", target = "Data mining for one gene")
                showTab(inputId = "tabs", target = "Data mining for a combination of gene")
                showTab(inputId = "tabs", target = "Extract Information")
                showTab(inputId = "tabs", target ="Subclustering")
                showTab(inputId = "tabs", target = "Add annotation")
                showTab(inputId = "tabs", target = "Add automatic annotation")
                for(i in c("downloadUMAPaftersbt", "labelsAfterSubset")){
                  shinyjs::hide(i)
                }
                for(i in c("downloadUMAPcellKept","cellsHighlightSize","cellsSize")){
                  shinyjs::hide(i)
                }
                
                for(i in c("downloadUMAP_post_annotation","labelsPostAnnot")){
                  shinyjs::hide(i)
                }
                
                ###########################################################################################################################################
                ##########################                                 DATA MINING                              #######################################
                ###########################################################################################################################################
                available_resolution <- grep(colnames(SeuratObjsubset@meta.data), pattern = "*_snn_res.*", value = TRUE)
                list_Color_Label <- list()
                updateSelectizeInput(session, "resolution_TreePage", choices=available_resolution, server=TRUE)
                s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm7","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Polr1b","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Mrpl36","E2f8")
                g2m.genes <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Pimreg","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
                tryCatch(
                    expr = SeuratObjsubset <- CellCycleScoring(SeuratObjsubset, s.features = s.genes, g2m.features = g2m.genes),
                    error = function(c){
                        shinyalert("Warning", "Cannot calculate phase because the following feature lists do not have enough features present in the object: S.Score exiting...The following feature lists do not have enough features present in the object: G2M.Score exiting...", type = "warning")
                    }
                )
                toremove <- c(grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
                varAvail <- names(SeuratObjsubset@meta.data)
                varAvail <- varAvail[-which(varAvail %in% toremove)]
                updateSelectizeInput(session,"InfoToPlotFilter",choices = varAvail, server = TRUE)
                output$AfterFiltering <- renderPlot({
                  MakeComplicatedVln(SeuratObjsubset,var = "orig.ident", "QC after subsetting datasets", minGenePerCell = input$minGenePerCells, maxGenePerCell = input$maxGenePerCells, maxCountPerCells = input$maxCountPerCells,percentMito = input$percentMito)
                })
                output$downloadAfterFiltering<- downloadHandler(
                  filename="VlnPlot_after_filt.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(MakeComplicatedVln(SeuratObjsubset,var = "orig.ident", "QC after subsetting datasets", minGenePerCell = input$minGenePerCells, maxGenePerCell = input$maxGenePerCells, maxCountPerCells = input$maxCountPerCells,percentMito = input$percentMito))
                    dev.off()
                  })
                output$clusteringPlot <- renderPlot({
                    vizu_UMAP(SeuratObjsubset, var = NULL, dlabel = input$showlabelcluster)
                })
                output$download_clustering_plot<- downloadHandler(
                  filename="UMAP_orig_ident.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(SeuratObjsubset, var=NULL, dlabel = input$showlabelcluster))
                    dev.off()
                  })
                
                output$barplotBeforeAfter <- renderPlot({
                    Plotnbcellsbeforeafter(objBefore = seuratObj,objAfter = SeuratObjsubset)
                })
                output$download_barplot_before_after<- downloadHandler(
                  filename="plot_before_after.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(Plotnbcellsbeforeafter(objBefore = seuratObj ,objAfter = SeuratObjsubset))
                    dev.off()
                  })
                output$cellcycle <- renderPlot({
                    validate(
                        need(try(SeuratObjsubset[["Phase"]]!=""),message = "No phase available in this seurat object")
                    )
                    vizu_UMAP(SeuratObjsubset, var = "Phase", dlabel = input$showlabelcc,color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                })
                output$download_cell_cycle<- downloadHandler(
                  filename="UMAP_cell_cycle.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(SeuratObjsubset, var="Phase", dlabel = input$showlabelcc, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                    dev.off()
                  })
                output$nbcells_by_datasets <- renderPlot({
                  nbCellsbydt(SeuratObjsubset)
                })
                output$download_nb_cells_dt<- downloadHandler(
                  filename="Plot_nb_cells_dt.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(nbCellsbydt(SeuratObjsubset))
                    dev.off()
                  })
                output$geneExpByCell <- renderPlot({
                  UmapNbGenesCell(SeuratObjsubset, sizePoint = input$pointSizeNbCell)
                })
                output$download_nb_cells_dt<- downloadHandler(
                  filename="UMAP_nb_genes_by_cells.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(UmapNbGenesCell(SeuratObjsubset, sizePoint = input$pointSizeNbCell))
                    dev.off()
                  })
                output$BatchVizu <- renderPlot({
                  vizu_UMAP(SeuratObjsubset, var = "file", dlabel = input$showlabelBatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                })
                output$downloadbatchVizu<- downloadHandler(
                  filename="UMAP_batch_visu.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(SeuratObjsubset, var = "file", dlabel = input$showlabelBatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                    dev.off()
                  })
                
                
                
                ### Color changer recursively  
                listenColors <- reactive({list(input$BooleanColorsFilter, input$InfoToPlotFilter, input$submitAnnot)})
                observeEvent(listenColors(),{
                  Idents(SeuratObjsubset) <- input$InfoToPlotFilter
                  orderLevel <- sort(levels(SeuratObjsubset))
                  if(input$BooleanColorsFilter == TRUE){
                    shinyjs::show("myPanelFilter")
                    cols <- reactive({
                      lapply(seq_along(levels(SeuratObjsubset)), function(i) {
                        colourpicker::colourInput(paste("col", i, sep="_"), paste("Choose colour for",orderLevel[i]), "black")
                      })
                    })
                    output$myPanelFilter <- renderUI({cols()})
                    # Put all the input in a vector
                    colors <- reactive({
                      lapply(seq_along(levels(SeuratObjsubset)), function(i) {
                        input[[paste("col", i, sep="_")]]
                      })
                    })
                    render_colors <- function(){
                      if (is.null(input$col_1)) {
                        cols <- rep("#000000", length(levels(SeuratObjsubset)))
                      } else {
                        cols <- unlist(colors())
                        labels_object <- levels(SeuratObjsubset)
                        list_Color_Label[[input$InfoToPlotFilter]]$color <<- cols
                        list_Color_Label[[input$InfoToPlotFilter]]$label <<- labels_object
                      }
                      return(cols)
                    }
                    output$downloadColorFileFilter <- downloadHandler(
                      filename = "ColorFile.csv",
                      content=function(file){
                        write_color_file(color_list = list_Color_Label, name_file = file)
                      }
                    )
                    output$plotchoiceFilter <- renderPlot({
                      cols <- render_colors()
                      vizu_UMAP(SeuratObjsubset,var =input$InfoToPlotFilter,  dlabel = input$showlabelfilter, color = cols, sizePoint = input$pointSizeBatch)
                    })
                    output$download_choice_filter<- downloadHandler(
                      filename="UMAP_choice.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        cols <- render_colors()
                        print(vizu_UMAP(SeuratObjsubset,var=input$InfoToPlotFilter,  dlabel = input$showlabelfilter, color = cols, sizePoint = input$pointSizeBatch))
                        dev.off()
                      })
                  }else{
                    shinyjs::hide("myPanelFilter")
                    output$plotchoiceFilter <- renderPlot({
                      vizu_UMAP(SeuratObjsubset,var =input$InfoToPlotFilter,  dlabel = input$showlabelfilter, sizePoint = input$pointSizeBatch)
                    })
                    output$download_choice_filter<- downloadHandler(
                      filename="UMAP_choice.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(SeuratObjsubset,var=input$InfoToPlotFilter,  dlabel = input$showlabelfilter, sizePoint = input$pointSizeBatch))
                        dev.off()
                      })
                  }
                })
                ####################################
                ######## cluster tree Page #########
                ####################################
                output$ClusterTree <- renderPlot({
                    validate(
                        need(length(grep(colnames(SeuratObjsubset@meta.data),pattern ="*_snn_res.*")) > 1, "Cannot create cluster tree with only one cluster")
                    )
                    clustree(SeuratObjsubset,prefix = paste0(SeuratObjsubset@active.assay,"_snn_res."))+theme(legend.position = "bottom")
                })
                output$DimplotTreePage <- renderPlot({
                    vizu_UMAP(SeuratObjsubset, var = input$resolution_TreePage, dlabel = input$addlabels_tree, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                })
                output$downloadUMAP_Tree_page<- downloadHandler(
                  filename="UMAP_resolution_tree_page.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(SeuratObjsubset, var = input$resolution_TreePage,dlabel = input$addlabels_tree, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                    dev.off()
                  })
                output$nbCellsByClust_Tree_Page <- DT::renderDataTable({
                    count_by_cluster <- data.frame(table(SeuratObjsubset[[input$resolution_TreePage]]))
                    colnames(count_by_cluster) <- c("Cluster",'Number of cells')
                    DT::datatable(count_by_cluster,rownames = F)
                })
                ###########################
                ### DE between clusters ###
                ###########################
                updateSelectizeInput(session, "resolutionDE", choices=available_resolution, server=TRUE, selected = available_resolution[1])
                reac <- reactiveValues(test = available_resolution[1])
                output$Markers <- renderDataTable(NULL)
                observeEvent(input$resolutionDE,{
                    reac$test <- input$resolutionDE
                    Idents(SeuratObjsubset) <- reac$test
                    Cluster_list <- levels(SeuratObjsubset)
                    updateSelectizeInput(session,"Cluster1",choices = Cluster_list, server =TRUE)
                    cluster_All <- c(Cluster_list,"All")
                    updateSelectizeInput(session,"Cluster2", choices = cluster_All ,server = TRUE)
                    output$DimplotMarkers <- renderPlot({
                        vizu_UMAP(SeuratObjsubset,var = input$resolutionDE, dlabel = input$addlabels_DE, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                    })
                    output$downloadUMAP_DE_page<- downloadHandler(
                      filename="UMAP_resolution_DE_page.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(SeuratObjsubset, var = input$resolution_DE,dlabel = input$addlabels_DE, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                        dev.off()
                      })
                    
                })
                output$Markers <- renderDataTable(
                  validate(
                    need(input$submit,"Provide the markers table for one cluster versus another or all. Possibility to adjust the minimum of average of log FC as a threshold for investigate in markers. you can choose to look for only the positive markers or to keep also the negative ones")
                  ))
                observeEvent(input$submit, {
                  Idents(SeuratObjsubset) <- reac$test
                  withProgress(message = "Differential expression", value = 0,{
                    output$Markers <- renderDataTable({
                      mList <- RunFindMarkers(SeuratObjsubset, cluster1 = input$Cluster1, cluster2 = input$Cluster2 ,posMarkers = input$onlyPos ,minimum.percent = input$percent, logthr = input$LogThreshold)
                    },server = FALSE ,extensions = c('Buttons'),
                    options = list(
                      autowidth = TRUE,
                      dom = 'frtipB',
                      buttons = list(list(extend ='csv',filename='TableMarkers'))
                    ),
                    filter = list(position = "top", clear = TRUE),)
                  })
                })
                
                ####################################
                ####### Data Mining one gene #######
                ####################################
                updateSelectizeInput(session, "ResolutionDataMining", choices=available_resolution, server=TRUE)
                toremove <- c(grep(names(SeuratObjsubset@meta.data),pattern ="*_snn_res.*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
                varAvail <- names(SeuratObjsubset@meta.data)
                varAvail <- varAvail[-which(varAvail %in% toremove)]
                updateSelectizeInput(session, "AnnotationDataMining", choices = varAvail, server = TRUE)
                gene <- SeuratObjsubset@assays$RNA@counts@Dimnames[[1]]
                updateSelectizeInput(session,"GeneList", choices=gene, server = TRUE)
                updateSelectizeInput(session, "AnnotAgainst", choices=varAvail, server = T)
                toListen <- reactive({list(input$ResolutionDataMining, input$AnnotationDataMining,input$format)})
                observeEvent(toListen(),{
                  if(input$format =="SVG"){
                      shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than tiff", "warning")
                  }
                
                  if(input$ResOrAnnot == "Resolution"){
                    Idents(SeuratObjsubset) <- input$ResolutionDataMining
                    ident_obj <- input$ResolutionDataMining
                  }else{
                    Idents(SeuratObjsubset) <- input$AnnotationDataMining
                    ident_obj <- input$AnnotationDataMining
                  }
                  observe({
                    for(i in c("downloadfeatureplot","VPdownload","HMdownload","DPdownload","downloadRawCount","downloadMatrix")){
                      shinyjs::toggleState(i,!is.null(input$GeneList))
                    }
                  })
                  ### Render FeaturePlot + saving them
                  lapply(seq(10),function(x) # allow to create the 10 graphs
                  {
                    output[[paste0('FeaturePlotMultGenes',x)]] = renderPlot({FP(SeuratObjsubset,input$GeneList[x], BoolColScales = input$colorScale, color = input$colorPickCol)})
                  })
                  output$downloadfeatureplot<- downloadHandler(
                    filename="featureplot.tiff",
                    content=function(file){
                      tiff(file, width = 900 , height = 600,res = 100)
                      print(FeaturePlot(SeuratObjsubset, features = input$GeneList))
                      dev.off()
                    }
                  )
                    
                    #Heatmap Generator 
                    output$Heatmap <- renderPlot({
                        Heatmap(SeuratObjsubset, genes = input$GeneList, var = ident_obj, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter, BoolColScales = input$colorScale, colorScale = input$colorPickCol)
                    })
                    if(input$format == "tiff"){
                        output$HMdownload <- downloadHandler(
                            filename = "heatmap.tiff",
                            content = function(file){
                                tiff(file, width = 900 , height = 600,res = 100)
                                print(Heatmap(SeuratObjsubset, genes = input$GeneList,var = ident_obj, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter,BoolColScales = input$colorScale, colorScale = input$colorPickCol))
                                dev.off()
                            }
                        )
                    }else{
                        output$HMdownload <- downloadHandler(
                            filename = "heatmap.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(Heatmap(SeuratObjsubset, genes = input$GeneList, var = ident_obj, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter, BoolColScales = input$colorScale, colorScale = input$colorPickCol))
                                dev.off()
                            }
                        )
                    }
                    
                    #Dotplot generator
                    output$Dotplot <- renderPlot({
                        Dotplot(SeuratObjsubset, genes = input$GeneList, BoolColScales = input$colorScale, colorScale = input$colorPickCol)
                    })
                    if(input$format == "tiff"){
                        output$DPdownload <- downloadHandler(
                            filename = "dotplot.tiff",
                            content = function(file){
                                tiff(file, width = 900 , height = 600,res = 100)
                                print(Dotplot(SeuratObjsubset, genes = input$GeneList, BoolColScales = input$colorScale, colorScale = input$colorPickCol))
                                dev.off()
                            }
                        )
                    }else{
                        output$DPdownload <- downloadHandler(
                            filename = "dotplot.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(Dotplot(SeuratObjsubset, genes = input$GeneList, BoolColScales = input$colorScale, colorScale = input$colorPickCol))
                                dev.off()
                            }
                        )
                    }
                    if(input$ResOrAnnot == "Resolution"){
                      dataTableCat <- reactive(aggregate(FetchData(SeuratObjsubset, vars=input$GeneList),data.frame(SeuratObjsubset[[input$ResolutionDataMining]]),mean))
                    }else{
                      dataTableCat <- reactive(aggregate(FetchData(SeuratObjsubset, vars=c(input$GeneList)),data.frame(SeuratObjsubset[[input$AnnotationDataMining]]),mean))
                    }
                    output$geneExp <- DT::renderDataTable({
                        validate(
                            need(input$GeneList != "","You have to choose at least one gene and here you will have for the mean of the expression for each cluster.")
                        )
                        DT::datatable(dataTableCat(), options = list(scrollX=TRUE), rownames = FALSE)
                    })
                    output$downloadRawCount<- downloadHandler(
                        filename = function(){
                            paste('data_',Sys.Date(),'.csv',sep = '')
                        },
                        content=function(file){
                            write.csv(dataTableCat(),file)
                        }
                    )
                    output$clusterOutput<- renderPlot({
                        vizu_UMAP(SeuratObjsubset,var= ident_obj,  dlabel = input$addlabels_ODM, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                    })
                    output$downloadUMAP_Cluster_resolution<- downloadHandler(
                      filename="UMAP_resolution_data_mining_page.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(SeuratObjsubset, var = ident_obj, dlabel = input$addlabels_ODM, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                        dev.off()
                    })
                    dataMatrixEntire <- reactive(FetchData(SeuratObjsubset, vars=input$GeneList))
                    output$downloadMatrix <- downloadHandler(
                      filename = function(){
                        paste('Matrix_',Sys.Date(),'.csv', sep ='')
                      },
                      content = function(file){
                        write.csv(dataMatrixEntire(), file)
                      }
                    )
                })

                ListenForViolin<- reactive({list(input$AnnotAgainst, input$SplitVln,input$ResolutionDataMining, input$AnnotationDataMining)})
                observeEvent(ListenForViolin(),{
                  Violins <- list()
                  if(input$ResOrAnnot == "Resolution"){
                    Idents(SeuratObjsubset) <- input$ResolutionDataMining
                    ident_obj_vln <- input$ResolutionDataMining
                    
                    
                  }else{
                    Idents(SeuratObjsubset) <- input$AnnotationDataMining
                    ident_obj_vln <- input$AnnotationDataMining
                    
                  }
                  updateSelectizeInput(session, "labelsToKeep1", choices = levels(SeuratObjsubset) , server =T)
                  
                  ### Render ViolinPlot + saving them
                  lapply(seq(10),function(x) # allow to create the 10 graphs
                  {
                    if(input$SplitVln == F){
                      output[[paste0('ViolinPlot',x)]] = renderPlot({SimpleViolinPlot(SeuratObjsubset,input$GeneList[x], var = ident_obj_vln, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)})
                    }else{
                      Idents(SeuratObjsubset) <- input$AnnotAgainst
                      updateSelectizeInput(session, "labelsToKeep2", choices = levels(SeuratObjsubset), server = T)
                      output[[paste0('ViolinPlot',x)]] = renderPlot({
                        gene_graph <- violinPlotSplited(SeuratObjsubset,input$GeneList[x],var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2=input$labelsToKeep2, color_list = list_Color_Label)
                        Violins[[x]] <- gene_graph
                        plot(gene_graph)
                      })
                    }
                  })
                  observeEvent(input$format, {
                    if(input$SplitVln == F && input$format== "tiff"){
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.tiff",
                        content = function(file){
                          tiff(file, width = 900 , height = 600,res = 100)
                          print(SimpleViolinPlot(SeuratObjsubset,input$GeneList, var = ident_obj_vln, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                          dev.off()
                        }
                      )
                    }else if (input$SplitVln == F && input$format == "SVG"){
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.svg",
                        content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(SimpleViolinPlot(SeuratObjsubset,input$GeneList, var = ident_obj_vln, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                          dev.off()
                        }
                      ) 
                    }else if(input$SplitVln == T && input$format =="tiff"){
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.tiff",
                        content = function(file){
                          tiff(file, width = 900 , height = 600,res = 100)
                          print(ViolinForDl(SeuratObjsubset, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2, color_list = list_Color_Label))
                          dev.off()
                        }
                      )
                    }else{
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.svg",
                        content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(ViolinForDl(SeuratObjsubset, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2, color_list = list_Color_Label))
                          dev.off()
                        }
                      )
                    }
                  })
                })
                ######################################################
                ############# Multiple genes data mining #############
                ######################################################
                updateSelectizeInput(session,"clusterwatch",choices =available_resolution,server =TRUE)
                toremove <- c(grep(names(SeuratObjsubset@meta.data),pattern ="*_snn_res.*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
                varAvail <- names(SeuratObjsubset@meta.data)
                varAvail <- varAvail[-which(varAvail %in% toremove)]
                updateSelectizeInput(session, "annotationwatch", choices = varAvail, server = TRUE)
                updateSelectizeInput(session,"GeneListPool", choices=gene, server = TRUE)
                coord <- Embeddings(SeuratObjsubset[["umap"]])[,1:2]
                listen_page_mining_mult_genes <- reactive({list(input$ResOrAnnotMult,input$annotationwatch,input$clusterwatch,input$SumorMeans,input$GeneListPool,input$format2)})
                #reactive graph from sum or mean button 
                observeEvent(listen_page_mining_mult_genes(),{
                    if(input$format2 =="SVG"){
                        shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than tiff", "warning")
                    }
                    if(input$ResOrAnnotMult == "Resolution"){
                      Idents(SeuratObjsubset) <- input$clusterwatch
                      ident_obj_min <- input$clusterwatch
                    }else{
                      Idents(SeuratObjsubset) <- input$annotationwatch
                      ident_obj_min <- input$annotationwatch
                    }
                  observe({
                    for(i in c("downloadMatrixMultList","downloadfeatureplotMultGenesPooled","VPdownloadMultGenesPooled","dotplotdownloadMultGenesPooled","downloadSumGeneExp")){
                      shinyjs::toggleState(i,!is.null(input$GeneListPool))
                    }
                  })
                  observe({
                    for(i in c("downloadMatrixMultFile","downloadfeatureplotMultGenesPooledFilesInput","downloadViolinPlotMultGenesPooledFilesInput","downloadDotPlotMultGenesPooledFilesInput","downloadmeanGeneExpPooledfromFile")){
                      shinyjs::toggleState(i,!is.null(input$fileGeneList))
                    }
                  })
                    #coord <- Embeddings(SeuratObjsubset[["umap"]])[,1:2]
                    output$dimplotSeurat4page <- renderPlot({
                        vizu_UMAP(SeuratObjsubset, var = ident_obj_min, dlabel = input$addlabels_MDM, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                    })
                    output$downloadUMAP_resolution_page_5<- downloadHandler(
                      filename="UMAP_resolution_MDM_page.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(SeuratObjsubset,var = ident_obj_min, dlabel = input$addlabels_MDM, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                        dev.off()
                      })
                    PoolExpTable <- reactive({
                      if(input$ResOrAnnotMult == "Resolution"){
                        out <- aggregate(FetchData(SeuratObjsubset, vars=c(input$GeneListPool)),data.frame(SeuratObjsubset[[input$clusterwatch]]),mean)
                      }else{
                        out <- aggregate(FetchData(SeuratObjsubset, vars=c(input$GeneListPool)),data.frame(SeuratObjsubset[[input$annotationwatch]]),mean)
                      }
                      out
                    })
                    output$sumGeneexp <- DT::renderDataTable({
                      validate(
                        need(input$GeneListPool !="", " ")
                      )
                      DT::datatable(PoolExpTable(), options = list(scrollX = TRUE), rownames = F)
                    })
                    output$downloadSumGeneExp<- downloadHandler(
                      filename = function(){
                        paste('data_',Sys.Date(),'.csv',sep = '')
                      },
                      content=function(file){
                        write.csv(PoolExpTable(),file)
                      }
                    )
                    
                    createEntireMatrix <- reactive({
                      if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                        if (input$SumorMeans == "Sum"){
                          if(input$ResOrAnnotMult == "Resolution"){
                            expressionMatrixEntire <- cbind(geneList = rowSums(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$clusterwatch]])
                          }else{
                            expressionMatrixEntire <- cbind(geneList = rowSums(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$annotationwatch]])
                          }
                        }else{
                          if(input$ResOrAnnotMult == "Resolution"){
                            expressionMatrixEntire <- cbind(geneList = rowMeans(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$clusterwatch]])
                          }else{
                            expressionMatrixEntire <- cbind(geneList = rowMeans(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$annotationwatch]])
                          }
                        }
                      }else{
                        sbt <- AddModuleScore(SeuratObjsubset, features =list(input$GeneListPool), name = "Gene_list")
                        if(input$ResOrAnnotMult == "Resolution"){
                          expressionMatrixEntire <- cbind(sbt[["Gene_list1"]],sbt[[input$clusterwatch]])
                        }else{
                          expressionMatrixEntire <- cbind(sbt[["Gene_list1"]],sbt[[input$annotationwatch]])
                        }
                      }
                      expressionMatrixEntire
                    })
                    output$downloadMatrixMultList <- downloadHandler(
                      filename = function(){
                        paste('Entire_matrix_',Sys.Date(),'.csv',sep = '')
                      },
                      content=function(file){
                        write.csv(createEntireMatrix(),file)
                      }
                    )
          
                  
                    output$FeaturePlotMultGenesPooled <- renderPlot({
                      validate(
                        need(input$GeneListPool != "","Choose at least one gene to display projection. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                      )
                      poolGene(SeuratObjsubset,gene = input$GeneListPool,typeOfNorm = input$SumorMeans, coord = coord, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_pool)
                    })
                    
                    output$ViolinPlotMultGenesPooled <- renderPlot({
                      validate(
                        need(input$GeneListPool != "", "Choose at least one gene to display violin plot. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                      )
                      if(input$ResOrAnnotMult == "Resolution"){
                        VlnPlotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                      }else{
                        VlnPlotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                      }
                      
                    })
                    
                    output$DotPlotMultGenePooled <- renderPlot({
                      validate(
                        need(input$GeneListPool != "", "Choose at least one gene to display dot plot. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                      )
                      if(input$ResOrAnnotMult == "Resolution"){
                        DotplotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled)
                      }else{
                        DotplotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled)
                      }
                    })
                    
                    if(input$format == "tiff"){
                        output$VPdownloadMultGenesPooled <- downloadHandler(
                            filename = "ViolinPlotMultGenesPooled.tiff",
                            content = function(file){
                                tiff(file, width = 900 , height = 600,res = 100)
                                print(if(input$ResOrAnnotMult == "Resolution"){
                                  VlnPlotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch,color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                                }else{
                                  VlnPlotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                                })
                                dev.off()
                            }
                        )
                    }else{
                        output$VPdownloadMultGenesPooled <- downloadHandler(
                            filename = "ViolinPlotMultGenesPooled.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(
                                  if(input$ResOrAnnotMult == "Resolution"){
                                    VlnPlotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                                  }else{
                                    VlnPlotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                                  }
                                )
                                dev.off()
                            }
                        )
                    }
                    
                    
                    
                      output$downloadfeatureplotMultGenesPooled<- downloadHandler(
                          filename = "featurePlotMultGenesPooled.tiff",
                          content = function(file){
                              tiff(file, width = 900 , height = 600,res = 100)
                              print(poolGene(SeuratObjsubset,gene = input$GeneListPool,typeOfNorm = input$SumorMeans, coord = coord, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_pool))
                              dev.off()
                          }
                      )
                    
                    
                    if(input$format2 == "tiff"){
                        output$dotplotdownloadMultGenesPooled<- downloadHandler(
                            filename = "DotPlotMultGenesPooled.tiff",
                            content = function(file){
                                tiff(file, width = 900 , height = 600,res = 100)
                                print(
                                  if(input$ResOrAnnotMult == "Resolution"){
                                    DotplotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled)
                                  }else{
                                    DotplotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled)
                                  }
                                )
                                dev.off()
                            }
                        )  
                    }else{
                        output$dotplotdownloadMultGenesPooled<- downloadHandler(
                            filename = "DotPlotMultGenesPooled.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(
                                  if(input$ResOrAnnotMult == "Resolution"){
                                    DotplotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled)
                                  }else{
                                    DotplotPooled(SeuratObjsubset,gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled)
                                  }
                                )
                                dev.off()
                            }
                        )
                    }   
                })
                ######################################### Gene list from file ################################
                #First of all we checked that every genes in the gene list are contained in the seurat object, if not we keep only those that are in the seurat object
                observeEvent(input$fileGeneList,{
                  CheckGenes <- reactive({
                    validate(
                      need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.")
                    )
                    list_gene_file <- read.csv(input$fileGeneList$datapath, header=TRUE)
                    colnames(list_gene_file) <- "cluster"
                    list_gene_file <- list_gene_file$cluster
                    before <- length(list_gene_file)
                    list_gene_file <- list_gene_file[list_gene_file %in% rownames(seuratObj)]
                    if(length(list_gene_file) != before){
                      shinyalert("Some genes cannot be found in this object, they have been removed", type = "warning")
                    }
                  })
                  CheckGenes()
                })
                #coord <- Embeddings(SeuratObjsubset[["umap"]])[,1:2]
                listen_page_mining_mult_genes_from_file <- reactive({list(input$ResOrAnnotMult,input$clusterwatch,input$annotationwatch,input$SumorMeans,input$fileGeneList, input$format2)})
                observeEvent(listen_page_mining_mult_genes_from_file(),{
                  # read file to and pass the good list for further analysis
                  Readfile <- function(){
                    validate(
                      need(input$fileGeneList, " ")
                    )
                    list_gene_file <- read.csv(input$fileGeneList$datapath, header=TRUE)
                    colnames(list_gene_file) <- "cluster"
                    list_gene_file <- list_gene_file$cluster
                    list_gene_file <- list_gene_file[list_gene_file %in% rownames(SeuratObjsubset)]
                    list_gene_file
                  }
                  exprPoolFromFile <- reactive({
                    list_gene_file <- Readfile()
                    if(input$ResOrAnnotMult == "Resolution"){
                      outFile <- aggregate(FetchData(SeuratObjsubset, vars=c(list_gene_file)), data.frame(SeuratObjsubset[[input$clusterwatch]]), mean)
                    }else{
                      outFile <- aggregate(FetchData(SeuratObjsubset, vars=c(list_gene_file)), data.frame(SeuratObjsubset[[input$annotationwatch]]), mean)
                    }
                    outFile
                  })
                  output$meanGeneExpPooledfromFile <- DT::renderDataTable({
                    DT::datatable(exprPoolFromFile(), options = list(scrollX =TRUE), rownames=F)
                  })
                  
                  output$downloadmeanGeneExpPooledfromFile<- downloadHandler(
                    filename = function(){
                      paste('data_',Sys.Date(),'.csv',sep = '')
                    },
                    content=function(file){
                      write.csv(exprPoolFromFile(),file)
                    }
                  )
                  
                  finaleExpressionMatrixFromFile <- reactive({
                    list_gene_file <- Readfile() 
                    if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                      if (input$SumorMeans == "Sum"){
                        if(input$ResOrAnnotMult == "Resolution"){
                          expressionMatrixEntireFile <- cbind(geneList = rowSums(FetchData(SeuratObjsubset, vars =c(list_gene_file))),SeuratObjsubset[[input$clusterwatch]])
                        }else{
                          expressionMatrixEntireFile <- cbind(geneList = rowSums(FetchData(SeuratObjsubset, vars =c(list_gene_file))),SeuratObjsubset[[input$annotationwatch]])
                        }
                      }else{
                        if(input$ResOrAnnotMult == "Resolution"){
                          expressionMatrixEntireFile <- cbind(geneList =rowMeans(FetchData(SeuratObjsubset, vars =c(list_gene_file))),SeuratObjsubset[[input$clusterwatch]])
                        }else{
                          expressionMatrixEntireFile <- cbind(geneList =rowMeans(FetchData(SeuratObjsubset, vars =c(list_gene_file))),SeuratObjsubset[[input$annotationwatch]])
                        }
                      }
                    }else{
                      sbt1 <- AddModuleScore(SeuratObjsubset, features =list_gene_file, name = "Gene_list")
                      if(input$ResOrAnnotMult == "Resolution"){
                        expressionMatrixEntireFile <- cbind(sbt1[["Gene_list1"]],sbt1[[input$clusterwatch]])
                      }else{
                        expressionMatrixEntireFile <- cbind(sbt1[["Gene_list1"]],sbt1[[input$annotationwatch]])
                      }
                    }
                    expressionMatrixEntireFile
                  })
                  
                  output$downloadMatrixMultFile <- downloadHandler(
                    filename = function(){
                      paste('Entire_matrix_',Sys.Date(),'.csv',sep = '')
                    },
                    content=function(file){
                      write.csv(finaleExpressionMatrixFromFile(),file)
                    }
                  )
                  
                  output$FeaturePlotMultGenesPooledFilesInput <- renderPlot({
                    validate(
                      need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.This plot will show the pool of each gene expression for each gene id contained in the file. \nExample of a file :\nHeader\nMlxipl\nSox9\nEtc...  ")
                    )
                    poolGene(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, coord = coord, file = input$fileGeneList, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_file)
                  })
                  
                  output$ViolinPlotMultGenesPooledFilesInput <- renderPlot({
                    validate(
                      need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.")
                    )
                    if(input$ResOrAnnotMult == "Resolution"){
                      VlnPlotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, cluster = input$clusterwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter) 
                    }else{
                      VlnPlotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, annotation = input$annotationwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter) 
                    }
                  })
                    
                  
                  output$DotPlotMultGenePooledFromFile <- renderPlot({
                    validate(
                      need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.")
                    )
                    if(input$ResOrAnnotMult == "Resolution"){
                      DotplotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, cluster = input$clusterwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled) 
                    }else{
                      DotplotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, annotation = input$annotationwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled) 
                    }
                  })
                  
                  if(input$format2 == "tiff"){
                      output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                          filename = "ViolinPlotMultGenesPooled.tiff",
                          content = function(file){
                              tiff(file, width = 900 , height = 600,res = 100)
                              print(if(input$ResOrAnnotMult == "Resolution"){
                                VlnPlotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, cluster = input$clusterwatch, BoolCol = input$BooleanColorsFilter ) 
                              }else{
                                VlnPlotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, annotation = input$annotationwatch, BoolCol = input$BooleanColorsFilter) 
                              })
                              dev.off()
                          }
                      )
                  }else{
                      output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                          filename = "ViolinPlotMultGenesPooled.svg",
                          content = function(file){
                              svg(file, width = 14 , height = 7)
                              print(
                                if(input$ResOrAnnotMult == "Resolution"){
                                  VlnPlotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, cluster = input$clusterwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter ) 
                                }else{
                                  VlnPlotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, annotation = input$annotationwatch, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter) 
                                }
                              )
                              dev.off()
                          }
                      )
                  }
                  
               
                    output$downloadfeatureplotMultGenesPooledFilesInput <- downloadHandler(
                        filename = "featurePlotFilesInput.tiff",
                        content = function(file){
                            tiff(file, width = 900 , height = 600,res = 100)
                            print(poolGene(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, coord = coord, file = input$fileGeneList, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_file))
                            dev.off()
                        }
                    ) 
                  
                  
                  if(input$format2 == "tiff"){
                      output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                          filename = "DotPlotMultGenesPooled.tiff",
                          content = function(file){
                              tiff(file, width = 900 , height = 600,res = 100)
                              print(
                                if(input$ResOrAnnotMult == "Resolution"){
                                  DotplotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, cluster = input$clusterwatch ) 
                                }else{
                                  DotplotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, annotation = input$annotationwatch) 
                                }
                              )
                              dev.off()
                          }
                      ) 
                  }else{
                      output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                          filename = "DotPlotMultGenesPooled.svg",
                          content = function(file){
                              svg(file, width = 14 , height = 7)
                              print(
                                if(input$ResOrAnnotMult == "Resolution"){
                                  DotplotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, cluster = input$clusterwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled) 
                                }else{
                                  DotplotPooled(obj = SeuratObjsubset, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, file = input$fileGeneList, annotation = input$annotationwatch, BoolColScales = input$colorScalePooled, colorScale = input$colorPickColPooled) 
                                }
                              )
                              dev.off()
                          }
                      )
                  }
              })
                
                
                ###########################################
                ########  Information extraction #########
                ###########################################
                updateSelectizeInput(session,"clusterNumber",choices = available_resolution,server =TRUE)
                toremove <- c(grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
                varAvail <- names(SeuratObjsubset@meta.data)
                varAvail <- varAvail[-which(varAvail %in% toremove)]
                updateSelectizeInput(session,"chooseVar2Plot", choices=varAvail, server=TRUE)
                
                cluster <- reactiveValues(corresponding_res = "*_snn_res.0") # Need a reactive value to pass through the observe event 
                observeEvent(input$clusterNumber,{
                    Idents(SeuratObjsubset) <- input$clusterNumber
                    DMclusters<- c(levels(SeuratObjsubset), "All")
                    updateSelectizeInput(session,"clusterInfo",choices = DMclusters, server =TRUE)
                    output$Dimplot_cluster <- renderPlot({
                        vizu_UMAP(SeuratObjsubset,var = input$clusterNumber, dlabel = input$addlabels_info, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                    })
                    output$downloadUMAP_info<- downloadHandler(
                      filename = "UMAP_info.tiff",
                      content = function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(SeuratObjsubset,var = input$clusterNumber, dlabel = input$addlabels_info, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                        dev.off()
                      }
                    )
                    
                    toListen <- reactive({list(input$clusterInfo,input$chooseVar2Plot, input$stackedBP, input$FreqOrVal)})
                    observeEvent(toListen(), {  
                      output$ggplot_information <- renderPlot({
                        p <- Plot2Render(obj = SeuratObjsubset,cluster = input$clusterInfo, StackedBarPlot = input$stackedBP,VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, freqOrValues = input$FreqOrVal, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                        p
                      })
                      output$legend_information <- renderPlot({
                        p <- Plot2Render(obj = SeuratObjsubset,cluster = input$clusterInfo, StackedBarPlot = input$stackedBP,VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, freqOrValues = input$FreqOrVal, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter )
                        grid.newpage()
                        legend <- cowplot::get_legend(p)
                        grid.draw(legend)
                      })
                      output$cluster_percent <- renderPlot({
                          p <- percentRender(obj = SeuratObjsubset, VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                          p
                      })
                      output$downloadInformationPlot2 <- downloadHandler(
                        filename = "Percent_Plot.tiff",
                        content = function(file){
                          tiff(file, width = 900 , height = 600,res = 100)
                          print(percentRender(obj = SeuratObjsubset, VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                          dev.off()
                        }
                      )
                      
                      output$get_info <- DT::renderDataTable({
                          Table2Render(SeuratObjsubset,VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, cluster = input$clusterInfo)
                      })
                      output$downloadInformationPlot <- downloadHandler(
                          filename = "Information_Plot.tiff",
                          content = function(file){
                              tiff(file, width = 900 , height = 600,res = 100)
                              print(Plot2Render(obj = SeuratObjsubset,cluster = input$clusterInfo, StackedBarPlot = input$stackedBP,VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, freqOrValues = input$FreqOrVal, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                              dev.off()
                          }
                      )
                      output$downloadMatrixInformation <- downloadHandler(
                          filename = "Information_Matrix.csv",
                          content = function(file){
                              write.csv(Table2Render(SeuratObjsubset,VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, cluster = input$clusterInfo),file)
                          }
                      )
                    })
                })
                #######################################
                ########### Add Annotation ############
                #######################################
                removeInfo <- c(grep(names(SeuratObjsubset@meta.data),pattern ="*_snn_res.*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
                list.metadata <- colnames(SeuratObjsubset@meta.data)
                list.metadata <- list.metadata[-which(list.metadata %in% removeInfo)]
                updateSelectizeInput(session,"Annot2Complete", choices =list.metadata)
                updateSelectizeInput(session,"res4Annot",choices = available_resolution,server =TRUE)
                resChoice <- NULL
                observeEvent(input$res4Annot,{
                    resChoice <<- input$res4Annot
                    Idents(SeuratObjsubset) <- input$res4Annot
                    clusterSbt <- levels(SeuratObjsubset)
                    updateSelectizeInput(session,"cluster4Annot",choices = clusterSbt, server = TRUE)
                    output$dimplot_res_annot <- renderPlot({
                        vizu_UMAP(SeuratObjsubset,var = input$res4Annot, dlabel = input$labelsAnnot,color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                    })
                    output$downloadUMAP_resolution_annot_page<- downloadHandler(
                      filename="UMAP_resolution_annot_page.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(SeuratObjsubset,var = input$res4Annot, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                        dev.off()
                      })
                })
                output$condPlot <- renderPlot({
                    vizu_UMAP(SeuratObjsubset, var = input$Annot2Complete, dlabel = input$labelsconditional, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                })
                output$downloadUMAP_Conditional<- downloadHandler(
                  filename="UMAP_conditional_page.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(SeuratObjsubset,var = input$Annot2Complete, dlabel = input$labelsconditional, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                    dev.off()
                  })
                observeEvent(input$submitAnnot,{
                    if(input$choicesCreateOrComplete == "Create"){
                        if(input$AnnotName == "" | input$clusterName == "" | length(input$cluster4Annot) == 0){
                            shinyalert("Oops", "The name of the column and the name of the annotation you want to give has to be filled", type = "error")
                        }else{
                          for(i in c("downloadUMAP_post_annotation","labelsPostAnnot")){
                            shinyjs::show(i)
                          }
                          if(is.element(input$clusterName, names(SeuratObjsubset@meta.data))){
                            shinyalert("Warning","The name that you give to your annotation is already existing in your dataset. Do you want to continue ?","warning", showConfirmButton = TRUE, showCancelButton = TRUE,callbackR = function(x){
                              if(x == TRUE){
                                annot <- Annotation_existent(SeuratObjsubset, nameVar = input$clusterName, annotation = input$AnnotName, resolution = resChoice,clusterToAnnotate = input$cluster4Annot)
                                SeuratObjsubset[[input$clusterName]] <<- annot
                                shinyalert("Done","You have correctly annotated your cluster", type = "success")
                                output$postAnnotation <- renderPlot({
                                  vizu_UMAP(SeuratObjsubset, var = input$clusterName, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                                }) 
                              }
                            } 
                            )
                          }else{
                            annot <- Annotation_unexistent(SeuratObjsubset, annotation = input$AnnotName, resolution = resChoice,clusterToAnnotate = input$cluster4Annot)
                            SeuratObjsubset[[input$clusterName]] <<- annot
                            shinyalert("Done","You have correctly annotated your cluster", type = "success")
                            output$postAnnotation <- renderPlot({
                              vizu_UMAP(SeuratObjsubset, var = input$clusterName, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                            }) 
                          }
                            output$downloadUMAP_post_annotation<- downloadHandler(
                              filename="UMAP_post_annot.tiff",
                              content=function(file){
                                tiff(file, width = 900 , height = 600,res = 100)
                                print(vizu_UMAP(SeuratObjsubset,var = input$clusterName, var = input$labelsPostAnnot,color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                                dev.off()
                              })
                            toremove <- c(grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T), grep(names(SeuratObjsubset@meta.data),pattern ="snn_res.*", value =T),"percent.mt","S.Score", "G2M.Score")
                            choices <- names(SeuratObjsubset@meta.data)
                            choicesForAnnot <- choices[-which(choices %in% toremove)]
                            updateSelectizeInput(session,"chooseVar2Plot", choices=choices)
                            updateSelectizeInput(session,"whichAnnot", choices=choicesForAnnot)
                            updateSelectizeInput(session,"Annot2Complete",choices = choicesForAnnot)
                            updateSelectizeInput(session,"AnnotationDataMining",choices = choicesForAnnot)
                            updateSelectizeInput(session,"annotationwatch",choices = choicesForAnnot)
                            updateSelectizeInput(session,"AnnotAgainst", choices = choicesForAnnot)
                            
                            
                        }
                        
                    }else{
                        if(input$AnnotName == "" | length(input$cluster4Annot) == 0){
                            shinyalert("Oops", "The name of the column and the name of the annotation you want to give has to be filled", type = "error")
                        }else{
                            f.clusters <- SeuratObjsubset[[input$Annot2Complete]]
                            f.clusters[SeuratObjsubset[[resChoice]] == input$cluster4Annot[1]] <- input$AnnotName
                            if(length(input$cluster4Annot) > 1){
                                for(i in 2:length(input$cluster4Annot)){
                                    f.clusters[SeuratObjsubset[[resChoice]] == input$cluster4Annot[i]] <- input$AnnotName 
                                }
                            }
                        }
                        output$postAnnotation <- renderPlot({
                            vizu_UMAP(SeuratObjsubset, var = input$Annot2Complete, dlabel = ,input$labelsPostAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                        })
                        output$downloadUMAP_post_annotation<- downloadHandler(
                          filename="UMAP_post_annot.tiff",
                          content=function(file){
                            tiff(file, width = 900 , height = 600,res = 100)
                            print(vizu_UMAP(SeuratObjsubset,var = input$Annot2Complete, input$labelsPostAnnot,color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                            dev.off()
                          })
                        SeuratObjsubset[[input$Annot2Complete]] <<- f.clusters
                        shinyalert("Done","You have correctly annotated your cluster", type = "success")    
                    }
                    
                })
                output$downloadRDSwithNewAnnot <- downloadHandler(
                    filename = "SeuratObj.rds",
                    content = function(file){
                        saveRDS(SeuratObjsubset,file)
                    }
                )
                
                ######################################
                ######## Automatic annotation ########
                ######################################
                shinyjs::disable("downloadProbabilities")
                gene <- SeuratObjsubset@assays$RNA@counts@Dimnames[[1]]
                listbutton <- reactive({list(input$NbCellType)})
                observe({
                  shinyjs::toggleState("AutomaticAnnotation",input$slotNameAnnotScina != "")
                })
                observeEvent(listbutton(),{
                  output$testPanel <- renderUI({
                    tagList(
                      lapply(seq(1,input$NbCellType), function(i){
                        tagList(
                          textInput(paste("cellType", i, sep="_"), "Give a name to your cell type"),
                          selectizeInput(paste("markerCellType",i,sep="_"), "Give some markers (between 2 and 20 is advised)",choices=character(0), multiple = T)
                        )
                      }) 
                    ) 
                  })
                  lapply(seq(1,input$NbCellType), function(i){
                    updateSelectizeInput(session,paste("markerCellType",i,sep="_"),choices=gene, server = TRUE)
                  })
                })
                output$FileExplanation <- renderText(
                  HTML(paste("If you want to pass a file here is an example of file :","Cell_type 1;Cell_type 2","Marker1_CT1;Marker1_CT2","Marker2_CT1;Marker2_CT2","Marker3_CT1;Marker3_CT2","Marker4_CT1;","Marker5_CT1;","<b>Important notice</b> :", "Do not use space to label your cell type","The list can have different numbers of markers", sep ="<br/>"))
                )
                observeEvent(input$AutomaticAnnotation,{
                  
                  shinyjs::disable(id = input$AutomaticAnnotation)
                  if(input$typeOfList == "List"){
                    
                    cell_type <- reactive({
                      lapply(seq(1,input$NbCellType), function(i) {
                        input[[paste("cellType", i, sep="_")]]
                      })
                    })
                    marker <- reactive({
                      lapply(seq(1,input$NbCellType), function(i) {
                        input[[paste("markerCellType", i, sep="_")]]
                      })
                    })
                    scina_list <- list()
                    marker_list <- marker()
                    celltype_list <- cell_type()
                    for(i in 1:length(marker_list)){
                      tmp <- unlist(celltype_list[i])
                      scina_list[tmp] <- marker_list[i]
                    }
                  }else{
                    scina_list <- read.csv(input$InputMarkerFile$datapath,sep = ";")
                  }
                  scina_results<- useSCINA(object = SeuratObjsubset, list_annotation = scina_list,allowUnknown = TRUE, allowOverlap = TRUE)
                  SeuratObjsubset <<- Fill_proba(object = SeuratObjsubset, slotName = input$slotNameAnnotScina, scina = scina_results)
                  shinyjs::enable(input$AutomaticAnnotation)
                  proba_avail <- c(grep(names(SeuratObjsubset@meta.data),pattern ="*_probabilities", value =T))
                  updateSelectizeInput(session,"Proba_selection", choices =proba_avail)
                  shinyjs::enable("downloadProbabilities")
                  shinyjs::show("Scina_annot_DL")
                  shinyjs::show("Probabilities_Dl")
                  shinyjs::show("Proba_selection")
                  output$Results_annotation <- renderPlot(
                    Dimplot_with_ggplot(object = SeuratObjsubset, group = input$slotNameAnnotScina)
                  )
                  output$Probabilities <- renderPlot(
                    FeaturePlot_proba(SeuratObjsubset, proba_col = input$Proba_selection, slotNameAnnot = input$slotNameAnnotScina)
                  )
                  output$Scina_annot_DL <- downloadHandler(
                    filename = "Scina_annotation_results.tiff",
                    content=function(file){
                      tiff(file, width = 900 , height = 600,res = 100)
                      print(Dimplot_with_ggplot(object = SeuratObjsubset, group = input$slotNameAnnotScina))
                      dev.off()
                    }
                  )
                  
                  output$Probabilities_Dl <- downloadHandler(
                    filename = "Scina_annotation_probabilities.tiff",
                    content=function(file){
                      tiff(file, width = 900 , height = 600,res = 100)
                      print(FeaturePlot_proba(SeuratObjsubset, proba_col = input$Proba_selection, slotNameAnnot = input$slotNameAnnotScina))
                      dev.off()
                    }
                  )
                  output$downloadProbabilities <- downloadHandler(
                    filename = "Scina_annotation_probabilities.csv",
                    content=function(file){
                      write.csv(t(scina_results$probabilities),file)
                    }
                  )
                  gc()
                  toremove <- c(grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="snn_res.*", value =T),"percent.mt","S.Score", "G2M.Score" ,grep(names(SeuratObjsubset@meta.data),pattern ="*_probabilities", value =T))
                  choices <- names(SeuratObjsubset@meta.data)
                  choices_annot <- choices[-which(choices %in% toremove )]
                  updateSelectizeInput(session,"InfoToPlotFilter", choices=choices_annot)
                  updateSelectizeInput(session,"whichAnnot", choices=choices_annot)
                  updateSelectizeInput(session,"Annot2Complete", choices =choices_annot)
                  updateSelectizeInput(session, "AnnotationDataMining", choices = choices_annot)
                  updateSelectizeInput(session,"annotationwatch",choices = choices_annot)
                  updateSelectizeInput(session,"AnnotAgainst", choices = choices_annot)
                  updateSelectizeInput(session,"InfoToPlot", choices = choices_annot)  
                })
                ######################################
                ########### Subclustering ############
                ######################################
                removeInfo <- c(grep(names(SeuratObjsubset@meta.data),pattern ="*_snn_res.*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nCount_*", value =T),grep(names(SeuratObjsubset@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score","seurat_clusters")
                annotation_sub <- colnames(SeuratObjsubset@meta.data)
                annotation_sub <- annotation_sub[-which(annotation_sub %in% removeInfo)]
                updateSelectizeInput(session,"clusterResPage7",choices = available_resolution,server =TRUE)
                updateSelectizeInput(session,"whichAnnot", choices = annotation_sub, server = TRUE)
                listenbutton <- reactive({list(input$clusterResPage7, input$whichAnnot,input$subcluster)})
                shinyjs::hide("NameObject")
                shinyjs::hide("downloadSubRds")
                shinyjs::hide("downloadLog")
                subsetSeuratObj <- NULL
                observeEvent(listenbutton(),{
                    test <- NULL
                    if(input$subcluster == "Cluster"){
                        Idents(SeuratObjsubset) <- input$clusterResPage7
                        clusterSbt <- levels(SeuratObjsubset)
                        updateSelectizeInput(session,"subclusterize",choices = clusterSbt, server = TRUE )
                        output$Dimplot_subclustering <- renderPlot({
                            vizu_UMAP(SeuratObjsubset,var = input$clusterResPage7, dlabel = input$labelsSubset, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                        })
                        output$downloadUMAP_sbt<- downloadHandler(
                          filename="UMAP_sbt.tiff",
                          content=function(file){
                            tiff(file, width = 900 , height = 600,res = 100)
                            print(vizu_UMAP(SeuratObjsubset,var = input$clusterResPage7, input$labelsSubset,color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                            dev.off()
                          })
                        observeEvent(input$subclusterize,{
                            test <- NULL
                            for (i in c("downloadUMAPcellKept","cellsHighlightSize","cellsSize")){
                              shinyjs::show(i)
                            }
                            if(length(input$subclusterize)==length(clusterSbt)){
                                test <- colnames(SeuratObjsubset)
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.tiff",
                                  content=function(file){
                                    tiff(file, width = 900 , height = 600,res = 100)
                                    print(MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                    dev.off()
                                  })
                            }else{
                              test <- colnames(SeuratObjsubset)[which(SeuratObjsubset[[]][input$clusterResPage7] == input$subclusterize[1])]
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test ,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.tiff",
                                  content=function(file){
                                    tiff(file, width = 900 , height = 600,res = 100)
                                    print(MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                    dev.off()
                                  })
                            }
                            if(length(input$subclusterize) > 1){
                                if(length(input$subclusterize) == length(clusterSbt)){
                                    test <- colnames(SeuratObjsubset)
                                    output$Dimplot_kept <- renderPlot({
                                        MakeUMAP_Red(SeuratObjsubset,highliht_cells = test ,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                    })
                                    output$downloadUMAPcellKept<- downloadHandler(
                                      filename="UMAP_cell_kept.tiff",
                                      content=function(file){
                                        tiff(file, width = 900 , height = 600,res = 100)
                                        print(MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                        dev.off()
                                      })
                                }else{
                                    for(i in 2:length(input$subclusterize)){
                                        test <- c(test, colnames(SeuratObjsubset)[which(SeuratObjsubset[[]][input$clusterResPage7] == input$subclusterize[i])])
                                        output$Dimplot_kept <- renderPlot({
                                            MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                        })
                                        output$downloadUMAPcellKept<- downloadHandler(
                                          filename="UMAP_cell_kept.tiff",
                                          content=function(file){
                                            tiff(file, width = 900 , height = 600,res = 100)
                                            print(MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                            dev.off()
                                          })
                                    }
                                }  
                            }
                            subsetSeuratObj <<- subset(SeuratObjsubset, cells = test)
                            
                        })
                    }else{
                        Idents(SeuratObjsubset) <- input$whichAnnot
                        annotSbt <- levels(SeuratObjsubset)
                        updateSelectInput(session, "subannot", choices = annotSbt)
                        output$Dimplot_subclustering <- renderPlot({
                            vizu_UMAP(SeuratObjsubset,var = input$whichAnnot, dlabel = input$labelsSubset, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter)
                        })
                        output$downloadUMAP_sbt<- downloadHandler(
                          filename="UMAP_sbt.tiff",
                          content=function(file){
                            tiff(file, width = 900 , height = 600,res = 100)
                            print(vizu_UMAP(SeuratObjsubset,var = input$whichAnnot, input$labelsSubset, color_list = list_Color_Label, BoolCol = input$BooleanColorsFilter))
                            dev.off()
                          })
                        observeEvent(input$subannot,{
                          for (i in c("downloadUMAPcellKept","cellsHighlightSize","cellsSize")){
                            shinyjs::show(i)
                          }
                            test <- NULL
                            if(length(input$subannot) == length(annotSbt)){
                                test <- colnames(SeuratObjsubset)
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.tiff",
                                  content=function(file){
                                    tiff(file, width = 900 , height = 600,res = 100)
                                    print(MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                    dev.off()
                                  })
                            }else{
                              test <- colnames(SeuratObjsubset)[which(SeuratObjsubset[[]][input$whichAnnot] == input$subannot[1])]
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.tiff",
                                  content=function(file){
                                    tiff(file, width = 900 , height = 600,res = 100)
                                    print(MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                    dev.off()
                                  })
                            }
                            if(length(input$subannot) > 1){
                                if(length(input$subannot) == length(annotSbt)){
                                    test <- colnames(SeuratObjsubset)
                                    output$Dimplot_kept <- renderPlot({
                                        MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                    })
                                    output$downloadUMAPcellKept<- downloadHandler(
                                      filename="UMAP_cell_kept.tiff",
                                      content=function(file){
                                        tiff(file, width = 900 , height = 600,res = 100)
                                        print(MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                        dev.off()
                                      })
                                }else{
                                    for(i in 2:length(input$subannot)){
                                        test <- c(test, colnames(SeuratObjsubset)[which(SeuratObjsubset[[]][input$whichAnnot] == input$subannot[i])])
                                        output$Dimplot_kept <- renderPlot({
                                            MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                        })
                                        output$downloadUMAPcellKept<- downloadHandler(
                                          filename="UMAP_cell_kept.tiff",
                                          content=function(file){
                                            tiff(file, width = 900 , height = 600,res = 100)
                                            print(MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                            dev.off()
                                          })
                                    }
                                }  
                            }
                            subsetSeuratObj <<- subset(SeuratObjsubset, cells = test) # the double <<- is used to fill a value wich is in observeEvent but that we want to keep for after
                            
                        })
                    }
                    
                })
                observe({
                    shinyjs::toggleState("dosubcluster",input$subcluster == "Annotation" && input$subannot != "" || input$subcluster == "Cluster" && input$subclusterize !="")
                })
                observeEvent(input$dosubcluster,{ #we do not normalize data once the subset is done
                  shinyjs::disable("dosubcluster")
                  withProgress(message = "Renormalizing data", value = 0,{
                    incProgress(1/4,message ="Running PCA")
                    subsetSeuratObj <<- RunPCA(subsetSeuratObj)
                    incProgress(1/4,message ="Finish PCA")
                    subsetSeuratObj <<- FindNeighbors(subsetSeuratObj)
                    incProgress(1/4,message ="Running UMAP")
                    subsetSeuratObj <<- RunUMAP(subsetSeuratObj, dims = 1:30)
                    incProgress(1/4,message ="Running UMAP")
                  })
                  withProgress(message ="Find cluster", value=0,{
                    for(i in seq(0,1,0.1)){
                      incProgress(i*10/10,message =paste0("res ",i*10,"/10"))
                      subsetSeuratObj <<- FindClusters(subsetSeuratObj, res =i)
                    }
                    
                  })
                  for(i in c("downloadUMAPaftersbt", "labelsAfterSubset")){
                    shinyjs::show(i)
                  }
                  
                  output$dimplot_aftersbt <- renderPlot({
                    vizu_UMAP(subsetSeuratObj, var = "seurat_clusters" ,dlabel = input$labelsAfterSubset)
                  })
                  output$downloadUMAPaftersbt<- downloadHandler(
                    filename="UMAP_After_sbt.tiff",
                    content=function(file){
                      tiff(file, width = 900 , height = 600,res = 100)
                      print(vizu_UMAP(subsetSeuratObj, var = "seurat_clusters", dlabel = input$labelsAfterSubset))
                      dev.off()
                    })
                  shinyjs::enable("dosubcluster")
                  shinyjs::show("NameObject")
                  shinyjs::show("downloadSubRds")
                  shinyjs::show("downloadLog")
                })
                observe({
                  shinyjs::toggleState("downloadSubRds",input$NameObject != "")
                })
                output$downloadSubRds <- downloadHandler(
                  filename = function(){
                    paste0(input$NameObject,".rds")
                  },
                  content = function(file){
                    shinyalert("warning", "Seurat object can be heavy, and it can take time, please be patient, don't click multiple times", "warning")
                    saveRDS(subsetSeuratObj,file)
                  }
                  
                )
                observe({
                  shinyjs::toggleState("downloadLog",input$NameObject)
                })
                output$downloadLog <- downloadHandler(
                  filename = function(){
                    paste0(input$NameObject,"_log.html")
                  },
                  content = function(file){
                    tempReport <- file.path(tempdir(),"logFile.Rmd")
                    file.copy("logFile.Rmd",tempReport ,overwrite = TRUE)
                    command <- list()
                    line <- paste0("SeuratObjsubset <- subset(seuratObj, subset = nFeature_RNA > ",input$minGenePerCells," & nFeature_RNA < ",input$maxGenePerCells," & nCount_RNA < ",input$maxCountPerCells," & percent.mt < ",input$percentMito,")")
                    command <- append(command,line)
                    if(input$normalization == "SCTransform"){
                      line <- paste0("SeuratObjsubset <- SCTransform(SeuratObjsubset)")
                      command <- append(command,line)
                    }else{
                      line <- paste0("SeuratObjsubset <- NormalizeData(SeuratObjsubset)")
                      command <- append(command,line)
                      line <- paste0("SeuratObjsubset <- FindVariableFeatures(SeuratObjsubset)")
                      command <- append(command,line)
                      line <- paste0("SeuratObjsubset <- ScaleData(SeuratObjsubset)")
                      command <- append(command,line)
                    }
                    command <- append(command,"SeuratObjsubset <- RunPCA(SeuratObjsubset)")
                    command <- append(command,"SeuratObjsubset <- FindNeighbors(SeuratObjsubset)")
                    command <- append(command,"SeuratObjsubset <- RunUMAP(SeuratObjsubset, dims = 1:30)")
                    for(i in seq(input$Resolution[1],input$Resolution[2], input$ResolutionStep)){
                      command <- append(command,paste0("SeuratObjsubset <- FindClusters(SeuratObjsubset, res =", i,")"))
                    }
                    if(input$subcluster == "Cluster"){
                      command <- append(command,paste0("subsetSeuratObj <- subset(SeuratObjsubset, subset =",input$clusterResPage7," == ",list(input$subclusterize),")"))
                    }else{
                      command <- append(command,paste0("subsetSeuratObj <- subset(SeuratObjsubset, subset = ",input$whichAnnot," == ",list(input$subannot),")"))
                    }
                    command <- append(command,("subsetSeuratObj <- RunPCA(subsetSeuratObj)"))
                    command <- append(command,("subsetSeuratObj <- FindNeighbors(subsetSeuratObj)"))
                    command <- append(command,("subsetSeuratObj <- RunUMAP(subsetSeuratObj, dims = 1:30)"))
                    for(i in seq(0,1,0.1)){
                      command<- append(command,paste0("subsetSeuratObj <- findClusters(subsetSeuratObj, res = ",i,")"))
                    }
                    params <- list(use = command)
                    rmarkdown::render(tempReport,output_file = file, params = params, envir = new.env(parent = globalenv()))
                  }
                )
            })
            #########################################################################################################################################################################
            ################################################################## If you have a RDS file ###############################################################################
            #########################################################################################################################################################################            
        }else{
            ####################################
            ############ QC metrics ############
            ####################################
            available_metadata <- colnames(seuratObj@meta.data)
            available_metadata <- available_metadata[-which(available_metadata %in% c("nCount_RNA","nFeature_RNA","percent.mt","S.Score", "G2M.Score","nCount_SCT", "nFeature_SCT"))]
            updateSelectizeInput(session, "InfoToPlot", choices= available_metadata, server=TRUE)
            observe({
              shinyjs::toggleState("downloadColorFile", input$BooleanColors != FALSE)
            })
            ## Violin plot on data
            output$featureQC <- renderPlot({
                makeVlnGreatAgain(seuratObj,  grouping = "orig.ident",var = c("nFeature_RNA","nCount_RNA","percent.mt"), col = 3) 
            })
            output$downloadVln<- downloadHandler(
              filename="Violin_QC.tiff",
              content=function(file){
                tiff(file, width = 900 , height = 600,res = 100)
                print(makeVlnGreatAgain(seuratObj, grouping = "orig.ident", var=c("nFeature_RNA","nCount_RNA","percent.mt"), col = 3))
                dev.off()
              })
            
            ## QC resolution
            output$plotClusterQC <- renderPlot({
                vizu_UMAP(seuratObj,var= "seurat_clusters",dlabel = input$addlabels_res)
            })
            output$downloadUMAP_resolution_qc<- downloadHandler(
              filename="UMAP_Cluster_QC.tiff",
              content=function(file){
                tiff(file, width = 900 , height = 600,res = 100)
                print(vizu_UMAP(seuratObj,var=NULL, dlabel = input$addlabels_res))
                dev.off()
              })
            
            ## cell cycle
            output$cellcycleQC <- renderPlot({
                validate(
                    need(try(seuratObj[["Phase"]]!=""),message = "No phase available in this seurat object")
                )
                vizu_UMAP(seuratObj,var = "Phase", dlabel = input$showlabelccQC)
            })
            output$downloadUMAP_CC<- downloadHandler(
              filename="UMAP_Cell_cycle.tiff",
              content=function(file){
                tiff(file, width = 900 , height = 600,res = 100)
                print(vizu_UMAP(seuratObj,var="Phase", dlabel = input$showlabelccQC))
                dev.off()
              })
            
            ## Hist nb cell by dataset
            output$nbcells_by_datasetsQC <- renderPlot({
                nbCellsbydt(seuratObj)
            })
            output$downloadHistNbCells<- downloadHandler(
              filename="Hist_nb_cells_dt.tiff",
              content=function(file){
                tiff(file, width = 800 , height = 600,res = 100)
                print(nbCellsbydt(seuratObj))
                dev.off()
              })
            
            ## Nb gene by cell
            output$geneExpByCellQC <- renderPlot({
                UmapNbGenesCell(seuratObj, sizePoint = input$pointSizeNbCellQC)
            })
            output$downloadUMAP_gene_by_cell<- downloadHandler(
              filename="UMAP_Nb_gene_by_cell.tiff",
              content=function(file){
                tiff(file, width = 800 , height = 600,res = 100)
                print(UmapNbGenesCell(seuratObj, sizePoint = input$pointSizeNbCellQC))
                dev.off()
              })
            
            list_Color_Label <- list()
            ### Color changer recursively  
            listenColors <- reactive({list(input$BooleanColors, input$InfoToPlot, input$submitAnnot)})
            observeEvent(listenColors(),{
              Idents(seuratObj) <- input$InfoToPlot
              orderLevel <- sort(levels(seuratObj))
              if(input$BooleanColors == TRUE){
                shinyjs::show("myPanel")
                cols <- reactive({
                  lapply(seq_along(levels(seuratObj)), function(i) {
                    colourpicker::colourInput(paste("col", i, sep="_"), paste("Choose colour for",orderLevel[i]), "black")
                  })
                })
                output$myPanel <- renderUI({cols()})
                # Put all the input in a vector
                colors <- reactive({
                  lapply(seq_along(levels(seuratObj)), function(i) {
                    input[[paste("col", i, sep="_")]]
                  })
                })
                render_colors <- function(){
                  if (is.null(input$col_1)) {
                    cols <- rep("#000000", length(levels(seuratObj)))
                  } else {
                    cols <- unlist(colors())
                    labels_object <- levels(seuratObj)
                    list_Color_Label[[input$InfoToPlot]]$color <<- cols
                    list_Color_Label[[input$InfoToPlot]]$label <<- labels_object
                  }
                  return(cols)
                }
                
                output$plotChoice <- renderPlot({
                  cols <- render_colors()
                  vizu_UMAP(seuratObj,var =input$InfoToPlot,  dlabel = input$addlabels_choice, color = cols, sizePoint = input$sizePoint_choice)
                })
                output$downloadUMAP_choice<- downloadHandler(
                  filename="UMAP_choice.tiff",
                  content=function(file){
                    tiff(file, width = 600 , height = 600,res = 100)
                    cols <- render_colors()
                    print(vizu_UMAP(seuratObj,var=input$InfoToPlot,  dlabel = input$addlabels_choice, color = cols, sizePoint = input$sizePoint_choice))
                    dev.off()
                })
              }else{
                shinyjs::hide("myPanel")
                output$plotChoice <- renderPlot({
                  vizu_UMAP(seuratObj,var =input$InfoToPlot,  dlabel = input$addlabels_choice, sizePoint = input$sizePoint_choice)
                })
                output$downloadUMAP_choice<- downloadHandler(
                  filename="UMAP_choice.tiff",
                  content=function(file){
                    tiff(file, width = 600 , height = 500, res = 100 )
                    print(vizu_UMAP(seuratObj,var=input$InfoToPlot,  dlabel = input$addlabels_choice, sizePoint = input$sizePoint_choice))
                    dev.off()
                })
              }
            })
            
            output$downloadColorFile<- downloadHandler(
              filename="Color_file.csv",
              content=function(file){
                write_color_file(color_list = list_Color_Label, name_file = file)
              })
            
            
            #############################
            ##### Cluster Tree page #####
            #############################
            available_resolution <- grep(colnames(seuratObj@meta.data), pattern = "*_snn_res.*", value = TRUE)
            updateSelectizeInput(session, "resolution_TreePage", choices=available_resolution, server=TRUE)
            output$ClusterTree <- renderPlot({
                validate(
                    need(length(grep(colnames(seuratObj@meta.data),pattern ="*_snn_res.*")) > 1, "Cannot create cluster tree with only one cluster")
                )
                clustree(seuratObj,prefix = paste0(seuratObj@active.assay,"_snn_res."))+theme(legend.position = "bottom")
            })
            
            #Umap tree page
            output$DimplotTreePage <- renderPlot({
                vizu_UMAP(seuratObj, var = input$resolution_TreePage, dlabel = input$addlabels_tree, color_list = list_Color_Label,BoolCol = input$BooleanColors)
            })
            output$downloadUMAP_Tree_page<- downloadHandler(
              filename="UMAP_From_tree_page.tiff",
              content=function(file){
                tiff(file, width = 800 , height = 600,res = 100)
                print(vizu_UMAP(seuratObj,var=input$resolution_TreePage, dlabel = input$addlabels_tree, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                dev.off()
            })
            
            output$nbCellsByClust_Tree_Page <- DT::renderDataTable({
                count_by_cluster <- data.frame(table(seuratObj[[input$resolution_TreePage]]))
                colnames(count_by_cluster) <- c("Cluster",'Number of cells')
                DT::datatable(count_by_cluster,rownames = F)
            })
            #################################
            ##### DE between clusters  ######
            #################################
            updateSelectizeInput(session, "resolutionDE", choices=available_resolution, server=TRUE, selected = available_resolution[1])
            reac <- reactiveValues(test = available_resolution[1])
            observeEvent(input$resolutionDE,{
                reac$test <- input$resolutionDE
                Idents(seuratObj) <- reac$test
                Cluster_list <- levels(seuratObj)
                updateSelectizeInput(session,"Cluster1",choices = Cluster_list, server =TRUE)
                cluster_All <- c(Cluster_list,"All")
                updateSelectizeInput(session,"Cluster2", choices = cluster_All ,server = TRUE)
                
                ### Umap differential expression
                output$DimplotMarkers <- renderPlot({
                    vizu_UMAP(seuratObj,var = input$resolutionDE, dlabel = input$addlabels_DE)
                })
                output$downloadUMAP_DE_page<- downloadHandler(
                  filename="UMAP_DE_page.tiff",
                  content=function(file){
                    tiff(file, width = 800 , height = 600,res = 100)
                    print(vizu_UMAP(seuratObj,var=input$resolution_TreePage, dlabel = input$addlabels_DE))
                    dev.off()
                })
                
            })
            output$Markers <- renderDataTable(
              validate(
                need(input$submit,"Provide the markers table for one cluster versus another or all. Possibility to adjust the minimum of average of log FC as a threshold for investigate in markers. you can choose to look for only the positive markers or to keep also the negative ones")
              ))
            observeEvent(input$submit, {
                Idents(seuratObj) <- reac$test
                withProgress(message = "Differential expression", value = 0,{
                    output$Markers <- renderDataTable({
                        mList <- RunFindMarkers(seuratObj, cluster1 = input$Cluster1, cluster2 = input$Cluster2 ,posMarkers = input$onlyPos ,minimum.percent = input$percent, logthr = input$LogThreshold)
                    },server = FALSE ,extensions = c('Buttons'),
                    options = list(
                        autowidth = TRUE,
                        dom = 'frtipB',
                        buttons = list(list(extend ='csv',filename='TableMarkers'))
                    ),
                    filter = list(position = "top", clear = TRUE),)
                })
            })
            
            #####################################
            ###### Data mining of one gene ######
            #####################################
            updateSelectizeInput(session, "ResolutionDataMining", choices=available_resolution, server=TRUE)
            toremove <- c(grep(names(seuratObj@meta.data),pattern ="*_snn_res.*", value =T),grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
            varAvail <- names(seuratObj@meta.data)
            varAvail <- varAvail[-which(varAvail %in% toremove)]
            updateSelectizeInput(session, "AnnotationDataMining", choices = varAvail, server = TRUE)
            updateSelectizeInput(session, "AnnotAgainst", choices = varAvail, server = TRUE)
            gene <- seuratObj@assays$RNA@counts@Dimnames[[1]]
            updateSelectizeInput(session,"GeneList", choices=gene, server = TRUE)
            toListen <- reactive({list(input$ResolutionDataMining, input$AnnotationDataMining, input$format)})
            observeEvent(toListen(),{
              if(input$format =="SVG"){
                shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than tiff", "warning")
              }
              if(input$ResOrAnnot == "Resolution"){
                Idents(seuratObj) <- input$ResolutionDataMining
                ident_obj <- input$ResolutionDataMining
                
              }else{
                Idents(seuratObj) <- input$AnnotationDataMining
                ident_obj <- input$AnnotationDataMining
              }
              observe({
                for(i in c("downloadfeatureplot","VPdownload","HMdownload","DPdownload","downloadRawCount","downloadMatrix")){
                  shinyjs::toggleState(i,!is.null(input$GeneList))
                }
              })
              
                
              ## Create and save the featurePlot
              lapply(seq(10),function(x) # allow to create the 10 graphd
              {
                output[[paste0('FeaturePlotMultGenes',x)]] = renderPlot({FP(seuratObj,input$GeneList[x],color = input$colorPickCol, BoolColScales = input$colorScale)})
              })
              lapply(seq(length(input$GeneList)), function(x){
                output$downloadfeatureplot<- downloadHandler(
                  filename="featureplot.tiff",
                  content=function(file){
                    tiff(file, width = 800 , height = 600,res = 100)
                    print(FeaturePlot(seuratObj, features = input$GeneList))
                    dev.off()
                  })
              })
                
              ## Get heatmap from list of genes and save it 
              output$Heatmap <- renderPlot({
                  Heatmap(seuratObj, genes = input$GeneList, color_list = list_Color_Label, var = ident_obj, BoolCol = input$BooleanColors, BoolColScales = input$colorScale, colorScale = input$colorPickCol )
              })
              if(input$format =="tiff"){
                  output$HMdownload <- downloadHandler(
                      filename = "heatmap.tiff",
                      content = function(file){
                          tiff(file, width = 800 , height = 600,res = 100)
                          print(Heatmap(seuratObj, genes = input$GeneList, color_list = list_Color_Label, var = ident_obj, BoolCol = input$BooleanColors, BoolColScales = input$colorScale, colorScale = input$colorPickCol))
                          dev.off()
                      }
                  )
              }else{
                  output$HMdownload <- downloadHandler(
                      filename = "heatmap.svg",
                      content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(Heatmap(seuratObj, genes = input$GeneList,color_list = list_Color_Label, var = ident_obj, BoolCol = input$BooleanColors, BoolColScales = input$colorScale, colorScale = input$colorPickCol))
                          dev.off()
                      }
                  )
              }
              
              
              #Dotplot generator
              
              output$Dotplot <- renderPlot({
                  Dotplot(seuratObj, genes = input$GeneList, BoolColScales = input$colorScale, color = input$colorPickCol)
              })
              if(input$format =="tiff"){
                  output$DPdownload <- downloadHandler(
                      filename = "dotplot.tiff",
                      content = function(file){
                          tiff(file, width = 900 , height = 600,res = 100)
                          print(Dotplot(seuratObj, genes = input$GeneList, BoolColScales = input$colorScale, color = input$colorPickCol))
                          dev.off()
                      }
                  ) 
              }else{
                  output$DPdownload <- downloadHandler(
                      filename = "dotplot.svg",
                      content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(Dotplot(seuratObj, genes = input$GeneList, BoolColScales = input$colorScale, color = input$colorPickCol))
                          dev.off()
                      }
                  ) 
              }
              
              output$clusterOutput<- renderPlot({
                  vizu_UMAP(seuratObj, var = ident_obj, dlabel = input$addlabels_ODM, color_list = list_Color_Label,BoolCol = input$BooleanColors)
              })
              output$downloadUMAP_Cluster_resolution<- downloadHandler(
                filename="UMAP_cluster_or_annot.tiff",
                content=function(file){
                  tiff(file, width = 900 , height = 600,res = 100)
                  print(vizu_UMAP(seuratObj,var=ident_obj, dlabel = input$addlabels_ODM, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                  dev.off()
                })
              
              if(input$ResOrAnnot == "Resolution"){
                dataTableCat <- reactive(aggregate(FetchData(seuratObj, vars=input$GeneList),data.frame(seuratObj[[input$ResolutionDataMining]]),mean))
                dataMatrixEntire <- reactive(cbind(FetchData(seuratObj, vars=input$GeneList),seuratObj[[input$ResolutionDataMiningu]]))                
              }else{
                dataTableCat <- reactive(aggregate(FetchData(seuratObj, vars=input$GeneList),data.frame(seuratObj[[input$AnnotationDataMining]]),mean))
                dataMatrixEntire <- reactive(cbind(FetchData(seuratObj, vars=input$GeneList),seuratObj[[input$AnnotationDataMining]]))
              }
              output$geneExp <- DT::renderDataTable({
                validate(
                  need(input$GeneList != "","You have to choose at least one gene and here you will have for the mean of the expression for each cluster.")
                )
                DT::datatable(dataTableCat(), options = list(scrollX=TRUE), rownames = FALSE)
              })
              output$downloadRawCount<- downloadHandler(
                filename = function(){
                  paste('data_',Sys.Date(),'.csv',sep = '')
                },
                content=function(file){
                  write.csv(dataTableCat(),file)
                }
              )
              #dataMatrixEntire <- reactive(FetchData(seuratObj, vars=input$GeneList))
              output$downloadMatrix <- downloadHandler(
                filename = function(){
                  paste('Matrix_',Sys.Date(),'.csv', sep ='')
                },
                content = function(file){
                  write.csv(dataMatrixEntire(), file)
                }
              )
            })
            
           
            ListenForViolin<- reactive({list(input$AnnotAgainst, input$SplitVln,input$ResolutionDataMining, input$AnnotationDataMining)})
            observeEvent(ListenForViolin(),{
              Violins <- list()
              if(input$ResOrAnnot == "Resolution"){
                Idents(seuratObj) <- input$ResolutionDataMining
                ident_obj_vln <- input$ResolutionDataMining
                
                
              }else{
                Idents(seuratObj) <- input$AnnotationDataMining
                ident_obj_vln <- input$AnnotationDataMining
                
              }
              updateSelectizeInput(session, "labelsToKeep1", choices = levels(seuratObj) , server =T)
              ### Render ViolinPlot + saving them
              lapply(seq(10),function(x) # allow to create the 10 graphs
              {
                if(input$SplitVln == F){
                  output[[paste0('ViolinPlot',x)]] = renderPlot({SimpleViolinPlot(seuratObj,input$GeneList[x],var = ident_obj_vln, color_list = list_Color_Label, BoolCol = input$BooleanColors)})
                }else{
                  Idents(seuratObj) <- input$AnnotAgainst
                  updateSelectizeInput(session, "labelsToKeep2", choices = levels(seuratObj), server = T)
                  output[[paste0('ViolinPlot',x)]] = renderPlot({
                    gene_graph <- violinPlotSplited(seuratObj,input$GeneList[x],var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2=input$labelsToKeep2, color_list = list_Color_Label)
                    Violins[[x]] <- gene_graph
                    plot(gene_graph)
                  })
                }
              })
              observeEvent(input$format, {
                if(input$SplitVln == F && input$format== "tiff"){
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.tiff",
                    content = function(file){
                      tiff(file, width = 900 , height = 600,res = 100)
                      print(SimpleViolinPlot(seuratObj, gene = input$GeneList, var = ident_obj_vln, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                      dev.off()
                    }
                  )
                }else if (input$SplitVln == F && input$format == "SVG"){
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.svg",
                    content = function(file){
                      svg(file, width = 14 , height = 7)
                      print(SimpleViolinPlot(seuratObj, gene = input$GeneList, var = ident_obj_vln, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                      dev.off()
                    }
                  ) 
                }else if(input$SplitVln == T && input$format =="tiff"){
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.tiff",
                    content = function(file){
                      tiff(file, width = 900 , height = 600,res = 100)
                      print(ViolinForDl(seuratObj, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2, color_list = list_Color_Label))
                      dev.off()
                    }
                  )
                }else{
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.svg",
                    content = function(file){
                      svg(file, width = 14 , height = 7)
                      print(ViolinForDl(seuratObj, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2, color_list = list_Color_Label))
                      dev.off()
                    }
                  )
                }
              })
            })
        
            ######################################################
            ############# Multiple genes data mining #############
            ######################################################
            updateSelectizeInput(session,"clusterwatch",choices =available_resolution,server =TRUE)
            toremove <- c(grep(names(seuratObj@meta.data),pattern ="*_snn_res.*", value =T),grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
            varAvail <- names(seuratObj@meta.data)
            varAvail <- varAvail[-which(varAvail %in% toremove)]
            updateSelectizeInput(session, "annotationwatch", choices = varAvail, server = TRUE)
            updateSelectizeInput(session,"GeneListPool", choices=gene, server = TRUE)
            listen_page_mining_mult_genes <- reactive({list(input$ResOrAnnotMult,input$clusterwatch,input$annotationwatch,input$SumorMeans, input$GeneListPool,input$format2)})
            #reactive graph from sum or mean button
            
            observeEvent(listen_page_mining_mult_genes(),{
              if(input$format2 =="SVG"){
                  shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than tiff", "warning")
              }
              if(input$ResOrAnnotMult == "Resolution"){
                Idents(seuratObj) <- input$clusterwatch
                ident_obj_mult <- input$clusterwatch
              }else{
                Idents(seuratObj) <- input$annotationwatch
                ident_obj_mult <- input$annotationwatch
              }
              coord <- Embeddings(seuratObj[["umap"]])[,1:2]
              observe({
                for(i in c("downloadMatrixMultList","downloadfeatureplotMultGenesPooled","VPdownloadMultGenesPooled","dotplotdownloadMultGenesPooled","downloadSumGeneExp")){
                  shinyjs::toggleState(i,!is.null(input$GeneListPool))
                }
              })
              observe({
                for(i in c("downloadMatrixMultFile","downloadfeatureplotMultGenesPooledFilesInput","downloadViolinPlotMultGenesPooledFilesInput","downloadDotPlotMultGenesPooledFilesInput","downloadmeanGeneExpPooledfromFile")){
                  shinyjs::toggleState(i,!is.null(input$fileGeneList))
                }
              })
              # print and save umap resolution
              output$dimplotSeurat4page <- renderPlot({
                  vizu_UMAP(seuratObj, var = ident_obj_mult, dlabel = input$addlabels_MDM, color_list = list_Color_Label,BoolCol = input$BooleanColors)
              })
              output$downloadUMAP_resolution_page_5<- downloadHandler(
                filename="UMAP_resolution_p5.tiff",
                content=function(file){
                  tiff(file, width = 900 , height = 600,res = 100)
                  print(vizu_UMAP(seuratObj,var=ident_obj_mult, dlabel = input$addlabels_MDM, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                  dev.off()
                })
              
              poolexpTable <- reactive({
                if(input$ResOrAnnotMult == "Resolution"){
                  out <- aggregate(FetchData(seuratObj, vars=c(input$GeneListPool)),data.frame(seuratObj[[input$clusterwatch]]),mean)
                }else{
                  out <- aggregate(FetchData(seuratObj, vars=c(input$GeneListPool)),data.frame(seuratObj[[input$annotationwatch]]),mean)
                }
                out
              })
              
              output$sumGeneexp <- DT::renderDataTable({
                validate(
                  need(input$GeneListPool !="", " ")
                )
                DT::datatable(poolexpTable(), options = list(scrollX = TRUE), rownames = F)
              })
              
              output$downloadSumGeneExp<- downloadHandler(
                  filename = function(){
                      paste('data_',Sys.Date(),'.csv',sep = '')
                  },
                  content=function(file){
                    write.csv(poolexpTable(),file)
                  }
              )
              
              createEntireMatrix <- reactive({
                if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                  if (input$SumorMeans == "Sum"){
                    if(input$ResOrAnnotMult == "Resolution"){
                      expressionMatrixEntire <- cbind(geneList = rowSums(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$clusterwatch]])
                    }else{
                      expressionMatrixEntire <- cbind(geneList = rowSums(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$annotationwatch]])
                    }
                  }else{
                    if(input$ResOrAnnotMult == "Resolution"){
                      expressionMatrixEntire <- cbind(geneList = rowMeans(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$clusterwatch]])
                    }else{
                      expressionMatrixEntire <- cbind(geneList = rowMeans(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$annotationwatch]])
                    }
                  }
                }else{
                  sbt <- AddModuleScore(SeuratObjsubset, features =list(input$GeneListPool), name = "Gene_list")
                  if(input$ResOrAnnotMult == "Resolution"){
                    expressionMatrixEntire <- cbind(sbt[["Gene_list1"]],sbt[[input$clusterwatch]])
                  }else{
                    expressionMatrixEntire <- cbind(sbt[["Gene_list1"]],sbt[[input$annotationwatch]])
                  }
                }
                expressionMatrixEntire
              })
              
              ### Create FeaturePlot pooled 
              output$FeaturePlotMultGenesPooled <- renderPlot({
                validate(
                  need(!is.null(input$GeneListPool),"Choose at least one gene to display projection. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                )
                poolGene(seuratObj,gene = input$GeneListPool,typeOfNorm = input$SumorMeans, coord = coord, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_pool)
              })
              
              
              
              ## Create Violin plot pooled
              output$ViolinPlotMultGenesPooled <- renderPlot({
                validate(
                  need(!is.null(input$GeneListPool),"Choose at least one gene to display projection. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                )
                if(input$ResOrAnnotMult == "Resolution"){
                  VlnPlotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, color_list = list_Color_Label)
                }else{
                  VlnPlotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, color_list = list_Color_Label)
                }
              })
              
              
              ## Create Dotplot pooled
              output$DotPlotMultGenePooled <- renderPlot({
                validate(
                  need(!is.null(input$GeneListPool),"Choose at least one gene to display projection. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                )
                if(input$ResOrAnnotMult == "Resolution"){
                  DotplotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch,BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                }else{
                  DotplotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                }
              })
              
              
              ### Download plot ###
              if(input$format2 =="tiff"){
                  output$VPdownloadMultGenesPooled <- downloadHandler(
                      filename = "ViolinPlotMultGenesPooled.tiff",
                      content = function(file){
                          tiff(file, width = 900 , height = 600,res = 100)
                          print(
                            if(input$ResOrAnnotMult == "Resolution"){
                             VlnPlotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, color_list = list_Color_Label)
                            }else{
                              VlnPlotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, color_list = list_Color_Label)
                            }
                          )
                          dev.off()
                      }
                  )
              }else{
                  output$VPdownloadMultGenesPooled <- downloadHandler(
                      filename = "ViolinPlotMultGenesPooled.svg",
                      content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(
                            if(input$ResOrAnnotMult == "Resolution"){
                              VlnPlotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, color_list = list_Color_Label)
                            }else{
                              VlnPlotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, color_list = list_Color_Label)
                            }
                          )
                          dev.off()
                      }
                  )
              }
              
              
              
                output$downloadfeatureplotMultGenesPooled<- downloadHandler(
                    filename = "featurePlotMultGenesPooled.tiff",
                    content = function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(
                          poolGene(seuratObj,gene = input$GeneListPool,typeOfNorm = input$SumorMeans, coord = coord, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_pool)
                        )
                        dev.off()
                    }
                )  
              
              
              if(input$format2=="tiff"){
                  output$dotplotdownloadMultGenesPooled<- downloadHandler(
                      filename = "DotPlotMultGenesPooled.tiff",
                      content = function(file){
                          tiff(file, width = 900 , height = 600,res = 100)
                          print(
                            if(input$ResOrAnnotMult == "Resolution"){
                              DotplotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                            }else{
                              DotplotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                            }
                          )
                          dev.off()
                      }
                  )
              }else{
                  output$dotplotdownloadMultGenesPooled<- downloadHandler(
                      filename = "DotPlotMultGenesPooled.svg",
                      content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(
                            if(input$ResOrAnnotMult == "Resolution"){
                              DotplotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                            }else{
                              DotplotPooled(seuratObj, gene = input$GeneListPool, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                            }
                          )
                          dev.off()
                      }
                  )
              }    
            })
            
            ######################################### Gene list from file ################################
            
            # Observe event useful in order to check if all genes from gene list passed in parameters are contained in the seurat object
            observeEvent(input$fileGeneList,{
              CheckGenes <- reactive({
                validate(
                  need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.")
                )
                list_gene_file <- read.csv(input$fileGeneList$datapath, header=TRUE)
                colnames(list_gene_file) <- "cluster"
                list_gene_file <- list_gene_file$cluster
                before <- length(list_gene_file)
                list_gene_file <- list_gene_file[list_gene_file %in% rownames(seuratObj)]
                if(length(list_gene_file) != before){
                  shinyalert("Some genes cannot be found in this object, they have been removed", type = "warning")
                }
              })
              CheckGenes()
            })
            
            listen_page_mining_mult_genes_from_file <- reactive({list(input$ResOrAnnotMult,input$clusterwatch,input$annotationwatch,input$SumorMeans,input$fileGeneList, input$format2)})
            observeEvent(listen_page_mining_mult_genes_from_file(),{
              
              
              
              exprPoolFromFile <- reactive({
                list_gene_file <- Readfile(seuratObj, filePath = input$fileGeneList)
                if(input$ResOrAnnotMult == "Resolution"){
                  outFile <- aggregate(FetchData(seuratObj, vars=c(list_gene_file)), data.frame(seuratObj[[input$clusterwatch]]), mean)
                }else{
                  outFile <- aggregate(FetchData(seuratObj, vars=c(list_gene_file)), data.frame(seuratObj[[input$annotationwatch]]), mean)
                }
                outFile
              })
              
              output$meanGeneExpPooledfromFile <- DT::renderDataTable({
                DT::datatable(exprPoolFromFile(), options = list(scrollX =TRUE), rownames=F)
              })
              
              output$downloadmeanGeneExpPooledfromFile<- downloadHandler(
                filename = function(){
                  paste('data_',Sys.Date(),'.csv',sep = '')
                },
                content=function(file){
                  write.csv(exprPoolFromFile(),file)
              })
              
              
              finaleExpressionMatrixFromFile <- reactive({
                list_gene_file <- Readfile(obj = seuratObj, filePath = input$fileGeneList)
                sbt1 <- AddModuleScore(seuratObj, features =list(list_gene_file), name = "Gene_list")
                if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                  if (input$SumorMeans == "Sum"){
                    if(input$ResOrAnnotMult == "Resolution"){
                      expressionMatrixEntireFile <- cbind(geneList = rowSums(FetchData(seuratObj, vars =c(list_gene_file))),seuratObj[[input$clusterwatch]])
                    }else{
                      expressionMatrixEntireFile <- cbind(geneList = rowSums(FetchData(seuratObj, vars =c(list_gene_file))),seuratObj[[input$annotationwatch]])
                    }
                  }else{
                    if(input$ResOrAnnotMult == "Resolution"){
                      expressionMatrixEntireFile <- cbind(geneList =rowMeans(FetchData(seuratObj, vars =c(list_gene_file))),seuratObj[[input$clusterwatch]])
                    }else{
                      expressionMatrixEntireFile <- cbind(geneList =rowMeans(FetchData(seuratObj, vars =c(list_gene_file))),seuratObj[[input$annotationwatch]])
                    }
                  }
                }else{
                  if(input$ResOrAnnotMult == "Resolution"){
                    expressionMatrixEntireFile <- cbind(sbt1[["Gene_list1"]],sbt1[[input$clusterwatch]])
                  }else{
                    expressionMatrixEntireFile <- cbind(sbt1[["Gene_list1"]],sbt1[[input$annotationwatch]])
                  }
                }
                expressionMatrixEntireFile
              })
              
              output$downloadMatrixMultFile <- downloadHandler(
                filename = function(){
                  paste('Entire_matrix_',Sys.Date(),'.csv',sep = '')
                },
                content=function(file){
                  write.csv(finaleExpressionMatrixFromFile(),file)
                }
              )
              
              
              
              coord <- Embeddings(seuratObj[["umap"]])[,1:2]
              output$FeaturePlotMultGenesPooledFilesInput <- renderPlot({
                validate(
                  need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.This plot will show the pool of each gene expression for each gene id contained in the file. \nExample of a file :\nHeader\nMlxipl\nSox9\nEtc...  ")
                )
                poolGene(seuratObj, typeOfNorm = input$SumorMeans, coord = coord, file = input$fileGeneList , BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_file)
              })
              
              output$ViolinPlotMultGenesPooledFilesInput <- renderPlot({
                validate(
                  need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.")
                )
                gene_list <- Readfile(obj = seuratObj, filePath = input$fileGeneList)
                if(input$ResOrAnnotMult == "Resolution"){
                  VlnPlotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, file = input$fileGeneList, color_list = list_Color_Label)
                }else{
                  VlnPlotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, file = input$fileGeneList, color_list = list_Color_Label)
                }
              })
              
              
              output$DotPlotMultGenePooledFromFile <- renderPlot({
                validate(
                  need(input$fileGeneList, "Choose a file that contain a list of genes to display the projection.")
                )
                if(input$ResOrAnnotMult == "Resolution"){
                  DotplotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, file = input$fileGeneList, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                }else{
                  DotplotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, file = input$fileGeneList, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                }
              })
              
              
              
              if(input$format2=="tiff"){
                output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                  filename = "ViolinPlotMultGenesPooled.tiff",
                  content = function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(
                      if(input$ResOrAnnotMult == "Resolution"){
                        VlnPlotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, file = input$fileGeneList, color_list = list_Color_Label)
                      }else{
                        VlnPlotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, file = input$fileGeneList, color_list = list_Color_Label)
                      }
                    )
                    dev.off()
                  }
                ) 
              }else{
                output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                  filename = "ViolinPlotMultGenesPooled.svg",
                  content = function(file){
                    svg(file, width = 14 , height = 7)
                    print(
                      if(input$ResOrAnnotMult == "Resolution"){
                        VlnPlotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, file = input$fileGeneList, color_list = list_Color_Label)
                      }else{
                        VlnPlotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, file = input$fileGeneList, color_list = list_Color_Label)
                      }
                    )
                    dev.off()
                  }
                )
              }
              
            
              output$downloadfeatureplotMultGenesPooledFilesInput <- downloadHandler(
                filename = "featurePlotFilesInput.tiff",
                content = function(file){
                  tiff(file, width = 900 , height = 600,res = 100)
                  print(poolGene(seuratObj, typeOfNorm = input$SumorMeans, coord = coord, file = input$fileGeneList, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled, size.points = input$sizePoint_file))
                  dev.off()
                }
              )  
              
              
              if(input$format2 == "tiff"){
                output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                  filename = "DotPlotMultGenesPooled.tiff",
                  content = function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(
                      if(input$ResOrAnnotMult == "Resolution"){
                        DotplotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, file = input$fileGeneList, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                      }else{
                        DotplotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, file = input$fileGeneList, BoolColScale = input$colorScalePooled, colorScale = input$colorPickColPooled)
                      }
                    )
                    dev.off()
                  }
                )
              }else{
                output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                  filename = "DotPlotMultGenesPooled.svg",
                  content = function(file){
                    svg(file, width = 14 , height = 7)
                    print(
                      if(input$ResOrAnnotMult == "Resolution"){
                        DotplotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, cluster = input$clusterwatch, file = input$fileGeneList)
                      }else{
                        DotplotPooled(seuratObj, typeOfNorm = input$SumorMeans, annotOrRes = input$ResOrAnnotMult, annotation = input$annotationwatch, file = input$fileGeneList)
                    })
                    dev.off()
                  }
                )
              }
              
            })
            ###########################################
            ########  Information extraction #########
            ###########################################
            updateSelectizeInput(session,"clusterNumber",choices = available_resolution,server =TRUE)
            toremove <- c(grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
            varAvail <- names(seuratObj@meta.data)
            varAvail <- varAvail[-which(varAvail %in% toremove)]
            updateSelectizeInput(session,"chooseVar2Plot", choices=varAvail, server=TRUE)
            
            cluster <- reactiveValues(corresponding_res = "*_snn_res.0") # Need a reactive value to pass through the observe event 
            observeEvent(input$clusterNumber,{
                Idents(seuratObj) <- input$clusterNumber
                DMclusters<- c(levels(seuratObj), "All")
                updateSelectizeInput(session,"clusterInfo",choices = DMclusters, server =TRUE)
                output$Dimplot_cluster <- renderPlot({
                    vizu_UMAP(seuratObj,var = input$clusterNumber, dlabel = input$addlabels_info, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                })
                output$downloadUMAP_info<- downloadHandler(
                  filename = "UMAP_info.tiff",
                  content = function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(seuratObj,var = input$clusterNumber, dlabel = input$addlabels_info, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                    dev.off()
                  }
                )
                toListen <- reactive({list(input$clusterInfo,input$chooseVar2Plot, input$stackedBP, input$FreqOrVal)})
                observeEvent(toListen(), {
                  
                    output$ggplot_information <- renderPlot({
                        p <- Plot2Render(seuratObj, cluster = input$clusterInfo, StackedBarPlot = input$stackedBP, VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, freqOrValues = input$FreqOrVal, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                        p +theme(legend.position = "none",axis.text.x = element_text(angle =60,hjust=1))
                    })
                    output$legend_information <- renderPlot({
                        p <- Plot2Render(seuratObj, cluster = input$clusterInfo, StackedBarPlot = input$stackedBP, VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, freqOrValues = input$FreqOrVal, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                        grid.newpage()
                        legend <- cowplot::get_legend(p)
                        grid.draw(legend)
                    })
                    
                    output$cluster_percent <- renderPlot({
                      g <- percentRender(seuratObj,VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                      g
                    })
                    output$downloadInformationPlot2 <- downloadHandler(
                      filename = "Percent_Plot.tiff",
                      content = function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(percentRender(seuratObj,VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                        dev.off()
                      }
                    )
                    
                    output$get_info <- DT::renderDataTable({
                        Table2Render(seuratObj, cluster = input$clusterInfo, VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber)
                    })
                    output$downloadInformationPlot <- downloadHandler(
                        filename = "Information_Plot.tiff",
                        content = function(file){
                            tiff(file, width = 900 , height = 600,res = 100)
                            print(Plot2Render(seuratObj, cluster = input$clusterInfo, StackedBarPlot = input$stackedBP, VariablePlot = input$chooseVar2Plot, NbCluster = input$clusterNumber, freqOrValues = input$FreqOrVal, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                            dev.off()
                        }
                    )
                    output$downloadMatrixInformation <- downloadHandler(
                        filename = "Information_Matrix.csv",
                        content = function(file){
                            write.csv(Table2Render(),file)
                        }
                    )
                })
                
                  
                
            })
            #######################################
            ########### Add Annotation ############
            #######################################
            removeInfo <- c(grep(names(seuratObj@meta.data),pattern ="*_snn_res.*", value =T),grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score")
            list.metadata <- colnames(seuratObj@meta.data)
            list.metadata <- list.metadata[-which(list.metadata %in% removeInfo)]
            updateSelectizeInput(session,"Annot2Complete", choices =list.metadata)
            updateSelectizeInput(session,"res4Annot",choices = available_resolution,server =TRUE)
            resChoice <- NULL
            observeEvent(input$res4Annot,{
                resChoice <<- input$res4Annot
                Idents(seuratObj) <- input$res4Annot
                clusterSbt <- levels(seuratObj)
                updateSelectizeInput(session,"cluster4Annot",choices = clusterSbt, server = TRUE)
                output$dimplot_res_annot <- renderPlot({
                    vizu_UMAP(seuratObj,var = input$res4Annot, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                })
                output$downloadUMAP_resolution_annot_page<- downloadHandler(
                  filename="UMAP_resolution_annot_page.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(seuratObj,var=input$res4Annot, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                    dev.off()
                  })
                
            })
            
            ##render plot displaying variable of interest when use update instead of create
            output$condPlot <- renderPlot({
                vizu_UMAP(seuratObj, var = input$Annot2Complete, dlabel = input$labelsconditional, color_list = list_Color_Label, BoolCol = input$BooleanColors)
            })
            output$downloadUMAP_Conditional<- downloadHandler(
              filename="UMAP_conditional.tiff",
              content=function(file){
                tiff(file, width = 900 , height = 600,res = 100)
                print(vizu_UMAP(seuratObj,var=input$Annot2Complete, dlabel = input$labelsconditional, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                dev.off()
              })
            
            observeEvent(input$submitAnnot,{
                if(input$choicesCreateOrComplete =="Create"){
                    if(input$AnnotName == "" | input$clusterName == "" | length(input$cluster4Annot) == 0){
                        shinyalert("Oops", "The name of the column and the name of the annotation you want to give has to be filled", type = "error")
                    }else{
                      for(i in c("downloadUMAP_post_annotation","labelsPostAnnot")){
                        shinyjs::show(i)
                      }
                      if(is.element(input$clusterName, names(seuratObj@meta.data))){
                        shinyalert("Warning","The name that you give to your annotation is already existing in your dataset. Do you want to continue ?", type = "warning", showConfirmButton = TRUE, showCancelButton = TRUE,callbackR = function(x){
                          if(x == TRUE){
                            annot <- Annotation_existent(seuratObj, nameVar = input$clusterName, annotation = input$AnnotName, resolution = resChoice,clusterToAnnotate = input$cluster4Annot)
                            seuratObj[[input$clusterName]] <<- annot
                            shinyalert("Done","You have correctly annotated your cluster", type = "success")
                            output$postAnnotation <- renderPlot({
                              vizu_UMAP(seuratObj, var = input$clusterName, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                            }) 
                            }
                          } 
                        )
                      }else{
                        annot <- Annotation_unexistent(seuratObj, annotation = input$AnnotName, resolution = resChoice,clusterToAnnotate = input$cluster4Annot)
                        seuratObj[[input$clusterName]] <<- annot
                        shinyalert("Done","You have correctly annotated your cluster", type = "success")
                        output$postAnnotation <- renderPlot({
                          vizu_UMAP(seuratObj, var = input$clusterName, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                        }) 
                      }
                      output$downloadUMAP_post_annotation<- downloadHandler(
                        filename="UMAP_post_annot.tiff",
                        content=function(file){
                          tiff(file, width = 900 , height = 600,res = 100)
                          print(vizu_UMAP(seuratObj,var=input$Annot2Complete, dlabel = input$labelsAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                          dev.off()
                      })
                      toremove <- c(grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),grep(names(seuratObj@meta.data),pattern ="snn_res.*", value =T),"percent.mt","S.Score", "G2M.Score")
                      choices <- names(seuratObj@meta.data)
                      choices_annot <- choices[-which(choices %in% toremove )]
                      updateSelectizeInput(session,"chooseVar2Plot", choices=choices_annot)
                      updateSelectizeInput(session,"whichAnnot", choices=choices_annot)
                      updateSelectizeInput(session,"Annot2Complete", choices =choices_annot)
                      updateSelectizeInput(session, "AnnotationDataMining", choices = choices_annot)
                      updateSelectizeInput(session,"annotationwatch",choices = choices_annot)
                      updateSelectizeInput(session,"AnnotAgainst", choices = choices_annot)
                      updateSelectizeInput(session,"InfoToPlot", choices = choices_annot)  
                    }
                }else{
                    if(input$AnnotName == "" | length(input$cluster4Annot) == 0){
                        shinyalert("Oops", "The name of the annotation and the cluster(s) you want to give has to be filled", type = "error")
                    }else{
                      annot <- Annotation_existent(seuratObj,nameVar = input$Annot2Complete, annotation = input$AnnotName, resolution = resChoice,clusterToAnnotate = input$cluster4Annot)
                      seuratObj[[input$Annot2Complete]] <<- annot
                      shinyalert("Done","You have correctly annotated your cluster", type = "success")
                    }
                    output$postAnnotation <- renderPlot({
                        vizu_UMAP(seuratObj, var = input$Annot2Complete, dlabel = input$labelsPostAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                    })
                    output$downloadUMAP_post_annotation<- downloadHandler(
                      filename="UMAP_post_annot.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(seuratObj,var=input$Annot2Complete, dlabel = input$labelsPostAnnot, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                        dev.off()
                    })
                }
                
            })
            output$downloadRDSwithNewAnnot <- downloadHandler(
                filename = "SeuratObj.rds",
                content = function(file){
                    saveRDS(seuratObj,file)
                }
            )
            ######################################
            ######## Automatic annotation ########
            ######################################
            #shinyjs::disable("downloadProbabilities")
            gene <- seuratObj@assays$RNA@counts@Dimnames[[1]]
            observe({
              shinyjs::toggleState("AutomaticAnnotation",input$slotNameAnnotScina != "")
            })
            listbutton <- reactive({list(input$NbCellType)})
            observeEvent(listbutton(),{
              output$testPanel <- renderUI({
                tagList(
                  lapply(seq(1,input$NbCellType), function(i){
                    tagList(
                      textInput(paste("cellType", i, sep="_"), "Give a name to your cell type"),
                      selectizeInput(paste("markerCellType",i,sep="_"), "Give some markers (between 2 and 20 is advised)",choices=character(0), multiple = T)
                    )
                  }) 
                )
                
              })
              updateSelectizeInput(session,paste("markerCellType",1,sep="_"),choices=gene, server = TRUE)
              lapply(seq(2,input$NbCellType), function(i){
                updateSelectizeInput(session,paste("markerCellType",i,sep="_"),choices=gene, server = TRUE)
              })
            })
            output$FileExplanation <- renderText(
             HTML(paste("If you want to pass a file here is an example of file :","Cell_type 1;Cell_type 2","Marker1_CT1;Marker1_CT2","Marker2_CT1;Marker2_CT2","Marker3_CT1;Marker3_CT2","Marker4_CT1;","Marker5_CT1;","<b>Important notice</b> :", "Do not use space to label your cell type","The list can have different numbers of markers", sep ="<br/>"))
            )
            observeEvent(input$AutomaticAnnotation,{
              shinyjs::disable(id = input$AutomaticAnnotation)
              if(input$typeOfList == "List"){
                
                cell_type <- reactive({
                  lapply(seq(1,input$NbCellType), function(i) {
                    input[[paste("cellType", i, sep="_")]]
                  })
                })
                marker <- reactive({
                  lapply(seq(1,input$NbCellType), function(i) {
                    input[[paste("markerCellType", i, sep="_")]]
                  })
                })
                scina_list <- list()
                marker_list <- marker()
                celltype_list <- cell_type()
                for(i in 1:length(marker_list)){
                  tmp <- unlist(celltype_list[i])
                  scina_list[tmp] <- marker_list[i]
                }
              }else{
                scina_list <- read.csv(input$InputMarkerFile$datapath,sep = ";")
              }
              scina_results<- useSCINA(object = seuratObj, list_annotation = scina_list,allowUnknown = TRUE, allowOverlap = TRUE)
              seuratObj <<- Fill_proba(object = seuratObj, slotName = input$slotNameAnnotScina, scina = scina_results)
              proba_avail <- c(grep(names(seuratObj@meta.data),pattern ="*_probabilities", value =T))
              shinyjs::enable(input$AutomaticAnnotation)
              shinyjs::show("Scina_annot_DL")
              shinyjs::show("Probabilities_Dl")
              shinyjs::show("Proba_selection")
              shinyjs::enable("downloadProbabilities")
              updateSelectizeInput(session,"Proba_selection", choices =proba_avail)
              output$Results_annotation <- renderPlot(
                Dimplot_with_ggplot(object = seuratObj, group = input$slotNameAnnotScina)
              )
              output$Probabilities <- renderPlot(
                FeaturePlot_proba(seuratObj, proba_col = input$Proba_selection, slotNameAnnot = input$slotNameAnnotScina)
              )
              output$Scina_annot_DL <- downloadHandler(
                filename = "Scina_annotation_results.tiff",
                content=function(file){
                  tiff(file, width = 900 , height = 600,res = 100)
                  print(Dimplot_with_ggplot(object = seuratObj, group = input$slotNameAnnotScina))
                  dev.off()
                }
              )
              
              output$Probabilities_Dl <- downloadHandler(
                filename = "Scina_annotation_probabilities.tiff",
                content=function(file){
                  tiff(file, width = 900 , height = 600,res = 100)
                  print(FeaturePlot_proba(seuratObj, proba_col = input$Proba_selection, slotNameAnnot = input$slotNameAnnotScina))
                  dev.off()
                }
              )
              output$downloadProbabilities <- downloadHandler(
                filename = "Scina_annotation_probabilities.csv",
                content=function(file){
                  write.csv(t(scina_results$probabilities),file)
                }
              )
              gc()
              toremove <- c(grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),grep(names(seuratObj@meta.data),pattern ="snn_res.*", value =T),"percent.mt","S.Score", "G2M.Score" ,grep(names(seuratObj@meta.data),pattern ="*_probabilities", value =T))
              choices <- names(seuratObj@meta.data)
              choices_annot <- choices[-which(choices %in% toremove )]
              updateSelectizeInput(session,"chooseVar2Plot", choices=choices_annot)
              updateSelectizeInput(session,"whichAnnot", choices=choices_annot)
              updateSelectizeInput(session,"Annot2Complete", choices =choices_annot)
              updateSelectizeInput(session, "AnnotationDataMining", choices = choices_annot)
              updateSelectizeInput(session,"annotationwatch",choices = choices_annot)
              updateSelectizeInput(session,"AnnotAgainst", choices = choices_annot)
              updateSelectizeInput(session,"InfoToPlot", choices = choices_annot)
            })
            ######################################
            ########### Subclustering ############
            ######################################
            removeInfo <- c(grep(names(seuratObj@meta.data),pattern ="*_snn_res.*", value =T),grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),"percent.mt","S.Score", "G2M.Score","seurat_clusters")
            annotation_sub <- colnames(seuratObj@meta.data)
            annotation_sub <- annotation_sub[-which(annotation_sub %in% removeInfo)]
            updateSelectizeInput(session,"clusterResPage7",choices = available_resolution,server =TRUE)
            updateSelectizeInput(session,"whichAnnot", choices = annotation_sub, server = TRUE)
            listenbutton <- reactive({list(input$clusterResPage7, input$whichAnnot,input$subcluster)})
            subsetSeuratObj <- NULL
            shinyjs::hide("NameObject")
            shinyjs::hide("downloadSubRds")
            shinyjs::hide("downloadLog")
            observeEvent(listenbutton(),{
                test <- NULL #test will contain the different cells that are highlighted with the clusters/annotations that we want to keep
                if(input$subcluster == "Cluster"){
                    Idents(seuratObj) <- input$clusterResPage7
                    clusterSbt <- levels(seuratObj)
                    updateSelectizeInput(session,"subclusterize",choices = clusterSbt, server = TRUE )
                    output$Dimplot_subclustering <- renderPlot({
                        vizu_UMAP(seuratObj,var = input$clusterResPage7, dlabel=input$labelsSubset, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                    })
                    output$downloadUMAP_sbt<- downloadHandler(
                      filename="UMAP_resolution_for_sbt.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(seuratObj,var=input$clusterResPage7, dlabel = input$labelsSubset, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                        dev.off()
                      })
                    observeEvent(input$subclusterize,{
                      for (i in c("downloadUMAPcellKept","cellsHighlightSize","cellsSize")){
                        shinyjs::show(i)
                      }
                        test <- NULL
                        if(length(input$subclusterize) == length(clusterSbt)){
                            output$Dimplot_kept <- renderPlot({
                              MakeUMAP_Red(seuratObj, highliht_cells = colnames(seuratObj),sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                            })
                            output$downloadUMAPcellKept<- downloadHandler(
                              filename="UMAP_cell_kept.tiff",
                              content=function(file){
                                tiff(file, width = 900 , height = 600,res = 100)
                                print(MakeUMAP_Red(seuratObj,test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                dev.off()
                            })
                        }else{
                            test <- colnames(seuratObj)[which(seuratObj[[]][input$clusterResPage7] == input$subclusterize[1])]
                            output$Dimplot_kept <- renderPlot({
                                MakeUMAPhighlight(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                            })
                            output$downloadUMAPcellKept<- downloadHandler(
                              filename="UMAP_cell_kept.tiff",
                              content=function(file){
                                tiff(file, width = 900 , height = 600,res = 100)
                                print(MakeUMAPhighlight(seuratObj,test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                dev.off()
                              })
                        }
                        if(length(input$subclusterize) > 1){
                            if(length(input$subclusterize) == length(clusterSbt)){
                                test <- colnames(seuratObj)
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAP_Red(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.tiff",
                                  content=function(file){
                                    tiff(file, width = 900 , height = 600,res = 100)
                                    print(MakeUMAP_Red(seuratObj,test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                    dev.off()
                                })
                            }
                            else{
                                for(i in 2:length(input$subclusterize)){
                                    test <- c(test, colnames(seuratObj)[which(seuratObj[[]][input$clusterResPage7] == input$subclusterize[i])])
                                    output$Dimplot_kept <- renderPlot({
                                      MakeUMAPhighlight(seuratObj, highliht_cells =  test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                    })
                                    output$downloadUMAPcellKept<- downloadHandler(
                                      filename="UMAP_cell_kept.tiff",
                                      content=function(file){
                                        tiff(file, width = 900 , height = 600,res = 100)
                                        print(MakeUMAPhighlight(seuratObj,test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                        dev.off()
                                    })
                                }
                            }
                        }
                        
                        subsetSeuratObj <<- subset(seuratObj, cells = test)
                    })
                }else{
                    Idents(seuratObj) <- input$whichAnnot
                    annotSbt <- levels(seuratObj)
                    updateSelectInput(session, "subannot", choices = annotSbt)
                    output$Dimplot_subclustering <- renderPlot({
                        vizu_UMAP(seuratObj,var = input$whichAnnot, dlabel = input$labelsSubset, color_list = list_Color_Label, BoolCol = input$BooleanColors)
                    })
                    
                    output$downloadUMAP_sbt<- downloadHandler(
                      filename="UMAP_resolution_for_sbt.tiff",
                      content=function(file){
                        tiff(file, width = 900 , height = 600,res = 100)
                        print(vizu_UMAP(seuratObj,var=input$whichAnnot, dlabel = input$labelsSubset, color_list = list_Color_Label, BoolCol = input$BooleanColors))
                        dev.off()
                      })
                    observeEvent(input$subannot,{
                      for (i in c("downloadUMAPcellKept","cellsHighlightSize","cellsSize")){
                        shinyjs::show(i)
                      }
                      test <- NULL
                      if(length(input$subannot) == length(annotSbt)){
                          test <- colnames(seuratObj)
                          output$Dimplot_kept <- renderPlot({
                            MakeUMAP_Red(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                          })
                          output$downloadUMAPcellKept<- downloadHandler(
                            filename="UMAP_cell_kept.tiff",
                            content=function(file){
                              tiff(file, width = 900 , height = 600,res = 100)
                              print(MakeUMAP_Red(seuratObj,test, sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                              dev.off()
                            })
                          
                      }else{
                          test <- colnames(seuratObj)[which(seuratObj[[]][input$whichAnnot] == input$subannot[1])]
                          output$Dimplot_kept <- renderPlot({
                            MakeUMAPhighlight(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                          })
                          output$downloadUMAPcellKept<- downloadHandler(
                            filename="UMAP_cell_kept.tiff",
                            content=function(file){
                              tiff(file, width = 900 , height = 600,res = 100)
                              print(MakeUMAPhighlight(seuratObj,test))
                              dev.off()
                            })
                      }
                      if(length(input$subannot) > 1){
                          if(length(input$subannot) == length(annotSbt)){
                              test <- colnames(seuratObj)
                              output$Dimplot_kept <- renderPlot({
                                  MakeUMAP_Red(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                              })
                              output$downloadUMAPcellKept<- downloadHandler(
                                filename="UMAP_cell_kept.tiff",
                                content=function(file){
                                  tiff(file, width = 900 , height = 600,res = 100)
                                  print(MakeUMAP_Red(seuratObj,test))
                                  dev.off()
                                })
                          }else{
                              for(i in 2:length(input$subannot)){
                                  test <- c(test, colnames(seuratObj)[which(seuratObj[[]][input$whichAnnot] == input$subannot[i])])
                                  output$Dimplot_kept <- renderPlot({
                                      MakeUMAPhighlight(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                  })
                                  output$downloadUMAPcellKept<- downloadHandler(
                                    filename="UMAP_cell_kept.tiff",
                                    content=function(file){
                                      tiff(file, width = 900 , height = 600,res = 100)
                                      print(MakeUMAPhighlight(seuratObj,test))
                                      dev.off()
                                    })
                              }
                          }
                      }
                      subsetSeuratObj <<- subset(seuratObj, cells = test) # the double <<- is used to fill a value which is in observeEvent but that we want to keep for after
                      
                    })
                }
                
            })
            observe({
                shinyjs::toggleState("dosubcluster",input$subcluster == "Annotation" && input$subannot != "" || input$subcluster == "Cluster"  && input$subclusterize !="")
            })
            observeEvent(input$dosubcluster,{
                shinyjs::disable("dosubcluster")
                withProgress(message = "Renormalizing data", value = 0,{
                    incProgress(1/4,message ="Running PCA")
                    if(length(grep("SCT",names(subsetSeuratObj@assays))) == 0){
                        tryCatch(subsetSeuratObj <<- RunPCA(subsetSeuratObj),#Need tryCatch to avoid error when there is no findvariablefeatures 
                                 error = function(c){
                                     shinyalert("warning","Can't find variable features in your RDS need to run FindVariableFeatures()", type = "warning")
                                     subsetSeuratObj <<- FindVariableFeatures(subsetSeuratObj)
                                     subsetSeuratObj <<- ScaleData(subsetSeuratObj)
                                      
                                 },
                                 finally = subsetSeuratObj <<- RunPCA(subsetSeuratObj) )
                    }else{
                        subsetSeuratObj <<- RunPCA(subsetSeuratObj)
                    }
                    incProgress(1/4,message ="Finish running PCA")
                    subsetSeuratObj <<- FindNeighbors(subsetSeuratObj)
                    incProgress(1/4,message ="Running UMAP")
                    subsetSeuratObj <<- RunUMAP(subsetSeuratObj, dims = 1:30)
                    incProgress(1/4,message ="Finish running UMAP")
                })
                withProgress(message ="Find cluster", value=0,{
                    for(i in seq(0,1,0.1)){
                        incProgress(i*10/10,message =paste0("res ",i*10,"/10"))
                        subsetSeuratObj <<- FindClusters(subsetSeuratObj, res =i)
                    }
                    
                })
                output$dimplot_aftersbt <- renderPlot({
                    vizu_UMAP(subsetSeuratObj, var = "orig.ident" ,dlabel = input$labelsAfterSubset)
                })
                
                for(i in c("downloadUMAPaftersbt", "labelsAfterSubset")){
                  shinyjs::show(i)
                }
                output$downloadUMAPaftersbt<- downloadHandler(
                  filename="UMAP_after_sbt.tiff",
                  content=function(file){
                    tiff(file, width = 900 , height = 600,res = 100)
                    print(vizu_UMAP(subsetSeuratObj,var=NULL, dlabel = input$labelsAfterSubset))
                    dev.off()
                })
                shinyjs::enable("dosubcluster")
                shinyjs::show("NameObject")
                shinyjs::show("downloadSubRds")
                shinyjs::show("downloadLog")
            })
            observe({
                shinyjs::toggleState("downloadSubRds",input$NameObject != "")
            })
            output$downloadSubRds <- downloadHandler(
                filename = function(){
                    paste0(input$NameObject,".rds")
                },
                content = function(file){
                    shinyalert("warning", "Seurat object can be heavy, and it can take time, please be patient, don't click multiple times", "warning")
                    saveRDS(subsetSeuratObj,file)
                }
                
            )
            observe({
                shinyjs::toggleState("downloadLog",input$NameObject != "")
            })
            
            output$downloadLog <- downloadHandler(
              filename = function(){
                paste0(input$NameObject,"_log.html")
              },
              content = function(file){
                write_rmarkdown_report_subclustering(name_file = file, ClustOrAnnot = input$subcluster, cluster = input$clusterResPage7, subsetCluster = input$subclusterize, Annot = input$whichAnnot, subsetAnnot = input$subannot)
              }
            )
        }
    })
}