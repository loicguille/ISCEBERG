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
    shinyjs::hide("Refresh")
    ##### read survey to give the rds available on the application ######
    output$Isceberg <- renderUI({
      img(src = "iceberg2.png", height = "90%", width ="90%")
    })
    output$memoryCons <- renderUI({
        img(src = "nb_giga.jpeg", height = "90%", width ="90%")
        })
    output$timeCons <- renderUI({
        img(src = "time.jpeg", height = "100%", width = "100%")
        })
    #output$Help_file <- renderUI({
    #  includeHTML("Help.html")
    #})
    observe({
        shinyjs::toggleState("createSeurat", !is.null(input$InputFile) && input$InputFile != "" ||  !is.null(input$MatrixFile) && input$MatrixFile != "" ||  !is.null(input$listFiles) && input$listFiles!= "")#allow to grey a button
    })
    observeEvent(input$file,{
        if(input$file =="RDS"){
            shinyalert("Information", "If you choose to give your proper RDS dataset, this has to be an output from seurat analyzed by yourself", type = "info")
        }
    })
    observeEvent(input$Refresh,{
        shinyalert("warning", "Do this operation will reset all what you have done before", type = "warning", showCancelButton = TRUE, callbackR = function(x){if(x == TRUE){shinyjs::refresh()}})
        
    })

    #### Little function
    FP <- function(object, gene){ # function that allow the construction of featureplot recursively
      validate(
        need(gene != "","Select at least one gene")
      )
      FeaturePlot(object, features = gene)
      
    }
    
    #Violin plot generator
    vnPlot <- function(object,gene){ #Take in entry seurat object and one gene, function used in data mining with one gene
      validate(
        need(gene != "","Choose at least one gene and it will show you the violin plot of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
      )
      VlnPlot(object, features = gene, pt.size = 0)
    }
    
    vizu_UMAP <- function(obj, var = NULL, dlabel){
      DimPlot(obj,group.by = var, label = dlabel, label.size = 6)
    }
    
    nbCellsbydt <- function(obj){
      cells_by_dt <- data.frame(table(obj$orig.ident))
      ggplot(cells_by_dt,aes(Var1,Freq, fill=Var1))+geom_bar(stat = "identity")+ggtitle("Number of cells by datasets")+geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25)
    }
    
    UmapNbGenesCell <- function(obj){
      coord_Umap <- Embeddings(obj[["umap"]])[,1:2]
      expressed_cells <- cbind.data.frame(counts=colSums(obj@assays$RNA@data > 0),coord_Umap)
      ggplot(expressed_cells, aes(UMAP_1,UMAP_2, color = counts))+geom_point(shape =20, size=0.25)+scale_color_gradient(low="lightgrey", high="blue")+theme_classic()+ggtitle("Number of expressed genes by cells")
    }
    
    makeVlnGreatAgain <- function(obj,var, grouping, col = 1){
      VlnPlot(object = obj, features = var, group.by = grouping, ncol = col, pt.size = 0 )
    }
    MakeUMAPhighlight <- function( obj, highliht_cells , sizePoint = 0.5, sizeHighlight = 0.25){
      DimPlot(obj, cells.highlight = highliht_cells, order = T,pt.size = sizePoint,sizes.highlight = sizeHighlight)+scale_color_manual(labels = c("Not kept","Kept"), values = c("gray","red"))+theme(legend.text = element_text(size = 10))
    }
    MakeUMAP_Red <- function( obj, highliht_cells , sizePoint = 0.5, sizeHighlight = 0.25){
      DimPlot(obj, cells.highlight = highliht_cells, order = T,pt.size = sizePoint,sizes.highlight = sizeHighlight)+scale_color_manual(labels = c("Kept"), values = c("red"))+theme(legend.text = element_text(size = 10))
    }
    MakeComplicatedVln <- function(obj, Giventitle, var){
      nfeat <- VlnPlot(obj,group.by = var, features = c("nFeature_RNA"), pt.size = 0)+theme(legend.position ="none")+geom_hline(yintercept = c(input$minGenePerCells,input$maxGenePerCells), color = c("red", "red"))
      ncount <- VlnPlot(obj,group.by = var, features = c("nCount_RNA"), pt.size = 0)+theme(legend.position ="none")+geom_hline(yintercept = input$maxCountPerCells, color = "red ")
      pmt <- VlnPlot(obj,group.by = var, features = c("percent.mt"), pt.size = 0)+theme(legend.position ="none")+geom_hline(yintercept = input$percentMito, color = "red")
      grid.arrange(nfeat,ncount,pmt,nrow=1, top = textGrob(Giventitle, gp =gpar(col="black", fontsize =20, fontface ="bold"))) ## arrange all the plot between them
    }
    Plotnbcellsbeforeafter <- function(objBefore, objAfter){
      dfNbcells <- data.frame(time = c("before_filtering", "after_filtering"),nb_cells = c(dim(objBefore)[2],dim(objAfter)[2]))
      ggplot(dfNbcells,aes(x= time, y = nb_cells, fill = time))+geom_bar(stat = "identity")+scale_x_discrete(limits=c("before_filtering", "after_filtering"))+ggtitle("Impact of filtering on the number of cells")+geom_text(aes(label=nb_cells), position=position_dodge(width=0.9), vjust=-0.25)
    }
    
    violinPlot <- function(seuratObj,genes,var1,var2, conditionVar2, conditionVar1){
      validate(need(var1 != var2 & conditionVar2 != "" & conditionVar1 != "" & genes != "","Need that the two variable your are looking for are different and that the list of the condition you want to look at aren't empty"))
      unfiltered_table <- FetchData(seuratObj, vars = c(genes, var1, var2))
      filtered_table <- data.frame()
      for(i in 1:length(conditionVar1)){
        for(j in 1:length(conditionVar2)){
          filtered_table <- rbind(filtered_table,unfiltered_table %>% filter(unfiltered_table[var1] == conditionVar1[i] & unfiltered_table[var2] == conditionVar2[j]))
        }
        j=1
      }
      
      colnames(filtered_table) <- c("gene","variable1", "variable2")
      graph <- ggplot(filtered_table)+geom_violin(aes(variable1,gene, fill = variable2), scale = "width") + xlab("Identity")+ylab("Expression Level") +guides(fill = guide_legend(title=var2))+theme_classic()+ggtitle(genes)+theme(plot.title = element_text(face="bold",hjust = 0.5))
      return(graph)
    }
    
    ViolinForDl <- function(seuratObject, geneList, var1, var2, conditionVar1, conditionVar2 ){## encapsulation of downloadable plot for split object in rendering Violinplot
      Violins <- list()
      for(i in 1:length(geneList)){
        Violins[[i]] <- violinPlot(seuratObj = seuratObject, genes = geneList[i],var1 = var1, var2 = var2, conditionVar1 = conditionVar1, conditionVar2 = conditionVar2)
      }
      do.call(grid.arrange,Violins)
    }
    #####################################################################################################################################################
    ############################### Create or Upload Seurat object ######################################################################################
    #####################################################################################################################################################
    observeEvent(input$createSeurat ,{
        shinyjs::hide("createSeurat") #Hide all buttons from first page so the user will need to reload the page if he want to see another file
        shinyjs::show("Refresh")
        shinyjs::disable("Refresh")
        shinyjs::hide("choiceFile")
        shinyjs::hide("file")
        shinyjs::hide("listFiles")
        shinyjs::hide("header")
        shinyjs::hide("separator")
        shinyjs::hide("GeneFile")
        shinyjs::hide("CellsFile")
        shinyjs::hide("MatrixFile")
        shinyjs::hide("InputFile")
        shinyjs::hide("Metadata")
        shinyjs::hide("MinGeneByCells")
        shinyjs::hide("MinCells")
        shinyjs::hide("download_barplot_before_after")
        tmp <- list()
        metadata <- list()
        matrix <- list()
        feature <- list()
        barcode <- list()
        
        ##Loading data and create seurat object
        rdsbool <- FALSE
        withProgress(message = "Reading file",value=0,{
            if(input$file == "RDS"){
                rdsbool <- TRUE
                seuratObj <- readRDS(input$InputFile$datapath)
            }
            if(input$file == "CSV" ){
                if(input$header=="Yes"){
                    headFile <- TRUE
                }else{
                    headFile <- FALSE
                }
                if(input$separator=="tab"){
                    sep <- "\t"
                }else{
                    sep <- input$separator
                }
                if(is.null(input$Metadata)){
                    for(i in 1:length(input$InputFile[,1])){
                        tmp[[i]] <- read.csv(input$InputFile[[i,'datapath']],sep = sep, header = headFile)
                        incProgress(1/(i+1), message = "Creating Seurat object")
                        tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = input$MinGeneByCells, min.features = input$MinCells)
                        tmp[[i]][["file"]] <- input$InputFile[[i,'name']]
                    }
                    seuratObj <- tmp[[1]]
                    if(length(tmp) > 1){
                        for(i in 2:length(input$InputFile[,1])){
                            seuratObj <- merge(seuratObj,tmp[[i]])
                        } 
                    }
                    seuratObj[["Project"]] <- "SeuratProject"
                }else{
                    if(length(input$InputFile$datapath) != length(input$Metadata$datapath)){
                        shinyalert("Oops", "The number of metadata files has to be the same as the number of input files", type = "error")
                    }else{
                        for(i in 1:length(input$InputFile[,1])){
                            tmp[[i]] <- read.csv(input$InputFile[[i,'datapath']],sep = sep, header = headFile)
                            metadata[[i]] <- read.csv(input$Metadata[[i,'datapath']])
                            incProgress(1/(i+1), message = "Creating Seurat object")
                            tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = input$MinGeneByCells, min.features = input$MinCells, meta.data = metadata[[i]])
                            tmp[[i]][["file"]] <- input$InputFile[[i,'name']]
                        }
                        seuratObj <- tmp[[1]]
                        if(length(tmp) > 1){
                            for(i in 2:length(input$InputFile[,1])){
                                seuratObj <- merge(seuratObj,tmp[[i]])
                            } 
                        }
                        seuratObj[["Project"]] <- "SeuratProject"
                    }
                }
            }
            if(input$file == "Mtx" ){
                print(input$MatrixFile)
                if(length(input$MatrixFile$datapath) != length(input$GeneFile$datapath) || length(input$MatrixFile$datapath) != length(input$CellsFile$datapath)){
                    shinyalert("Oops", "The number of matrix has to be the same as the number of cells files and feature files", type = "error")
                }else{
                    for(i in 1:length(input$MatrixFile[,1])){
                        tmp[[i]] <- readMM(input$MatrixFile[[i,'datapath']])
                        barcode.names <- read.delim(input$CellsFile[[i,'datapath']], header = FALSE, stringsAsFactors = FALSE)
                        feature.names <- read.delim(input$GeneFile[[i,'datapath']], header = FALSE, stringsAsFactors = FALSE)
                        colnames(tmp[[i]]) <- barcode.names$V1
                        rownames(tmp[[i]]) <- feature.names$V1
                        incProgress(1/(i+1), message = "Creating Seurat object")
                        tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = input$MinGeneByCells, min.features = input$MinCells)
                        tmp[[i]][["file"]] <- input$MatrixFile[[i,'name']]
                    } #Can't read tmp file with function ReadMtx 
                    seuratObj <- tmp[[1]]
                    if(length(tmp) > 1){
                        for(i in 2:length(input$MatrixFile[,1])){
                            seuratObj <- merge(seuratObj,tmp[[i]])
                        }   
                    }
                    seuratObj[["Project"]] <- "SeuratProject"
                }
                
            }
            if(input$file == "H5"){
                if(is.null(input$Metadata)){
                    for(i in 1:length(input$InputFile[,1])){
                        tmp[[i]] <- Read10X_h5(input$InputFile[[i,'datapath']])
                        incProgress(1/(i+1), message = "Creating Seurat object")
                        tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = input$MinGeneByCells, min.features = input$MinCells)
                        tmp[[i]][["file"]] <- input$InputFile[[i,'name']]
                    }
                    seuratObj <- tmp[[1]]
                    if(length(tmp) > 1){
                        for(i in 2:length(input$InputFile[,1])){
                            seuratObj <- merge(seuratObj,tmp[[i]])
                        }
                    }
                    seuratObj[["Project"]] <- "SeuratProject"
                }else{
                    if(length(input$InputFile$datapath) != length(input$Metadata$datapath)){
                        shinyalert("Oops", "The number of metadata files has to be the same as the number of input files", type = "error")
                    }else{
                        for(i in 1:length(input$InputFile[,1])){
                            tmp[[i]] <- read.csv(input$InputFile[[i,'datapath']])
                            metadata[[i]] <- read.csv(input$Metadata[[i,'datapath']])
                            incProgress(1/(i+1), message = "Creating Seurat object")
                            tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = input$MinGeneByCells, min.features = input$MinCells, meta.data = metadata[[i]])
                            tmp[[i]][["file"]] <- input$InputFile[[i,'name']]
                        }
                        seuratObj <- tmp[[1]]
                        if(length(tmp) > 1){
                            for(i in 2:length(input$InputFile[,1])){
                                seuratObj <- merge(seuratObj,tmp[[i]])
                            } 
                        }
                        seuratObj[["Project"]] <- "SeuratProject"
                    }
                }
            }
            if(input$file == "Txt"){
                if(input$header=="Yes"){
                    headFile <- TRUE
                }else{
                    headFile <- FALSE
                }
                if(input$separator=="tab"){
                    sep <- "\t"
                }else{
                    sep <- input$separator
                }
                if(is.null(input$Metadata)){
                    for(i in 1:length(input$InputFile[,1])){
                        tmp[[i]] <- read.table(input$InputFile[[i,'datapath']], sep = sep, header = headFile)
                        incProgress(1/(i+1), message = "Creating Seurat object")
                        tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = input$MinGeneByCells, min.features = input$MinCells)
                        tmp[[i]][["file"]] <- input$InputFile[[i,'name']]
                    }
                    seuratObj <- tmp[[1]]
                    if(length(tmp) > 1){
                        for(i in 2:length(input$InputFile[,1])){
                            seuratObj <- merge(seuratObj,tmp[[i]])
                        }   
                    }
                    seuratObj[["Project"]] <- "SeuratProject"
                }else{
                    if(length(input$InputFile$datapath) != length(input$Metadata$datapath)){
                        shinyalert("Oops", "The number of metadata files has to be the same as the number of input files", type = "error")
                    }else{
                        for(i in 1:length(input$InputFile[,1])){
                            tmp[[i]] <- read.csv(input$InputFile[[i,'datapath']], sep = sep, header = headFile)
                            metadata[[i]] <- read.csv(input$Metadata[[i,'datapath']])
                            incProgress(1/(i+1), message = "Creating Seurat object")
                            tmp[[i]] <- CreateSeuratObject(tmp[[i]], min.cells = input$MinGeneByCells, min.features = input$MinCells, meta.data = metadata[[i]])
                            tmp[[i]][["file"]] <- input$InputFile[[i,'name']]
                        }
                        seuratObj <- tmp[[1]]
                        if(length(tmp) > 1){
                            for(i in 2:length(input$InputFile[,1])){
                                seuratObj <- merge(seuratObj,tmp[[i]])
                            } 
                        }
                        seuratObj[["Project"]] <- "SeuratProject"
                    }
                }
            }
            seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^mt-")
             output$PreProcessPlot <- renderPlot({
                makeVlnGreatAgain(seuratObj, grouping = "orig.ident",var = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
            })
            output$downloadVln<- downloadHandler(
              filename="VlnPlot_QC.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(makeVlnGreatAgain(seuratObj,  grouping = "orig.ident",var = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
                dev.off()
              })
            output$BeforeFiltering <- renderPlot({
              MakeComplicatedVln(seuratObj,var = "orig.ident", "QC before subsetting datasets")
            })
            output$downloadBeforeFiltering<- downloadHandler(
              filename="VlnPlot_before_filt.png",
              content=function(file){
                png(file, width = 14 , height = 7)
                print(MakeComplicatedVln(seuratObj,var = "orig.ident", "QC before subsetting datasets"))
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
                showTab(inputId = "tabs", target = "Add annotation")
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
                shinyjs::show("downloadSeurat")
                shinyjs::show("downloadLogFile")
                shinyjs::show("NameSeuratLog")
                shinyjs::show("download_barplot_before_after")
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
                        file.copy("logFile.Rmd", tempReport,overwrite = T)
                        line <- paste0("\tSeuratObjsubset <- subset(seuratObj, subset = nFeature_RNA > ",input$minGenePerCells," & nFeature_RNA < ",input$maxGenePerCells," & nCount_RNA < ",input$maxCountPerCells," & percent.mt < ",input$percentMito,")\n")
                        write(line,tempReport, append =TRUE)
                        if(input$normalization == "SCTransform"){
                            line <- paste0("\tSeuratObjsubset <- SCTransform(SeuratObjsubset)\n")
                            write(line,tempReport, append =TRUE)
                        }else{
                            line <- paste0("\tSeuratObjsubset <- NormalizeData(SeuratObjsubset)\n")
                            write(line,tempReport, append =TRUE)
                            line <- paste0("\tSeuratObjsubset <- FindVariableFeatures(SeuratObjsubset)\n")
                            write(line,tempReport, append =TRUE)
                            line <- paste0("\tSeuratObjsubset <- ScaleData(SeuratObjsubset)\n")
                            write(line,tempReport, append =TRUE)
                        }
                        write("\tSeuratObjsubset <- RunPCA(SeuratObjsubset)\n",tempReport, append =TRUE)
                        write("\tSeuratObjsubset <- FindNeighbors(SeuratObjsubset)\n", tempReport, append = TRUE)
                        write("\tSeuratObjsubset <- RunUMAP(SeuratObjsubset, dims = 1:30)\n", tempReport, append = TRUE)
                        for(i in seq(input$Resolution[1],input$Resolution[2], input$ResolutionStep)){
                            write(paste0("\tSeuratObjsubset <- FindClusters(SeuratObjsubset, res =", i,")\n"),tempReport, append = TRUE)
                        }
                        rmarkdown::render(tempReport,output_file = file, params = command, envir = new.env(parent = globalenv()))
                    }
                )
                #)
                shinyjs::enable("runFiltNorm")
                showTab(inputId = "tabs", target = "Cluster tree")
                showTab(inputId = "tabs", target = "DE between clusters")
                showTab(inputId = "tabs", target = "Data mining for one gene")
                showTab(inputId = "tabs", target = "Data mining for a combination of gene")
                showTab(inputId = "tabs", target = "Extract Information")
                showTab(inputId = "tabs", target ="Subclustering")
                showTab(inputId = "tabs", target = "Add annotation")
                
                ###########################################################################################################################################
                ##########################                                 DATA MINING                              #######################################
                ###########################################################################################################################################
                available_resolution <- grep(colnames(SeuratObjsubset@meta.data), pattern = "*_snn_res.*", value = TRUE)
                updateSelectizeInput(session, "resolution_TreePage", choices=available_resolution, server=TRUE)
                s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm7","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7","Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Polr1b","Nasp","Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b","Mrpl36","E2f8")
                g2m.genes <- c("Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Pimreg","Smc4","Ccnb2","Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
                tryCatch(
                    expr = SeuratObjsubset <- CellCycleScoring(SeuratObjsubset, s.features = s.genes, g2m.features = g2m.genes),
                    error = function(c){
                        shinyalert("Warning", "Cannot calculate phase because the following feature lists do not have enough features present in the object: S.Score exiting...The following feature lists do not have enough features present in the object: G2M.Score exiting...", type = "warning")
                    }
                )
                output$AfterFiltering <- renderPlot({
                  MakeComplicatedVln(SeuratObjsubset,var = "orig.ident", "QC after subsetting datasets")
                })
                output$downloadAfterFiltering<- downloadHandler(
                  filename="VlnPlot_after_filt.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(MakeComplicatedVln(SeuratObjsubset,var = "orig.ident", "QC after subsetting datasets"))
                    dev.off()
                  })
                output$clusteringPlot <- renderPlot({
                    vizu_UMAP(SeuratObjsubset, var = NULL, dlabel = input$showlabelcluster)
                })
                output$download_clustering_plot<- downloadHandler(
                  filename="UMAP_orig_ident.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(SeuratObjsubset, var=NULL, dlabel = input$showlabelcluster))
                    dev.off()
                  })
                
                output$barplotBeforeAfter <- renderPlot({
                    Plotnbcellsbeforeafter(objBefore = seuratObj,objAfter = SeuratObjsubset)
                })
                output$download_barplot_before_after<- downloadHandler(
                  filename="plot_before_after.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(Plotnbcellsbeforeafter(objBefore = seuratObj ,objAfter = SeuratObjsubset))
                    dev.off()
                  })
                output$cellcycle <- renderPlot({
                    validate(
                        need(try(SeuratObjsubset[["Phase"]]!=""),message = "No phase available in this seurat object")
                    )
                    vizu_UMAP(SeuratObjsubset, var = "Phase", dlabel = input$showlabelcc)
                })
                output$download_cell_cycle<- downloadHandler(
                  filename="UMAP_cell_cycle.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(SeuratObjsubset, var="Phase", dlabel = input$showlabelcc))
                    dev.off()
                  })
                output$nbcells_by_datasets <- renderPlot({
                  nbCellsbydt(SeuratObjsubset)
                })
                output$download_nb_cells_dt<- downloadHandler(
                  filename="Plot_nb_cells_dt.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(nbCellsbydt(SeuratObjsubset))
                    dev.off()
                  })
                output$geneExpByCell <- renderPlot({
                  UmapNbGenesCell(SeuratObjsubset)
                })
                output$download_nb_cells_dt<- downloadHandler(
                  filename="UMAP_nb_genes_by_cells.png",
                  content=function(file){
                    png(file,width = 1500 , height = 1100,res = 150)
                    print(UmapNbGenesCell(SeuratObjsubset))
                    dev.off()
                  })
                output$BatchVizu <- renderPlot({
                  vizu_UMAP(SeuratObjsubset, var = "file", dlabel = input$showlabelBatch)
                })
                output$downloadbatchVizu<- downloadHandler(
                  filename="UMAP_batch_visu.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(SeuratObjsubset, var = "file", dlabel = input$showlabelBatch))
                    dev.off()
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
                    vizu_UMAP(SeuratObjsubset, var = input$resolution_TreePage, dlabel = input$addlabels_tree)
                })
                output$downloadUMAP_Tree_page<- downloadHandler(
                  filename="UMAP_resolution_tree_page.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(SeuratObjsubset, var = input$resolution_TreePage,dlabel = input$addlabels_tree))
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
                        vizu_UMAP(SeuratObjsubset,var = input$resolutionDE, dlabel = input$addlabels_DE)
                    })
                    output$downloadUMAP_DE_page<- downloadHandler(
                      filename="UMAP_resolution_DE_page.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(SeuratObjsubset, var = input$resolution_DE,dlabel = input$addlabels_DE))
                        dev.off()
                      })
                    
                })
                
                observeEvent(input$submit, {
                    Idents(SeuratObjsubset) <- reac$test
                    withProgress(message = "Differential expression", value = 0,{
                        if(input$onlyPos == "Yes"){
                            posmarkers <- "TRUE"
                        }
                        if(input$onlyPos =="No"){
                            posmarkers  <- "FALSE"
                        }
                        if(input$Cluster2 == "All"){
                            clusterComp <- ""
                        }else{
                            clusterComp <- input$Cluster2
                        }
                        cluster <- input$Cluster1
                        minimum.percent <- input$percent
                        logthr <- input$LogThreshold
                        makeMarkersList <- function(){
                            validate(
                                need(input$submit,"Here you can get the markers table for one cluster versus another or all. You can adjust the minimum of average of log FC to consider that one genes is a marker of one cluster. you can choose to look for only the positive markers or to keep also the negative markers")
                            )
                            if(clusterComp != ""){
                                FindMarkers(SeuratObjsubset, ident.1 = cluster, ident.2 = clusterComp ,only.pos = posmarkers ,min.pct = minimum.percent, logfc.threshold = logthr)
                            }else{
                                FindMarkers(SeuratObjsubset, ident.1 = cluster ,only.pos = posmarkers ,min.pct = minimum.percent, logfc.threshold = logthr)
                            }
                        }
                        output$Markers <- renderDataTable({
                            mList <- makeMarkersList()
                        },server = FALSE ,extensions = c('Buttons'),
                        options = list(
                            autowidth =TRUE,
                            dom = 'frtipB',
                            buttons = list(list(extend ='csv',filename='TableMarkers'))
                        ),
                        filter = list(position = "top", clear = TRUE))
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
                      shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than PNG", "warning")
                  }
                
                  if(input$ResOrAnnot == "Resolution"){
                    Idents(SeuratObjsubset) <- input$ResolutionDataMining
                    ident_obj <- input$ResolutionDataMining
                  }else{
                    Idents(SeuratObjsubset) <- input$AnnotationDataMining
                    ident_obj <- input$AnnotationDataMining
                  }
                  ### Render FeaturePlot + saving them
                  lapply(seq(10),function(x) # allow to create the 10 graphs
                  {
                    output[[paste0('FeaturePlotMultGenes',x)]] = renderPlot({FP(SeuratObjsubset,input$GeneList[x])})
                  })
                  output$downloadfeatureplot<- downloadHandler(
                    filename="featureplot.png",
                    content=function(file){
                      png(file, width = 1500 , height = 1100,res = 150)
                      print(FeaturePlot(SeuratObjsubset, features = input$GeneList))
                      dev.off()
                    }
                  )
                    
                    #Heatmap Generator 
                    Heatmap <- function(){
                        validate(
                            need(input$GeneList != "","Choose at least one gene and it will show you the heatmap of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
                        )
                        DoHeatmap(SeuratObjsubset,features=input$GeneList, slot= "data")+scale_fill_gradient(low = 'ivory2', high= 'navy')
                    }
                    output$Heatmap <- renderPlot({
                        Heatmap()
                    })
                    if(input$format == "PNG"){
                        output$HMdownload <- downloadHandler(
                            filename = "heatmap.png",
                            content = function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(Heatmap())
                                dev.off()
                            }
                        )
                    }else{
                        output$HMdownload <- downloadHandler(
                            filename = "heatmap.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(Heatmap())
                                dev.off()
                            }
                        )
                    }
                    
                    #Dotplot generator
                    Dotplot <- function(){
                        validate(
                            need(input$GeneList != "","Choose at least one gene and it will show you the dotplot of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
                        )
                        DotPlot(SeuratObjsubset, features=input$GeneList)
                    }
                    output$Dotplot <- renderPlot({
                        Dotplot()
                    })
                    if(input$format == "PNG"){
                        output$DPdownload <- downloadHandler(
                            filename = "dotplot.png",
                            content = function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(Dotplot())
                                dev.off()
                            }
                        )
                    }else{
                        output$DPdownload <- downloadHandler(
                            filename = "dotplot.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(Dotplot())
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
                        vizu_UMAP(SeuratObjsubset,var= NULL,  dlabel = input$addlabels_ODM)
                    })
                    output$downloadUMAP_Cluster_resolution<- downloadHandler(
                      filename="UMAP_resolution_data_mining_page.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(SeuratObjsubset,var = NULL, dlabel = input$addlabels_ODM))
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
                      output[[paste0('ViolinPlot',x)]] = renderPlot({vnPlot(SeuratObjsubset,input$GeneList[x])})
                    }else{
                      Idents(SeuratObjsubset) <- input$AnnotAgainst
                      updateSelectizeInput(session, "labelsToKeep2", choices = levels(SeuratObjsubset), server = T)
                      output[[paste0('ViolinPlot',x)]] = renderPlot({
                        gene_graph <- violinPlot(SeuratObjsubset,input$GeneList[x],var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2=input$labelsToKeep2)
                        Violins[[x]] <- gene_graph
                        plot(gene_graph)
                      })
                    }
                  })
                  observeEvent(input$format, {
                    if(input$SplitVln == F && input$format== "PNG"){
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.png",
                        content = function(file){
                          png(file, width = 1500 , height = 1100,res = 150)
                          print(VlnPlot(SeuratObjsubset, features = input$GeneList, pt.size = 0))
                          dev.off()
                        }
                      )
                    }else if (input$SplitVln == F && input$format == "SVG"){
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.svg",
                        content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(VlnPlot(SeuratObjsubset, features = input$GeneList, pt.size = 0))
                          dev.off()
                        }
                      ) 
                    }else if(input$SplitVln == T && input$format =="PNG"){
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.png",
                        content = function(file){
                          png(file, width = 1500 , height = 1100,res = 150)
                          print(ViolinForDl(SeuratObjsubset, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2))
                          dev.off()
                        }
                      )
                    }else{
                      output$VPdownload <- downloadHandler(
                        filename = "violinplot.svg",
                        content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(ViolinForDl(SeuratObjsubset, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2))
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
                listen_page_mining_mult_genes <- reactive({list(input$ResOrAnnotMult,input$annotationwatch,input$clusterwatch,input$SumorMeans,input$GeneListPool,input$format2)})
                #reactive graph from sum or mean button 
                observeEvent(listen_page_mining_mult_genes(),{
                    if(input$format2 =="SVG"){
                        shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than PNG", "warning")
                    }
                    if(input$ResOrAnnotMult == "Resolution"){
                      Idents(SeuratObjsubset) <- input$clusterwatch 
                    }else{
                      Idents(SeuratObjsubset) <- input$annotationwatch
                    }
                    coord <- Embeddings(SeuratObjsubset[["umap"]])[,1:2]
                    output$dimplotSeurat4page <- renderPlot({
                        vizu_UMAP(SeuratObjsubset, var = NULL, dlabel = input$addlabels_MDM)
                    })
                    output$downloadUMAP_resolution_page_5<- downloadHandler(
                      filename="UMAP_resolution_MDM_page.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(SeuratObjsubset,var = NULL, dlabel = input$addlabels_MDM))
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
          
                    poolGene <- function(){
                        validate(
                            need(input$GeneListPool != "","Choose at least one gene to display projection. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                        )
                        if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                            if (input$SumorMeans == "Sum"){
                                pool_genes <- rowSums(FetchData(SeuratObjsubset, vars =c(input$GeneListPool)))  
                            }else if (input$SumorMeans == "Mean"){
                                pool_genes <- rowMeans(FetchData(SeuratObjsubset, vars =c(input$GeneListPool)))
                            }
                            pool_genes <- cbind.data.frame(coord,counts = pool_genes)
                            ggplot(pool_genes)+geom_point(aes(UMAP_1,UMAP_2, color = counts),shape=20, size=0.25)+ scale_color_gradient(low="lightgrey",high="blue") + theme_classic() 
                        }else{
                            sbt <- AddModuleScore(SeuratObjsubset, features =list(input$GeneListPool), name = "Gene_list")
                            FeaturePlot(sbt, features = "Gene_list1")
                        }
                        
                    }
                    output$FeaturePlotMultGenesPooled <- renderPlot({
                        poolGene()
                    })
                    VlnPlotPooled <- function(){
                        validate(
                            need(input$GeneListPool != "", "Choose at least one gene to display violin plot")
                        )
                        if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                            if(input$SumorMeans == "Sum"){
                                if(input$ResOrAnnotMult == "Resolution"){
                                  pool_genesVln <- cbind(SeuratObjsubset[[input$clusterwatch]],rowSums(FetchData(SeuratObjsubset, vars = c(input$GeneListPool))))
                                }else{
                                  pool_genesVln <- cbind(SeuratObjsubset[[input$annotationwatch]],rowSums(FetchData(SeuratObjsubset, vars = c(input$GeneListPool))))
                                }
                                
                            }else{
                                if(input$ResOrAnnotMult == "Resolution"){
                                  pool_genesVln <- cbind(SeuratObjsubset[[input$clusterwatch]],rowMeans(FetchData(SeuratObjsubset, vars = c(input$GeneListPool))))
                                }else{
                                  pool_genesVln <- cbind(SeuratObjsubset[[input$annotationwatch]],rowMeans(FetchData(SeuratObjsubset, vars = c(input$GeneListPool))))
                                }
                                
                            }
                            colnames(pool_genesVln) <- c("resolution","GeneList")
                            ggplot(pool_genesVln, aes(x=resolution,y=GeneList, fill = resolution))+geom_violin()+scale_fill_hue()
                        }else{
                            sbt <- AddModuleScore(SeuratObjsubset, features =list(input$GeneListPool), name = "Gene_list")
                            VlnPlot(sbt, features = "Gene_list1")
                        }
                    }
                    output$ViolinPlotMultGenesPooled <- renderPlot({
                        VlnPlotPooled()
                    })
                    DotplotPooled <- function(){
                        validate(
                            need(input$GeneListPool != "", "Choose at least one gene to display dot plot")
                        )
                        if(input$SumorMeans == "Sum"){
                            if(input$ResOrAnnotMult == "Resolution"){
                              Idents(SeuratObjsubset) <- input$clusterwatch 
                            }else{
                              Idents(SeuratObjsubset) <- input$annotationwatch
                            }
                            test <- as.data.frame(cbind(rowSums(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$clusterwatch]]))
                            colnames(test)<- c("V1","V2")
                            percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                            colnames(percentExpressed) <- c("cluster", "logical","number")
                            percentExpressed <- filter(percentExpressed, logical == "TRUE")
                            totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                            percent <- percentExpressed$number/totalCellsbyCluster$n
                            percent <- as.data.frame(cbind(cluster = levels(SeuratObjsubset),RatioCellsExpressingGeneList = percent))
                            avgExpression <-test %>% group_by(V2) %>%summarise(avg =sum(V1))
                            percent <- cbind(percent, SumExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                            percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                            ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = SumExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                            
                        }else if(input$SumorMeans=="Mean"){
                            if(input$ResOrAnnotMult == "Resolution"){
                              Idents(SeuratObjsubset) <- input$clusterwatch 
                            }else{
                              Idents(SeuratObjsubset) <- input$annotationwatch
                            }
                            test <- as.data.frame(cbind(rowMeans(FetchData(SeuratObjsubset, vars =c(input$GeneListPool))),SeuratObjsubset[[input$clusterwatch]]))
                            colnames(test)<- c("V1","V2")
                            percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                            colnames(percentExpressed) <- c("cluster", "logical","number")
                            percentExpressed <- filter(percentExpressed, logical == "TRUE")
                            totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                            percent <- percentExpressed$number/totalCellsbyCluster$n
                            percent <- as.data.frame(cbind(cluster = levels(SeuratObjsubset),RatioCellsExpressingGeneList = percent))
                            avgExpression <-test %>% group_by(V2) %>%summarise(avg =mean(V1))
                            percent <- cbind(percent, AvgExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                            percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                            ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = AvgExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                        }else{
                            sbt <- AddModuleScore(SeuratObjsubset, features =list(input$GeneListPool), name = "Gene_list")
                            DotPlot(sbt, features = "Gene_list1")
                        }
                        
                    }
                    output$DotPlotMultGenePooled <- renderPlot({
                        DotplotPooled()
                    })
                    if(input$format == "PNG"){
                        output$VPdownloadMultGenesPooled <- downloadHandler(
                            filename = "ViolinPlotMultGenesPooled.png",
                            content = function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(VlnPlotPooled())
                                dev.off()
                            }
                        )
                    }else{
                        output$VPdownloadMultGenesPooled <- downloadHandler(
                            filename = "ViolinPlotMultGenesPooled.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(VlnPlotPooled())
                                dev.off()
                            }
                        )
                    }
                    
                    
                    if(input$format2 =="PNG"){
                        output$downloadfeatureplotMultGenesPooled<- downloadHandler(
                            filename = "featurePlotMultGenesPooled.png",
                            content = function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(poolGene())
                                dev.off()
                            }
                        )
                    }else{
                        output$downloadfeatureplotMultGenesPooled<- downloadHandler(
                            filename = "featurePlotMultGenesPooled.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(poolGene())
                                dev.off()
                            }
                        )
                    }
                    
                    if(input$format2 == "PNG"){
                        output$dotplotdownloadMultGenesPooled<- downloadHandler(
                            filename = "DotPlotMultGenesPooled.png",
                            content = function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(poolGene())
                                dev.off()
                            }
                        )  
                    }else{
                        output$dotplotdownloadMultGenesPooled<- downloadHandler(
                            filename = "DotPlotMultGenesPooled.svg",
                            content = function(file){
                                svg(file, width = 14 , height = 7)
                                print(poolGene())
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
                  
                  poolGeneFromList <- function(){
                    validate(
                      need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.This plot will show the pool of each gene expression for each gene id contained in the file. \nExample of a file :\nHeader\nMlxipl\nSox9\nEtc...  ")
                    )
                    coord <- Embeddings(SeuratObjsubset[["umap"]])[,1:2]
                    gene_list <- Readfile()
                    if(input$ResOrAnnotMult == "Resolution"){
                      Idents(SeuratObjsubset) <- input$clusterwatch
                    }else{
                      Idents(SeuratObjsubset) <- input$annotationwatch
                    }
                    if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                      if(input$SumorMeans =="Sum"){
                        gene_list_pooled <- rowSums(FetchData(SeuratObjsubset, vars = gene_list))
                      }
                      if(input$SumorMeans =="Mean"){
                        gene_list_pooled <- rowMeans(FetchData(SeuratObjsubset, vars = gene_list))
                      }
                      gene_list_pooled <- cbind.data.frame(coord, counts = gene_list_pooled)
                      ggplot(gene_list_pooled)+geom_point(aes(UMAP_1,UMAP_2, color = counts),shape=20, size=0.25)+ scale_color_gradient(low="lightgrey",high="blue") + theme_classic()
                    }else{
                      sbt1 <- AddModuleScore(SeuratObjsubset, features =gene_list, name = "Gene_list")
                      FeaturePlot(sbt1, features = "Gene_list1")
                    }
                  }
                  
                  output$FeaturePlotMultGenesPooledFilesInput <- renderPlot({
                    poolGeneFromList()
                  })
                  poolGenefromListVP <- function(){
                      validate(
                          need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.")
                      )
                      gene_list <- Readfile()
                      if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                          if(input$SumorMeans =="Sum"){
                              if(input$ResOrAnnotMult == "Resolution"){
                                pool_genesVln <- cbind(SeuratObjsubset[[input$clusterwatch]],rowSums(FetchData(SeuratObjsubset, vars = gene_list)))
                              }else{
                                pool_genesVln <- cbind(SeuratObjsubset[[input$annotationwatch]],rowSums(FetchData(SeuratObjsubset, vars = gene_list)))
                              }
                              
                          }
                          if(input$SumorMeans =="Mean"){
                              if(input$ResOrAnnotMult == "Resolution"){
                                pool_genesVln <- cbind(SeuratObjsubset[[input$clusterwatch]],rowMeans(FetchData(SeuratObjsubset, vars = gene_list)))
                              }else{
                                pool_genesVln <- cbind(SeuratObjsubset[[input$annotationwatch]],rowMeans(FetchData(SeuratObjsubset, vars = gene_list)))
                              }
                              
                          }
                          colnames(pool_genesVln) <- c("resolution","GeneList")
                          ggplot(pool_genesVln, aes(x=resolution,y=GeneList, fill = resolution))+geom_violin()+scale_fill_hue()
                      }else{
                          if(input$ResOrAnnotMult == "Resolution"){
                            Idents(SeuratObjsubset) <- input$clusterwatch 
                          }else{
                            Idents(SeuratObjsubset) <- input$annotationwatch
                          }
                          sbt1 <- AddModuleScore(SeuratObjsubset, features =list(gene_list), name = "Gene_list")
                          VlnPlot(sbt1, features = "Gene_list1")
                      }
                      
                  }
                  output$ViolinPlotMultGenesPooledFilesInput <- renderPlot({
                      poolGenefromListVP()
                  })
                    
                  DotplotPooledFromFile <- function(){
                      validate(
                          need(input$fileGeneList, "Choose a file that contain a list of genes to display the projection.")
                      )
                      gene_list <- Readfile()
                      if(input$SumorMeans == "Sum"){
                          if(input$ResOrAnnotMult == "Resolution"){
                          Idents(SeuratObjsubset) <- input$clusterwatch 
                          }else{
                            Idents(SeuratObjsubset) <- input$annotationwatch
                          }
                          test <- as.data.frame(cbind(rowSums(FetchData(SeuratObjsubset, vars =gene_list)),SeuratObjsubset[[input$clusterwatch]]))
                          colnames(test)<- c("V1","V2")
                          percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                          colnames(percentExpressed) <- c("cluster", "logical","number")
                          percentExpressed <- filter(percentExpressed, logical == "TRUE")
                          totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                          percent <- percentExpressed$number/totalCellsbyCluster$n
                          percent <- as.data.frame(cbind(cluster = levels(SeuratObjsubset),RatioCellsExpressingGeneList = percent))
                          avgExpression <-test %>% group_by(V2) %>%summarise(avg =sum(V1))
                          percent <- cbind(percent, SumExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                          percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                          ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = SumExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                          
                      }else if (input$SumorMeans == "Mean"){
                        if(input$ResOrAnnotMult == "Resolution"){
                          Idents(SeuratObjsubset) <- input$clusterwatch 
                          }else{
                            Idents(SeuratObjsubset) <- input$annotationwatch
                          }
                          test <- as.data.frame(cbind(rowMeans(FetchData(SeuratObjsubset, vars =gene_list)),SeuratObjsubset[[input$clusterwatch]]))
                          colnames(test)<- c("V1","V2")
                          percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                          colnames(percentExpressed) <- c("cluster", "logical","number")
                          percentExpressed <- filter(percentExpressed, logical == "TRUE")
                          totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                          percent <- percentExpressed$number/totalCellsbyCluster$n
                          percent <- as.data.frame(cbind(cluster = levels(SeuratObjsubset),RatioCellsExpressingGeneList = percent))
                          avgExpression <-test %>% group_by(V2) %>%summarise(avg =mean(V1))
                          percent <- cbind(percent, AvgExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                          percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                          ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = AvgExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                      }else{
                          if(input$ResOrAnnotMult == "Resolution"){
                            Idents(SeuratObjsubset) <- input$clusterwatch 
                          }else{
                            Idents(SeuratObjsubset) <- input$annotationwatch
                          }
                          sbt1 <- AddModuleScore(SeuratObjsubset, features =list(gene_list), name = "Gene_list")
                          DotPlot(sbt1, features = "Gene_list1")
                      }
                      
                  }
                  output$DotPlotMultGenePooledFromFile <- renderPlot({
                      DotplotPooledFromFile()
                  })
                  if(input$format2 == "PNG"){
                      output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                          filename = "ViolinPlotMultGenesPooled.png",
                          content = function(file){
                              png(file, width = 1500 , height = 1100,res = 150)
                              print(poolGenefromListVP())
                              dev.off()
                          }
                      )
                  }else{
                      output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                          filename = "ViolinPlotMultGenesPooled.svg",
                          content = function(file){
                              svg(file, width = 14 , height = 7)
                              print(poolGenefromListVP())
                              dev.off()
                          }
                      )
                  }
                  
                  if(input$format2 =="PNG"){
                      output$downloadfeatureplotMultGenesPooledFilesInput <- downloadHandler(
                          filename = "featurePlotFilesInput.png",
                          content = function(file){
                              png(file, width = 1500 , height = 1100,res = 150)
                              print(poolGeneFromList())
                              dev.off()
                          }
                      ) 
                  }else{
                      output$downloadfeatureplotMultGenesPooledFilesInput <- downloadHandler(
                          filename = "featurePlotFilesInput.svg",
                          content = function(file){
                              svg(file, width = 14 , height = 7)
                              print(poolGeneFromList())
                              dev.off()
                          }
                      ) 
                  }
                  
                  if(input$format2 == "PNG"){
                      output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                          filename = "DotPlotMultGenesPooled.png",
                          content = function(file){
                              png(file, width = 1500 , height = 1100,res = 150)
                              print(DotplotPooledFromFile())
                              dev.off()
                          }
                      ) 
                  }else{
                      output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                          filename = "DotPlotMultGenesPooled.svg",
                          content = function(file){
                              svg(file, width = 14 , height = 7)
                              print(DotplotPooledFromFile())
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
                        vizu_UMAP(SeuratObjsubset,var = input$clusterNumber, dlabel = input$addlabels_info)
                    })
                    output$downloadUMAP_info<- downloadHandler(
                      filename = "UMAP_info.png",
                      content = function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(SeuratObjsubset,var = input$clusterNumber, dlabel = input$addlabels_info))
                        dev.off()
                      }
                    ) 
                    toListen <- reactive({list(input$clusterInfo,input$chooseVar2Plot, input$stackedBP, input$FreqOrVal)})
                    observeEvent(toListen(), {
                        Plot2Render <- function(){
                            if(input$clusterInfo == "All"){
                                if(input$stackedBP == "No"){
                                    table_tmp2 <-as.data.frame(table(c(SeuratObjsubset[[input$chooseVar2Plot]], SeuratObjsubset[[input$clusterNumber]])))
                                    colnames(table_tmp2) <- c("condition","cluster", "freq")
                                    test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
                                    g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = position_dodge())+scale_fill_manual(values = test)
                                }else{
                                    if(input$FreqOrVal =="Frequence"){
                                        table_tmp2 <-as.data.frame(table(c(SeuratObjsubset[[input$chooseVar2Plot]], SeuratObjsubset[[input$clusterNumber]])))
                                        colnames(table_tmp2) <- c("condition","cluster", "freq")
                                        test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
                                        g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "fill")+scale_fill_manual(values = test)
                                    }else{
                                        table_tmp2 <-as.data.frame(table(c(SeuratObjsubset[[input$chooseVar2Plot]], SeuratObjsubset[[input$clusterNumber]])))
                                        colnames(table_tmp2) <- c("condition","cluster", "freq")
                                        test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
                                        g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = test)
                                    }
                                }
                            }else{
                                table_tmp <- as.data.frame(table(c(SeuratObjsubset[[input$chooseVar2Plot]],SeuratObjsubset[[input$clusterNumber]])))
                                colnames(table_tmp) <- c("condition","cluster", "freq")
                                table_tmp <- subset(table_tmp, table_tmp$cluster == input$clusterInfo) # allow to fix the color with the good number of variable ^
                                test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp$condition)))
                                g <- ggplot(table_tmp, aes(x = condition,freq, fill = condition ))+geom_bar(stat = "identity")+scale_fill_manual(values = test)
                            }
                        }
                        percentRender <- function(){
                            table_tmp2 <-as.data.frame(table(c(SeuratObjsubset[[input$chooseVar2Plot]], SeuratObjsubset[[input$clusterNumber]])))
                            colnames(table_tmp2) <- c("condition","cluster", "freq")
                            tmp <- as.data.frame(rep(table(SeuratObjsubset[[input$chooseVar2Plot]]),length(unique(SeuratObjsubset[[input$clusterNumber]]))))
                            colnames(tmp)<- "total"
                            table_tmp2 <- cbind(table_tmp2,tmp)
                            table_tmp2$percent <- table_tmp2$freq/table_tmp2$total*100
                            test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$cluster)))
                            g <- ggplot(table_tmp2, aes(x= condition, percent, fill = cluster))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = test)
                            ggplotly(g+theme(axis.text.x = element_text(angle =60,hjust=1)))
                        }
                        output$ggplot_information <- renderPlotly({
                            p <- Plot2Render()
                            ggplotly(p +theme(legend.position = "none",axis.text.x = element_text(angle =60,hjust=1)))
                        })
                        output$legend_information <- renderPlot({
                            p <- Plot2Render()
                            grid.newpage()
                            legend <- cowplot::get_legend(p)
                            grid.draw(legend)
                        })
                        output$cluster_percent <- renderPlotly({
                            percentRender()
                        })
                        output$downloadInformationPlot2 <- downloadHandler(
                          filename = "Percent_Plot.png",
                          content = function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(percentRender())
                            dev.off()
                          }
                        )
                        Table2Render <- function(){
                            if(input$clusterInfo == "All"){
                                as.data.frame(table(c(SeuratObjsubset[[input$chooseVar2Plot]], SeuratObjsubset[[input$clusterNumber]])))
                            }else{
                                as.data.frame(table(SeuratObjsubset[[input$chooseVar2Plot]][which(SeuratObjsubset[[input$clusterNumber]]==input$clusterInfo),]))
                            }
                        }
                        output$get_info <- DT::renderDataTable({
                            Table2Render()
                        })
                        output$downloadInformationPlot <- downloadHandler(
                            filename = "Information_Plot.png",
                            content = function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(Plot2Render())
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
                        vizu_UMAP(SeuratObjsubset,var = input$res4Annot, dlabel = input$labelsAnnot)
                    })
                    output$downloadUMAP_resolution_annot_page<- downloadHandler(
                      filename="UMAP_resolution_annot_page.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(SeuratObjsubset,var = input$res4Annot, dlabel = input$labelsAnnot))
                        dev.off()
                      })
                })
                output$condPlot <- renderPlot({
                    vizu_UMAP(SeuratObjsubset, var = input$Annot2Complete, dlabel = input$labelsconditional)
                })
                output$downloadUMAP_Conditional<- downloadHandler(
                  filename="UMAP_conditional_page.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(SeuratObjsubset,var = input$Annot2Complete, dlabel = input$labelsconditional))
                    dev.off()
                  })
                observeEvent(input$submitAnnot,{
                    if(input$choicesCreateOrComplete == "Create"){
                        if(input$AnnotName == "" | input$clusterName == "" | length(input$cluster4Annot) == 0){
                            shinyalert("Oops", "The name of the column and the name of the annotation you want to give has to be filled", type = "error")
                        }else{
                            if(is.element(input$clusterName, names(SeuratObjsubset@meta.data))){
                                f.clusters <- SeuratObjsubset[[input$clusterName]]
                                f.clusters[SeuratObjsubset[[resChoice]] == input$cluster4Annot[1]] <- input$AnnotName
                                if(length(input$cluster4Annot) > 1){
                                    for(i in 2:length(input$cluster4Annot)){
                                        f.clusters[SeuratObjsubset[[resChoice]] == input$cluster4Annot[i]] <- input$AnnotName 
                                    }
                                }
                                SeuratObjsubset[[input$clusterName]] <<- f.clusters
                            }else{
                                f.clusters <- rep("unknown", dim(SeuratObjsubset)[2])
                                names(f.clusters) <- colnames(SeuratObjsubset)
                                f.clusters[SeuratObjsubset[[resChoice]] == input$cluster4Annot[1]] <- input$AnnotName
                                if(length(input$cluster4Annot) > 1){
                                    for(i in 2:length(input$cluster4Annot)){
                                        f.clusters[SeuratObjsubset[[resChoice]] == input$cluster4Annot[i]] <- input$AnnotName 
                                    }
                                }
                                SeuratObjsubset[[input$clusterName]] <<- f.clusters
                            }
                            shinyalert("Done","You have correctly annotated your cluster", type = "success")
                            output$postAnnotation <- renderPlot({
                                vizu_UMAP(SeuratObjsubset, var = input$clusterName, dlabel = input$labelsPostAnnot)
                            })
                            output$downloadUMAP_post_annotation<- downloadHandler(
                              filename="UMAP_post_annot.png",
                              content=function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(vizu_UMAP(SeuratObjsubset,var = input$clusterName, var = input$labelsPostAnnot))
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
                            vizu_UMAP(SeuratObjsubset, var = input$Annot2Complete, dlabel = ,input$labelsPostAnnot)
                        })
                        output$downloadUMAP_post_annotation<- downloadHandler(
                          filename="UMAP_post_annot.png",
                          content=function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(vizu_UMAP(SeuratObjsubset,var = input$Annot2Complete, input$labelsPostAnnot))
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
                            vizu_UMAP(SeuratObjsubset,var = input$clusterResPage7, dlabel = input$labelsSubset)
                        })
                        output$downloadUMAP_sbt<- downloadHandler(
                          filename="UMAP_sbt.png",
                          content=function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(vizu_UMAP(SeuratObjsubset,var = input$clusterResPage7, input$labelsSubset))
                            dev.off()
                          })
                        observeEvent(input$subclusterize,{
                            test <- NULL
                            
                            if(length(input$subclusterize)==length(clusterSbt)){
                                test <- colnames(SeuratObjsubset)
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.png",
                                  content=function(file){
                                    png(file, width = 1500 , height = 1100,res = 150)
                                    print(MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                    dev.off()
                                  })
                            }else{
                              test <- colnames(SeuratObjsubset)[which(SeuratObjsubset[[]][input$clusterResPage7] == input$subclusterize[1])]
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test ,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.png",
                                  content=function(file){
                                    png(file, width = 1500 , height = 1100,res = 150)
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
                                      filename="UMAP_cell_kept.png",
                                      content=function(file){
                                        png(file, width = 1500 , height = 1100,res = 150)
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
                                          filename="UMAP_cell_kept.png",
                                          content=function(file){
                                            png(file, width = 1500 , height = 1100,res = 150)
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
                            vizu_UMAP(SeuratObjsubset,var = input$whichAnnot, dlabel = input$labelsSubset)
                        })
                        output$downloadUMAP_sbt<- downloadHandler(
                          filename="UMAP_sbt.png",
                          content=function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(vizu_UMAP(SeuratObjsubset,var = input$whichAnnot, input$labelsSubset))
                            dev.off()
                          })
                        observeEvent(input$subannot,{
                            test <- NULL
                            if(length(input$subannot) == length(annotSbt)){
                                test <- colnames(SeuratObjsubset)
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.png",
                                  content=function(file){
                                    png(file, width = 1500 , height = 1100,res = 150)
                                    print(MakeUMAP_Red(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                    dev.off()
                                  })
                            }else{
                              test <- colnames(SeuratObjsubset)[which(SeuratObjsubset[[]][input$whichAnnot] == input$subannot[1])]
                                output$Dimplot_kept <- renderPlot({
                                    MakeUMAPhighlight(SeuratObjsubset, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                                })
                                output$downloadUMAPcellKept<- downloadHandler(
                                  filename="UMAP_cell_kept.png",
                                  content=function(file){
                                    png(file, width = 1500 , height = 1100,res = 150)
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
                                      filename="UMAP_cell_kept.png",
                                      content=function(file){
                                        png(file, width = 1500 , height = 1100,res = 150)
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
                                          filename="UMAP_cell_kept.png",
                                          content=function(file){
                                            png(file, width = 1500 , height = 1100,res = 150)
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
                    output$dimplot_aftersbt <- renderPlot({
                        vizu_UMAP(subsetSeuratObj, dlabel = input$labelsAfterSubset)
                    })
                    output$downloadUMAPaftersbt<- downloadHandler(
                      filename="UMAP_After_sbt.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(subsetSeuratObj, dlabel = input$labelsAfterSubset))
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
                        file.copy("logFile.Rmd", tempReport,overwrite = T)
                        line <- paste0("\tSeuratObjsubset <- subset(seuratObj, subset = nFeature_RNA > ",input$minGenePerCells," & nFeature_RNA < ",input$maxGenePerCells," & percent.mt < ",input$percentMito,")")
                        write(line,tempReport, append =TRUE)
                        if(input$normalization == "SCTransform"){
                            line <- paste0("\tSeuratObjsubset <- SCTransform(SeuratObjsubset)")
                            write(line,tempReport, append =TRUE)
                        }else{
                            line <- paste0("\tSeuratObjsubset <- NormalizeData(SeuratObjsubset)")
                            write(line,tempReport, append =TRUE)
                            line <- paste0("\tSeuratObjsubset <- FindVariableFeatures(SeuratObjsubset)")
                            write(line,tempReport, append =TRUE)
                            line <- paste0("\tSeuratObjsubset <- ScaleData(SeuratObjsubset)")
                            write(line,tempReport, append =TRUE)
                        }
                        write("\tSeuratObjsubset <- RunPCA(SeuratObjsubset)",tempReport, append =TRUE)
                        write("\tSeuratObjsubset <- FindNeighbors(SeuratObjsubset)", tempReport, append = TRUE)
                        write("\tSeuratObjsubset <- RunUMAP(SeuratObjsubset, dims = 1:30)\n", tempReport, append = TRUE)
                        for(i in seq(input$Resolution[1],input$Resolution[2], input$ResolutionStep)){
                            write(paste0("\tSeuratObjsubset <- FindClusters(SeuratObjsubset, res =", i,")"),tempReport, append = TRUE)
                        }
                        
                        if(input$subcluster == "Cluster"){
                            line <- paste0("\tsubsetSeuratObj <- subset(SeuratObjsubset, subset =",input$clusterResPage7," == ",list(input$subclusterize),")")
                        }else{
                            line <- paste0("\tsubsetSeuratObj <- subset(SeuratObjsubset, subset = ",input$whichAnnot," == ",list(input$subannot),")")
                        }
                        write(line,tempReport, append =TRUE)
                        write("\tsubsetSeuratObj <- RunPCA(subsetSeuratObj)",tempReport, append =TRUE)
                        write("\tsubsetSeuratObj <- FindNeighbors(subsetSeuratObj)", tempReport, append = TRUE)
                        write("\tsubsetSeuratObj <- RunUMAP(subsetSeuratObj, dims = 1:30)\n", tempReport, append = TRUE)
                        write("\tfor(i in seq(0,1,0.1)){",tempReport, append = T)
                        write("\t\tsubsetSeuratObj <- findClusters(subsetSeuratObj, res = i)",tempReport, append = T)
                        write("\t}\n",tempReport, append = T)
                        rmarkdown::render(tempReport,output_file = file, params = command, envir = new.env(parent = globalenv()))
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
            
            ## Violin plot on data
            output$featureQC <- renderPlot({
                makeVlnGreatAgain(seuratObj,  grouping = "orig.ident",var = c("nFeature_RNA","nCount_RNA","percent.mt"), col = 3) 
            })
            output$downloadVln<- downloadHandler(
              filename="Violin_QC.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(makeVlnGreatAgain(seuratObj, grouping = "orig.ident", var=c("nFeature_RNA","nCount_RNA","percent.mt"), col = 3))
                dev.off()
              })
            
            ## QC resolution
            output$plotClusterQC <- renderPlot({
                vizu_UMAP(seuratObj,var= NULL,dlabel = input$addlabels_res)
            })
            output$downloadUMAP_resolution_qc<- downloadHandler(
              filename="UMAP_Cluster_QC.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
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
              filename="UMAP_Cell_cycle.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(vizu_UMAP(seuratObj,var="Phase", dlabel = input$showlabelccQC))
                dev.off()
              })
            
            ## Hist nb cell by dataset
            output$nbcells_by_datasetsQC <- renderPlot({
                nbCellsbydt(seuratObj)
            })
            output$downloadHistNbCells<- downloadHandler(
              filename="Hist_nb_cells_dt.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(nbCellsbydt(seuratObj))
                dev.off()
              })
            
            ## Nb gene by cell
            output$geneExpByCellQC <- renderPlot({
                UmapNbGenesCell(seuratObj)
            })
            output$downloadUMAP_gene_by_cell<- downloadHandler(
              filename="UMAP_Nb_gene_by_cell.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(UmapNbGenesCell(seuratObj))
                dev.off()
              })
            
            ## choice plot 
            output$plotChoice <- renderPlot({
                vizu_UMAP(seuratObj,var =input$InfoToPlot,  dlabel = input$addlabels_choice)
            })
            output$downloadUMAP_choice<- downloadHandler(
              filename="UMAP_choice.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(vizu_UMAP(seuratObj,var=input$InfoToPlot,  dlabel = input$addlabels_choice))
                dev.off()
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
                vizu_UMAP(seuratObj, var = input$resolution_TreePage, dlabel = input$addlabels_tree)
            })
            output$downloadUMAP_Tree_page<- downloadHandler(
              filename="UMAP_From_tree_page.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(vizu_UMAP(seuratObj,var=input$resolution_TreePage, dlabel = input$addlabels_tree))
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
                  filename="UMAP_DE_page.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(seuratObj,var=input$resolution_TreePage, dlabel = input$addlabels_DE))
                    dev.off()
                })
                
            })
            output$Markers <- renderDataTable(NULL)
            observeEvent(input$submit, {
                Idents(seuratObj) <- reac$test
                withProgress(message = "Differential expression", value = 0,{
                    if(input$onlyPos == "Yes"){
                        posmarkers <- "TRUE"
                    }
                    if(input$onlyPos =="No"){
                        posmarkers  <- "FALSE"
                    }
                    if(input$Cluster2 == "All"){
                        clusterComp <- ""
                    }else{
                        clusterComp <- input$Cluster2
                    }
                    cluster <- input$Cluster1
                    minimum.percent <- input$percent
                    logthr <- input$LogThreshold
                    makeMarkersList <- function(){
                        validate(
                            need(input$submit,"Here you can get the markers table for one cluster versus another or all. You can adjust the minimum of average of log FC to consider that one genes is a marker of one cluster. you can choose to look for only the positive markers or to keep also the negative markers")
                        )
                        if(clusterComp != ""){
                            FindMarkers(seuratObj, ident.1 = cluster, ident.2 = clusterComp ,only.pos = posmarkers ,min.pct = minimum.percent, logfc.threshold = logthr)
                        }else{
                            FindMarkers(seuratObj, ident.1 = cluster ,only.pos = posmarkers ,min.pct = minimum.percent, logfc.threshold = logthr)
                        }
                    }
                    output$Markers <- renderDataTable({
                        mList <- makeMarkersList()
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
                shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than PNG", "warning")
              }
              if(input$ResOrAnnot == "Resolution"){
                Idents(seuratObj) <- input$ResolutionDataMining
                ident_obj <- input$ResolutionDataMining
                
                
              }else{
                Idents(seuratObj) <- input$AnnotationDataMining
                ident_obj <- input$AnnotationDataMining
                
              }
                
                
                ## Create and save the featurePlot
                lapply(seq(10),function(x) # allow to create the 10 graphd
                {
                  output[[paste0('FeaturePlotMultGenes',x)]] = renderPlot({FP(seuratObj,input$GeneList[x])})
                })
                lapply(seq(length(input$GeneList)), function(x){
                  output$downloadfeatureplot<- downloadHandler(
                    filename="featureplot.png",
                    content=function(file){
                      png(file, width = 1500 , height = 1100,res = 150)
                      print(FeaturePlot(seuratObj, features = input$GeneList))
                      dev.off()
                    })
                })
                
                
                
                
                
                #Heatmap Generator 
                Heatmap <- function(){
                    validate(
                        need(input$GeneList != "","Choose at least one gene and it will show you the heatmap of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
                    )
                    DoHeatmap(seuratObj,features=input$GeneList, slot = "data")+scale_fill_gradient(low = 'ivory2', high= 'navy')
                }
                output$Heatmap <- renderPlot({
                    Heatmap()
                })
                if(input$format =="PNG"){
                    output$HMdownload <- downloadHandler(
                        filename = "heatmap.png",
                        content = function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(Heatmap())
                            dev.off()
                        }
                    )
                }else{
                    output$HMdownload <- downloadHandler(
                        filename = "heatmap.svg",
                        content = function(file){
                            svg(file, width = 14 , height = 7)
                            print(Heatmap())
                            dev.off()
                        }
                    )
                }
                #Dotplot generator
                Dotplot <- function(){
                    validate(
                        need(input$GeneList != "","Choose at least one gene and it will show you the dotplot of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
                    )
                    DotPlot(seuratObj, features=input$GeneList)
                }
                output$Dotplot <- renderPlot({
                    Dotplot()
                })
                if(input$format =="PNG"){
                    output$DPdownload <- downloadHandler(
                        filename = "dotplot.png",
                        content = function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(Dotplot())
                            dev.off()
                        }
                    ) 
                }else{
                    output$DPdownload <- downloadHandler(
                        filename = "dotplot.svg",
                        content = function(file){
                            svg(file, width = 14 , height = 7)
                            print(Dotplot())
                            dev.off()
                        }
                    ) 
                }
                
                output$clusterOutput<- renderPlot({
                    vizu_UMAP(seuratObj, var = ident_obj, dlabel = input$addlabels_ODM)
                })
                output$downloadUMAP_Cluster_resolution<- downloadHandler(
                  filename="UMAP_cluster_or_annot.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(seuratObj,var=ident_obj, dlabel = input$addlabels_ODM))
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
                  output[[paste0('ViolinPlot',x)]] = renderPlot({vnPlot(seuratObj,input$GeneList[x])})
                }else{
                  Idents(seuratObj) <- input$AnnotAgainst
                  updateSelectizeInput(session, "labelsToKeep2", choices = levels(seuratObj), server = T)
                  output[[paste0('ViolinPlot',x)]] = renderPlot({
                    gene_graph <- violinPlot(seuratObj,input$GeneList[x],var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2=input$labelsToKeep2)
                    Violins[[x]] <- gene_graph
                    plot(gene_graph)
                  })
                }
              })
              observeEvent(input$format, {
                if(input$SplitVln == F && input$format== "PNG"){
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.png",
                    content = function(file){
                      png(file, width = 1500 , height = 1100,res = 150)
                      print(VlnPlot(seuratObj, features = input$GeneList, pt.size = 0))
                      dev.off()
                    }
                  )
                }else if (input$SplitVln == F && input$format == "SVG"){
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.svg",
                    content = function(file){
                      svg(file, width = 14 , height = 7)
                      print(VlnPlot(seuratObj, features = input$GeneList, pt.size = 0))
                      dev.off()
                    }
                  ) 
                }else if(input$SplitVln == T && input$format =="PNG"){
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.png",
                    content = function(file){
                      png(file, width = 1500 , height = 1100,res = 150)
                      print(ViolinForDl(seuratObj, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2))
                      dev.off()
                    }
                  )
                }else{
                  output$VPdownload <- downloadHandler(
                    filename = "violinplot.svg",
                    content = function(file){
                      svg(file, width = 14 , height = 7)
                      print(ViolinForDl(seuratObj, geneList = input$GeneList, var1 = ident_obj_vln, var2 = input$AnnotAgainst, conditionVar1 = input$labelsToKeep1, conditionVar2 = input$labelsToKeep2))
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
                  shinyalert("warning", "SVG format is a good choice of graph for your paper because it can be easily formated but is heavier than PNG", "warning")
              }
              if(input$ResOrAnnotMult == "Resolution"){
                Idents(seuratObj) <- input$clusterwatch
                ident_obj_mult <- input$clusterwatch
              }else{
                Idents(seuratObj) <- input$annotationwatch
                ident_obj_mult <- input$annotationwatch
              }
              coord <- Embeddings(seuratObj[["umap"]])[,1:2]
              
              # print and save umap resolution
              output$dimplotSeurat4page <- renderPlot({
                  vizu_UMAP(seuratObj, var = ident_obj_mult, dlabel = input$addlabels_MDM)
              })
              output$downloadUMAP_resolution_page_5<- downloadHandler(
                filename="UMAP_resolution_p5.png",
                content=function(file){
                  png(file, width = 1500 , height = 1100,res = 150)
                  print(vizu_UMAP(seuratObj,var=ident_obj_mult, dlabel = input$addlabels_MDM))
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
                      expressionMatrixEntire <- cbind(geneList = rowSums(FetchData(seuratObj, vars =c(input$GeneListPool))),seuratObj[[input$clusterwatch]])
                    }else{
                      expressionMatrixEntire <- cbind(geneList = rowSums(FetchData(seuratObj, vars =c(input$GeneListPool))),seuratObj[[input$annotationwatch]])
                    }
                  }else{
                    if(input$ResOrAnnotMult == "Resolution"){
                      expressionMatrixEntire <- cbind(geneList = rowMeans(FetchData(seuratObj, vars =c(input$GeneListPool))),seuratObj[[input$clusterwatch]])
                    }else{
                      expressionMatrixEntire <- cbind(geneList = rowMeans(FetchData(seuratObj, vars =c(input$GeneListPool))),seuratObj[[input$clusterwatch]])
                    }
                  }
                }else{
                  sbt <- AddModuleScore(seuratObj, features =list(input$GeneListPool), name = "Gene_list")
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
              poolGene <- function(){
                  validate(
                      need(input$GeneListPool != "","Choose at least one gene to display projection. This plot can pool some gene and get the sum of the normalized counts. \nAddModuleScore function will calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.")
                  )
                  if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                      if (input$SumorMeans == "Sum"){
                          pool_genes <- rowSums(FetchData(seuratObj, vars =c(input$GeneListPool)))  
                      }else if (input$SumorMeans == "Mean"){
                          pool_genes <- rowMeans(FetchData(seuratObj, vars =c(input$GeneListPool)))
                      }
                      pool_genes <- cbind.data.frame(coord,counts = pool_genes)
                      ggplot(pool_genes)+geom_point(aes(UMAP_1,UMAP_2, color = counts),shape=20, size=0.25)+ scale_color_gradient(low="lightgrey",high="blue") + theme_classic() 
                  }else{
                      sbt <- AddModuleScore(seuratObj, features =list(input$GeneListPool), name = "Gene_list")
                      FeaturePlot(sbt, features = "Gene_list1")
                  }
                  
              }
              output$FeaturePlotMultGenesPooled <- renderPlot({
                  poolGene()
              })
              
              VlnPlotPooled <- function(){
                  validate(
                      need(input$GeneListPool != "", "Choose at least one gene to display violin plot")
                  )
                  if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                      if(input$SumorMeans == "Sum"){
                        if(input$ResOrAnnotMult == "Resolution"){
                          pool_genesVln <- cbind(seuratObj[[input$clusterwatch]],rowSums(FetchData(seuratObj, vars = input$GeneListPool)))
                        }else{
                          pool_genesVln <- cbind(seuratObj[[input$annotationwatch]],rowSums(FetchData(seuratObj, vars = input$GeneListPool)))
                        }
                      }else{
                        if(input$ResOrAnnotMult == "Resolution"){
                          pool_genesVln <- cbind(seuratObj[[input$clusterwatch]],rowMeans(FetchData(seuratObj, vars = input$GeneListPool)))
                        }else{
                          pool_genesVln <- cbind(seuratObj[[input$annotationwatch]],rowMeans(FetchData(seuratObj, vars = input$GeneListPool)))
                        }
                      }
                      colnames(pool_genesVln) <- c("resolution","GeneList")
                      ggplot(pool_genesVln, aes(x=resolution,y=GeneList, fill = resolution))+geom_violin()+scale_fill_hue()
                  }else{
                      sbt <- AddModuleScore(seuratObj, features =list(input$GeneListPool), name = "Gene_list")
                      VlnPlot(sbt, features = "Gene_list1")
                  }
              }
              output$ViolinPlotMultGenesPooled <- renderPlot({
                  VlnPlotPooled()
              })
              
              
              DotplotPooled <- function(){
                  validate(
                      need(input$GeneListPool != "", "Choose at least one gene to display dot plot")
                  )
                  if(input$SumorMeans == "Sum"){
                      if(input$ResOrAnnotMult == "Resolution"){
                        Idents(seuratObj) <- input$clusterwatch
                      }else{
                        Idents(seuratObj) <- input$annotationwatch
                      }
                      test <- as.data.frame(cbind(rowSums(FetchData(seuratObj, vars =c(input$GeneListPool))),seuratObj[[input$clusterwatch]]))
                      colnames(test)<- c("V1","V2")
                      percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                      colnames(percentExpressed) <- c("cluster", "logical","number")
                      percentExpressed <- filter(percentExpressed, logical == "TRUE")
                      totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                      percent <- percentExpressed$number/totalCellsbyCluster$n
                      percent <- as.data.frame(cbind(cluster = levels(seuratObj),RatioCellsExpressingGeneList = percent))
                      avgExpression <-test %>% group_by(V2) %>%summarise(avg =sum(V1))
                      percent <- cbind(percent, SumExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                      percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                      ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = SumExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                      
                  }else if(input$SumorMeans=="Mean"){
                      if(input$ResOrAnnotMult == "Resolution"){
                        Idents(seuratObj) <- input$clusterwatch
                      }else{
                        Idents(seuratObj) <- input$annotationwatch
                      }
                      test <- as.data.frame(cbind(rowMeans(FetchData(seuratObj, vars =c(input$GeneListPool))),seuratObj[[input$clusterwatch]]))
                      colnames(test)<- c("V1","V2")
                      percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                      colnames(percentExpressed) <- c("cluster", "logical","number")
                      percentExpressed <- filter(percentExpressed, logical == "TRUE")
                      totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                      percent <- percentExpressed$number/totalCellsbyCluster$n
                      percent <- as.data.frame(cbind(cluster = levels(seuratObj),RatioCellsExpressingGeneList = percent))
                      avgExpression <-test %>% group_by(V2) %>%summarise(avg =mean(V1))
                      percent <- cbind(percent, AvgExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                      percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                      ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = AvgExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                  }else{
                      sbt <- AddModuleScore(seuratObj, features =list(input$GeneListPool), name = "Gene_list")
                      DotPlot(sbt, features = "Gene_list1")
                  }
                  
              }
              output$DotPlotMultGenePooled <- renderPlot({
                  DotplotPooled()
              })
              ### Download plot ###
              if(input$format2 =="PNG"){
                  output$VPdownloadMultGenesPooled <- downloadHandler(
                      filename = "ViolinPlotMultGenesPooled.png",
                      content = function(file){
                          png(file, width = 1500 , height = 1100,res = 150)
                          print(VlnPlotPooled())
                          dev.off()
                      }
                  )
              }else{
                  output$VPdownloadMultGenesPooled <- downloadHandler(
                      filename = "ViolinPlotMultGenesPooled.svg",
                      content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(VlnPlotPooled())
                          dev.off()
                      }
                  )
              }
              
              
              if(input$format2 =="PNG"){
                  output$downloadfeatureplotMultGenesPooled<- downloadHandler(
                      filename = "featurePlotMultGenesPooled.png",
                      content = function(file){
                          png(file, width = 1500 , height = 1100,res = 150)
                          print(poolGene())
                          dev.off()
                      }
                  )  
              }else{
                  output$downloadfeatureplotMultGenesPooled<- downloadHandler(
                      filename = "featurePlotMultGenesPooled.svg",
                      content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(poolGene())
                          dev.off()
                      }
                  )
              }
              
              if(input$format2=="PNG"){
                  output$dotplotdownloadMultGenesPooled<- downloadHandler(
                      filename = "DotPlotMultGenesPooled.png",
                      content = function(file){
                          png(file, width = 1500 , height = 1100,res = 150)
                          print(poolGene())
                          dev.off()
                      }
                  )
              }else{
                  output$dotplotdownloadMultGenesPooled<- downloadHandler(
                      filename = "DotPlotMultGenesPooled.svg",
                      content = function(file){
                          svg(file, width = 14 , height = 7)
                          print(poolGene())
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
              Readfile <- function(){
                validate(
                  need(input$fileGeneList, " ")
                )
                list_gene_file <- read.csv(input$fileGeneList$datapath, header=TRUE)
                colnames(list_gene_file) <- "cluster"
                list_gene_file <- list_gene_file$cluster
                list_gene_file <- list_gene_file[list_gene_file %in% rownames(seuratObj)]
                list_gene_file
              }
              
              
              exprPoolFromFile <- reactive({
                list_gene_file <- Readfile()
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
                list_gene_file <- Readfile()
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
              
              
              poolGeneFromList <- function(){
                validate(
                  need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.This plot will show the pool of each gene expression for each gene id contained in the file. \nExample of a file :\nHeader\nMlxipl\nSox9\nEtc...  ")
                )
                coord <- Embeddings(seuratObj[["umap"]])[,1:2]
                gene_list <- Readfile()
                if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                  if(input$SumorMeans =="Sum"){
                    gene_list_pooled <- rowSums(FetchData(seuratObj, vars = gene_list))
                  }
                  if(input$SumorMeans =="Mean"){
                    gene_list_pooled <- rowMeans(FetchData(seuratObj, vars = gene_list))
                  }
                  gene_list_pooled <- cbind.data.frame(coord, counts = gene_list_pooled)
                  ggplot(gene_list_pooled)+geom_point(aes(UMAP_1,UMAP_2, color = counts),shape=20, size=0.25)+ scale_color_gradient(low="lightgrey",high="blue") + theme_classic()
                }else{
                  sbt1 <- AddModuleScore(seuratObj, features =gene_list, name = "Gene_list")
                  FeaturePlot(sbt1, features = "Gene_list1")
                }
              }
              
              output$FeaturePlotMultGenesPooledFilesInput <- renderPlot({
                poolGeneFromList()
              })
              
              poolGenefromListVP <- function(){
                validate(
                  need(input$fileGeneList,"Choose a file that contain a list of genes to display the projection.")
                )
                gene_list <- Readfile()
                if(input$SumorMeans =="Sum" | input$SumorMeans == "Mean"){
                  if(input$SumorMeans =="Sum"){
                    if(input$ResOrAnnotMult == "Resolution"){
                      pool_genesVln <- cbind(seuratObj[[input$clusterwatch]],rowSums(FetchData(seuratObj, vars = gene_list)))
                    }else{
                      pool_genesVln <- cbind(seuratObj[[input$annotationwatch]],rowSums(FetchData(seuratObj, vars = gene_list)))
                    }
                  }
                  if(input$SumorMeans =="Mean"){
                    if(input$ResOrAnnotMult == "Resolution"){
                      pool_genesVln <- cbind(seuratObj[[input$clusterwatch]],rowMeans(FetchData(seuratObj, vars = gene_list)))
                    }else{
                      pool_genesVln <- cbind(seuratObj[[input$annotationwatch]],rowMeans(FetchData(seuratObj, vars = gene_list)))
                    }
                  }
                  colnames(pool_genesVln) <- c("resolution","GeneList")
                  ggplot(pool_genesVln, aes(x=resolution,y=GeneList, fill = resolution))+geom_violin()+scale_fill_hue()
                }else{
                  if(input$ResOrAnnotMult == "Resolution"){
                    Idents(seuratObj) <- input$clusterwatch
                  }else{
                    Idents(seuratObj) <- input$annotationwatch
                  }
                  sbt1 <- AddModuleScore(seuratObj, features =gene_list, name = "Gene_list")
                  VlnPlot(sbt1, features = "Gene_list1")
                }
                
              }
              output$ViolinPlotMultGenesPooledFilesInput <- renderPlot({
                poolGenefromListVP()
              })
              
              DotplotPooledFromFile <- function(){
                validate(
                  need(input$fileGeneList, "Choose a file that contain a list of genes to display the projection.")
                )
                gene_list <- Readfile()
                if(input$SumorMeans == "Sum"){
                  if(input$ResOrAnnotMult == "Resolution"){
                    Idents(seuratObj) <- input$clusterwatch
                  }else{
                    Idents(seuratObj) <- input$annotationwatch
                  }
                  test <- as.data.frame(cbind(rowSums(FetchData(seuratObj, vars =gene_list)),seuratObj[[input$clusterwatch]]))
                  colnames(test)<- c("V1","V2")
                  percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                  colnames(percentExpressed) <- c("cluster", "logical","number")
                  percentExpressed <- filter(percentExpressed, logical == "TRUE")
                  totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                  percent <- percentExpressed$number/totalCellsbyCluster$n
                  percent <- as.data.frame(cbind(cluster = levels(seuratObj),RatioCellsExpressingGeneList = percent))
                  avgExpression <-test %>% group_by(V2) %>%summarise(avg =sum(V1))
                  percent <- cbind(percent, SumExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                  percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                  ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = SumExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                  
                }else if (input$SumorMeans == "Mean"){
                  if(input$ResOrAnnotMult == "Resolution"){
                    Idents(seuratObj) <- input$clusterwatch
                  }else{
                    Idents(seuratObj) <- input$annotationwatch
                  }
                  test <- as.data.frame(cbind(rowMeans(FetchData(seuratObj, vars =gene_list)),seuratObj[[input$clusterwatch]]))
                  colnames(test)<- c("V1","V2")
                  percentExpressed <- test %>% group_by(V2) %>%count(V1 > 0)
                  colnames(percentExpressed) <- c("cluster", "logical","number")
                  percentExpressed <- filter(percentExpressed, logical == "TRUE")
                  totalCellsbyCluster <- test %>% group_by(V2) %>% tally()
                  percent <- percentExpressed$number/totalCellsbyCluster$n
                  percent <- as.data.frame(cbind(cluster = levels(seuratObj),RatioCellsExpressingGeneList = percent))
                  avgExpression <-test %>% group_by(V2) %>%summarise(avg =mean(V1))
                  percent <- cbind(percent, AvgExp =avgExpression$avg, GeneList = rep("Gene_list",length(percent$cluster)))
                  percent$RatioCellsExpressingGeneList <- as.numeric(percent$RatioCellsExpressingGeneList)
                  ggplot(percent , aes(x= GeneList, y  = cluster , size = RatioCellsExpressingGeneList, color = AvgExp))+geom_point()+scale_color_gradient(low = "lightgrey", high= "blue")
                }else{
                  if(input$ResOrAnnotMult == "Resolution"){
                    Idents(seuratObj) <- input$clusterwatch
                  }else{
                    Idents(seuratObj) <- input$annotationwatch
                  }
                  sbt1 <- AddModuleScore(seuratObj, features =gene_list, name = "Gene_list")
                  DotPlot(sbt1, features = "Gene_list1")
                }
                
              }
              output$DotPlotMultGenePooledFromFile <- renderPlot({
                DotplotPooledFromFile()
              })
              
              
              
              if(input$format2=="PNG"){
                output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                  filename = "ViolinPlotMultGenesPooled.png",
                  content = function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(poolGenefromListVP())
                    dev.off()
                  }
                ) 
              }else{
                output$downloadViolinPlotMultGenesPooledFilesInput <- downloadHandler(
                  filename = "ViolinPlotMultGenesPooled.svg",
                  content = function(file){
                    svg(file, width = 14 , height = 7)
                    print(poolGenefromListVP())
                    dev.off()
                  }
                )
              }
              
              if(input$format2=="PNG"){
                output$downloadfeatureplotMultGenesPooledFilesInput <- downloadHandler(
                  filename = "featurePlotFilesInput.png",
                  content = function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(poolGeneFromList())
                    dev.off()
                  }
                )  
              }else{
                output$downloadfeatureplotMultGenesPooledFilesInput <- downloadHandler(
                  filename = "featurePlotFilesInput.svg",
                  content = function(file){
                    svg(file, width = 14 , height = 7)
                    print(poolGeneFromList())
                    dev.off()
                  }
                )
              }
              
              if(input$format2 == "PNG"){
                output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                  filename = "DotPlotMultGenesPooled.png",
                  content = function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(DotplotPooledFromFile())
                    dev.off()
                  }
                )
              }else{
                output$downloadDotPlotMultGenesPooledFilesInput<- downloadHandler(
                  filename = "DotPlotMultGenesPooled.svg",
                  content = function(file){
                    svg(file, width = 14 , height = 7)
                    print(DotplotPooledFromFile())
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
                    vizu_UMAP(seuratObj,var = input$clusterNumber, dlabel = input$addlabels_info)
                })
                output$downloadUMAP_info<- downloadHandler(
                  filename = "UMAP_info.png",
                  content = function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(seuratObj,var = input$clusterNumber, dlabel = input$addlabels_info))
                    dev.off()
                  }
                )
                toListen <- reactive({list(input$clusterInfo,input$chooseVar2Plot, input$stackedBP, input$FreqOrVal)})
                observeEvent(toListen(), {
                    Plot2Render <- function(){
                        if(input$clusterInfo == "All"){
                            if(input$stackedBP == "No"){
                                table_tmp2 <-as.data.frame(table(c(seuratObj[[input$chooseVar2Plot]], seuratObj[[input$clusterNumber]])))
                                colnames(table_tmp2) <- c("condition","cluster", "freq")
                                test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
                                g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = position_dodge())+scale_fill_manual(values = test)
                            }else{
                                if(input$FreqOrVal =="Frequence"){
                                    table_tmp2 <-as.data.frame(table(c(seuratObj[[input$chooseVar2Plot]], seuratObj[[input$clusterNumber]])))
                                    colnames(table_tmp2) <- c("condition","cluster", "freq")
                                    test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
                                    g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "fill")+scale_fill_manual(values = test)
                                }else{
                                    table_tmp2 <-as.data.frame(table(c(seuratObj[[input$chooseVar2Plot]], seuratObj[[input$clusterNumber]])))
                                    colnames(table_tmp2) <- c("condition","cluster", "freq")
                                    test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
                                    g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = test)
                                }
                            }
                        }else{
                            table_tmp <- as.data.frame(table(c(seuratObj[[input$chooseVar2Plot]],seuratObj[[input$clusterNumber]])))
                            colnames(table_tmp) <- c("condition","cluster", "freq")
                            table_tmp <- subset(table_tmp, table_tmp$cluster == input$clusterInfo)
                            test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp$condition)))
                            g <- ggplot(table_tmp, aes(x = condition,freq, fill = condition ))+geom_bar(stat = "identity")+scale_fill_manual(values = test)
                        }
                    }
                    percentRender <- function(){
                        table_tmp2 <-as.data.frame(table(c(seuratObj[[input$chooseVar2Plot]], seuratObj[[input$clusterNumber]])))
                        colnames(table_tmp2) <- c("condition","cluster", "freq")
                        tmp <- as.data.frame(rep(table(seuratObj[[input$chooseVar2Plot]]),length(unique(seuratObj[[input$clusterNumber]]))))
                        colnames(tmp)<- "total"
                        table_tmp2 <- cbind(table_tmp2,tmp)
                        table_tmp2$percent <- table_tmp2$freq/table_tmp2$total*100
                        test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$cluster)))
                        g <- ggplot(table_tmp2, aes(x= condition, percent, fill = cluster))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = test)+theme(axis.text.x = element_text(angle =60,hjust=1))
                    }
                    output$ggplot_information <- renderPlotly({
                        p <- Plot2Render()
                        ggplotly(p +theme(legend.position = "none",axis.text.x = element_text(angle =60,hjust=1)))
                    })
                    output$legend_information <- renderPlot({
                        p <- Plot2Render()
                        grid.newpage()
                        legend <- cowplot::get_legend(p)
                        grid.draw(legend)
                    })
                    output$cluster_percent <- renderPlotly({
                        g <- percentRender()
                        ggplotly(g)
                    })
                    output$downloadInformationPlot2 <- downloadHandler(
                      filename = "Percent_Plot.png",
                      content = function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(percentRender())
                        dev.off()
                      }
                    )
                    Table2Render <- function(){
                        if(input$clusterInfo == "All"){
                            as.data.frame(table(c(seuratObj[[input$chooseVar2Plot]], seuratObj[[input$clusterNumber]])))
                        }else{
                            as.data.frame(table(seuratObj[[input$chooseVar2Plot]][which(seuratObj[[input$clusterNumber]]==input$clusterInfo),]))
                        }
                    }
                    output$get_info <- DT::renderDataTable({
                        Table2Render()
                    })
                    output$downloadInformationPlot <- downloadHandler(
                        filename = "Information_Plot.png",
                        content = function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(Plot2Render())
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
                    vizu_UMAP(seuratObj,var = input$res4Annot, dlabel = input$labelsAnnot)
                })
                output$downloadUMAP_resolution_annot_page<- downloadHandler(
                  filename="UMAP_resolution_annot_page.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
                    print(vizu_UMAP(seuratObj,var=input$res4Annot, dlabel = input$labelsAnnot))
                    dev.off()
                  })
                
            })
            
            ##render plot displaying variable of interest when use update instead of create
            output$condPlot <- renderPlot({
                vizu_UMAP(seuratObj, var = input$Annot2Complete, dlabel = input$labelsconditional)
            })
            output$downloadUMAP_Conditional<- downloadHandler(
              filename="UMAP_conditional.png",
              content=function(file){
                png(file, width = 1500 , height = 1100,res = 150)
                print(vizu_UMAP(seuratObj,var=input$Annot2Complete, dlabel = input$labelsconditional))
                dev.off()
              })
            
            observeEvent(input$submitAnnot,{
                if(input$choicesCreateOrComplete =="Create"){
                    if(input$AnnotName == "" | input$clusterName == "" | length(input$cluster4Annot) == 0){
                        shinyalert("Oops", "The name of the column and the name of the annotation you want to give has to be filled", type = "error")
                    }else{
                        if(is.element(input$clusterName, names(seuratObj@meta.data))){
                            f.clusters <- seuratObj[[input$clusterName]]
                            f.clusters[seuratObj[[resChoice]] == input$cluster4Annot[1]] <- input$AnnotName
                            if(length(input$cluster4Annot) > 1){
                                for(i in 2:length(input$cluster4Annot)){
                                    f.clusters[seuratObj[[resChoice]] == input$cluster4Annot[i]] <- input$AnnotName 
                                }
                            }
                            seuratObj[[input$clusterName]] <<- f.clusters
                        }else{
                            f.clusters <- rep("unknown", dim(seuratObj)[2])
                            names(f.clusters) <- colnames(seuratObj)
                            f.clusters[seuratObj[[resChoice]] == input$cluster4Annot[1]] <- input$AnnotName
                            if(length(input$cluster4Annot) > 1){
                                for(i in 2:length(input$cluster4Annot)){
                                    f.clusters[seuratObj[[resChoice]] == input$cluster4Annot[i]] <- input$AnnotName 
                                }
                            }
                            seuratObj[[input$clusterName]] <<- f.clusters
                        }
                        shinyalert("Done","You have correctly annotated your cluster", type = "success")
                        output$postAnnotation <- renderPlot({
                            vizu_UMAP(seuratObj, var = input$clusterName, dlabel = T)
                        }) 
                        output$downloadUMAP_post_annotation<- downloadHandler(
                          filename="UMAP_post_annot.png",
                          content=function(file){
                            png(file, width = 1500 , height = 1100,res = 150)
                            print(vizu_UMAP(seuratObj,var=input$Annot2Complete, dlabel = T))
                            dev.off()
                          })
                        toremove <- c(grep(names(seuratObj@meta.data),pattern ="nCount_*", value =T),grep(names(seuratObj@meta.data),pattern ="nFeature_*", value =T),grep(names(seuratObj@meta.data),pattern ="snn_res.*", value =T),"percent.mt","S.Score", "G2M.Score")
                        choices <- names(seuratObj@meta.data)
                        choices_annot <- choices[-which(choices %in% toremove )]
                        updateSelectizeInput(session,"chooseVar2Plot", choices=choices)
                        updateSelectizeInput(session,"whichAnnot", choices=choices_annot)
                        updateSelectizeInput(session,"Annot2Complete", choices =choices_annot)
                        updateSelectizeInput(session, "AnnotationDataMining", choices = choices_annot)
                        updateSelectizeInput(session,"annotationwatch",choices = choices_annot)
                        updateSelectizeInput(session,"AnnotAgainst", choices = choices_annot)
                        
                    }
                }else{
                    if(input$AnnotName == "" | length(input$cluster4Annot) == 0){
                        shinyalert("Oops", "The name of the annotation and the cluster(s) you want to give has to be filled", type = "error")
                    }else{
                        f.clusters <- seuratObj[[input$Annot2Complete]]
                        f.clusters[seuratObj[[resChoice]] == input$cluster4Annot[1]] <- input$AnnotName
                        if(length(input$cluster4Annot) > 1){
                            for(i in 2:length(input$cluster4Annot)){
                                f.clusters[seuratObj[[resChoice]] == input$cluster4Annot[i]] <- input$AnnotName 
                            }
                        }
                        seuratObj[[input$Annot2Complete]] <<- f.clusters
                    }
                    output$postAnnotation <- renderPlot({
                        vizu_UMAP(seuratObj, var = input$Annot2Complete, dlabel = input$labelsPostAnnot)
                    })
                    
                    
                    shinyalert("Done","You have correctly annotated your cluster", type = "success")
                    output$downloadUMAP_post_annotation<- downloadHandler(
                      filename="UMAP_post_annot.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(seuratObj,var=input$Annot2Complete, dlabel = input$labelsPostAnnot))
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
                        vizu_UMAP(seuratObj,var = input$clusterResPage7, dlabel=input$labelsSubset)
                    })
                    output$downloadUMAP_sbt<- downloadHandler(
                      filename="UMAP_resolution_for_sbt.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(seuratObj,var=input$clusterResPage7, dlabel = input$labelsSubset))
                        dev.off()
                      })
                    observeEvent(input$subclusterize,{
                        test <- NULL
                        if(length(input$subclusterize) == length(clusterSbt)){
                            test <- colnames(seuratObj)
                            output$Dimplot_kept <- renderPlot({
                              MakeUMAP_Red(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                            })
                            output$downloadUMAPcellKept<- downloadHandler(
                              filename="UMAP_cell_kept.png",
                              content=function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(MakeUMAP_Red(seuratObj,test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                dev.off()
                            })
                        }else{
                            test <- colnames(seuratObj)[which(seuratObj[[]][input$clusterResPage7] == input$subclusterize[1])]
                            output$Dimplot_kept <- renderPlot({
                                MakeUMAPhighlight(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                            })
                            output$downloadUMAPcellKept<- downloadHandler(
                              filename="UMAP_cell_kept.png",
                              content=function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
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
                                  filename="UMAP_cell_kept.png",
                                  content=function(file){
                                    png(file, width = 1500 , height = 1100,res = 150)
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
                                      filename="UMAP_cell_kept.png",
                                      content=function(file){
                                        png(file, width = 1500 , height = 1100,res = 150)
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
                        vizu_UMAP(seuratObj,var = input$whichAnnot, dlabel = input$labelsSubset)
                    })
                    output$downloadUMAP_sbt<- downloadHandler(
                      filename="UMAP_resolution_for_sbt.png",
                      content=function(file){
                        png(file, width = 1500 , height = 1100,res = 150)
                        print(vizu_UMAP(seuratObj,var=input$whichAnnot, dlabel = input$labelsSubset))
                        dev.off()
                      })
                    observeEvent(input$subannot,{
                        test <- NULL
                        if(length(input$subannot) == length(annotSbt)){
                            test <- colnames(seuratObj)
                            output$Dimplot_kept <- renderPlot({
                              MakeUMAP_Red(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                            })
                            output$downloadUMAPcellKept<- downloadHandler(
                              filename="UMAP_cell_kept.png",
                              content=function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
                                print(MakeUMAP_Red(seuratObj,test, sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize))
                                dev.off()
                              })
                            
                        }else{
                            test <- colnames(seuratObj)[which(seuratObj[[]][input$whichAnnot] == input$subannot[1])]
                            output$Dimplot_kept <- renderPlot({
                              MakeUMAPhighlight(seuratObj, highliht_cells = test,sizePoint = input$cellsSize, sizeHighlight = input$cellsHighlightSize)
                            })
                            output$downloadUMAPcellKept<- downloadHandler(
                              filename="UMAP_cell_kept.png",
                              content=function(file){
                                png(file, width = 1500 , height = 1100,res = 150)
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
                                  filename="UMAP_cell_kept.png",
                                  content=function(file){
                                    png(file, width = 1500 , height = 1100,res = 150)
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
                                      filename="UMAP_cell_kept.png",
                                      content=function(file){
                                        png(file, width = 1500 , height = 1100,res = 150)
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
                    vizu_UMAP(subsetSeuratObj, var = NULL, dlabel = input$labelsAfterSubset)
                })
                output$downloadUMAPaftersbt<- downloadHandler(
                  filename="UMAP_after_sbt.png",
                  content=function(file){
                    png(file, width = 1500 , height = 1100,res = 150)
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
                    tempReport <- file.path(tempdir(),"logFile.Rmd")
                    file.copy("logFile.Rmd", tempReport,overwrite = T)
                    if(input$subcluster == "Cluster"){
                        line <- paste0("\tsubsetSeuratObj <- subset(seuratObj, subset =",input$clusterResPage7," == ",list(input$subclusterize),")")
                    }else{
                        line <- paste0("\tsubsetSeuratObj <- subset(seuratObj, subset = ",input$whichAnnot," == ",list(input$subannot),")")
                    }
                    write(line,tempReport, append =TRUE)
                    write("\tsubsetSeuratObj <- RunPCA(subsetSeuratObj)",tempReport, append =TRUE)
                    write("\tsubsetSeuratObj <- FindNeighbors(subsetSeuratObj)", tempReport, append = TRUE)
                    write("\tsubsetSeuratObj <- RunUMAP(subsetSeuratObj, dims = 1:30)\n", tempReport, append = TRUE)
                    write("\tfor(i in seq(0,1,0.1)){",tempReport, append = T)
                    write("\t\tsubsetSeuratObj <- findClusters(subsetSeuratObj, res = i)",tempReport, append = T)
                    write("\t}\n",tempReport, append = T)
                    rmarkdown::render(tempReport,output_file = file, params = command, envir = new.env(parent = globalenv()))
                }
            )
        }
    })
}