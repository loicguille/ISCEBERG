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
library(rmarkdown)

options(shiny.maxRequestSize =100000*1024^2)
# Define UI for application that draws a histogram
ui <- fluidPage(
    theme = bslib::bs_theme(bootswatch = "flatly"),
    useShinyjs(),
    # Application title
    titlePanel(title="ISCEBERG"),
    
    navbarPage("",id = "tabs",
               tabPanel("Pre-processing", icon = icon("upload"),
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("file", "Type of file", choices = c("Mtx","CSV","H5","RDS","Txt"), selected = NULL),
                                conditionalPanel(condition = "input.file != 'Mtx'",fileInput("InputFile", "Upload file (don't forget to load the data once this is over)", multiple = TRUE), fileInput("Metadata", "Would you add metadata ? (optional)", multiple = TRUE)),
                                conditionalPanel(condition = "input.file == 'Mtx' ", fileInput("MatrixFile", "Upload Matrix", multiple = TRUE), fileInput("GeneFile", "Upload gene or feature file", multiple = TRUE), fileInput("CellsFile", "Upload barcode", multiple = TRUE)),
                                conditionalPanel(condition = "input.file == 'CSV'  | input.file == 'Txt' ", radioButtons("header", "Is the file(s) containing header(s) ?", choices =c("Yes", "No")),radioButtons("separator", "What is the separator of your file ?", choices = c(",", ";","tab"))),
                                br(),
                                conditionalPanel(condition = "input.file != 'RDS' ", numericInput("MinGeneByCells", "How many cells a gene has to be found in to be kept", min = 0,max = 100 , value = 3),numericInput("MinCells", "How many genes a cells has to express to be kept", value = 100, min = 0, max = 1000, step = 10)),
                                actionButton("createSeurat","Run Pre-processing"),
                                actionButton("Refresh","Refresh data"),
                                width = 2
                            ),
                            mainPanel(
                                column(6,uiOutput('Isceberg')),
                                fluidPage(
                                  htmltools::tags$iframe(src = "Help.html", width = '100%',  height = 1000,  style = "border:none;"))
                            ,
                                fluidRow(
                                  titlePanel("Ressource consumption"),
                                  column(6,
                                         uiOutput("memoryCons")),
                                  column(6,
                                         uiOutput("timeCons")))
                            )
                        ),
                        
               ),
               tabPanel("filtering", icon = icon("filter"),
                        sidebarLayout(
                            sidebarPanel(
                                numericInput("percentMito", "Percent-mito max per cells : ", min = 0, max = 100, value = 25),
                                numericInput("minGenePerCells","Minimum number of genes expressed per cells : ", min = 0, max = 2500, value = 500, step = 100),
                                numericInput("maxGenePerCells","Maximum number of genes expressed per cells : ", min = 1000, value = 4500, step = 100),
                                numericInput("maxCountPerCells","Maximum number of reads by cells : ", min = 20000, value = 40000, step = 2000),
                                radioButtons("normalization", "What kind of normalization do you want to use on your data : ", choices = c("SCTransform","LogNormalisation")),
                                sliderInput("Resolution","Minimum and maximum resolution for clustering : ", min = 0, max = 5, value= c(0,1),step = 0.05),
                                numericInput("ResolutionStep", "Resolution step : ", min = 0.05, max = 1, value = 0.1, step = 0.05),
                                actionButton("runFiltNorm", "Run analysis"),
                                br(),
                                textInput("NameSeuratLog", "Give a name to your seurat object (mandatory)"),
                                br(),
                                downloadButton("downloadSeurat","Download seurat object"),
                                br(),
                                br(),
                                downloadButton("downloadLogFile","Download log file"),
                                width = 2
                            ),
                            mainPanel(
                                fluidRow(
                                    column(6,
                                           plotOutput("BeforeFiltering"),
                                           downloadButton("downloadBeforeFiltering", "Download plot"),
                                           plotOutput("barplotBeforeAfter"),
                                           downloadButton("download_barplot_before_after", "Download plot"),
                                           plotOutput("cellcycle"),
                                           fluidRow(column(3,downloadButton("download_cell_cycle", "Download UMAP")), column(3, checkboxInput("showlabelcc","Add label to plot", value = TRUE))),
                                           plotOutput("clusteringPlot"),
                                           fluidRow(column(3,downloadButton("download_clustering_plot", "Download UMAP")), column(3, checkboxInput("showlabelcluster","Add label to plot", value = TRUE)))
                                           ),
                                    column(6,
                                           plotOutput("AfterFiltering"),
                                           downloadButton("downloadAfterFiltering", "Download plot"),
                                           plotOutput("nbcells_by_datasets"),
                                           downloadButton("download_nb_cells_dt", "Download plot"),
                                           plotOutput("geneExpByCell"),
                                           downloadButton("download_UMAP_nb_genes_cells", "Download UMAP"),
                                           plotOutput("BatchVizu"),
                                           fluidRow(column(3,downloadButton("downloadbatchVizu", "Download plot")), column(3,checkboxInput("showlabelBatch","Add label to plot", value = TRUE)))
                                           
                                    )
                                
                            )
                        )
                        )
               ),
               tabPanel("QC", icon = icon("ruler"),
                        sidebarLayout(
                            sidebarPanel(id = "sidebar",
                                         selectizeInput("InfoToPlot","What type of informations do you want to project ?", choices = character(0)),width = 2),
                            mainPanel(
                                fluidRow(
                                    column(6,
                                        plotOutput("featureQC"),
                                        downloadButton('downloadVln', "Download plot"),
                                        plotOutput("cellcycleQC"),
                                        fluidRow(column(3,downloadButton('downloadUMAP_CC', "Download UMAP")),column(3,checkboxInput("showlabelccQC","Add label to plot", value = TRUE))),
                                        plotOutput("plotClusterQC"),
                                        fluidRow(column(3,downloadButton('downloadUMAP_resolution_qc', "Download UMAP")), column(3, checkboxInput("addlabels_res", " Add labels to plot", value = TRUE)))
                                    ),
                                    column(6,
                                       plotOutput("nbcells_by_datasetsQC"),
                                       downloadButton('downloadHistNbCells', "Download plot"),
                                       plotOutput("geneExpByCellQC"),
                                      downloadButton('downloadUMAP_gene_by_cell', "Download UMAP"),
                                       plotOutput("plotChoice"),
                                       fluidRow(column(3,downloadButton('downloadUMAP_choice', "Download UMAP")), column(3, checkboxInput("addlabels_choice", " Add labels to plot", value = TRUE)))))
                                )
                        )
               ),
               tabPanel("Cluster tree", icon = icon("sitemap"),
                        sidebarLayout(
                            sidebarPanel(
                                selectizeInput("resolution_TreePage","What resolution do you want to apply on the dimplot ?" ,choices = character(0)),
                                width = 2
                            ),
                            mainPanel(
                                plotOutput("ClusterTree"),
                                fluidRow(
                                    column(6,
                                           plotOutput("DimplotTreePage"),
                                           fluidRow(column(3,downloadButton('downloadUMAP_Tree_page', "Download UMAP")), column(3, checkboxInput("addlabels_tree", " Add labels to plot", value = TRUE)))),
                                    column(6,
                                           DTOutput("nbCellsByClust_Tree_Page"))
                                )
                            )
                            
                        )
               ),
               tabPanel("DE between clusters", icon = icon("chart-pie"),
                        sidebarLayout(
                            sidebarPanel(
                                selectizeInput("resolutionDE","What resolution do you want to have ?", choices= character(0)),
                                selectizeInput("Cluster1","From which cluster do you want the markers",choices=character(0)),
                                selectizeInput("Cluster2","Versus which group do you want to find the markers", choices=character(0)),
                                radioButtons("onlyPos", "Do you want only the positive markers", choices = c("Yes","No"), selected = "Yes"),
                                numericInput("LogThreshold","Value of minimum threshold difference (Log_scale)", value = 0.25, min = 0.1, max = 3,step = 0.05),
                                numericInput("percent","Minimal fraction of the cells that have to contain a gene in order to test it (default 0.2 = 20%)", value = 0.2, min = 0,max = 1, step = 0.05),
                                actionButton("submit","Load"),
                                width = 2
                            ),
                            mainPanel(
                                fluidRow(
                                    column(6,plotOutput("DimplotMarkers"),
                                           fluidRow(column(3,downloadButton('downloadUMAP_DE_page', "Download UMAP")),column(3,checkboxInput("addlabels_DE", " Add labels to plot", value = TRUE))))),
                                withSpinner(DT::dataTableOutput("Markers")) 
                            )
                        )
               ),
               tabPanel("Data mining for one gene",icon = icon("search"),
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("ResOrAnnot", "Data you want to use", choices =c("Resolution","Annotation")),
                                conditionalPanel(condition = "input.ResOrAnnot == 'Resolution'" ,selectizeInput("ResolutionDataMining","Select a resolution to use to calculate the number of clusters", choices = character(0))),
                                conditionalPanel(condition = "input.ResOrAnnot == 'Annotation'" ,selectizeInput("AnnotationDataMining","Select an annotation to use for vizualization", choices = character(0))),
                                selectizeInput("GeneList", "Select a list of genes that you want to plot (max 10)", choices=character(0), multiple = T,  options = list(maxItems=10)),
                                checkboxInput("SplitVln","Split plot (only available for violin plot) ", value = F),
                                conditionalPanel(condition = " input.SplitVln == true", selectizeInput("AnnotAgainst","Group for split", choices = character(0))),
                                conditionalPanel(condition = " input.SplitVln == true", selectizeInput("labelsToKeep1","Labels to keep", multiple = T, choices = character(0))),
                                conditionalPanel(condition = " input.SplitVln == true", selectizeInput("labelsToKeep2","Labels to keep", multiple = T, choices = character(0))),
                                radioButtons("format", "Download plot as : ", choices=c("PNG", "SVG")),
                                downloadButton('downloadMatrix',"Download entire matrix"),
                                width = 2
                            ),
                            mainPanel(
                              fluidRow(column(6,DT::dataTableOutput("geneExp"),
                                              downloadButton('downloadRawCount', 'Download')),
                                       column(6,plotOutput('clusterOutput'),
                                              fluidRow(column(3,downloadButton('downloadUMAP_Cluster_resolution', "Download UMAP")), column(3,checkboxInput("addlabels_ODM", " Add labels to plot", value = TRUE))))
                              ),
                                tabsetPanel(
                                    tabPanel("Feature plot",
                                             fluidRow(
                                               column(6,plotOutput('FeaturePlotMultGenes1')),
                                               column(6,plotOutput('FeaturePlotMultGenes2'))
                                             ),
                                             fluidRow(
                                               column(6,plotOutput('FeaturePlotMultGenes3')),
                                               column(6,plotOutput('FeaturePlotMultGenes4'))
                                             ),
                                             fluidRow(
                                               column(6,plotOutput('FeaturePlotMultGenes5')),
                                               column(6,plotOutput('FeaturePlotMultGenes6'))
                                             ),         
                                             fluidRow(
                                               column(6,plotOutput('FeaturePlotMultGenes7')),
                                               column(6,plotOutput('FeaturePlotMultGenes8'))
                                             ),
                                             fluidRow(
                                               column(6,plotOutput('FeaturePlotMultGenes9')),
                                               column(6,plotOutput('FeaturePlotMultGenes10'))
                                             ),
                                             downloadButton('downloadfeatureplot',"Download plot")
                                    ),
                                    tabPanel("Violin plot",
                                             fluidPage(
                                               fluidRow(
                                                 column(6,plotOutput('ViolinPlot1')),
                                                 column(6,plotOutput('ViolinPlot2'))
                                               ),
                                               fluidRow(
                                                 column(6,plotOutput('ViolinPlot3')),
                                                 column(6,plotOutput('ViolinPlot4'))
                                               ),
                                               fluidRow(
                                                 column(6,plotOutput('ViolinPlot5')),
                                                 column(6,plotOutput('ViolinPlot6'))
                                               ),         
                                               fluidRow(
                                                 column(6,plotOutput('ViolinPlot7')),
                                                 column(6,plotOutput('ViolinPlot8'))
                                               ),
                                               fluidRow(
                                                 column(6,plotOutput('ViolinPlot9')),
                                                 column(6,plotOutput('ViolinPlot10'))
                                               )),
                                             #plotOutput('ViolinPlot')),
                                             downloadButton("VPdownload","Download plot")
                                             
                                    ),
                                    tabPanel("Heatmap",
                                             plotOutput('Heatmap'),
                                             downloadButton("HMdownload","Download plot")
                                    ),
                                    tabPanel("Dotplot",
                                             plotOutput('Dotplot'),
                                             downloadButton("DPdownload","Download plot")
                                    )
                                )
                                  
                            )
                        )
               ),
               tabPanel(
                   "Data mining for a combination of gene", icon = icon("blender"),
                   sidebarLayout(
                       sidebarPanel(
                           h1("Gene List"),
                           radioButtons("ResOrAnnotMult", "Data you want to use", choices =c("Resolution","Annotation")),
                           conditionalPanel(condition = "input.ResOrAnnotMult == 'Resolution'" ,selectizeInput("clusterwatch","Select a resolution to use to calculate the number of clusters", choices = character(0))),
                           conditionalPanel(condition = "input.ResOrAnnotMult == 'Annotation'" ,selectizeInput("annotationwatch","Select an annotation to use for vizualization", choices = character(0))),
                           selectizeInput("GeneListPool", "Gene list you want to pool together", choices= character(0), multiple=T),
                           fileInput("fileGeneList","Upload gene list file"),
                           radioButtons("SumorMeans", "Do you want that the UMAP projection represent the sum or the mean of the expression of the genes list ?", choices=c("Sum","Mean","AddModuleScore Function"), selected = "Sum"),
                           radioButtons("format2", "Download plot as : ", choices=c("PNG", "SVG")),
                           downloadButton('downloadMatrixMultList',"Download entire matrix from list"),
                           br(),
                           br(),
                           downloadButton('downloadMatrixMultFile',"Download entire matrix from file"),
                           width=2),
                       mainPanel(
                           fluidRow(
                               column(6,titlePanel("Pool of a gene list"),
                                      tabsetPanel(
                                          tabPanel("Feature plot",
                                                   plotOutput("FeaturePlotMultGenesPooled"),
                                                   downloadButton('downloadfeatureplotMultGenesPooled',"Download plot")
                                          ),
                                          tabPanel("Violin plot",
                                                   plotOutput('ViolinPlotMultGenesPooled'),
                                                   downloadButton("VPdownloadMultGenesPooled","Download plot")
                                          ),
                                          tabPanel("Dotplot",
                                                   plotOutput("DotPlotMultGenePooled"),
                                                   downloadButton("dotplotdownloadMultGenesPooled", "Download plot"))
                                      ),
                                      DT::dataTableOutput("sumGeneexp"),
                                      downloadButton('downloadSumGeneExp', "Download matrix"),
                                      titlePanel("Pool of a gene list from an input file"),
                                        tabsetPanel(
                                            tabPanel("Feature plot",
                                                  plotOutput("FeaturePlotMultGenesPooledFilesInput"),
                                                  downloadButton('downloadfeatureplotMultGenesPooledFilesInput',"Download plot"),
                                            ),
                                            tabPanel("Violin plot",
                                               plotOutput("ViolinPlotMultGenesPooledFilesInput"),
                                               downloadButton('downloadViolinPlotMultGenesPooledFilesInput',"Download plot")),
                                            tabPanel("Dot plot",
                                                     plotOutput("DotPlotMultGenePooledFromFile"),
                                                     downloadButton('downloadDotPlotMultGenesPooledFilesInput',"Download plot"))
                                        ),
                                      DT::dataTableOutput("meanGeneExpPooledfromFile"),
                                      downloadButton('downloadmeanGeneExpPooledfromFile', "Download matrix"),),
                               column(6,plotOutput("dimplotSeurat4page"),
                                      fluidRow(column(3,downloadButton('downloadUMAP_resolution_page_5', "Download UMAP")),(column(3, checkboxInput("addlabels_MDM", " Add labels to plot", value = TRUE)))))
                           )
                           
                       )
                   )
               ),
               ####### PAGE 6 Information Plot
               tabPanel(
                   "Extract Information",icon = icon("atom"),
                   sidebarLayout(
                       sidebarPanel(
                           h1("Get information"),
                           selectizeInput("clusterNumber", "Select a resolution ", choices = character(0)),
                           selectizeInput("clusterInfo", "From which cluster do you want to provide informations ?", choices = character(0)),
                           conditionalPanel(condition = "input.clusterInfo == 'All'", radioButtons("stackedBP", "Do you want a stack graph ?", choices = c("Yes", "No"), selected = "No")),
                           conditionalPanel(condition = "input.stackedBP == 'Yes'", radioButtons("FreqOrVal", "Do you want frequence or values ?", choices = c("Frequence", "Values"), selected = "Frequence")),
                           selectizeInput("chooseVar2Plot", "Which metadata has to be shown ? ", choices=character(0)),
                           width=2),
                       mainPanel(
                           fluidRow(
                               column(6,
                                      plotOutput("Dimplot_cluster"),
                                      fluidRow(column(3,downloadButton('downloadUMAP_info', "Download UMAP")),(column(3, checkboxInput("addlabels_info", " Add labels to plot", value = TRUE)))),
                                      plotlyOutput("ggplot_information"),
                                      downloadButton('downloadInformationPlot', "Download plot"),
                                      plotlyOutput("cluster_percent"),
                                      downloadButton('downloadInformationPlot2', "Download plot")
                                    
                               ),
                               column(6,
                                      DT::dataTableOutput("get_info"),
                                      downloadButton('downloadMatrixInformation', "Download matrix"),
                                      plotOutput("legend_information"))
                           )
                       )
                   )
               ),
               ######### PAGE 7 write Annotation
               tabPanel("Add annotation",icon = icon("pen"),
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("choicesCreateOrComplete", "Do you want to create a new annotation or complete a former one ? ", choices = c("Create","Update")),
                                conditionalPanel(condition = "input.choicesCreateOrComplete == 'Update'",selectizeInput("Annot2Complete","Which annotation do you want to complete ?", choices = character(0))),
                                conditionalPanel(condition="input.choicesCreateOrComplete == 'Create'",textInput("clusterName","Type of the annotation ? Name it without space and - (ex: Cell_type, ...)")),
                                selectizeInput("res4Annot", "Select a resolution ", choices = character(0)),
                                selectizeInput("cluster4Annot", "Which cluster do you want to annotate ?", choices = character(0), multiple=T),
                                textInput("AnnotName","Labels of the annotation ? Name it without space (ex: Hepato, ...)"),
                                actionButton("submitAnnot", "Annotate"),
                                br(),
                                br(),
                                downloadButton("downloadRDSwithNewAnnot","Download Seurat Object"),
                                width=2    
                            ),
                        
                            mainPanel(
                                fluidRow(
                                    column(6,
                                           plotOutput("dimplot_res_annot"),
                                           fluidRow(column(3,downloadButton('downloadUMAP_resolution_annot_page', "Download UMAP")),column(3,checkboxInput("labelsAnnot","Add labels to plot", value = T)))
                                           
                                          ),
                                    column(6,
                                           plotOutput("postAnnotation"),
                                           fluidRow(column(3,downloadButton('downloadUMAP_post_annotation', "Download UMAP")),column(3,checkboxInput("labelsPostAnnot","Add labels to plot", value = T))),
                                           conditionalPanel(condition = "input.choicesCreateOrComplete == 'Update'", 
                                                            plotOutput("condPlot"),
                                                            fluidRow(column(3,downloadButton('downloadUMAP_Conditional', "Download UMAP")),column(3,checkboxInput("labelsconditional","Add labels to plot", value = T))))
                                           )
                                )
                            )
                        )
                ),
               ######### PAGE 8 subclustering
               tabPanel("Subclustering", icon = icon("cut"),
                        sidebarLayout(
                            sidebarPanel(
                                radioButtons("subcluster", "From what do you want to subcluster your data ?", choices =c("Cluster", "Annotation")),
                                conditionalPanel(condition = "input.subcluster == 'Cluster'", selectizeInput("clusterResPage7", "Select a resolution", choices = character(0))),
                                conditionalPanel(condition = "input.subcluster == 'Cluster'", selectizeInput("subclusterize","Which clusters do you want to get fro the subclustering ? ",multiple =T, choices = character(0))),
                                conditionalPanel(condition = "input.subcluster == 'Annotation'", selectizeInput("whichAnnot", "Select an Annotation", choices = character(0))),
                                conditionalPanel(condition = "input.subcluster == 'Annotation'", selectizeInput("subannot", "Which part of the annotation do you want to get for the subclustering ? ",multiple=T, choices = character(0))),
                                actionButton("dosubcluster", "Subclusterize"),
                                br(),
                                br(),
                                textInput("NameObject", "Name for your seurat object"),
                                downloadButton('downloadSubRds', "Download Subset Object"),
                                br(),
                                br(),
                                downloadButton('downloadLog', "Download log report"),width=2),
                            mainPanel(
                                fluidRow(
                                    column(6, 
                                           plotOutput("Dimplot_subclustering"),
                                           fluidRow(column(3,downloadButton('downloadUMAP_sbt', "Download UMAP")),column(3,checkboxInput("labelsSubset","Add labels to plot", value = T))),
                                           titlePanel("After Subclustering"),
                                           plotOutput("dimplot_aftersbt"),
                                           fluidRow(column(3,downloadButton('downloadUMAPaftersbt', "Download UMAP")),column(3,checkboxInput("labelsAfterSubset","Add labels to plot", value = T)))
                                           ),
                                    column(6, 
                                           plotOutput("Dimplot_kept"),
                                           fluidRow(column(3,downloadButton('downloadUMAPcellKept', "Download UMAP")),column(3,numericInput("cellsHighlightSize","Highlited cells size", min= 0, max =5,step = 0.25, value = 0.25)), column(3,numericInput("cellsSize","Cells size", min= 0, max =5,step = 0.25, value = 0.5))))
                                )
                            )
                        )
                )
               
               
    )
    
)

