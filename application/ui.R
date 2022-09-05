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
                                includeMarkdown("../Documentation/Help.md"),
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
                                           plotOutput("barplotBeforeAfter"),
                                           plotOutput("cellcycle"),
                                           plotOutput("clusteringPlot")),
                                    column(6,
                                           plotOutput("AfterFiltering"),
                                           plotOutput("nbcells_by_datasets"),
                                           plotOutput("geneExpByCell")
                                           
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
                                        plotOutput("cellcycleQC"),
                                        plotOutput("plotClusterQC")
                                    ),
                                    column(6,
                                           plotOutput("nbcells_by_datasetsQC"),
                                           plotOutput("geneExpByCellQC"),
                                           plotOutput("plotChoice"))))
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
                                           plotOutput("DimplotTreePage")),
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
                                    column(6,plotOutput("DimplotMarkers"))),
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
                                radioButtons("keepscale","Do you want to show the feature plot with the same scale ? ", choices=c("Yes", "No")),
                                radioButtons("format", "Download plot as : ", choices=c("PNG", "SVG")),
                                downloadButton('downloadMatrix',"Download entire matrix"),
                                width = 2
                            ),
                            mainPanel(

                                tabsetPanel(
                                    tabPanel("Feature plot",
                                             plotOutput("FeaturePlotMultGenes"),
                                             downloadButton('downloadfeatureplot',"Download plot")
                                    ),
                                    tabPanel("Violin plot",
                                             plotOutput('ViolinPlot'),
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
                                ),
                                fluidRow(column(6,DT::dataTableOutput("geneExp"),
                                                downloadButton('downloadRawCount', 'Download')),
                                         column(6,plotOutput('clusterOutput'))
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
                               column(6,plotOutput("dimplotSeurat4page"))
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
                           selectizeInput("chooseVar2Plot", "Which metadata has to be shown ? ", choices=character(0))
                           ,width=2),
                       mainPanel(
                           fluidRow(
                               column(6,
                                      plotOutput("Dimplot_cluster"),
                                      plotlyOutput("ggplot_information"),
                                      downloadButton('downloadInformationPlot', "Download plot"),
                                      plotlyOutput("cluster_percent")
                                    
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
                                conditionalPanel(condition="input.choicesCreateOrComplete == 'Create'",textInput("clusterName","Type of the annotation ? Name it without space (ex: Cell_type, ...)")),
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
                                           
                                    ),
                                    column(6,
                                           conditionalPanel(condition = "input.choicesCreateOrComplete == 'Update'", plotOutput("condPlot")),
                                           plotOutput("postAnnotation"))
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
                                           titlePanel("After Subclustering"),
                                           plotOutput("dimplot_aftersbt")),
                                    column(6, 
                                           plotOutput("Dimplot_kept"))
                                )
                            )
                        )
                )
               
               
    )
    
)

