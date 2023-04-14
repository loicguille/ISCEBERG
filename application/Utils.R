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




vizu_UMAP <- function(obj, var, dlabel = TRUE, color = NULL, color_list = NULL, BoolCol = FALSE, sizePoint = NULL){ # Get vizualization with or without labels for download
  if(!is.null(color_list[[var]]$color) && BoolCol == TRUE){
    color <- unlist(color_list[[var]]$color)
  }
  DimPlot(obj,group.by = var, label = dlabel, label.size = 6, cols = color, pt.size = sizePoint)
}



MakeUMAPhighlight <- function( obj, highliht_cells , sizePoint = 0.5, sizeHighlight = 0.25){# Create UMAP highlight for keeping an eye on what will be subsetted, the size of the point can be change. 
  DimPlot(obj, cells.highlight = highliht_cells, order = T,pt.size = sizePoint,sizes.highlight = sizeHighlight)+scale_color_manual(labels = c("Not kept","Kept"), values = c("gray","red"))+theme(legend.text = element_text(size = 10))
}

MakeUMAP_Red <- function( obj, highliht_cells , sizePoint = 0.5, sizeHighlight = 0.25){# same as before but for the case when all the datasets as been kept. 
  DimPlot(obj, cells.highlight = highliht_cells, order = T,pt.size = sizePoint,sizes.highlight = sizeHighlight)+scale_color_manual(labels = c("Kept"), values = c("red"))+theme(legend.text = element_text(size = 10))
}


MakeComplicatedVln <- function(obj, Giventitle, var, minGenePerCell, maxGenePerCell, maxCountPerCells,percentMito){ #Violin for the first graph on QC page and the two first graph of the filtering page.  
  nfeat <- VlnPlot(obj,group.by = var, features = c("nFeature_RNA"), pt.size = 0)+theme(legend.position ="none")+geom_hline(yintercept = c(minGenePerCell,maxGenePerCell), color = c("red", "red"))
  ncount <- VlnPlot(obj,group.by = var, features = c("nCount_RNA"), pt.size = 0)+theme(legend.position ="none")+geom_hline(yintercept = maxCountPerCells, color = "red ")
  pmt <- VlnPlot(obj,group.by = var, features = c("percent.mt"), pt.size = 0)+theme(legend.position ="none")+geom_hline(yintercept = percentMito, color = "red")
  grid.arrange(nfeat,ncount,pmt,nrow=1, top = textGrob(Giventitle, gp =gpar(col="black", fontsize =20, fontface ="bold"))) ## arrange all the plot between them
}

# Plotnbcellsbeforeafter <- function(objBefore, objAfter){#Get the number of cells before and after filtering of the data
#   dfNbcells <- data.frame(time = c("before_filtering", "after_filtering"),nb_cells = c(dim(objBefore)[2],dim(objAfter)[2]))
#   ggplot(dfNbcells,aes(x= time, y = nb_cells, fill = time))+geom_bar(stat = "identity")+scale_x_discrete(limits=c("before_filtering", "after_filtering"))+ggtitle("Impact of filtering on the number of cells")+geom_text(aes(label=nb_cells), position=position_dodge(width=0.9), vjust=-0.25)+theme_cowplot()
# }


### Encapsulation of findMarkers function for the Differential expression page, all the parameter will be passed to this function. 
RunFindMarkers <- function(obj, cluster1, cluster2 = NULL, posMarkers, minimum.percent, logthr){
  if(posMarkers == "Yes"){
    resultPos <- TRUE
  }else{
    resultPos <- FALSE
  }
  if(cluster2 == "All"){
    FindMarkers(object = obj, ident.1 = cluster1 ,only.pos = resultPos ,min.pct = minimum.percent, logfc.threshold = logthr)
  }else{
    FindMarkers(object = obj, ident.1 = cluster1, ident.2 = cluster2 ,only.pos = resultPos ,min.pct = minimum.percent, logfc.threshold = logthr)
  }
  
}



Readfile <- function(obj, filePath){
    validate(
      need(filePath, " ")
    )
    list_gene_file <- read.csv(filePath$datapath, header=TRUE)
    colnames(list_gene_file) <- "cluster"
    list_gene_file <- list_gene_file$cluster
    list_gene_file <- list_gene_file[list_gene_file %in% rownames(obj)]
    return(list_gene_file)
}

