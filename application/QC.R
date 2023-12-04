# QC Plot function
library(Seurat)
library(Matrix)
library(bslib)
library(ggplot2)
library(DT)
library(cowplot)
library(RColorBrewer)
library(plotly)
library(grid)
library(dplyr)
library(gridExtra)


nbCellsbydt <- function(obj){#Get graph of number of cells by datasets
  cells_by_dt <- data.frame(table(obj$orig.ident))
  ggplot(cells_by_dt,aes(Var1,Freq, fill=Var1))+geom_bar(stat = "identity")+ggtitle("Number of cells by datasets")+geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25)+theme_cowplot()+theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))
}

# UmapNbGenesCell <- function(obj, sizePoint){#Get umap with number of genes by cells. 
#   coord_Umap <- Embeddings(obj[["umap"]])[,1:2]
#   expressed_cells <- cbind.data.frame(counts=colSums(obj@assays$RNA@data > 0),coord_Umap)
#   ggplot(expressed_cells, aes(UMAP_1,UMAP_2, color = counts))+geom_point(shape =20, size=sizePoint)+scale_color_gradient(low="lightgrey", high="blue")+theme_cowplot()+ggtitle("Number of expressed genes by cells")
# }

makeVlnGreatAgain <- function(obj,var, grouping, col = 1){ #Create violin plot for features and group by the annotation/resolution needed
  VlnPlot(object = obj, features = var, group.by = grouping, ncol = col, pt.size = 0 )
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

VlnWithThresholdAppearance <- function(obj, var, threshold){
  VlnPlot(obj, group.by = "file" ,features =  var, pt.size = 0)+theme(legend.position ="none")+geom_hline(yintercept = threshold, color = "red")
}

Plotnbcellsbeforeafter <- function(objBefore, objAfter){#Get the number of cells before and after filtering of the data
  dfNbcells <- data.frame(time = c("before_filtering", "after_filtering"),nb_cells = c(dim(objBefore)[2],dim(objAfter)[2]))
  ggplot(dfNbcells,aes(x= time, y = nb_cells, fill = time))+geom_bar(stat = "identity")+scale_x_discrete(limits=c("before_filtering", "after_filtering"))+ggtitle("Impact of filtering on the number of cells")+geom_text(aes(label=nb_cells), position=position_dodge(width=0.9), vjust=-0.25)+theme_cowplot()
}

DataTableGenebyResorAnnot <- function(obj, metadata){
  expressed_cells <- cbind.data.frame(counts=colSums(obj@assays$RNA@data > 0),obj[[metadata]])
  colnames(expressed_cells) <- c("counts","value")
  dataTable <- aggregate(expressed_cells$counts, list(expressed_cells$value), FUN=mean)
  colnames(dataTable) <- c(metadata,"mean_expressed_nb_gene")
  dataTable$mean_expressed_nb_gene <- round(dataTable$mean_expressed_nb_gene)
  dataTable <- as.data.frame(dataTable)
}
