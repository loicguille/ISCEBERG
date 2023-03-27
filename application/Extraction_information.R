### function for information extraction page
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

Plot2Render <- function(obj,cluster,StackedBarPlot,VariablePlot, NbCluster, freqOrValues, color = NULL, color_list = NULL, BoolCol = FALSE){
  if(!is.null(color_list[[VariablePlot]]$color) && BoolCol == TRUE){
    color <- unlist(color_list[[VariablePlot]]$color)
  }
  if(cluster == "All"){
    if(StackedBarPlot == "No"){
      table_tmp2 <-as.data.frame(table(c(obj[[VariablePlot]], obj[[NbCluster]])))
      colnames(table_tmp2) <- c("condition","cluster", "freq")
      test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
      if(!is.null(color)){
        g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = position_dodge())+scale_fill_manual(values = color)+theme_cowplot()
      }else{
        g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = position_dodge())+scale_fill_manual(values = test)+theme_cowplot()
      }
      
    }else{
      if(freqOrValues =="Frequence"){
        table_tmp2 <-as.data.frame(table(c(obj[[VariablePlot]], obj[[NbCluster]])))
        colnames(table_tmp2) <- c("condition","cluster", "freq")
        test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
        if(!is.null(color)){
          g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "fill")+scale_fill_manual(values = color)+theme_cowplot()
        }else{
          g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "fill")+scale_fill_manual(values = test)+theme_cowplot()
        }
      }else{
        table_tmp2 <-as.data.frame(table(c(obj[[VariablePlot]], obj[[NbCluster]])))
        colnames(table_tmp2) <- c("condition","cluster", "freq")
        test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$condition)))
        if(!is.null(color)){
          g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = color)+theme_cowplot()
        }else{
          g <- ggplot(table_tmp2, aes(x= cluster, freq, fill = condition))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = test)+theme_cowplot()
        }
      }
    }
  }else{
    table_tmp <- as.data.frame(table(c(obj[[VariablePlot]],obj[[NbCluster]])))
    colnames(table_tmp) <- c("condition","cluster", "freq")
    table_tmp <- table_tmp[which(table_tmp$cluster == cluster),]
    test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp$condition)))
    if(!is.null(color)){
      g <- ggplot(table_tmp, aes(x= condition, freq, fill = condition))+geom_bar(stat="identity")+scale_fill_manual(values = color)+theme_cowplot()
    }else{
      g <- ggplot(table_tmp, aes(x= condition, freq, fill = condition))+geom_bar(stat="identity")+scale_fill_manual(values = test)+theme_cowplot()
    }
  }
}


percentRender <- function(obj, VariablePlot, NbCluster, color = NULL, color_list = NULL, BoolCol = FALSE){
  if(!is.null(color_list[[NbCluster]]$color) && BoolCol == TRUE){
    color <- unlist(color_list[[NbCluster]]$color)
  }
  table_tmp2 <-as.data.frame(table(c(obj[[VariablePlot]], obj[[NbCluster]])))
  colnames(table_tmp2) <- c("condition","cluster", "freq")
  tmp <- as.data.frame(rep(table(obj[[VariablePlot]]),length(unique(obj[[NbCluster]]))))
  colnames(tmp)<- "total"
  table_tmp2 <- cbind(table_tmp2,tmp)
  table_tmp2$percent <- table_tmp2$freq/table_tmp2$total*100
  test <- colorRampPalette(brewer.pal(9,"Set1"))(length(unique(table_tmp2$cluster)))
  if(!is.null(color)){
    g <- ggplot(table_tmp2, aes(x= condition, percent, fill = cluster))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = color)+theme(axis.text.x = element_text(angle =60,hjust=1))+theme_cowplot()
  }else{
    g <- ggplot(table_tmp2, aes(x= condition, percent, fill = cluster))+geom_bar(stat="identity", position = "stack")+scale_fill_manual(values = test)+theme(axis.text.x = element_text(angle =60,hjust=1))+theme_cowplot()
    
  }
}

Table2Render <- function(obj,VariablePlot, NbCluster, cluster){
  if(cluster == "All"){
    as.data.frame(table(c(obj[[VariablePlot]], obj[[NbCluster]])))
  }else{
    as.data.frame(table(obj[[VariablePlot]][which(obj[[NbCluster]]== cluster),]))
  }
}
