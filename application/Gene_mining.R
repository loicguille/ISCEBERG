#################################################################################################################################################################################################
############################################################################### DATA MINING PAGE FUNCTION #######################################################################################
#################################################################################################################################################################################################
#Create featurePlot 
FP <- function(object, gene, color = NULL, BoolColScales = FALSE){ # function that allow the construction of featureplot recursively
  validate(
    need(gene != "","Select at least one gene")
  )
  if(!is.null(color) && BoolColScales == TRUE){
    colorList <- c("lightgrey", color)
    FeaturePlot(object, features = gene, cols = colorList)
  }else{
    FeaturePlot(object, features = gene)
  }
}


#Create heatmap from genes
Heatmap <- function(obj , genes, color = NULL, color_list = NULL, var = NULL, BoolCol = FALSE, BoolColScales = FALSE, colorScale = NULL){
  validate(
    need(genes != "","Choose at least one gene and it will show you the heatmap of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
  )
  if(!is.null(color_list[[var]]$color) && BoolCol == TRUE){
    color <- unlist(color_list[[var]]$color)
  }
  if(BoolColScales == TRUE){
    DoHeatmap(object = obj, features=genes, group.colors = color,slot = "data")+scale_fill_gradient(low = 'lightgrey', high= colorScale)
  }else{
    DoHeatmap(object = obj, features=genes, group.colors = color,slot = "data")+scale_fill_gradient(low = 'lightgrey', high= "blue")
  }
  
}


#Create Dotplot from genes
Dotplot <- function(obj, genes, colorScale = NULL, BoolColScales = FALSE){
  validate(
    need(genes != "","Choose at least one gene and it will show you the dotplot of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
  )
  if(BoolColScales == TRUE){
    color <- c("lightgrey",colorScale)
    DotPlot(obj, features=genes, cols = color)
  }else{
    DotPlot(obj, features=genes)
  }
  
}

#Get the violin plot from when they are splited in data mining for one gene page.
violinPlotSplited <- function(seuratObj,genes,var1,var2, conditionVar2, conditionVar1, color = NULL, color_list = NULL){  
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
  if(!is.null(color_list[[var2]]$color)){
    color <- unlist(color_list[[var2]]$color)
    graph <- ggplot(filtered_table)+geom_violin(aes(variable1,gene, fill = variable2), scale = "width") + xlab("Identity")+ylab("Expression Level") +guides(fill = guide_legend(title=var2))+theme_classic()+ggtitle(genes)+theme(plot.title = element_text(face="bold",hjust = 0.5))+scale_fill_manual(values = color)
  }else{
    graph <- ggplot(filtered_table)+geom_violin(aes(variable1,gene, fill = variable2), scale = "width") + xlab("Identity")+ylab("Expression Level") +guides(fill = guide_legend(title=var2))+theme_classic()+ggtitle(genes)+theme(plot.title = element_text(face="bold",hjust = 0.5))
  }
  return(graph)
}

## encapsulation of downloadable plot for split object in rendering Violinplot
ViolinForDl <- function(seuratObject, geneList, var1, var2, conditionVar1, conditionVar2, color_list = NULL ){
  Violins <- list()
  for(i in 1:length(geneList)){
    Violins[[i]] <- violinPlotSplited(seuratObj = seuratObject, genes = geneList[i],var1 = var1, var2 = var2, conditionVar1 = conditionVar1, conditionVar2 = conditionVar2, color_list = color_list)
  }
  do.call(grid.arrange,Violins)
}

#Simple Violin plot generator
SimpleViolinPlot <- function(object,gene, var , color = NULL, color_list = NULL, BoolCol= FALSE){ #Take in entry seurat object and one gene, function used in data mining with one gene
  validate(
    need(gene != "","Choose at least one gene and it will show you the violin plot of this gene in function of the different clusters, you can choose 10 different gene ids to plot")
  )
  if(!is.null(color_list[[var]]$color) &&  BoolCol == TRUE){
    color <- unlist(color_list[[var]]$color)
  }
  VlnPlot(object, group.by = var,features = gene, pt.size = 0, cols = color)
}