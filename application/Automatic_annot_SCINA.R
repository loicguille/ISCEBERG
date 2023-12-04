#####################################################################################################################################################################################################
############################################################################## AUTOMATIC ANNOTATION #################################################################################################
#####################################################################################################################################################################################################

#object = seurat object, list_annotation = marker list, allowUnknown = boolean, allowOverlap = Boolean, group = annotation/resolution, colors = list of colors associated to each label
#scina = results of scina, slotName = slot of annotation 

source("New_SCINA.R")
useSCINA <- function(object,list_annotation, allowUnknown = TRUE , allowOverlap = TRUE){
  if(object@active.assay == "RNA"){
    scina <- SCINA_New(object@assays$RNA@counts[,],list_annotation, allow_unknown = allowUnknown, rm_overlap = allowOverlap)
  }else{
    scina <- SCINA_New(object@assays$SCT@counts[,],list_annotation, allow_unknown = allowUnknown, rm_overlap = allowOverlap)
  }
  return(scina)
}


Dimplot_with_ggplot <- function(object, group, CellsSize, colors = NULL){
  coord_Umap <- Embeddings(object[["umap"]])[,1:2]
  Proba_plot <- data.frame(coord_Umap,object[[group]])
  ggplot(Proba_plot, aes(UMAP_1,UMAP_2, color = Proba_plot[[group]]))+geom_point(size = CellsSize)+theme_cowplot()+scale_color_manual(name = group ,values = c(hue_pal()(length(levels(factor(Proba_plot[[group]]))))))
}

FeaturePlot_proba <- function(object, proba_col,slotNameAnnot, CellsSize, colors = NULL, BoolColScales = FALSE){
  validate(
    need(proba_col != "", "Wait for the automatic annotation")
  )
  coord_Umap <- Embeddings(object[["umap"]])[,1:2]
  if(proba_col == "All_probabilities"){
    Proba_plot <- data.frame(coord_Umap,object[[slotNameAnnot]], proba = object$All_probabilities)
    ggplot(Proba_plot, aes(UMAP_1,UMAP_2, color = Proba_plot[[slotNameAnnot]]))+geom_point(aes(alpha = proba), size = CellsSize)+theme_cowplot()+scale_color_manual(name = slotNameAnnot, values = c(hue_pal()(length(levels(factor(Proba_plot[[slotNameAnnot]])))-1),"lightgrey"))
  }else{
    prob_test <- data.frame(coord_Umap, object[[proba_col]])
    ggplot(prob_test, aes(UMAP_1,UMAP_2, color = prob_test[[proba_col]]))+geom_point(size = CellsSize)+theme_cowplot()+scale_color_gradientn(name = proba_col ,colors = c("lightgrey","blue"), limits = c(min(prob_test[[proba_col]]),max(prob_test[[proba_col]])))
  }
}

Fill_proba <- function(object, scina, slotName){
  for(i in rownames(scina$probabilities)){
    object[[paste(i,"probabilities", sep = "_")]] <- scina$probabilities[i,]
  }
  object[[slotName]] <- scina$cell_labels
  object$All_probabilities <- ifelse(scina$probabilities[1,] > scina$probabilities[2,],scina$probabilities[1,],scina$probabilities[2,])
  object$All_probabilities[which(object[[slotName]] == "unknown")] <- NA
  return(object)
}
