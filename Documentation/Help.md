# ISCEBERG 

<img src= ../application/www/iceberg2.png height="200">

## Pre-processing

Several options can be chosen : 

- Mtx : for a matrix file you will have to include at least one matrix, one feature file , and one barcode file. You can upload more than one of these files but the number has to be the same between matrix, feature, and barcode files.

- CSV or Txt : The CSV and TXT file have the same behaviour. You will have to choose the separator and whether your file has a header or not. If you have a  metadata file you can include it and it will be added into your vizualisation. If your dataset cannot be automatically analyzed please contact us. 

- H5 : the H5 format is the major output for single-cell sequencing data and you will just have to load file(s) in this format to create an output. You can also add metadata.

- RDS : if you choose RDS format you have to load preprocessed data with Seurat and save in RDS format. The object opened has to be a Seurat object otherwise it will be impossible to open it.

Once you have chosen the type of file, you can change the values we used and recommend by default for the preprocessing concerning the minimal number of time that a gene has to be expressed to be kept and the number of genes a cells has to express to be kept. 

## Filtering page (if you have chosen H5, CSV, MTX or TXT)

Once your Seurat object is created you have to filter on several values, in general we apply a filter on the maximum percent of mitotic  genes expressed by  cells (which is fixed by default to 25%). Seurat also recommend to filter on the cells that expressed a low number of genes which can be uninformative, the values of this filter is fixed by default to 500 genes by cells minimum, but can be easily changed in function of your datasets. Seurat also considers that a maximum number of genes expressed by cells has to be applied on the data because cell doublets or multiplets may exhibit an aberrantly high gene count.

A good indication for the filter that you can use are the three plots that are present on the right of the filter. The nFeature_RNA corresponds to the number of genes expressed by cells in the current state, so you can adapt the filter. The percent.mt corresponds to the percent of mitotic genes expressed by cells.

After that you will have to choose the method of normalization that you want to apply on your data : 

- SCTransform is the method that Seurat recommends to use, but this method takes much time to run than lognormalize and can fail if you have too much cells in your dataset (recommended method).

- LogNormalize : Can be used if you have a very high number of cells as an alternative. 

Then you will have to choose the resolution you want to achieve with your clusterization. First you have to choose a minimum and a maximum and then you will have to choose a step. Then you can click on "Run analysis". 

Once the analysis is done, several QC plot will appear in order to check how your data are composed and if they don't have any problem: 

- The graphic concerning the impact of filtering on the number of cells is a good way to see how many cells have been excluded with your filter. If you think that you have filtered too much/few cells,  you can readapt your filter and click again on "Run analysis". 

- The Phase graph is important to see if there is one step of the cell cycle that determines your clustering.

- The number of expressed genes by cell is also a good way to see if the number of expressed genes by cell determines or biases your clustering. Normally the SCTransform normalization tends to avoid this bias. 

Several buttons will also appear once the filters have been run on the dataset: one to save the dataset you have filtered, and one to download the log file containing the version of the package used to analyze the data and the command with the filter you have applied on the dataset. In order to download this two files, you will have to provide a name.

## QC page (if you have chosen RDS)
 
This page allows the user to check the same plot as in the results of the filtering page (Phase, number of expressed cells, clustering, ...)

## Cluster tree page

This page presents you how your cells have been split up in function of the resolution you have chosen. In order to build this tree, the dataset has to contain at least 2 different resolutions. A button will allow you to choose which resolution you want to represent on your UMAP. A higher resolution includes more cluster. This implies that the cell composition of the different cluster will change in function of the resolution. 
Another things to notice is that the cluster number is conversly correlated to the number of cells. 

## Differential expression between clusters page

In order to do differential expression between clusters :

1. Choose the resolution at which you want to do the differential expression.

2. Choose the first cluster.

3. Choose the second cluster (if you choose "All" the differential will be made against all clusters grouped, i.e. their mean).

4. Choose if you want only the positive markers or not (in case of "No" both positive and negative markers determined).

5. Choose a minimum value of threshold (Log2 scaled) for comparison fold expression.

6. Fill the minimal fraction of the cells that have to contain a gene to be tested (if you choose 0.2 at least 20% of the first cluster has to contain the gene to be tested).

7. Click on "Load".

Once the differential has been made, a table will appear containing all the differentially expressed genes. This table contain the p-value, the adjusted p-value the average Log2 fold change and the percent of cells that expressed this genes in the first and the second cluster. You can apply filters on the results for each column and save the table by clicking on "CSV".

## Data mining for one gene 

The data mining page will allow you to visualize gene expression through different types of representations such as a feature plot (a projection of the cells that expressed the gene) a violin plot (compute expression level for each cluster), an heatmap (each bars correspond to a cell and the cell are grouped by cluster) and a dotplot. 

10 genes maximum can be visualized together on the same scale or not (option only for featurePlot). Two different formats for downloading the plot are available (PNG and SVG). SVG graphs are heavier than png but they are also much more scalable/formattable for publication. The table shows the mean expression of the genes you have choosen by cluster.

## Data mining for a combination of genes

The idea here is to combine the expression of a set of genes that is given as entry. Three different methods can be used in order to pool the genes together :
- the sum, where all the expression counts are summerized for each cell and projected on the different plots available. 

- the mean, where the mean of the expression of the selected set of genes is mean for each cells and projected on the different plot available.

- the addmodulescore function calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.

Three different methods of vizualization are also available, featurePlot, violin Plot or Dotplot which can be save either in SVG or in PNG. There are Two different ways to fill the query list, you can pass gene by gene or by providing a file as entry. 

## Information extraction

Information extraction is used to know the composition of the cluster with regard to certain variables contained in the object used. For example if your dataset is composed of cells from condition control and treated you can have the number or percent of treated/control cells for each cluster. The plots are interative and can be downloaded in PNG format. 

## Add annotation page 

Once you have identify a cluster of interest based on the differential expressed genes or else, the browser will allow you to add an annotation on this cluster. In order to do so, you can choose *Create*, and put a name **without space**, once you have done that you have to choose the resolution based on the cluster you want to annotate. Finally you can write your annotation and click on Annotate. If you have correctly done all the previous steps a pop-up will be displayed telling you that you have correctly annotated your object, and the projection of what you have annotated. 
If you want to annotate a metadata that is already present in your dataset, you can choose the option to *Update* a previous annotation. This will display the UMAP projection of the annotation you want to annotate. Then you have to fill the form with the same informations as before and click on "annotate". 
Once the annotation has been done you have the possibility to download the Seurat object. This object will contains all the annotation you have done before. 

## Subclustering page 

If you are not interested in the entire object that you have loaded, but only certain conditions/clusters, this page will allow you to subclusterize your data. If you choose to subcluster from cluster, you will have to give the resolution you want and the cluster number you want to extract (the cluster you want to extract will be in red in the projection on the right). If you choose "annotation" you can choose based on which annotation and the label you want to select cells.  

The subclustering will calculate de novo a clustering with multiple choices of resolution. When this is over, you will have a projection of the cluster as well as the possibility to download the Seurat object (this can be long, don't push the button several times or the browser might crash !!!) and the log file that recapitulates all the operations you have done before.

