# ISCEBERG 

## Pre-processing

Several options can be choosed : 

- Mtx : for a matrix file you will have to pass at least one matrix, one feature file , and one barcode file. You can give more of these files but the number has to be the same between matrix, feature file, and barcode file.

- CSV or Txt : The csv and TXT file have the same behaviour you will have to choose the separator and if your file has a header or not. If you have metadata you can give it and it will be added into your vizualisation. If your dataset cannot be automatically analyzed please contact us. 

- H5 : the H5 format is the major output for single cell and you will just have to pass file(s) into this format to create an output. You can also add metadata.

- RDS : if you choose RDS format you have to pass preprocessed data with seurat and save into an rds format. So the object opened has to be a seurat object otherwise it will be impossible to open it.

Once you have choosen the type of file, you can change the values we used and recommend by default for the preprocessing for values for the minimal number of time that a gene has to be expressed to be keeped and the number of genes a cells has to expressed to be keeped 

## Filtering page (if you have choose H5, CSV, MTX or TXT)

Once your seurat object is created you have to filter it with some values, in general we apply filter on the percent max of gene mitotic by cells (which is fixed by default to 25%), Seurat also recommend to filter on the cells that expressed a low number of genes which can be uninformative, the values of this filter is fixed by default to 500 genes by cells minimum but can be easily change in function of your datasets. Seurat also consider that a maximum number of genes expressed by cells have to be apply on the data because cell doublets or multiplets may exhibit an aberrantly high gene count.

A good indication for the filter that you can use are the three plots that are present on the right of the filter. The nFeature_RNA correspond to the number of genes expressed by cells for the moment so you can adapt the filter. The percent.mt correspond to the percent of mitotic gene in by cells. 

After that you will have to choose the method of normalization that you want to apply on your data : 

- SCTransform is the method that Seurat recommend to use because , this method take much time to run than lognormalize and can fail if you have too much cells in your dataset

- LogNormalize : Can be used if you have a very high number of cells. 

Then you will have to choose the resolution you want to make to do your clusterization first you have to choose a minimum and a maximum and then you will have to choose a step. Then you can click on "Run analysis". 

Once the analysis is done several QC plot will appear in order to check how your data are composed and if they don't have any problem: 

- The graphic concerning the impact of filtering on the number of cells is a good way to see how many cells have been excluded with your filter. If you think that you have filtered too much/few cells you can readapt your filter and click again on "Run analysis". 

- The Phase graph is a important to see if there is one step of the cell cycle that cluster our data.

- The number of expressed genes by cells is also a good way to see if the number of expressed genes by cells cluster our cells, normally the SCTransform normalization tend to disturb this bias. 

Some button will also appear once the filter have been run on the dataset, one to save the dataset you have filtered, and one to download the log file containing the version of the package used to analyze the data and the command with the filter you have made on the dataset. In order to download this two files you will have to provide a name.

## QC page (if you have choose available data or RDS)
 
This page all the user to check the same plot ads the results of the filtering page (Phase, nmber of expressed cells, clustering, ...)

## Cluster tree page

This page resent you how your cells have been split up in function of the resolution you have. In order to build this tree, the dataset have to contain at least 2 different resolution. a button will allow you to choose which resolution you want to represent on your UMAP. A higher resolution include more cluster. This implies that the composition in cells of the different cluster will change in function of the resolution. 
One other things to notice is that the cluster number is conversly correlated to the number of cells. 

## Differential expression between clusters page

In order to do differential expression between clusters :

1. Choose the resolution where you want to do the differential expression

2. Choose the first cluster

3. Choose the second cluster (if you choose All the differential will be made against all clusters grouped)

4. Choose if you want only the positive markers or not (in case of No the two ways will be made)

5. Choose a minimum value of threshold (Log2 scaled)

6. fill the minimal fraction of the cells that have to contains a gene to be tested (if you choose 0.2 at least 20% of the first cluster has to contain a gene to be tested)

7. Click on "Load"

Once the differential has been made a table will appear containing all the differentially expressed genes. This table contain the pvalue, the adjusted pvalue the average Log2 fold change and the percent of cell that expressed this genes in the first and the second cluster. You can apply filter on the results you have for each column and save it by clicking on "CSV".

## Data mining for one gene 

Data mining page will allow you to check different gene over different types of representation such as feature plot (a projection of the cells that expressed the gene) the violin plot (compute expression level for each cluster), a heatmap (each bars correspond to a cell and the cell are grouped by cluster) and a dotplot. 

10 genes maximum can be put together on the same scale or not (for featurePlot). Two different format for downloading the plot are available png and svg, svg graphs are heavier than png but they are also much more formattable for publication. The table show the mean expression of the genes you have choosen by cluster.

## Data mining for a combination of genes

The idea here is to combine the expression of every genes that is given as entry. Three different methods can be used in order to pool the genes together :
- the sum where all the expression counts are summerized for each cells and projecting on the different plot available. 

- the mean where the expression of genes is mean for each cells and projecting on the different plot available.

- the addmodulescore function calculate the average expression levels of each program (cluster) on single cell level, subtracted by the aggregated expression of control feature sets. All analyzed features are binned based on averaged expression, and the control features are randomly selected from each bin.

Three different method of vizualization are also available, featurePlot, violin Plot or Dotplot which can be save either in SVG or in PNG. Two different ways to fill list, you can pass gene by gene or by giving a file as entry. 

## Information extraction

Information extraction is used to know the composition of the cluster in some variable contain in the object used. For example if your dataset is composed of cells from condition control and KO you can have the number or percent of KO/control cells for each cluster. The plot are interative and can be downloaded into png format. 

## Add annotation page 

Once you have identify a cluster of interest based on the differential expressed genes or something, the browser will allow you to put an annotation on this cluster. In order to do so you can choose *Create*, and put a name **without space**, once you have done that you have to choose the resolution based on the cluster you want to annotate. Finally you can write your annotation and click on Annotate. If you have correctly done all the previous step a pop-up will be display telling you that you have correctly annotated your object, and the projection of what you have annotated. 
If you want to annotate a metadata that is already present in your dataset, you can choose the option to *Update* a previous annotation. This will display the umap projection of the annotation you want to annotate. Then you have to fill the form with the same informations as before and click on annotate. 
Once the annotation has been down you have the possibility to download the seurat object. This object will contains all the annotation you have done before. 

## Subclustering page 

if you are not interested in all the object that you have loaded but just about some conditions/clusters, this page will allow you to subclusterize your data. If you choose to subcluster from cluster you will have to give the resolution you want and the cluster number you want to extract (the cluster you want to extract will be in red in the projection on the right). If you choose annotation you can choose from which annotation and the label you want to get.  

The subclustering will calculate de novo a clustering with multiple choices of resolution and when this is over you will have a projection of the cluster as well as the possibility to download the seurat object (this can be long, don't push the button several time !!!) and the log file that recapitulate all the operations you have done before.

