# ISCEBERG
## Interactive Single Cell Expression Browser for Exploration of RNAseq data using Graphics 

<img src=application/www/iceberg.png height="200">

## Overview of ISCEBERG

Here is our single cell data explorer in shiny this has been build in order to analyse/vizualize and extract informations from single cell datasets. This application will allow you to analyze you proper single cell datasets (RDS, H5, MTX, TXT or CSV) or to vizualize one of the several dataset that have already been analyze by our team. The purpose of this application is to allow the user to download dataframe or plot or even subset of the dataset of interest.

## Citation

## Details 

All the functionnality are listed here, for more details see Documentation section :

- Pre-processing (read, create and apply a first filtering on your data)
- filtering (only if you have chosen  Mtx, H5, Txt, or Csv, apply a filtering on your dataset and show some QC plot)
- QC (only if you have chosen RDS, show some QC plot)
-  Cluster tree (create a cluster tree from different resolution in the object)
-  DE between cluster
-  Data mining (Plot expression of one or several genes with different vizualisation method)
-  Data mining for a combination of gene (Plot expression of a list of genes passed in parameter or with a file)
-  Extract information
-  Add annotation 
-  Subclustering (subcluster data based on annotation or cluster)

## Installation 

Here is the procedure in order to install application :

`gitclone ISCEBERG`

`cd ISCEBERG`

In order to create a docker image run the command

`docker build -t image_name .`(if you are already in the directory containing the dockerfile)

Once the image is create you can run

`docker run -p 3838:3838 image_name`

Then type the adress in your web browser :

`localhost:3838`

## Documentation

Help is available in Documentation directory in the document [Help.md](https://github.com/loicguille/ISCEBERG/blob/master/Documentation/Help.md). 
