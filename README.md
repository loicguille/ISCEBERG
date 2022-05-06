# ISCEBERG
## Interactive Single Cell Expression Browser for Exploration of RNAseq data using Graphics 

<img src=application/www/iceberg2.png height="200">

## Overview of ISCEBERG

Here is our single cell data explorer shiny,  this has been developped in order to analyze/vizualize and extract informations from single cell datasets. This application can transform raw counts into single cell dataset from different input format (H5, MTX, TXT or CSV) or allow you to give as input a preprocess RDS file that contain already a single cell dataset in format Seurat. The purpose of this application is to go much deeper into your single cell datasets without using code lines.

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
