# ISCEBERG
## Interactive Single Cell Expression Browser for Exploration of RNAseq data using Graphics 

<img src=application/www/iceberg2.png height="200">

## Overview of ISCEBERG

Here we describe our single cell data explorer. This shiny application has been developped to analyze, vizualize and extract informations from single-cell sequencing datasets. ISCEBERG can transform raw counts into filtered, normalized counts and computes clustering, UMAP to provide a Seurat single-cell sequencing dataset from various input formats (H5, MTX, TXT or CSV). ISCEBERG allows the use of preprocessed RDS files as input containing already a Seurat single-cell sequencing dataset. The purpose of this application is to explore much deeper your single-cell datasets without using R code lines. Graphics and associated tables are downloadable. To ensure Reproducibility and Traceability, preprocessed and subclusterized data are provided as Seurat objects (to be re-analyzed or re-loaded later) with a report of all command lines used to generate them. 

## Citation

## Details 

All the functionnalities are listed here. For more details see the Documentation section :

- Pre-processing (read, create and apply a first filtering on your data)
- Filtering (only if you have chosen  Mtx, H5, Txt, or Csv, apply a filtering on your dataset about number of genes expressed per cell, mitochondrial DNA percentage)
- QC (QC plots about cells distribution among studied samples, projection of metadata on UMAP, cell cycle phase score computing)
- Cluster tree (create a cluster tree from different resolutions present in the object)
- Differential Expression between clusters
- Data mining (expression plot  of one or several genes with different vizualisation methods)
- Data mining for a combination of genes (expression plot of a list of genes given in parameter or with a file)
- Extract information (ex : distribution of metadata groups across clusters)
- Add annotation 
- Subclustering (subcluster data based on annotation or cluster)

## Installation 

Here is the procedure to install our application :

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
