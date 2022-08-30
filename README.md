[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6563734.svg)](https://doi.org/10.5281/zenodo.6563734)

# ISCEBERG
## Interactive Single Cell Expression Browser for Exploration of RNAseq data using Graphics 

<img src=application/www/iceberg2.png height="200">

## Overview of ISCEBERG

Here we describe our single cell data explorer. This Shiny application has been developped to analyze, visualize and extract informations from single-cell sequencing datasets. ISCEBERG can transform raw counts into filtered, normalized counts and computes clustering, UMAP projections to provide a Seurat single-cell sequencing dataset from various input formats (H5, MTX, TXT or CSV). ISCEBERG allows the use of preprocessed RDS files as input containing already a Seurat single-cell sequencing dataset. The purpose of this application is to explore much deeper and easily visualize your single-cell datasets without using R code lines. Graphics and associated tables are downloadable. To ensure Reproducibility and Traceability, preprocessed and subclusterized data are provided as Seurat objects (to be re-analyzed or re-loaded later) with a report of all command lines used to generate them. 

## Citation

Loïc Guille, Manuel Johanns, Francesco Zummo, Bart Staels, Philippe Lefebvre, Jérôme Eeckhoute, & Julie Dubois-Chevalier. (2022). ISCEBERG : Interactive Single Cell Expression Browser for Exploration of RNAseq data using Graphics (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.6563734

## Details 

All the functionnalities are listed here. For more details see the Documentation section :

- Pre-processing (read, create and apply a first filtering on your data)
- Filtering (only if you have chosen  Mtx, H5, Txt, or Csv, apply a filtering on your dataset on number of genes expressed per cell, mitochondrial DNA percentage)
- QC (QC plots about cells distribution among studied samples, projection of metadata on UMAP, cell cycle phase score computing)
- Cluster tree (create a cluster tree from different resolutions present in the object)
- Differential gene expression between clusters
- Data mining (expression plots  of one or several genes with different vizualisation methods)
- Data mining for a combination of genes (expression plots of a list of genes given in parameter or with a file)
- Extract information (ex : distribution of metadata groups across clusters)
- Add annotations
- Subclustering (subcluster data based on annotation or cluster)

## Installation 

Here is the procedure to install our application :

`gitclone https://github.com/loicguille/ISCEBERG.git`

`cd ISCEBERG`

In order to create a docker image run the command

`docker build -t image_name .`(if you are already in the directory containing the dockerfile)

Once the image has been created you can run

`docker run -p 3838:3838 image_name`

Then type the adress in your web browser :

`localhost:3838`

## Documentation

Help is available in Documentation directory in the document [Help.md](https://github.com/loicguille/ISCEBERG/blob/master/Documentation/Help.md). 

## Authors

Loïc Guille;  Manuel Johanns;  Francesco Zummo;  Bart Staels;  Philippe Lefebvre;  Jérôme Eeckhoute;  Julie Dubois-Chevalier

## Grants

This work was supported by the Agence Nationale de la Recherche (ANR) grants “HSCreg” (ANR-21-CE14-0032-01) , “European Genomic Institute for Diabetes” E.G.I.D (ANR-10-849 LABX-0046), a French State fund managed by ANR under the frame program Investissements d’Avenir I-SITE ULNE / ANR-16-IDEX-0004 ULNE,  by grants from the Fondation pour la Recherche Médicale (FRM : EQU202203014645) and by European Commission

Agence Nationale de la Recherche:
    
    • EGID - EGID Diabetes Pole (10-LABX-0046) 

European Commission:
    
    • ImmunoBile - Bile acid, immune-metabolism, lipid and glucose homeostasis (694717) 
