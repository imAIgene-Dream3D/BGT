# Behavioral Guided Transcriptomics (BGT) Pipeline

This README is part of the extended BGT pipeline documentation, which builds upon the initial [BEHAV3D repository](https://github.com/AlievaRios/BEHAV3D). BEHAV3D is a multispectral imaging and analytical pipeline to 1) classify T cell based on behavior; 2) study tumor/organoid death dynamics; 3) backprojecting identified behavior in the imaging dataset. This pipeline runs with the output of [BEHAV3D T cell classification module](https://github.com/AlievaRios/BEHAV3D?tab=readme-ov-file#2-t-cell-behavior-classification-module)  and includes additional features able to integrate this phenotypic T cell behavioral information wiht single cell transcriptomics of T cells. 

## Overview

The BGT introduces new functionalities that integrates imaging and transcriptomic data, allowing for a more comprehensive analysis of the molecular drivers of immune cells behavior upon tumor organoid targeting. 

## Input data

### Dataset of classified T cells based on behavior. 
This data is obtained by performing multispectral time-lapse imaging of T cells incubated with tumor organoids and processed with BEHAV3D. Refer to [BEHAV3D](https://github.com/AlievaRios/BEHAV3D) and our detailed Protocol to obtain this data[link to NaPo BEHAV3D] 
### single cell sequencing data from T cells incubated with tumor organoids
This data is obtained from running SORT seq on T cells incubated with tumor organoids.

## Scripts

[Timepoint graph generation](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Timepoint_graph.R)
To determine the optimal time-point for the separation of highly engaging T cells, execute a dynamic graphical analysis focusing on a specified cluster of interest. This cluster is indicative of highly engaged T cells, often referred to as super-engagers. This analysis marks the time window during which the super-engager population is most discernible, aiding in the accurate prediction of the timepoint for T cell separation. The classified_tcell_track_data.rds file, encompassing the classified behavioral data for T cells, alongside the determined starting point of imaging (imaging time), are critical inputs for this analysis. These parameters should be incorporated into the provided BGT_config template within the simulation settings section. It is imperative that these are accurate and representative of your specific experimental setup.

[Population_seperation simulation](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Mutiple_timepoint_population_seperation_simulation.R)
Utilize population_seperation_simulation for analyzing the frequency and distribution of T-cell behavior over time in simulated co-culture environment. L

### scRNA-seq Analysis Pipeline
To use the scRNA-seq analysis pipeline you have to download the sequencing data deposited in the GEO depository with the accession number of (ADD ACCESSION NUMBER). Afterwards, proceed to run the following scripts. 
Please note that you either need to run all parts of the analysis one after the previous one, or use the intermediate outputs of the previous part of the pipeline by downloading them from Zenodo depository (ADD LINK TO ZENODO).
1. [Pre-processing of raw SORTseq data (QC, filtering)](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part1_SORTseq_QC.ipynb)
2. [Subset Identification](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part2_Subset_Identification.ipynb)
3. [Pseudotime Analysis](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part3_Pseudotime.ipynb)
4. [Dynamic Gene Clustering](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part4_Dynamic_Genes.ipynb)
5. [Comparison to In vivo Data](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part5_In_vivo_Comparison.ipynb)

[Behavioral Guided Transcriptomics](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Behavioral-guided_transcriptomics.R)
we perform integrative analysis of scRNA -seq Data with Behavioral Clusters derived from Imaging. Here we utilize myriads of  computational methods to cluster scRNA-seq data, annotate these clusters with behaviorally relevant data, and infer the probabilistic behavioral landscape of individual immune ccells. The goal is to correlate gene expression profiles with phenotypic behavioral patterns to unravel cellular mechanisms in a spatial and temporal context.

### What additional types of data does the extended BEHAV3D work with?
1. Classified T cell Track Data (.rds) file
2. CD8/CD4 Engagement frequency (.csv) file
3. single cell Sequencing data

- [//]: # (Commented instructions: Describe any new data types that the extended pipeline can process, such as additional imaging modalities or integration with other types of 'omics' data.)

### What additional outputs does the extended BEHAV3D provide?
1. Classified T cell Track Data (.rds) file
- [//]: # (Commented instructions: List new outputs that the extended pipeline will generate, such as advanced visualizations, integration with genomic data, etc.)


## Extended Software and Hardware Requirements
As the BGT pipeline is an extended implementation to the BEHAV3D protocol, the input data required to apply it is generated from BEHAV3D pipeline, which has been extensively described [REF NP BEHAV3D]. 
BEHAV3D demands high-end workstations (RAM 64GB, CPU 3.7 HZ or more with at least 8 cores, GPU with 16GB VRAM or more) for data processing with Imaris ensuring analysis efficiency. For applying the BGT data analysis pipeline, a general workstation will be sufficient. 
In terms of user expertise, a basic understanding of R programming is required to run the data analysis pipeline.

- [//]: # (Commented instructions: Mention any new software, hardware, or other system requirements that are specific to the extended pipeline.)

## Installations

BGT uses the following R libraries (R version 4.3.2 (2023-10-31)) :

| Package        | Version |
|----------------|---------|
| abind          | 1.4-5   |
| dplyr          | 1.1.4   |
| dtwclust       | 5.5.12  |
| fs             | 1.6.3   |
| future         | 1.33.1  |
| furrr          | 0.3.1   |
| ggplot2        | 3.4.4   |
| gplots         | 3.1.3   |
| MESS           | 0.5.9   |
| optparse       | 1.7.3   |
| parallel       | 4.3.0   |
| patchwork      | 1.2.0   |
| pheatmap       | 1.0.12  |
| plyr           | 1.8.8   |
| randomForest   | 4.7-1.1 |
| readr          | 2.1.4   |
| reshape2       | 1.4.4   |
| scales         | 1.3.0   |
| Seurat         | 5.0.1   |
| SeuratObject   | 5.0.1   |
| spatstat       | 3.0-6   |
| sp             | 1.6-1   |
| stats          | 4.3.0   |
| tibble         | 3.2.1   |
| tidyr          | 1.3.0   |
| tidyverse      | 2.0.0   |
| umap           | 0.2.10.0|
| viridis        | 0.6.4   |
| viridisLite    | 0.4.2   |
| xlsx           | 0.6.5   |
| yaml           | 2.3.8   |
| zoo            | 1.8-12  |


## Example Datasets for Extended Features

In this repository we provide example datasets consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids from the BEHAV3D original paper. Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: Either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).

- [//]: # (Commented instructions: If applicable, provide links to example datasets that are specifically used for demonstrating the new features of the extended pipeline.)

## Repository Structure 
This repository contains a collection of scripts and example datasets enabling the following dowstream analysis. Follow the structure in the script folder for each module and each analysis type. Introduce the corresponding folder/ file direction on your own computer where required (note that to specify directory paths in R (/) forward slash is recommended):

- [//]: # (Commented instructions: Describe the organization of the repository, including where users can find scripts, configurations, and data related to the new features.)

## Set-up 
BEHAV3D uses 2 specific fiels to customize the analysis:\

BEHAV3D config
Contains all experiment-specific settings and paths to data for all modules in BGT
An example version can be found in ...BGT/configs/config_template.yml
Explanation on what each variable changes is commented in that template

- [//]: # (Commented instructions: Detail any additional setup or configuration steps required to use the new features.)

## Demos 
You can run BGT on demo data to see examples of the results. This will take <10 minutes

There are 2 demos:

T cell population seperation simulation (For 'in-silico T cell and Tumor organoid co-culture)
Behavioral guided transcriptomics (For 'Behavioral integration with single cell Sequencing)
>Step 1 To set up the demo on you local PC, run ...BGT/demos/xyz.R
This sets up the paths in the BEHAV3D config file for the demo, then run the different modules on the demo (look below).

- [//]: # (Commented instructions: If you have created demos to showcase the new features, provide instructions on how users can run these demos.)

## Extended Modules

### [Identification of Seperation Time]

- [//]: # (Commented instructions: Provide a brief description of each new module, how to run it, and what outputs it generates.)

### [Population Seperation Simulation]

- [//]: # (Commented instructions: Provide a brief description of each new module, how to run it, and what outputs it generates.)

### [scRNAseq data preprocessing]

- [//]: # (Commented instructions: Provide a brief description of each new module, how to run it, and what outputs it generates. @Peter, put here your part and links to code uploaded in the scripts folder)

### [Pseudotime trajectory inference]

- [//]: # (Commented instructions: Provide a brief description of each new module, how to run it, and what outputs it generates. @Farid, put here your part and links to code uploaded in the scripts folder)

### [Behavioral Guided Transcriptomics]

- [//]: # (Commented instructions: Provide a brief description of each new module, how to run it, and what outputs it generates. @Farid, put here your part and links to code uploaded in the scripts folder)


## Further Reading and Documentation
### Online UCSC Browser
Link to [online UCSC browser](https://cells.ucsc.edu/?ds=behavior-guided-tx) to see how to use it go to [this page](https://cellbrowser.readthedocs.io/en/master/interface.html) 
<img src="https://cellbrowser.readthedocs.io/en/master/_images/cellbrowser-BasicUiFeatures.converted.jpg" width="70%">


Currently you can view it by :
- Cell
- Oganoid
- Experiment
- Plate
- Pseudotime
- Experimental Condition
- PT Cluster
- [//]: # (Commented instructions: If there's further documentation, such as a wiki, or publications related to the extended features, provide links here.)


