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

[Timepoint graph generation](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Timepoint_graph.R), [Multiple_timepoint_simulation_for population_seperation](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Mutiple_timepoint_population_seperation_simulation.R)

### scRNA-seq Analysis Pipeline
To use the scRNA-seq analysis pipeline you have to download the sequencing data deposited in the GEO depository with the accession number of (ADD ACCESSION NUMBER). Afterwards, proceed to run the following scripts. 
Please note that you either need to run all parts of the analysis one after the previous one, or use the intermediate outputs of the previous part of the pipeline by downloading them from Zenodo depository (ADD LINK TO ZENODO).
1. [Pre-processing of raw SORTseq data (QC, filtering)](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part1_SORTseq_QC.ipynb)
2. [Subset Identification](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part2_Subset_Identification.ipynb)
3. [Pseudotime Analysis](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part3_Pseudotime.ipynb)
4. [Dynamic Gene Clustering](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part4_Dynamic_Genes.ipynb)
5. [Comparison to In vivo Data](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part5_In_vivo_Comparison.ipynb)

## Additional Data Types and Outputs

### What additional types of data does the extended BEHAV3D work with?

- [//]: # (Commented instructions: Describe any new data types that the extended pipeline can process, such as additional imaging modalities or integration with other types of 'omics' data.)

### What additional outputs does the extended BEHAV3D provide?

- [//]: # (Commented instructions: List new outputs that the extended pipeline will generate, such as advanced visualizations, integration with genomic data, etc.)

## Citing the Extended Pipeline

- [//]: # (Commented instructions: Provide a reference format for users to cite the extended pipeline if it has been published or is available as a preprint.)

## Extended Software and Hardware Requirements

- [//]: # (Commented instructions: Mention any new software, hardware, or other system requirements that are specific to the extended pipeline.)

## Installation of New Features

- [//]: # (Commented instructions: Outline the steps for installing new components of the extended pipeline, including any additional libraries or dependencies.)

## Input Data for Extended Features

- [//]: # (Commented instructions: Explain how users should prepare and organize their data to make use of the new features in the extended pipeline.)

## Example Datasets for Extended Features

- [//]: # (Commented instructions: If applicable, provide links to example datasets that are specifically used for demonstrating the new features of the extended pipeline.)

## Repository Structure for Extended Features

- [//]: # (Commented instructions: Describe the organization of the repository, including where users can find scripts, configurations, and data related to the new features.)

## Set-up for Extended Features

- [//]: # (Commented instructions: Detail any additional setup or configuration steps required to use the new features.)

## Demos for Extended Features

- [//]: # (Commented instructions: If you have created demos to showcase the new features, provide instructions on how users can run these demos.)

## Extended Modules

### [New Module Name]

- [//]: # (Commented instructions: Provide a brief description of each new module, how to run it, and what outputs it generates.)

### [scRNAseq data preprocessing]

- [//]: # (Commented instructions: Provide a brief description of each new module, how to run it, and what outputs it generates. @Peter, put here your part and links to code uploaded in the scripts folder)

### [Pseudotime trajectory inference]

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


