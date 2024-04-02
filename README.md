# Behavioral Guided Transcriptomics (BGT) Pipeline

BGT is an imaging-guided transcriptomic experimental and analytical pipeline for co-cultures of Patient Derived Organoides with (engeneered) T cells described in [Dekkers&Alieva et al Nat. Biotech (2022)](https://www.nature.com/articles/s41587-022-01397-w). This repository compiles the analytical procedures required to integrate data from single-cell RNA sequencing of T cells with single-cell imaging data, providing for each sequenced cell probabilities of exhibiting different functional behaviors.

## Input data

### Imaging-derived dataset of behavior-classified T cells. 
This data is obtained by performing multispectral time-lapse imaging of T cells incubated with tumor organoids and processed with BEHAV3D. Refer to [BEHAV3D](https://github.com/AlievaRios/BEHAV3D) repository and our detailed Protocol to obtain this data [Alieva et al, Nat Protoc, 2024](https://www.nature.com/articles/s41596-024-00972-6) $${\color{red}Avi}$$ what is the name of the file that is the output of BEHAV3D and input in BGT, has to be the same name
### Single cell sequencing data from T cells incubated with patient derived organoids.
This data is obtained from running SORT seq on T cells incubated with tumor organoids.
$${\color{red}Miguel}$$ What is the name of the input file for scRNA seq that we provide as a demo
The scRNA seq file we provide for the demo is 10T_master.rds 

## Extended Software and Hardware Requirements
As the BGT pipeline is an extended implementation to the BEHAV3D protocol, the input data required to apply it is generated from BEHAV3D pipeline, which has been extensively described [Alieva et al, Nat Protoc, 2024](https://www.nature.com/articles/s41596-024-00972-6) 
BEHAV3D demands high-end workstations (RAM 64GB, CPU 3.7 HZ or more with at least 8 cores, GPU with 16GB VRAM or more) for data processing with Imaris ensuring analysis efficiency. For applying the BGT data analysis pipeline, a general workstation will be sufficient. 
In terms of user expertise, a basic understanding of R programming is required to run the data analysis pipeline.

- [//]: # (Commented instructions: Mention any new software, hardware, or other system requirements that are specific to the extended pipeline.)

## Installations 
$${\color{red}Avi/and/Miguel}$$ we need to have an integrated version for installation of all the packages related to both parts. and ideally provide a docker image with the package for both parts
BGT uses the following R libraries (R version 4.3.2 (2023-10-31)) :

| Package        | Version  |
|----------------|----------|
| abind          | 1.4-5    |
| dplyr          | 1.1.4    |
| dtwclust       | 5.5.12   |
| fs             | 1.6.3    |
| future         | 1.33.1   |
| furrr          | 0.3.1    |
| ggplot2        | 3.4.4    |
| gplots         | 3.1.3    |
| MESS           | 0.5.9    |
| optparse       | 1.7.3    |
| parallel       | 4.3.0    |
| patchwork      | 1.2.0    |
| pheatmap       | 1.0.12   |
| plyr           | 1.8.8    |
| randomForest   | 4.7-1.1  |
| readr          | 2.1.4    |
| reshape2       | 1.4.4    |
| scales         | 1.3.0    |
| Seurat         | 5.0.1    |
| SeuratObject   | 5.0.1    |
| spatstat       | 3.0-6    |
| sp             | 1.6-1    |
| stats          | 4.3.0    |
| tibble         | 3.2.1    |
| tidyr          | 1.3.0    |
| tidyverse      | 2.0.0    |
| umap           | 0.2.10.0 |
| viridis        | 0.6.4    |
| viridisLite    | 0.4.2    |
| xlsx           | 0.6.5    |
| yaml           | 2.3.8    |
| zoo            | 1.8-12   |
| BiocManager    | 3.18     |
| ggthemes       | (latest) |
| gridExtra      | (latest) |
| openxlsx       | (latest) |
| devtools       | (latest) |
| monocle3       | (latest) |

**Note:** For packages listed with "(latest)" as their version, you should install the most recent version available at the time of installation. This can be done using the `install.packages()` function for CRAN packages and the `BiocManager::install()` function for Bioconductor packages, without specifying a version number.

## Demo and intermediate files 

In this repository we provide example datasets as well as intermediate files to run the different steps of the pipeline:
$${\color{red}Avi/and/Miguel}$$ 
- [//]: # (Commented instructions: list here the datasets that are provided in the demo folder including also the intermediary datasets that have.)

RNAseq
- 10T_master.rds
- All_TEGs_MNN_Seurat.rds
- CD4_IL17RBNeg_Dynamic_Genes.rds
- CD4_IL17RBNeg_MNN_Seurat.rds
- CD4_IL17RBPos_Dynamic_Genes.rds
- CD4_IL17RBPos_MNN_Seurat.rds
- CD8_Dynamic_Genes.rds
- CD8_MNN_Seurat.rds
- Savas_CD8_Trm_MarkerGenes.rds
intermediate demo objects to test the pipeline partially or full

## Repository Structure 
This repository contains a collection of scripts and example datasets enabling the following dowstream analysis. Follow the structure in the script folder for each module and each analysis type. Introduce the corresponding folder/ file direction on your own computer where required (note that to specify directory paths in R (/) forward slash is recommended):

## Set-up 
BGT uses 2 specific fiels to customize the analysis:\

BGT_config
Contains all experiment-specific settings and paths to data for all modules in BGT
An example version can be found in ...BGT/configs/BGT_config_template.yml
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

##  Modules
![BGT workflow](https://github.com/AlievaRios/BGT/blob/dev_avi/BGT%20analysis%20workflow.jpg)
### (1) Evaluation of super-engager population dynamics in co-culture

This module focuses on analyzing T-cell engagement dynamics, particularly highlighting the activity within cluster 9—representative of super-engagers. It enables a detailed comparison between CD4 and CD8 T-cells' engagement over time in co-culture experiments, visualizing the engagement percentage of T-cells in the super-engager state across various time points. [This analysis](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Timepoint_graph.R) marks the time window during which the super-engager population is most discernible, aiding in the accurate prediction of the timepoint for T cell separation. The classified_tcell_track_data.rds file, encompassing the classified behavioral data for T cells, alongside the determined starting point of imaging (imaging time), are critical inputs for this analysis. These parameters should be incorporated into the provided BGT_config template within the simulation settings section. It is imperative that these are accurate and representative of your specific experimental setup.

#### Configuration
Before running the module, ensure the config file (`config_template.yml`) is correctly set up with your specific parameters, including:
- `classified_tcell_track_data_filepath_rds`: Path to the classified T cell track data file.
- `imaging_time`: Time (in hours) when imaging starts in the experiment; this module will focus on time points after this.

***To run from the command line:***
```bash
Rscript /path/to/TimepointGraph/timepoint_graph.R -c </Path/to/BGT/config_template.yml> -f
```
The `-f` flag forces the re-import and processing of data even if the output files already exist. Optionally, use `-t` to specify an RDS file containing processed T-cell track data, if not already included in the BGT config.

***To run from RStudio:***

**Step 1:** For a demo run, use the [T cell behavioral dynamics](/scripts/TimepointGraph/timepoint_graph.R). Adjust the path to your **BGT config** on [line 18](/scripts/TimepointGraph/timepoint_graph.R#L18) if you're using a different data folder or config file.

***Output Files***

Outputs are saved in the specified `output_dir` within the `Results/Seperation_time_graph/` directory:
- **T-cell_engagementVsTime.png**: Visualization of the percentage of T cells in cluster 9 (super-engagers) over time.
- **T-cell_engagementVsTime.csv**: Raw data file containing the time points and percentages of T cells in cluster 9, ready for further analysis.

#### Key Features
- **Temporal Dynamics Insight**: The module reveals critical engagement dynamics, identifying when super-engager T cells are most prevalent.
- **Cluster 9 Focus**: By concentrating on cluster 9, it provides an in-depth look at the behavior of T cells demonstrating high engagement levels.
- **CD4 vs. CD8 Comparison**: Offers valuable insights into distinct engagement patterns between CD4 and CD8 T cells, enhancing the understanding of immune mechanisms.

For additional analysis or adjustments to visualization parameters, consult the [script](/scripts/TimepointGraph/timepoint_graph.R). This might include altering smoothing functions or the analysis time window, depending on the specifics of your experimental setup.

----
### (2) scRNAseq data preprocessing
This module contains scripts to process and analyze single cell sequencing data of T cells.
To test this pipeline you can either use the provided demos datasets or download an example sequencing data deposited in the GEO depository with the accession number [GSE172325](https://www-ncbi-nlm-nih-gov.ezproxy.u-pec.fr/geo/query/acc.cgi?acc=GSE172325). 
This module is structured in the following way:
1. [Pre-processing of raw SORTseq data (QC, filtering)](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part1_SORTseq_QC.ipynb)
2. [Subset Identification](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part2_Subset_Identification.ipynb)
3. [Pseudotime Analysis](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part3_Pseudotime.ipynb)
4. [Dynamic Gene Clustering](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part4_Dynamic_Genes.ipynb)
5. [Comparison to In vivo Data](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Part5_In_vivo_Comparison.ipynb)

----

### (3)  Population separation in silico simulation  

[This module](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Mutiple_timepoint_population_seperation_simulation.R) simulates T-cell population separation into engagers and non-engagers over multiple time points, providing insights into T-cell dynamics in a co-culture environment. It utilizes user-defined parameters from the configuration file to analyze and visualize the engagement behavior of CD4 and CD8 T-cells.

#### Configuration
Set up the `config_template.yml` with the necessary parameters for this module, including:
- `n_timepoints`: The number of separation time points to include in the simulation.
- `first_timepoint_interval_starts` and `first_timepoint_interval_ends`: The starting and ending minutes of the first time interval for separation.
- `reoccuring_time`: Time difference (in minutes) for subsequent time intervals.

***To run from the command line:***
```bash
Rscript /path/to/PopulationSeparationSimulation/population_separation_simulation.R -c </Path/to/BGT/config_template.yml> -f
```
Use the `-f` flag to force re-import and processing of data if output files already exist. The `-t` option allows specifying an RDS file with processed T-cell track data, if not included in the BGT config.

***To run from RStudio:***

**Step 2:** For a demo, navigate to the [population_separation_simulation script](/scripts/PopulationSeparationSimulation/population_separation_simulation.R). Adjust the **BGT config** file path on [line 18](/scripts/PopulationSeparationSimulation/population_separation_simulation.R#L18) if working with new data in a different folder.

***Output Files***

Generated outputs in the `Results/Population_seperation_simulation/` directory include:
- Plots for CD8 and CD4 T-cell engagement at specified intervals, indicating the percentage of cells in different clusters.
  - `plot_cd8_timepoint_X.png` and `plot_cd4_timepoint_X.png` where X denotes the specific time interval.
- CSV files detailing the engagement frequency analysis for CD8 and CD4 T-cells:
  - `CD8_engagement_behavior_freq.csv` and `CD4_engagement_behavior_freq.csv`, providing a quantitative measure of T-cell engagement over time.

#### Key Features
- **Dynamic Analysis**: Enables a comprehensive view of T-cell behavior across several time points, facilitating the understanding of engagement dynamics.
- **Flexible Time Points**: Leverages user-defined intervals to tailor the analysis to specific experimental setups, enhancing relevance to your research.
- **Detailed Visualization**: Offers clear, visual representations of T-cell engagement, aiding in the identification of patterns and trends in cellular behavior.

Adjustments to the analysis parameters or visualization aspects can be made by modifying the [script](/scripts/PopulationSeparationSimulation/population_separation_simulation.R) or the config file, depending on the specifics of the experimental design or analysis goals.

----

### (4) Behavioral probability mapping module

This module uses the principle of probability transitivity to infer for each sequenced T cell a probability to belonging to a particular T cell beahvioral class.  
[Behavioral probability mapping module](https://github.com/AlievaRios/BGT/blob/dev_avi/scripts/Behavioral-guided_transcriptomics.R) uses as input processes Seurat object from scRNA analysis module along with Tcell_engagement_behavior_freq.csv from the  Population separation in silico simulation module.


#### Configuration
Ensure the `config_template.yml` is correctly configured for this module:
- `scRNA_seq_dataset`: Path to the single-cell RNA sequencing dataset file.
- Other parameters relevant to this module are set directly within the script or derived from the behavioral analysis outputs.

***To run from the command line:***
```bash
Rscript /path/to/BehavioralGuidedTranscriptomics/behavioral_guided_transcriptomics.R -c </Path/to/BGT/config_template.yml> -f
```
Use `-f` to force re-import and processing of data if output files already exist. Optionally, `-t` can specify an RDS file with processed T-cell track data not included in the BGT config.

***To run from RStudio:***

**Step 3:** For demonstration purposes, utilize the [Behavioral probability mapping script](/scripts/BehavioralGuidedTranscriptomics/behavioral_guided_transcriptomics.R). If using new data or a different configuration, update the **BGT config** file path on [line 18](/scripts/BehavioralGuidedTranscriptomics/behavioral_guided_transcriptomics.R#L18).

***Output Files***

The module outputs include:
- UMAP plots illustrating the distribution of behavioral signatures across the scRNA-seq dataset.
- Probability maps and CSV files detailing the behavioral signature composition for each analyzed cluster.
- Heatmaps and PDFs visualizing the transition of T-cell behaviors along the pseudotime trajectory.

Outputs are stored in `Results/Behavioral_guided_transcriptomics/` and its subdirectories.

#### Key Features
- **Behavioral Signature Mapping**: Elucidates the transcriptomic underpinnings of T-cell behaviors, offering a nuanced view of immune cell dynamics.
- **Probability Maps**: Constructs detailed maps showing the likelihood of each cell belonging to different behavioral categories, enabling a deep dive into cellular function and heterogeneity.
- **Pseudotime Analysis**: Integrates pseudotime analysis to explore the evolution of T-cell states and behaviors over time, revealing insights into cellular progression and differentiation in response to tumor exposure.

For modifications or further analysis, refer to the [script](/scripts/BehavioralGuidedTranscriptomics/behavioral_guided_transcriptomics.R) and adjust parameters or visualization settings as needed based on your specific research questions or experimental design.


## Example output dataset
To enable users to delve into an illustrative dataset generated with BGT, we offer an example showcasing CD8+ TEGs co-cultured with Breast Cancer Organoids, accessible via an [online UCSC browser](https://cells.ucsc.edu/?ds=behavior-guided-tx). For information on how to use it go to [this page](https://cellbrowser.readthedocs.io/en/master/interface.html) 
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


