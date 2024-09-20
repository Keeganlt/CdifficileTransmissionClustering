# Source code for producing the results and figures

The code is divided into scripts for each figure or analysis. Most scripts only rely on [Load_Data](Load_Data.R) and [GlobalFunctions](GlobalFunctions.R) however, some functions rely on other scripts, those are noted below. 

* [Load_Data](Load_Data.R): loads in all of the data objects, colors used, and dependencies (see below). This code will only work for people who are on the IRB and have access to the data.
* [GlobalFunctions](GlobalFunctions.R): is where all of the functions that I use throughout multiple scripts are called. Descriptions of the functions and what they do are commented in this file.

The rest of the scripts contain individual analyses or figure generating code except for [PariwiseDistDat](PariwiseDistDat.R) and [ClusteringData](ClusteringData.R). These two generate data that is used in multiple other scripts. 
* [PariwiseDistDat](PariwiseDistDat.R): calculates the pairwise genomic distance beteween all isolates in the study
* [ClusteringData](ClusteringData.R): 

1. [TimeOnWard](TimeOnWard.R): calculates the time on ward for patients both with and without *C. difficile*
2. [Figure1](Figure1.R): calculates the prevalence and contains the code to make Figure 1. This script is dependent on [ClusteringData](ClusteringData.R)
3. [SamplingStats](SamplingStats.R): calculates the stats about the number of samples taken for patients
4. [ToxinStats](ToxinStats.R): calculates the stats about patients with and without toxigenic *C. difficile*
5. [Figure2](Figure2.R): contains the code to make Figure 2. This script is dependent on [PariwiseDistDat](PariwiseDistDat.R)
6. [Figure3](Figure3.R): contains the code to make Figure 3
7. [OddsRecovery](OddsRecovery.R): calculates the odds ratio of recovering a positive isolate given characteristics about the patient, their environment, and healthcare worker hands
8. [ClusteringStats](ClusteringStats.R): calculates stats around the probability of being in a cluster given characteristics about the patient, their environment, and healthcare worker hands
9. [Figure4](Figure4.R): contains the code to make Figure 4. This script is dependent on [ClusteringData](ClusteringData.R)
10. [DataDescriptive](DataDescriptive.R): contains code to make supplemental Figure 1
11. [Table1](Table1.R): contains code to make supplemental Table 1
12. [FigureSI2](FigureSI2.R): contains the code to make supplemental Figure 2
13. [PhyloTree](PhyloTree.R): code to make the phylogenetic tree with tips colored by patient ID, shaped by sampling location, indicating toxin status and clade
14. [FigureSI5](FigureSI5.R): contains the code to make supplemental Figure 5
15. [FigureSI6](FigureSI6.R): contains the code to make supplemental Figure 6. This script is dependent on [PariwiseDistDat](PariwiseDistDat.R)
16. [FigureSI7](FigureSI7.R): contains the code to make supplemental Figure 7
17. [FigureSI8](FigureSI8.R): contains the code to make supplemental Figure 8

## Dependencies

The following R packages are required for this project:

- [`ape`](https://cran.r-project.org/package=ape): for analysis of phylogenetics and evolution.
- [`cowplot`](https://cran.r-project.org/package=cowplot): for publication-quality plots with a custom theme.
- [`dplyr`](https://cran.r-project.org/package=dplyr): for data manipulation.
- [`GGally`](https://cran.r-project.org/package=GGally): for extended functionalities for `ggplot2`.
- [`gdata`](https://cran.r-project.org/package=gdata): for various data manipulation tasks.
- [`ggnewscale`](https://cran.r-project.org/package=ggnewscale): for adding multiple color scales to a `ggplot2` plot.
- [`ggplot2`](https://cran.r-project.org/package=ggplot2): for data visualization.
- [`ggraph`](https://cran.r-project.org/package=ggraph): for creating network visualizations.
- [`ggtext`](https://cran.r-project.org/package=ggtext): for improved text rendering in `ggplot2`.
- [`ggtree`](https://bioconductor.org/packages/release/bioc/html/ggtree.html): for visualization and annotation of phylogenetic trees.
- [`grid`](https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid-package.html): for low-level grid graphics functions.
- [`gridExtra`](https://cran.r-project.org/package=gridExtra): for arranging multiple grid-based plots.
- [`igraph`](https://cran.r-project.org/package=igraph): for network analysis and visualization.
- [`lubridate`](https://cran.r-project.org/package=lubridate): for working with dates and times.
- [`network`](https://cran.r-project.org/package=network): for creating and modifying network objects.
- [`png`](https://cran.r-project.org/package=png): for reading and writing PNG images.
- [`Polychrome`](https://cran.r-project.org/package=Polychrome): for generating qualitative color palettes.
- [`RColorBrewer`](https://cran.r-project.org/package=RColorBrewer): for color schemes for maps and other graphics.
- [`readxl`](https://cran.r-project.org/package=readxl): for reading Excel files.
- [`sna`](https://cran.r-project.org/package=sna): for social network analysis.
- [`stringr`](https://cran.r-project.org/package=stringr): for string manipulation.
- [`tidygraph`](https://cran.r-project.org/package=tidygraph): for tidy API for graph manipulation.
- [`tidyverse`](https://cran.r-project.org/package=tidyverse): a collection of R packages for data science.
- [`wesanderson`](https://cran.r-project.org/package=wesanderson): for color palettes inspired by Wes Anderson films.
