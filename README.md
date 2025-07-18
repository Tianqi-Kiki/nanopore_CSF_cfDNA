# Nanopore sequencing of cell free DNA identifies methylation and fragmentation profiles from cerebral spinal fluid from lung cancer brain metastases 

This project analyzed the cfDNA fragmentation, methylation (5mC) and hydroxymethylation (5hmC) profiles from cerebrospinal fluid (CSF) in non-small cell lung cancer (NSCLC) patients with brain metastases and non-cancer controls. We found distinct nucleosome ratio patterns between cancer and control, which is unique to CSF compared with plasma samples. Distinct 5mC and 5hmC patterns were also observed between cancer and control CSF. Pathway analysis revealed cancer related genes and protein markers in NSCLC. The scripts below were used to generate the main figures, including the sequencing coverage, cfDNA fragment size distributions (in both CSF and plasma), methylation analysis and pathway analysis. 

## Contents
### CSF_cfDNA_sequencing_coverage.R 
This is a R notebook-style script used in Rstudio to calculate the sequencing coverage of CSF-cfDNA using nanopore sequencing. This has been tested with the following:
* RStudio 1.4.1717
* R 4.1.0
* data.table 1.14
* dplyr 1.0.7
* ggplot2 3.3.5
* tidyverse 1.3.1
