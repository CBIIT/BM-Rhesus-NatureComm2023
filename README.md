# BM-Rhesus-NatureComm2023


BubblePlot foler: generate_bubble_plot.R and required example inputs

## Description

Example code to generate bubble plot. 
This script generates bubble plot based on TOA (Time Of Acquisition) data, ATAC-Seq accessibility data and annotation data of CREB1.

## Getting Started

The version indicates the tested working version of R libraries.

### Prerequisite when R v4.1.2 is used

* rlang_1.0.2
* rtracklayer_1.54.0
* GenomicRanges_1.46.1
* scales_1.2.1
* wesanderson_0.3.6
* viridis_0.6.2
* ggplot2_3.3.5

### Prerequisite when R v4.2.2 is used

* rlang_1.0.6
* rtracklayer_1.56.0 or tracklayer_1.58.0
* GenomicRanges_1.50.1 or GenomicRanges_1.50.2
* scales_1.2.1
* wesanderson_0.3.6
* viridis_0.6.2
* ggplot2_3.4.0


### Required input data (provided for example gene, CREB1)

#### Processed example ATAC-Seq data (full data set available at GEO; GSE188879, GSE189032)

* atac_110098_peak_list.csv
* atac.norm.counts.basal.chr12.csv
* atac.norm.counts.post.chr12.csv
* example_subset_CREB1_only_Macaca_mulatta.Mmul_10.103.chr.gtf

### Expected output

* bubble-test.pdf will be generated under the current working directory.

## Workflow
* source('generate_bubble_plot.R')
