# Notes on ChIP-seq and other-seq-related tools

Issues with suggestions and pull requests are welcome!

# Table of content

* [ChIP-seq](#chip-seq)
  * [ChIP-seq pipelines](#chip-seq-pipelines)
  * [Peaks](#peaks)
  * [Visualization](#visualization)
  * [Motifs](#motifs)
  * [Differential peak detection](#differential-peak-detection)
  * [Misc](#misc)
* [ATAC-seq](#atac-seq)
  * [ATAC-seq pipelines](#atac-seq-pipelines)
* [Histone-seq](#histone-seq)
* [Data](#data)
  * [Motif DBs](#motif-dbs)


## ChIP-seq

### ChIP-seq pipelines

ChIP-seq peak calling using MACS2: `macs2 callpeak -f BAMPE -g hs -B --SPMR --verbose 3 --cutoff- analysis --call-summits -q 0.01`

- `AQUAS` ChIP-seq processing pipeline - The AQUAS pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications, https://github.com/kundajelab/chipseq_pipeline

- `Crunch` - Completely Automated Analysis of ChIP-seq Data, http://crunch.unibas.ch/, https://www.biorxiv.org/content/early/2016/03/09/042903

- `ChiLin` - QC, peak calling, motif analysis for ChIP-seq and DNAse-seq data used by CistromeDb. References to other tools. https://github.com/cfce/chilin
    - Qin, Qian, Shenglin Mei, Qiu Wu, Hanfei Sun, Lewyn Li, Len Taing, Sujun Chen, et al. “ChiLin: A Comprehensive ChIP-Seq and DNase-Seq Quality Control and Analysis Pipeline.” BMC Bioinformatics 17, no. 1 (October 3, 2016): 404. https://doi.org/10.1186/s12859-016-1274-4.

- `ChIPLine` - a pipeline for ChIP-seq analysis, https://github.com/ay-lab/ChIPLine

### Peaks

- `mosaics` - This package provides functions for fitting MOSAiCS and MOSAiCS-HMM, a statistical framework to analyze one-sample or two-sample ChIP-seq data of transcription factor binding and histone modification. https://bioconductor.org/packages/release/bioc/html/mosaics.html

- `triform` - finds enriched regions (peaks) in transcription factor ChIP-sequencing data. https://bioconductor.org/packages/release/bioc/html/triform.html


### Visualization

- `ChAsE` - Chromatin Analysis & Exploration Tool. http://chase.cs.univie.ac.at/overview

- `ChromHeatMap` - Heat map plotting by genome coordinate. https://bioconductor.org/packages/release/bioc/html/ChromHeatMap.html

- `BAM2WIG` - a flexible tool to generate read coverage profile (WIG file) from a BAM file. http://www.epigenomes.ca/tools-and-software

- `UROPA` - Universal RObustPeak Annotator. http://loosolab.mpi-bn.mpg.de/

- `Zerone` - combine multiple ChIP-seq profiles into one discretized profile. HMM with zero-inflated negative multinomial emissions across windowed genome. QC using SVM trained on ENCODE data to distinguish good from bad samples. Requires two negative controls. Compared against peaks called by MACS, BayesPeak, JAMM. https://github.com/nanakiksc/zerone
    - Cuscó, Pol, and Guillaume J. Filion. “Zerone: A ChIP-Seq Discretizer for Multiple Replicates with Built-in Quality Control.” Bioinformatics 32, no. 19 (October 1, 2016): 2896–2902. https://doi.org/10.1093/bioinformatics/btw336.


### Motifs

- `rGADEM` - de novo motif discovery, https://bioconductor.org/packages/release/bioc/html/rGADEM.html

- `gimmemotifs` - framework for TF motif analysis using an ensemble of motif predictors. `maelstrom` tool to detect differential motif activity between multiple different conditions. Includes manually curated database of motifs. Benchmark of 14 motif detection tools - Homer, MEME, BioProspector are among the top performing. Extensive analysis results. https://github.com/vanheeringen-lab/gimmemotifs, and documentation https://gimmemotifs.readthedocs.io
    - Bruse, Niklas, and Simon J. van Heeringen. “GimmeMotifs: An Analysis Framework for Transcription Factor Motif Analysis,” November 20, 2018. https://doi.org/10.1101/474403.

- `motifStack` - Plot stacked logos for single or multiple DNA, RNA and amino acid sequence, https://bioconductor.org/packages/release/bioc/html/motifStack.html

- `Logolas` - R package for Enrichment Depletion Logos (EDLogos) and String Logos, https://github.com/kkdey/Logolas


### Differential peak detection

- `csaw` - Detection of differentially bound regions in ChIP-seq data with sliding windows, with methods for normalization and proper FDR control. https://bioconductor.org/packages/release/bioc/html/csaw.html



### Misc





## ATAC-seq

### ATAC-seq pipelines

ATAC-seq peak calling using MACS2: `macs2 callpeak --nomodel --nolambda -- keep-dup all --call-summits -f BAMPE -g hs`

- `ATACProc` - ATAC-seq processing pipeline, https://github.com/ay-lab/ATACProc

- `pepatac` - A modular, containerized pipeline for ATAC-seq data processing. https://github.com/databio/pepatac, examples and documentation at http://code.databio.org/PEPATAC/


## Histone-seq

Homer program ‘findPeaks’ with the style ‘histone’. Peaks within 1 kb were merged into a single peak. Broad peaks in H3K36me3, H3K27me3 and H3K9me3 were called using the Homer program ‘findPeaks’ with the options ‘-region –size 1000 –minDist 2500’. When Homer runs with these options, the initial sets of peaks were 1 kb wide and peaks within 2.5 kb were merged. 

- `DEScan2` - Differential Enrichment Scan 2 R package, Integrated peak and differential caller, specifically designed for broad epigenomic signals. http://bioconductor.org/packages/release/bioc/html/DEScan2.html

- `RSEG` - ChIP-seq analysis for identifying genomic regions and their boundaries marked by diffusive histone modification markers, such as H3K36me3 and H3K27me3, http://smithlabresearch.org/software/rseg/


## Data

- `Cistrome DB` - ChIP-seq peaks for TFs, histone modifications, DNAse/ATAC. Downloadable cell type-specific, hg38 BED files. http://cistrome.org/db/#/
    - [Toolkit](http://dbtoolkit.cistrome.org/) to answer questions like "What factors regulate your gene of interest?", "What factors bind in your interval?", "What factors have a significant binding overlap with your peak set?"

- `ChIP-Atlas` is an integrative and comprehensive database for visualizing and making use of public ChIP-seq data. ChIP-Atlas covers almost all public ChIP-seq data submitted to the SRA (Sequence Read Archives) in NCBI, DDBJ, or ENA, and is based on over 78,000 experiments. Besides "Peak browser" tool, includes "Target Genes", "Colocalization", "Enrichment Analysis" tools. http://chip-atlas.org/, code on GetHub, https://github.com/inutano/chip-atlas

- `ReMap 2018` - database of public and ENCODE human ChIP-seq datasets in BED format, in different coordinate systems, downloadable. http://pedagogix-tagc.univ-mrs.fr/remap/index.php

### Motif DBs

- `HOCOMOCO` - TF binding models from ChIP-seq data. http://hocomoco11.autosome.ru/

