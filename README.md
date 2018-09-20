# Notes on ChIP-seq and other-seq-related tools

Issues with suggestions and pull requests are welcome!

# Table of content

* [ChIP-seq](#chip-seq)
  * [ChIP-seq pipelines](#chip-seq-pipelines)
  * [Visualization](#visualization)
  * [Motif plots](#motif-plots)
* [ATAC-seq](#atac-seq)
  * [ATAC-seq pipelines](#atac-seq-pipelines)
* [Histone-seq](#histone-seq)

## ChIP-seq

### ChIP-seq pipelines

ChIP-seq peak calling using MAXS2: `macs2 callpeak -f BAMPE -g hs -B --SPMR --verbose 3 --cutoff- analysis --call-summits -q 0.01`

- `AQUAS` ChIP-seq processing pipeline - The AQUAS pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications, https://github.com/kundajelab/chipseq_pipeline

- `Crunch` - Completely Automated Analysis of ChIP-seq Data, http://crunch.unibas.ch/, https://www.biorxiv.org/content/early/2016/03/09/042903

- `ChIPLine` - a pipeline for ChIP-seq analysis, https://github.com/ay-lab/ChIPLine

### Visualization

- `ChAsE` - Chromatin Analysis & Exploration Tool. http://chase.cs.univie.ac.at/overview

- `ChromHeatMap` - Heat map plotting by genome coordinate. https://bioconductor.org/packages/release/bioc/html/ChromHeatMap.html

- `BAM2WIG` - a flexible tool to generate read coverage profile (WIG file) from a BAM file. http://www.epigenomes.ca/tools-and-software

### Motif plots

- `motifStack` - Plot stacked logos for single or multiple DNA, RNA and amino acid sequence, https://bioconductor.org/packages/release/bioc/html/motifStack.html

- `Logolas` - R package for Enrichment Depletion Logos (EDLogos) and String Logos, https://github.com/kkdey/Logolas

### Misc

- `csaw` - Detection of differentially bound regions in ChIP-seq data with sliding windows, with methods for normalization and proper FDR control. https://bioconductor.org/packages/release/bioc/html/csaw.html

- `mosaics` - This package provides functions for fitting MOSAiCS and MOSAiCS-HMM, a statistical framework to analyze one-sample or two-sample ChIP-seq data of transcription factor binding and histone modification. https://bioconductor.org/packages/release/bioc/html/mosaics.html

- `triform` - finds enriched regions (peaks) in transcription factor ChIP-sequencing data. https://bioconductor.org/packages/release/bioc/html/triform.html

- `rGADEM` - de novo motif discovery, https://bioconductor.org/packages/release/bioc/html/rGADEM.html



## ATAC-seq

### ATAC-seq pipelines

ATAC-seq peak calling using MACS2: `macs2 callpeak --nomodel --nolambda -- keep-dup all --call-summits -f BAMPE -g hs`

- `ATACProc` - ATAC-seq processing pipeline, https://github.com/ay-lab/ATACProc


## Histone-seq

Homer program ‘findPeaks’ with the style ‘histone’. Peaks within 1 kb were merged into a single peak. Broad peaks in H3K36me3, H3K27me3 and H3K9me3 were called using the Homer program ‘findPeaks’ with the options ‘-region –size 1000 –minDist 2500’. When Homer runs with these options, the initial sets of peaks were 1 kb wide and peaks within 2.5 kb were merged. 

- `RSEG` - ChIP-seq analysis for identifying genomic regions and their boundaries marked by diffusive histone modification markers, such as H3K36me3 and H3K27me3, http://smithlabresearch.org/software/rseg/






