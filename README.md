# Notes on ChIP-seq and other-seq-related tools

These notes are not intended to be comprehensive. They include notes about methods, packages and tools I would like to explore. For a comprehensive overview of the subject, consider [other bioinformatics resources](https://github.com/mdozmorov/Bioinformatics_notes) and [collections of links to various resources](https://github.com/mdozmorov/MDmisc_notes). Issues with suggestions and pull requests are welcome!

# Table of content

* [ChIP-seq](#chip-seq)
  * [ChIP-seq pipelines](#chip-seq-pipelines)
  * [Peaks](#peaks)
  * [Visualization](#visualization)
  * [Motifs](#motifs)
  * [Differential peak detection](#differential-peak-detection)
  * [Interpretation](#interpretation)
  * [Misc](#misc)
* [DNAse-seq](#dnase-seq)
* [ATAC-seq](#atac-seq)
  * [ATAC-seq pipelines](#atac-seq-pipelines)
* [Histone-seq](#histone-seq)
* [Data](#data)
  * [Motif DBs](#motif-dbs)


## ChIP-seq

### ChIP-seq pipelines

- Regulatory Genomics Toolbox: Python library and set of tools for the integrative analysis of high throughput regulatory genomics data. http://www.regulatory-genomics.org, https://github.com/CostaLab/reg-gen
    - HINT (Hmm-based IdeNtification of Transcription factor footprints) is a framework for detection of DNA footprints from DNase-Seq and histone modification ChIP-Seq data.
    - Motif Analysis tools allows the search of motifs with binding sites enriched in particular genomic regions
    - ODIN and THOR are HMM-based approaches to detect and analyse differential peaks in pairs of ChIP-seq data. 
    - RGT-Viz is a collection of tests for association analysis and tools for visualizaiton of genomic data such as files in BED and BAM format
    - Triplex Domain Finder (TDF) statistically characterizes the triple helix potential of RNA and DNA regions.

- ENCODE3 pipeline v1 specifications, https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq

ChIP-seq peak calling using MACS2: `macs2 callpeak -f BAMPE -g hs -B --SPMR --verbose 3 --cutoff- analysis --call-summits -q 0.01`

- `AQUAS` ChIP-seq processing pipeline - The AQUAS pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications, https://github.com/kundajelab/chipseq_pipeline

- `Crunch` - Completely Automated Analysis of ChIP-seq Data, http://crunch.unibas.ch/, https://www.biorxiv.org/content/early/2016/03/09/042903

- `ChiLin` - QC, peak calling, motif analysis for ChIP-seq and DNAse-seq data used by CistromeDb. References to other tools. https://github.com/cfce/chilin
    - Qin, Qian, Shenglin Mei, Qiu Wu, Hanfei Sun, Lewyn Li, Len Taing, Sujun Chen, et al. “ChiLin: A Comprehensive ChIP-Seq and DNase-Seq Quality Control and Analysis Pipeline.” BMC Bioinformatics 17, no. 1 (October 3, 2016): 404. https://doi.org/10.1186/s12859-016-1274-4.

- `ChIP-eat` - a pipeline for aligning reads, calling peaks, predicting TFBSs. https://bitbucket.org/CBGR/chip-eat/src/master/
    - Gheorghe, Marius, Geir Kjetil Sandve, Aziz Khan, Jeanne Chèneby, Benoit Ballester, and Anthony Mathelier. “A Map of Direct TF–DNA Interactions in the Human Genome.” Nucleic Acids Research 47, no. 4 (February 28, 2019): e21–e21. https://doi.org/10.1093/nar/gky1210.

- `ChIPLine` - a pipeline for ChIP-seq analysis, https://github.com/ay-lab/ChIPLine

- `S3norm` - Chip-seq normalization to sequencing depth AND signal-to-noise ratio to the common reference. Negative Binomial for modeling background, convert counts to -log10(p-values), use monotonic nonlinear model to match the means of the common peaks and backgrounds in two datasets. https://github.com/guanjue/S3norm
    - Xiang, Guanjue, Cheryl Keller, Belinda Giardine, Lin An, Ross Hardison, and Yu Zhang. “S3norm: Simultaneous Normalization of Sequencing Depth and Signal-to-Noise Ratio in Epigenomic Data.” BioRxiv, January 1, 2018, 506634. https://doi.org/10.1101/506634.


### Peaks

- `Genrich` - Detecting sites of genomic enrichment in ChIP-seq and ATAC-seq. https://github.com/jsh58/Genrich, unpublished but highly tested and recommented, https://informatics.fas.harvard.edu/atac-seq-guidelines.html

- `mosaics` - This package provides functions for fitting MOSAiCS and MOSAiCS-HMM, a statistical framework to analyze one-sample or two-sample ChIP-seq data of transcription factor binding and histone modification. https://bioconductor.org/packages/release/bioc/html/mosaics.html

- `ROSE` - rank-ordering of super-enhancers. https://bitbucket.org/young_computation/rose

- `RSEG` - ChIP-seq broad domain analysis. http://smithlabresearch.org/software/rseg/
    - Song, Qiang, and Andrew D Smith. “Identifying Dispersed Epigenomic Domains from ChIP-Seq Data.” Bioinformatics 27, no. 6 (2011): 870–871.

- `triform` - finds enriched regions (peaks) in transcription factor ChIP-sequencing data. https://bioconductor.org/packages/release/bioc/html/triform.html


### Visualization

- Visualizations of ChIP-Seq data using Heatmaps, https://www.biostars.org/p/180314/

- `ChAsE` - Chromatin Analysis & Exploration Tool. http://chase.cs.univie.ac.at/overview

- `ChromHeatMap` - Heat map plotting by genome coordinate. https://bioconductor.org/packages/release/bioc/html/ChromHeatMap.html

- `BAM2WIG` - a flexible tool to generate read coverage profile (WIG file) from a BAM file. http://www.epigenomes.ca/tools-and-software

- `EaSeq` - peak calling (MACS), visualization, and analysis of ChIP-seq experiments. GUI, Windows-based, stand-alone. Figure 1, 3 - range of functionality, compared with other tools. https://easeq.net/downloadeaseq/. Description of tools: http://easeq.net/tools.pdf, Visualization examples: http://easeq.net/plots.pdf, Workflow examples: http://easeq.net/examples.pdf
    - Lerdrup, Mads, Jens Vilstrup Johansen, Shuchi Agrawal-Singh, and Klaus Hansen. “An Interactive Environment for Agile Analysis and Visualization of ChIP-Sequencing Data.” Nature Structural & Molecular Biology 23, no. 4 (April 2016): 349–57. https://doi.org/10.1038/nsmb.3180.

- `UROPA` - Universal RObustPeak Annotator. http://loosolab.mpi-bn.mpg.de/

- `Zerone` - combine multiple ChIP-seq profiles into one discretized profile. HMM with zero-inflated negative multinomial emissions across windowed genome. QC using SVM trained on ENCODE data to distinguish good from bad samples. Requires two negative controls. Compared against peaks called by MACS, BayesPeak, JAMM. https://github.com/nanakiksc/zerone
    - Cuscó, Pol, and Guillaume J. Filion. “Zerone: A ChIP-Seq Discretizer for Multiple Replicates with Built-in Quality Control.” Bioinformatics 32, no. 19 (October 1, 2016): 2896–2902. https://doi.org/10.1093/bioinformatics/btw336.


### Motifs

- `rGADEM` - de novo motif discovery, https://bioconductor.org/packages/release/bioc/html/rGADEM.html

- `gimmemotifs` - framework for TF motif analysis using an ensemble of motif predictors. `maelstrom` tool to detect differential motif activity between multiple different conditions. Includes manually curated database of motifs. Benchmark of 14 motif detection tools - Homer, MEME, BioProspector are among the top performing. Extensive analysis results. https://github.com/vanheeringen-lab/gimmemotifs, and documentation https://gimmemotifs.readthedocs.io
    - Bruse, Niklas, and Simon J. van Heeringen. “GimmeMotifs: An Analysis Framework for Transcription Factor Motif Analysis,” November 20, 2018. https://doi.org/10.1101/474403.

- `motifStack` - Plot stacked logos for single or multiple DNA, RNA and amino acid sequence, https://bioconductor.org/packages/release/bioc/html/motifStack.html

- `Logolas` - R package for Enrichment Depletion Logos (EDLogos) and String Logos, https://github.com/kkdey/Logolas

- `homerkit` - Read HOMER motif analysis output in R. https://github.com/slowkow/homerkit


### Differential peak detection

- `csaw` - Detection of differentially bound regions in ChIP-seq data with sliding windows, with methods for normalization and proper FDR control. https://bioconductor.org/packages/release/bioc/html/csaw.html

- `DiffBind` - Differential Binding Analysis of ChIP-Seq Peak Data. https://bioconductor.org/packages/release/bioc/html/DiffBind.html

### Interpretation

- `Toolkit for Cistrome Data Browser` - online tool to answer questions like:
    - What factors regulate your gene of interest?
    - What factors bind in your interval?
    - What factors have a significant binding overlap with your peak set?

- Chongzhi Zhang software page, http://faculty.virginia.edu/zanglab/software.htm
    - `BART` (Binding Analysis for Regulation of Transcription), a bioinformatics tool for predicting functional transcription factors (TFs) that bind at genomic cis-regulatory regions to regulate gene expression in the human or mouse genomes, given a query gene set or a ChIP-seq dataset as input. http://bartweb.uvasomrc.io/
    - `MARGE` (Model-based Analysis of Regulation of Gene Expression), a comprehensive computational method for inference of cis-regulation of gene expression leveraging public H3K27ac genomic profiles in human or mouse. http://cistrome.org/MARGE/
    - `MANCIE` (Matrix Analysis and Normalization by Concordant Information Enhancement), a computational method for high-dimensional genomic data integration. https://cran.r-project.org/web/packages/MANCIE/index.html
    - `SICER` (Spatial-clustering Identification of ChIP-Enriched Regions), a ChIP-Seq data analysis method. https://home.gwu.edu/~wpeng/Software.htm


### Misc

- mm10 blacklist - http://usevision.org/data/blacklist/blacklist.full.bed


## DNAse-seq

- DNAse-seq analysis guide. Tools for QC, peak calling, analysis, footprint detection, motif analysis, visualization, all-in-one tools ([Table 2](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby057/5053117#118754375))
    - Liu, Yongjing, Liangyu Fu, Kerstin Kaufmann, Dijun Chen, and Ming Chen. “A Practical Guide for DNase-Seq Data Analysis: From Data Management to Common Applications.” Briefings in Bioinformatics, July 12, 2018. https://doi.org/10.1093/bib/bby057.


## ATAC-seq

- ATAC-seq peak-calling and differential analysis pipeline. https://github.com/nf-core/atacseq, http://nf-co.re

- Dimensionality Reduction for scATAC-seq Data, http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html

### ATAC-seq pipelines

ATAC-seq peak calling using MACS2: `macs2 callpeak --nomodel --nolambda -- keep-dup all --call-summits -f BAMPE -g hs`

- `ATACProc` - ATAC-seq processing pipeline, https://github.com/ay-lab/ATACProc

- `HMMRATAC` - hidden Markov model for ATAC-seq to identify open chromatin regions. Parametric modeling of nucleosome-free regions and three nucleosomal reatures (mono-, di-, and tri-nucleosomes). First, train on 1000 auto-selected regions, then predict. Tested on "active promoters" and "strong enhancers" chromatin states (positive examples), and "heterochromatin" (negative examples). Compared with MACS2, F-seq. https://github.com/LiuLabUB/HMMRATAC
    - Tarbell, Evan D, and Tao Liu. “HMMRATAC: A Hidden Markov ModeleR for ATAC-Seq.” Nucleic Acids Research, June 14, 2019, gkz533. https://doi.org/10.1093/nar/gkz533.

- `pepatac` - A modular, containerized pipeline for ATAC-seq data processing. https://github.com/databio/pepatac, examples and documentation at http://code.databio.org/PEPATAC/


## Histone-seq

Homer program ‘findPeaks’ with the style ‘histone’. Peaks within 1 kb were merged into a single peak. Broad peaks in H3K36me3, H3K27me3 and H3K9me3 were called using the Homer program ‘findPeaks’ with the options ‘-region –size 1000 –minDist 2500’. When Homer runs with these options, the initial sets of peaks were 1 kb wide and peaks within 2.5 kb were merged. 

- `DEScan2` - broad peak (histone, ATAC, DNAse) analysis (peak caller, peak filtering and alignment across replicates, creation of a count matrix). Peak caller uses a moving window and calculated a Poisson likelihood of a peak as compared to a region outside the window. https://bioconductor.org/packages/release/bioc/html/DEScan2.html
    - Righelli, Dario, John Koberstein, Nancy Zhang, Claudia Angelini, Lucia Peixoto, and Davide Risso. “Differential Enriched Scan 2 (DEScan2): A Fast Pipeline for Broad Peak Analysis.” PeerJ Preprints, 2018.

- `HMCan` and `HMCan-diff` - histone ChIP-seq peak caller (and differential) that accounts for CNV, also for CG bias. Hidden Markov Model to detect peak signal. Control-FREEC to detect CNV in ChIP-seq data. Outperforms others, CCAT second best. https://www.cbrc.kaust.edu.sa/hmcan/
    - Ashoor, Haitham, Aurélie Hérault, Aurélie Kamoun, François Radvanyi, Vladimir B. Bajic, Emmanuel Barillot, and Valentina Boeva. “HMCan: A Method for Detecting Chromatin Modifications in Cancer Samples Using ChIP-Seq Data.” Bioinformatics (Oxford, England) 29, no. 23 (December 1, 2013): 2979–86. https://doi.org/10.1093/bioinformatics/btt524.
    - Ashoor, Haitham, Caroline Louis-Brennetot, Isabelle Janoueix-Lerosey, Vladimir B. Bajic, and Valentina Boeva. “HMCan-Diff: A Method to Detect Changes in Histone Modifications in Cells with Different Genetic Characteristics.” Nucleic Acids Research 45, no. 8 (05 2017): e58. https://doi.org/10.1093/nar/gkw1319.

- `RSEG` - ChIP-seq analysis for identifying genomic regions and their boundaries marked by diffusive histone modification markers, such as H3K36me3 and H3K27me3, http://smithlabresearch.org/software/rseg/


## Data

- `Cistrome DB` - ChIP-seq peaks for TFs, histone modifications, DNAse/ATAC. Downloadable cell type-specific, hg38 BED files. http://cistrome.org/db/#/
    - [Toolkit](http://dbtoolkit.cistrome.org/) to answer questions like "What factors regulate your gene of interest?", "What factors bind in your interval?", "What factors have a significant binding overlap with your peak set?"

- `ChIP-Atlas` is an integrative and comprehensive database for visualizing and making use of public ChIP-seq data. ChIP-Atlas covers almost all public ChIP-seq data submitted to the SRA (Sequence Read Archives) in NCBI, DDBJ, or ENA, and is based on over 78,000 experiments. Besides "Peak browser" tool, includes "Target Genes", "Colocalization", "Enrichment Analysis" tools. http://chip-atlas.org/, code on GetHub, https://github.com/inutano/chip-atlas

- `CODEX ChIP-seq` - CODEX provides access to processed and curated NGS experiments, including ChIP-Seq (transcription factors and histones), RNA-Seq and DNase-Seq. Human, mouse. Download tracks, analyze correlations, motifs, compare between organisms, more. http://codex.stemcells.cam.ac.uk/

- `RAEdb` - enhancer database. Enhancers identified from STARR-seq and MPRA studies. Epromoters - promoters containing enhancers. Human (hg38)/mouse (mm10) data, select cell lines. BED/FASTQ download. Links to EnhancerAtlas, VISTA, SuperEnhancer databases. http://www.computationalbiology.cn/RAEdb/index.php
    - Cai, Zena, Ya Cui, Zhiying Tan, Gaihua Zhang, Zhongyang Tan, Xinlei Zhang, and Yousong Peng. “RAEdb: A Database of Enhancers Identified by High-Throughput Reporter Assays.” Database: The Journal of Biological Databases and Curation 2019 (January 1, 2019). https://doi.org/10.1093/database/bay140.

- `ReMap 2018` - database of public and ENCODE human ChIP-seq datasets in BED format, in different coordinate systems, downloadable. http://pedagogix-tagc.univ-mrs.fr/remap/index.php

- `UniBind` database of robustly predicted TF binding sites. Using ChIP-seq data from 1983 studies, PWMs, binding energy, and many other parameters. HOT/XOT regions are likely artifacts, not TFBSs. Downloadable (hg38) database.  https://unibind.uio.no/
    - Gheorghe, Marius, Geir Kjetil Sandve, Aziz Khan, Jeanne Chèneby, Benoit Ballester, and Anthony Mathelier. “A Map of Direct TF–DNA Interactions in the Human Genome.” Nucleic Acids Research 47, no. 4 (February 28, 2019): e21–e21. https://doi.org/10.1093/nar/gky1210.

### Motif DBs

- `HOCOMOCO` - TF binding models from ChIP-seq data. http://hocomoco11.autosome.ru/

