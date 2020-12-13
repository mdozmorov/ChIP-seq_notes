# Notes on ChIP-seq and other-seq-related tools

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

ChIP-seq, ATAC-seq related tools and genomics data analysis resources. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [ChIP-seq](#chip-seq)
  - [ChIP-seq pipelines](#chip-seq-pipelines)
    - [CUT&RUN](#cutrun)
  - [Quality control](#quality-control)
  - [Peaks](#peaks)
  - [Visualization](#visualization)
  - [Motif search](#motif-search)
  - [Differential peak detection](#differential-peak-detection)
  - [Enrichment](#enrichment)
  - [Interpretation](#interpretation)
  - [Blacklisted](#blacklisted)
- [DNAse-seq](#dnase-seq)
- [ATAC-seq](#atac-seq)
  - [ATAC-seq pipelines](#atac-seq-pipelines)
- [Histone-seq](#histone-seq)
- [Misc](#misc)
- [Data](#data)
  - [Motif DBs](#motif-dbs)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## ChIP-seq

- [ChIP-seq-analysis](https://github.com/crazyhottommy/ChIP-seq-analysis) - ChIP-seq analysis notes from Ming Tang


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

#### CUT&RUN

- [CUT&Tag Data Processing and Analysis Tutorial](https://yezhengstat.github.io/CUTTag_tutorial/)
- [Methods with detailed commands of CUT&RUN data analysis](https://www.biorxiv.org/content/10.1101/2020.08.31.272856v1.full-text) - from Divya S. Vinjamur et al. "[ZNF410 represses fetal globin by devoted control of CHD4/NuRD](https://doi.org/10.1101/2020.08.31.272856)," bioRxiv, August 31, 2020.

- CUT&TAG technology Cleavage Under Target and Tagmentation. Compared with CUT&RUN that uses MNase, it uses Tn5 transposase, reactions performed within intact cells, performed on a solid support (tethered). Better suited for low cell numbers, low cost. Tested on H3K27me3 and RNAPII profiling in K562, compared with the same CUT&RUN data. Sharper peaks, nearly 20X more that ChIP-seq. Compared with ATAC-seq in K562, H3K4me2, better signal-to-noise ratio, even at low sequencing depth. Tested using NPAT and CTCF transcription factors. Methods - alignment (bowtie2) and peak calling (MACS2) settings.
    - Kaya-Okur, Hatice S., Steven J. Wu, Christine A. Codomo, Erica S. Pledger, Terri D. Bryson, Jorja G. Henikoff, Kami Ahmad, and Steven Henikoff. “[CUT&Tag for Efficient Epigenomic Profiling of Small Samples and Single Cells](https://doi.org/10.1038/s41467-019-09982-5).” Nature Communications 10, no. 1 (December 2019)

- CUT&RUN technology, chromatin profiling strategy, antibody-targeted controlled cleavage by micrococcal nuclease. Cost-efficient, low input requirements, easier.
    - Skene, Peter J, and Steven Henikoff. “[An Efficient Targeted Nuclease Strategy for High-Resolution Mapping of DNA Binding Sites](https://elifesciences.org/articles/21856).” Genes and Chromosomes

- [SEARC](https://seacr.fredhutch.org/) (Sparse Enrichment Analysis for CUT&RUN) peak caller for CUT&RUN data. Data-driven, peaks with respect to global background or IgG control. Compared to MACS2 and HOMER, more precise and maintains true positive rate at low read depth. Better call wide peaks. Input - bedGraph, output - BED. Command line and [web server](https://seacr.fredhutch.org/)
    - Meers, MP, Tenenbaum, D and Henikoff S (2019). "[Peak calling by sparse enrichment analysis for CUT&RUN chromatin profiling](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4)". Epigenetics & Chromatin 2019 12:42.

- [CUT&RUNTools](https://bitbucket.org/qzhudfci/cutruntools/src/master/) - a pipeline to fully process CUT&RUN data and identify protein binding and genomic footprinting from antibody-targeted primary cleavage data. Implemented in R, Python, Bach, runs under the SLURM job submission. At the core, creates a cut matrix of from enzyme cleavage data. Compared with Atactk and Centipede. (Tested, didn't work)
    - Zhu, Qian. “[CUT&RUNTools: A Flexible Pipeline for CUT&RUN Processing and Footprint Analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4),” 2019, 12. 

### Quality control

- `phantompeakqualtools` - This package computes informative enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data. It can also be used to obtain robust estimates of the predominant fragment length or characteristic tag shift values in these assays. https://github.com/kundajelab/phantompeakqualtools

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


### Motif search

- [Non-redundant TF motif matches genome-wide](https://www.vierstra.org/resources/motif_clustering#downloads) - Clustering of 2179 motif models. hg38/mm10 BED files download with coordinates

- `marge` - API for HOMER in R for Genomic Analysis using Tidy Conventions https://robertamezquita.github.io/marge/, https://github.com/robertamezquita/marge

- `LISA` - epigenetic Landscape In silico Subtraction analysis, enriched TFs and chromatin regulators in a list of genes. http://lisa.cistrome.org

- `AME-MEME` - enrichment of known user-provided motifs in a sequence, http://meme-suite.org/doc/ame.html, http://meme-suite.org/tools/ame

- `rGADEM` - de novo motif discovery, https://bioconductor.org/packages/release/bioc/html/rGADEM.html

- [gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs) - framework for TF motif analysis using an ensemble of motif predictors. `maelstrom` tool to detect differential motif activity between multiple different conditions. Includes manually curated database of motifs. Benchmark of 14 motif detection tools - Homer, MEME, BioProspector are among the top performing. Extensive analysis results. [Documentation](https://gimmemotifs.readthedocs.io). [Tweet with updates](https://twitter.com/svheeringen/status/1313778158999666688?s=20)
    - Bruse, Niklas, and Simon J. van Heeringen. “[GimmeMotifs: An Analysis Framework for Transcription Factor Motif Analysis](https://doi.org/10.1101/474403),” November 20, 2018

- `DECOD` - Differential motif finder. k-mer-based. http://gene.ml.cmu.edu/DECOD/
    - Huggins, Peter, Shan Zhong, Idit Shiff, Rachel Beckerman, Oleg Laptenko, Carol Prives, Marcel H. Schulz, Itamar Simon, and Ziv Bar-Joseph. “DECOD: Fast and Accurate Discriminative DNA Motif Finding.” Bioinformatics 27, no. 17 (September 1, 2011): 2361–67. https://doi.org/10.1093/bioinformatics/btr412.

- `motifStack` - Plot stacked logos for single or multiple DNA, RNA and amino acid sequence, https://bioconductor.org/packages/release/bioc/html/motifStack.html

- `Logolas` - R package for Enrichment Depletion Logos (EDLogos) and String Logos, https://github.com/kkdey/Logolas

- `homerkit` - Read HOMER motif analysis output in R. https://github.com/slowkow/homerkit


### Differential peak detection

- `csaw` - Detection of differentially bound regions in ChIP-seq data with sliding windows, with methods for normalization and proper FDR control. https://bioconductor.org/packages/release/bioc/html/csaw.html

- `DiffBind` - Differential Binding Analysis of ChIP-Seq Peak Data. https://bioconductor.org/packages/release/bioc/html/DiffBind.html

### Enrichment

- UniBind Enrichment Analysis, also differential enrichment. Input - BED file in hg38 version. [LOLA](https://bioconductor.org/packages/release/bioc/html/LOLA.html) as an enrichment and database engine.

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


### Blacklisted

- [Repetitive centromeric, telomeric and satellite regions known to have low sequencing confidence - blacklisted regions defined by the ENCODE project](http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz) - from Upton et al., “Epigenomic Profiling of Neuroblastoma Cell Lines.”

- [Manually annotated GRCh38 blacklisted regions](https://www.encodeproject.org/files/ENCFF356LFX/) on ENCODE data portal. [Tweet by Anshul Kundaje](https://twitter.com/anshulkundaje/status/1263546023151992832?s=20)

- Blacklisted regions, how they were created. https://github.com/Boyle-Lab/Blacklist/, https://www.encodeproject.org/annotations/ENCSR636HFF/
    - Amemiya, Haley M., Anshul Kundaje, and Alan P. Boyle. “The ENCODE Blacklist: Identification of Problematic Regions of the Genome.” Scientific Reports 9, no. 1 (December 2019): 9354. https://doi.org/10.1038/s41598-019-45839-z.

- hg19 blacklist - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz, http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz

- hg38 blacklist - http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz

- mm10 blacklist - http://usevision.org/data/blacklist/blacklist.full.bed

- [UCSC unusual regions on assembly structure](http://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=problematic&hgta_table=comments&hgta_doSchema=describe+table+schema), [Tweet](https://twitter.com/GenomeBrowser/status/1260693767125778434?s=20)


## DNAse-seq

- DNAse-seq analysis guide. Tools for QC, peak calling, analysis, footprint detection, motif analysis, visualization, all-in-one tools ([Table 2](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby057/5053117#118754375))
    - Liu, Yongjing, Liangyu Fu, Kerstin Kaufmann, Dijun Chen, and Ming Chen. “A Practical Guide for DNase-Seq Data Analysis: From Data Management to Common Applications.” Briefings in Bioinformatics, July 12, 2018. https://doi.org/10.1093/bib/bby057.


## ATAC-seq

- [Benchmarking ATAC-seq peak calling](https://bigmonty12.github.io/peak-calling-benchmark) by Austin Montgomery

- Awesome ATAC-seq analysis by Nathan Sheffield, https://github.com/databio/awesome-atac-analysis

- `TOBIAS` - Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal. https://github.com/loosolab/TOBIAS

- ATAC-seq peak-calling and differential analysis pipeline. https://github.com/nf-core/atacseq, http://nf-co.re

- Dimensionality Reduction for scATAC-seq Data, http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html

### ATAC-seq pipelines

- [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) - ATAC-seq and DNase-seq processing pipeline by Anshul Kundaje

- `scATAC-pro` - A comprehensive tool for processing, analyzing and visulizing single cell chromatin accessibility sequencing data. https://github.com/wbaopaul/scATAC-pro

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

## Misc

- `covtobed` - a tool to generate BED coverage tracks from BAM files. https://github.com/telatin/covtobed

## Data

- `Catchitt` - method for predicting TFBSs, meader of ENCODE-DREAM challenge. Other methods - table in supplementary. AUPRC to benchmark performance. DNAse-seq is the best predictor, RNA-seq and sequence-based features are not informative. Java implementation http://jstacs.de/index.php/Catchitt, predicted peaks for 32 transcription factors in 22 primary cell types and tissues (682 total) BED hg19 files, conservative and relaxed predictions https://www.synapse.org/#!Synapse:syn11526239/wiki/497341
    - Keilwagen, Jens, Stefan Posch, and Jan Grau. “Accurate Prediction of Cell Type-Specific Transcription Factor Binding.” Genome Biology 20, no. 1 (December 2019). https://doi.org/10.1186/s13059-018-1614-y.

- `hTFtarget` - database of TF-gene target regulations from >7K human ChIP-seq experiments. http://bioinfo.life.hust.edu.cn/hTFtarget/

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

