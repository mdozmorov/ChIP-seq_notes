# Notes on ChIP-seq and other-seq-related tools

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

ChIP-seq, ATAC-seq related tools and genomics data analysis resources. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Databases](#databases)
  - [Motif DBs](#motif-dbs)
- [ChIP-seq](#chip-seq)
  - [ChIP-seq pipelines](#chip-seq-pipelines)
    - [Normalization](#normalization)
    - [CUT&RUN](#cutrun)
  - [Quality control](#quality-control)
  - [Peaks](#peaks)
  - [Enhancers](#enhancers)
  - [Visualization](#visualization)
  - [Motif enrichment](#motif-enrichment)
  - [Differential peak detection](#differential-peak-detection)
  - [Enrichment](#enrichment)
  - [Interpretation](#interpretation)
  - [Excludable](#excludable)
- [DNAse-seq](#dnase-seq)
- [ATAC-seq](#atac-seq)
  - [ATAC-seq pipelines](#atac-seq-pipelines)
- [Histone-seq](#histone-seq)
  - [Broad peak analysis](#broad-peak-analysis)
- [Technology](#technology)
- [Machine learning](#machine-learning)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Databases

- [UniBind database](https://unibind.uio.no/) - TFBS predictions of approx. 56 million TFBSs with experimental and computational support for direct TF-DNA interactions for 644 TFs in > 1000 cell lines and tissues. Processed approx. 10,000 public ChIP-seq datasets from nine species using [ChIP-eat](https://bitbucket.org/CBGR/chip-eat/src/master/). ChIP-eat combines both computational (high PWM score) and experimental (centrality to ChIP-seq peak summit) support to find high-confidence direct TF-DNA interactions in a ChIP-seq experiment-specific manner, uses the DAMO tool. Input data - ReMap 2018 and GTRD. Robust and permissive collections. Over 197,000 Cis-regulatory modules. [Downloads](https://unibind.uio.no/downloads/) of BED, FASTA, PWMs, [Tracks for the UCSC GenomeBrowser](https://unibind.uio.no/genome-tracks/), [API](https://unibind.uio.no/api/), [Enrichment analysis, online](https://unibind.uio.no/enrichment/) with or without background, differential enrichment. [UniBind Enrichment BitBucket](https://bitbucket.org/CBGR/unibind_enrichment/src/master/). <details>
    <summary>Paper</summary>
    Puig, Rafael Riudavets, Paul Boddie, Aziz Khan, Jaime Abraham Castro-Mondragon, and Anthony Mathelier. “UniBind: Maps of High-Confidence Direct TF-DNA Interactions across Nine Species” BMC Genomics, (December 2021) https://doi.org/10.1186/s12864-021-07760-6

    Gheorghe, Marius, Geir Kjetil Sandve, Aziz Khan, Jeanne Chèneby, Benoit Ballester, and Anthony Mathelier. “A Map of Direct TF–DNA Interactions in the Human Genome.” Nucleic Acids Research 47, no. 4 (February 28, 2019): e21–e21. https://doi.org/10.1093/nar/gky1210
</details>

- [ReMap](https://remap.univ-amu.fr/) is an integrative analysis of Homo sapiens, Mus musculus and Arabidopsis thaliana transcriptional regulators from DNA-binding experiments such as ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA). Human hg38 and Arabidopsis TAOR10. All peaks, non-redundant peaks, cis-Regulatory Modules. [GitHub](https://github.com/remap-cisreg). [Download](https://remap.univ-amu.fr/download_page) genomic coordinates. <details>
    <summary>Paper</summary>
     Chèneby, Jeanne, Zacharie Ménétrier, Martin Mestdagh, Thomas Rosnet, Allyssa Douida, Wassim Rhalloussi, Aurélie Bergon, Fabrice Lopez, and Benoit Ballester. “[ReMap 2020: A Database of Regulatory Regions from an Integrative Analysis of Human and Arabidopsis DNA-Binding Sequencing Experiments](https://doi.org/10.1093/nar/gkz945).” Nucleic Acids Research, October 29, 2019

     Hammal, Fayrouz, Pierre de Langen, Aurélie Bergon, Fabrice Lopez, and Benoit Ballester. “ReMap 2022: A Database of Human, Mouse, Drosophila and Arabidopsis Regulatory Regions from an Integrative Analysis of DNA-Binding Sequencing Experiments.” Nucleic Acids Research, November 9, 2021, gkab996. https://doi.org/10.1093/nar/gkab996.
</details>

- [ADASTRA](https://adastra.autosome.org) -  the database of Allelic Dosage-corrected Allele-Specific human Transcription factor binding sites (over 500K sites across 1073 human TFs and 649 cell types, reprocessed data from [GTRD](#gtrd), pipeline at [GitHub](https://github.com/autosome-ru/ADASTRA-pipeline)) at nearly 270K SNPs. Background Allele Dosage (BAD) maps. Many SNPs overlap eQTLs. <details>
    <summary>Paper</summary>
    Abramov, Sergey, Alexandr Boytsov, Daria Bykova, Dmitry D. Penzar, Ivan Yevshin, Semyon K. Kolmykov, Marina V. Fridman, et al. “Landscape of Allele-Specific Transcription Factor Binding in the Human Genome.” Nature Communications 12, no. 1 (December 2021): 2751. https://doi.org/10.1038/s41467-021-23007-0.
</details>

- [ANANASTRA](https://ananastra.autosome.org/) - ANnotation and enrichment ANalysis of Allele-Specific TRAnscription factor binding at SNPs. Annotates a given list of SNPs with allele-specific binding events across a wide range of transcription factors and cell types using [ADASTRA](#adastra). Enrichment analysis of SNPs in cell type-specific TFBSs (Fisher's exact, one-sided). [API](https://ananastra.autosome.org/api/v4/). <details>
    <summary>Paper</summary>
    Boytsov, Alexandr, Sergey Abramov, Ariuna Z Aiusheeva, Alexandra M Kasianova, Eugene Baulin, Ivan A Kuznetsov, Yurii S Aulchenko, et al. “ANANASTRA: Annotation and Enrichment Analysis of Allele-Specific Transcription Factor Binding at SNPs.” Nucleic Acids Research, April 21, 2022, gkac262. https://doi.org/10.1093/nar/gkac262.
</details>

- [Catchitt](http://jstacs.de/index.php/Catchitt) - method for predicting TFBSs, leader of ENCODE-DREAM challenge. Other methods - table in supplementary. AUPRC to benchmark performance. DNAse-seq is the best predictor, RNA-seq and sequence-based features are not informative. Java implementation, predicted peaks for 32 transcription factors in 22 primary cell types and tissues (682 total) BED hg19 files, conservative and relaxed predictions, [download](https://www.synapse.org/#!Synapse:syn11526239/wiki/497341). <details>
    <summary>Paper</summary>
    Keilwagen, Jens, Stefan Posch, and Jan Grau. “Accurate Prediction of Cell Type-Specific Transcription Factor Binding.” Genome Biology 20, no. 1 (December 2019). https://doi.org/10.1186/s13059-018-1614-y
</details>

- [RAEdb](http://www.computationalbiology.cn/RAEdb/index.php) - enhancer database. Enhancers identified from STARR-seq and MPRA studies. Epromoters - promoters containing enhancers. Human (hg38)/mouse (mm10) data, select cell lines. BED/FASTQ download. Links to EnhancerAtlas, VISTA, SuperEnhancer databases. <details>
    <summary>Paper</summary>
    Cai, Zena, Ya Cui, Zhiying Tan, Gaihua Zhang, Zhongyang Tan, Xinlei Zhang, and Yousong Peng. “RAEdb: A Database of Enhancers Identified by High-Throughput Reporter Assays.” Database: The Journal of Biological Databases and Curation 2019 (January 1, 2019). https://doi.org/10.1093/database/bay140.
</details>

- [ChIP-Atlas](http://chip-atlas.org/) - a large database and analysis suite of public ChIP-seq and DNAse-seq experiments (Over 76K experiments, SRA uniformly processed data). Analyses: **Visualization** of peaks in IGV browser, BED file download, **Target genes** identification, **Colocalization** of factors (antigens), **Enrichment analysis** - permutation enrichment of BED regions, with custom background possible. [GitHub](https://github.com/shinyaoki/chipatlas/tree/master/sh), [Documentation](https://github.com/inutano/chip-atlas/wiki). <details>
    <summary>Paper</summary>
    Oki, Shinya, Tazro Ohta, Go Shioi, Hideki Hatanaka, Osamu Ogasawara, Yoshihiro Okuda, Hideya Kawaji, Ryo Nakaki, Jun Sese, and Chikara Meno. “ChIP‐Atlas: A Data‐mining Suite Powered by Full Integration of Public ChIP‐seq Data” EMBO Reports, (December 2018) https://doi.org/10.15252/embr.201846255
</details>

- [TRRUST](https://www.grnpedia.org/trrust/) database (Transcriptional Regulatory Relationships Unraveled by Sentence-based Text mining). Over 8K regulatory interactions for 800 TFs in human, and over 6K interactions for 828 mouse TFs. Mouse and human TF regulatory networks overlap, complement each other. More information than in PAZAR, TFactS, TRED, TFe databases. [Download](https://www.grnpedia.org/trrust/downloadnetwork.php), TSV format. [Tools](https://www.grnpedia.org/trrust/Network_search_form.php): 1. Search a gene, 2. Enrichment of key regulators for query genes. <details>
    <summary>Paper</summary>
    Han, Heonjong, Jae-Won Cho, Sangyoung Lee, Ayoung Yun, Hyojin Kim, Dasom Bae, Sunmo Yang, et al. “TRRUST v2: An Expanded Reference Database of Human and Mouse Transcriptional Regulatory Interactions.” Nucleic Acids Research 46, no. D1 (January 4, 2018): D380–86. https://doi.org/10.1093/nar/gkx1013.
</details>

- <a name="gtrd">[GTRD](http://gtrd.biouml.org)</a> - transcription factor binding sites and data (ChIP-seq, ChIP-seo, DNAse-seq, MNase-seq, ATAC-seq, RNA-seq), uniformly processed, over 35K experiments. Seven species, TFs linked to [CIS-BP](#cisbp). All cell types are assigned onthology. Experiment search, processed data/peaks download (BED, bigBed, bigWig). <details>
    <summary>Paper</summary>
    Yevshin, Ivan, Ruslan Sharipov, Tagir Valeev, Alexander Kel, and Fedor Kolpakov. “GTRD: A Database of Transcription Factor Binding Sites Identified by ChIP-Seq Experiments.” Nucleic Acids Research 45, no. D1 (January 4, 2017): D61–67. https://doi.org/10.1093/nar/gkw951.

    Kolmykov, Semyon, Ivan Yevshin, Mikhail Kulyashov, Ruslan Sharipov, Yury Kondrakhin, Vsevolod J Makeev, Ivan V Kulakovskiy, Alexander Kel, and Fedor Kolpakov. “GTRD: An Integrated View of Transcription Regulation.” Nucleic Acids Research 49, no. D1 (January 8, 2021): D104–11. https://doi.org/10.1093/nar/gkaa1057.
</details>

- [Cistrome DB](http://cistrome.org/db/#/) - ChIP-seq peaks for TFs, histone modifications, DNAse/ATAC. Downloadable cell type-specific, hg38 BED files. [Toolkit](http://dbtoolkit.cistrome.org/) to answer questions like "What factors regulate your gene of interest?", "What factors bind in your interval?", "What factors have a significant binding overlap with your peak set?" <details>
    <summary>Paper</summary>
    Zheng R, Wan C, Mei S, Qin Q, Wu Q, Sun H, Chen CH, Brown M, Zhang X, Meyer CA, Liu XS. Cistrome Data Browser: expanded datasets and new tools for gene regulatory analysis. Nucleic Acids Res, 2018 Nov 20. https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky1094/5193328
    
    Mei S, Qin Q, Wu Q, Sun H, Zheng R, Zang C, Zhu M, Wu J, Shi X, Taing L, Liu T, Brown M, Meyer CA, Liu XS. Cistrome data browser: a data portal for ChIP-Seq and chromatin accessibility data in human and mouse. Nucleic Acids Res, 2017 Jan 4;45(D1):D658-D662. https://academic.oup.com/nar/article/45/D1/D658/2333932
</details>

- [CODEX ChIP-seq](http://codex.stemcells.cam.ac.uk/) - CODEX provides access to processed and curated NGS experiments, including ChIP-Seq (transcription factors and histones), RNA-Seq and DNase-Seq. Human, mouse. Download tracks, analyze correlations, motifs, compare between organisms, more. <details>
    <summary>Paper</summary>
    Sánchez-Castillo, Manuel and Ruau, David and Wilkinson, Adam C. and Ng, Felicia S.L. and Hannah, Rebecca and Diamanti, Evangelia and Lombard, Patrick and Wilson, Nicola K. and Gottgens, Berthold. "CODEX: a next-generation sequencing experiment database for the haematopoietic and embryonic stem cell communities" Nucleic Acids Research, Database Issue, September 2014 https://doi.org/10.1093/nar/gku895
</details>

- [hTFtarget](http://bioinfo.life.hust.edu.cn/hTFtarget/) - database of TF-gene target regulations from >7K human ChIP-seq experiments.

### Motif DBs

- <a name="cisbp">[CIS-BP](http://cisbp.ccbr.utoronto.ca/index.php)</a> (The Catalog of Inferred Sequence Binding Preferences) - database of inferred sequence binding preferences. DNA sequence preferences for >1,000 TFs encompassing 54 different DBD classes from 131 diverse eukaryotes. PBM microarray assays to analyze TF binding preferences. Closely related DBDs (70% Amino Acid identity) almost always have very similar DNA sequence preferences, enabling inference of motifs for approx. 34% of the 70,000 known or predicted eukaryotic TFs. Tools to scan single sequence for TF binding, two sequences for differential TF binding (including SNP effect scan), protein scan, motif scan. Bulk download of PWMs, protein sequences, TF information, logos. <details>
    <summary>Paper</summary>
    Weirauch, Matthew T., Ally Yang, Mihai Albu, Atina G. Cote, Alejandro Montenegro-Montero, Philipp Drewe, Hamed S. Najafabadi, et al. “Determination and Inference of Eukaryotic Transcription Factor Sequence Specificity.” Cell 158, no. 6 (September 2014): 1431–43. https://doi.org/10.1016/j.cell.2014.08.009.
</details>

- <a name="hocomoco">[HOCOMOCO](https://hocomoco11.autosome.org/)</a> (Homo sapiens comprehensive model collection) - TFBS models and PWMs. Human- and mouse-specific models. HOCOMOCO v11 contains binding models for 453 mouse and 680 human transcription factors and includes 1302 mononucleotide and 576 dinucleotide position weight matrices. Uniformly processed data from [GTRD](#gtrd), peaks called with four peak callers (). Used [ChIPMunk](https://autosome.org/ChIPMunk/) in four computational models, including using DNA shape. Added [MoLoTool](https://molotool.autosome.org/), a web app to scan DNA sequences for TFBSs with PWMs. One model per TF is manually selected. Twice as many models as in JASPAR. <details>
    <summary>Paper</summary>
    Kulakovskiy, Ivan V., Yulia A. Medvedeva, Ulf Schaefer, Artem S. Kasianov, Ilya E. Vorontsov, Vladimir B. Bajic, and Vsevolod J. Makeev. “HOCOMOCO: A Comprehensive Collection of Human Transcription Factor Binding Sites Models.” Nucleic Acids Research 41, no. D1 (January 1, 2013): D195–202. https://doi.org/10.1093/nar/gks1089.

    Kulakovskiy, Ivan V., Ilya E. Vorontsov, Ivan S. Yevshin, Anastasiia V. Soboleva, Artem S. Kasianov, Haitham Ashoor, Wail Ba-alawi, et al. “HOCOMOCO: Expansion and Enhancement of the Collection of Transcription Factor Binding Sites Models.” Nucleic Acids Research 44, no. D1 (January 4, 2016): D116–25. https://doi.org/10.1093/nar/gkv1249.

    Kulakovskiy, Ivan V, Ilya E Vorontsov, Ivan S Yevshin, Ruslan N Sharipov, Alla D Fedorova, Eugene I Rumynskiy, Yulia A Medvedeva, et al. “HOCOMOCO: Towards a Complete Collection of Transcription Factor Binding Models for Human and Mouse via Large-Scale ChIP-Seq Analysis.” Nucleic Acids Research 46, no. D1 (January 4, 2018): D252–59. https://doi.org/10.1093/nar/gkx1106.
</details>

- [SwissRegulon](https://swissregulon.unibas.ch/sr/) - a database of regulatory motifs (PWMs) across model organisms (prokaryots, eukaryots). Data partly comes from JASPAR and TRANSFAC, reprocessing of ChIP-seq experiments. GBrowse for browsing TFBSs. Other tools. <details>
    <summary>Paper</summary>
    Pachkov, Mikhail, Piotr J. Balwierz, Phil Arnold, Evgeniy Ozonov, and Erik van Nimwegen. “SwissRegulon, a Database of Genome-Wide Annotations of Regulatory Sites: Recent Updates.” Nucleic Acids Research 41, no. D1 (November 23, 2012): D214–20. https://doi.org/10.1093/nar/gks1145.
</details>

## ChIP-seq

- [ChIP-seq-analysis](https://github.com/crazyhottommy/ChIP-seq-analysis) - ChIP-seq analysis notes from Ming Tang

### ChIP-seq pipelines

- [ChIP-AP](https://github.com/JSuryatenggara/ChIP-AP) - ChIP-seq analysis pipeline integrating multiple tools and peak callers (FastQC, Clumpify and BBDuk from the BBMap Suite, Trimmomatic, BWA, Samtools, deepTools, MACS2, GEM, SICER2, HOMER, Genrich, IDR, and the MEME-Suite). QC, cleanup, alignment, peak-calling, pathway analysis. High-confidence peaks based on overlaps by different peak callers. Input - single- or paired-end FASTQ files or aligned BAM files. Conda installable. Command line and GUI. [Documentation](https://github.com/JSuryatenggara/ChIP-AP/wiki/ChIP-AP-Guide). <details>
    <summary>Paper</summary>
    Suryatenggara, Jeremiah, Kol Jia Yong, Danielle E. Tenen, Daniel G. Tenen, and Mahmoud A. Bassal. "ChIP-AP: an integrated analysis pipeline for unbiased ChIP-seq analysis." Briefings in Bioinformatics 23, no. 1 (January 2022) https://doi.org/10.1093/bib/bbab537
</details>

- Regulatory Genomics Toolbox: Python library and set of tools for the integrative analysis of high throughput regulatory genomics data. http://www.regulatory-genomics.org, https://github.com/CostaLab/reg-gen
    - HINT (Hmm-based IdeNtification of Transcription factor footprints) is a framework for detection of DNA footprints from DNase-Seq and histone modification ChIP-Seq data.
    - Motif Analysis tools allows the search of motifs with binding sites enriched in particular genomic regions
    - ODIN and THOR are HMM-based approaches to detect and analyse differential peaks in pairs of ChIP-seq data. 
    - RGT-Viz is a collection of tests for association analysis and tools for visualizaiton of genomic data such as files in BED and BAM format
    - Triplex Domain Finder (TDF) statistically characterizes the triple helix potential of RNA and DNA regions.

- ENCODE3 pipeline v1 specifications, https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq

- [CHIPS](https://github.com/liulab-dfci/chips) - A Snakemake pipeline for quality control and reproducible processing of chromatin profiling data (ChIP-seq, ATAC-seq). Alignment, extensive QC, peak calling, downstream analysis (annotation, motif finding, putative targets). Generates an HTML report, plotly interactive plots. Distributed as a Conda recipe.
    - Taing, Len, Gali Bai, Clara Cousins, Paloma Cejas, Xintao Qiu, Zachary T. Herbert, Myles Brown, et al. “[CHIPS: A Snakemake Pipeline for Quality Control and Reproducible Processing of Chromatin Profiling Data](https://doi.org/10.12688/f1000research.52878.1).” F1000Research, (June 30, 2021)

- `AQUAS` ChIP-seq processing pipeline - The AQUAS pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline specifications, https://github.com/kundajelab/chipseq_pipeline

- `Crunch` - Completely Automated Analysis of ChIP-seq Data, http://crunch.unibas.ch/, https://www.biorxiv.org/content/early/2016/03/09/042903

- `ChiLin` - QC, peak calling, motif analysis for ChIP-seq and DNAse-seq data used by CistromeDb. References to other tools. https://github.com/cfce/chilin
    - Qin, Qian, Shenglin Mei, Qiu Wu, Hanfei Sun, Lewyn Li, Len Taing, Sujun Chen, et al. “ChiLin: A Comprehensive ChIP-Seq and DNase-Seq Quality Control and Analysis Pipeline.” BMC Bioinformatics 17, no. 1 (October 3, 2016): 404. https://doi.org/10.1186/s12859-016-1274-4.

- `ChIP-eat` - a pipeline for aligning reads, calling peaks, predicting TFBSs. https://bitbucket.org/CBGR/chip-eat/src/master/
    - Gheorghe, Marius, Geir Kjetil Sandve, Aziz Khan, Jeanne Chèneby, Benoit Ballester, and Anthony Mathelier. “A Map of Direct TF–DNA Interactions in the Human Genome.” Nucleic Acids Research 47, no. 4 (February 28, 2019): e21–e21. https://doi.org/10.1093/nar/gky1210.

- `ChIPLine` - a pipeline for ChIP-seq analysis, https://github.com/ay-lab/ChIPLine

#### Normalization

- [CHIPIN](https://github.com/BoevaLab/CHIPIN) - ChIP-seq Intersample Normalization using gene expression. Assumption - non-differential genes should have non-differential peaks.

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

- [CUT&RUNTools 2.0](https://github.com/fl-yu/CUT-RUNTools-2.0) - extended functionality to handle single-cell data, data normalization, peak calling (MACS2, SEACR), dimensionality reduction (Latent Semantic Indexing), downstream functional analysis.
    - Yu, Fulong, Vijay G Sankaran, and Guo-Cheng Yuan. “[CUT&RUNTools 2.0: A Pipeline for Single-Cell and Bulk-Level CUT&RUN and CUT&Tag Data Analysis](https://doi.org/10.1093/bioinformatics/btab507),” Bioinformatics, 09 July 2021

- [CUT&RUNTools](https://bitbucket.org/qzhudfci/cutruntools/src/master/) - a pipeline to fully process CUT&RUN data and identify protein binding and genomic footprinting from antibody-targeted primary cleavage data. Implemented in R, Python, Bach, runs under the SLURM job submission. At the core, creates a cut matrix of from enzyme cleavage data. Compared with Atactk and Centipede. (Tested, didn't work)
    - Zhu, Qian. “[CUT&RUNTools: A Flexible Pipeline for CUT&RUN Processing and Footprint Analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4),” 2019, 12. 

### Quality control

- `phantompeakqualtools` - This package computes informative enrichment and quality measures for ChIP-seq/DNase-seq/FAIRE-seq/MNase-seq data. It can also be used to obtain robust estimates of the predominant fragment length or characteristic tag shift values in these assays. https://github.com/kundajelab/phantompeakqualtools

### Peaks

- [LanceOtron](https://github.com/LHentges/LanceOtron) - deep learning-based peak caller from TF and histone ChIP-seq, ATAC-seq, DNAse-seq. Input - bigWig coverage file (+input, if available). Image recognition using wide and deep model (logistic regression producing enrichment scores, CNN, multilayer perceptron, Fig. 1c, Methods). Trained on hand-labeled data. Outperforms MACS2. Visualization using [MLV genome visualization](https://github.com/Hughes-Genome-Group/mlv) software. [Website](https://lanceotron.molbiol.ox.ac.uk/) with videos, documentation. <details>
    <summary>Paper</summary>
    Hentges, Lance D., Martin J. Sergeant, Damien J. Downes, Jim R. Hughes, and Stephen Taylor. "LanceOtron: a deep learning peak caller for ATAC-seq, ChIP-seq, and DNase-seq." bioRxiv (2021). https://doi.org/10.1101/2021.01.25.428108
</details>

- [epic2](https://github.com/biocore-ntnu/epic2) - diffuse ChIP-seq peak caller, Cython reimplementation of SICER, 30X times faster, 1/7 memory use. Available on Conda and [GitHub](https://github.com/biocore-ntnu/epic2)
    - Stovner, Endre Bakken. “[Epic2 Efficiently Finds Diffuse Domains in ChIP-Seq Data](https://doi.org/10.1093/bioinformatics/btz232),” Bioinformatics. 2019 Nov 1

- `Genrich` - Detecting sites of genomic enrichment in ChIP-seq and ATAC-seq. https://github.com/jsh58/Genrich, unpublished but highly tested and recommented, https://informatics.fas.harvard.edu/atac-seq-guidelines.html

- `mosaics` - This package provides functions for fitting MOSAiCS and MOSAiCS-HMM, a statistical framework to analyze one-sample or two-sample ChIP-seq data of transcription factor binding and histone modification. https://bioconductor.org/packages/release/bioc/html/mosaics.html

- `RSEG` - ChIP-seq broad domain analysis. http://smithlabresearch.org/software/rseg/
    - Song, Qiang, and Andrew D Smith. “Identifying Dispersed Epigenomic Domains from ChIP-Seq Data.” Bioinformatics 27, no. 6 (2011): 870–871.

- `triform` - finds enriched regions (peaks) in transcription factor ChIP-sequencing data. https://bioconductor.org/packages/release/bioc/html/triform.html

### Enhancers

- [ROSE](https://bitbucket.org/young_computation/rose) - rank-ordering of super-enhancers using H3K27ac ChIP-seq data, by the Young lab.

- [LILI](https://github.com/BoevaLab/LILY) - a pipeline by Boeva lab for detection of super-enhancers using H3K27ac ChIP-seq data, which includes explicit correction for copy number variation inherent to cancer samples. The pipeline is based on the ROSE algorithm originally developed by the Young lab.

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


### Motif enrichment

- [memes](https://github.com/snystrom/memes/) - an R package interfacing MEME suite (DREME, ME, FIMO, TOMTOM). Using universalmotif_df R/Bioconductor object to make results compatible across tools. De novo motif discovery, differential motifs, known motif enrichment analysis. Visualization capabilities. Case example on ChIP-seq peaks in Drosophila wing development. Requires installation of MEME suite. Docker container with RStudio and everything configured. [Bioconductor](https://bioconductor.org/packages/memes/), [pkgdown website](https://snystrom.github.io/memes-manual/index.html).  <details>
    <summary>Paper</summary>
    Nystrom, Spencer L, and Daniel J McKay. “Memes: A Motif Analysis Environment in R Using Tools from the MEME Suite.” PLOS COMPUTATIONAL BIOLOGY, n.d., 14.
</details>

- [ChEA3](https://maayanlab.cloud/chea3/) - transcription factor enrichment in gene lists. Six reference libraries of TF regulatory signatures (ARCHS4, ENCODE, GTeX, ReMap, Enrichr, Literature). Fisher's exact test. Outperform VIPER, DoRothEA, BART, TFEA.ChIP, oPOSSUM, MAGICACT. <details>
    <summary>Paper</summary>
    Keenan, Alexandra B, Denis Torre, Alexander Lachmann, Ariel K Leong, Megan L Wojciechowicz, Vivian Utti, Kathleen M Jagodnik, Eryk Kropiwnicki, Zichen Wang, and Avi Ma’ayan. “[ChEA3: Transcription Factor Enrichment Analysis by Orthogonal Omics Integration](https://doi.org/10.1093/nar/gkz446).” Nucleic Acids Research, (July 2, 2019)
</details>

- [TFEA.ChIP](https://bioconductor.org/packages/TFEA.ChIP/) - R package for transcription factor enrichment of gene lists (hypergeometric and GSEA) using experimental ChIP-seq datasets (ENCODE, GEO). Tested on known signatures, compared with two PWM-based and ChIP-based, performs comparably or better.
    - Puente-Santamaria, Laura, Wyeth W. Wasserman, and Luis Del Peso. "[TFEA.ChIP: a tool kit for transcription factor binding site enrichment analysis capitalizing on ChIP-seq datasets](https://doi.org/10.1093/bioinformatics/btz573)." Bioinformatics, (2019)

- <a name="pwmscan">[PWMScan](https://ccg.epfl.ch/pwmtools/pwmscan.php)</a> - web tool for scanning entire genomes with a position-specific weight matrix. Multiple genomes and assemblies hosted on the server. Multiple PWM collections for Eukaryotic DNA (JASPAR, [HOCOMOCO](#hocomoco), SwissRegulon, UniPROBE, [CIS-BP](#cisbp), from Jomla, Isakova publications) matrix_scan C program for matching PWMs. Compared with other motif scanning tools (PoSSuMseqrch, Patser, RSAT, STORM, HOMER), overlap >99%. Output - [BEDdetail](https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7) format. [Code](https://sourceforge.net/projects/pwmscan/). <details>
    <summary>Paper</summary>
    Ambrosini, Giovanna, Romain Groux, and Philipp Bucher. “PWMScan: A Fast Tool for Scanning Entire Genomes with a Position-Specific Weight Matrix.” Edited by John Hancock. Bioinformatics 34, no. 14 (July 15, 2018): 2483–84. https://doi.org/10.1093/bioinformatics/bty127.
</details>

- [gimmemotifs](https://github.com/vanheeringen-lab/gimmemotifs) - framework for TF motif analysis using an ensemble of motif predictors. `maelstrom` tool to detect differential motif activity between multiple different conditions. Includes manually curated database of motifs. Benchmark of 14 motif detection tools - Homer, MEME, BioProspector are among the top performing. Extensive analysis results. [Documentation](https://gimmemotifs.readthedocs.io). [Tweet with updates](https://twitter.com/svheeringen/status/1313778158999666688?s=20)
    - Bruse, Niklas, and Simon J. van Heeringen. “[GimmeMotifs: An Analysis Framework for Transcription Factor Motif Analysis](https://doi.org/10.1101/474403),” November 20, 2018

- `DECOD` - Differential motif finder. k-mer-based. http://gene.ml.cmu.edu/DECOD/
    - Huggins, Peter, Shan Zhong, Idit Shiff, Rachel Beckerman, Oleg Laptenko, Carol Prives, Marcel H. Schulz, Itamar Simon, and Ziv Bar-Joseph. “DECOD: Fast and Accurate Discriminative DNA Motif Finding.” Bioinformatics 27, no. 17 (September 1, 2011): 2361–67. https://doi.org/10.1093/bioinformatics/btr412.

- [Non-redundant TF motif matches genome-wide](https://www.vierstra.org/resources/motif_clustering#downloads) - Clustering of 2179 motif models. hg38/mm10 BED files download with coordinates

- [homerkit](https://github.com/slowkow/homerkit) - Read HOMER motif analysis output in R. 

- [LISA](http://lisa.cistrome.org/) - epigenetic Landscape In silico Subtraction analysis, enriched TFs and chromatin regulators in a list of genes.

- [Logolas](https://github.com/kkdey/Logolas) - R package for Enrichment Depletion Logos (EDLogos) and String Logos.

- [marge](https://robertamezquita.github.io/marge/) - API for HOMER in R for Genomic Analysis using Tidy Conventions, [GitHub](https://github.com/robertamezquita/marge)

- [motifStack](https://bioconductor.org/packages/motifStack/) - R package for plotting stacked logos for single or multiple DNA, RNA and amino acid sequence.

- [pyjaspar](https://github.com/asntech/pyjaspar) - A Pythonic interface to query and access JASPAR transcription factor motifs

- [RcisTarget](https://bioconductor.org/packages/RcisTarget/) - R package to identify transcription factor binding motifs enriched on a list of genes or genomic regions.

- [rGADEM](https://bioconductor.org/packages/rGADEM/) - R package for de novo motif discovery.


### Differential peak detection

- [normR](https://bioconductor.org/packages/normr/) - a Bioconductor R package, data-driven normalization and difference calling approach for ChIP-seq data. Models ChIP- and control read counts by binomial mixture model. One component models the background, the other models the signal. Can work without control.  
    - Helmuth, Johannes, et al. "[normR: Regime enrichment calling for ChIP-seq data](https://doi.org/10.1101/082263)." BioRxiv (2016)

- `csaw` - Detection of differentially bound regions in ChIP-seq data with sliding windows, with methods for normalization and proper FDR control. https://bioconductor.org/packages/release/bioc/html/csaw.html

- `DiffBind` - Differential Binding Analysis of ChIP-Seq Peak Data. https://bioconductor.org/packages/release/bioc/html/DiffBind.html

### Enrichment

- UniBind Enrichment Analysis, also differential enrichment. Input - BED file in hg38 version. [LOLA](https://bioconductor.org/packages/release/bioc/html/LOLA.html) as an enrichment and database engine.

### Interpretation

- [Cistrome-GO](http://go.cistrome.org/) - functional enrichment analysis of genes regulated by TFs in human and mouse. Solo mode (ChIP-seq peaks only) or ensemble mode (integrates ChIP-seq peaks and RNA-seq differentially expressed genes). Implementation of BETA method. MACS2 peaks, DESeq2 output. Gene-centric regulatory potential (RP) score (exponentially weighted by distance sum of peaks). Human (hg19/hg38), Mouse (mm9/mm10). <details>
    <summary>Paper</summary>
    Li, Shaojuan, Changxin Wan, Rongbin Zheng, Jingyu Fan, Xin Dong, Clifford A. Meyer, and X. Shirley Liu. “Cistrome-GO: A Web Server for Functional Enrichment Analysis of Transcription Factor ChIP-Seq Peaks.” Nucleic Acids Research 47, no. W1 (July 2, 2019): W206–11. https://doi.org/10.1093/nar/gkz332.
</details>

- `Toolkit for Cistrome Data Browser` - online tool to answer questions like:
    - What factors regulate your gene of interest?
    - What factors bind in your interval?
    - What factors have a significant binding overlap with your peak set?

- Chongzhi Zhang software page, http://faculty.virginia.edu/zanglab/software.htm
    - `BART` (Binding Analysis for Regulation of Transcription), a bioinformatics tool for predicting functional transcription factors (TFs) that bind at genomic cis-regulatory regions to regulate gene expression in the human or mouse genomes, given a query gene set or a ChIP-seq dataset as input. http://bartweb.uvasomrc.io/
    - `MARGE` (Model-based Analysis of Regulation of Gene Expression), a comprehensive computational method for inference of cis-regulation of gene expression leveraging public H3K27ac genomic profiles in human or mouse. http://cistrome.org/MARGE/
    - `MANCIE` (Matrix Analysis and Normalization by Concordant Information Enhancement), a computational method for high-dimensional genomic data integration. https://cran.r-project.org/web/packages/MANCIE/index.html
    - `SICER` (Spatial-clustering Identification of ChIP-Enriched Regions), a ChIP-Seq data analysis method. https://home.gwu.edu/~wpeng/Software.htm


### Excludable

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

- [awesome-atac-analysis](https://github.com/databio/awesome-atac-analysis) Awesome ATAC-seq analysis by Nathan Sheffield.

- [Benchmarking ATAC-seq peak calling](https://bigmonty12.github.io/peak-calling-benchmark) by Austin Montgomery

- ATAC-seq analysis considerations. Considering multiple workflows, settling on csaw-based. Normalization by library complexity (downsampling) is important. [Workflow](https://github.com/reskejak/ATAC-seq/blob/master/ATACseq_workflow.txt) and [GitHub](https://github.com/reskejak/ATAC-seq) with all scripts. <details>
    <summary>Paper</summary>
    Reske, Jake J., Mike R. Wilson, and Ronald L. Chandler. “ATAC-Seq Normalization Method Can Significantly Affect Differential Accessibility Analysis and Interpretation.” Epigenetics & Chromatin 13, no. 1 (December 2020): 22. https://doi.org/10.1186/s13072-020-00342-y.
</details>

### ATAC-seq pipelines

- [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline) - ATAC-seq and DNase-seq processing pipeline by Anshul Kundaje

- [HINT-ATAC](http://www.regulatory-genomics.org/hint/introduction/) - a footprinting method considering ATAC-seq protocol biases. Uses a position dependency model (PDM) to learn the cleavage preferences (Methods). Compared against three footprinting methods, DNase2TF, PIQ, Wellington. PDMs are crucial for correction of cleavage bias for ATAC-seq for all methods. Also improves correction for DNAse-seq data. Comparison of protocols, Omni-ATAC (best performance), Fast-ATAC. Part of [RGT, Regulatory Genomics Toolbox](https://github.com/CostaLab/reg-gen). [Tutorial](https://www.regulatory-genomics.org/hint/tutorial-differential-footprints-on-scatac-seq/). <details>
    <summary>Paper</summary>
    Li, Zhijian, Marcel H. Schulz, Thomas Look, Matthias Begemann, Martin Zenke, and Ivan G. Costa. “Identification of Transcription Factor Binding Sites Using ATAC-Seq.” Genome Biology, (December 2019). https://doi.org/10.1186/s13059-019-1642-2
</details>

- [HMMRATAC](https://github.com/LiuLabUB/HMMRATAC) - hidden Markov model for ATAC-seq to identify open chromatin regions. Parametric modeling of nucleosome-free regions and three nucleosomal reatures (mono-, di-, and tri-nucleosomes). First, train on 1000 auto-selected regions, then predict. Tested on "active promoters" and "strong enhancers" chromatin states (positive examples), and "heterochromatin" (negative examples). Compared with MACS2, F-seq. <details>
    <summary>Paper</summary>
    Tarbell, Evan D, and Tao Liu. “HMMRATAC: A Hidden Markov ModeleR for ATAC-Seq.” Nucleic Acids Research, June 14, 2019, gkz533. https://doi.org/10.1093/nar/gkz533
</details>

- ATAC-seq peak calling using MACS2: `macs2 callpeak --nomodel --nolambda -- keep-dup all --call-summits -f BAMPE -g hs`

- [ATACProc](https://github.com/ay-lab/ATACProc) - ATAC-seq processing pipeline

- [atacseq](https://github.com/nf-core/atacseq) - nf-core ATAC-seq peak-calling and differential analysis pipeline.

- [pepatac](https://github.com/databio/pepatac) - A modular, containerized pipeline for ATAC-seq data processing. [Examples and documentation](http://code.databio.org/PEPATAC/)

- [TOBIAS](https://github.com/loosolab/TOBIAS) - Transcription factor Occupancy prediction By Investigation of ATAC-seq Signal. 

## Histone-seq

Homer program ‘findPeaks’ with the style ‘histone’. Peaks within 1 kb were merged into a single peak. Broad peaks in H3K36me3, H3K27me3 and H3K9me3 were called using the Homer program ‘findPeaks’ with the options ‘-region –size 1000 –minDist 2500’. When Homer runs with these options, the initial sets of peaks were 1 kb wide and peaks within 2.5 kb were merged. 

- `DEScan2` - broad peak (histone, ATAC, DNAse) analysis (peak caller, peak filtering and alignment across replicates, creation of a count matrix). Peak caller uses a moving window and calculated a Poisson likelihood of a peak as compared to a region outside the window. https://bioconductor.org/packages/release/bioc/html/DEScan2.html
    - Righelli, Dario, John Koberstein, Nancy Zhang, Claudia Angelini, Lucia Peixoto, and Davide Risso. “Differential Enriched Scan 2 (DEScan2): A Fast Pipeline for Broad Peak Analysis.” PeerJ Preprints, 2018.

- `HMCan` and `HMCan-diff` - histone ChIP-seq peak caller (and differential) that accounts for CNV, also for CG bias. Hidden Markov Model to detect peak signal. Control-FREEC to detect CNV in ChIP-seq data. Outperforms others, CCAT second best. https://www.cbrc.kaust.edu.sa/hmcan/
    - Ashoor, Haitham, Aurélie Hérault, Aurélie Kamoun, François Radvanyi, Vladimir B. Bajic, Emmanuel Barillot, and Valentina Boeva. “HMCan: A Method for Detecting Chromatin Modifications in Cancer Samples Using ChIP-Seq Data.” Bioinformatics (Oxford, England) 29, no. 23 (December 1, 2013): 2979–86. https://doi.org/10.1093/bioinformatics/btt524.
    - Ashoor, Haitham, Caroline Louis-Brennetot, Isabelle Janoueix-Lerosey, Vladimir B. Bajic, and Valentina Boeva. “HMCan-Diff: A Method to Detect Changes in Histone Modifications in Cells with Different Genetic Characteristics.” Nucleic Acids Research 45, no. 8 (05 2017): e58. https://doi.org/10.1093/nar/gkw1319.

- `RSEG` - ChIP-seq analysis for identifying genomic regions and their boundaries marked by diffusive histone modification markers, such as H3K36me3 and H3K27me3, http://smithlabresearch.org/software/rseg/

### Broad peak analysis

- [EDD](https://github.com/CollasLab/edd) - Enriched Domain Detector, a ChIP-seq peak caller for detection of megabase domains of enrichment.

- [epic2](https://github.com/biocore-ntnu/epic2) - an ultraperformant reimplementation of SICER. It focuses on speed, low memory overhead and ease of use.
    - Stovner, Endre Bakken, and Pål Sætrom. "[epic2 efficiently finds diffuse domains in ChIP-seq data](https://doi.org/10.1093/bioinformatics/btz232)." Bioinformatics, (2019)

- [DEScan2](https://bioconductor.org/packages/DEScan2/) - Integrated peak and differential caller, specifically designed for broad epigenomic signals, R package.

## Technology

- **CLIP-seq** (cross-linking and immunoprecipitation) technology, detects sites bound by a protein to RNAs.Figure 1 - technology overview, Figure 2 - details of HITS-CLIP/iCLIP/irCLIP/eCLIP/PAR-CLIP/Proximity-CLIP. Computational analysis, Table 3 - peak detection software. Databases ([doRiNA](https://dorina.mdc-berlin.de/), [ENCORI](https://starbase.sysu.edu.cn/), [POSTAR3](http://111.198.139.65/)). <details>
    <summary>Paper</summary>
    Hafner, Markus, Maria Katsantoni, Tino Köster, James Marks, Joyita Mukherjee, Dorothee Staiger, Jernej Ule, and Mihaela Zavolan. "CLIP and complementary methods." Nature Reviews Methods Primers 1, no. 1 (2021): 1-23. https://doi.org/10.1038/s43586-021-00018-1
</details>

## Machine learning

- [maxATAC](https://github.com/MiraldiLab/maxATAC) - TFBS prediction from ATAC-seq (bulk and pseudobulk) in any cell type (whole genome, chromosome, or region). Deep dilated convolutional neural networks, bigWig and BED predictions of TFBSs. [Models avaliable for 127 human TFs](https://github.com/MiraldiLab/maxATAC_data) (h5 files). Outperforms baseline (average ChIP-seq signal, motif scanning) for most TFs and cell lines. AUPR is similar to the top performer in the ENCODE-DREAM in vivo TFBS prediction challenge (0.4). OMNI-ATAC-seq data for three cell lines, to be available. ATAC-seq scaling to signal per replicate to 20 million mapped reads (RP20M) and min-max normalized to 99th percentile signals. Python, separate functions  for each step (prepare, average, normalize, train, predict, benchmark, peaks, variants). [Tweet 1](https://twitter.com/tareian_it_up/status/1487614524492505090?s=20&t=1dQPuanBrvUUlP_g-Uo9jQ), [Tweet 2](https://twitter.com/EmilyMiraldi/status/1494414950848253953?s=20&t=1dQPuanBrvUUlP_g-Uo9jQ). <details>
    <summary>Paper</summary>
    Cazares, Tareian A, Faiz W Rizvi, Balaji Iyer, Xiaoting Chen, Michael Kotliar, Joseph A Wayman, Anthony Bejjani, et al. “MaxATAC: Genome-Scale Transcription-Factor Binding Prediction from ATAC-Seq with Deep Neural Networks.” Preprint. Bioinformatics, January 29, 2022. https://doi.org/10.1101/2022.01.28.478235.
</details>

- Segmentation and genome annotation (SAGA) algorithms review. Methods and tools for finding patterns from multiple ChIP-seq, histone-seq, etc. measures (Table 1). Hidden Markov Model (HMM), Dynamic Bayesian Network (DBN) algorithms. HMM intuition, math, solution algorithms. Visualization. Future work, challenges. <details>
    <summary>Paper</summary>
    Libbrecht, Maxwell W., Rachel C. W. Chan, and Michael M. Hoffman. “Segmentation and Genome Annotation Algorithms for Identifying Chromatin State and Other Genomic Patterns.” Edited by Tamar Schlick. PLOS Computational Biology 17, no. 10 (October 14, 2021): e1009423. https://doi.org/10.1371/journal.pcbi.1009423.
</details>


## Misc

- `covtobed` - a tool to generate BED coverage tracks from BAM files. https://github.com/telatin/covtobed

- [UCSC Genome Browser API](http://genome.ucsc.edu/goldenPath/help/api.html) to retrieve DNA sequence from coordinates. 
    - https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chrM;start=4321;end=5678
    - https://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:4336341,4336599