# Codes for Manuscript

Reproducible codes used for conducting analysis and generating figures and tables of the manuscript [Normalization and de-noising of single-cell Hi-C data with BandNorm and scVI-3D. BioRxiv (2021). Ye Zheng*, Siqi Shen* and Sündüz Keleş. * contribute equally.](https://www.biorxiv.org/content/10.1101/2021.03.10.434870v1)

The software source codes are deposited at

- **BandNorm** (R package): https://github.com/sshen82/BandNorm

- **scVI-3D** (Python tool): https://github.com/yezhengSTAT/scVI-3D

The curated scHi-C data used in this manuscript are cleaned and organized at https://pages.stat.wisc.edu/~sshen82/bandnorm/ for easy access and direct usage in BandNorm and scVI-3D. This online data repository contains five scHi-C data sets, all at 1Mb. The dataset name, reference and raw data source are summarized below:

- **Ramani2017**: GEO accession number GSE84920. [Ramani, Vijay, et al. "Massively multiplex single-cell Hi-C." Nature methods 14.3 (2017): 263-266.](https://www.nature.com/articles/nmeth.4155)

- **Lee2019**: GEO accession number GSE130711. [Lee, Dong-Sung, et al. "Simultaneous profiling of 3D genome structure and DNA methylation in single human cells." Nature methods 16.10 (2019): 999-1006.](https://www.nature.com/articles/s41592-019-0547-z)

- **Li2019**: GEO accession number GSE119171. [Li, Guoqiang, et al. "Joint profiling of DNA methylation and chromatin architecture in single cells." Nature methods 16.10 (2019): 991-993.](https://www.nature.com/articles/s41592-019-0502-z)

- **Kim2020**: https://noble.gs.washington.edu/proj/schic-topic-model/. [Kim, Hyeon-Jin, et al. "Capturing cell type-specific chromatin compartment patterns by applying topic modeling to single-cell Hi-C data." PLoS computational biology 16.9 (2020): e1008173.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008173)

- Tan2021**: GEO accession number GSE162511. [Tan, Longzhi, et al. "Changes in genome architecture and transcriptional dynamics progress independently of sensory experience during post-natal brain development." Cell 184.3 (2021): 741-758.](https://www.sciencedirect.com/science/article/pii/S0092867420317542)

The metadata summaries that contains the cell barcode, cell type, sequencing depth and sparsity for each data set are also provided at https://pages.stat.wisc.edu/~sshen82/bandnorm/Summary/.

In this repository, the scripts are classified into three subfolders based on the purpose. 

- **RunAlgorithms** folder contains scripts to run BandNorm, scVI-3D, scHi-C normalization baseline methods that we created for comparison and other publically available tools including [scHiCluster](https://github.com/zhoujt1994/scHiCluster), [scHiC Topics](https://github.com/khj3017/schic-topic-model) and [Higashi](https://github.com/ma-compbio/Higashi). The outputs of those scripts provide lower-dimensional embeddings and normalized scHi-C counts, which will be used for the downstream analysis and generating figures.

- **AnalysisCodes** folder contains analysis related scripts to implement evaluation, comparison and biological discoveries.

- **FigureGeneration** folder contains scripts that generate the figures in the manuscript.

For questions, please email Siqi Shen (sshen82@wisc.edu) or submit new issue at this repository.

This BandNorm_and_scVI-3D_manuscript repository is protected by GNU GPL 3 license. Please refer to LICENSE.md for using codes in this repository for your own purposes.
