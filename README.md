#  Integrative Multi-Omics Analysis of Chemical Exposure in *Daphnia magna*

##  Overview

This project investigates the molecular effects of river-derived polar organic pollutants on *Daphnia magna* using an integrative multi-omics approach. Transcriptomics, metabolomics (positive and negative ionisation modes), and environmental chemical exposure data were analysed to identify biomarkers, biological pathways, and molecular responses associated with pollutant exposure.

---

##  Objectives

- Identify biomarkers associated with river organic pollutants.
- Compare molecular responses across different sampling locations and exposure concentrations.
- Integrate transcriptomics and metabolomics to uncover coordinated biological responses.
- Perform pathway enrichment and network-based analyses to understand affected biological processes.

---

##  Experimental Design

- **Model organism:** *Daphnia magna*
- **Sampling locations:** 12 river sites (D01–D12)
- **Exposure levels:** Control, 1× environmental concentration, and 10× environmental concentration
- **Omics data:**
  - RNA-seq
  - Positive-mode metabolomics
  - Negative-mode metabolomics
  - Chemical exposure profiles

---

##  Analysis Workflow

1. Data preprocessing
2. Feature selection
3. Iterative Random Forest (iRF)
4. Multi-omics integration
5. Gene–metabolite correlation analysis
6. WGCNA network analysis
7. Functional annotation
8. Pathway enrichment
9. Biological interpretation

---

##  Computational Methods

- Iterative Random Forest (iRF)
- Weighted Gene Co-expression Network Analysis (WGCNA)
- RGCCA
- Pearson Correlation Analysis
- KEGG Annotation
- HMDB Annotation
- g:Profiler Pathway Enrichment
- Cytoscape Network Analysis

---

##  Key Results

- Identified candidate biomarkers associated with environmental pollutants.
- Integrated transcriptomics and metabolomics to reveal coordinated molecular responses.
- Constructed co-expression networks highlighting biologically relevant modules.
- Performed functional annotation and pathway enrichment to identify detoxification and stress-response pathways.

---

##  Technologies

- R
- Python
- WGCNA
- RGCCA
- gprofiler2
- Cytoscape
- DESeq2
- KEGG
- HMDB

---

##  Repository Contents

| File | Description |
|------|-------------|
| `annotation_enrichment_colab_final.R` | Annotation and pathway enrichment pipeline |
| `combined_top_features.csv` | Top features identified during analysis |
| `gene_metabolite_correlations.csv` | Gene–metabolite correlation results |
| `DATA_ACCESS.md` | Information regarding dataset availability |

---

##  Data Availability

The original transcriptomics, metabolomics, chemical exposure, and annotation datasets are **not included** in this repository because they are subject to intellectual property restrictions.

Only analysis code, processed summary results, and project documentation are provided.

---

##  Future Improvements

- Interactive visualisation dashboard
- Additional multi-omics integration approaches
- Validation using external datasets
- Expanded biological pathway analysis

---

##  Author

**Prithvi Athreya Jagadish**

MSc Bioinformatics  
University of Birmingham

---

##  License

This repository is intended for educational and research purposes.
