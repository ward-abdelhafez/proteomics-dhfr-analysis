# Proteomics Data Analysis: DHFR Knockout Study

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebook-orange.svg)](https://jupyter.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🔬 Project Overview

A comprehensive computational proteomics workflow demonstrating **quality control**, **peptide identification validation**, and **search parameter extraction** from mass spectrometry data.

**Biological Context:** DHFR (Dihydrofolate Reductase) knockout study in human cells  
**Dataset:** [PXD044513](https://www.ebi.ac.uk/pride/archive/projects/PXD044513) from PRIDE Archive  
**Instrument:** Orbitrap Tribrid MS (Thermo Fisher Scientific)  
**Search Platform:** Proteome Discoverer 2.2.0.388 (Sequest HT + Percolator)

---

## 📊 Key Findings

| Metric | Value | Assessment |
|--------|-------|------------|
| MS2 Spectra Acquired | 48,537 | Appropriate depth |
| High-Confidence PSMs | 34,094 | FDR < 1% |
| Unique Peptides | 21,720 | Rich coverage |
| **Identification Rate** | **70.2%** | **Excellent** ✅ |
| Mean Peptide Length | 14.2 aa | Optimal tryptic |
| Charge +2 Dominance | 75.3% | Expected profile |
| RT-Hydrophobicity Correlation | r = 0.45 | Good LC separation |

---

## 🚀 Quick Start

### Prerequisites

```bash
# Conda environment (recommended)
conda create -n proteomics python=3.9
conda activate proteomics

### Installation
# Clone repository
git clone https://github.com/ward-abdelhafez/proteomics-dhfr-analysis.git
cd proteomics-dhfr-analysis

# Install dependencies
pip install -r requirements.txt

### Running Analysis
jupyter lab
# Open notebooks in order:
# 1. notebooks/01_Data_Loading_QC.ipynb
# 2. notebooks/02_Parameter_Extraction.ipynb

### Repository Structure
proteomics-dhfr-analysis/
├── notebooks/          # Jupyter notebooks
├── data/              # Data directory
├── results/           # Analysis outputs
├── docs/              # Documentation
└── scripts/           # Utility scripts

### Technologies Used
pyteomics - MS data parsing
pandas - Data manipulation
numpy - Numerical operations
matplotlib/seaborn - Visualization
scipy - Statistical analysis

### Documentation
Full report: docs/Proteomics_Analysis_Complete_Updated.md

### License
MIT License - see LICENSE file.
