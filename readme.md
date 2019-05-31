# Liver Cell Metabolism Simulation (Liver cancer cells excrete glutamate of cytosolic but not mitochondrial origin)

- Models and scripts for simulations of growth of HepG2 cells.

- Abstract:
Altered amino acid metabolism is an emerging characteristic of cancer. We measured 
the exchange fluxes of two liver cancer cell lines (HepG2 and Huh7) under different
nutrient conditions. The amino acids glutamate, cysteine, arginine, phenylalanine
and the branched chain amino acids (BCAA) were consumed at rates exceeding the 
requirement for protein synthesis, whilst alanine and glutamate were excreted. We
investigated the intracellular fluxes using a genome scale metabolic model, and 
found that mitochondrial glutamate is consumed whilst cytosolic glutamate is
excreted. Cytoplasmic glutamate originates primarily from deamination of BCAA and
transamination of glutamine for biosynthesis of nucleotides. We find that inhibiting
the glutamate transporter reduces growth of liver cancer cells. Based on our analysis
we discuss a number of other potential drug targets.

- KeyWords:

**Utilisation:** maximising growth, predictive simulation, experimental data integration; **Model Source:** HMR2.00; **Taxonomy:** Homo sapiens  **Condition:** HepG2/Huh7 different glucose concentrations 

- Reference: Manuscript

- Pubmed ID: NaN

- Last update: 2018-09-27




This repository is administered by name @avlant, Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology


## Installation

### Required Software:

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  * You need a functional Matlab installation of **Matlab_R_2015_b** (MATLAB 7.3 and higher)
  * The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)
  * Optional for graphviz plots (e.g figure 2C). Python 3 with the graphviz package.

### Run
To reproduce the figures of the paper run corresponding files (e.g. figure1a.m for figure 1A)