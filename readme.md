# Liver Cell Metabolism Simulation (Liver cancer cells excrete glutamate of cytosolic but not mitochondrial origin)

- Models and scripts for simulations of growth of HepG2 cells.

- Abstract:
Many cancer cells consume glutamine at high rates; counterintuitively, they simultaneously excrete glutamate, the first intermediate in glutamine metabolism. Glutamine consumption has been linked to replenishment of TCA intermediates and synthesis of ATP, but the reason for glutamate excretion is unclear. Here we dynamically profile the uptake and excretion fluxes of a liver cancer cell line (HepG2) and use genome-scale metabolic modeling for in-depth analysis. We find that up to 30% of the glutamine is metabolized in the cytosol, primarily for nucleotide synthesis, producing cytosolic glutamate. We hypothesize that excreting glutamate helps the cell to increase the nucleotide synthesis rate to sustain growth. Indeed, we show experimentally that partial inhibition of glutamate excretion reduces cell growth. Our integrative approach thus links glutamine addiction to glutamate excretion in cancer and points towards potential drug targets.

- KeyWords:

**Utilisation:** maximising growth, predictive simulation, experimental data integration; **Model Source:** HMR2.00; **Taxonomy:** Homo sapiens  **Condition:** HepG2 different glucose concentrations 

- Reference: Manuscript

- Pubmed ID: NaN

- Last update: 2020-03-05




This repository is administered by @avlant, Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology


## Installation

### Required Software:

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  * You need a functional Matlab installation of **Matlab_R_2015_b** (MATLAB 7.3 and higher)
  * The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)
  * Optional for graphviz plots (e.g figure 2C). Python 3 with the graphviz package.

### Run
To reproduce the figures of the paper run corresponding files (e.g. figure1a.m for figure 1A)
