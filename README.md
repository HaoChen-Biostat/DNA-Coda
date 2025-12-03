# DNA-Coda
Differential Network Analysis for Compositional Data with Applications to Metagenomics and Single Cell Sequencing Studies

DNA_Coda The code implements the proposed method (DNA_Coda) in "Differential Network Analysis for
Compositional Data with Applications to Metagenomics and Single Cell Sequencing Studies" by Chencheng Ma, Hao Chen, Yijuan Hu, and Xiang Zhan.

Compositional data are frequently encountered in metagenomics and single cell sequencing experiments.  Understanding condition-specific variations among compositional component-component interactions—represented as differential network matrix—is a fundamental task in these studies.  This article introduces a rigorous specification of the compositional differential network matrix and its basis counterpart, with the latter proven to be asymptotically identifiable under appropriate regularity conditions. Leveraging this connection, we propose DNA-Coda: a direct estimation method targeting the basis differential network matrix without intermediate estimation of individual precision matrices. We derive convergence rates for the DNA-Coda estimator and provide theoretical guarantees for support recovery consistency. The advantages of DNA-Coda over existing methods are illustrated by both simulation studies and application examples to metagenomics and single cell sequencing data analysis.
Keywords: Cell-cell communication; Compositional data analysis/Coda; Differential network analysis; Metagenomics; Precision matrix.

## Descriptions

The `simulation` folder contains files for reproducing the simulation studies:

The `DNA_Coda` folder: Files implementing the **proposed method (DNA_Coda)** and its oracle estimator
- `dpm.c`: C code for the implementation of DNA_Coda.(Linux: dpm.so, Windows: dpm.dll)
- `DNA_Coda.R`,: R wrapper function for dpm.c  to implement differential network matrix estimation for compositional data
- `powerlaw data generation.R`: functions to generate differential network matrix(Power-law graph)
- `Summary_roc_norm.R`: functions to summarize the results of support recovery(TPR, FPR)  and estimation accuracy( matrix norms)
- `example.R`: example usage, including DNA_Coda and its oracle

