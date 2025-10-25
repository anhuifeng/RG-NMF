# RG-NMF
Identifying pancreatic cell type via robust graph-regularized non-negative matrix factorization


Precise profiling of highly heterogeneous cell types and their interactions is critical for understanding pancreatic physiology and pathological mechanisms. However, the inherent high-dimensionality, high noise, and dropout events in single-cell RNA sequencing data significantly impede accurate cell-type annotation, presenting a key bottleneck in identifying pancreatic cell types. To address this challenge, this paper proposed a robust graph-regularized non-negative matrix factorization (RG-NMF) model. A hybrid loss function applying mean squared error to non-zero values and Cauchy loss to zero-expression values was presented, mitigating dropout-induced perturbations during reconstruction. Through dual graph regularization, RG-NMF simultaneously preserved local manifold structures in both cell and gene spaces. Experimental validation across three pancreatic scRNA-seq datasets demonstrated that RG-NMF significantly enhances cell-type identification accuracy while maintaining robust cross-dataset stability.

Keywords: Cell type identification, Single-cell RNA sequencing, Non-negative matrix factorization, Graph-regularization

Matlab version: R2018b and later versions

