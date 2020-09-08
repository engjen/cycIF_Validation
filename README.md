# cycIF_Validation
*Code and data for validation of CyCIF multiplex imaging platform*

Multiplex imaging technologies are increasingly used for single-cell phenotyping and spatial characterization of tissues; however, quantitative, reproducible analysis is a technical and computational challenge. We developed an open-source python-based image analysis tool, mplex-image, to achieve fully-reproducible multiplex image visualization and analysis. We deploy this tool in the accompanying Jupyter notebooks to validate specificity, sensitivity, reproducibility and normalization of the multiplex imaging platform cyclic immunofluorescence (CyCIF). 

Through our work, we learned general principles of antibody staining performance, signal removal and background removal, summarized below:

**1. Antibody staining**
 - CyCIF method shows differences in dynamic range, but [similar signal-to-background as standard IF.](https://github.com/engjen/cycIF_Validation/blob/master/Fig1_SinglevsCyclic_analysis_44290.ipynb)
 - Antibodies applied ealier in the panel are brighter, while [later antibodies can show decreased signal.](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_EarlyvsLate_K157_analysis.ipynb)
 - Antibodies applied earlier may show non-specific staining; these [non-specific staining artifacts  are abrogated by applying antibodies later.](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_OrderOptimization_K154vsK175_analysis.ipynb)
 - Technical replicates show differences in dynamic range, but [replicates show similar signal-to-background ratios.](https://github.com/engjen/cycIF_Validation/blob/master/Fig3_TMAReplicates_analysis.ipynb)
 - Given a stable signal-to-background ratio, [batch normalization](https://github.com/engjen/cycIF_Validation/blob/master/Fig3_Normalize_JE-TMA-reps_cluster_analysis.ipynb) can be performed by dividing each marker's signal by its background. [We show this improves unsupervised cluster analysis](https://github.com/engjen/cycIF_Validation/blob/master/Fig3_Normalize_JE-TMA-reps_plot.ipynb)
 
 **2. Signal removal using hydrogen peroxide**
 - [Inceased concentrations](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_Quenching_analysis.ipynb) above 3% hydrogen peroxide do not improve speed of signal removal.
 - [Increased incubation times](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_Quenching_analysis.ipynb) of improve signal removal somewhat, but do not result in complete signal removal.
 - [Increased heat](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_Quenching_analysis.ipynb) of quenching solution results in complete signal removal but must be balanced with increased tissue loss.
 
 **3. Background autofluorescence removal**
 - A single quenching of 3% H2O2 applied for 15 - 30 minutes dramatically reduced tissue autofluorescence and this ["pre-quenching" step should be performed before staining.](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_Quenching_Single_Cell.ipynb)
 - Additional rounds of quenching [decrease the bright autofluorescent cells linearly](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_Quenching_Single_Cell.ipynb), while most cells show little additional decrease.
 - Therefore, we recommend taking a blank image after pre-quenching and and blank image after all stainiing has completed, and [combining these two images (weighted by round) for autofluorescence subtraction.](https://github.com/engjen/cycIF_Validation/blob/master/Fig2_44290-146_subtractAF_scale_by_round.ipynb)
