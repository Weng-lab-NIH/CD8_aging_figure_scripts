# Single-cell RNA analysis of CD8 T cell aging in healthy individuals
These R notebooks contains all the analysis done for the single-cell data, starting from CellRanger output and going upto the final generation of each single-cell figure. 

The **preprocessing** notebook contains code for reading cellranger output, filtering, normalizing, integrating, and clustering.

**Notebook 1** contains code for differential expression and produces figures for the tSNE projection and differential expression heatmap in Figure 1.

**Notebook 2** contains code for generating linear models examining how the 9 CD8 subpopulations change over time, creating an expression heatmap for the 50 most-changed genes by age, and creating a network plot of functional changes associated with age and associated genes for Figure 2.

**Notebook 3_4** contains code for statistically testing whether genes are changing in number of expressing cells and/or average expression level per cell, and plotting representative examples. It also contains code for gene set enrichment analysis (GSEA), and the corresponding bubble plots for Figures 3 and 4.

**Notebook 5** contains code for ******both****** training and plotting the results of our age-prediction machine learning model on our training data, our validation data, and an external test dataset for Figure 5. 

**Notebook 6** compares the predicted ages of each CD8 T cell subset with chronological age, as well as mutation burden. Mutation burden was calculated using the [USCMD](https://github.com/Weng-lab-NIH/USCMD) method. 

**Notebook 7** Calculates the results of applying the age-prediction model to CD8+ T cells from HIV patients as well as CAR T Cells. 
