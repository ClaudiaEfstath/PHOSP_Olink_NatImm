Code relating to Liew et al., Nature Immunology 2023 

Penalized logistic regression

This script creates a penalised logistic regression (PLR) with LASSO regression, forest plots of results, and stability testing of results. The PLR is a logistic regression technique that allows analysis of high demensional data by imposing a pentaly function to shrink the coeffcient variables, this allows only the most signficant variable to be included in the final model. 

Network_nat_imm_revs - Network analysis

To further describe the results of the the PLR, the relationships of the top proteins identifed by the PLR were analysed with the qgraph package.
Data was restricted to the NPX results generated through Olink Explore 384 Inflammation panel. Script used to generate Fig. 2.


This code has been edited to remove senstive data. The full dataset can be accessed contacting phosp@leicester.ac.uk 
