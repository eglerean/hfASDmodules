# hfASDmodules
Code used in **Glerean et al. 2015 "Reorganization of functionally connected brain subnetworks in high-functioning autism"** (submitted)

This is part of the code used for the article mentioned above. 

## Preprocessing and head motion quality control
The BraMiLa Matlab tools were used for further preprocessing and head motion quality control
https://git.becs.aalto.fi/bml/bramila/
specifically
- bramila_clean_signal.m (to further clean the fMRI data as described in Power et al. 2014 http://www.sciencedirect.com/science/article/pii/S1053811913009117
- bramila_diagnostics.m (for quality control)

## Graph-theoretical analysis (macro-level, mesoscopic, micro-level)
All graph-theoretical analysis was done in Python 2.7 using custom tools available at:
https://git.becs.aalto.fi/rmkujala/brainnets

The scripts are using the extensive library developed by the Complex Networks group at the Neuroscience and Biomedical Engineering department of Aalto University
https://git.becs.aalto.fi/complex-networks/verkko/tree/master
http://becs.aalto.fi/en/research/complex_networks/

## ABIDE dataset
The subfolder ABIDE contains
- 


## Statistics
The following scripts were used for statistics
- bramila_ttest_np.m 
-- Used to efficiently compute difference of the means using permutations. It uses Matlab parallel computing toolbox.
- bramila_mantel.m
-- Used to perform Mantel test (correlation between two distance/similarity matrices). 
- %% add here micro-level stuff


