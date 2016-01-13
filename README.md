# hfASDmodules #
Code used in **Glerean et al. 2015 "Reorganization of functionally connected brain subnetworks in high-functioning autism"** (in press: http://onlinelibrary.wiley.com/doi/10.1002/hbm.23084/abstract)

This is part of the code used for the article mentioned above. 

## Preprocessing and head motion quality control ##
The BraMiLa Matlab tools were used for further preprocessing and head motion quality control
https://git.becs.aalto.fi/bml/bramila/
specifically
- bramila_clean_signal.m (to further clean the fMRI data as described in Power et al. 2014 http://www.sciencedirect.com/science/article/pii/S1053811913009117)
- bramila_diagnostics.m (for quality control)

## Graph-theoretical analysis (macro-level, mesoscopic, micro-level) ##
Tools for graph-theoretical analysis in Python 2.7 are available at:
https://git.becs.aalto.fi/rmkujala/brainnets

The scripts are using the extensive library developed by the Complex Networks group at the Neuroscience and Biomedical Engineering department of Aalto University
https://git.becs.aalto.fi/complex-networks/verkko/tree/master
http://becs.aalto.fi/en/research/complex_networks/

## ABIDE dataset ##
Please refere to the readme of the subfolder ABIDE 

## Statistics ##
The following MATLAB functions and toolboxes were used for permutation based statistics:
* bramila_ttest_np.m from https://git.becs.aalto.fi/bml/bramila/
..* Used to efficiently compute difference of the means using permutations. It uses Matlab parallel computing toolbox.
* bramila_mantel.m from https://git.becs.aalto.fi/bml/bramila/ (a copy available also in the ABIDE subfolder)
..* Used to perform Mantel test (correlation between two distance/similarity matrices). 
* Micro-level statistics were computed with http://bia.korea.ac.kr/people/~cheolhan/software/ 


