# Intersubject similarity result for ABIDE reproducibility dataset #

Files in this folder:

1. ABIDE_results.mat - output from abide_mantel.m
2. Phenotypic_V1_0b.csv - Original phenotipic file from the ABIDE project
3. abide_mantel.m - computes mantel test using scaled inclusivity values for the selected ABIDE subjects
4. bramila_mantel.m - the mantel test function, part of BraMiLa tools: https://git.becs.aalto.fi/bml/bramila/
5. cbrewer - External toolbox for Cbrewer colormaps - http://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab
6. fd_mean.mat - Mean framewise displacement for each subject
7. pheno_details.m - Script to parse Phenotypic_V1_0b.csv and automatically extract the subjects 
8. scores.mat - Output from pheno_details.m, used by abide_mantel.m
9. si_median.mat - Pairwise similarity matrices based on median scaled inclusivity for the 12 NT consensus subnetworks. Scaled inclusivity was computed from adjacency matrices using https://git.becs.aalto.fi/rmkujala/brainnets/
 
