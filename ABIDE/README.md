# Intersubject similarity result for ABIDE reproducibility dataset #

Files in this folder:

1. ABIDE_results.mat - output from abide_mantel.m
2. Phenotypic_V1_0b.csv - Original phenotipic file (Composite Phenotypic File) from the ABIDE project http://fcon_1000.projects.nitrc.org/indi/abide/ 
3. abide_mantel.m - computes mantel test using scaled inclusivity values for the selected ABIDE subjects
4. abide_mantel_v2.m - revised version of abide_mantel.m with added site location regression and further head motion checks
5. bramila_mantel.m - the mantel test function, part of BraMiLa tools: https://git.becs.aalto.fi/bml/bramila/
6. bramila_ttest2_np - a permutation based two sample t-test, part of BraMiLa tools: https://git.becs.aalto.fi/bml/bramila/
7. correlation_stability_analysis_sites_regressed.m - [Assessing the stability of results from ABIDE dataset, as described by Felix D. Schönbrodt, Marco Perugini, At what sample size do correlations stabilize?, Journal of Research in Personality, Volume 47, Issue 5, October 2013, Pages 609-612, ISSN 0092-6566, http://dx.doi.org/10.1016/j.jrp.2013.05.009.] (http://www.sciencedirect.com/science/article/pii/S0092656613000858)
8. cbrewer - External toolbox for Cbrewer colormaps - http://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab
9. fd_mean.mat - Mean framewise displacement for each subject
10. fd_mean_NTAS.mat - Mean framewise displacement for each subject, included the NTs from ABIDE
11. headmotion_vs_modularity.m - compares individual values of modularity with mean framewise displacement
12. pheno_details.m - Script to parse Phenotypic_V1_0b.csv and automatically extract the subjects 
13.  scores.mat - Output from pheno_details.m, used by abide_mantel.m
14. si_median.mat - Pairwise similarity matrices based on median scaled inclusivity for the 12 NT consensus subnetworks. Scaled inclusivity was computed from adjacency matrices using https://git.becs.aalto.fi/rmkujala/brainnets/
15. si_median_NTAS.mat - Same as point 14, with added NT controls from ABIDE
 
