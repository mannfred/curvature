[Note: Thanks to K. Thompson, lead author of "Patterns, predictors, and consequences of dominance in hybrids" (2021, The American Naturalist), from which the structure of this Dryad repository is emulating.]

This readme accompanies the paper "Plant-pollinator specialization: Origin and measurement of curvature" (Boehm MMA, Jankowski JE, Cronk QCB, 2021). The following is an ordered, annotated list of .R files briefly describing their contribution to the manuscript. If reproducing the analyses, I recommend opening the .Rproj file in RStudio and using the package `here`; this way, all scripts will run without needing to modify any of the working directories. Scripts are ordered alpha-numerically. Scripts that share the same prefix number are related, but seemed too long to put into one .R file. 

Within the RStudio Project, the "data" folder contains two subfolders: the original raw data used in this study ("raw_data"), and several saved datasets that allow you to reproduce some individual results without having to re-do all previous steps in R ("derived_data"). Data that are not the products of any other scripts are in the "raw_data' subfolder. All variable descriptions are outlined in the metadata.txt file. 

Please feel free to direct any questions to mannfred.boehm@ubc.ca. 



# -------------------


/Rscripts

1_epimedium_garden_data.R
	This script imports "epimedium_growth_data.csv", the raw data collected from our Epimedium development study, and cleans and reorganizes it. It then runs some test to determine if the stages we assigned are statistically distinct, and groups similar stages together. This script exports "epimedium_growth_data_pivot_redefined_stages.rds"


2_seq_alignment_koreanum.R
2_seq_alignment_violaceum.R
	These scripts import "epimedium_growth_data_pivot_redefined_stages.rds" and use the newly defined stage data to align censored developmental sequences. It then runs tests to determine how many days each (consensus) stage lasts. Exports "stages_days_sort_koreanum.rds" and "stages_days_sort_violaceum.rds" for script #4. 


3_curvature_splines_dorsal.R
	This script imports the .tps files created in tpsUtil and tpsDig (landmarked Epimedium specimens). It subsets the landmarks to the dorsal curve, rotates them so the chord is parallel to the xaxis, and fits circles to them. It then exports "circles_fitby_pracma.csv", and "circle_rms_error_pracma.csv". It then computes total curvature and arclength and exports the dataframe "spline_curvature_tbl_dorsal.rds". 


3_geomorph_analysis.R
	This script imports the size-class data "epimedium_curv_size_data_nogran.csv" and .tps files (see above), both which are needed for morphometric trajectory analysis. Exports "proc_data.rds" and "pca_data.rds". 


4_plotting_stage_vs_curv.R
	Imports the total curvature estimates ("spline_curvature_tbl_dorsal.rds") and size-class data (see above) to visualize and test how curvature varies with developmental stage.


4_plotting_stage_vs_time_and_size.R
	Imports "stages_days_sort_koreanum.rds" and "stages_days_sort_violaceum.rds" from script #2 and "epimedium_growth_data_pivot_redefined_stages.rds" from script #1 to visualize and test differences in duration of stages between the two taxa. 


5_ordinations
	Imports the the three .rds files from the two "#3 scripts" to test for correlations between the PCs of shape space and total curvature. Some exploratory analyses (e.g. RDA) are included here, but were not a part of the manuscript.


6_measure_altmetrics
 	Imports the .tps files and makes linear measurements from the landmarks. Exports "alternative_metrics.csv".


7_compare_altmetrics
	Imports the metrics from script #6, the total curvature estimates from script #3, and the fitted circles and rms error estimates (also from script #3). Runs pairwise regression for all curvature metrics and extracts the residuals from the model regressing the inverse radius and total curvature metrics. It then fits another linear model to rms error and the residuals from the previous model. 







	
