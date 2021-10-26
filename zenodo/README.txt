This readme accompanies the paper "Plant-pollinator specialization: Origin and measurement of curvature" (Boehm MMA, Jankowski JE, Cronk QCB, 2021). The following is an ordered, annotated list of .R files briefly describing their contribution to the manuscript. If reproducing the analyses, I recommend opening the .Rproj file in RStudio and using the package `here`; this way, all scripts will run without needing to modify any of the working directories. Scripts are ordered alpha-numerically. Scripts that share the same prefix number are related, but seemed too long to put into one .R file. 

Following the R script descriptions is information describing the data used the analysis. Within the RStudio Project, the "data" folder contains two subfolders: the original raw data used in this study ("raw_data"), and several saved datasets that allow you to reproduce some individual results without having to re-do all previous steps in R ("derived_data"). Data that are not the products of any other scripts are in the "raw_data' subfolder. All variable descriptions are outlined (following the R script descriptions) in the metadata.txt file. 

Please feel free to direct any questions to mannfred.boehm@ubc.ca. 



# -------------------

/Rscripts

1_epimedium_garden_data.R
	This script imports "epimedium_growth_data.csv", the raw data collected from our Epimedium development study, and cleans and reorganizes it. It then runs some test to determine if the stages we assigned are statistically distinct, and groups similar stages together. This script exports "epimedium_growth_data_pivot_redefined_stages.rds"


2_seq_alignment_koreanum.R
2_seq_alignment_violaceum.R
	These scripts import "epimedium_growth_data_pivot_redefined_stages.rds" and use the newly defined stage data to align censored developmental sequences. It then runs tests to determine how many days each (consensus) stage lasts. Exports "stages_days_sort_koreanum.rds" and "stages_days_sort_violaceum.rds" for script #4. 


3_circle_fitting.R
	This script imports the .tps files created in tpsUtil and tpsDig (landmarked Epimedium specimens). It subsets the landmarks to the dorsal curve, rotates them so the chord is parallel to the xaxis, and fits circles to them. It then exports "circles_fitby_pracma.csv", and "circle_rms_error_pracma.csv". 


3_curvature_splines_dorsal.R
	This script imports the same .tps files as above, and then computes total curvature and arclength. This file exports the dataframe "spline_curvature_tbl_dorsal.rds". 


3_pointwise_k.R
	Imports the same .tps files as above, and computes point-curvature at each of the 16 landmarks. Also infers the constant curvature rate from a fitted circle and estimates the normalized root mean square distance from point-curvature. 


3_geomorph_analysis.R
	This script imports the size-class data "epimedium_curv_size_data_nogran.csv" and .tps files (see above), both which are needed for morphometric trajectory analysis. Exports "proc_data.rds" and "pca_data.rds". 


4_plotting_stage_vs_curv.R
	Imports the total curvature estimates ("spline_curvature_tbl_dorsal.rds") and size-class data (see above) to visualize and test how curvature varies with developmental stage.


4_plotting_stage_vs_time_and_size.R
	Imports "stages_days_sort_koreanum.rds" and "stages_days_sort_violaceum.rds" from script #2 and "epimedium_growth_data_pivot_redefined_stages.rds" from script #1 to visualize and test differences in duration of stages between the two taxa. 


5_ordinations.R
	Imports the the three .rds files from the two "#3 scripts" to test for correlations between the PCs of shape space and total curvature. Some exploratory analyses (e.g. RDA) are included here, but were not a part of the manuscript.


6_measure_altmetrics.R
 	Imports the .tps files and makes linear measurements from the landmarks. Exports "alternative_metrics.csv".


7_compare_altmetrics.R
	Imports the metrics from script #6, the total curvature estimates from script #3, and the fitted circles and rms error estimates (also from script #3). Runs pairwise regression for all curvature metrics and extracts the residuals from the model regressing the inverse radius and total curvature metrics. It then fits another linear model to rms error and the residuals from the previous model. 


8_allometric_curves.R
	This script runs tests to determine if a set of allometrically scaled curves produce curvature estimates correlated with size.


8_isometric_curves.R
	This script does the same as above, but for a set of isometrically scaled curves. 





# ----------------------------------
The following contains information about variable names in the '/data/raw_data' and '/data/derived_data' directories.


# /data/raw_data 

/epimedium_photos/koreanum/koreanum_appended.tps
/epimedium_photos/koreanum/violaceum_appended.tps
	LM = number of landmarks
	IMAGE = image ID from associated .jpg file
	ID = tpsDig ID
	SCALE = pixels to mm conversion

/epimedium_curv_size_data.csv
/epimedium_curv_size_data_nogran.csv 
	species_individual_panicle_flower = species ID, individual plant ID, inflorescence ID, flower ID
	sepal_size_mm = length of the sepal in millimeters
	stage = developmental stage
	size_class = developmental stage  
	species = species epithet
	indiv = individual ID

/epimedium_growth_data.csv
	Species_Individual_Panicle_Flower = species ID, individual plant ID, inflorescence ID, flower ID
	April_* = date, where * is the day of the month

/epimedium_original_flower_stages.csv
	Original_Stage_Designations = stage labels assigned during the experiment
	New_Stage_Designations = stage labels assigned after grouping statistically indistinguishable groups
	Definition = qualitative description of each stage

/match_matrix.txt
	rows and columns 1-24 = IUPAC amino acid abbreviation
	
# ----------------------------------
# /data/derived_data 
	
/RDS_files/epimedium_growth_data_pivot_redefined_stages.rds
	Species_Individual_Panicle_Flower = species ID, individual plant ID, inflorescence ID, flower ID
	date = date observation made (year-month-day)
	size = sepal size in mm
	stage = developmental stage
	Species_epithet = specific epithet only
	days = age of flower in days
	identity = unique identifier for each flower in study
	flower_ID = flower identifier
	panicle_ID = inflorescence identifier
	indiv_ID = individual identifier
	species_ID = species identifier
	spp_ind_ID = previous two identifiers concatenated
	new_stage = re-binned developmental stages (see variable 'stage')

/RDS_files/geomorph_data_frame.rds
	coords = landmark coordinates in xy plane
	groups = Epimedium species
	colors = hexidecimal RGB triplet for colouring graphs
	new_stage = see previous RDS file
	Csize = centroid size of specimen
	numerical_stage = developmental stage expressed as a number instead of a factor

/RDS_files/pca_data.rds
	pc.summary = summary of each pc axis from the PCA of shape space (see: '3_geomorph_analysis.R' in README)
	pc.shapes = coordinates for end member specimens for each pc of shape space

/RDS_files/proc_data.rds
	coords = xy coordinates for procrustes aligned specimens
	points.VCV = variance-covariance matrix of coordinates between each procrustes aligned specimen
	points.var = variance around each xy landmark coordinate of the consensus configuration
	consensus = xy coordinates of the consensus configuration

/RDS_files/size_data_koreanum.rds
	predicted_days_elapsed = number of days elapsed to reach the current developmental stage
	all other variables defined in the first .rds file

/RDS_files/size_data_violaceum.rds
	all variables defined in the above .rds file

/RDS_files/spline_curvature_tbl_dorsal.rds
	total_K = total curvature (degrees)
	name = species ID, individual plant ID, inflorescence ID, flower ID

/RDS_files/stages_days_sort_koreanum.rds
	ID = species ID, individual plant ID, inflorescence ID, flower ID
	elapsed_days = number of days elapsed to reach the current developmental stage
	stage = developmental stage

/RDS_files/stages_days_sort_violaceum.rds
	all variables defined in the above .rds file

/alternative_metrics.csv
	chord = chord length (mm)
	depth = versine length (mm)
	theta = angle of deflection (degrees)

/circle_rms_error.csv
	rms_error = root mean square error from each circle fit (see: '3_circle_fitting.R' in README)

/circles_fitby_pracma.csv
	radius = radius (mm)
	arcs = arc length (mm)
	csize = centroid size

/epimedium_growth_data_pivot.csv
	all variables defined in first .rds file (this is the original .csv version)	

/msa_koreanum.txt
	column1 = species ID, individual plant ID, inflorescence ID, flower ID
	column2 = multiple sequence alignment

/msa_violaceum.txt
	all variables defined in the above .txt file
	
