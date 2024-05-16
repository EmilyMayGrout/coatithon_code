steps for the data analysis for the low res coati subgrouping paper

# this paper has analysis on three coati groups (Presedente, Galaxy, and Trago). Because Trago didn't show fission-fusion
behaviours, they have less code than Galaxy and Presedente.


# install packages:

install.packages(c("dbscan", "rgdal", "lubridate", "stringr", "fields", "viridis", "tidyverse", "hms", "dplyr", "tidyr", "ggthemes", "vioplot", "sf"))

#data for code review are in the coati subgrouping data folder on Google Drive. These data are put in the correct format in the preprocessing stage:

This folder contains a folder for each coati group which contain the gps data, id data, and metadata stored as .RData and .csv files. 
These gps data are filtered to one gps point every 10 minutes and these data are used in the analysis scripts (they are the output of the preprocessing scripts). 
There are 2 gps files per group, one is the Eastings and Northings (xy) and the other is Latitudes and Longitudes (latlon). 
The CoatiTrioMLmatrix.csv is the genetics matrix for all groups data combined (these data are subsetted accordingly for each group in the genetics script)
All analyses use the xy data.


data_dir is the file directory where all the data is stored (the coati id file and the associated groups gps file 
e.g. for Galaxy the files are 'galaxy_xy_10min_level1.RData' (a list of 'xs','ys','ts') and 'galaxy_coati_ids.RData')
code_dir is the file directory where the coati_function_library_V1.R is stored
plot_dir is the file directory where you want to put the plots - I have a different directory for each group
gps_file is the filename of the data for the specific group 
id_file is the filename of the coati ids with their age/sex class

#----------------------------------------------------------------------------------------------------------
preprocessing code (5 files): 

# the preprocessing code puts the raw data into the correct formats for running the analysis (level0)

# for some of the groups, data were removed for specific individuals (e.g. if a collar had fallen off). 
This data removal were done in the preprocess_lowres_remove_wrong_data_level1_[group].R scripts
              
	       Galaxy -> 'preprocess_gps_lowres_galaxy.R' and 'preprocess_lowres_remove_wrong_data_level1_galaxy.R'
	       Presedente -> 'preprocess_gps_lowres_presedente.R' and 'preprocess_lowres_remove_wrong_data_level1_presedente.R'
	       Trago -> 'preprocess_Trago.R'

output files for analysis and plotting:
	       Galaxy -> 'galaxy_xy_10min_level1.RData' (a list of 'xs','ys','ts') and 'galaxy_coati_ids.RData'
	       Presedente -> 'presedente_xy_10min_level1.RData' (a list of 'xs','ys','ts') and 'presedente_coati_ids.RData'
	       Trago -> 'trago_xy_10min_level0.RData' (a list of 'xs','ys','ts') and 'coati_ids.RData'

#-----------------------------------------------------------------------------------------------------------

functions code:
	     'coati_function_library_V1.R'

fission-fusion analysis code:
	     'fission_fusion_galaxy_V1.R'
	     'fission_fusion_presedente_V1.R'
             'fission_fusion_trago_V1.R'

split duration code:
	     'merge_analysis_galaxy_V1.R'
	     'merge_analysis_presedente_V1.R'

genetics code for MRQAP analysis:
	     'genetics_galaxy_V3.R' (most up to date versions of the scripts)
	     'genetics_presedente_V3.R'
	     'genetics_trago_V1.R' 

plot scripts:
	for Figure 1: 'figure1_fissionfusion_plot.R' 
	#for the google maps image in Figure1a, you need to find your own API for this code to work
	for Figure A11: 'figureS4_subgrouping_overtime_combined.R'

#---------------------------------------------------------------------------------------------------------------
Preprocess data:

	a. before R: 
		
		-Galaxy group + Trago group: gps data (on logger.bin file downloaded from eObs collar) are parsed 
		from an eObs decoder to a csv file which contains date, time, latitude, longitude (and some other 
		information which is ignored here)
		
		- For Presedente group: the gps data was downloaded from movebank because there were too many 
		logger.bin files with the data to join them manually in R
	
	b. in R: preprocess_gps_lowres_[group].R script (for level0) - Some extra data wrangling was needed to put the data 
		for all individuals of each group into the right format for the preprocessing:
		
		- For Galaxy group: the gps data for each individual were on 2 files (because the memory card got 
		full a few days before the data were all collected, so I bound them in the preprocess_gps_lowres_galaxy.R
		to make one txt file per individual with "name", "lon", "lat", and "datetime" script (and saved this 
		file into a new txt file, so each individual had 1 file)
		
		- For Trago group: all gps data were on the same logger.bin file, so did not need any extra restructuring 
		of the data before running the rest of preprocess_Trago.R code to put these data in the 
		'trago_xy_10min_level0.RData' and 'coati_ids.RData'
		
		- For Presedente group: the csv file containing all individuals data were put into the same txt file 
		structure as Galaxy group for further preprocessing (each individual had a txt file containing 
		"name", "lon", "lat", and "datetime"

        c. in R: preprocess_gps_lowres_[group].R script (for level0) - for all groups, subsampled the gps to 1 gps point every 
		10 mins and convert gps positions to utm. 
		
	d. use level0 output files for cleaning data in preprocess_lowres_remove_wrong_data_level1_[group].R to make "level1"
	this makes a list with xs, ys, ts in the file called "[group]_xs_10min_level1.RData
	
        e. level1 files are used for Galaxy and Presedente analysis, Trago analysis use level0 because no additional cleaning was needed

#---------------------------------------------------------------------------------------------------------------

1. calculate split durations:
	-before calculating split durations, the fission_fusion_galaxy_V1.R and fission_fusion_presedente_V1.R code
	must be run so the splits_df dataframe is saved and can be loaded into the merge_analysis scripts.

	-To calculate the duration of splits and merges for both groups, the merge_analysis_galaxy_V1.R must be run 
	before the merge_analysis_presedente_V1.R. This is because the calculations for the durations for both groups 
	is in merge_analysis_presedente_V1.R

#---------------------------------------------------------------------------------------------------------------

2. MRQAP analysis:

	-To run the MRQAP analysis, you need to save the "ffnet_reorder" matrices for each of the groups which are made in the:
	'fission_fusion_galaxy_V1.R', 'fission_fusion_presedente_V1.R', and 'fission_fusion_trago_V1.R' scripts for each group respectively.
	-You will need to load the age and sex homophily matrices (which were manually created in excel for all groups separately), 
	and the genetics matrices which is named 'CoatiTrioMLmatrix.csv' 
	-Must make sure the order of ids for each matrices are the same for the mrqap analysis!! - this is organised in the script

#---------------------------------------------------------------------------------------------------------------

4. plotting:

Figure 1 - (a-c) in figure1_fissionfusion_plot.R script
	   (e) in fission_fusion_galaxy_V1.R	
Figure 2 - (a-c) in fission_fusion_galaxy_V1.R
	   (d-f) in fission_fusion_presedente_V1.R
Figure 3 - (a) in fission_fusion_galaxy_V1.R
	   (b) in genetics_galaxy_V3.R
	   (c) in fission_fusion_presedente_V1.R
	   (d) in genetics_presedente_V3.R
Figure A1 - (a) in fission_fusion_galaxy_V1.R
	    (b) in fission_fusion_presedente_V1.R
Figure A2 - (a-b) in fission_fusion_galaxy_V1.R
	    (c-d) in fission_fusion_presedente_V1.R
Figure A3 - (a-b) in fission_fusion_galaxy_V1.R
	    (c-d) in fission_fusion_presedente_V1.R
Figure A4 - (a) in fission_fusion_galaxy_V1.R
	    (b) in fission_fusion_presedente_V1.R
Figure A5 - (a) in fission_fusion_galaxy_V1.R
	    (b) in fission_fusion_presedente_V1.R
Figure A6 - (a) in fission_fusion_galaxy_V1.R
	    (b) in fission_fusion_presedente_V1.R
Figure A7 - (a) in fission_fusion_galaxy_V1.R
	    (b) in fission_fusion_presedente_V1.R
Figure A8 - (a) in fission_fusion_galaxy_V1.R
	    (b) in fission_fusion_presedente_V1.R
Figure A9 - (a) in fission_fusion_galaxy_V1.R
	    (b) in fission_fusion_presedente_V1.R
Figure A10 - (b) in genetics_trago_V1.R
	     (a) in fission_fusion_trago_V1.R
Figure A11 - in figureA11_subgrouping_overtime_combined.R
Figure A12 - (a) in fission_fusion_galaxy_V1.R
	     (b) in fission_fusion_presedente_V1.R
Figure A13 - (a) in genetics_galaxy_V3.R
	     (b) in genetics_presedente_V3.R

MRQAP Tables:
Table 1, A1, and A2 - Galaxy results from genetics_galaxy_V3 and Presidente results from genetics_presedente_V3.R
