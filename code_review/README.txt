steps for the data analysis for the low res fission fusion paper

# TODO: clean code for durations of splits is in merge_analysis_presedente
## DO NOT PUT DATA INTO THIS FOLDER, OTHERWISE IT WILL BE PUBLIC ON GITHUB


#this paper has analysis on three coati groups (Presedente, Galaxy, and Trago). Because Trago didn't show fission-fusion
behaviours, they have less code than Galaxy and Presedente for review. Below is the information for the code that needs
to be reviewed and the names of the data files for each group where the plots were made from.


#install packages:

install.packages(c("dbscan", "rgdal", "lubridate", "stringr", "fields", "viridis", "tidyverse", "hms", "dplyr", "tidyr", "ggthemes", "vioplot"

#files needed for review:

preprocessing code (5 files): 
               Galaxy -> 'preprocess_gps_lowres.R' and 'preprocess_lowres_remove_wrong_data_level1.R'
	       Presedente -> 'preprocess_gps_lowres_presedente.R' and 'preprocess_lowres_remove_wrong_data_level1_presedente.R'
	       Trago -> 'preprocess_Trago.R'

output files for analysis and plotting:
	       Galaxy -> 'galaxy_xy_10min_level1.RData' (a list of 'xs','ys','ts') and 'galaxy_coati_ids.RDS'
	       Presedente -> 'presedente_xy_10min_level1.RData' (a list of 'xs','ys','ts') and 'presedente_coati_ids.RDS'
	       Trago -> 'trago_xy_10min_level0.RData' (a list of 'xs','ys','ts') and 'coati_ids.RData'

functions code:
	     'coati_function_library_V1.R'

plot code:
	     'fission_fusion_galaxy_V1.R'
	     'fission_fusion_presedente_V1.R'
             'fission_fusion_trago_V1.R'

1. Preprocess data:

	a. before R: 
		
		-Galaxy group + Trago group: gps data (on logger.bin file downloaded from eObs collar) are parsed 
		from an eObs decoder to a csv file which contains date, time, latitude, longitude (and some other 
		information which is ignored here)
		
		- For Presedente group: the gps data was downloaded from movebank because there were too many 
		logger.bin files with the data to join them manually in R
	
	b. in R: preprocess_gps_lowres script (for level0) - Some extra data wrangling was needed to put the data 
		for all individuals of each group into the right format for the preprocessing:
		
		- For Galaxy group: the gps data for each individual were on 2 files (because the memory card got 
		full a few days before the data were all collected, so I binded them in the preprocess_gps_lowres 
		to make one txt file per individual with "name", "lon", "lat", and "datetime" script (and saved this 
		file into a new txt file, so each individual had 1 file)
		
		- For Trago group: all gps data were on the same logger.bin file, so did not need any extra restructuring 
		of the data before running the rest of preprocess_gps_lowres code to put these data in the 
		'trago_xy_10min_level0.RData' and 'coati_ids.RData'
		
		- For Presedente group: the csv file containing all individuals data were put into the same txt file 
		structure as Galaxy group for further preprocessing (each individual had a txt file containing 
		"name", "lon", "lat", and "datetime"

        c. in R: preprocess_gps_lowres script (for level0) - for all groups, subsampled the gps to 1 gps point every 
		10 mins and convert gps positions to utm. 
		
	d. use level0 output files for cleaning data in preprocess_lowres_remove_wrong_data_level1.R to make "level1"
	this makes a list with xs, ys, ts in the file called "galaxy_xs_10min_level1.RData
	
        e. level1 data are used for all analysis

2. plots:

fission_fusion_galaxy_V1.R - for code review I have put most plot code for galaxy group into the same script
Figure 1 - ### Visualisation of fission-fusion movements (in diff script)
Figure 2abc - L54
Figure 3a - L94
Figure 3b - L363
Figure 4a - L217
Figure 5 - L386
Figure 6a - L124
*Figure S2a - L74
*Figure S3a - L165
*Figure S4ab - L205

fission_fusion_trago_V1
*Figure S1 - L80

fission_fusion_presidente_V1
Figure 2def - L72
Figure 3c - L92
Figure 3d - L349
Figure 4b - L214
Figure 6b -L127
*Figure S2b - L53
*Figure S3b - L162
*Figure S4cd -L202

#durations of splits is in merge_analysis_presedente.
#figure 1 visualisation in coatithon_notgithub/Galaxy_fission_fusion/split_visual_gal.R

* asterix are supplementary plots