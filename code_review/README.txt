steps for the data analysis for the low res fission fusion paper

# TODO: get all code for all of the first paper into the code_review folder in dropbox/coatithon_code
## DO NOT PUT DATA INTO THIS FOLDER, OTHERWISE IT WILL BE PUBLIC ON GITHUB


#this paper has analysis on three coati groups (Presedente, Galaxy, and Trago). Because Trago didn't show fission-fusion
behaviours, they have less code than Galaxy and Presedente for review. Below is the imformation for the code that needs
to be reviewed and the names of the data files for each group where the plots were made from.

#files needed for review:

preprocessing code: 
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
	     'fission_fusion_V1.R'
	     'fission_fusion_presidente_V1.R'
             'fission_fusion_trago_V1.R'

1. Preprocess data:
Galaxy group
	1. before R: gps data (on logger.bin file) are parsed from a decoder to latitude and longitude
	2. in R: preprocess_gps_lowres (for level0) which is subsampling to 1 gps point every 10 mins 
		and convert gps positions to utm
	3. use level0 for cleaning data in preprocess_lowres_remove_wrong_data_level1.R to make "level1"
	this makes a list with xs, ys, ts in the file called "galaxy_xs_10min_level1.RData
	4. level1 data are used for all analysis

2. plots:

fission_fusion_galaxy_V1.R - for code review I have put most plot code
for galaxy group into the same script
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