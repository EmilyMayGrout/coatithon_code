
- The coati subgrouping data folder contains a folder for each coati group (Galaxy, Trago, and Presedente) that contain the gps data, id data, 
and metadata stored as .RData and .csv files. 
	Files per group:
	-coati_id.csv ==> raw csv file used in the preprocessing code to make the galaxy_coati_ids.RData file which is used for all analysis)
	-[group name]_age_matrix.csv ==> matrix manually made in excel for age holophily for group members (1 = same age catagory, 0 is not the same age catagory)
	-[group name]_sex_matrix.csv ==> matrix manually made in excel for sex holophily for group members (1 = same sex, 0 is not the same sex)
	-[group name]_coati_ids.RData ==> an .RData file which contains the individual name, collar ID, age, sex, age/sex color for plotting
	-[group name]_latlon_10min_level0.RData ==> the .RData file which contain 2 matrices and 1 vector:
		- one matrix contains the latitude (lats) and the other contains longitude (lons). Group members are the rows, the column is each time point
		- the vector (ts) gives the times for the column indexes for the lats and lons
	-[group name]_xy_10min_level0.RData ==> the .RData file which contain 2 matrices and 1 vector:
		- one matrix contains the Eastings (xs) and the other contains Northings (ys). Group members are the rows, the column is each time point
		- the vector (ts) gives the times for the column indexes for the xs and ys

-The CoatiTrioMLmatrix.csv is the genetics matrix for all groups (these data are subsetted accordingly for each group in the genetics script)

*** analysis and plotting scripts with a detailed README.txt file can be found on GitHub:
https://github.com/EmilyMayGrout/coatithon_code/tree/main/code_review  
