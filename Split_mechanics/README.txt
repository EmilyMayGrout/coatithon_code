In the Split_mechanics folder there are the following scripts:
-----------------------------------------------------------------------
- 1. identify_splits_and_merges.R
	INPUT FILES: xs, ys, ts, coati_ids
	This script uses the "sticky" dbscan method to identify splits and merges using the detect_fissions_and_fusions function from the coati_function_library (L876)
	This function makes a list called ff_data_50 where we have:
				1. events_detected: a dataframe with the tidx for each event, and a row for each event with the individuals in group A and B (same format as manual events)
				2. groups_list
				3. together: [n_inds x n_inds x n_times matrix] of whether individuals are connected (1) or not (0) or unknown (NA)
				4. changes: data frame containing all the subgroups membership changes (not just ones identified as fissions or fusions)
	OUTPUT FILE: "gal_events_detected.Rda" and "pres_events_detected.Rda" --> which is the ff_data_50$events_detected for Galaxy and Presidente
	PLOTTING: - look at one event with analyse_ff_event function
		  - also made a sequence of plots for a chosen event but code not working to save the animation of the event

-----------------------------------------------------------------------
- 2. charecterize_splits_and_merges.R
	- INPUT FILES: either manual events dataframe (processed/split_analysis_processed/group_manual_split_merge_clean.csv), named "events"
		     or the automated events dataframe is created with the detect_fissions_and_fusions function, 
			the events_detected dataframe from this output is extracted and named "events"
	- OUTPUT FILE: This script is extracting additional information about the events and saving it in:
			if auto: processed/2022/galaxy/galaxy_auto_ff_events_characterized.RData (also for presendente)
		 	if manual: processed/2022/galaxy/galaxy_manual_ff_events_characterized.RData (also for presendente)
        - additional info includes: distance travelled for group A and B after a split or before a merge, the full groups distance travelled during the event,
	  turn angles for each group, age of group members, number of adults/subadults/juveniles in each subgroup, subgroup sizes (155 obs x 38 variables)

	- this script also contains Ari's code to look at the overlap between the periods (overlaps in before, during, and after periods)

-----------------------------------------------------------------------
- calling_during_fissions_and_fusions.R
	- This script is for looking at how call rates relate to fission-fusion dynamics
	- Call rate for each individual for each event before, during and after the event (also some additional info, like the subgroup id and the distance moved)
	- INPUT FILES: 'galaxy_auto_ff_events_characterized.RData' made in charecterize_splits_and_merges
			'all_data_hms_synched.csv' --> made from coati_synch.R in coatithon_code using the all_data_hms dataframe made in 
						       coaticalls/call_code/calls_around_fission_fusions_galaxy_prelim.R on L96, 
						       using csv files per individual per day/events and binding them into one dataframe called all_data_hms
	- #still working on charecterising the different sorts of fissions and fusions

-----------------------------------------------------------------------
- compare_automated_manual_ffevents.R
	- making txt files for manual and automated events to compare outside of RStudio --> we found they matched up pretty well, with automated finding more events

-----------------------------------------------------------------------
- find_callrates_grouptogether.R
	- finding the baseline call rates when the group is together, to see if calls rates do change significantly during fissions and fusions

-----------------------------------------------------------------------
- split_mechanics_exploration.R
	- My code where there's some copied code from charecterize_splits_and_merges.R
	- Mostly code for making plots for presenting initial results 
	INPUT: manual events in processed/split_analysis_processed/',group,'_manual_split_merge_clean2.csv'
	       automated events in processed/2022/galaxy/galaxy_auto_ff_events_characterized.RData'
	OUTPUT: saved events in processed/split_analysis_processed/galaxy_manual_events_withinfo2.RData' (files saved for auto_events, and for Presidente group)
	PLOTTING: histogram of speeds before splits, distance travelled by the smaller and larger subgroup (for fissions and fusions), plotting age/sex class that moves

-----------------------------------------------------------------------
- split_mechanics_leadership_markdown.Rmd -> the markdown file is saved in coatithon/results/split_mechanics_leadership_markdown.html

------------------------------------------------------------------------
- split_merge_leadership.R
	- look at leadership in split and merge events
	- created 3 different methods (functions in this script) to assess leadership 
		-1. position along trajectory: FUNCTION -> ind_disp_along_group_path
		-2. group crossing time: FUNCTION -> ind_crossing_thresh_times_along_group_path (own_finish_line = F)
		-3. crossing time for each individuals own finish line) -> ind_crossing_thresh_times_along_group_path (own_finish_line = T)
	-computed entropy for leadership for real vs permuted data

	INPUT: #NEED TO ADD AS ITS A BIT COMPLICATED

	OUTPUT: Leader rank position saved in processed/2022/galaxy/galaxy_LeaderRank_position.RData
					      processed/2022/galaxy/galaxy_LeaderRank_crosstime.RData
                                              processed/2022/galaxy/galaxy_LeaderRank_crosstime_ownfinishline.RData
  						(and saved for Presidente)


