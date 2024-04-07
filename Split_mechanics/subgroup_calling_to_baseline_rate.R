

#load the changed_df made in the calling_during_fission_and_fusion script, called ind_fission_data_long_filt
load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/changed_df.RData")

#load the call_rates_together df from find_callrates_grouptogether
load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/call_rates_together_gal.RData")



subgroup_change <- ind_fission_data_long_filt

subgroup_change_with_ind_baseline<-data.frame()
for(i in unique(subgroup_change$event_idx)){
  sc<-subgroup_change[which(subgroup_change$event_idx == i),]
  # change_ids<-sc[sc$change_group == "change","ind_idx"]  
  # nchange_ids<-sc[sc$change_group == "not_change","ind_idx"]
  # 
  sc$base_calls_mean<-NA
  for(j in sc$ind_idx){
    sc[sc$ind_idx == j & sc$call == "contact_call_rate", "base_calls_mean"]<-mean(call_rates_together[which(call_rates_together$ind_idx == j & call_rates_together$call_type == "contact call" ),"call_rate"])
    sc[sc$ind_idx == j & sc$call == "agg_call_rate", "base_calls_mean"]<-mean(call_rates_together[which(call_rates_together$ind_idx == j & call_rates_together$call_type == "aggression call" ),"call_rate"])
  }
  subgroup_change_with_ind_baseline<-rbind(subgroup_change_with_ind_baseline,sc)
  
}
