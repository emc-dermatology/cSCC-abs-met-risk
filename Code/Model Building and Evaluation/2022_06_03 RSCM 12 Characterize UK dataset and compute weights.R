#---------------------------------------------------
# Aim: Preprocess NCC UK dataset, and compute the weights.
# Author: B. Rentroia Pacheco
# Input: UK dataset
# Output: Analysis and weight computation of the UK dataset
#---------------------------------------------------

# 0. Setup
#---------------------------------------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_02_17 RSCM_8 A1 Internal Validation MI_Boot.R"))
script_id = "RSCM12"

# Decide the sets for the numerator and denominator of the incidence: 
# 2013_onw: only consider patients who do not have any cSCC prior 2013
# 2013_diag: only consider patients who had a cSCC diagnosed in 2013
# all: all available data
weight_scheme = "2013_diag" 
#---------------------------------------------------

# 1. Load datasets
#---------------------------------------------------
#NCC dataset:
DF_UK = read.spss(paste0(dir_project,"Intermediate Data/UK cohort/UK_data_multimp_18052021_2.sav"), to.data.frame=TRUE, use.value.labels=T, convert.dates=T)
#Full cohort info:
uk_fullcohort = read.csv(paste0(dir_project,"Intermediate Data/UK cohort/RSCM 12_1 UK cohort merged.csv"))

#---------------------------------------------------
# 2. Preprocess the UK dataset
#---------------------------------------------------
DF_UK=DF_UK %>%
  mutate(Sex = relevel(Sex, ref="Female"),
         site_primary_cSCC_categories =  relevel(site_primary_cSCC_categories, ref="Trunk + extremities"),
         Ethnicity = relevel(DF_UK$Ethnicity, ref=2),
         mets_numeric = mets,
         tissue_involvement=relevel(tissue_involvement, ref="Dermis"),
         mets = factor(ifelse(mets==1,"Case","Control"),levels=c("Control","Case")),
         setID = factor(setID),
         transplant_new = factor(ifelse(transplant_new==1,"OTR","No OTR"),levels=c("No OTR","OTR")),
         haem_new = factor(ifelse(haem_new==1,"HM","No HM"),levels=c("No HM","HM")))%>%
  as.data.frame()

# Note: After talking with Zoe Venables, it seems that the fu_metastasis_years computed in the previous dataset needs to be updated, so we use the one Zoe provided.
ncc_ids = unique(DF_UK$pid)
uk_full_cohort_ncc_subset = uk_fullcohort[uk_fullcohort$anonpid%in%ncc_ids,]
merged_uk_subs = merge(uk_full_cohort_ncc_subset,DF_UK[,c("pid","Sex","age","fu_metastasis_yrs","fu_yrs","mets")],by.x="anonpid",by.y="pid")
merged_uk_subs$metsfreesurv =merged_uk_subs$metsfreesurv/365 
merged_uk_subs=merged_uk_subs %>% 
  mutate(age_x_minus_y=age.x-age.y,
         is_sex_equal=(ifelse(sex==1,"Male","Female")==Sex),
         also_mets = primary1315mets ==mets.y,
         fup_diff = metsfreesurv-fu_metastasis_yrs) %>%
  mutate(anonpid = factor(anonpid)) %>% group_by(anonpid) %>%
  slice_min(order_by = age_x_minus_y)

#Note: previous metastasis follow up is the same as the one sent in the full cohort for cases:
hist(merged_uk_subs$fup_diff[merged_uk_subs$mets.y=="Case"]) 
# For controls there are more differences:
hist(merged_uk_subs$fup_diff[merged_uk_subs$mets.y=="Control"]) 

# We replace the previous follow up with the one provided by Zoe Venables in the most recent version of the dataset:
DF_UK$fu_metastasis_yrs = merged_uk_subs$metsfreesurv[match(DF_UK$pid,merged_uk_subs$anonpid)]

# Add info on whether the cSCC belongs to a patient with records before 2013:
DF_UK$recs_before_2013 = merged_uk_subs$PRE13SCC[match(DF_UK$pid,merged_uk_subs$anonpid)]
DF_UK$cscc_diagnosed_2013 = merged_uk_subs$YEAR13[match(DF_UK$pid,merged_uk_subs$anonpid)]

# We remove all follow up columns from the dataset, because they are not double checked.
DF_UK[, c("fu_recurrence_yrs","fu_end_histology_yrs","fu_yrs")]=NULL

# Preprocessing uk full cohort:
uk_fullcohort=uk_fullcohort %>% 
  mutate(mets = ifelse(is.na(mets),0,mets),
         VS = ifelse(VS=="X",NA,VS),
         metsfreesurv = ifelse(metsfreesurv<0,0,metsfreesurv),
         primary1315mets =  ifelse(is.na(primary1315mets),0,primary1315mets)) %>%as.data.frame()

# If we are computing incidence on a subset of the cohort, we need to filter out some patients:
if(weight_scheme =="2013_onw"){
  uk_fullcohort = uk_fullcohort %>% filter(PRE13SCC==0)%>%as.data.frame()  #Only patients without a cSCC before 2013 will be considered, since data are incomplete before this year
}else if (weight_scheme =="2013_diag"){
  uk_fullcohort = uk_fullcohort %>% filter(YEAR13==1)%>%as.data.frame() # Only patients who had a cSCC diagnosed in 2013
}


#---------------------------------------------------
# 3. Patient characteristics
#---------------------------------------------------
# Create function to generate a gtsummary, as this is performed multiple times in the code:
generate_pt_table = function(dt, filename){
  pat_charc_UK= dt%>%
    sjlabelled::remove_all_labels() %>%
    tibble %>%
    tbl_summary(by = mets, type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"),include = -setID) %>%
    add_n() %>%
    add_p(list(all_continuous() ~ "paired.wilcox.test",all_categorical()~"mcnemar.test"),group=setID) %>% # test for a difference between groups
    # add column with total number of non-missing observations
    modify_header(label = "**Variable**") %>% # update the column header
    bold_labels()
  
  pat_charc_UK %>%
    as_flex_table() %>%
    flextable::save_as_docx(path=filename)
  
}

# Register patient characteristics before imputation
vars_pt_characteristics = c("mets","Sex","fu_metastasis_yrs","age","previous_cSCC_n_continuous","site_primary_cSCC_categories","diameter_num","thickness_num","tissue_involvement","Differentiation_cat","site_first_metastasis","Morphology","perineural_invasion","lymphovascular_invasion","immunonew","transplant_new","haem_new")

#Save patient characteristics table (all NCC patients)
filename_all = paste0(dir_int_results,"RSCM 12 Patient_characteristics stratified UK Model coding_all.docx")
generate_pt_table(DF_UK[,c("setID",vars_pt_characteristics)],filename_all)


## EXCLUDE PATIENTS ##
# Remove patients with metastasis at baseline:
cases_with_metastasis_at_baseline = intersect(which(DF_UK$mets=="Case"),which(DF_UK$fu_metastasis_yrs==0))
set_ids_mets_baseline = as.character(DF_UK$setID[cases_with_metastasis_at_baseline])
pids_mets_baseline =as.character(DF_UK$pid[cases_with_metastasis_at_baseline])

# Exclude controls with data problems:
exc_patients= c(636233,768354,955768)# These patients must be excluded because there is a data quality concern about the follow up (636233,768354).

DF_UK_excluded = DF_UK %>%filter(as.character(setID)%in%set_ids_mets_baseline|pid%in%exc_patients)%>%as.data.frame()
set_ids_exclude = as.character(DF_UK$setID[cases_with_metastasis_at_baseline])

DF_UK = DF_UK %>%filter(!as.character(setID)%in%set_ids_mets_baseline & (!pid%in%exc_patients))%>%as.data.frame()

filename_nobsmets = paste0(dir_int_results,"RSCM 12 Patient_characteristics stratified UK Model coding_all_no mets baseline.docx")
generate_pt_table(DF_UK[,c("setID",vars_pt_characteristics)],filename_nobsmets)

# Removing patients in the NCC, if incidence is meant to be computed in a subset of the original full cohort:
if(weight_scheme %in% c("2013_onw","2013_diag")){
  if(weight_scheme =="2013_onw"){
    #Remove patients with csccs before 2013
    DF_UK = DF_UK %>%filter(recs_before_2013==0)%>%as.data.frame()
    DF_UK_excluded = DF_UK_excluded%>%filter(recs_before_2013==0)%>%as.data.frame()
  }else if (weight_scheme =="2013_diag"){
    #Remove patients whose cSCC was not diagnosed in 2013
    DF_UK = DF_UK %>%filter(cscc_diagnosed_2013==1)%>%as.data.frame()
    DF_UK_excluded = DF_UK_excluded%>%filter(cscc_diagnosed_2013==1)%>%as.data.frame()
  }
  pat_characteristics_subset_name = paste0(dir_results,"RSCM 12 Patient_characteristics stratified UK Model coding",weight_scheme,".docx")
  
  generate_pt_table(DF_UK[,c("setID",vars_pt_characteristics)],pat_characteristics_subset_name)
}

# Due to the exclusion of cases with metastasis at baseline + possibly the exclusion of cases/controls due to the weight_scheme, we have to rearrange the case-control pairs, to maximize the number of pairs in the final dataset:
# The excluded patients can be replaced by a control pair of a patient with a metastasis at baseline:
set_ids_replace_control = setdiff(as.character(DF_UK_excluded$setID[DF_UK_excluded$pid%in%exc_patients]),set_ids_exclude)
# Controls with a cscc before 2013 can also be replaced by a control pair of a patient with a metastasis at baseline:
cases_without_a_pair = DF_UK %>% filter(setID%in%(names(table(DF_UK$setID))[table(DF_UK$setID)==1])) %>%
  filter(mets=="Case")%>%arrange(fu_metastasis_yrs)%>%pull(setID)%>%as.character()
set_ids_replace_control = unique(c(set_ids_replace_control,cases_without_a_pair))

# Gather all excluded controls together:
# Controls from cases with metastasis at baseline:
DF_UK_excluded_controls = DF_UK_excluded %>%filter(mets=="Control") %>%as.data.frame()
# Add controls without a pair to the dataframe above:
controls_without_a_pair = DF_UK %>% filter(setID%in%(names(table(DF_UK$setID))[table(DF_UK$setID)==1])) %>%
  filter(mets=="Control")%>%as.data.frame()
DF_UK_excluded_controls = rbind(DF_UK_excluded_controls,controls_without_a_pair)
# and remove them from the original dataframe
DF_UK =DF_UK %>% filter(!setID %in%controls_without_a_pair$setID)%>%as.data.frame()

# Rescue some controls back into the NCC dataset, by matching again on follow up. 
for (set_id in set_ids_replace_control){
  case_fup = DF_UK[which(DF_UK$setID==set_id),"fu_metastasis_yrs"]
  control_could_be_rescued = DF_UK_excluded_controls[which(DF_UK_excluded_controls$fu_metastasis_yrs>case_fup),]
  if(nrow(control_could_be_rescued)>0){
    rescued_control = control_could_be_rescued[which.min(control_could_be_rescued$fu_metastasis_yrs),]
    rescued_control[,"setID"] = set_id
    DF_UK = rbind(DF_UK,rescued_control)
    DF_UK_excluded_controls = DF_UK_excluded_controls %>%filter(!pid %in%rescued_control[["pid"]]) %>%as.data.frame()
  }
}

#Remove cases for which no pair could be found:
cases_without_a_pair = DF_UK %>% filter(setID%in%(names(table(DF_UK$setID))[table(DF_UK$setID)==1])) %>%
  filter(mets=="Case")%>%arrange(fu_metastasis_yrs)%>%pull(setID)%>%as.character()

DF_UK = DF_UK[!as.character(DF_UK$setID)%in%cases_without_a_pair,]

## Save patient characteristics and final dataframe ##
filename_after_rearrangement = paste0(dir_results,"RSCM 12 Patient_characteristics stratified UK",weight_scheme,"_final.docx")
generate_pt_table(DF_UK[,c("setID",vars_pt_characteristics)],filename_after_rearrangement)

DF_UK = droplevels(DF_UK)

save(DF_UK, file=paste0(dir_int_output,script_id,"UK_dataset",weight_scheme,".RData"))  

# Check median FUP time in the NCC:
quantile(prodlim(Hist(fu_metastasis_yrs, mets_numeric)~1, data=DF_UK, reverse=TRUE))

## Characteristics of the full cohort ##

# Number of patients: 
length(unique(uk_fullcohort$anonpid))
#Follow up time (first cSCC):
if (weight_scheme=="2013_diag"){
  # We do not have access to all cSCCs before 2013, just that the patients had a cSCC before 2013. Therefore we assume they only had one cSCC before 2013:
  pts_with_atleast_2 = sum(uk_fullcohort$PRE13SCC)
  pts_with_only_1 = nrow(uk_fullcohort) - pts_with_atleast_2
  #Summary number cSCCs:
  print(summary(as.numeric(c(rep(2,pts_with_atleast_2),rep(1,pts_with_only_1)))))
  
  print(quantile(prodlim(Hist(metsfreesurv, mets)~1, data=uk_fullcohort, reverse=TRUE)))
}
#Obtain Sex, Age and Follow up time for the first cSCC ocurrence:
uk_full_cohort_ordered_chrono = uk_fullcohort[order(uk_fullcohort$metsfreesurv),]
uk_full_cohort_ordered_chrono= uk_full_cohort_ordered_chrono[!duplicated(uk_full_cohort_ordered_chrono$anonpid),]

pat_charc_full_UK=uk_full_cohort_ordered_chrono[,c("sex","age","metsfreesurv")] %>%
  sjlabelled::remove_all_labels() %>%
  tibble %>%
  tbl_summary() %>%
  add_n() %>%
  # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()

pat_charc_full_UK %>%
  as_flex_table() %>%
  flextable::save_as_docx(path=paste0(dir_results,"RSCM 12 Patient_characteristics full UK cohort_",weight_scheme,".docx"))

#Number of cSCCs per patient: each row in the uk_fullcohort dataset is corresponds to a cSCC record
summary(as.numeric(table(as.character(uk_fullcohort$anonpid))))

#---------------------------------------------------
# 4. Compute weights
#---------------------------------------------------
# KM-like weights will be computed for the UK cohort, in a similar manner as they were computed for the Dutch dataset.
# One key difference: matching was only performed based on follow up time.
# Moreover, here we know that some metastatic cases (n=180 out of 1076) were not included in the NCC dataset due to missing variables. However they should be included in the weight computation.

# Note: depending on the variable weight_scheme, the weights are different:
# 2013_onw: only patients without cSCCs prior to 2013 are considered. This is because registry is quite incomplete before 2013.
# 2013_diag: only patients with a cSCC diagnosed in 2013 are considered. This is better because there was a bias on the controls chosen (all of them had a cSCC diagnosed in 2013)
# all: no filtering applied: all patients that developed metastasis/all patients with a cscc between 2013 and 2015

weights_UK = data.frame(anonpid = DF_UK$pid, mets = DF_UK$mets,fup = NA,sampling_prob = 1,weights = NA)

# Weights for the cases:
# Case weights are computed as in Zhou et al. 2022: they correspond to the inverse of the proportion of cases that are sampled for the NCC.
met_incidence_cases = setdiff(as.character(uk_fullcohort$anonpid[which(uk_fullcohort$mets==1)]),as.character(uk_fullcohort$anonpid[which(uk_fullcohort$metsfreesurv<=0)])) #Note: there is one patient (456412), who had a mestastatis at baseline but was not indicated as having a primary+metastasis within 2013. Either way, it should not be counted as an incident case

pi = length(which(DF_UK$mets=="Case"))/length(met_incidence_cases)
weights_UK=weights_UK %>%
  mutate(sampling_prob = replace(sampling_prob,mets=="Case",pi),
          weights =replace(weights,mets=="Case",1/pi))%>%
  as.data.frame()

# The expression for control weights is slightly more complicated.
# It corresponds to the inverse of the sampling probability. 
met_cases_ncc = DF_UK$pid[which(DF_UK$mets=="Case")]
for (m.case in met_cases_ncc){
  m.case_info = uk_fullcohort[which(uk_fullcohort$anonpid==m.case),]
  m.case_fup =  as.numeric(m.case_info["metsfreesurv"])
  control_could_be_sampled_for_mcase = intersect(which(as.numeric(uk_fullcohort$metsfreesurv) >m.case_fup),which(uk_fullcohort$mets!=1))
  pot_controls_df = uk_fullcohort[control_could_be_sampled_for_mcase,]
  samp_probs = pot_controls_df %>% group_by (anonpid) %>% dplyr::mutate(count = dplyr::n()) %>% ungroup() %>% dplyr::mutate(prop = count/dplyr::n())%>%select(anonpid,prop)
  samp_probs_subset = samp_probs %>%filter(anonpid%in%weights_UK$anonpid) %>%as.data.frame()
  #Here we multiply the complement of the sampling probability: (cumulative probability of not being sampled until m.case)*(1-prob_sampled for case m.case)
  weights_UK$sampling_prob[match(samp_probs_subset$anonpid,weights_UK$anonpid)]=weights_UK$sampling_prob[match(samp_probs_subset$anonpid,weights_UK$anonpid)]*(1-samp_probs_subset$prop)
}

weights_UK$sampling_prob[which(weights_UK$mets=="Control")]=1-weights_UK$sampling_prob[which(weights_UK$mets=="Control")]
weights_UK$weights = 1/weights_UK$sampling_prob
weights_UK$fup = uk_fullcohort$metsfreesurv[match(weights_UK$anonpid,uk_fullcohort$anonpid)]

nr_mets_baseline = length(intersect(which(uk_fullcohort$mets==1),which(uk_fullcohort$metsfreesurv==0)))
total_potential_mets = length(unique(uk_fullcohort$anonpid))-nr_mets_baseline
weights_UK$weights_rescaled = weights_UK$weights*total_potential_mets/sum(weights_UK$weights) # Rescaling so that the sum of the weights is equal to the total number of patients that can develop a metastasis.

# Sanity check: weights should increase with decreasing follow up time
plot(weights_UK$weights_rescaled[which(weights_UK$mets=="Control")],weights_UK$fup[which(weights_UK$mets=="Control")],xlab="Rescaled weights",ylab="Follow up time (days)")
write.csv(weights_UK,paste0(dir_int_output,"/",script_id," weights_",weight_scheme," UK cohort.csv"))
#---------------------------------------------------

