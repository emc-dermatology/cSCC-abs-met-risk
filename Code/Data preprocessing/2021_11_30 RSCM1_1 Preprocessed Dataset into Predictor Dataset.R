#---------------------------------------------------
# Aim: Generate datasets for model development
# Author: B. Rentroia Pacheco
# Input: Preprocessed Dataset
# Output: Input dataset for the clinical prediction model
#---------------------------------------------------


# The idea of this script is to add as many variables possible to the preprocessed dataset. The datasets should be ready to be used in the modelling part.

# 0. Setup
#-----------------------
mice_folder = paste0(dir_int_results,"/MICE/")

if(!dir.exists(mice_folder)){
dir.create(mice_folder)
}
#-----------------------

#-----------------------
# 1. Load the data
#-----------------------
load(paste0(dir_int_output,"RSCM1 Dutch Dataset Preprocessed_v3_0.RData"))
#-----------------------


# Dutch Dataset, with variables similar to Tokez 2022
#-----------------------
dutch_ds_baseline_model = dutch_ds

# Remove variables that are not meant to be used in this research/would require some extra work in the collection to be accurate.
rem_vars = c("Ratio_width_diameter","Recurrence_cSCC","Sebaceous_gland_carcinoma_possibility","Morphology_present_less30pct","Resection_margin_cont","Invasion_sum","Tumor_measurements","Number_of_cSCC_after_culprit","Morphology_subgroup","PNI","Lymphovascular_invasion","Verspringende_groei_tumor","Resection_margin_location","Resection_margin_cont_cm","Irradical_location","Extra_remarks_Antien","PNI_diameter","Metastasis_non_histologically_confirmed_date_of_primary")
dutch_ds_baseline_model = dutch_ds_baseline_model[,!colnames(dutch_ds_baseline_model)%in%rem_vars]

# Recode variables the same way that Tokez 2022 did, just for comparison terms with her results.
dutch_ds_baseline_model$Differentiation = factor(revalue(dutch_ds_baseline_model$Differentiation,c("Good (>75% well differentiated)"="Good/moderate","Moderate (25-75% well differentiated)"="Good/moderate","Poor (<25% well differentiated)"="Poor/undifferentiated")),levels=c("Good/moderate","Poor/undifferentiated"))
dutch_ds_baseline_model$Metastasis_numeric=ifelse(dutch_ds_baseline_model$Metastasis=="Case",1,0)
dutch_ds_baseline_model$Sex = factor(as.character(dutch_ds_baseline_model$Sex),levels=c("vrouw","man")) # To ensure reference level is the same
dutch_ds_baseline_model$Tumor_location_cats = factor(revalue(dutch_ds_baseline_model$Tumor_location_cats,c("Face"="Face","Trunk"="Trunk and extremities","Eyelid"="Face","Ear"="Face","Scalp/neck"="Scalp and neck","Lip"="Face","Upper limb"="Trunk and extremities","Lower limb"="Trunk and extremities")),levels = c("Trunk and extremities","Face","Scalp and neck"))

# Survival time should not be used directly in the imputation. Therefore we calculate the Nelson-Aalen estimator:
dutch_ds_baseline_model$H0 = nelsonaalen(data=dutch_ds_baseline_model,timevar = Vitfup_metastasis_years, statusvar = Metastasis_numeric)

#Perform Multiple imputation
# Some variables need to be excluded from the multiple imputation:
# Text variables:
id_vars = c("Key_NKR","Administratie_nr","Key_EID","PALGA_excerpt_id")
text_vars = c("Conclusion","Microscopy")
# Variables that will not be used in the model, because Tumor diameter factorized was a variable created by Selin and Tumor location has too many levels.
exc_vars = c("Tumor_diameter_factorized","Tumor_location",id_vars,text_vars)
# These variables are excluded from the imputation to avoid multicollinearity issues:
exc_vars_imp = c("Metastasis_numeric","Set_id","Tissue_involvement","Vitfup_metastasis_years","Metastasis_location_lymphnode","Metastasis_location_cutaneous","Immunocommpromised_at_time_of_cSCC","PNI_bin","Lymphovascular_invasion_bin")

dir_folder_results_mice = paste0(mice_folder,"Similar to UK dataset")
if(!dir.exists(dir_folder_results_mice)){dir.create(dir_folder_results_mice)}
performMultipleImputation(dutch_ds_baseline_model,exc_vars_ds =exc_vars ,exc_vars_imp = exc_vars_imp ,winsorization="manual",order_vars="monotone",nr_imp_ds=10,cor_pred_mat=0.0001,niter=50,seed=1,dir_folder_results_mice,"")

# Separate variable sets
uk_predictors = c("Sex","Age","Tumor_diameter","Number_of_cSCC_before_culprit","Tumor_location_cats","Breslow_thickness","Differentiation","PNI_or_LVI","Immunocommpromised_at_time_of_cSCC", "Tissue_involvement","Morphology_subgroup_bin")  
outcome_vars = c("Metastasis","Vitfup_metastasis_years")
obj_save = dutch_ds_baseline_model[,c(outcome_vars,uk_predictors)]
local({
  dutch_ds <- dutch_ds_baseline_model[,c(outcome_vars,uk_predictors)]
  save(dutch_ds, file=paste0(dir_int_output,"RSCM1_1 Dutch Dataset UKcoding_v3_0.RData"))  
})
#-----------------------

# Dutch dataset,with variable coding for model development
#-----------------------
dutch_ds_v1 = dutch_ds_baseline_model
dutch_ds_v1$isBiopsy = revalue(dutch_ds$Biopt_of_reex_of_shaveex,c("Re-excisie"="Excisie","Shave excisie"="Excisie","Normale excisie"="Excisie"))
dutch_ds_v1$Biopt_of_reex_of_shaveex = NULL
dutch_ds_v1$Tumor_location_cats = factor(revalue(dutch_ds$Tumor_location,c("Forehead"="Face","Eyelid"="Face","Coeur"="Trunk/Extremities","Arm"="Trunk/Extremities","Hand"="Trunk/Extremities","Leg"="Trunk/Extremities","Feet"="Trunk/Extremities","Buttocks"="Trunk/Extremities","Nose"="Face","Cheeks"="Face","Jaw"="Face","Lip"="Face","Vertex"="Scalp/neck","Neck"="Scalp/neck","Retroauricular"="Face","Ear"="Face","Trunk"="Trunk/Extremities")),
                                         levels = c("Trunk/Extremities","Scalp/neck","Face"))

dutch_ds_v1$Keratinization = factor(revalue(dutch_ds_baseline_model$Keratinization,c("Good (>30% keratinized)"="Good","Moderate (10-30% keratinized)"="Poor/moderate","Poor&undifferentiated (<10% keratinized)"="Poor/moderate")),levels=c("Good","Poor/moderate"))
dutch_ds_v1$Peritumoral_infiltration = factor(revalue(dutch_ds_baseline_model$Peritumoral_infiltration,c("Absent"="Absent/Moderate","Moderate"="Absent/Moderate","Abundant"="Abundant")),levels=c("Absent/Moderate","Abundant"))
dutch_ds_v1$Solar_elastosis = factor(revalue(dutch_ds_baseline_model$Solar_elastosis,c("No/nihil"="None/Moderate","Moderate"="None/Moderate","Extensive"="Extensive")),levels=c("None/Moderate","Extensive"))

dutch_ds_v1$logBreslowThickness = log(dutch_ds_v1$Breslow_thickness)
dutch_ds_v1$Breslow_thickness = NULL

dir_folder_results_mice = paste0(mice_folder,"v3_0 coding")
if(!dir.exists(dir_folder_results_mice)){dir.create(dir_folder_results_mice)}

performMultipleImputation(dutch_ds_v1,exc_vars_ds =exc_vars ,exc_vars_imp = exc_vars_imp ,winsorization="manual",order_vars="monotone",nr_imp_ds=10,cor_pred_mat=0.0001,niter=50,seed=1,dir_folder_results_mice,"")
local({
  dutch_ds <- dutch_ds_v1
  save(dutch_ds, file=paste0(dir_int_output,"RSCM1_1 Dutch Dataset coding_v3_0.RData"))  
})
#-----------------------

# Patient characteristics with new coding:
#-----------------------
table2 <- 
  dutch_ds_v1[,!colnames(dutch_ds_v1)%in%setdiff(c(id_vars,text_vars),"Set_id")] %>%
  sjlabelled::remove_all_labels() %>%
  #mutate_if(is.factor,fct_explicit_na,na_level = "Unknown") %>%
  tibble %>%
  tbl_summary(
    by = Metastasis, type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"),include = -Set_id) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(list(all_continuous() ~ "paired.wilcox.test",all_categorical()~"mcnemar.test"),group=Set_id) %>% # test for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()

table2 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path=paste0(dir_results,"RSCM1_1 Patient_characteristics stratified 390pts Model coding.docx"))
#-----------------------

