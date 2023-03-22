#---------------------------------------------------
# Aim: Reproduce results of absolute risk model publication
# Author: B. Rentroia Pacheco
# Input: Collected dataset
# Output: Weighted cox model + performance estimates
#---------------------------------------------------
#---------------------------------------------------
# 0. Setup
#---------------------------------------------------
dir_scripts = "" # Directory of the repository
dir_project = "" # Directory with all the data
dir_raw_data = paste0(dir_project,"Raw Data/")
dir_results = paste0(dir_project,"Analyses/Results_",Sys.Date(),"/")
dir.create(dir_results) 
dir_int_output = paste0(dir_results,"/Intermediate output/")
dir_int_results = paste0(dir_results,"/Intermediate results/")
dir.create(dir_int_output)
dir.create(dir_int_results)
source(paste0(dir_scripts,"/2022_09_21 Required libraries.R")) # Libraries are loaded in this script, although they need to be installed
#---------------------------------------------------

#---------------------------------------------------
# 1. Structure dataset
#---------------------------------------------------
# Note: The raw datasets for this publication are:
# 1. NKR data: obtained from NKR. It contains all patients that were diagnosed with a first cSCC in 2007/08 in the Netherlands
# 2. PALGA data: tumor records of all the patients identified in the NKR data 

#Preprocess dataset for model development
source(paste0(dir_scripts,"/Multiple Imputation/2021_11_30 RSCM3_1 Multiple Imputation MICE Functions.R"))
source(paste0(dir_scripts,"/Data preprocessing/2021_11_30 RSCM1_1 Preprocessed Dataset into Predictor Dataset.R"))
#---------------------------------------------------

#---------------------------------------------------
# 2. Investigation of competing risks
#---------------------------------------------------
source(paste0(dir_scripts,"/Data preprocessing/2022_09_21 RSCM1_2 Competing risks.R"))
#---------------------------------------------------

#---------------------------------------------------
# 3. Model fitting and recalibration (includes bootstrap validation): This script takes about 1h to run
#---------------------------------------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_02_08 RSCM 8 Model building.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_12_29 RSCM8_A2 Check assumptions final model.R"))
#---------------------------------------------------

#---------------------------------------------------
# 4. Evaluate the model (internal evaluation)
#---------------------------------------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_02_28 RSCM 9 Compute risk probabilities form final model.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_04_01 RSCM 11 Decision Curve Analysis.R"))
#---------------------------------------------------

#---------------------------------------------------
# 5. Evaluate the model on the UK dataset (external evaluation)
#---------------------------------------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_09_14 RSCM 12_1 Check full cohort UK.R")) 
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_06_03 RSCM 12 Characterize UK dataset and compute weights.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_09_25 RSCM 13 Validate Dutch model on the UK dataset.R"))
#---------------------------------------------------

#---------------------------------------------------
# 6. Analyze results of the survey
#---------------------------------------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_09_14 RSCM 11_1 Analyze results of the survey for decision cutoffs.R"))
#---------------------------------------------------

