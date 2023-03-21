#---------------------------------------------------
# Aim: Check the dataset from the UK, of the full cohort, so that weights can be computed.
# Author: B. Rentroia Pacheco
# Input: UK dataset provided by Zoe Venables on the 12th of November 2022, 9 Jan and 6 Feb
# Output: Sanity checks run on the dataset
#---------------------------------------------------

#---------------------------------------------------
# 0. Setup
#---------------------------------------------------
library(openxlsx)
library(dplyr)
script_id = "RSCM 12_1"
#---------------------------------------------------

#---------------------------------------------------
# 1. Load the data
#---------------------------------------------------
# Load information for the computation of the weights, based on the full cohort:
uk_fullcohort = read.xlsx(paste0(dir_project,"Intermediate Data/UK cohort/2022_09_14 UK full cohort.xlsx")) #This was the dataset that was sent by Zoe Venables, with both cases and potential controls
uk_info2013 = read.csv(paste0(dir_project,"Raw Data/UK PHE/ODR1819_225 (DARS- NIC- 656837) extra data B5.csv")) #This was the dataset sent by Zoe Venables with additional information regarding the controls (whether cSCC was diagnosed in 2013, whether patient had a cSCC before 2013)
uk_cases_info2013 = read.csv(paste0(dir_project,"Raw Data/UK PHE/FILE0170240_NIC656837_NDRS_CanReg_Latest.csv")) #This was the dataset sent by Zoe Venables with additional information regarding the cases (whether cSCC was diagnosed in 2013, whether patient had a cSCC before 2013)
# Load the Nested case-control dataset from the UK
DF_UK = read.spss(paste0(dir_project,"Intermediate Data/UK_data_multimp_18052021_2.sav"), to.data.frame=TRUE, use.value.labels=T, convert.dates=T) #This is the dataset with the cases and controls used in Venables et al 2022.
#---------------------------------------------------

#---------------------------------------------------
# 2. Preprocess the datasets
#---------------------------------------------------

# Merge the uk_fullcohort and the uk_info 2013 datasets. The latter contains 2 additional variables that were requested later. These variables indicate whether the patient had a cSCC before 2013 (PRE13SCC) and whether the tumor was diagnosed in 2013 (YEAR13). This allows the computation of the incidence in different scenarios.
# Note that the patient ids are not exactly the same due to some rounding issues. 
all_controls_info = uk_fullcohort[!uk_fullcohort$mets%in%1,]
uk_fullcohort_merged = cbind(all_controls_info[order(all_controls_info$anonpid),],uk_info2013[order(uk_info2013$ANONPID),])
stopifnot(round(uk_fullcohort_merged$anonpid,9) - uk_fullcohort_merged$ANONPID<0.0000001) #To check that the patient ids match (except for rounding errors)
uk_fullcohort_merged$ANONPID = NULL
uk_fullcohort_merged$mets = 0

mets_info = uk_fullcohort[which(uk_fullcohort$mets==1),]
uk_mets_merged = cbind(mets_info[order(mets_info$anonpid),],uk_cases_info2013[order(uk_cases_info2013$ANONID),])
stopifnot(round(uk_mets_merged$anonpid,9) - uk_mets_merged$ANONID<0.0000001) #To check that the patient ids match (except for rounding errors)
uk_mets_merged$ANONID = NULL
colnames(uk_mets_merged)[which(colnames(uk_mets_merged)=="PRE13")]="PRE13SCC"
uk_fullcohort = rbind(uk_mets_merged,uk_fullcohort_merged)

# Check that the number of metastasis is consistent with Venables 2022 publication:
table(uk_fullcohort$mets) # 1566 is the number of metastasis identified
table(uk_fullcohort$primary1315mets) # 1076 is the number of metastasis identified, with primary in 2013

# Get number of patients with first primary cSCC after 2013:
table(uk_fullcohort$PRE13SCC[!duplicated(uk_fullcohort$anonpid)],uk_fullcohort$mets[!duplicated(uk_fullcohort$anonpid)])

# Check that patients can be matched to the NCC dataset we had access previously: 
ncc_ids = unique(DF_UK$pid)
uk_full_cohort_ncc_subset = uk_fullcohort[uk_fullcohort$anonpid%in%ncc_ids,]
merged_uk_subs = merge(uk_full_cohort_ncc_subset,DF_UK[,c("pid","Sex","age","fu_metastasis_yrs","fu_yrs","mets")],by.x="anonpid",by.y="pid")

# Some preprocessing to check any inconsistencies on follow up:
merged_uk_subs$metsfreesurv =merged_uk_subs$metsfreesurv/365 
merged_uk_subs$mets_or_death_free_surv = merged_uk_subs$fu_metastasis_yrs
merged_uk_subs$mets_or_death_free_surv[which(merged_uk_subs$VS  =="D")] = merged_uk_subs$fu_yrs[which(merged_uk_subs$VS  =="D")]
merged_uk_subs$mets_or_death_free_surv[which(merged_uk_subs$mets.y==1)] = merged_uk_subs$fu_metastasis_yrs[which(merged_uk_subs$mets.y==1)]

merged_uk_subs=merged_uk_subs %>% 
  mutate(age_x_minus_y=age.x-age.y,
         fup_x_minus_y=metsfreesurv-mets_or_death_free_surv ,
         is_sex_equal=(ifelse(sex==1,"Male","Female")==Sex),
         also_mets = primary1315mets ==mets.y) %>%
  mutate(anonpid = factor(anonpid)) %>% group_by(anonpid) %>%
  slice_min(order_by = age_x_minus_y)

# Check how many patients have a mismatch in age
length(which(abs(merged_uk_subs$age_x_minus_y)>0.5)) # for 191 patients, age differs from the previous dataset by more than 0.5 years. This depends on the primary cSCC that was selected. Full cohort has age at first diagnosis. Hence the discrepancies
# Check how many patients have a mismatch in survival time
length(which(abs(merged_uk_subs$fup_x_minus_y)>0.2)) # Only 5 patients -> These are due to data updates
#---------------------------------------------------

# Check a few numbers, to create the flowchart:
#---------------------------------------------------
pts_w_csccs_bef2013 = unique(uk_fullcohort$anonpid[which(uk_fullcohort$PRE13SCC==1)])
pts_w_csccs_onw2013 = unique(uk_fullcohort$anonpid[which(uk_fullcohort$PRE13SCC==0)])
pts_w_diag2013 = unique(uk_fullcohort$anonpid[which(uk_fullcohort$YEAR13==1)])
length(intersect(pts_w_csccs_bef2013,pts_w_diag2013))
length(intersect(pts_w_csccs_onw2013,pts_w_diag2013))
#---------------------------------------------------


# Save merged dataset
#---------------------------------------------------
write.csv(uk_fullcohort,paste0(dir_project,"Intermediate Data/UK cohort/RSCM 12_1 UK cohort merged.csv"))
#---------------------------------------------------
