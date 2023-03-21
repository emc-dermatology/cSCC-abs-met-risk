#---------------------------------------------------
# Aim: Analyze results of the survey of the SCOUT consortium
# Author: B. Rentroia Pacheco
# Input: Results of the survey
# Output: Summary plots of the survey
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------
library(ggplot2)
library(plyr)
library(dplyr)
library(gtsummary)
library(gridExtra)
library(reshape2)
script_id = "RSCM 11_1"
dir_folder_results = paste0(dir_results,"/",script_id)
#-----------------------

#-----------------------
# 1. Load the data
#-----------------------
survey_results = read.xlsx(paste0(dir_project,"Intermediate Data/2022_09_14 Survey results/2022_09_29 Online_Interview_Questionnaire.xlsx"))
survey_results =survey_results %>% rename(Profession = "What's.your.profession?", 
                          Work_location = "Where.do.you.work?") %>%as.data.frame()
survey_results=survey_results %>% mutate(Cutoff_surveillance = as.numeric(Cutoff_surveillance),
                          Cutoff_adjuvant_treat = as.numeric(Cutoff_adjuvant_treat),
                          Cutoff_systemic = as.numeric(Cutoff_systemic)) %>%as.data.frame()
#-----------------------


#-----------------------
# 2. Explore each variable individually
#-----------------------
# Work location and Profession
par(mar = c(3, 16, 2, 2))
barplot(sort(table(survey_results$Profession)),las=2,horiz=TRUE)
dev.off()

par(mar = c(3, 10, 2, 2))
barplot(sort(table(survey_results$Work_location)),las=2,horiz=TRUE)
dev.off()
survey_results=survey_results %>% mutate(WorkLocation_cats = revalue(Work_location, c("Oceania"="Other","Asia"="Other","South America"="Other"))) %>%as.data.frame()

survey_results=survey_results %>% mutate(Profession_cats = factor(Profession),WorkLocation_cats = factor(WorkLocation_cats,levels=c("Europe","North America","Other")))%>%as.data.frame()
cutoffs = c("Cutoff_surveillance","Cutoff_adjuvant_treat","Cutoff_systemic")

table_pc <- 
  survey_results[,c("Profession_cats","WorkLocation_cats",cutoffs)] %>%
  sjlabelled::remove_all_labels() %>%
  tibble %>%
  tbl_summary( type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous")) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 
table_pc %>%
  as_flex_table() %>%
  flextable::save_as_docx(path=paste0(dir_folder_results,"Survey results all.docx")) # This was the result that was published in the supplementary information

#We also looked at dermatologists only, as these are the ones that provide routine follow up to cscc patients:
table_derm <- 
  survey_results[survey_results$Profession_cats%in%c("Dermatologist","Dermatologist, Mohs surgeon "),c("Profession_cats","WorkLocation_cats",cutoffs)] %>%
  sjlabelled::remove_all_labels() %>%
  tibble %>%
  tbl_summary( type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous")) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() 
table_derm %>%
  as_flex_table() %>%
  flextable::save_as_docx(path=paste0(dir_folder_results,"Survey results DERM.docx"))

#-----------------------

#-----------------------
# 3. Explore cutoffs per profession and geographical area
#-----------------------
survey_results 
xlabs_cutoffs =  c("Metastatic risk threshold\n for surveillance","Metastatic risk threshold\n for adjuvant therapy","Metastatic risk threshold\n for systemic treatment")
plots_cutoffs= vector(mode = "list", length = 9)
for(i in c(1:length(cutoffs))){
  cutoff_i = cutoffs[i]
  plots_cutoffs[[i]]=ggplot(survey_results,aes_string(x=cutoff_i))+geom_histogram()+theme_bw()+xlim(c(0,50))+ylim(c(0,17))+xlab(xlabs_cutoffs[i])
  plots_cutoffs[[i+3]]=ggplot(survey_results,aes_string(y=cutoff_i,x="Profession_cats",fill="Profession_cats"))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_x_discrete(guide = guide_axis(n.dodge = 2))+ylim(c(0,60))+ylab(xlabs_cutoffs[i])
  plots_cutoffs[[i+6]]=ggplot(survey_results,aes_string(y=cutoff_i,x="WorkLocation_cats",fill="WorkLocation_cats"))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_x_discrete(guide = guide_axis(n.dodge = 2))+ylim(c(0,60))+ylab(xlabs_cutoffs[i])
}
g_cat=grid.arrange(grobs=plots_cutoffs,nrow=3,ncol=3) 
ggsave(paste0(dir_folder_results, "distribution answers cutoffs.pdf"),width=25,height=9.1,g_cat)

survey_results.m = melt(survey_results[,cutoffs])
survey_results.m$variable = factor(revalue(survey_results.m$variable,c("Cutoff_surveillance" = "Surveillance cutoff","Cutoff_adjuvant_treat"="Adjuvant treatment cutoff","Cutoff_systemic"= "Systemic treatment cutoff")))
p_cutoff = ggplot(survey_results.m,aes(y=value,x = variable,fill = variable))+geom_boxplot()+theme_bw()+ylab("Cutoff (%)")+xlab("")
ggsave(paste0(dir_folder_results, "distribution answers overall.pdf"),width=9,height=7,p_cutoff)
#-----------------------
