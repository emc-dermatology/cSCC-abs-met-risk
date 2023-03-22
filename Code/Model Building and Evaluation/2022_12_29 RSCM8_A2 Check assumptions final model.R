#---------------------------------------------------
# Aim: Check assumptions of final Cox model
# Author: B. Rentroia Pacheco
# Input: Imputed dataset
# Output: Diagnostic plots
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------
folder_name = paste0(dir_int_results,"MICE/v3_0 coding/_ws_","manual","_ovr_","monotone","_nimp_10_cpm_","1e-04")
load(paste0(folder_name,"/RSCM3_1 MICE object.RData")) # multiple imputed dutch dataset (dutch_ds), can be fecthed with dutch_ds_mimp object.
load(paste0(dir_int_output,"RSCM1_1 Dutch Dataset coding_v3_0.RData"))

dir.create(paste0(dir_int_results,"Cox Assumptions/"))
dir_results_assum = paste0(dir_int_results,"Cox Assumptions/RSCM 8A2_")

#-----------------------
# 1. Check proportionality assumption
#-----------------------
test.ph_all=data.frame()
for(i in c(1:10)){
  final_fit_simp = coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric) ~ Sex + Age+Tumor_location_cats + Number_of_cSCC_before_culprit + Tumor_diameter+ Differentiation + PNI_or_LVI +Tissue_involvement,data=complete(dutch_ds_mimp,i),weights = weights)

  # Cox zph
  test.ph = as.data.frame(cox.zph(final_fit_simp)$table) # All okay!
  test.ph$imp = i
      
  test.ph_all = rbind(test.ph_all,test.ph)
  
  
}
test.ph_all = test.ph_all[-1,]
write.csv(test.ph_all,paste0(dir_results_assum,"PH assumption.csv"))

#-----------------------
# 2. Check whether splines are required
#-----------------------
anova_tests = as.data.frame(matrix(NA,ncol=1,nrow=0))
for (i in c(1:10)){
  ds=complete(dutch_ds_mimp,i)
  ds =  ds[rep(c(1:nrow(ds)),round(weights)),]
  rknots=3
  fit_null=coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric) ~Age+ Number_of_cSCC_before_culprit + Tumor_diameter,data=ds)
  fit_rcsage=coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric) ~rcs(Age,rknots)+ Number_of_cSCC_before_culprit + Tumor_diameter,data=ds)
  fit_rcsnrcscc=coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric) ~Age+ rcs(Number_of_cSCC_before_culprit,rknots+1) + Tumor_diameter,data=ds)
  fit_rcstdiam=coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric) ~Age+ Number_of_cSCC_before_culprit + rcs(Tumor_diameter,rknots),data=ds)
  
  anova_fit = anova(fit_null, fit_rcsage) 
  anova_tests [1,] = anova_fit[,4][2] #This gives the p-value of the chi-squared test
  
  anova_fit = anova(fit_null, fit_rcsnrcscc) 
  anova_tests [2,] = anova_fit[,4][2] #This gives the p-value of the chi-squared test
  
  anova_fit = anova(fit_null, fit_rcstdiam)
  anova_tests [3,] = anova_fit[,4][2] #This gives the p-value of the chi-squared test
  
}
write.csv(anova_tests,paste0(dir_results_assum,"splines_investigation.csv"))

