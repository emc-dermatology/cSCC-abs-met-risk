#---------------------------------------------------
# Aim: Compute risk probabilities based on the model built in 2022_02_08 RSCM 8
# Author: B. Rentroia Pacheco
# Input: weighted Cox Regression Model
# Output: Risk probabilities
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------


folder_name = paste0(dir_int_results,"MICE/v3_0 coding/_ws_","manual","_ovr_","monotone","_nimp_10_cpm_","1e-04")
load(paste0(folder_name,"/RSCM3_1 MICE object.RData")) # multiple imputed dutch dataset (dutch_ds), can be fecthed with dutch_ds_mimp object.
load(paste0(dir_int_output,"RSCM1_1 Dutch Dataset coding_v3_0.RData"))
source(paste0(dir_scripts,"/Multiple Imputation/2021_11_30 RSCM3_1 Multiple Imputation MICE Functions.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_01_19 RSCM_A 6 Weighted Performance Metrics functions.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_02_17 RSCM_8 A1 Internal Validation MI_Boot.R"))

script_id = "RSCM9"
ds_id = "v3_0"
#-----------------------

#-----------------------
# 1. Load final model
#-----------------------
weighted_cox_table_coefs_recalib=read.csv(paste0(dir_results,"RSCM8_wghcox_model_recalib",ds_id,".csv"))
coefs_final_model = log(as.numeric(gsub("\\s.*","",weighted_cox_table_coefs_recalib[,"Estimate"]))) #Log because we need the coefficients for the linear predictor, not the Hazard ratios!
important_quantities = read.csv(paste0(paste0(dir_results,"RSCM8_wghcox_model_recalib_importantqntsv3_0.csv")))
mean_age = important_quantities$x[which(important_quantities$X=="Mean age")]
mean_ncscc = important_quantities$x[which(important_quantities$X=="Mean number of cSCC")]
mean_tdiam = important_quantities$x[which(important_quantities$X=="Mean tumor diameter")]
mean_lp = important_quantities$x[which(important_quantities$X =="MeanLP_recalib")]
recalibrated_baseline_survival = round(important_quantities$x[which(important_quantities$X =="Recalibrated baseline survival")],3)

#-----------------------

#-----------------------
# 2. Compute risk probabilities
#-----------------------
model_variables = c("Age","Differentiation","Number_of_cSCC_before_culprit","PNI_or_LVI","Sex","Tissue_involvement","Tumor_diameter","Tumor_location_cats")

#Note: probabilities for each patient are going to be averaged:
surv_probs = matrix(NA,ncol=10,nrow=nrow(dutch_ds))
for (m in 1:dutch_ds_mimp$m){
  ds = complete(dutch_ds_mimp,m)
  dutch_ds_subset_model_form = model.matrix.lm(~.,ds[,c(model_variables,"Metastasis")])[,-1]
  dutch_ds_subset_model_form=dutch_ds_subset_model_form %>% as.data.frame()%>% mutate(Age = Age-mean_age,
                                                                                      Number_of_cSCC_before_culprit = Number_of_cSCC_before_culprit-mean_ncscc,
                                                                                      Tumor_diameter = Tumor_diameter - mean_tdiam)
  lps = as.vector(t(coefs_final_model)%*% t(dutch_ds_subset_model_form[,1:10]))
  surv_probs[,m] = recalibrated_baseline_survival^exp(lps-mean_lp)
}

surv_probs_avg_across_mi = rowMeans(surv_probs)

fake_pt = data.frame("Age"=70-mean_age,"Differentiation"=1,"Number_of_cSCC_before_culprit"=1-mean_ncscc,"PNI_or_LVI"=0,"Sex"=1,"Tissue_involvementSubcutaneous fat"=0,"Tissue_involvementBeyond subcutaneous fat"=0,"Tumor_diameter"=2-mean_tdiam,"Tumor_location_catsScalp/neck"=0 ,"Tumor_location_catsFace"=1)
lp_pt=sum(as.numeric(coefs_final_model*fake_pt))
1-recalibrated_baseline_survival^exp(lp_pt-mean_lp) # Note: the mean lp was subtracted because we had to recalculate the baseline survival. See the github issue that BRP raised in the survival package: https://github.com/therneau/survival/issues/187
#-----------------------

#-----------------------
# 3. Plot risk probabilities
#-----------------------
df_plot = data.frame(SurvProb = surv_probs_avg_across_mi)
df_plot$Metastasis = dutch_ds$Metastasis
df_plot$PALGA_id = dutch_ds$PALGA_excerpt_id
write.csv(df_plot,paste0(dir_int_output,"RSCM 9 weighted cox model probabilities on 390 dataset_updated_centered_LP.csv"))

#Sanity check: AUC should be around 0.78:
roc(df_plot$Metastasis,df_plot$SurvProb)

# 3 options for showing probabilities
df_plot$sign = ifelse(df_plot$Metastasis=="Case",1,-1)
df_plot %>%ggplot(aes(x=(1-SurvProb)*100,fill = Metastasis))+geom_density(alpha=0.4, position = 'identity')+theme_bw()+xlab("Probability of metastasis within 5 years (%)")+ylab("Density")+scale_x_continuous(breaks=(seq(0,100,5)))+scale_fill_manual(values=c("#00BFC4","#F8766D"))
p2=df_plot %>%ggplot(aes(x=(1-SurvProb)*100,fill = Metastasis))+geom_histogram(alpha=0.4, position = 'identity')+theme_bw()+xlab("Probability of metastasis within 5 years (%)")+ylab("Number of patients")+scale_x_continuous(breaks=(seq(0,100,5)))+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"RSCM 9 probability histogram.pdf"),p2,height=5,width=8)
df_plot %>%ggplot(aes(y=(1-SurvProb)*100,x = Metastasis,fill = Metastasis))+geom_boxplot(alpha=0.4)+theme_bw()+ylab("Probability of metastasis within 5 years (%)")+scale_y_continuous(breaks=(seq(0,100,5)))

# Alternative histogram
p4=ggplot(df_plot,aes(x=(1-SurvProb)*100,fill = Metastasis)) + geom_histogram(data=subset(df_plot,Metastasis == 'Control'), alpha = 0.5,aes( y = -..count..))+geom_histogram(data=subset(df_plot,Metastasis == 'Case'),aes( y = ..count..), alpha = 0.5)+theme_bw()+xlab("Probability of metastasis within 5 years (%)")+ylab("Number of patients")+scale_x_continuous(breaks=(seq(0,100,5)))+ylim(c(-125,125))+scale_fill_manual(values=c("#F8766D","#00BFC4"))
ggsave(paste0(dir_int_results,"RSCM 9 probability histogram posneg.pdf"),p4,height=5,width=8)

#-----------------------

#-----------------------
# 4. Sanity check on the computed probabilities
#-----------------------

#Weights:
psamp = read.csv(paste0(dir_int_output,"RSCM 5_3 Sampling weights VitFUP and PA.csv"))
psamp_rordered = psamp[match(dutch_ds$Administratie_nr,psamp$pat.ids),]
stopifnot(length(which(psamp_rordered$Metastasis!=dutch_ds$Metastasis_numeric))==0)
weights = psamp_rordered$weights*12203/sum(psamp_rordered$weights) #Weight rescaling. We use the sum of potential controls + cases = 12203 in the full cohort

# Fit Model:
fm_model="Sex+Age+Tumor_location_cats+Number_of_cSCC_before_culprit+Tumor_diameter+logBreslowThickness+Differentiation+Immunocommpromised_at_time_of_cSCC+PNI_or_LVI+Tissue_involvement+Morphology_subgroup_bin"
fit_cox_weighted = build_pooled_model(dutch_ds_mimp,"wghcox",fm_model,weights)
summary_model_cox_weighted = summary(pool(fit_cox_weighted), conf.int = TRUE)
fit_shrunk_coefs = fit_cox_weighted$analyses[[1]]
fit_shrunk_coefs$coefficients = summary_model_cox_weighted[,c("estimate")]
orig_bas_survival = pool_baseline_survival(dutch_ds_mimp,fit_cox_weighted)

probs_i = matrix(NA,ncol=10,nrow=nrow(dutch_ds)) # probabilities with the survival package
probs_i_manually_computed = matrix(NA,ncol=10,nrow=nrow(dutch_ds)) # probabilities computed with formula
probs_with_pooled_model  = matrix(NA,ncol=10,nrow=nrow(dutch_ds)) # probabilities with the pooled model
for(i in c(1:10)){
  ds_interest = complete(dutch_ds_mimp,i)
  ds_interest$Vitfup_metastasis_years = 5
  dutch_ds_subset_model_form = model.matrix.lm(~.,ds_interest[,c(model_variables,"Metastasis")])[,-1]
  dutch_ds_subset_model_form=dutch_ds_subset_model_form %>% as.data.frame()%>% mutate(Age = Age-fit_cox_weighted$analyses[[i]]$means[["Age"]],
                                                                                      Number_of_cSCC_before_culprit = Number_of_cSCC_before_culprit-fit_cox_weighted$analyses[[i]]$means[["Number_of_cSCC_before_culprit"]],
                                                                                      Tumor_diameter = Tumor_diameter - fit_cox_weighted$analyses[[i]]$means[["Tumor_diameter"]])
  
  probs_i[,i] = predict(fit_cox_weighted$analyses[[i]],ds_interest,type="survival")
  bs_5y = summary(survfit(fit_cox_weighted$analyses[[i]]),times = 5)$surv
  lp = as.vector(t(fit_cox_weighted$analyses[[i]]$coefficients)%*% t(dutch_ds_subset_model_form[,1:10]))
  probs_i_manually_computed[,i] = bs_5y^exp(lp)
  
  #Pooled model:
  dutch_ds_subset_model_form = model.matrix.lm(~.,ds_interest[,c(model_variables,"Metastasis")])[,-1]
  dutch_ds_subset_model_form=dutch_ds_subset_model_form %>% as.data.frame()%>% mutate(Age = Age-mean_age,
                                                                                      Number_of_cSCC_before_culprit = Number_of_cSCC_before_culprit-mean_ncscc,
                                                                                      Tumor_diameter = Tumor_diameter - mean_tdiam)
  lp = as.vector(t(fit_shrunk_coefs$coefficients)%*% t(dutch_ds_subset_model_form[,1:10]))
  probs_with_pooled_model[,i] = orig_bas_survival^exp(lp)
}

plot(rowMeans(probs_i),rowMeans(probs_i_manually_computed))
plot(rowMeans(probs_i),rowMeans(probs_with_pooled_model),xlab="Survival probabilities at 5 years: survival package",ylab= "Survival probabilities computed with the formula pooled model")
text(0.8, 0.2, as.character(round(cor(rowMeans(probs_i),rowMeans(probs_with_pooled_model)),5)),
     cex=0.65, pos=3,col="black") 

# Effect of recalibration:
plot(rowMeans(probs_with_pooled_model),surv_probs_avg_across_mi,xlab="Survival probabilities from the uncalibrated model",ylab="Survival probabilities from the recalibrated model")
abline(0,1)
boxplot(rowMeans(probs_with_pooled_model),surv_probs_avg_across_mi)

# Effect of probabilities before and after fixing bug:
#df_plot_before = read.csv(paste0(dir_results,"/RSCM 9 weighted cox model probabilities on 390 dataset.csv"))
#plot(df_plot_before$SurvProb,surv_probs_avg_across_mi,xlab="Survival probabilities with bug",ylab="Survival probabilities without the bug")

