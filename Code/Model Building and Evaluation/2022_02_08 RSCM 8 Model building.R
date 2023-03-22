#---------------------------------------------------
# Aim: Build an absolute risk model using weighted cox regression and backward selection
# Author: B. Rentroia Pacheco
# Input: Imputed datasets
# Output: weighted Cox Regression Model
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------
folder_name = paste0(dir_int_results,"MICE/v3_0 coding/_ws_","manual","_ovr_","monotone","_nimp_10_cpm_","1e-04")
load(paste0(folder_name,"/RSCM3_1 MICE object.RData")) # multiple imputed dutch dataset (dutch_ds), can be fecthed with dutch_ds_mimp object.
load(paste0(dir_int_output,"RSCM1_1 Dutch Dataset coding_v3_0.RData"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_01_06 RSCM_A 5 Recalibration Functions.R")) #These functions are not used for the final model, but are needed for this script to run
source(paste0(dir_scripts,"/Multiple Imputation/2021_11_30 RSCM3_1 Multiple Imputation MICE Functions.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_01_19 RSCM_A 6 Weighted Performance Metrics functions.R"))
script_id = "RSCM8"
ds_id = "v3_0"

dir.create(paste0(dir_int_results,"Model Performance/"))
#-----------------------

#-----------------------
# 1. Useful functions
#-----------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_02_17 RSCM_8 A1 Internal Validation MI_Boot.R"))
#-----------------------
# 2. Modeling
#-----------------------
# Routinely available variables:
fm_model="Sex+Age+Tumor_location_cats+Number_of_cSCC_before_culprit+Tumor_diameter+logBreslowThickness+Differentiation+Immunocommpromised_at_time_of_cSCC+PNI_or_LVI+Tissue_involvement+Morphology_subgroup_bin"

#Weights for development and evaluation:
psamp = read.csv(paste0(dir_int_output,"RSCM 5_3 Sampling weights VitFUP and PA.csv"))
psamp_rordered = psamp[match(dutch_ds$Administratie_nr,psamp$pat.ids),]
stopifnot(length(which(psamp_rordered$Metastasis!=dutch_ds$Metastasis_numeric))==0)
weights = psamp_rordered$weights*12203/sum(psamp_rordered$weights) #Weight rescaling. We use the sum of potential controls + cases = 12203 in the full cohort

# 2.1 Unweighted Cox regression with backward selection: for comparison with final model, as a sanity check
fit_cox_uns = build_pooled_model(dutch_ds_mimp,"uncox",fm_model,weights)
summary_model_cox_uns = summary(pool(fit_cox_uns), conf.int = TRUE)
exp_estimates=round(exp(summary_model_cox_uns[,c("estimate","2.5 %","97.5 %")]),3)
write.csv(data.frame(Variable=summary_model_cox_uns$term,Estimate=paste0(exp_estimates$estimate,"(",exp_estimates$`2.5 %`,"-",exp_estimates$`97.5 %`,")"),Pvalue=round(summary_model_cox_uns$p.value,3)),paste0(dir_int_results,script_id,"_uncox_model_",ds_id,".csv"))

# 2.2 Unweighted Cox regression with LASSO: for comparison with final model -> Not used
expr <- expression( df  <-  data.frame(Sex,Age,Tumor_location_cats,Number_of_cSCC_before_culprit,Tumor_diameter,logBreslowThickness,Differentiation,Immunocommpromised_at_time_of_cSCC,PNI_or_LVI,Tissue_involvement,Morphology_subgroup_bin),
                    x <- model.matrix(~.,data = df),
                    df  <-  data.frame(Metastasis_numeric,Vitfup_metastasis_years,Sex,Age,Tumor_location_cats,Number_of_cSCC_before_culprit,Tumor_diameter,logBreslowThickness,Differentiation,Immunocommpromised_at_time_of_cSCC,PNI_or_LVI,Tissue_involvement,Morphology_subgroup_bin),
                    df <- truncateFUP(df,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5),
                    Metastasis_numeric <- df[,"Metastasis_numeric"],
                    Vitfup_metastasis_years <- df[,"Vitfup_metastasis_years"],
                    y <- Surv(Vitfup_metastasis_years,Metastasis_numeric),
                    seed=123,
                    fit <- cv.glmnet(x[,-1],y,family="cox",alpha=1,standardize=TRUE))

fit_cox_uns_lasso = with(dutch_ds_mimp,expr)
df_fit_cox_uns_lasso = sapply(1:dutch_ds_mimp$m,function(x) as.matrix(coef(fit_cox_uns_lasso$analyses[[x]],s="lambda.1se")))
keep_vars=rowSums(t(apply(df_fit_cox_uns_lasso,1,function(x) x!=0)))>=5
df_fit_cox_uns_lasso[!keep_vars,]=0
write.csv(data.frame(Variable=rownames(coef(fit_cox_uns_lasso$analyses[[1]],s="lambda.1se")),Estimate=round(exp(rowMeans(df_fit_cox_uns_lasso)),3),Pvalue=NA),paste0(dir_int_results,script_id,"_uncoxLASSO_model_",ds_id,".csv"))

# 2.3a Weighted Cox regression with backward selection: final model

fit_cox_weighted = build_pooled_model(dutch_ds_mimp,"wghcox",fm_model,weights)
summary_model_cox_weighted = summary(pool(fit_cox_weighted), conf.int = TRUE)
exp_estimates_weighted=round(exp(summary_model_cox_weighted[,c("estimate","2.5 %","97.5 %")]),3)
weighted_cox_table_coefs = data.frame(Variable=summary_model_cox_weighted$term,Estimate=paste0(exp_estimates_weighted$estimate,"(",exp_estimates_weighted$`2.5 %`,"-",exp_estimates_weighted$`97.5 %`,")"),Pvalue=round(summary_model_cox_weighted$p.value,3))
write.csv(weighted_cox_table_coefs,paste0(dir_results,script_id,"_wghcox_model_",ds_id,".csv"))
orig_bas_survival = pool_baseline_survival(dutch_ds_mimp,fit_cox_weighted)
# Coefficients across multiple imputed datasets:
coefs_across_imp_datasets = as.data.frame(matrix(NA,ncol=10,nrow=10))
for (i in c(1:10)){
  coefs_across_imp_datasets[,i]=fit_cox_weighted$analyses[[i]]$coefficients
}
rownames(coefs_across_imp_datasets) = names(fit_cox_weighted$analyses[[i]]$coefficients)
df=truncateFUP(dutch_ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5)
fit_complete_weighted  = coxph(Surv(df$Vitfup_metastasis_years,df$Metastasis_numeric)~ Age+Differentiation + Number_of_cSCC_before_culprit + PNI_or_LVI+Sex+Tissue_involvement+Tumor_diameter+Tumor_location_cats,df,weights = weights)
coefs_across_imp_datasets[,"Complete"] = fit_complete_weighted$coefficients
write.csv(round(exp(coefs_across_imp_datasets),3),paste0(dir_results,script_id,"_wghcox_model_all_imp",ds_id,".csv"))

# 3. Evaluation Final Model
# Evaluation is performed using discrimination and calibration measurements.
# Note: Unweighted cox and Recalibrated unweighted cox models were built and evaluated, just to have a reference of what kind of performance could be achieved with these models. They are not the appropriate methodology to analyze these type of cohorts.
#-----------------------
list_perf_models = vector(mode = "list", length = 3)
names(list_perf_models)=c("uncox","wghcox","uncox_recalib")
#Truncation:
dutch_ds = truncateFUP(dutch_ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5)

for (model in names(list_perf_models)){
  print(model)
  if (model == "clogit"){
    fm_model = "Metastasis_numeric ~Sex + Age + Tumor_location_cats + Number_of_cSCC_before_culprit + Tumor_diameter + logBreslowThickness + Differentiation + Immunocommpromised_at_time_of_cSCC + PNI_or_LVI +Tissue_involvement+Morphology_subgroup_bin+strata(Set_id)"
  }else if (model %in%c("uncox","uncoxnet","wghcox","wghcoxnet","uncox_recalib","wghcox_recalibSurv","wghcoxnet_recalibSurv")){
    fm_model = "Surv(Vitfup_metastasis_years,Metastasis_numeric) ~Sex + Age + Tumor_location_cats + Number_of_cSCC_before_culprit + Tumor_diameter+logBreslowThickness + Differentiation + Immunocommpromised_at_time_of_cSCC + PNI_or_LVI +Tissue_involvement+Morphology_subgroup_bin"
  }
  
  if (model %in%c("uncoxnet","wghcoxnet","wghcoxnet_recalibSurv")){
    bw = FALSE
  }else{
    bw=TRUE
  }
  list_perf_models[[model]] = mi_boot(fm_model,model,dutch_ds,bw,weights,100,dutch_ds_mimp,seed=42,5)
}

#3.2 Visualize summary (boxplots): Variability across multiple imputations
summary_internal_validation_folds = as.data.frame(matrix(NA,nrow=0,ncol=25)) #n rows are folds, columns are metrics
for(model in names(list_perf_models)){
  temp_df = cbind(list_perf_models[[model]][[1]],list_perf_models[[model]][[2]])
  temp_df$Model=model
  summary_internal_validation_folds = rbind(summary_internal_validation_folds,temp_df)
}
colnames(summary_internal_validation_folds)=c(colnames(list_perf_models[[model]][[2]]),paste0("TR_",colnames(list_perf_models[[model]][[2]])),"Model")
sum_int_val_folds.m=melt(summary_internal_validation_folds,id="Model")
sum_int_val_folds.m$Split = factor(ifelse(grepl("TR_",sum_int_val_folds.m$variable),"Apparent","Optimism-corrected"),levels=c("Apparent","Optimism-corrected"))
sum_int_val_folds.m$Weighted = ifelse(grepl("^w|^TR_w",sum_int_val_folds.m$variable),"Weighted","Unweighted")
sum_int_val_folds.m$Metric =revalue(gsub("^w|\\.C Index","",gsub("^TR_","", sum_int_val_folds.m$variable)),c("auc"="AUC","cind"= "Cind","itcp_binary"="LR_intercept","slope_binary"="LR_slope","itcp_Surv"="Cox_intercept","slope_Surv" ="Cox_slope"))
sum_int_val_folds.m$BestValue = ifelse(sum_int_val_folds.m$Metric%in%c("AUC","Cind","Cox_slope","LR_slope"),1,0)
sum_int_val_folds.m$Metric=factor(sum_int_val_folds.m$Metric,levels=c("AUC","Cind","LR_intercept","LR_slope","Cox_intercept","Cox_slope"))
gp1=ggplot(sum_int_val_folds.m,aes(x=Model,y=value,alpha=Split,fill=Model))+ geom_boxplot()+facet_grid(Metric~Weighted ,scales="free")+theme_bw()+scale_alpha_ordinal(range = c(0.5, 1))+geom_hline(data=sum_int_val_folds.m,aes(yintercept=BestValue),col="pink",linetype="dashed",lwd=1.2)+scale_fill_manual(values = c("darkgreen","darkolivegreen3","dodgerblue4","deepskyblue3"))+scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave(paste0(dir_int_results,"Model Performance/",script_id," performance estimates variability across MI",ds_id,".pdf"),gp1,width=8,height=8)

# 3.2 Visualize summary (boxplots): Variability across multiple imputations and bootstrap
summary_internal_validation_folds = as.data.frame(matrix(NA,nrow=0,ncol=25)) #n rows are folds, columns are metrics
for(model in names(list_perf_models)){
  temp_df = do.call("rbind", (list_perf_models[[model]][[6]]))
  temp_df$Model=model
  summary_internal_validation_folds = rbind(summary_internal_validation_folds,temp_df)
}
colnames(summary_internal_validation_folds)=c(colnames(list_perf_models[[model]][[2]]),"Model")
sum_int_val_folds.m=melt(summary_internal_validation_folds,id="Model")
sum_int_val_folds.m$Weighted = ifelse(grepl("^w|^TR_w",sum_int_val_folds.m$variable),"Weighted","Unweighted")
sum_int_val_folds.m$Metric =revalue(gsub("^w|\\.C Index","",gsub("^TR_","", sum_int_val_folds.m$variable)),c("auc"="AUC","cind"= "Cind","itcp_binary"="LR_intercept","slope_binary"="LR_slope","itcp_Surv"="Cox_intercept","slope_Surv" ="Cox_slope"))
sum_int_val_folds.m$BestValue = ifelse(sum_int_val_folds.m$Metric%in%c("AUC","Cind","Cox_slope","LR_slope"),1,0)
sum_int_val_folds.m$Metric=factor(sum_int_val_folds.m$Metric,levels=c("AUC","Cind","LR_intercept","LR_slope","Cox_intercept","Cox_slope"))
gp2=ggplot(sum_int_val_folds.m,aes(x=Model,y=value,fill=Model))+ geom_boxplot()+facet_grid(Metric~Weighted ,scales="free")+theme_bw()+scale_alpha_ordinal(range = c(0.5, 1))+geom_hline(data=sum_int_val_folds.m,aes(yintercept=BestValue),col="pink",linetype="dashed",lwd=1.2)+scale_fill_manual(values = c("darkgreen","darkolivegreen3","dodgerblue4","deepskyblue3"))+scale_x_discrete(guide = guide_axis(n.dodge = 2))
ggsave(paste0(dir_int_results,"Model Performance/",script_id," performance estimates variability across MI and Boot",ds_id,".pdf"),gp2,width=8,height=8)


## Add performance estimates at 1 and 3 years for the weighted cox model:
fm_model = "Surv(Vitfup_metastasis_years,Metastasis_numeric) ~Sex + Age+Tumor_location_cats + Number_of_cSCC_before_culprit + Tumor_diameter+logBreslowThickness + Differentiation + Immunocommpromised_at_time_of_cSCC + PNI_or_LVI +Tissue_involvement+Morphology_subgroup_bin"
list_perf_models[["wghcox_1y"]] = mi_boot(fm_model,"wghcox",dutch_ds,TRUE,weights,100,dutch_ds_mimp,seed=42,1)
list_perf_models[["wghcox_3y"]] = mi_boot(fm_model,"wghcox",dutch_ds,TRUE,weights,100,dutch_ds_mimp,seed=42,3)


# Confidence intervals across multiple imputed datasets. The following code combined performance metrics across imputations
# Function pool_auc estimates pooled estimates. It requires the estimated metrics and associated standard errors
metrics_names = colnames(list_perf_models[[1]][["Optimism-corrected of pooled model"]])
performance_with_ci = list()
for (model in names(list_perf_models)){
  temp_df = as.data.frame(matrix(NA,nrow = length(metrics_names),ncol=3))
  rownames(temp_df) = metrics_names
  for (metric in metrics_names){
    if (metric %in%c("auc","wauc","cind.C Index","wcind")){
      logit_tr = TRUE
    }else{
      logit_tr = FALSE
    }
    temp_df[metric,]=as.data.frame(pool_auc_adapted(list_perf_models[[model]][["Optimism-corrected of pooled model"]][,metric],list_perf_models[[model]][["Optimism se"]][,metric],dutch_ds_mimp$m,log_auc = logit_tr))
  }
  colnames(temp_df)=c("lower_CI","Estimate","upper_CI")
  performance_with_ci[[model]]=temp_df
}
sink(paste0(dir_results,script_id," performance estimates 5plus variables across MI",ds_id,".txt"))
print(performance_with_ci)
closeAllConnections()

# Recalibrate model, using the average calibration slope as the shrinkage factor.
exp_estimates_weighted_recalib = round(exp(summary_model_cox_weighted[,c("estimate","2.5 %","97.5 %")]*performance_with_ci[["wghcox"]]["wslope_Surv",2]),2)
weighted_cox_table_coefs_recalib = data.frame(Variable=summary_model_cox_weighted$term,Estimate=paste0(exp_estimates_weighted_recalib$estimate," (",exp_estimates_weighted_recalib$`2.5 %`,"-",exp_estimates_weighted_recalib$`97.5 %`,")"),Pvalue=NA)
write.csv(weighted_cox_table_coefs_recalib,paste0(dir_results,script_id,"_wghcox_model_recalib",ds_id,"rounded.csv"))

exp_estimates_weighted_recalib = round(exp(summary_model_cox_weighted[,c("estimate","2.5 %","97.5 %")]*performance_with_ci[["wghcox"]]["wslope_Surv",2]),7)
weighted_cox_table_coefs_recalib = data.frame(Variable=summary_model_cox_weighted$term,Estimate=paste0(exp_estimates_weighted_recalib$estimate," (",exp_estimates_weighted_recalib$`2.5 %`,"-",exp_estimates_weighted_recalib$`97.5 %`,")"),Pvalue=NA)
write.csv(weighted_cox_table_coefs_recalib,paste0(dir_results,script_id,"_wghcox_model_recalib",ds_id,".csv"))

#Intercept also needs to be recalibrated.
#For this, we obtain the linear predictor with the new coefficients:
fit_shrunk_coefs = fit_cox_weighted$analyses[[1]]
fit_shrunk_coefs$coefficients = summary_model_cox_weighted[,c("estimate")]*performance_with_ci[["wghcox"]]["wslope_Surv",2]
log_log_bsurvival = vector("numeric",length=10)
avg_tumor_diameter = vector("numeric",length=10)
avg_lp = vector("numeric",length=10)
for (m in 1:10){
  ds= complete(dutch_ds_mimp,m)
  lp=get_cox_linear_predictor(fit_shrunk_coefs,ds,ds,NULL)
  #And we fit a model with these
  fit_recalib = coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric)~offset(lp),data=ds,weights=weights)
  #print(summary(survfit(fit_recalib),times = 5)$surv)
  #Pool baseline survival
  log_log_bsurvival[m]=log(-log(summary(survfit(fit_recalib),times = 5)$surv)) # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
  avg_tumor_diameter [m]= fit_cox_weighted$analyses[[m]]$means["Tumor_diameter"]
  avg_lp[m]=mean(lp)
}
recalibrated_baseline_survival = exp(-exp(mean(log_log_bsurvival)))

# Get recalibrated baseline survival for 1,3 and 5 years (slightly different procedure, but yields the same results as above for 5y)
log_log_bsurvival = as.data.frame(matrix(NA,ncol=3,nrow=10))
for (m in 1:10){
  ds= complete(dutch_ds_mimp,m)
  prob=as.vector(get_prob_survival_model(fit_shrunk_coefs,ds,ds,5))
  log_log_pred = log(-log(prob))
  #And we fit a model with these
  fit_recalib = coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric)~offset(log_log_pred),data=ds,weights=weights)
  print(summary(survfit(fit_recalib),times = 5)$surv)
  #Pool baseline survival
  log_log_bsurvival[m,1]=log(-log(summary(survfit(fit_recalib),times = 1)$surv)) # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
  log_log_bsurvival[m,2]=log(-log(summary(survfit(fit_recalib),times = 3)$surv)) # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
  log_log_bsurvival[m,3]=log(-log(summary(survfit(fit_recalib),times = 5)$surv)) # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
}
recalibrated_baseline_survival = exp(-exp(mean(log_log_bsurvival[,3]))) # 5 years
recalibrated_baseline_survival_3y = exp(-exp(mean(log_log_bsurvival[,2]))) # 3 years
recalibrated_baseline_survival_1y = exp(-exp(mean(log_log_bsurvival[,1]))) # 1 years


# Save mean of continuous predictors, so they can be used in the final equation, as well as different recalibrated survival estimates at different timepoints:
important_quantities = c(recalibrated_baseline_survival,means_fit[["Age"]],means_fit[["Number_of_cSCC_before_culprit"]],mean(avg_tumor_diameter),mean(avg_lp),recalibrated_baseline_survival_1y,recalibrated_baseline_survival_3y)
names(important_quantities) = c("Recalibrated baseline survival","Mean age","Mean number of cSCC","Mean tumor diameter","MeanLP_recalib","Baselinesurv_1y","Baselinesurv_3y")
write.csv(important_quantities,paste0(dir_results,script_id,"_wghcox_model_recalib_importantqnts",ds_id,".csv"))

