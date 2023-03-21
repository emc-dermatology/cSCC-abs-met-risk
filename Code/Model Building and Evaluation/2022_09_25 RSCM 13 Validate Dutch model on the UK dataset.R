#---------------------------------------------------
# Aim: Validate Dutch model on the UK dataset
# Author: B. Rentroia Pacheco
# Input: UK dataset
# Output: Performance metrics on the Dutch dataset
#---------------------------------------------------

# 0. Setup
#---------------------------------------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_02_17 RSCM_8 A1 Internal Validation MI_Boot.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_01_19 RSCM_A 6 Weighted Performance Metrics functions.R"))
#---------------------------------------------------

# 1. Load datasets
#---------------------------------------------------
# Decide the sets for the numerator and denominator of the incidence: 
# 2013_onw: only consider patients who do not have any cSCC prior 2013 -> This introduces a bias, as most controls have only 1 prior cSCC, because they are from 2013, and cases are more likely to have multiple prior cSCCs.
# 2013_diag: only consider patients who had a cSCC diagnosed in 2013 -> This was the decision that made more sense, eventhough the cohort is smaller, and therefore CIs are larger.
# all: all available data
weights_scheme = "2013_diag" 
script_id = paste0("RSCM13 ",weights_scheme," ")

# NCC dataset:
load(paste0(dir_int_output,"RSCM12UK_dataset",weights_scheme,".RData"))  
# Weights:
weights_UK = read.csv(paste0(dir_int_output,"/RSCM12 weights_",weights_scheme," UK cohort.csv"))
weights_vec = weights_UK$weights_rescaled
stopifnot(DF_UK$pid==weights_UK$anonpid)
#---------------------------------------------------

#-----------------------
# 2. Load final model built on the Dutch dataset
#-----------------------
weighted_cox_table_coefs_recalib=read.csv(paste0(dir_results,"RSCM8_wghcox_model_recalibv3_0.csv"))
coefs_final_model = log(as.numeric(gsub("\\s.*","",weighted_cox_table_coefs_recalib[,"Estimate"]))) #Log because we need the coefficients for the linear predictor, not the Hazard ratios!
coefs_final_model[1] = coefs_final_model[1]*10 # Coefficent of age is divided by 10, because we decided later that age should be provided in decades.
coefs_final_model=round(coefs_final_model,2)
important_quantities = read.csv(paste0(dir_results,"RSCM8_wghcox_model_recalib_importantqntsv3_0.csv"))
mean_age = round(important_quantities$x[which(important_quantities$X=="Mean age")]/10,2) #/10 because we want effect per decade
mean_ncscc = round(important_quantities$x[which(important_quantities$X=="Mean number of cSCC")],2)
mean_tdiam = round(important_quantities$x[which(important_quantities$X=="Mean tumor diameter")],2)
mean_lp = round(important_quantities$x[which(important_quantities$X =="MeanLP_recalib")],2)
recalibrated_baseline_survival = round(important_quantities$x[which(important_quantities$X =="Baselinesurv_3y")],3)
#-----------------------

#---------------------------------------------------
# 3. Multiple imputation
#---------------------------------------------------
# Only keeping the relevant variables for imputation and application of the clinico-pathological model
subset_vars_UK = c("pid","setID", "site_primary_cSCC_categories","age","Sex","mets","Ethnicity","previous_cSCC_n_continuous","deprivation_IMD", "diameter_num","thickness_num","Differentiation_cat","invasion_bone","tissue_involvement","perineural_invasion","perineural_lymfo_bin","fu_metastasis_yrs","Morphology","immunonew")
DF_UK <- DF_UK[,subset_vars_UK]

# In multiple imputation, follow up should not be used directly, the cumulative hazard rate should be used instead.
DF_UK$mets_num = ifelse(DF_UK$mets=="Case",1,0)
DF_UK$H0 = nelsonaalen(data=DF_UK,timevar = fu_metastasis_yrs, statusvar = mets_num)
DF_UK$mets_num = NULL

#Multiple imputation of missing values.
# setup:
ini_UK <- mice(DF_UK, maxit=0, print=F)
pred_UK <- ini_UK$pred
# Variables that should not be used in the imputation because they are identifiers or for the reason above (follow up time)
pred_UK[ ,"setID"] <- 0
pred_UK[ ,"pid"] <- 0
pred_UK[ ,"fu_metastasis_yrs"] <- 0

mimp = 10
multimp_cat <- mice(DF_UK, m=mimp, maxit=50, meth='pmm', seed = 100, pred=pred_UK)
plot(multimp_cat) # checked convergence, all variables seem fine
#---------------------------------------------------

#---------------------------------------------------
# 4. Compute probability and performance metrics.
#---------------------------------------------------
model_variables = c("Age","Differentiation","Number_of_cSCC_before_culprit","PNI_or_LVI","Sex","Tissue_involvement","Tumor_diameter","Tumor_location_cats")
recalibrated_baseline_survivals = round(important_quantities[,2],3)[c(6,7,1)]
stopifnot(important_quantities[,1][c(6,7,1)]==c("Baselinesurv_1y","Baselinesurv_3y","Recalibrated baseline survival"))

# Auxiliary function for the performance metrics loop:
compute_perf_metrics_validation = function(ds_j,weights_j,model_estimate,b,df,perf_metric,horizon){
  df[b,"AUC"] = roc(ds_j$Metastasis, model_estimate, ci = TRUE,direction=">")$auc
  df[b,"C-ind"] = rcorr.cens(model_estimate,Surv(ds_j[,"fu_metastasis_yrs"],ds_j[,"Metastasis_numeric"]))[[1]]
  df[b,"wAUC"]=WeightedAUC(WeightedROC(model_estimate,ds_j$Metastasis, weight =weights_j))
  df[b,"wC-ind"] =intsurv::cIndex(ds_j[,"fu_metastasis_yrs"],event=ds_j[,"Metastasis_numeric"],model_estimate,weight=weights_j)[[1]]
  
  
  if(perf_metric %in%c("AJCC")){
    df[b,"AUC"] =1-df[b,"AUC"]
    df[b,"C-ind"] = 1-df[b,"C-ind"]
  }else{
    df[b,"wAUC"] =1-df[b,"wAUC"]
    df[b,"wC-ind"] = 1-df[b,"wC-ind"]
  }
  
  
  if(perf_metric %in% c("AbRM","BWH")){
    # Observed and expected ratio # based on https://github.com/danielegiardiello/Prediction_performance_survival/blob/617603162e8fe31d121c99b5c533f053f7923cc4/03_predsurv_extended.md:
    horizon = ifelse(horizon==3,2.95,horizon)
    obj = summary(survfit(Surv(fu_metastasis_yrs, Metastasis_numeric) ~ 1, 
                          data = ds_j), 
                  times = horizon)
    
    df[b,"OE"] = (1 - obj$surv) / (1-mean(model_estimate))
    wobj = summary(survfit(Surv(fu_metastasis_yrs, Metastasis_numeric) ~ 1, 
                           data = ds_j,weights=weights_j), 
                   times = horizon)
    #print(1 - wobj$surv)
    df[b,"wOE"] = (1 - wobj$surv) / (1-weighted.mean(model_estimate,weights_j))
    
    # Calibration slope at fixed time point
    log_log_prob = log(-log(model_estimate))
    df[b,"Cslope"]=coef(coxph(Surv(fu_metastasis_yrs, Metastasis_numeric) ~ log_log_prob,data=ds_j))[1]   
    df[b,"wCslope"]=coef(coxph(Surv(fu_metastasis_yrs, Metastasis_numeric) ~ log_log_prob,data=ds_j,weights = weights_j))[1]  
  }
  return(df)
}

# Computation of performance metrics for complete cases (cc) or the entire dataset "", at different timepoints (1y or 3y)
for(imp in c("cc","")){
  # Performance metrics setup:
  timepoints = c(1,3)
  performance_pooled = as.data.frame(matrix(NA,ncol=17,nrow=3))
  colnames(performance_pooled)=c("N",paste0(rep(c("AUC","C-ind","wAUC","wC-ind","Cslope","wCslope","OE","wOE"),each=2),rep(paste0("-",as.character(timepoints),"y"),times=2)))
  rownames(performance_pooled)=c("AbRM","BWH","AJCC")
  
  
  for(tm in timepoints){
    if(tm ==1){
      baseline_survival = recalibrated_baseline_survivals[1]
    }else if (tm==3){
      baseline_survival = recalibrated_baseline_survivals[2]
    }
    
    perf_metrcs_empydf = as.data.frame(matrix(NA,ncol=2,nrow=multimp_cat$m))
    list_m_metrics = list("AbRM"=list("AUC"= perf_metrcs_empydf,"wAUC"= perf_metrcs_empydf,"C-ind"= perf_metrcs_empydf,"wC-ind"= perf_metrcs_empydf,"Cslope" =perf_metrcs_empydf,"wCslope" =perf_metrcs_empydf,"OE" =perf_metrcs_empydf,"wOE" =perf_metrcs_empydf),
                          "BWH"=list("AUC"= perf_metrcs_empydf,"wAUC"= perf_metrcs_empydf,"C-ind"= perf_metrcs_empydf,"wC-ind"= perf_metrcs_empydf,"Cslope" =perf_metrcs_empydf,"wCslope" =perf_metrcs_empydf,"OE" =perf_metrcs_empydf,"wOE" =perf_metrcs_empydf),
                          "AJCC"=list("AUC"= perf_metrcs_empydf,"wAUC"= perf_metrcs_empydf,"C-ind"= perf_metrcs_empydf,"wC-ind"= perf_metrcs_empydf))
    
    perf_metrcs_empydf_withinstages = as.data.frame(matrix(NA,ncol=5,nrow=multimp_cat$m))
    colnames(perf_metrcs_empydf_withinstages) = c("EntireCohort","T1","T2a","T2b","T3")
    list_m_metrics_withinbwhstages = list("AUC"= perf_metrcs_empydf_withinstages,"wAUC"= perf_metrcs_empydf_withinstages,"C-ind"= perf_metrcs_empydf_withinstages,"wC-ind"= perf_metrcs_empydf_withinstages)
    
    bwh_staging = list()
    ajcc_staging = list()
    
    surv_probs = surv_probs_BWH = matrix(NA,ncol=multimp_cat$m,nrow=nrow(DF_UK))
    cdiffs_comparisons = cvars_comparisons = matrix(NA,ncol=multimp_cat$m,nrow=2)
    rownames(cdiffs_comparisons) = rownames(cvars_comparisons) = c("BWH_vs_ARM","AJCC_vs_ARM")
    
    # Subset complete cases if requested:
    i.incompletecase = is.na(DF_UK$diameter_num)|is.na(DF_UK$perineural_invasion)|is.na(DF_UK$tissue_involvement)|is.na(DF_UK$invasion_bone)|is.na(DF_UK$Differentiation_cat)
    incompletepairs = DF_UK$setID[i.incompletecase] #This is needed for correct bootstrapping with pairs
    i.incompletecase = i.incompletecase|as.character(DF_UK$setID)%in%as.character(incompletepairs) 
    if(imp =="cc"){
      surv_probs = surv_probs[!i.incompletecase,]
      surv_probs_BWH = surv_probs_BWH[!i.incompletecase,]
    }
    
    # Compute performance metrics for each multiple imputed dataset:
    for (m in 1:multimp_cat$m){
      ds = complete(multimp_cat,m)
      
      #Subset dataset if we only want complete cases:
      if(imp =="cc"){
        ds = ds[!i.incompletecase,]
        weights_ds = weights_vec[!i.incompletecase]
      }else{
        weights_ds = weights_vec
      }
      
      
      performance_pooled$N = nrow(ds)
      
      # Code variables in a similar way to Dutch dataset to compute survival probabilities:
      ds = ds %>%dplyr::rename(Age = age,
               Differentiation = Differentiation_cat,
               Number_of_cSCC_before_culprit = previous_cSCC_n_continuous,
               PNI_or_LVI=perineural_lymfo_bin,
               Tissue_involvement=tissue_involvement,
               Tumor_diameter = diameter_num,
               Tumor_location_cats=site_primary_cSCC_categories,
               Metastasis=mets,
               PNI_bin = perineural_invasion,
               Invasion_of_bones = invasion_bone)%>% 
        dplyr::mutate(Age = Age/10,
              PNI_or_LVI = factor(ifelse(as.character(PNI_or_LVI)=="PNI and/or lympho","Yes","No"),levels=c("No","Yes")),
               Tissue_involvement = factor(revalue(Tissue_involvement,c("Invasion in SC fat"="Subcutaneous fat" ,"Invasion beyond SC fat, in underlying tissue"="Beyond subcutaneous fat")),levels=c("Dermis","Subcutaneous fat", "Beyond subcutaneous fat")),
               Tumor_diameter = Tumor_diameter/10,
               Tumor_location_cats = factor(revalue(Tumor_location_cats,c("Trunk + extremities"="Trunk/Extremities")),levels=c("Trunk/Extremities", "Scalp/neck", "Face")),
               Number_of_cSCC_before_culprit=Number_of_cSCC_before_culprit+1,
               PNI_bin = factor(ifelse(PNI_bin==1,"Yes","No"),levels=c("No","Yes")),
               Invasion_of_bones=factor(ifelse(Invasion_of_bones==1,"Yes","No"),levels=c("No","Yes")))
      
      ds$Tumor_diameter=ifelse(ds$Tumor_diameter>4,4,ds$Tumor_diameter)
      ds$Number_of_cSCC_before_culprit=ifelse(ds$Number_of_cSCC_before_culprit>5,5,ds$Number_of_cSCC_before_culprit)
      
      ds_subset_model_form = model.matrix.lm(~.,ds[,c(model_variables,"Metastasis")])[,-1]
      ds_subset_model_form=ds_subset_model_form %>% as.data.frame()%>% mutate(Age = Age-mean_age,
                                                                                          Number_of_cSCC_before_culprit = Number_of_cSCC_before_culprit-mean_ncscc,
                                                                                          Tumor_diameter = Tumor_diameter - mean_tdiam)
      lps = as.vector(t(coefs_final_model)%*% t(ds_subset_model_form[,1:10]))
      surv_probs[,m] = baseline_survival^exp(lps-mean_lp)
      
      
      # Compute BWH stages
      ds=ds %>% mutate(HR_Tumordiam = ifelse(Tumor_diameter>=2,1,0), 
                         HR_Diff = ifelse(Differentiation=="Poor/undifferentiated",1,0),
                         HR_PNI = ifelse(PNI_bin=="Yes",1,0),
                         HR_Tinv = ifelse(Tissue_involvement=="Beyond subcutaneous fat",1,0),
                         Sum_HR = HR_Tumordiam+HR_Diff+HR_PNI+HR_Tinv,
                         Sum_HR_bones = ifelse(as.character(Invasion_of_bones)=="Yes",4,Sum_HR),
                         Decision_BWH = ifelse(Sum_HR>1,1,0)) #Division at T2a
        
        ds$BWH_staging = ifelse(ds$Sum_HR_bones==0,"T1",NA)
        ds$BWH_staging = ifelse(ds$Sum_HR_bones==1,"T2a",ds$BWH_staging)
        ds$BWH_staging = ifelse(ds$Sum_HR_bones%in%c(2,3),"T2b",ds$BWH_staging)
        ds$BWH_staging = ifelse(ds$Sum_HR_bones>=4,"T3",ds$BWH_staging)
        
        bwh_staging[[m]]=ds$BWH_staging
        surv_probs_BWH[,m] = 1-as.numeric(as.character(revalue(ds$BWH_staging,c("T1"="0.007","T2a"="0.045","T2b"="0.367","T3"="0.5"))))
      
      # Compute AJCC  stages
      ds=ds %>% mutate(T1 = ifelse(Tumor_diameter<2 &((PNI_bin=="No"|(PNI_bin=="Yes"&Tissue_involvement=="Dermis"))&Tissue_involvement!="Beyond subcutaneous fat"&thickness_num<=6),1,0),
                         T2 = ifelse(Tumor_diameter>=2 & Tumor_diameter<4 & ((PNI_bin=="No"|(PNI_bin=="Yes"&Tissue_involvement=="Dermis"))&Tissue_involvement!="Beyond subcutaneous fat"&thickness_num<=6),1,0),
                         T3 = ifelse((Tumor_diameter>=4|(PNI_bin=="Yes"&Tissue_involvement!="Dermis")|Tissue_involvement=="Beyond subcutaneous fat"|thickness_num>6)&Invasion_of_bones=="No",1,0),
                         T4 = ifelse(Invasion_of_bones=="Yes",1,0),
                         Decision_AJCC= ifelse(T1==1,0,1))#Division at T1/T2
        
        ds$AJCC_staging = ifelse(ds$T1==1,"T1",NA)
        ds$AJCC_staging = ifelse(ds$T2==1,"T2",ds$AJCC_staging)
        ds$AJCC_staging = ifelse(ds$T3==1,"T3",ds$AJCC_staging)
        ds$AJCC_staging = ifelse(ds$T4==1,"T4",ds$AJCC_staging)
        
        ajcc_staging[[m]]=ds$AJCC_staging
        
        ds$Metastasis_numeric= ifelse(ds$Metastasis=="Case",1,0)
        ds = truncateFUP(ds,"Metastasis_numeric","fu_metastasis_yrs","Metastasis",tm)
        
      #Compare C-index:
        cdiffs_comparisons[1,m] = compareC(ds$fu_metastasis_yrs,ds$Metastasis_numeric,surv_probs_BWH[,m],surv_probs[,m])$est.diff_c
        cdiffs_comparisons[2,m] = compareC(ds$fu_metastasis_yrs,ds$Metastasis_numeric,1-as.numeric(as.character(revalue(ajcc_staging[[m]],c("T1"="0","T2"="1","T3"="2","T4"="3")))),surv_probs[,m])$est.diff_c

        cvars_comparisons[1,m] = compareC(ds$fu_metastasis_yrs,ds$Metastasis_numeric,surv_probs_BWH[,m],surv_probs[,m])$est.vardiff_c
        cvars_comparisons[2,m] = compareC(ds$fu_metastasis_yrs,ds$Metastasis_numeric,1-as.numeric(as.character(revalue(ajcc_staging[[m]],c("T1"="0","T2"="1","T3"="2","T4"="3")))),surv_probs[,m])$est.vardiff_c
        
      # Confidence intervals calculation:
        brep=100
        
        bootstrap_estimates = as.data.frame(matrix(NA,ncol=8,nrow=brep))
        colnames(bootstrap_estimates)=c("AUC", "C-ind","wAUC","wC-ind","Cslope","wCslope","OE","wOE")
        
        list_b_metrics = list("AbRM"=list("Boot_estimates"= bootstrap_estimates,"Boot_sdev"= bootstrap_estimates),
                              "BWH"=list("Boot_estimates"= bootstrap_estimates,"Boot_sdev"= bootstrap_estimates),
                              "AJCC"=list("Boot_estimates"= bootstrap_estimates,"Boot_sdev"= bootstrap_estimates))

        set.seed(123)
        for(b in c(1:brep)){
          # bootstrap sample: note we sample based on pairs, not individual tumors
          #Sample set ids
          set_ids_sampled = sample(unique(as.numeric(as.character(ds$setID))),replace=T)
          #Find indexes of sampled pairs
          j = unlist(purrr::map(set_ids_sampled, function(x) which(ds$setID==x)))
          
          ds_j = ds[j,]
          weights_j = weights_ds[j]
          
          # compute performance metrics
          list_b_metrics[["AbRM"]][["Boot_estimates"]]  = compute_perf_metrics_validation(ds_j,weights_j,surv_probs[j,m],b,list_b_metrics[["AbRM"]][["Boot_estimates"]],"AbRM",tm)
          bwh_model_estimates=1-as.numeric(as.character(revalue(ds_j$BWH_staging,c("T1"="0.007","T2a"="0.045","T2b"="0.367","T3"="0.5"))))
          list_b_metrics[["BWH"]][["Boot_estimates"]]  = compute_perf_metrics_validation(ds_j,weights_j,bwh_model_estimates,b,list_b_metrics[["BWH"]][["Boot_estimates"]],"BWH",tm)
          ajcc_model_estimates=as.numeric(as.character(revalue(ds_j$AJCC_staging,c("T1"="0","T2"="1","T3"="2","T4"="3"))))
          list_b_metrics[["AJCC"]][["Boot_estimates"]]  = compute_perf_metrics_validation(ds_j,weights_j,ajcc_model_estimates,b,list_b_metrics[["AJCC"]][["Boot_estimates"]],"AJCC",tm)
        }
        
        # Means and standard deviations
        for(model_name in names(list_b_metrics)){
          bootstrap_means = apply(list_b_metrics[[model_name]][["Boot_estimates"]],2,mean)
          bootstrap_serrors = apply(list_b_metrics[[model_name]][["Boot_estimates"]],2,sd)
          
          for(metric_name in names(list_m_metrics[[model_name]])){
            list_m_metrics[[model_name]][[metric_name]][m,] = c(bootstrap_means[[metric_name ]],bootstrap_serrors[[metric_name ]])
          }
        }
        
        # Performance of the model within T-stages:
        for(st in c("T1","T2a","T2b","T3")){
          i.st = which(ds$BWH_staging==st)
          list_m_metrics_withinbwhstages[["AUC"]][m,st] = tryCatch(expr = {roc(ds$Metastasis[i.st], surv_probs[i.st,m], ci = TRUE,direction=">")$auc},error = function(e){NA})
          list_m_metrics_withinbwhstages[["C-ind"]][m,st] = rcorr.cens(surv_probs[i.st,m],Surv(ds[i.st,"fu_metastasis_yrs"],ds[i.st,"Metastasis_numeric"]))[[1]]
          
        }
        
    }
    
    # Pool estimates of the performance metrics across multiple imputed datasets:
    
    pdf(paste0(dir_int_results,"/",script_id," distribution perf metrics across mimp UK",imp,".pdf"),width=20,height=8)
    par(mfrow =c(3,8))
    for(model_name in names(list_b_metrics)){
      for(metric_name in names(list_m_metrics[[model_name]])){
        if(metric_name %in%c("AUC","wAUC","C-ind","wC-ind")){
          pooled_estimate = round(pool_auc_adapted(list_m_metrics[[model_name]][[metric_name]][,1],list_m_metrics[[model_name]][[metric_name]][,2],mimp,log_auc = TRUE),2)
          hist(log(list_m_metrics[[model_name]][[metric_name]][,1]),xlab = metric_name,main="")
        }else{
          # Note: Calibration slopes are regression coefficients, so no transformation is performed.
          pooled_estimate = round(pool_auc_adapted(list_m_metrics[[model_name]][[metric_name]][,1],list_m_metrics[[model_name]][[metric_name]][,2],mimp,log_auc = FALSE),2)
          hist(list_m_metrics[[model_name]][[metric_name]][,1],xlab = metric_name,main="")
        }
        performance_pooled[model_name,paste0(metric_name,"-",as.character(tm),"y")] = paste0(pooled_estimate[2]," (",pooled_estimate[1],"-",pooled_estimate[3],")")
      }
    }
    dev.off()
    
    # Save survival probabilities for calibration plots:
    if(tm==1){
      surv_probs_1y=surv_probs
    }else if (tm==3){
      surv_probs_3y = surv_probs
    }
    rownames(surv_probs) = ds$pid
    write.csv(surv_probs,paste0(dir_results,"/",script_id," survival_probs_UK",imp,"_",tm,"y",".csv"))
  }
  write.csv(performance_pooled,paste0(dir_results,"/",script_id," discrimination_performance_UK",imp,".csv"))
  
  write.csv(pnorm(rowMeans(cdiffs_comparisons)/sqrt(rowMeans(cvars_comparisons)))*2,paste0(dir_results,"/",script_id," discrimination_performance_UK_pvalues",imp,".csv"))

}

# Calibration plots:
ds = complete(multimp_cat,1)
ds$Metastasis_numeric= ifelse(ds$mets=="Case",1,0)
ds$"Vitfup_metastasis_years" = ds$fu_metastasis_yrs

# Calibration plots: different flavours
pdf(paste0(paste0(dir_results,"/",script_id," calibrationplots_UK.pdf")),width=20,height=7)
par(mfrow=c(1,3))
n_pts_bp = round(nrow(surv_probs_3y)/10) # This ensures that we have about 10 groups in the calibration plot.
calib_plot(rowMeans(surv_probs_3y),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0,1),TRUE,"survival",weights_vec,"wghcox",new_plot = TRUE,dots_col ="deepskyblue3")
calib_plot(rowMeans(surv_probs_BWH),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0,1),TRUE,"survival",weights_vec,"wghcox",new_plot = FALSE,dots_col = "mediumvioletred")

calib_plot(rowMeans(surv_probs_3y),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0,0.5),TRUE,"survival",weights_vec,"wghcox",surv_or_met = "met",new_plot = TRUE,dots_col ="deepskyblue3")
calib_plot(rowMeans(surv_probs_BWH),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0,0.5),TRUE,"survival",weights_vec,"wghcox",surv_or_met = "met",new_plot = FALSE,dots_col = "mediumvioletred")

calib_plot(rowMeans(surv_probs_3y),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0.0001,1),TRUE,"survival",weights_vec,"wghcox",log_plot = TRUE,new_plot = TRUE,dots_col ="deepskyblue3")
calib_plot(rowMeans(surv_probs_BWH),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0.0001,1),TRUE,"survival",weights_vec,"wghcox",log_plot = TRUE,new_plot = FALSE,dots_col = "mediumvioletred")

dev.off()

# Separate plot for publication:
pdf(paste0(paste0(dir_results,"/",script_id," calibrationplots_UK_1plot.pdf")),width=7,height=7)
par(mfrow=c(1,1))
calib_plot(rowMeans(surv_probs_3y),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0,0.5),TRUE,"survival",weights_vec,"wghcox",surv_or_met = "met",new_plot = TRUE,dots_col ="deepskyblue3")
calib_plot(rowMeans(surv_probs_BWH),ds$Metastasis_numeric,ds$fu_metastasis_yrs,ds,3,n_pts_bp,c(0,0.5),TRUE,"survival",weights_vec,"wghcox",surv_or_met = "met",new_plot = FALSE,dots_col = "mediumvioletred")
dev.off()

df_plot = data.frame(SurvProb = rowMeans(surv_probs_3y))
df_plot$Metastasis = ds$mets
p2=df_plot %>%ggplot(aes(x=(1-SurvProb)*100,fill = Metastasis ))+geom_histogram(alpha=0.4, position = 'identity')+theme_bw()+xlab("Probability of metastasis within 5 years (%)")+ylab("Number of patients")+scale_x_continuous(breaks=(seq(0,100,5)))+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,script_id," probability histogram.pdf"),p2,height=5,width=8)
#---------------------------------------------------

# 5. Compute most frequent BWH and AJCC stages across multiple imputed datasets
#---------------------------------------------------
#BWH:
bwh_staging_df = t(as.data.frame(matrix(unlist(bwh_staging),nrow=length(bwh_staging),byrow=TRUE)))
#readapted from : https://stackoverflow.com/questions/28094468/how-to-calculate-the-mode-of-all-rows-or-columns-from-a-dataframe
modefunc <- function(x){
  set.seed(123) # In case we have the half of the iterations have the same category and the other have another
  tabresult <- table(x)
  themode <- names(tabresult)[which(tabresult == max(tabresult))]
  if(sum(tabresult == max(tabresult))>1) themode <- sample(themode,1)
  return(themode)
}
bwh_staging_avg = apply(bwh_staging_df,1,modefunc)
bwh_staging_avg = factor(as.character(bwh_staging_avg),levels=c("T1","T2a","T2b","T3"))

#AJCC:
ajcc_staging_df = t(as.data.frame(matrix(unlist(ajcc_staging),nrow=length(ajcc_staging),byrow=TRUE)))
ajcc_staging_avg = apply(ajcc_staging_df,1,modefunc)
ajcc_staging_avg = factor(as.character(ajcc_staging_avg),levels=c("T1","T2","T3","T4"))

# Save Stages in the UK cohort:
i.incompletecase_BWH = is.na(DF_UK$diameter_num)|is.na(DF_UK$perineural_invasion)|is.na(DF_UK$tissue_involvement)|is.na(DF_UK$Differentiation_cat)
DF_UK$BWH_staging = bwh_staging_avg
DF_UK$BWH_staging[i.incompletecase_BWH]=NA 

i.incompletecase_AJCC = is.na(DF_UK$diameter_num)|is.na(DF_UK$perineural_invasion)|is.na(DF_UK$tissue_involvement)|is.na(DF_UK$Differentiation_cat)|is.na(DF_UK$thickness_num)
DF_UK$AJCC_staging = ajcc_staging_avg
DF_UK$AJCC_staging[i.incompletecase_AJCC ]=NA 

st_syst_UK=DF_UK[,c("AJCC_staging","BWH_staging","mets","setID")] %>%
  sjlabelled::remove_all_labels() %>%
  tibble %>%
  tbl_summary(
    by = mets, type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"),include = -setID) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(list(all_continuous() ~ "paired.wilcox.test",all_categorical()~"mcnemar.test"),group=setID) %>% # test for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()

st_syst_UK %>%
  as_flex_table() %>%
  flextable::save_as_docx(path=paste0(dir_results,script_id," Staging systems stratified UK no mets baseline.docx"))

#---------------------------------------------------

# 6. Plot distribution of probabilities within BWH stages and within AJCC
#---------------------------------------------------
ds = complete(multimp_cat,1)
df_plot = data.frame(Met_risk = 1 - rowMeans(surv_probs),"BWH_Staging" = bwh_staging_avg,Metastasis = ds$mets,"AJCC_Staging"=ajcc_staging_avg)
plot_boxplot_UK=ggplot(df_plot,aes(y=100*Met_risk,x=BWH_Staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("BWH stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/",script_id," estimated prob by BWH staging UK dataset.pdf"),plot_boxplot_UK,width=12,height=7)
plot_boxplot_UK_cc=ggplot(df_plot[!i.incompletecase,],aes(y=100*Met_risk,x=BWH_Staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("BWH stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/",script_id," estimated prob by BWH staging UK dataset CC.pdf"),plot_boxplot_UK,width=12,height=7)

plot_boxplot_UK_AJCC=ggplot(df_plot,aes(y=100*Met_risk,x=AJCC_Staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("AJCC8 stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/",script_id," estimated prob by AJCC staging UK dataset.pdf"),plot_boxplot_UK_AJCC,width=12,height=7)
plot_boxplot_UK_AJCC_cc=ggplot(df_plot[!i.incompletecase,],aes(y=100*Met_risk,x=AJCC_Staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("AJCC8 stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/",script_id," estimated prob by AJCC staging UK dataset CC.pdf"),plot_boxplot_UK_AJCC_cc,width=12,height=7)

#---------------------------------------------------


# 7. Decision curve analysis for different prevalence:
#---------------------------------------------------
# Net benefit: All models with cutoffs:
df_plot$Metastasis_numeric = ifelse(ds$mets=="Case",1,0)
df_plot$fup_metastasis_years = ds$fu_metastasis_yrs
ds =df_plot
ds = ds %>% mutate(BWH_risk_estimates = as.numeric(as.character(revalue(BWH_Staging,c("T1"=0.007,"T2a"=0.045,"T2b"=0.367,"T3"=0.5)))),
                   BWH_T1vsrest = ifelse(BWH_Staging=="T1",0,1),
                   BWH_T2avsT2b = ifelse(BWH_Staging%in%c("T1","T2a"),0,1),
                   Decision_AJCC = ifelse(AJCC_Staging=="T1",0,1))

vars=c("Met_risk","BWH_risk_estimates","Decision_AJCC")
col_vars = c("deepskyblue3","mediumvioletred","seagreen")
var_labels = c("Absolute risk model","BWH","AJCC:T1 vs rest")
names(col_vars)=names(var_labels)=vars

# Net benefit: All models with probability:
stopifnot(baseline_survival==round(important_quantities[which(important_quantities[,1]=="Baselinesurv_3y"),2],3))

for( tm in timepoints){
  if(tm==1){
    ds$Met_risk = 1-rowMeans(surv_probs_1y)
  }else if (tm==3){
    ds$Met_risk = 1-rowMeans(surv_probs_3y)
  }
  
  for (imp in c("","cc")){
    if(imp ==""){
      ds_dca = ds
      weights_ds = weights_vec
    }else if (imp=="cc"){
      ds_dca = ds[!i.incompletecase,]
      weights_ds = weights_vec[!i.incompletecase]
    }
    
    ds_dca = truncateFUP(ds_dca,"Metastasis_numeric","fup_metastasis_years","Metastasis",tm)
    #DCA with assumed prevalence, this is to to have a reference of the decision curve analysis made with an R package.
    # However, this does not take censoring and sampling procedure into account.
    p_list=list()
    seq_prevs = c(0.02,0.05)
    for(i in c(1:length(seq_prevs))){
      p_list[[i]]=dca(Metastasis_numeric ~ Met_risk+BWH_risk_estimates+Decision_AJCC, 
                      data = ds_dca, 
                      thresholds = seq(0, 0.20, by = 0.01),
                      label = list(Met_risk = "Absolute risk model",BWH_risk_estimates = "BWH",AJCC_decision = "AJCC T1 vs rest"),
                      prevalence=seq_prevs[i]) %>%
        plot()+ ggtitle(paste0("Assumed prevalence = ",as.character(round(seq_prevs[i],3)))) 
    }
    grid_dca_prob =  grid.arrange(arrangeGrob(grobs=p_list,ncol=2,nrow=1))
    ggsave(paste0(dir_results,"/",script_id," DCA_UK_assumed prevalence",imp,"_",tm,"y",".pdf"),width=16,height=5,grid_dca_prob)
  
    #Probabilities:
    pdf(paste0(dir_results,"/",script_id," DCA_UK_surv_weighted ",imp,"_",tm,"y",".pdf"),width=7,height=6)
    if(imp ==""){
      plot_dca_surv_wgh(ds_dca,vars,col_vars,ds_dca$Metastasis_numeric,ds_dca$fup_metastasis_years,c(-0.01,0.02),c(0,20),5,weights_ds,var_labels)
    }else if(imp =="cc"){
      plot_dca_surv_wgh(ds_dca,vars,col_vars,ds_dca$Metastasis_numeric,ds_dca$fup_metastasis_years,c(-0.01,0.02),c(0,20),5,weights_ds,var_labels)
    }
    dev.off()
  }
}
#---------------------------------------------------

# 7. Comparison different cutoffs:
#---------------------------------------------------
th_performance = as.data.frame(matrix(NA,ncol=6,nrow=8))
colnames(th_performance) = c("Sensitivity","Specificity","NPV","PPV","NCase","NControl")
rownames(th_performance) = c("BWH_T1vsrest","BWH_T2avsT2b","Decision_AJCC","Abs_1p","Abs_2p","Abs_3p","Abs_4p","Abs_5p")

#Make new rows for different cutoffs
for (i in seq(1,5,1)){
  ds[,paste0("Abs_",as.character(i),"p")]=ifelse(ds$Met_risk>i/100,1,0)
}

for(i in rownames(th_performance)){
  th_performance[i,"Sensitivity"]=sensitivity(factor(ds[,i],levels=c("0","1")),factor(ds$Metastasis_numeric,levels=c("0","1")),positive = "1")
  th_performance[i,"Specificity"]=specificity(factor(ds[,i],levels=c("0","1")),factor(ds$Metastasis_numeric,levels=c("0","1")),negative = "0")
  th_performance[i,"NPV"]=negPredValue(factor(ds[,i],levels=c("0","1")),factor(ds$Metastasis_numeric,levels=c("0","1")),negative = "0",prevalence = 0.03)
  th_performance[i,"PPV"]=posPredValue(factor(ds[,i],levels=c("0","1")),factor(ds$Metastasis_numeric,levels=c("0","1")),positive = "1",prevalence = 0.03)
  freq_pts = table(factor(ds[,i],levels=c("0","1")))
  th_performance[i,"NCase"] = freq_pts["1"]
  th_performance[i,"NControl"] = freq_pts["0"]
}
th_performance$PctHRiskclass = th_performance$NCase/(nrow(DF_UK))
df_plot_pm = melt(as.matrix(th_performance))

p_pm = ggplot(df_plot_pm[!df_plot_pm$Var2%in%c("NCase","NControl"),], aes(x = Var2, y = value*100,col=Var1)) +scale_color_manual(values=c("deepskyblue","deepskyblue4","orange","olivedrab1","olivedrab2","olivedrab3","limegreen","olivedrab4")) +
  geom_point(size = 5) +
  geom_line(aes(group =  Var1))+theme_bw()+xlab("Metrics")+ylab("Value (%)")
ggsave(paste0(dir_int_results,"/",script_id," Perf Metrics UK.pdf"),width=16,height=5,p_pm)

#---------------------------------------------------
