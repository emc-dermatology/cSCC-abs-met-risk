#Function bootstrap
truncateFUP = function(df,outcome_num,outcome_fup,outcome,th){
  i.truncated = which(df[,outcome_fup]>th)
  df[i.truncated,outcome_fup]=th
  df[i.truncated,outcome_num]=0
  df[i.truncated,outcome]="Control"
  return(df)
}

mi_boot = function(fm_model,model_method,ds,bw,weights,B,ds_mimp,seed,tm_pt){
  set.seed(seed) #for reproducibility
  #1. Build pooled model
  mice_models = build_pooled_model(ds_mimp,model_method,gsub(".*~","",fm_model),weights)
  
  #Pool baseline survival
  log_log_bsurvival = vector("numeric",length=ds_mimp$m)
  
  if(!model_method%in%c("wghcoxnet_median","wghcoxnet_average","wghcoxnet_recalibscore")){
    mice_models_pooled = mice::pool(mice_models)
    #Pool coefficients
    fit_train = mice_models$analyses[[1]] 
    fit_train$coefficients =mice_models_pooled$pooled$estimate
    
    #Extract each baseline survival
    for(m in c(1:ds_mimp$m)){
      log_log_bsurvival[m]=log(-log(summary(survfit(mice_models$analyses[[m]]),times = tm_pt)$surv)) # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
    }
  }else{
    ds <- truncateFUP(ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5)
    Metastasis_numeric <- ds[,"Metastasis_numeric"]
    Vitfup_metastasis_years <- ds[,"Vitfup_metastasis_years"]
    
    fake_fit = coxph(as.formula(fm_model),data=ds,weights=weights) #This is just for us to have an object we can use for functions get_prob_survival and get_performance_metrics
    
    #This fake fit is populated according to the pooled coefficients:
    stopifnot(row.names(coef(mice_models[[1]]))==names(coef(fake_fit)))
    
    df_fit_cox_wgh_lasso=sapply(1:length(mice_models),function(x) as.matrix(coef(mice_models[[x]],s="lambda.1se")))
    
    if(model_method =="wghcoxnet_median"){
      fake_fit$coefficients = apply(df_fit_cox_wgh_lasso,1,median)
    }else if (model_method =="wghcoxnet_average"){
      fake_fit$coefficients = apply(df_fit_cox_wgh_lasso,1,mean)
    }
    
    for(i in c(1:ds_mimp$m)){
      ds_i = complete(ds_mimp,i)
      ds_i$Mitotic_Rate_ws <- ifelse(ds_i$Mitotic_rate>2,2,ds_i$Mitotic_rate)
      ds_i$Mitotic_rate<-NULL
      ds_x = ds_i[,strsplit(gsub("\\+",",",gsub(".*~|\\s","",fm_model)),",")[[1]]]
      x <- model.matrix(~.,data = ds_x,na.action="na.pass") #na.action is needed, otherwise rows with NAs are removed
      #lp = predict(fit[[i]],newx=data.matrix(x[,-1]),s="lambda.1se",type="link") #This gives the linear predictor (product of coefficients and values)
      #fit_recalib = coxph(Surv(Vitfup_metastasis_years,Metastasis_numeric) ~lp,x=TRUE,weights=weights)
      
      #Extract each baseline survival
      fit_m=mice_models[[i]]
      log_log_bsurvival[i]=log(-log(summary(survfit(fit_m,x =x[,-1],y=Surv(Vitfup_metastasis_years,Metastasis_numeric),weights=weights),times = tm_pt)$surv))
      
    }
    fit_train=fake_fit
  }
  
  baseline_survival_pooled = exp(-exp(mean(log_log_bsurvival)))
  
  #2. Compute optimism corrected metrics:
  optimism_estimates = optimism_estimates_se=  bootstrap_corrected_estimates = perf_estimates  = as.data.frame(matrix(NA,ncol=12,nrow=ds_mimp$m))
  colnames(perf_estimates) = colnames(optimism_estimates_se)=colnames(optimism_estimates)=colnames(bootstrap_corrected_estimates ) = c("auc","wauc","cind.C Index","wcind","itcp_binary","slope_binary","itcp_Surv","slope_Surv","witcp_binary","wslope_binary","witcp_Surv","wslope_Surv")
  all_corrected_estimates = list()
  selected_vars = list()
  for(m in c(1:ds_mimp$m)){
    print(paste0("Computing performance estimates for: imputation",m))
    complete_i = complete(ds_mimp,m)
    # We will add winsorization of mitotic rate here, because this way we don't change the multiple imputation
    complete_i$Mitotic_Rate_ws=ifelse(complete_i$Mitotic_rate>2,2,complete_i$Mitotic_rate)
    prob_i = get_prob_survival_model(fit_train,complete_i,complete_i,5,baseline_survival_pooled,model_method)[1,] #Note: the baseline survival is the one that influences the survival probabilities computed, so that's why tm_pt is not entered as a parameter.
    # baseline_survival_pooled is used, instead of the survival of model i, because we use the pooled baseline survival in our final model, so we wanted to know whether this can for instance, lead to a bad miscalibration of the model. The differences are minimal though.
    if(model_method == "uncox_recalib"){
      recalib.utils = recalibrate_surv_probs(prob_i,complete_i ,weights)
      predictions_i_recalib=recalib.utils [["predictions_tr"]]
      # Assess performance
      perf_estimates[m,] = get_performance_measures(predictions_i_recalib,recalib.utils[["log_log_predictions_tr"]],"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",recalib.utils[["fit"]],model_method,weights,tm_pt)
    }else{
      # Assess performance
      perf_estimates[m,] = get_performance_measures(prob_i,complete_i,"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",fit_train,model_method,weights,tm_pt)
      
    }
    #compute optmisim
    bootstrap_estimate = bootstrap_validation(fm_model,model_method,complete_i,bw,weights,B,ds_mimp,tm_pt)
    optimism_estimates[m,] = bootstrap_estimate[["Optimism"]]
    bootstrap_corrected_estimates [m,] = bootstrap_estimate[["Optimism-Corrected"]]
    allbootstrap_estimates = bootstrap_estimate[["All bootstrap optimism"]]
    optimism_estimates_se[m,] = apply(allbootstrap_estimates,2,function(x)sqrt(var(x))) # standard error
    all_corrected_estimates[[m]]=as.data.frame(t(t(rep(1, B))) %*% as.numeric(perf_estimates[m,])) - allbootstrap_estimates
    selected_vars[[m]]=bootstrap_estimate[["Selected_Features"]]
  }
  ret_list = list("Optimism-corrected of pooled model" = perf_estimates-optimism_estimates,"Apparent Performance of pooled model"=perf_estimates,"Optimism over imputed datasets" = optimism_estimates,"Corrected estimates over imputed datasets"=bootstrap_corrected_estimates,"Optimism se"=optimism_estimates_se,"ALL corrected estimates"=all_corrected_estimates,"baseline_hazard_pooled"=baseline_survival_pooled,"Selected variables"=selected_vars)
  return(ret_list)
}


# this function performs bootstrap validation of a model fitted on data ds. Backward selection can be applied.
bootstrap_validation = function(fm_model,model_method,ds,bw,weights,B,ds_mimp,tm_pt){
  #0. Setup
  #set.seed(1)
  #fm_train <- as.formula(fm_model)
  ds <- truncateFUP(ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5)
  
  #1. Construct a model in the entire dataset ds
  fit_train <-build_model(fm_model,ds,model_method,weights,bw,ds_mimp,tm_pt=tm_pt)
  if(model_method %in% c("wghcoxnet_median","wghcoxnet_average")){
    predictions_model_orig = as.numeric(get_prob_survival_model(fit_train[["model"]],ds,ds,5,fit_train[["baseline_hazard"]],model_method)) #Note: the baseline survival is the one that influences the survival probabilities computed, so that's why tm_pt is not entered as a parameter.
  }else{
    predictions_model_orig = as.numeric(get_prob_survival_model(fit_train[["model"]],ds,ds,tm_pt))
  }
  if(model_method == "uncox_recalib"){
    recalib.utils = recalibrate_surv_probs(predictions_model_orig,ds,weights)
    predictions_model_orig=recalib.utils [["predictions_tr"]]
    # Assess performance
    perf_orig_on_orig = get_performance_measures(predictions_model_orig,recalib.utils[["log_log_predictions_tr"]],"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",recalib.utils[["fit"]],model_method,weights,tm_pt)
  }else{
    # Assess performance
    perf_orig_on_orig = get_performance_measures(predictions_model_orig,ds,"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",fit_train[["model"]],model_method,weights,tm_pt)
    
  }
  
  optimism_boostraps = as.data.frame(matrix(NA,ncol=12,nrow=B))
  total_feats = strsplit(gsub(".*~|\\s","",fm_model),"\\+")[[1]]
  selected_feats = as.data.frame(matrix(0,ncol=B,nrow=length(total_feats)))
  rownames(selected_feats)=total_feats
  #Check that dataframe has set id information. This will be used when sampling pairs for the bootstrapped samples
  stopifnot("Set_id"%in%colnames(ds))
  for (brep in c(1:B)){
    # Boostrap repetitions:
    #2. Draw a bootstrap sample with replacement from the sample
    # j = sample(nrow(ds),replace=T) #this would be for individual tumors
    #Sample set ids
    set_ids_sampled = sample(unique(as.numeric(as.character(ds[,"Set_id"]))),replace=T)
    #Find indexes of sampled pairs
    j = unlist(purrr::map(set_ids_sampled, function(x) which(ds[,"Set_id"]==x)))
    
    # TODO: Perhaps I should have stratified this. (bootstrap within cases and controls)
    #3. Build a model in the bootstrap sample
    fit_bootstrap = build_model(fm_model,ds[j,],model_method,weights[j],bw,ds_mimp,j,tm_pt=tm_pt)
    selected_feats[fit_bootstrap[["bck_feats"]],brep]=1
    #4. Apply model in the bootstrap sample
    if(model_method %in% c("wghcoxnet_median","wghcoxnet_average")){
      predictions_model= as.numeric(get_prob_survival_model(fit_bootstrap[["model"]],ds[j,],ds[j,],5,fit_bootstrap[["baseline_hazard"]],model_method))
    }else{
      predictions_model= as.numeric(get_prob_survival_model(fit_bootstrap[["model"]],ds[j,],ds[j,],tm_pt))}
    if(model_method == "uncox_recalib"){
      recalib.utils = recalibrate_surv_probs(predictions_model,ds[j,],weights[j])
      predictions_model_recalib=recalib.utils [["predictions_tr"]]
      #4.1 get all prediction measures
      perf_bsmodel_on_bs = get_performance_measures(predictions_model_recalib,recalib.utils[["log_log_predictions_tr"]],"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",recalib.utils[["fit"]],model_method,weights[j],tm_pt)
    }else{
      
      perf_bsmodel_on_bs = get_performance_measures(predictions_model,ds[j,],"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",fit_bootstrap[["model"]],model_method,weights[j],tm_pt)
      
      
    }
    if(model_method %in% c("wghcoxnet_median","wghcoxnet_average")){
      #5. Apply model to the original sample
      predictions_bs_model_on_orig = as.numeric(get_prob_survival_model(fit_bootstrap[["model"]],ds[j,],ds,5,fit_bootstrap[["baseline_hazard"]],model_method))
    }else{
      predictions_bs_model_on_orig = as.numeric(get_prob_survival_model(fit_bootstrap[["model"]],ds[j,],ds,tm_pt))
      
    }
    
    if(model_method == "uncox_recalib"){
      recalib.utils = recalibrate_surv_probs(predictions_model,ds[j,],weights[j],predictions_bs_model_on_orig,ds)
      perf_bsmodel_on_orig = get_performance_measures(recalib.utils[["predictions_test"]],recalib.utils[["log_log_predictions_tst"]],"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",recalib.utils[["fit"]],model_method,weights,tm_pt)
    }else{
      perf_bsmodel_on_orig = get_performance_measures(predictions_bs_model_on_orig,ds,"Metastasis","Metastasis_numeric","Vitfup_metastasis_years",fit_bootstrap[["model"]],model_method,weights,tm_pt)
    }
    #6. Compute optimism   
    optimism = unlist(perf_bsmodel_on_bs)-unlist(perf_bsmodel_on_orig)
    optimism_boostraps[brep,]=optimism
  }
  
  #7. Optimism corrected metrics
  optimism_estimates = colMeans(optimism_boostraps)
  int_val_estimates = unlist(perf_orig_on_orig)-optimism_estimates 
  return(list("Optimism-Corrected"=int_val_estimates,"Optimism" = optimism_estimates,"All bootstrap optimism" = optimism_boostraps,"Selected_Features"=selected_feats))
  
}

build_model = function(fm_train,ds,model_method,weights,bw,ds_mimp=NULL,i.subset=NULL,tm_pt=5){
  #NOTE: i.subset is needed only for ds_mimp
  bck_feats = NA
  fm_train = as.formula(as.character(fm_train))
  if(bw){
    if (is.null(ds_mimp)){
      #if(!is.null(i.subset)){ds=ds[i.subset,] } Note: this was not a bug: i subset is only for the multiple imputed dataset
      if(model_method =="wghcox"){
        bw_select = fastbw(cph(fm_train,x=TRUE,y=TRUE,data=ds,weights = weights),rule= "aic")
      }else if (model_method %in%c("uncox","uncox_recalib")){
        bw_select = fastbw(cph(fm_train,x=TRUE,y=TRUE,data=ds),rule= "aic")
      }
      bck_feats = bw_select$names.kept
    }else{
      bck_feats = list()
      for (m in 1:ds_mimp$m){
        data_compl <- as.data.frame(complete(ds_mimp, m))
        # We will add winsorization of mitotic rate here, because this way we don't change the multiple imputation
        data_compl$Mitotic_Rate_ws=ifelse(data_compl$Mitotic_rate>2,2,data_compl$Mitotic_rate)
        
        if(!is.null(i.subset)){data_compl=data_compl[i.subset,] }
        
        if(model_method =="wghcox"){
          bw_select = fastbw(cph(fm_train,x=TRUE,y=TRUE,data=data_compl,weights = weights),rule= "aic")
        }else if (model_method %in%c("uncox","uncox_recalib")){
          bw_select = fastbw(cph(fm_train,x=TRUE,y=TRUE,data=data_compl),rule= "aic")
        }
        bck_feats[[m]] =bw_select$names.kept
        
      }
      chosen_bck = table(unlist(bck_feats))
      bck_feats = names(chosen_bck)[which(chosen_bck > ds_mimp$m/2)]
      
    }
    #update formula
    fm_model = paste0(gsub("~.*","",fm_train)[2],"~",paste(bck_feats,collapse = "+")) 
    fm_train = as.formula(fm_model) # Note: the estimate that are recomputed are essentially the same. The pvalues/standard errors are the ones that differ, predictions are the same.
  }
  if (model_method =="wghcox"){
    if (!bw){fm_train = as.formula(fm_train)}
    fit_train = coxph(fm_train, data = ds,weights = weights)
    baseline_hazard = summary(survfit(fit_train),times = tm_pt)$surv
  }
  if (model_method %in%c("uncox","uncox_recalib")){
    if (!bw){fm_train = as.formula(fm_train)}
    fit_train = coxph(fm_train, data = ds)
    baseline_hazard = summary(survfit(fit_train),times = tm_pt)$surv
  }
  if( model_method %in% c("wghcoxnet_median","wghcoxnet_average","wghcoxnet_recalibvar")){
    # Fit the model to select features:
    # I got stuck in how to access the column names within the expression environment without hard coding them, so I just coded the loop across the imputed datasets
    mice_models=list()
    log_log_bsurvival = vector("numeric",length=ds_mimp$m)
    for (i in c(1:ds_mimp$m)){
      complete_i = as.data.frame(complete(ds_mimp,i))
      complete_i$Mitotic_Rate_ws <- ifelse(complete_i$Mitotic_rate>2,2,complete_i$Mitotic_rate)
      complete_i$Mitotic_rate<-NULL
      if(!is.null(i.subset)){complete_i=complete_i[i.subset,] }
      predictors = strsplit(gsub("\\+",",",gsub(".*~|\\s","",fm_train))[[3]], ",")[[1]]
      fit_m=train_wghcoxnet(complete_i,10,predictors,weights,tm_pt=tm_pt)
      mice_models[[i]] = fit_m[["fit_train"]]
      log_log_bsurvival[i]=log(-log(fit_m[["baseline_survival"]]))
    }
    # Pool model:
    
    if (model_method %in%c("wghcoxnet_median","wghcoxnet_average")){
      fake_fit = coxph(fm_train,data=ds,weights=weights) #This is just for us to have an object we can use for functions get_prob_survival and get_performance_metrics
      
      #This fake fit is populated according to the pooled coefficients:
      stopifnot(row.names(coef(mice_models[[1]]))==names(coef(fake_fit)))
      
      df_fit_cox_wgh_lasso=sapply(1:length(mice_models),function(x) as.matrix(coef(mice_models[[x]])))
      
      if(model_method =="wghcoxnet_median"){
        fake_fit$coefficients = apply(df_fit_cox_wgh_lasso,1,median)
      }else if (model_method =="wghcoxnet_average"){
        fake_fit$coefficients = apply(df_fit_cox_wgh_lasso,1,mean)
      }
      fit_train = fake_fit
      
      baseline_hazard = exp(-exp(mean(log_log_bsurvival)))
      #print(names(fake_fit$coefficients))
      bck_feats = row.names(coef(mice_models[[1]]))[which(fake_fit$coefficients!=0)]
      
      bck_feats = unique(gsub("man|Scalp/neck|Face|Poor/.*|Abundant|Yes|Extensive|R1/R2.*|Beyond.*|Subcutane.*|No/.*","",bck_feats)) #Otherwise, variables with multiple selected levels will be counted multiple times.
      
    }else if(model_method %in%"wghcoxnet_recalibvar"){
      bck_feats = get_most_selected_bck(mice_models,"wghcoxnet_recalibvar")[[2]]
      #update formula
      fm_model = paste0(gsub("~.*","",fm_train)[2],"~",paste(bck_feats,collapse = "+")) 
      fm_train = as.formula(fm_model) # Note: the estimate that are recomputed are essentially the same. The pvalues/standard errors are the ones that differ, predictions are the same.
      fit_train = coxph(fm_train, data = ds,weights = weights)
      baseline_hazard = summary(survfit(fit_train),times = tm_pt)$surv
      
    }
  }
  
  
  return(list("model"=fit_train,"baseline_hazard" = baseline_hazard,"bck_feats"=bck_feats))
}

get_performance_measures=function(surv_probabilities,ds,outcome,outcome_num,outcome_FU,model_fit,model_method,weights,tm_pt){
  
  ds = truncateFUP(ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",tm_pt)
  
  # Compute metrics
  auc_model = pROC::roc(ds[,outcome], surv_probabilities, quiet = TRUE,direction = ">")$auc 
  wauc_model = WeightedAUC(WeightedROC(1-surv_probabilities, ds[,outcome],weight=weights))
  
  cind_model = rcorr.cens(surv_probabilities,Surv(ds[,outcome_FU],ds[,outcome_num]))[1] 
  wcind_model = intsurv::cIndex(ds[,outcome_FU],event=ds[,outcome_num],1-surv_probabilities,weight=weights)[[1]]
  
  res_tr=calib_plot(surv_probabilities, ds[,outcome_num],ds[,outcome_FU],ds,tm_pt,20,c(0,1),FALSE,outcome_type="survival",weights=NULL,model_method) 
  wres_tr=calib_plot(surv_probabilities, ds[,outcome_num],ds[,outcome_FU],ds,tm_pt,20,c(0,1),FALSE,outcome_type="survival",weights=weights,model_method) 
  
  
  #Store metrics
  return(list(auc = auc_model,
              wauc = wauc_model,
              cind = cind_model,
              wcind = wcind_model,
              itcp_binary = res_tr[[1]],
              slope_binary = res_tr[[2]],
              itcp_Surv =res_tr[[3]] ,
              slope_Surv =res_tr[[4]] ,
              witcp_binary = wres_tr[[1]],
              wslope_binary =wres_tr[[2]],
              witcp_Surv = wres_tr[[3]],
              wslope_Surv = wres_tr[[4]]
  ))
}

remove_least_freq_selected_vars = function(fit,ds_mimp,sel_model,weights=NULL,fm_model){
  bck_feats_sel = get_most_selected_bck(fit,sel_model)
  
  # Note: the bck_feats_sel variable has the names of the variables. If we want to get splines, we need to update the model:
  # Mitotic rate
  if(grepl("rcs\\(Mitotic_rate,3\\)",fm_model)){
    if("Mitotic_rate"%in%bck_feats_sel[[2]]){
      bck_feats_sel[[2]][which(bck_feats_sel[[2]]=="Mitotic_rate")]="rcs(Mitotic_rate,3)"
    }
  }
  # Resection margin
  if(grepl("rcs\\(Resection_margin_cont_cm,3\\)",fm_model)){
    if("Resection_margin_cont_cm"%in%bck_feats_sel[[2]]){
      bck_feats_sel[[2]][which(bck_feats_sel[[2]]=="Resection_margin_cont_cm")]="rcs(Resection_margin_cont_cm,3)"
    }
  }
  if(sel_model =="clogit"){
    fit.without <- with(ds_mimp, expression(f1 = clogit(as.formula(paste0("Metastasis_numeric ~",paste(bck_feats_sel[[2]],collapse = "+"))),x=TRUE)))
  }else if(sel_model %in%c("uncox","uncox_recalib")){
    fit.without <- with(ds_mimp, expression(df  <-  data.frame(Metastasis_numeric,Vitfup_metastasis_years,Sex,Age,Tumor_location_cats,Number_of_cSCC_before_culprit,Tumor_diameter,logBreslowThickness,Differentiation,Immunocommpromised_at_time_of_cSCC,PNI_or_LVI,Tissue_involvement,Morphology_subgroup_bin),
                                            df <- truncateFUP(df,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5),
                                            Metastasis_numeric <- df[,"Metastasis_numeric"],
                                            Vitfup_metastasis_years <- df[,"Vitfup_metastasis_years"],
                                            Mitotic_Rate_ws <- ifelse(Mitotic_rate>2,2,Mitotic_rate),
                                            f1 = coxph(as.formula(paste0("Surv(Vitfup_metastasis_years,Metastasis_numeric) ~",paste(bck_feats_sel[[2]],collapse = "+"))),x=TRUE)))
  }else if (sel_model %in% c("wghcox","wghcoxnet_recalibvar")){
    fit.without <- with(ds_mimp, expression(df  <-  data.frame(Metastasis_numeric,Vitfup_metastasis_years,Sex,Age,Tumor_location_cats,Number_of_cSCC_before_culprit,Tumor_diameter,logBreslowThickness,Differentiation,Immunocommpromised_at_time_of_cSCC,PNI_or_LVI,Tissue_involvement,Morphology_subgroup_bin),
                                            df <- truncateFUP(df,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5),
                                            Metastasis_numeric <- df[,"Metastasis_numeric"],
                                            Vitfup_metastasis_years <- df[,"Vitfup_metastasis_years"],
                                            Mitotic_Rate_ws <- ifelse(Mitotic_rate>2,2,Mitotic_rate),
                                            f1 = coxph(as.formula(paste0("Surv(Vitfup_metastasis_years,Metastasis_numeric) ~",paste(bck_feats_sel[[2]],collapse = "+"))),x=TRUE,weights=weights)))
  }else if (sel_model %in%c("wghcoxnet_median","wghcoxnet_average","wghcoxnet_recalibscore")){
    fit.without  = fit # because we do not refit a model
  }
  return(fit.without)
}


build_mice_models = function(ds_mimp,model_method,fm_model,weights){
  if(model_method =="wghcox"){
    expr <- expression(df  <-  data.frame(Metastasis_numeric,Vitfup_metastasis_years,Sex,Age,Tumor_location_cats,Number_of_cSCC_before_culprit,Tumor_diameter,logBreslowThickness,Differentiation,Immunocommpromised_at_time_of_cSCC,PNI_or_LVI,Tissue_involvement,Morphology_subgroup_bin),
                       df <- truncateFUP(df,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5),
                       Metastasis_numeric <- df[,"Metastasis_numeric"],
                       Vitfup_metastasis_years <- df[,"Vitfup_metastasis_years"],
                       Mitotic_Rate_ws <- ifelse(Mitotic_rate>2,2,Mitotic_rate),
                       f1 <- cph(as.formula(paste0("Surv(Vitfup_metastasis_years,Metastasis_numeric) ~",fm_model)),weights=weights),
                       f2 <- fastbw(f1, rule="aic"))
  }else if (model_method %in%c("uncox","uncox_recalib")){
    expr <- expression(df  <-  data.frame(Metastasis_numeric,Vitfup_metastasis_years,Sex,Age,Tumor_location_cats,Number_of_cSCC_before_culprit,Tumor_diameter,logBreslowThickness,Differentiation,Immunocommpromised_at_time_of_cSCC,PNI_or_LVI,Tissue_involvement,Morphology_subgroup_bin),
                       df <- truncateFUP(df,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5),
                       Metastasis_numeric <- df[,"Metastasis_numeric"],
                       Vitfup_metastasis_years <- df[,"Vitfup_metastasis_years"],
                       Mitotic_Rate_ws <- ifelse(Mitotic_rate>2,2,Mitotic_rate),
                       f1 <- cph(as.formula(paste0("Surv(Vitfup_metastasis_years,Metastasis_numeric) ~",fm_model))),
                       f2 <- fastbw(f1, rule="aic"))
  }else if (model_method %in% c("wghcoxnet_median","wghcoxnet_average","wghcoxnet_recalibscore","wghcoxnet_recalibvar")){
    
    #I got stuck in how to access the column names within the expression environment without hard coding them, so I just coded the loop across the imputed datasets
    mi_models=list()
    for (i in c(1:ds_mimp$m)){
      ds = complete(ds_mimp,i)
      ds$Mitotic_Rate_ws <- ifelse(ds$Mitotic_rate>2,2,ds$Mitotic_rate)
      ds$Mitotic_rate<-NULL
      predictors = strsplit(gsub("\\+",",",gsub("\\s","",fm_model)), ",")[[1]]
      fit= train_wghcoxnet(ds,10,predictors,weights)
      #ds_x = ds[,strsplit(gsub("\\+",",",gsub("\\s","",fm_model)), ",")[[1]]]
      #x <- model.matrix(~.,data = ds_x)
      #ds=truncateFUP(ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5)
      #Metastasis_numeric <- ds[,"Metastasis_numeric"]
      #Vitfup_metastasis_years <- ds[,"Vitfup_metastasis_years"]
      #y <- Surv(Vitfup_metastasis_years,Metastasis_numeric)
      #fit <- cv.glmnet(x[,-1],y,family="cox",alpha=0.1,standardize=TRUE,weights = weights,seed=123)
      mi_models[[i]] = fit[["fit_train"]]
    }
  }
  
  if(!model_method %in%c("wghcoxnet_median","wghcoxnet_average","wghcoxnet_recalibscore","wghcoxnet_recalibvar")){
    mi_models = with(ds_mimp,expr)
  }
  
  return(mi_models)
} 

get_most_selected_bck = function(fit,model_name){
  if (model_name %in% c("wghcox","uncox","uncox_recalib","clogit")){
    votes=unlist(lapply(fit$analyses,function(x) x$names.kept))
    print(table(votes))
    # Keep only the features that were consistently chosen across datasets: 
    bck_feats = names(table(votes))
    consistently_chosen_bck_feats = names(table(votes))[which(table(votes)>5)] # improv: half of the imputed datasets was hard coded, because we always used 10 imputed datasets. if code is ever going to be more flexible, this should be changed
  }else if(model_name %in% c("wghcoxnet_median","wghcoxnet_average","wghcoxnet_recalibscore","wghcoxnet_recalibvar")){
    for (i in c(1:length(fit))){
      df_fit_cox_wgh_lasso=sapply(1:length(fit),function(x) as.matrix(coef(fit[[x]],s="lambda.1se")))
      
      bck_feats=NA
      #Obtain variables that were selected at least half of the time:
      var_selected=rownames(coef(fit[[1]],s="lambda.1se"))[rowSums(t(apply(df_fit_cox_wgh_lasso,1,function(x) x!=0)))>=5] # improv: half of the imputed datasets was hard coded, because we always used 10 imputed datasets.  if code is ever going to be more flexible, this should be changed
      consistently_chosen_bck_feats= gsub("man|Scalp/neck|Face|Poor/.*|Abundant|Yes|Extensive|R1/R2.*|Beyond.*|Subcutane.*|No/.*","",var_selected)
    }
  }
  return(list(bck_feats,consistently_chosen_bck_feats))
}
build_pooled_model = function(ds_mimp,model_method,fm_model,weights){
  mi_models = build_mice_models(ds_mimp,model_method,fm_model,weights)
  mi_models_selection = remove_least_freq_selected_vars(mi_models,ds_mimp,model_method,weights,fm_model)
  return(mi_models_selection)
}

pool_baseline_survival = function(ds_mimp,mice_models){
  #Pool baseline survival
  log_log_bsurvival = vector("numeric",length=ds_mimp$m)
  for(m in c(1:ds_mimp$m)){
    log_log_bsurvival[m]=log(-log(summary(survfit(mice_models$analyses[[m]]),times = 5)$surv)) # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-015-0048-4
    print(summary(survfit(mice_models$analyses[[m]]),times = 5)$surv)
  }
  baseline_survival_pooled = exp(-exp(mean(log_log_bsurvival)))
  
  return(baseline_survival_pooled)
  
}

# Note: This function was copied from the psfmi function pool_auc_adapted, without truncating the confidence intervals to the interval 0 to 1:
pool_auc_adapted <- function(est_auc, 
                             est_se, 
                             nimp = 5, 
                             log_auc=TRUE)
{
  
  RR_se <- function(est, se, nimp){
    m <- nimp
    w_auc <-
      mean(se^2) # within variance
    b_auc <-
      var(est) # between variance
    tv_auc <-
      w_auc + (1 + (1/m)) * b_auc # total variance
    se_total <-
      sqrt(tv_auc)
    r <- (1 + 1 / m) * (b_auc / w_auc)
    v <- (m - 1) * (1 + (1/r))^2
    t <- qt(0.975, v)
    res <- c(se_total, t)
    return(res)
  }
  
  est_auc <-
    unlist(est_auc)
  est_auc_se <-
    unlist(est_se)
  if(length(est_auc) != nimp)
    stop("Include c-statistic value for each imputed dataset")
  
  if(log_auc){
    est_auc_log <-
      log(est_auc/(1-est_auc))
    est_auc_se_log <-
      est_auc_se / (est_auc * (1-est_auc))
    
    se_total <-
      RR_se(est_auc_log, est_auc_se_log, nimp=nimp) # pooled se
    
    # Backtransform
    inv.auc <- exp(mean(est_auc_log)) /
      (1 + exp(mean(est_auc_log)))
    inv.auc.u <- exp(mean(est_auc_log) + (se_total[2]*se_total[1])) /
      (1 + exp(mean(est_auc_log) + (se_total[2]*se_total[1])))
    inv.auc.l <- exp(mean(est_auc_log) - (se_total[2]*se_total[1])) /
      (1 + exp(mean(est_auc_log) - (se_total[2]*se_total[1])))
    auc_res <- round(matrix(c(inv.auc.l, inv.auc, inv.auc.u),
                            1, 3, byrow = T), 4)
    dimnames(auc_res) <- list(c("C-statistic (logit)"),
                              c("95% Low", "C-statistic", "95% Up"))
  } else {
    mean_auc <-
      mean(est_auc)
    se_total <-
      RR_se(est_auc, est_auc_se, nimp=nimp)
    auc_u <-
      mean(est_auc) + (se_total[2]*se_total[1])
    #if(auc_u > 1) auc_u <- 1.00
    auc_l <- mean(est_auc) - (se_total[2]*se_total[1])
    #if(auc_l < 0) auc_l <- 0.00
    auc_res <-
      round(matrix(c(auc_l, mean_auc, auc_u),
                   1, 3, byrow = T), 4)
    dimnames(auc_res) <-
      list(c("C-statistic"), c("95% Low", "C-statistic", "95% Up"))
  }
  return(auc_res)
}

recalibrate_surv_probs = function(predictions,ds,weights,predictions_test=NULL,ds_test=NULL){
  predictions_test_recalib=NULL
  log_log_predictions_recalib_test = NULL
  log_log_predictions_recalib = as.data.frame(log(-log(predictions)))
  colnames(log_log_predictions_recalib) ="log_log_preds"
  log_log_predictions_recalib[,c("outcome_num","outcome_fup","Metastasis","Vitfup_metastasis_years","Metastasis_numeric")] = data.frame("outcome_num" =ds[,"Metastasis_numeric"],"outcome_fup"=ds[,"Vitfup_metastasis_years"],"Metastasis"=ds[,"Metastasis"],"Vitfup_metastasis_years"=ds[,"Vitfup_metastasis_years"],"Metastasis_numeric" =ds[,"Metastasis_numeric"])
  
  if(!is.null(predictions_test)){
    log_log_predictions_recalib_test = as.data.frame(log(-log(predictions_test)))
    colnames(log_log_predictions_recalib_test) ="log_log_preds"
    log_log_predictions_recalib_test[,c("outcome_num","outcome_fup","Metastasis","Vitfup_metastasis_years","Metastasis_numeric")] = data.frame("outcome_num" =ds_test[,"Metastasis_numeric"],"outcome_fup"=ds_test[,"Vitfup_metastasis_years"],"Metastasis"=ds_test[,"Metastasis"],"Vitfup_metastasis_years"=ds_test[,"Vitfup_metastasis_years"],"Metastasis_numeric" =ds_test[,"Metastasis_numeric"])
  }
  recalib.pred = recalibrate_probs_survival(predictions,ds[,"Vitfup_metastasis_years"],ds[,"Metastasis_numeric"],weights,5,predictions_test)
  predictions_recalib  = as.numeric(recalib.pred[["Predictions_train"]])
  fit_train_recalib = recalib.pred[["recalibrated_fit"]]
  predictions_test_recalib  = as.numeric(recalib.pred[["Predictions_tst"]])
  return(list("predictions_tr"=predictions_recalib,"fit"=fit_train_recalib,"predictions_test"=predictions_test_recalib,"log_log_predictions_tr"=log_log_predictions_recalib,"log_log_predictions_tst"=log_log_predictions_recalib_test))
}

combine_perf_metrics_across_imputations = function(list_perf_models,m){
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
      temp_df[metric,]=as.data.frame(pool_auc_adapted(list_perf_models[[model]][["Optimism-corrected of pooled model"]][,metric],list_perf_models[[model]][["Optimism se"]][,metric],m,log_auc = logit_tr))
    }
    colnames(temp_df)=c("lower_CI","Estimate","upper_CI")
    performance_with_ci[[model]]=temp_df
  }
  return(performance_with_ci)
}

train_wghcoxnet = function(ds,folds,predictors,weights,seed=123,cl_parallel=NULL,tm_pt=5){
  #Tidy up predictor matrix:
  ds_x = ds[,predictors]
  x <- model.matrix(~.,data = ds_x)
  ds=truncateFUP(ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",5)
  Metastasis_numeric <- ds[,"Metastasis_numeric"]
  Vitfup_metastasis_years <- ds[,"Vitfup_metastasis_years"]
  y <- Surv(Vitfup_metastasis_years,Metastasis_numeric)
  
  fit_train = optimize_wghcoxnet(x,y,weights,ds,folds,seed,cl_parallel)
  
  bs_y = summary(survfit(fit_train,x =x[,-1],y=y,weights=weights),times = tm_pt)$surv
  rm(.Random.seed, envir=.GlobalEnv) #This ensures that seed was only set inside this function
  return(list("fit_train"=fit_train,"baseline_survival"=bs_y,"best_alpha"=fit_train$best_alpha))
}

optimize_wghcoxnet =function(x,y,weights,ds,folds,seed,cl_parallel){
  chosen_vars_net = 0
  sd_extra=0
  alphas=c(0.1,0.5,1)
  while (chosen_vars_net ==0){
    #Set folds for reproducibility (will only be used by wghcoxnet and uncoxnet)
    set.seed(seed+sd_extra)
    #print(123+sd_extra)
    flds <- createFolds(ds[,"Metastasis"], k = folds, list = TRUE, returnTrain = FALSE)
    foldids = rep(1,length(ds[,"Metastasis"]))
    for (fdi in names(flds)[2:length(flds)]){
      foldids[flds[[fdi]]]=as.numeric(gsub("Fold","",fdi))
    }
    cv.fit = cva.glmnet(x[,-1],y,family="cox",alpha=alphas,standardize=TRUE,weights = weights,seed=seed,foldid = foldids,outerParallel=cl_parallel)
    
    #Find best alpha:
    # Compute c-indexes for each model
    i=1
    cind_aucs = c()
    for (aph in cv.fit$alpha){
      lp = predict(cv.fit$modlist[[i]],newx=data.matrix(x[,-1]),type="link",s="lambda.1se")
      cind_aucs[as.character(aph)] = 1-rcorr.cens(lp,Surv(ds[,"Vitfup_metastasis_years"],ds[,"Metastasis_numeric"]))[1] # rcorr.cens expects a survival numeric predictor (higher the less risk), hence we take the complement
      i=i+1
    }
    #Find best alpha:
    best_candidate = which.max(cind_aucs)
    lower_bound = cind_aucs[best_candidate]-0.02 # We only decrease alpha (and consequently increase predictors) if AUC is 2% higher.
    cind_aucs = cind_aucs[which(cind_aucs>lower_bound)]
    best_alpha = as.numeric(names(cind_aucs[length(cind_aucs)]))
    best_lambda = cv.fit$modlist[which(cv.fit$alpha==best_alpha)][[1]]$lambda.1se
    
    fit_train = glmnet(x[,-1],y,family="cox",alpha=best_alpha,standardize=TRUE,lambda = best_lambda,weights = weights,seed=seed)
    chosen_vars_net = length(which(coef(fit_train)>0.00000000001))
    sd_extra=sd_extra+1
    if(sd_extra==30){alphas=c(0.01,alphas)}
  }
  rm(.Random.seed, envir=.GlobalEnv) #This ensures that seed was only set inside this function
  fit_train$best_alpha = best_alpha
  return(fit_train)
}

get_prob_survival_model = function(cox_model_fit,tr_data,tst_data,times,baseline_survival=NULL,model_method="uncox",fm_model=NULL){
  #Note: this function was tested with the output of riskRegression's function: predictSurvProb, and it is the same for unweighted cox regression. This function does not work with strata!!
  #Get baseline survival
  if(is.null(baseline_survival)){
    summ_baseline_surv = summary(survfit(cox_model_fit),times = times)
    baseline_survival = summ_baseline_surv$surv
  }
  center_test = ifelse(model_method%in%c("wghcoxnet_median","wghcoxnet_mean","wghcoxnet"),FALSE,TRUE)
  lp = get_cox_linear_predictor(cox_model_fit,tr_data,tst_data,fm_model,center_test)
  surv_probs = matrix(NA,nrow=length(times),ncol =nrow(tst_data))
  for(tm in c(1:length(times))){
    surv_probs[tm,] = baseline_survival[tm]^exp(lp)#TO DO: check that baseline survival of coxnet is compatible with this formulation. At the moment, I'm assuming that baseline survival is compatible with the coefficients that are outputed by the model.
  }
  return(surv_probs)
  
}

get_cox_linear_predictor = function(cox_model_fit,tr_data,tst_data,fm_model,center_test = TRUE){
  # Compute Linear predictor
  #center_data = preProcess(tr_data,method = "center")
  #tst_data = predict(center_data,tst_data)
  numeric_columns_means = cox_model_fit$means[which(cox_model_fit$means!=0)] # For weighted cox regression, it is better to use the means of the variables from the model fit directly
  if(center_test){ #This was added because coefficients of coxnet are in the original units.
    tst_data[,names(numeric_columns_means)]= tst_data[,names(numeric_columns_means)]- matrix(numeric_columns_means, nrow=nrow(tst_data), ncol=length(numeric_columns_means), byrow=TRUE)
  }
  #lp = coxLP(cox_model_fit,tst_data,center=FALSE) #we don't use this function anymore because it does not work with coxnet
  if(is.null(fm_model)){
    if(ncol(tr_data)>1){
      lp = as.numeric(as.vector(t(as.matrix(coef(cox_model_fit))) %*% as.matrix(t(model.matrix.lm(~.,tst_data[,attr(cox_model_fit$terms,"term.labels")],na.action="na.pass")[,-1])))) #model.matrix.lm is similar to model.matrix, but allows NA values
    }else{
      if(is.numeric(tr_data[,1])){
        lp =  as.numeric(coef(cox_model_fit)*tst_data[,1])
      }
    }
  }else{
    vars = strsplit(gsub(" ","",gsub("\\+",",",gsub(".*~","",fm_model))),",")[[1]]
    coefs = cox_model_fit$coefficients
    x_matrix = model.matrix.lm(~.,rbind(tr_data,tst_data)[,vars],na.action="na.pass")
    x_matrix = x_matrix[(nrow(tr_data)+1):nrow(x_matrix),]
    if(all(rownames(coefs)==colnames(x_matrix)[-1])){x_matrix = x_matrix[,-1]}
    stopifnot(rownames(coefs)==colnames(x_matrix))
    lp = as.numeric(as.vector(t(as.matrix(coefs)) %*% as.matrix(t(x_matrix)))) #model.matrix.lm is similar to model.matrix, but allows NA values
    
  }
  return(lp)
}
