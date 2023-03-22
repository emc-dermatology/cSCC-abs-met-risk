#---------------------------------------------------
# Aim: Auxilliary functions for multiple imputation
# Author: B. Rentroia Pacheco
# Input: Preprocessed Dataset
# Output: Dataset for model imput
#---------------------------------------------------

#-----------------------
# 0. Setup
#-----------------------
dir_project = "" #fill with project directory
#-----------------------

# Winsorization of outliers 
# We have visualized and identified outliers. For each of the variables, we have checked that the outliers are at least biologically plausible and not likely data errors.
# For the remaining outliers, we will winsorize them.
# There are different ways of identifying the outliers. I chose to go with the following approach:
# Outliers are identified based on 3*MAD boundary
winsorize_outliers = function(x,type_var = "numeric"){
  
  if(type_var == "numeric"){
    #1. Identify the outliers according to:
    median_x = median(na.omit(x))
    mad_x = mad(na.omit(x))
    low_outliers = which(x<( median_x -3*mad_x))
    high_outliers = which(x>(median_x +3*mad_x))
    
    if (length(low_outliers)<0.01*length(x)){
      one_pct_percentile = quantile(na.omit(x),0.01)
    }else{
      one_pct_percentile = median_x -3*mad_x
    }
    
    if(length(high_outliers)<0.01*length(x)){
      nnine_pct_percentile = quantile(na.omit(x),0.99)
    }else{
      nnine_pct_percentile = median_x +3*mad_x
    }
    
  }else if (type_var == "count"){
    outs_poisson = aout.pois(data = x, param = median(x,na.rm=TRUE),alpha=0.01)
    high_outliers = which(outs_poisson$is.outlier==TRUE)
    high_outliers =setdiff(high_outliers,which(is.na(x)==TRUE))
    low_outliers = c()
    
    one_pct_percentile = 0
    nnine_pct_percentile = min(outs_poisson$data[high_outliers])
    
  }else if (type_var == "bimodal"){
    median_x = median(na.omit(x[which(x!=0)]))
    mad_x = mad(na.omit(x[which(x!=0)]))
    low_outliers = c()
    high_outliers = which(x>(median_x +3*mad_x))
    
    if(length(high_outliers)<0.01*length(x)){
      nnine_pct_percentile = quantile(na.omit(x),0.99)
    }else{
      nnine_pct_percentile = median_x +3*mad_x
    }
    
  }
  
  
  
  #If variable only contains integers, then we should replace the outliers with integer values.
  if(all(na.omit(x) == floor(na.omit(x)))){
    one_pct_percentile = floor(one_pct_percentile)
    nnine_pct_percentile = floor(nnine_pct_percentile)
  }
  
  #2. Winsorize outliers
  
  outliers = c()
  new_values = c()
  
  if(length(low_outliers)>0){
    if(all(one_pct_percentile>=x[low_outliers])){
      outliers = c(outliers,x[low_outliers])
      new_values= c(new_values,rep(one_pct_percentile,length(low_outliers)))
      
      x[low_outliers] = one_pct_percentile
      
    }else{
      print("1st percentile is lower than the outliers")
    }
  }
  
  if(length(high_outliers)>0){
    if(all(nnine_pct_percentile<=x[high_outliers])){
      outliers = c(outliers,x[high_outliers])
      new_values= c(new_values,rep(nnine_pct_percentile,length(high_outliers)))
      
      x[high_outliers] = nnine_pct_percentile
    }else{
      print("99th percentile is higher than the outliers")
    }
  }
  
  if(length(outliers)>0){
    df = as.data.frame(matrix(c(outliers,new_values),ncol=2))
  }else{
    df=as.data.frame(matrix(0,ncol=2,nrow=0))
  }
  return(list(x,df))
}


# Function that performs multiple imputation on a dataset.
#exc_vars_ds = c("Tumor_diameter_factorized","Ratio_width_diameter","Recurrence_cSCC","PNI","Lymphovascular_invasion","Sebaceous_gland_carcinoma_possibility")
#exc_vars_imp = c("Resection_margin_cont","Morphology_subgroup_bin","Morphology_present_less30pct","Invasion_sum","Metastasis_numeric","Set_id")




performMultipleImputation = function(ds_raw,exc_vars_ds,exc_vars_imp,winsorization,order_vars,nr_imp_ds,cor_pred_mat,niter,seed,dir_folder_results,analysis_id,save=TRUE){
  
  if (save){
    #0. Set folder:
    dir_subfolder_results = paste0(dir_folder_results,"/",analysis_id,"_ws_",winsorization,"_ovr_",order_vars,"_nimp_",nr_imp_ds,"_cpm_",cor_pred_mat)
    if(!dir.exists(dir_subfolder_results)){dir.create(dir_subfolder_results)}
  }
  # 1. Remove any variables that are not meant to go to imputation:
  ds = ds_raw[,setdiff(colnames(ds_raw),c(exc_vars_ds))]
  
  num_vars = colnames(ds)[unlist(lapply(ds, is.numeric))]
  cat_vars  = colnames(ds)[unlist(lapply(ds, is.factor))]
  
  # 2. Winsorization
  # Diagnostics of winsorizing process: Distribution before and after
  if(save){
    pdf(paste0(dir_subfolder_results,"/RSCM3_1 Winsorizing process before and after.pdf"),width=16,height=16)
    par(mfrow=c(2,3))
  }
  summary_thresholds_outliers = as.data.frame(matrix(0,ncol=3,nrow=0))
  thresholds_man = c(4,5,3,10,15,5,5,1,5)
  names(thresholds_man) = c("Tumor_diameter","Number_of_cSCC_before_culprit","Number_of_cSCC_after_culprit","Depth_of_Invasion","Breslow_thickness","Tumor_budding","Mitotic_rate","Resection_margin_cont_cm","Tumor_width")
  
  for (n_var in num_vars){
    #print(n_var)
    x = ds[,n_var]
    if(save){plot(density(na.omit(x)),main=n_var)}
    
    if(winsorization =="autom"){
      if(n_var %in%c("Number_of_cSCC_before_culprit","Number_of_cSCC_after_culprit","Tumor_budding","Mitotic_rate")){
        new_values = winsorize_outliers(x,"count")
      }else if (n_var %in% c("PNI_diameter","Resection_margin_cont_cm")){
        new_values = winsorize_outliers(x,"bimodal")
      }else{
        new_values = winsorize_outliers(x,"numeric")
      }
      ds[,n_var] = new_values[[1]]
      if(nrow(new_values[[2]])>1){
        thresholds_df = as.data.frame(new_values[[2]])
        thresholds_df$V3 = n_var
      }else{
        thresholds_df = as.data.frame(matrix(c(NA,NA,n_var),ncol=3,nrow=1))
      }
      summary_thresholds_outliers = rbind(summary_thresholds_outliers,thresholds_df )
      thresholds = unique(thresholds_df [,"V2"])
    }else if (winsorization =="manual"){
      thresholds = thresholds_man[n_var]
      i.outliers = which(ds[,n_var]>thresholds)
      ds[i.outliers,n_var]=thresholds
      
      
    }
    if(save){
      abline(v=thresholds,lty="dashed",col="red")
      lines(density(na.omit(ds[,n_var])),col="blue")}
  }
  
  if(save){dev.off()}
  
  if(save){
    if(winsorization =="autom"){
      write.csv(summary_thresholds_outliers,paste0(dir_subfolder_results,"/RSCM3 MAD thresholds.csv"))
    }
  }
  
  # 3. Multiple Imputation using MICE 
  
  #If any variable takes 1 value and has missing values, the remaining have to be imputed with that value:
  # This was added to avoid mice returning matrices with missing values after imputation
  i.impute_directly = which((apply(droplevels(ds),2,function (x)length(table(x)))==1 & apply(ds,2,function(x)sum(is.na(x))>0))==TRUE)
  if(length(i.impute_directly)>0){
    for (col in names(i.impute_directly)){
      print(paste0(col," was imputed directly"))
      ds[is.na(ds[,col]),col]=names(table(droplevels(ds)[,col]))
    }
  }
  # Keeping unused levels can leave to problems during mice
  ds=droplevels(ds)
  
  # Variables with near zero variance will not be used as predictors in imputation models, to avoid convergence issues during mice
  near_zero_vars=setdiff(colnames(ds)[nearZeroVar(ds)],c("PNI_bin","Lymphovascular_invasion_bin","Invasion_beyond_subcutaneous_fat","Invasion_of_muscle"))
  print(paste0("Variables", paste0(near_zero_vars,collapse=", ")," will not be used as predictors in imputation models."))
  exc_vars_imp = c(exc_vars_imp,near_zero_vars)
  
  # 3.1 Which predictors will we use in the imputation of each variable?
  #This can depend on the winsorization method
  pred_mt_mimp = quickpred(ds,include = c("Age","Sex","Metastasis"),exclude =exc_vars_imp, mincor=cor_pred_mat) # Initially we kept these variables to guarantee we would still have predictors if we had set a high threshold on the correlation, but later we decided to include all predictors.
  
  if(save){
    #Save plot
    pdf(paste0(dir_subfolder_results,"/RSCM3_1 predictors_imputation.pdf"),width=16,height=16)
    heatmap(pred_mt_mimp,col = c("white","darkblue"),margins=c(15,15),scale="none")
    dev.off()
  }
  #This matrix should be checked.
  
  # 3.2 Passive imputation and post processing:
  #Some variables are coded through passive imputation
  #Initial imputation. This is just to get access to the structure of inputs of mice.
  ini <- mice(ds, maxit=0, print=F,predictorMatrix = pred_mt_mimp,seed=123)
  meth<- ini$meth
  pred_mt_imp <- ini$pred
  post <- ini$post
  
  if ("Invasion_sum"%in%colnames(ds)){
    if("Invasion_of_other_bones_deep_structures" %in% colnames(ds)){
      pred_mt_imp[c("Invasion_of_muscle","Invasion_of_bones","Invasion_of_other_bones_deep_structures","Invasion_cartilage","Clark_level"),"Invasion_sum"]=0
      meth["Invasion_sum"]<- "~ I(as.numeric(Invasion_beyond_subcutaneous_fat=='Yes') + as.numeric(Invasion_of_muscle=='Yes') + as.numeric(Invasion_of_bones=='Yes')+ as.numeric(Invasion_of_other_bones_deep_structures=='Yes')+ as.numeric(Invasion_cartilage=='Yes'))"
    }else if ("Invasion_cartilage" %in% colnames(ds)){
      if(!all(ds[!is.na(ds$Invasion_cartilage),"Invasion_cartilage"]=="No")){
      pred_mt_imp[c("Invasion_of_muscle","Invasion_of_bones","Invasion_cartilage","Clark_level"),"Invasion_sum"]=0
      meth["Invasion_sum"]<- "~ I(as.numeric(Invasion_beyond_subcutaneous_fat=='Yes') + as.numeric(Invasion_of_muscle=='Yes') + as.numeric(Invasion_of_bones=='Yes'))"
      }else if(!all(ds[!is.na(ds$Invasion_of_bones),"Invasion_of_bones"]=="No")){
      pred_mt_imp[c("Invasion_of_muscle","Invasion_of_bones","Clark_level"),"Invasion_sum"]=0
      meth["Invasion_sum"]<- "~ I(as.numeric(Invasion_beyond_subcutaneous_fat=='Yes') + as.numeric(Invasion_of_muscle=='Yes') + as.numeric(Invasion_of_bones=='Yes'))"
      }else{
        pred_mt_imp[c("Invasion_of_muscle","Clark_level"),"Invasion_sum"]=0
        meth["Invasion_sum"]<- "~ I(as.numeric(Invasion_beyond_subcutaneous_fat=='Yes') + as.numeric(Invasion_of_muscle=='Yes'))"
        
      }
    }
  }
  
  if ("PNI_or_LVI"%in%colnames(ds)){
    pred_mt_imp[c("PNI_bin","Lymphovascular_invasion_bin"),"PNI_or_LVI"]=0
    meth["PNI_or_LVI"]<- "~ I(ifelse(PNI_bin=='Yes' |Lymphovascular_invasion_bin=='Yes', 'Yes','No'))"
  }
  
  if ("PNI_diameter"%in%colnames(ds)){
    post["PNI_diameter"] = "imp[[j]][, i] <- ifelse(data$PNI_bin[!r[, j]]== 'No',0,imp[[j]][, i])"
  }
  if("Resection_margin_cont_cm"%in% colnames(ds)){
    post["Resection_margin_cont_cm"] = "imp[[j]][, i] <- ifelse(data$Resection_margin_cat[!r[, j]] == 'R1/R2 (irradicaal)',0,imp[[j]][, i])"
  }
  
  # 3.3 Multiple imputation using MICE:
  #The order of the multiple imputation has to be decided (done with the vis parameter)
  dutch_ds_mimp=mice(ds,print = FALSE, maxit = niter, m = nr_imp_ds, seed = seed,vis = order_vars,predictorMatrix = pred_mt_mimp,meth = meth,post=post)
  
  if(save){
    # 3.4 Diagnostics:
    # The first thing to check is whether the imputation has converged over the iterations
    pdf(paste0(dir_subfolder_results,"/RSCM3_1 Mice convergence plots.pdf"),width=9,height=12)
    print(plot(dutch_ds_mimp))
    dev.off()
    
    # The second is to look at the distribution of the imputed and the original data (for continuous variables)
    pdf(paste0(dir_subfolder_results,"/RSCM3_1 Mice density plots after imputation.pdf"),width=9,height=6)
    print(densityplot(dutch_ds_mimp, layout = c(3, 2),lwd=2))
    dev.off()
    pdf(paste0(dir_subfolder_results,"/RSCM3_1 Mice boxplots after imputation.pdf"),width=9,height=6)
    print(mice::bwplot(dutch_ds_mimp, layout = c(3, 2)))
    dev.off()
    #stripplot(dutch_ds_mimp, layout = c(3, 2))
    
    # For categorical variables:
    plots_cat_vars = list()
    for (c_var in setdiff(cat_vars,"Set_id")){
      df_plot = data.frame(ds[,c_var],rep("Original",nrow(ds)))
      colnames(df_plot) = c(c_var,"dataset")
      if (sum(is.na(df_plot[,c_var]))>0){
        df_plot  = df_plot[!is.na(df_plot[,c_var]),]
      }
      for (i in c(1:nr_imp_ds)){
        temp_df = data.frame(complete(dutch_ds_mimp,i)[,c_var],rep(i,nrow(ds)))
        colnames(temp_df ) = c(c_var,"dataset")
        df_plot = rbind(df_plot,temp_df)
      }
      
      
      #ggplot
      plots_cat_vars[[c_var]]=df_plot %>% group_by(across(all_of(c_var)),dataset) %>% dplyr::count() %>% ungroup() %>% group_by(dataset) %>% mutate(percent = n/sum(n)) %>% ggplot(aes_string(x=c_var,fill="dataset",y="percent"))+geom_col(position="dodge",col="black")+theme_bw()+scale_fill_manual(values=c("#8A0000","#8B0000","#8C0000","#8D0000","#8E0000", "#900000", "#910000", "#920000", "#930000", "#940000","#00008B"))+ylab("Percentage")+xlab(c_var)+scale_x_discrete(guide = guide_axis(n.dodge = 2))
    }
    
    g=marrangeGrob(grobs=plots_cat_vars,nrow=3,ncol=3) 
    ggsave(paste0(dir_subfolder_results,"/RSCM3_1 bps missing imputation cat variables.pdf"),width=16.5,height=9.1,g)

    # 3.5 Check logged events as well:
    write.csv(dutch_ds_mimp$loggedEvents,paste0(dir_subfolder_results,"/RSCM3_1 MICE logged events.csv"))
    
    
    # 3.6 Sanity check: compare Tumor diameter factorized with Tumor diameter
    pdf (paste0(dir_subfolder_results,"/RSCM3_1 Tumor diam imput vs Selin factor.pdf"),width=16,height=12)
    par(mfrow=c(2,3))
    for(i in c(1:nr_imp_ds)){
      plot(complete(dutch_ds_mimp,i)[,"Tumor_diameter"][is.na(ds$Tumor_diameter)],ds_raw$Tumor_diameter_factorized[is.na(ds_raw$Tumor_diameter)],xlab = paste0("Tumor diameter in imputed dataset ",i),ylab= "Tumor diameter factorized") 
      cor_diameters= cor(complete(dutch_ds_mimp,i)[,"Tumor_diameter"][is.na(ds$Tumor_diameter)],ds_raw$Tumor_diameter_factorized[is.na(ds_raw$Tumor_diameter)],use="pairwise.complete.obs",method="spearman") 
      text(1,3,paste0("Separman's cor\n",round(cor_diameters,2)))
    }
    dev.off()
    
    # 4. Save objects:
    #-----------------------
    save(dutch_ds_mimp,file=paste0(dir_subfolder_results,"/RSCM3_1 MICE object.RData"))
    sink(file = paste0(dir_subfolder_results,"/RSCM3_1 MICE Specs.txt"))
    print("Imputation methods: ")
    print(dutch_ds_mimp$meth)
    paste0("Postprocessing: ")
    print(dutch_ds_mimp$post)
    paste0("Number of iterations: ")
    print(dutch_ds_mimp$iteration)
    closeAllConnections()
    #-----------------------
  }
  return(dutch_ds_mimp)
}