#---------------------------------------------------
# Aim: Show decision curve analysis for the developed model
# Author: B. Rentroia Pacheco
# Input: Recalibrated Model
# Output: Decision curve analysis
#---------------------------------------------------


#-----------------------
# 0. Setup
#-----------------------
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_01_19 RSCM_A 6 Weighted Performance Metrics functions.R"))
source(paste0(dir_scripts,"/Model Building and Evaluation/2022_02_17 RSCM_8 A1 Internal Validation MI_Boot.R"))
#-----------------------

#-----------------------
# 1. Get data:
#-----------------------
computed_probs = read.csv(paste0(dir_int_output,"RSCM 9 weighted cox model probabilities on 390 dataset_updated_centered_LP.csv"))
load(paste0(dir_int_output,"RSCM1_1 Dutch Dataset coding_v3_0.RData"))
psamp = read.csv(paste0(dir_int_output,"RSCM 5_3 Sampling weights VitFUP and PA.csv"))
psamp_rordered = psamp[match(dutch_ds$Administratie_nr,psamp$pat.ids),]
stopifnot(length(which(psamp_rordered$Metastasis!=dutch_ds$Metastasis_numeric))==0)
weights = psamp_rordered$weights*12203/sum(psamp_rordered$weights) #Weight rescaling. We use the sum of potential controls + cases = 12203 in the full cohort

#Multiple imputed datasets as well:
folder_name = paste0(dir_int_results,"MICE/v3_0 coding/_ws_","manual","_ovr_","monotone","_nimp_10_cpm_","1e-04")
load(paste0(folder_name,"/RSCM3_1 MICE object.RData")) # multiple imputed dutch dataset (dutch_ds), can be fecthed with dutch_ds_mimp object.
load(paste0(dir_int_output,"RSCM1_1 Dutch Dataset coding_v3_0.RData"))

#BWH:
# We compute the BWH across all imputed datasets:
# We used the supplementary information of the Venables 2022 publication to compute this staging: https://pubmed.ncbi.nlm.nih.gov/34862598/
# And average across all:
bwh_staging = list()
for (i in c(1:10)){
  ds = complete(dutch_ds_mimp,i)
  ds=ds %>% dplyr::mutate(HR_Tumordiam = ifelse(Tumor_diameter>=2,1,0), 
                               HR_Diff = ifelse( Differentiation=="Poor/undifferentiated",1,0),
                               HR_PNI = ifelse(PNI_bin=="Yes",1,0),
                               HR_Tinv = ifelse(Tissue_involvement=="Beyond subcutaneous fat",1,0),
                               Sum_HR = HR_Tumordiam+HR_Diff+HR_PNI+HR_Tinv,
                               Sum_HR = ifelse(!(is.na(Invasion_of_bones))&Invasion_of_bones=="Yes",4,Sum_HR),
                               Decision_BWH = ifelse(Sum_HR>1,1,0)) #Division at T2a
                                
  ds$BWH_staging = ifelse(ds$Sum_HR==0,"T1",NA)
  ds$BWH_staging = ifelse(ds$Sum_HR==1,"T2a",ds$BWH_staging)
  ds$BWH_staging = ifelse(ds$Sum_HR%in%c(2,3),"T2b",ds$BWH_staging)
  ds$BWH_staging = ifelse(ds$Sum_HR>=4,"T3",ds$BWH_staging)
  
  bwh_staging[[i]]=ds$BWH_staging
}

bwh_staging_df = t(as.data.frame(matrix(unlist(bwh_staging),nrow=length(bwh_staging),byrow=TRUE)))
#readapted from : https://stackoverflow.com/questions/28094468/how-to-calculate-the-mode-of-all-rows-or-columns-from-a-dataframe
modefunc <- function(x){
  set.seed(123) # In case we have the half of the iterations have the same category and the other have another
  tabresult <- table(x)
  themode <- names(tabresult)[which(tabresult == max(tabresult))]
  if(sum(tabresult == max(tabresult))>1) themode <- sample(themode,1)
  return(themode)
}
bwh_staging_avg= apply(bwh_staging_df,1,modefunc)


dutch_ds$BWH_staging = bwh_staging_avg
dutch_ds=dutch_ds %>% dplyr::mutate(Decision_BWH = ifelse(bwh_staging_avg%in%c("T1","T2a"),0,1),
                    BWH_risk_estimates = as.numeric(revalue(as.character(bwh_staging_avg),c("T1"=0.007,"T2a"=0.045,"T2b"=0.367,"T3"=0.5)))) # As reported in: Jambusaria-Pahlajani 2013: https://pubmed.ncbi.nlm.nih.gov/23325457/

# Add the information on BWH staging in the patient characteristics table
dutch_ds$BWH_NA = ifelse(is.na(dutch_ds$Differentiation)|is.na(dutch_ds$Invasion_of_bones)|is.na(dutch_ds$Tumor_diameter)|is.na(dutch_ds$Tissue_involvement)|is.na(dutch_ds$PNI_bin),1,0)
table(dutch_ds$Metastasis,dutch_ds$BWH_staging,dutch_ds$BWH_NA)
dutch_ds_test = dutch_ds
dutch_ds_test$BWH_staging = ifelse(dutch_ds_test$BWH_NA==1,NA,dutch_ds_test$BWH_staging)
dutch_ds_test$BWH_staging = factor(dutch_ds_test$BWH_staging,levels=c("T1","T2a","T2b","T3"))

# AJCC:
# There is an erratum of the AJCC8, according to Selin Tokez. It's about >=2 or >2. Everyone uses >2, so we will use that one. Results shouldnt change.
# We used the supplementary information of the Venables 2022 publication to compute this staging: https://pubmed.ncbi.nlm.nih.gov/34862598/
ajcc_staging = list()
for (i in c(1:10)){
  ds = complete(dutch_ds_mimp,i)
  ds=ds %>% mutate(T1 = ifelse(Tumor_diameter<2 &((PNI_bin=="No"|(PNI_bin=="Yes"&Tissue_involvement=="Dermis"))&Invasion_beyond_subcutaneous_fat=="No"&Depth_of_Invasion<=6),1,0),
                   T2 = ifelse(Tumor_diameter>=2 & Tumor_diameter<4 & ((PNI_bin=="No"|(PNI_bin=="Yes"&Tissue_involvement=="Dermis"))&Invasion_beyond_subcutaneous_fat=="No"&Depth_of_Invasion<=6),1,0),
                   T3 = ifelse((Tumor_diameter>=4|(PNI_bin=="Yes"&Tissue_involvement!="Dermis")|Invasion_beyond_subcutaneous_fat=="Yes"|Depth_of_Invasion>6)&Invasion_of_bones=="No",1,0),
                   T4 = ifelse(Invasion_of_bones=="Yes",1,0),
                   Decision_AJCC= ifelse(T1==1,0,1))#Division at T1/T2
  
  ds$AJCC_staging = ifelse(ds$T1==1,"T1",NA)
  ds$AJCC_staging = ifelse(ds$T2==1,"T2",ds$AJCC_staging)
  ds$AJCC_staging = ifelse(ds$T3==1,"T3",ds$AJCC_staging)
  ds$AJCC_staging = ifelse(ds$T4==1,"T4",ds$AJCC_staging)
  
  ajcc_staging[[i]]=ds$AJCC_staging
}
ajcc_staging_df = t(as.data.frame(matrix(unlist(ajcc_staging),nrow=length(ajcc_staging),byrow=TRUE)))
ajcc_staging_avg= apply(ajcc_staging_df,1,modefunc)
dutch_ds$AJCC_staging = ajcc_staging_avg
dutch_ds=dutch_ds %>% mutate(Decision_AJCC = as.numeric(ifelse(ajcc_staging_avg=="T1",0,1)))


dutch_ds_test$AJCC_staging =dutch_ds$AJCC_staging
dutch_ds_test$AJCC_NA = ifelse(is.na(dutch_ds$Tumor_diameter)|is.na(dutch_ds$Tissue_involvement)|is.na(dutch_ds$PNI_bin)|is.na(dutch_ds$Depth_of_Invasion)|is.na(dutch_ds$Invasion_of_bones),1,0)
dutch_ds_test$AJCC_staging = ifelse(dutch_ds_test$AJCC_NA==1,NA,dutch_ds_test$AJCC_staging)
dutch_ds_test$AJCC_staging = factor(dutch_ds_test$AJCC_staging,levels=c("T1","T2","T3","T4"))
dutch_ds$weights=weights

df_surv = merge(computed_probs,dutch_ds[,c("PALGA_excerpt_id","Metastasis_numeric","Vitfup_metastasis_years","Decision_BWH","BWH_risk_estimates","AJCC_staging","Decision_AJCC","BWH_staging","weights")],by.x="PALGA_id",by.y="PALGA_excerpt_id")
df_surv$Met_risk = 1-df_surv$SurvProb
df_surv$Met_risk_th = ifelse(df_surv$Met_risk>0.05,1,0)
#-----------------------

#-----------------------
# 1. Summary statistics 
#-----------------------
i.completecase = is.na(dutch_ds$Tumor_diameter)|is.na(dutch_ds$PNI_bin)|is.na(dutch_ds$Invasion_beyond_subcutaneous_fat)|is.na(dutch_ds$Depth_of_Invasion)|is.na(dutch_ds$Invasion_of_bones)|is.na(dutch_ds$Differentiation)
pdf(paste0(dir_int_results,"/RSCM 11 distribution AJCC and BWH imputed.pdf"),width=12,height=10)
par(mfrow=c(2,2))
  barplot(table(dutch_ds$AJCC_staging),xlab="AJCC staging",ylab="No.patients",main="with imputed variables")
  barplot(table(dutch_ds$BWH_staging),xlab="BWH staging",ylab = "No.patients",main="with imputed variables")
  barplot(table(dutch_ds$AJCC_staging[!i.completecase]),xlab="AJCC staging",ylab="No.patients",main="complete case")
  barplot(table(dutch_ds$BWH_staging[!i.completecase]),xlab="BWH staging",ylab = "No.patients",main="complete case")
dev.off()

p1=ggplot(dutch_ds,aes(x=AJCC_staging,fill=Metastasis))+geom_bar()+theme_bw()+xlab("AJCC staging")
p2=ggplot(dutch_ds,aes(x=BWH_staging,fill=Metastasis))+geom_bar()+theme_bw()+xlab("BWH staging")
p3=ggplot(dutch_ds[!i.completecase,],aes(x=AJCC_staging,fill=Metastasis))+geom_bar()+theme_bw()+xlab("AJCC staging")
p4=ggplot(dutch_ds[!i.completecase,],aes(x=BWH_staging,fill=Metastasis))+geom_bar()+theme_bw()+xlab("BWH staging")
p_list=list(p1,p2,p3,p4)
grid_ps = grid.arrange(arrangeGrob(grobs=p_list))
ggsave(paste0(dir_int_output,"/RSCM 11 case control distribution staging.pdf"),width=12,height=10,grid_ps)

table(dutch_ds$BWH_staging,dutch_ds$AJCC_staging)


#Add BWH and AJCC staging to the patient characteristics table:
table_stg=dutch_ds_test[,c("AJCC_staging","BWH_staging","Metastasis","Set_id")] %>%
  sjlabelled::remove_all_labels() %>%
  tibble %>%
  tbl_summary(
    by = Metastasis, type = list(where(is.integer) ~ "continuous",where(is.numeric) ~ "continuous"),include = -Set_id) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p(list(all_continuous() ~ "paired.wilcox.test",all_categorical()~"mcnemar.test"),group=Set_id) %>% # test for a difference between groups
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels()

table_stg %>%
  as_flex_table() %>%
  flextable::save_as_docx(path=paste0(dir_results,script_id," Patient_characteristics stratified 390pts Staging systems.docx"))

# Add p-values:
# Merge categories for computation of p-values:
dutch_ds_test$BWH_staging = revalue(dutch_ds_test$BWH_staging,c("T3"="T2b"))
dutch_ds_test$AJCC_staging = revalue(dutch_ds_test$AJCC_staging,c("T4"="T3"))
for(var in c("AJCC_staging","BWH_staging")){
  tab_var = merge(dutch_ds_test[which(dutch_ds_test$Metastasis=="Case"),c("Set_id",var)],dutch_ds_test[which(dutch_ds$Metastasis=="Control"),c("Set_id",var)],by="Set_id")
  tab_var = tab_var[complete.cases(tab_var),]
  mcnemar.test(tab_var[,2],tab_var[,3])
}

# Check probabilities within risk groups
df_surv$Metastasis = factor(df_surv$Metastasis,levels=c("Control","Case"))
#BWH:
plot_boxplot_bwh = ggplot(df_surv,aes(y=Met_risk*100,x=BWH_staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("BWH stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/RSCM 11 estimated prob by BWH staging.pdf"),plot_boxplot_bwh ,width=12,height=7)
plot_boxplot_bwh_cc = ggplot(df_surv[df_surv$PALGA_id%in%dutch_ds$PALGA_excerpt_id[!i.completecase],],aes(y=Met_risk*100,x=BWH_staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("BWH stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/RSCM 11 estimated prob by BWH staging CC.pdf"),plot_boxplot_bwh_cc ,width=12,height=7)

#AJCC:
plot_boxplot_ajcc = ggplot(df_surv,aes(y=Met_risk*100,x=AJCC_staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("AJCC stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/RSCM 11 estimated prob by AJCC staging.pdf"),plot_boxplot_ajcc ,width=12,height=7)
plot_boxplot_ajcc_cc = ggplot(df_surv[df_surv$PALGA_id%in%dutch_ds$PALGA_excerpt_id[!i.completecase],],aes(y=Met_risk*100,x=AJCC_staging,fill= Metastasis))+geom_boxplot(outlier.shape = NA)+theme_bw()+geom_point(aes(fill =Metastasis), shape = 21,position = position_jitterdodge()) +ylab("Predicted metastatic risk (%)")+xlab("AJCC stage")+scale_y_log10()+scale_fill_manual(values=c("#00BFC4","#F8766D"))
ggsave(paste0(dir_results,"/RSCM 11 estimated prob by AJCC staging CC.pdf"),plot_boxplot_ajcc_cc ,width=12,height=7)

#-----------------------
# 2. Discrimination:  
#-----------------------
#With and without imputed datasets
discriminative_performances = as.data.frame(matrix(NA,ncol=4,nrow=4))
colnames(discriminative_performances)= c("AUC","Cind","wAUC","wCind")
rownames(discriminative_performances)= c("AJCC","BWH","AJCC_cc","BWH_cc")

for( imp in c("","_cc")){
  if(imp ==""){
    ds = dutch_ds
    weights_ds = weights
  }else if(imp =="_cc"){
    ds = dutch_ds[!i.completecase,]
    weights_ds = weights[!i.completecase]
  }
  
  for (staging_system in c("AJCC","BWH")){
    if(staging_system == "AJCC"){
      staging_ranks =  as.numeric(revalue(ds$AJCC_staging,c("T1"=1,"T2"=2,"T3"=3,"T4"=4)))
    }else if (staging_system =="BWH"){
      staging_ranks = as.numeric(revalue(ds$BWH_staging,c("T1"=1,"T2a"=2,"T2b"=3,"T3"=4)))
    }
    
    # unweighted:
    discriminative_performances[paste0(staging_system,imp),"Cind"]=round(1-rcorr.cens(staging_ranks,Surv(ds[,"Vitfup_metastasis_years"],ds[,"Metastasis_numeric"]))[1],2)
    auc_estimate = pROC::roc(ds[,"Metastasis_numeric"], staging_ranks , quiet = TRUE,direction = "<",ci=TRUE)
    discriminative_performances[paste0(staging_system,imp),"AUC"]= paste0(round(auc_estimate$auc,2),"(",round(auc_estimate$ci[1],2),"-",round(auc_estimate$ci[3],2),")")
  
    # weighted:
    discriminative_performances[paste0(staging_system,imp),"wAUC"] = round(WeightedAUC(WeightedROC(staging_ranks, ds[,"Metastasis"],weight=weights_ds)),2)
    discriminative_performances[paste0(staging_system,imp),"wCind"]=round(intsurv::cIndex(ds[,"Vitfup_metastasis_years"],event=ds[,"Metastasis_numeric"],staging_ranks,weight=weights_ds)[[1]],2)
    }
}
#write.csv(discriminative_performances,paste0(dir_folder_results,"/RSCM 11 AUC_Cind_staging_systems.csv"))

all_performances_mimp_multipleTP = as.data.frame(matrix(NA,nrow=0,ncol = 4))

for (tm in c(1,3,5)){
    ds = dutch_ds
    ds = truncateFUP(ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",tm)
    discriminative_performances_mimp_mean= as.data.frame(matrix(NA,ncol=10,nrow=9))
    rownames(discriminative_performances_mimp_mean)= c("AJCCCind","AJCCAUC","AJCCwCind","AJCCwAUC","BWHCind","BWHAUC","BWHwCind","BWHwAUC","BWHwCslope")
    discriminative_performances_mimp_se= as.data.frame(matrix(NA,ncol=10,nrow=9))
    rownames(discriminative_performances_mimp_se)= c("AJCCCind","AJCCAUC","AJCCwCind","AJCCwAUC","BWHCind","BWHAUC","BWHwCind","BWHwAUC","BWHwCslope")
    
    for (i in c(1:10)){
      for (staging_system in c("AJCC","BWH")){
        if(staging_system == "AJCC"){
          staging_ranks =  as.numeric(revalue(ajcc_staging_df[,i],c("T1"=1,"T2"=2,"T3"=3,"T4"=4)))
        }else if (staging_system =="BWH"){
          staging_ranks = as.numeric(revalue(bwh_staging_df[,i],c("T1"=1,"T2a"=2,"T2b"=3,"T3"=4)))
          risk_probs_BWH = as.numeric(revalue(bwh_staging_df[,i],c("T1"=0.007,"T2a"=0.045,"T2b"=0.367,"T3"=0.5)))
          log_log_prob_BWH = log(-log(1-risk_probs_BWH))
        }
        
        # unweighted:
        cind = rcorr.cens(staging_ranks,Surv(ds[,"Vitfup_metastasis_years"],ds[,"Metastasis_numeric"]))
        discriminative_performances_mimp_mean[paste0(staging_system,"Cind"),i]=1- cind [1]
        discriminative_performances_mimp_se[paste0(staging_system,"Cind"),i]=cind['S.D.']/2 # https://stat.ethz.ch/pipermail/r-help/2011-May/278817.html
        auc_estimate = pROC::roc(ds[,"Metastasis_numeric"], staging_ranks , quiet = TRUE,direction = "<",ci=TRUE)
        discriminative_performances_mimp_mean[paste0(staging_system,"AUC"),i] = auc_estimate$auc
        discriminative_performances_mimp_se[paste0(staging_system,"AUC"),i] = sqrt(pROC::var(auc_estimate))
        
        # weighted: for this we have to perform 100 bootstrap to get the standard errors.
        bs_wauc = vector(mode="numeric",length=100)
        bs_wcind = vector(mode="numeric",length=100)
        bs_cslope = vector(mode="numeric",length=100)
        for (brep in c(1:100)){
          # Boostrap repetitions:
          #2. Draw a bootstrap sample with replacement from the sample
          j = sample(nrow(ds),replace=T)
          bs_wauc[brep]=WeightedAUC(WeightedROC(staging_ranks[j], ds[j,"Metastasis"],weight=weights[j]))
          bs_wcind[brep]=intsurv::cIndex(ds[j,"Vitfup_metastasis_years"],event=ds[j,"Metastasis_numeric"],staging_ranks[j],weight=weights[j])[[1]]
          if (staging_system =="BWH"){
            bs_cslope[brep]=coef(coxph(Surv(ds[j,"Vitfup_metastasis_years"], ds[j,"Metastasis_numeric"]) ~ log_log_prob_BWH[j],weights = weights[j]))[1]  
          }
          else{
            bs_cslope[brep]=NA
          }
               
        }
        
        discriminative_performances_mimp_mean[paste0(staging_system,"wAUC"),i] = mean(bs_wauc)
        discriminative_performances_mimp_mean[paste0(staging_system,"wCind"),i] = mean(bs_wcind)
        discriminative_performances_mimp_mean[paste0(staging_system,"wCslope"),i] = mean(bs_cslope)
        
        discriminative_performances_mimp_se[paste0(staging_system,"wAUC"),i] = sqrt(var(bs_wauc))
        discriminative_performances_mimp_se[paste0(staging_system,"wCind"),i] = sqrt(var(bs_wcind))
        discriminative_performances_mimp_se[paste0(staging_system,"wCslope"),i] = sqrt(var(bs_cslope))
        
        #discriminative_performances_mimp_mean[paste0(staging_system,"wAUC"),i] = WeightedAUC(WeightedROC(staging_ranks, dutch_ds[,"Metastasis"],weight=weights))
        #discriminative_performances_mimp_mean[paste0(staging_system,"wCind"),i]=intsurv::cIndex(dutch_ds[,"Vitfup_metastasis_years"],event=dutch_ds[,"Metastasis_numeric"],staging_ranks,weight=weights)[[1]]
      }  
    }
    
    # Pool
    discriminative_performances_mimp = as.data.frame(matrix(NA,nrow=2,ncol = 5))
    colnames(discriminative_performances_mimp) = c("AUC","Cind","wAUC","wCind","wCslope")
    rownames(discriminative_performances_mimp)=c("AJCC","BWH")
    for (rw in c(1:nrow(discriminative_performances_mimp_mean))){
      metric_name = gsub("AJCC|BWH","",rownames(discriminative_performances_mimp_mean)[rw])
      staging_system = gsub("wCind|AUC|wAUC|Cind|wCslope","",rownames(discriminative_performances_mimp_mean)[rw])
      pooled_metric = round(pool_auc(discriminative_performances_mimp_mean[rw,],discriminative_performances_mimp_se[rw,],10),3)
      discriminative_performances_mimp[staging_system,metric_name] =paste0(round(pooled_metric[2],2)," (",round(pooled_metric[1],2),"-",round(pooled_metric[3],2),")")
    }
    all_performances_mimp_multipleTP = rbind(all_performances_mimp_multipleTP,discriminative_performances_mimp) 
}
rownames(all_performances_mimp_multipleTP) = paste0(rep(rownames(discriminative_performances_mimp),3),rep(c("_1y","_3y","_5y"),each=2))
write.csv(all_performances_mimp_multipleTP,paste0(dir_results,"/RSCM 11 AUC_Cind_slope_staging_systems_wCI.csv"))

#-----------------------


#-----------------------
# 3. Decision curve analysis  
#-----------------------
# Net benefit: All models with cutoffs:

vars=c("Met_risk","BWH_risk_estimates","Decision_AJCC")
col_vars = c("deepskyblue3","mediumvioletred","seagreen")
var_labels = c("Absolute risk model","BWH","AJCC:T1 vs rest")
names(col_vars)=names(var_labels)=vars

for (imp in c("","_cc")){
  
  if(imp ==""){
    ds = df_surv
  }else if(imp =="_cc"){
    ds = df_surv[!i.completecase,]
  }
  
  #p_list=list()
  #seq_prevs = c(seq(0.02,0.05,0.01))
  #for(i in c(1:length(seq_prevs))){
  #  p_list[[i]]=dca(Metastasis_numeric ~ Met_risk+Decision_BWH+Decision_AJCC, 
  #      data = ds, 
  #      thresholds = seq(0, 0.2, by = 0.01),
  #      label = list(Met_risk = "Absolute risk model",Decision_BWH = "BWH T1/T2a vs T2b/T3",Decision_AJCC="AJCC T1 vs rest",Met_risk_extended = "Extended absolute risk model"),
  #      prevalence=seq_prevs[i]) %>%
  #    plot(smooth = TRUE)+ ggtitle(paste0("Assumed prevalence = ",as.character(round(seq_prevs[i],3)))) 
  #}
  #grid_dca_thresholds = grid.arrange(arrangeGrob(grobs=p_list))
  #ggsave(paste0(dir_folder_results,"/RSCM 11 DCA thresholds with extended",imp,".pdf"),width=16,height=10,grid_dca_thresholds )
  
  # Net benefit: All models with probability:
  p_list=list()
  #seq_prevs = c(seq(0.02,0.05,0.01))
  j=1
  for(xl in c(0.2,1)){
    seq_prevs = c(0.02,0.05)
    for(i in c(1:length(seq_prevs))){
      p_list[[j]]=dca(Metastasis_numeric ~ Met_risk+BWH_risk_estimates+Decision_AJCC, 
                      data = ds, 
                      thresholds = seq(0, xl, by = 0.01),
                      label = list(Met_risk = "Absolute risk model",BWH_risk_estimates = "BWH", Decision_AJCC = "AJCC: T1 vs rest"),
                      prevalence=seq_prevs[i]) %>%
        plot()+ ggtitle(paste0("Assumed prevalence = ",as.character(round(seq_prevs[i],3)))) 
    j=j+1
    }
  }
  grid_dca_prob = grid.arrange(arrangeGrob(grobs=p_list,ncol=2,nrow=2))
  ggsave(paste0(dir_results,"/RSCM 11 DCA probabilities",imp,".pdf"),width=16,height=12,grid_dca_prob)

  # With weighted decision curve analysis:
  pdf(paste0(dir_results,"/RSCM 11 DCA_surv_weighted probabilities ",imp,".pdf"),width=7,height=6)
  plot_dca_surv_wgh(ds,vars,col_vars,ds$Metastasis_numeric,ds$Vitfup_metastasis_years,c(-0.01,0.02),c(0,20),5,ds$weights,var_labels)
  dev.off()
}
#-----------------------

#-----------------------
# 4. Performance metrics: Auxilliary analysis, to check performance metrics at different cutoffs
#-----------------------
th_performance = as.data.frame(matrix(NA,ncol=6,nrow=8))
colnames(th_performance) = c("Sensitivity","Specificity","NPV","PPV","NCase","NControl")
rownames(th_performance) = c("BWH_T1vsrest","BWH_T2avsT2b","Decision_AJCC","Abs_1p","Abs_2p","Abs_3p","Abs_4p","Abs_5p")

#Make new rows for different cutoffs
ds = df_surv
ds = ds %>% mutate(BWH_T1vsrest = ifelse(BWH_staging=="T1",0,1),
                   BWH_T2avsT2b = ifelse(BWH_staging%in%c("T1","T2a"),0,1))
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
th_performance$PctHRiskclass = th_performance$NCase/(nrow(ds))
df_plot_pm = melt(as.matrix(th_performance))

p_pm = ggplot(df_plot_pm[!df_plot_pm$Var2%in%c("NCase","NControl"),], aes(x = Var2, y = value*100,col=Var1)) +scale_color_manual(values=c("deepskyblue","deepskyblue4","orange","olivedrab1","olivedrab2","olivedrab3","limegreen","olivedrab4")) +
  geom_point(size = 5) +
  geom_line(aes(group =  Var1))+theme_bw()+xlab("Metrics")+ylab("Value (%)")
ggsave(paste0(dir_int_results,"/RSCM 11 Perf Metrics.pdf"),width=16,height=5,p_pm)
