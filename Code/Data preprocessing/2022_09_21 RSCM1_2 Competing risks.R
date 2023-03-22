#---------------------------------------------------
# Aim: Investigate competing risks
# Author: B. Rentroia Pacheco
# Input: Dutch dataset and weights
# Output: Competing risk curves
#---------------------------------------------------

# 0. Setup
#-----------------------
# library(survival)
# library(cmprsk)

# 1. Load data
#-----------------------
NKR_data = read.csv(paste0(dir_project,"Intermediate Data/2022_06_02 Recomputed FUP by Loes/2022_06_08 NKR_PALGA_FU_data.csv"))
# This dataset has the old version of the follow up. We only load it to get the follow up of the metastasis that are not included in the NCC cohort.
ncr_palga_db_incorrect_fup =  fread(paste0(dir_project,"Intermediate Data/Intermediate datasets old follow up 2022 06 17/2021_12_15 RSCM A_3 NCR_PALGA PA and FU_zonderbetrouwbaarheid3.csv"))

#This contains all unique cSCC records for patients in the 2007/08 cohort:
full_cohort_records = read.csv(paste0(dir_int_output,"RSCM 5_2_1 2022_06_15 NCRPALGA Tumor records preprocessed one per potential control.csv"))

#This includes all the patients with metastasis in the 2007/08 cohort, even if they were not included in the NCC cohort:
full_cohort_including_metastasis = read.csv(paste0(dir_int_output,"RSCM 5_2_1 2022_06_15 NCRPALGA Tumor records preprocessed one per potential control including metastasis.csv"))

load(paste0(dir_int_output,"RSCM1 Dutch Dataset Preprocessed_v3_0.RData")) #dutch_ds


# 2. Investigate competing events
#----------------------
# Get death information
NKR_data$Vitstat = NKR_data$vitstat2020
dutch_ds$row_order = 1:nrow(dutch_ds)
dutch_ds=merge(dutch_ds,NKR_data[,c("key_NKR","Vitstat")],by.x = "Key_NKR",by.y = "key_NKR",all.x=TRUE)
dutch_ds$event = ifelse(dutch_ds$Metastasis=="Case",1,0)
dutch_ds$event = ifelse(dutch_ds$Metastasis=="Control" & dutch_ds$Vitstat==1,2,dutch_ds$event)
dutch_ds$event = factor(dutch_ds$event,0:2,labels = c("Censor","Metastasis","Death"))
dutch_ds = dutch_ds[order(dutch_ds$row_order),]
dutch_ds$row_order = NULL
# Plot competing events:
df_plot = melt(table(dutch_ds$event,dutch_ds$Metastasis))
colnames(df_plot)=c("Event","Metastasis","Frequency")
ggplot(df_plot,aes(x=Metastasis,y=Frequency,fill=Event))+geom_bar(stat="identity",position = position_dodge(0.45))+theme_bw()+geom_text(aes(label=Frequency),position=position_dodge(width=0.5),vjust=-0.25)

# Plot Cumulative Incidence Curves
# CIF estimates (code based : https://www.ahajournals.org/doi/suppl/10.1161/CIRCULATIONAHA.115.017719)
event_num = as.numeric(as.character(revalue(as.character(dutch_ds$event),c("Censor"=0,"Metastasis"=1,"Death"=2))))
cif_metastasis <- cuminc(ftime=dutch_ds$Vitfup_metastasis_years,fstatus=event_num,cencode=0)
pdf(paste0(dir_results,"RSCM1 CIF plot.pdf"),width=20,height=5)
par(mfrow=c(1,3))
plot(cif_metastasis$"1 1"$time,cif_metastasis$"1 1"$est,type="l",
     ylim=c(0,1),
     xlab="Survival time (years)",ylab="Probability",xlim=c(0,10),
     lty=1,col="red",cex.lab=1.7,cex.axis=2,cex.sub=2)
#title("Cumulative Incidence functions")
lines(cif_metastasis$"1 2"$time,cif_metastasis$"1 2"$est,type="l",
      lty=1,col="blue")
legend("topleft",legend=c("Metastasis (CIF)","Death (CIF)","Metastasis (1-KM)","Death (1-KM)"),lty = c(1,1,2,2), col = c("red","blue","red","blue"),bty="n",cex=2)
# KM estimates
# Metastasis
km.mets <- survfit(Surv(dutch_ds$Vitfup_metastasis_years,ifelse(dutch_ds$Metastasis=="Case",1,0)) ~ 1)
lines(km.mets$time,1-km.mets$surv,type="l",lty=2,col="red")
# Death
km.death <- survfit(Surv(dutch_ds$Vitfup_metastasis_years,ifelse(event_num==2,1,0)) ~ 1)
lines(km.death$time,1-km.death$surv,type="l",lty=2,col="blue")
# What's the difference between the cumulative incidence at 5 years?
cif_metastasis$"1 1"$est[344]-(1-km.mets$surv[220])

#Weighted Plot:
#Weights for development and evaluation:
psamp = read.csv(paste0(dir_int_output,"RSCM 5_3 Sampling weights VitFUP and PA.csv"))
psamp_rordered = psamp[match(dutch_ds$Administratie_nr,psamp$pat.ids),]
stopifnot(length(which(psamp_rordered$Metastasis!=dutch_ds$Metastasis_numeric))==0)
weights = psamp_rordered$weights*12203/sum(psamp_rordered$weights) #Weight rescaling. We use the sum of potential controls + cases = 12203 in the full cohort

# weighted KM estimates:
# Metastasis
km.mets <- survfit(Surv(dutch_ds$Vitfup_metastasis_years,ifelse(dutch_ds$Metastasis=="Case",1,0)) ~ 1,weights = weights)
plot(km.mets$time,1-km.mets$surv,type="l",lty=2,col="red",ylim=c(0,0.6),xlim=c(0,10),xlab="Survival time (years)",ylab="Probability",cex.lab=1.7,cex.axis=2,cex.sub=2)
# Death
km.death <- survfit(Surv(dutch_ds$Vitfup_metastasis_years,ifelse(event_num==2,1,0)) ~ 1,weights = weights)
lines(km.death$time,1-km.death$surv,type="l",lty=2,col="blue")
# weighted Finegray:
pdata <- finegray(Surv(dutch_ds$Vitfup_metastasis_years, as.factor(event_num)) ~ 1,weights=weights)
fgfit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1,
               weight=fgwt, data=pdata)
fg_surv = survfit(fgfit)
lines(fg_surv$time,1-fg_surv$surv,type="l",lty=1,col="red")
legend("topleft",legend=c("Metastasis (CIF)","Death (CIF)","Metastasis (1-KM)","Death (1-KM)"),lty = c(1,1,2,2), col = c("red","blue","red","blue"),bty="n",cex=2)

pdata <- finegray(Surv(dutch_ds$Vitfup_metastasis_years, as.factor(revalue(as.character(event_num),c("1"="2","2"="1")))) ~ 1,weights=weights)
fgfit_death <- coxph(Surv(fgstart, fgstop, fgstatus) ~ 1,
                     weight=fgwt, data=pdata)
fg_surv_death = survfit(fgfit_death)
lines(fg_surv_death$time,1-fg_surv_death$surv,type="l",lty=1,col="blue")

#Plot on the full cohort:
# We consider a random record of each patient
full_cohort_crisk = NKR_data
full_cohort_crisk$Metastasis = ifelse(full_cohort_crisk$key_NKR%in%unique(full_cohort_including_metastasis$KEY_NKR[which(full_cohort_including_metastasis$metastasis==1)]),1,0)
# Select a random record for each:
# randomize record order
full_cohort_records = full_cohort_records[sample(1:nrow(full_cohort_records),nrow(full_cohort_records),replace=FALSE),]
# pick the first per administration number
full_cohort_one_record = full_cohort_records[!duplicated(full_cohort_records$administratienummer),]
#This is the follow up counting from the first cSCC in the NKR dataset:
full_cohort_crisk$Vitfup = as.Date(full_cohort_crisk$end_fu_2020,"%d%B%Y") - as.Date(full_cohort_crisk$date_first_cscc_nkr,"%d%B%Y") 
#This is the follow up of a randomly selected cSCC record:
full_cohort_crisk$Vitfup_random_record =  full_cohort_crisk$Vitfup
full_cohort_crisk$Vitfup_random_record = full_cohort_one_record$VitFU_timeyears[match(full_cohort_crisk$administratienummer,full_cohort_one_record$administratienummer)]
#Now we replace the follow up by te follow up until metastasis, for patients that experienced a metastasis:
full_cohort_crisk$Vitfup_metastasis = ncr_palga_db_incorrect_fup$FU_metastasis_years[match(full_cohort_crisk$key_NKR,ncr_palga_db_incorrect_fup$KEY_NKR)]
full_cohort_crisk$Vitfup_metastasis[which(full_cohort_crisk$Metastasis==0)] = full_cohort_crisk$Vitfup_random_record[which(full_cohort_crisk$Metastasis==0)]

#Solve some discrepancies in the follow up:
#No negative follow up
full_cohort_crisk$Vitfup_metastasis[which(full_cohort_crisk$Vitfup_metastasis<0)]=0
#If patient only has records in NKR, use Vitfup variable. These are patients that had a record in NKR, but not in PALGA.
full_cohort_crisk[is.na(full_cohort_crisk$Vitfup_metastasis),"Vitfup_metastasis"] = full_cohort_crisk$Vitfup[is.na(full_cohort_crisk$Vitfup_metastasis)]/365
full_cohort_crisk = full_cohort_crisk[!is.na(full_cohort_crisk$Vitfup_metastasis),]

#Code events as follows: 0 is censored, 1 is metastasis, 2 is death
full_cohort_crisk$event = 2*full_cohort_crisk$vitstat2020
full_cohort_crisk$event[which(full_cohort_crisk$Metastasis==1)]=1

cif_metastasis <- cuminc(ftime=full_cohort_crisk$Vitfup_metastasis,fstatus=full_cohort_crisk$event,cencode=0)
plot(cif_metastasis$"1 1"$time,cif_metastasis$"1 1"$est,type="l",
     ylim=c(0,0.6),
     xlab="Survival time (years)",ylab="Probability",xlim=c(0,10),
     lty=1,col="red",cex.lab=1.7,cex.axis=2,cex.sub=2)
#title("Cumulative Incidence functions")
lines(cif_metastasis$"1 2"$time,cif_metastasis$"1 2"$est,type="l",
      lty=1,col="blue")
legend("topleft",legend=c("Metastasis (CIF)","Death (CIF)","Metastasis (1-KM)","Death (1-KM)"),lty = c(1,1,2,2), col = c("red","blue","red","blue"),bty="n",cex=2)
# KM estimates
# Metastasis
km.mets <- survfit(Surv(full_cohort_crisk$Vitfup_metastasis,full_cohort_crisk$Metastasis) ~ 1)
lines(km.mets$time,1-km.mets$surv,type="l",lty=2,col="red")
# Death
km.death <- survfit(Surv(full_cohort_crisk$Vitfup_metastasis,ifelse(full_cohort_crisk$event==2,1,0)) ~ 1)
lines(km.death$time,1-km.death$surv,type="l",lty=2,col="blue")
dev.off()


dutch_ds$Vitstat=NULL