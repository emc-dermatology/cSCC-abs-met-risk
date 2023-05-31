# Functions to compute weighted versions of the performance metrics. 
# These functions are adapted versions of the original ones

weighted_scaled_brier = function (obs,pred,weights=NULL){
  if(is.null(weights)){
    scaled_brier(obs,pred)
  }else{
    1 - (weighted.mean((obs - pred)^2,weights) / ((weighted.mean(obs,weights)) * (1 - weighted.mean(obs,weights))))
  }
}

# Function to produce the calibration plots
calib_plot =function(predictions,outcome_num,outcome_FUP,ds,u,bin_nr,lim,pl,outcome_type,weights=NULL,md,log.calib=TRUE,smooth=TRUE,log_plot=FALSE,main_plot="",surv_or_met = "surv",new_plot = TRUE,dots_col = "darkred"){
  
  #Truncate outcome at u time
  outcome_num_u = outcome_num
  i.greaterFUP = which(outcome_FUP>u)
  if(length(i.greaterFUP)>0){
    outcome_num_u[outcome_FUP>u]=0
  }
  
  outcome_num_u = ifelse(outcome_num_u==1,0,1)
  if(length(i.greaterFUP)>0){
    outcome_FUP[outcome_FUP>u] = u
  }
  # TODO: there is a bug with the truncation that needs to be fixed.
  
  # Get calibration slope:
  logit <- qlogis(predictions)
  i <- ! is.infinite(logit)
  logit_notinf = logit[i]
  outcome_not_inf = outcome_num_u[i]
  f.recal <- c(coef(glm(outcome_not_inf~logit_notinf,family="binomial",weights=weights[i])),NA,NA)
  names(f.recal)=c("Intercept_LR","Slope_LR","Intercept_surv","Slope_surv")
  #Note, if predictor is constant, then slope will be NA:
  if(is.na(f.recal["Slope_LR"])){f.recal["Slope_LR"]=0}
  
  if(outcome_type=="survival"){
    # Intercept:
    if(md %in%c("uncoxnet","wghcoxnet","wghcoxnet_median","wghcoxnet_average")){
      fit_intercept = NA
    }else{
      # We can compute model intercept at fixed time point or entire follow up time. 
      # We focused on fixed time point estimates, as we are interested in time point estimates at 3 or 5y.
      # For estimates for the entire follow up time, a poisson model should be used.
      #ds_intcp = truncateFUP(ds,"Metastasis_numeric","Vitfup_metastasis_years","Metastasis",u)
      #p = log(predict(fit_model, newdata=ds_intcp,type="expected",reference="zero"))
      #i.notzero = is.finite(p)
      #fit_intercept = coef(glm(outcome_num[i.notzero] ~ offset(p[i.notzero]), family="poisson",weights=weights[i.notzero]))
      if(is.null(weights)){
        obj = summary(survfit(Surv(outcome_FUP,outcome_num) ~ 1), 
                      times = u,extend=TRUE)
        fit_intercept = (1 - obj$surv) / (1-mean(predictions))
      }else{
        obj = summary(survfit(Surv(outcome_FUP,outcome_num) ~ 1,weights = weights), 
                      times = u,extend=TRUE)
        fit_intercept = (1 - obj$surv) / (1-weighted.mean(predictions,weights))
      }
      
    }
    # Slope
    if(md %in%c("uncoxnet","wghcoxnet","wghcoxnet_median","wghcoxnet_average")){
      #print("h")
      lp = get_cox_linear_predictor(fit_model,ds,ds,fm_model,FALSE) #We have checked that this gives the same result as predict function
      #lp = predict(fit_model,newx=ds,type="link",reference="zero")
      fit_slope = coef(coxph(Surv(outcome_FUP,outcome_num)~lp,weights=weights))
      if (is.na(fit_slope)){fit_slope=0} # This happens when coefficients are 0
    }else{
      log_log_pred = log(-log(predictions))
      #lp = predict(fit_model,newdata=ds,type="lp",reference="zero")
      #fit_slope=coef(coxph(Surv(outcome_FUP,outcome_num)~lp,data=ds,weights=weights))
      fit_slope=coef(coxph(Surv(outcome_FUP,outcome_num)~log_log_pred,data=ds,weights=weights))
    }
  }
  # Smooth calibration: 
  Sm <- lowess(predictions, outcome_num_u, iter=0)
  
  if(outcome_type=="survival"){
    if(log_plot|(surv_or_met=="met")){
      ylab=paste0("1-Kaplan-Meier survival ",u,' years') 
      xlab = "Predicted metastastic probability"
    }else{
      ylab=paste0("Kaplan-Meier survival ",u,' years') 
      xlab = "Predicted survival probability"
    }
  }else if (outcome_type=="binary"){
    ylab=paste0("Actual probability of survival( ",u,' years)')
    predictions= predictions
    xlab = "Predicted survival probability"
  }
  # Plot:
  if(pl){
    if(new_plot){
      if(log_plot){
        plot(.5, .5, xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab,log ="xy",main=main_plot,cex.lab=1.5, cex.axis=1.5)
      }else{
        plot(.5, .5, xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab,main=main_plot,cex.lab=1.5, cex.axis=1.5)
      }
      abline(0, 1, lwd=6, col=gray(.85))
    }
    
    lt <- 1; leg <- "Ideal"; marks <- -1; lwd <- 6; col <- gray(.85)
    # Logistic Regression
    if(log.calib& outcome_type=="binary"){
      lt <- c(lt, 1); leg <- c(leg, "Logistic calibration")
      lwd <- c(lwd, 1); col <- c(col, 'black')
      marks <- c(marks, -1)
      logit <- seq(-7, 7, length=200)
      prob  <- plogis(logit)
      pred.prob <- f.recal[1] + f.recal[2] * logit
      pred.prob <- plogis(pred.prob)
      lines(prob, pred.prob, lty=1)
    }
    #Non-parametric
    if(smooth & is.null(weights)){
      lines(Sm, lty=3)
      lt <- c(lt, 3)
      lwd <- c(lwd, 1); col <- c(col, 'black')
      marks <- c(marks, -1)
      leg <- c(leg, "Nonparametric")}
    logit <- seq(-7, 7, length=200)
    prob  <- plogis(logit)
    pred.prob <- f.recal[1] + f.recal[2] * logit
    pred.prob <- plogis(pred.prob)
    
    # Partition data
    bin <- cut2(predictions,m=bin_nr,levels.mean=TRUE,digits=7)
    means <- as.numeric(levels(bin))
    
    if (outcome_type=="survival"){
      Srv = Surv(outcome_FUP,outcome_num)
      # Adapted from rms's function groupKM, this could not be used directly because it didn't have weighted KM estimates
      q <- unclass(bin)  #need integer
      g <- length(levels(q))
      km      <- double(g)
      pred    <- km
      std.err <- km
      events  <- integer(g)
      numobs  <- events
      
      #Compute KM estimates for each group:
      for(i in 1:g){
        s <- q==i
        nobs <- sum(s); ne <- sum(outcome_num[s])
        
        pred[i] <- mean(predictions[s], na.rm=TRUE) #mean predictions of each group
        dummystrat <- as.factor(rep("1", nobs))
        if(is.null(weights)){
          #f <- survfitKM(dummystrat,Srv[s,]) 
          f = survfit(Srv[s,] ~ 1)
        }else{
          f = survfit(Srv[s,] ~ 1,weights =weights[s])
          #f <- survfitKM(dummystrat,Srv[s,],weights = weights[s]) 
        }
        
        #Note: if last u> last time, we take the estimates at the last time point:
        if(max(f$time)<u){
          print(paste0("Group",i," did not have any observation with survival time later than the end time point! The longest survival estimates were extended."))
          i.add = which.max(f$time)
          f$time = c(f$time,u)
          f$surv = c(f$surv,f$surv[i.add])
          f$std.err = c(f$std.err,f$std.err[i.add])
        }
        
        ##doesn't need conf.int since only need s.e.
        tt <- c(0, f$time)
        ss <- c(1, f$surv)
        se <- c(0, f$std.err)
        tm <- max((1:length(tt))[tt <= u+1e-6])
        km[i] <- ss[tm]
        std.err[i] <- se[tm]
        numobs[i]  <- nobs
        events[i]  <- ne
        n <- length(tt)
        if(u > tt[n]+1e-6 & ss[n]>0)
        {
          km[i] <- NA
          std.err[i] <- NA
        }
      }
      
      #For confidence intervals
      z <- cbind(x=pred, n=numobs, events=events, KM=km, 
                 std.err=std.err)
      #print(z)
      
      ciupper <- function(surv, d) ifelse(surv==0, 0, pmin(1, surv*exp(d)))
      cilower <- function(surv, d) ifelse(surv==0, 0, surv*exp(-d))
      prop = km
      names(prop)=pred
      
      conf.int=0.95
      zcrit <- qnorm((conf.int+1)/2)
      low <- cilower(km, zcrit*std.err)
      hi  <- ciupper(km, zcrit*std.err)
      if(conf.int){
        if (log_plot|(surv_or_met=="met")){
          hi = ifelse(hi==1,0.999,hi)
          errbar(1-pred, 1-km, 1-hi, 1-low, add=TRUE)
         }else{

             errbar(pred, km, hi, low, add=TRUE)

          
           
        }
            
      }
     }else if (outcome_type =="binary"){
      if(is.null(weights)){
        prop <- tapply(outcome_num_u, bin, function(x) mean(x, na.rm=TRUE))
      }else{
        prop = rep(0,length(levels(bin)))
        names(prop)=levels(bin)
        for(bin_uniq in unique(bin)){
          i.p= which(bin==bin_uniq)
          prop[as.character(bin_uniq)]=weighted.mean(outcome_num_u[i.p],w=weights[i.p])
        }
      }
    }
    
    if(log_plot|(surv_or_met =="met")){
      points(1-means, 1-prop, pch=16,col=dots_col)
      
    }else{

        points(means, prop, pch=16,col=dots_col)
      
      
    
    }
    #lines(means, prop,col="darkred"); lt <- c(lt, 1)
    leg <- c(leg, "Grouped observations")
    lt <- c(lt, 0)
    col <- c(col, dots_col); lwd <- c(lwd, 1)     
    marks <- c(marks, 16)
    
    #Add legend:
    #if(log_plot){
    #  legend(0.0001,0.8, leg, lty=lt, pch=marks, cex=1.1, lwd=lwd, col=col, bty="n")
    #}else{
    #  legend(0.01,1.05, leg, lty=lt, pch=marks, cex=1.1, lwd=lwd, col=col, bty="n")
    #}
    
    #Add probability distribution:
    if(log_plot|(surv_or_met =="met")){
      x<-1-predictions
    }else{
        x<-predictions

      
    }
   
    bins <- seq(lim[1], lim[2], length=101)
    x <- x[x >= lim[1] & x <= lim[2]]
    f <- table(cut(x, bins))
    j <- f > 0
    bins <- (bins[-101])[j]
    f <- f[j]
    #print(bins)
    #print(f)
 
    if(log_plot){ 
      #print(diff(lim))
      #print(max(f))
      #f <- lim[1] + .015 * diff(lim) * f / max(f)
      #print(f)
      #segments(bins, lim[1], bins, f) #This doesn't help because it seems like the probabilities are not concentrated below 0.1 when in fact they are
    }else{
      f <- lim[1] + .15 * diff(lim) * f / max(f)
      segments(bins, lim[1], bins, f,col =dots_col)
    }

  }
  if(outcome_type=="survival"){
    f.recal[c("Intercept_surv","Slope_surv")]=c(fit_intercept,fit_slope)
  }
  
  if (length(unique(predictions))==1){
    f.recal["Slope_LR"]=0
    f.recal["Slope_surv"]=0
  }
  return(f.recal)
}

# Function to compute net benefit:
# @param th: probability threshold
# @param ev: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param ev_tm: numeric vector with follow up time of subjects
# @param rp: numeric vector with model predictions
# @param tm: number indicating at which timepoint the observed survival should be computed.  It should be between 0 and max of outcome_FUP.
# @param wghs: if required output is weighted net benefit, weights is a numeric vector with weights of each observation. If required output is unweighted net benefit, weights should be NULL
compute_nb_surv_wgh = function(th,ev,ev_tm,rp,tm,wghs){
  
  total_n = length(ev)
  
  # Compute the proportion of subjects with risk probability larger than threshold th:
  p_x1 =  weighted.mean(ifelse(rp>th,1,0),wghs)
  hr_pts = which(rp>th)
  
  # Focus on high risk subjects:
  ev_tm_hr = ev_tm[hr_pts]
  ev_hr = ev[hr_pts]
  wghs_hr = wghs[hr_pts]

  # Compute observed survival for high-risk subjects
  if(length(hr_pts)>1){
    #fit <- coxph(Surv(ev_tm_hr, ev_hr) ~ 1,weights = wghs_hr) 
    expfit <- summary(survfit(Surv(ev_tm_hr, ev_hr) ~ 1,weights = wghs_hr),times=tm,extend = TRUE)$surv
  }else{
    expfit = ifelse(ev_hr[hr_pts]==1,0,1)
  }
  
  # Compute true positives and false positives. Note True positives = tp* total_n and False positives = fp*total_n. However, these would be divided by total_n in the net benefit computation. So the total_n would cancel out.
  tp = (1-expfit)*p_x1 # tp = (1-expfit)*p_x1*total_n
  fp = expfit*p_x1 # fp = expfit*p_x1*total_n
  
  # Compute net benefit at threshold th
  nb = tp-(fp)*(th/(1-th)) #nb = tp/total_n-(fp/total_n)*(th/(1-th))
  
  return(nb)
}

# Function to plot decision curve analysis:
# @param ds: dataframe containing risk probabilities to be evaluated with a decision curve. Each column in this dataframe should correspond to a different risk prediction method
# @param vars: character vector with column names of columns present in dataframe ds, which represent risk probabilities of different methods
# @param col_vars: character vector with the color id for each variable present in the vars vector
# @param ev: numeric vector with status indicator, 0 =alive, 1 = dead. Follows the same conventions as parameter event in the Surv function from the survival package.
# @param ev_time: numeric vector with follow up time of subjects
# @param ylim: vector of length two, with lower and upper limit for the y axis in the plot
# @param xlim: vector of length two, with lower and upper limit for the x axis in the plot
# @param tm: number indicating at which timepoint the data should be truncated at.  It should be between 0 and max of outcome_FUP.
# @param wghs: if required output is weighted decision curve, weights is a numeric vector with weights of each observation. If required output is unweighted decision curve, weights should be NULL
plot_dca_surv_wgh = function(ds,vars,col_vars,ev,ev_time,ylim,xlim,tm,wghs,var_labels=NULL){
  nb_th = seq(0,0.99,0.01)
  
  #Net benefit of treating all:
  nb_all=sapply(nb_th,compute_nb_surv_wgh,ev,ev_time,rep(1,length(ev)),tm,wghs)
  plot(nb_th*100,nb_all,type="l",ylim=ylim,xlim=xlim,ylab = "Net benefit",xlab="Threshold probability (%)",col="black",lwd=2)
  
  #Net benefit of methods of interest:
  for (var in vars){
    #Note: in the dcurves package, the model is refitted for the desired time
    outcome =Surv(time=ev_time,event=ev)
    variable=ds[,var]
    
    nb_var=sapply(nb_th,compute_nb_surv_wgh,ev,ev_time,ds[,var],tm,wghs)
    lines(nb_th*100,nb_var,lwd=2,col=col_vars[var])
  }
  
  #Net benefit of treating none:
  lines(nb_th*100,rep(0,length(nb_th)),col="darkgray",lwd=2)
  
  #Add a legend:
  if(is.null(var_labels)){
    legend("topright",legend=c("Treat all","Treat none",vars),lty = 1,lwd=2, col = c("black","darkgray",col_vars),bty="n")
  }else{
    legend("topright",legend=c("Treat all","Treat none",var_labels[vars]),lty = 1,lwd=2, col = c("black","darkgray",col_vars),bty="n")
  }
  
}