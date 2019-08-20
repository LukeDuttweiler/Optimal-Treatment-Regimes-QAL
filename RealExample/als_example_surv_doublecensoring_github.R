### created at 13/07/2019 the sensitivity analysis real example for ALS data set by modeing two types of censoring

#==============================#
#== step 1: data processing ===#
#==============================#

library(sas7bdat)
library(dplyr)
setwd("~/Desktop/HSun/Research Project/Causal inference/ALSResearch/CeftriaxoneTrial/Erin data/data/20170316")

allcensoreddata <- read.sas7bdat("alldata_censored.sas7bdat")
#alldata <- read.sas7bdat("alldata.sas7bdat")
basecensoreddata <- read.sas7bdat("baseline_censored.sas7bdat")
#alsf <- read.sas7bdat("D:/HSun/Research Project/Causal inference/ALSResearch/CeftriaxoneTrial/data/
#Ceftriaxone_deidentified_dataset/alsf.sas7bdat")

L = 120
K = 4
S = L/K
### select the variables: the first line, baseline variables; second line, time variable;
### third line, treatment and censoring variables; forth line, the others are time varying confounders.
work_dat = allcensoreddata %>% select(
  id, sex, Age, race, siteofonset, months_symp_to_diag, Ceftriaxone, Riluzole,bmi0,
  week, week2, week3,
  als_total, als_bulbar, als_breathing, als_gross_motor, als_fine_motor,
  maxfvc, bmi_pct_change_base, meangrip, alssqol_average_total) %>%
  mutate(QOL = 1)

work_dat$race[work_dat$race != "Cauciasian"] = "Other Race"
### organize the data into the required structure: gap 4 weeks
# the stage of peg insertion
last_week = basecensoreddata$last_week
g_tube = basecensoreddata$week_of_Gtube

last_week_l = pmin(last_week, L)
# bound at week 120
g_tube[g_tube >= L] = NaN
g_tube[g_tube <= 4] = 4
g_tube = round(g_tube %/% K)

# the censoring endpoint and round the censoring time
censor = pmax(basecensoreddata$Eventual_LTFU, basecensoreddata$Eventual_CensEOF)
cenLTFU = basecensoreddata$Eventual_LTFU
cenEOF = basecensoreddata$Eventual_CensEOF
censor[last_week >= L] = 0 # because the censoring occures after our upper bound, they are not censored before L
cenLTFU[last_week >= L] = 0 #
cenEOF[last_week >= L] = 0 #

qol_dat = c()
for(i in 1:481){
  subdat = work_dat[work_dat$id == i, ]
  refined_week = (censor[i]!=0) * round(last_week_l[i]/K)*K + (1 - (censor[i]!=0)) * last_week_l[i]
  s = ceiling(refined_week/K) + censor[i]
  
  if(censor[i] == 0){
    filter_subdat  = subdat %>% filter(week %% K == 0 & week <= L)
    if(refined_week %% 4 == 0) filter_subdat = filter_subdat[-dim(filter_subdat)[1],]
    #if(refined_week %% 4 > 0)  filter_subdat = rbind(filter_subdat, filter_subdat[dim(filter_subdat)[1],])
    
    filter_subdat$stage = seq(0, (dim(filter_subdat)[1]-1))
    filter_subdat$cen = filter_subdat$cenLTFU = filter_subdat$cenEOF = 0
    
    if(!is.nan(g_tube[i])){
      filter_subdat$trt = c(rep(0, g_tube[i]), rep(1, dim(filter_subdat)[1]-g_tube[i]))
    }
    else{
      filter_subdat$trt = 0
    }
    if(refined_week %% 4 > 0){
      filter_subdat$outcome = c(rep(4, refined_week %/% 4), refined_week %% 4)
    }
    else{
      filter_subdat$outcome = 4
    }
  } else if (cenLTFU[i] == 1){ ## for LTFU
    filter_subdat  = subdat %>% filter(week %% K == 0 & week <= L)
    filter_subdat$stage = seq(0, (dim(filter_subdat)[1]-1))
    if(refined_week > last_week_l[i])  filter_subdat = rbind(filter_subdat, filter_subdat[dim(filter_subdat)[1],])
    filter_subdat$cen = filter_subdat$cenLTFU = c(rep(0, dim(filter_subdat)[1]-1),1)
    filter_subdat$cenEOF = 0
    
    if(!is.nan(g_tube[i])){
      filter_subdat$trt = c(rep(0, g_tube[i]), rep(1, dim(filter_subdat)[1]-g_tube[i]))
    }
    else{
      filter_subdat$trt = 0
    }
    filter_subdat$outcome = c(rep(4, refined_week %/% 4), refined_week %% 4)
  } else if (cenEOF[i] == 1){ ## for EOF
    filter_subdat  = subdat %>% filter(week %% K == 0 & week <= L)
    filter_subdat$stage = seq(0, (dim(filter_subdat)[1]-1))
    if(refined_week > last_week_l[i])  filter_subdat = rbind(filter_subdat, filter_subdat[dim(filter_subdat)[1],])
    filter_subdat$cen = filter_subdat$cenEOF = c(rep(0, dim(filter_subdat)[1]-1),1)
    filter_subdat$cenLTFU = 0
    
    if(!is.nan(g_tube[i])){
      filter_subdat$trt = c(rep(0, g_tube[i]), rep(1, dim(filter_subdat)[1]-g_tube[i]))
    }
    else{
      filter_subdat$trt = 0
    }
    filter_subdat$outcome = c(rep(4, refined_week %/% 4), refined_week %% 4)
  }
  
  qol_dat = rbind(qol_dat, filter_subdat)
}
qol_dat$maxfvc_square = qol_dat$maxfvc^2
qol_dat$maxfvc_cubic = qol_dat$maxfvc^3
qol_dat$time = qol_dat$stage + 1
qol_dat$stage2 = qol_dat$stage^2
qol_dat$stage3 = qol_dat$stage^3

#=========================#
#== step 2: fit ps model =#
#=========================#

### this step using HAL is pretty arbitrary and needs a lot of time due to computation complexity
### load the help functions
setwd("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL/code_github")
source("support_github.r")
### check some basic features of the data

cen = qol_dat %>% group_by(id) %>% summarise(cen = sum(cen))
cat("the censoring rate is", mean(cen$cen))
obs_surv = qol_dat %>% group_by(id) %>% summarise(outcome = sum(outcome))
cat("the mean observed survival time is", mean(obs_surv$outcome))
quantile(obs_surv$outcome)
hist(obs_surv$outcome)
obs_qal = qol_dat %>% group_by(id) %>% summarise(qal = sum(QOL * outcome))
cat("the mean observed qal is", mean(obs_qal$qal))
#a = dat %>% group_by(ID) %>% summarise(trt = sum(trt) > 0)
#cat("the mean observed treatment rate is", mean(a$trt))

upp = quantile(obs_qal$qal, prob = 0.95)
# define the upper bound for QAL

### logistics regression:
# for the usage of thsi functions, check the 'ipw' package. I modify something based on the proposed method, but
# the structure is not changed.
# here, data = qol_dat[qol_dat$time > 1, ] is because we assume no peg is inserted at time 1.
trt <- ipwtm_qal(exposure = trt, family = "binomial", link = "logit",
                 numerator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
                   Ceftriaxone + Riluzole + bmi0,
                 denominator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
                   Ceftriaxone + Riluzole + bmi0 + als_total + als_bulbar + als_breathing + als_gross_motor + maxfvc +
                   bmi_pct_change_base + meangrip + alssqol_average_total,
                 id = id, tstart = stage, timevar = time, type = "first",
                 data = qol_dat[qol_dat$time > 1, ])

cen_LTFU <- ipwtm_qal(exposure = cenLTFU, family = "binomial", link = "logit",
                      numerator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
                        Ceftriaxone + Riluzole + bmi0,
                      denominator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
                        Ceftriaxone + Riluzole + bmi0 + als_total + als_bulbar + als_breathing + als_gross_motor + maxfvc +
                        bmi_pct_change_base + meangrip + alssqol_average_total,
                      id = id, tstart = stage, timevar = time, type = "first",
                      data = qol_dat[qol_dat$time > 1, ])

# cen_EOF<- ipwtm_qal(exposure = cenEOF, family = "binomial", link = "logit",
#                       numerator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
#                         Ceftriaxone + Riluzole + bmi0,
#                       denominator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
#                         Ceftriaxone + Riluzole + bmi0 + als_total + als_bulbar + als_breathing + als_gross_motor + maxfvc +
#                         bmi_pct_change_base + meangrip + alssqol_average_total,
#                       id = id, tstart = stage, timevar = time, type = "first",
#                       data = qol_dat[qol_dat$time > 1, ])

library(survival)
cen_EOF = Surv(time = obs_surv$outcome, event = cenEOF)
cen_EOF_surv = survfit(cen_EOF~1)
cen_EOF_surv = cbind(cen_EOF_surv$time, cen_EOF_surv$surv)
### record the surv prob at each landmark
cen_EOF_surv = cen_EOF_surv[cen_EOF_surv[,1] %% K == 0,]
cen_EOF_surv[,1] = cen_EOF_surv[,1] / K
EOF_surv_weight = data.frame(stage = 1:S, weight = 1)
for(s in 1:S){
  if(s == 1) EOF_surv_weight$weight[s] = 1
  else{
    if(s %in% cen_EOF_surv[,1]) EOF_surv_weight$weight[s] = 1/cen_EOF_surv[cen_EOF_surv[,1] == s,2]
    else{
      EOF_surv_weight$weight[s] = EOF_surv_weight$weight[s-1]
    }
  }
}


### construct the data for the main function.
# here, x1, x2, ... is user-specified. You can add/create the covariates you are interested in. These covariates can be
# considered in the optimal regimes.
# notice that x1 and x2 are the covariates used in the real example section of chapter 3.
clean_dat = qol_dat %>% select(ID = id, stage = time, x1 = maxfvc, x2 = bmi_pct_change_base,
                               x3 = als_bulbar, x1_square = maxfvc_square, x1_cubic = maxfvc_cubic,
                               qol = QOL, trt, outcome, censor = cen, start = stage)
### assign the weight function.
clean_dat$weight = 1
clean_dat$weight[clean_dat$stage > 1] = EOF_surv_weight$weight[clean_dat$stage[clean_dat$stage > 1]] * cen_LTFU$ipw.weights * trt$ipw.weights
clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame

### the test time points when computing RMST
a = c(0, sort(unique(pmin(obs_qal$qal, upp))))

### the idx of subjects censored or not
final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
uncen_idx = which(final_cen$final_cen == 0)


#=========================#
# step 3: optimal regime =#
#=========================#

### load the optimizer
library(rgenoud)

### the sample size
n = 481

#=====================================================================================================#
#=========== the example when only bmi and fev are considered ========================================#
#=========== this is the case used in chapter 3 ======================================================#
#=====================================================================================================#

temp1 <- genoud(fn=indicator_opt1_als, index = c(3, 4), nvars=3, default.domains=1, starting.values=rep(0, 3),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod1 = loss_indicator2_als(n = 481, eta = temp1$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))

temp2 <- genoud(fn=smooth_opt2_als, index = c(3, 4), nvars=3, default.domains=1, starting.values=rep(0, 3),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod2 = loss_indicator2_als(n = 481, eta = temp2$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))

# these are the benchmar regimes
### no peg
loss_indicator2_als(n = 481, eta = c(-1, 0, 0), clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))
### immediate peg
loss_indicator2_als(n = 481, eta = c(1, 0, 0), clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))
### fvc less than 50
loss_indicator2_als(n = 481, eta = c(50, -1, 0), clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))
### bmi reduce more than 10%
loss_indicator2_als(n = 481, eta = c(-10, 0, -1), clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))


#=====================================================================================================#
#=========== the example when bmi, fvc and bulbar are considered =====================================#
#=========== this is the case suggested by emory =====================================================#
#=====================================================================================================#

temp3 <- genoud(fn=indicator_opt1_als, index = c(3, 4, 5), nvars=4, default.domains=1, starting.values=rep(0, 4),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod3 = loss_indicator2_als(n = 481, eta = temp3$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))

temp4 <- genoud(fn=smooth_opt2_als, index = c(3, 4, 5), nvars=4, default.domains=1, starting.values=rep(0, 4),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod4 = loss_indicator2_als(n = 481, eta = temp4$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))

# # these are the benchmar regimes
# ### no peg
# loss_indicator2_als(n = 481, eta = c(-1, 0, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))
# ### immediate peg
# loss_indicator2_als(n = 481, eta = c(1, 0, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))
# ### fvc less than 50
# loss_indicator2_als(n = 481, eta = c(50, -1, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))
# ### bmi reduce more than 10%
# loss_indicator2_als(n = 481, eta = c(-10, 0, -1, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))

#=====================================================================================================#
#=========== the example when bmi, fvc^2 are considered ===================================#
#=====================================================================================================#

temp5 <- genoud(fn=indicator_opt1_als, index = c(3, 4, 6), nvars=4, default.domains=1, starting.values=rep(0, 4),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod5 = loss_indicator2_als(n = 481, eta = temp5$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 6))

temp6 <- genoud(fn=smooth_opt2_als, index = c(3, 4, 6), nvars=4, default.domains=1, starting.values=rep(0, 4),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod6 = loss_indicator2_als(n = 481, eta = temp6$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 6))


#=====================================================================================================#
#=========== the example when bmi, fvc^2 and bulbar are considered ===================================#
#=====================================================================================================#

temp7 <- genoud(fn=indicator_opt1_als, index = c(3, 4, 5, 6), nvars=5, default.domains=1, starting.values=rep(0, 5),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod7 = loss_indicator2_als(n = 481, eta = temp7$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5, 6))

temp8 <- genoud(fn=smooth_opt2_als, index = c(3, 4, 5, 6), nvars=5, default.domains=1, starting.values=rep(0, 5),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod8 = loss_indicator2_als(n = 481, eta = temp8$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5, 6))

#=====================================================================================================#
#=========== the example when bmi, fvc^3 and bulbar are considered ===================================#
#=====================================================================================================#

temp9 <- genoud(fn=indicator_opt1_als, index = c(3, 4, 6, 7), nvars=5, default.domains=1, starting.values=rep(0, 5),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod9 = loss_indicator2_als(n = 481, eta = temp9$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 6, 7))

temp10 <- genoud(fn=smooth_opt2_als, index = c(3, 4, 6, 7), nvars=5, default.domains=1, starting.values=rep(0, 5),
                 max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod10 = loss_indicator2_als(n = 481, eta = temp10$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 6, 7))

#=====================================================================================================#
#=========== the example when bmi, fvc^3 and bulbar are considered ===================================#
#=====================================================================================================#

temp11 <- genoud(fn=indicator_opt1_als, index = c(3, 4, 5, 6, 7), nvars=6, default.domains=1, starting.values=rep(0, 6),
                 max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod11 = loss_indicator2_als(n = 481, eta = temp11$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5, 6, 7))

temp12 <- genoud(fn=smooth_opt2_als, index = c(3, 4, 5, 6, 7), nvars=6, default.domains=1, starting.values=rep(0, 6),
                 max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod12 = loss_indicator2_als(n = 481, eta = temp12$par, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5, 6, 7))


### save the optimation results
real_logit_surv_double = list(temp1 = temp1, temp2 = temp2, temp3 = temp3, temp4 = temp4 ,
                      temp5 = temp5, temp6 = temp6, temp7 = temp7, temp8 = temp8, temp9 = temp9, temp10 = temp10, 
                      temp11 = temp11, temp12 = temp12)
setwd("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL/code_github")
save(real_logit_surv_double, file = "real_logit_surv_double.rda")





#=====================================================================================================#
#================ generate the survival function =====================================================#
#=====================================================================================================#

index = c(3,4)
load("real_logit_surv_double.rda")


### IPW estimator: 
surv_curve = function(eta, index){
  n_dat = clean_dat
  
  pp = c()
  ### construct the follow indicator: need to be complex
  id = clean_dat$ID
  
  ### construct the design matrix for optiml regime (here, we include the intercept and interested time-varying covariates)
  regime_design = cbind(1, as.matrix(n_dat[, index]))
  ### construct the follow indicator: need to be complex
  id = clean_dat$ID
  trt_eta_raw  = (regime_design %*% as.matrix(eta)) > 0
  
  #trt_eta_raw  = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0
  trt_eta = unlist(lapply(1:length(unique(id)), function(x) pmin(cumsum(trt_eta_raw[id == x]), 1)))
  
  n_dat$follow = n_dat$trt * trt_eta_raw + (1-n_dat$trt) * (1 - trt_eta_raw)
  #n_dat$follow = n_dat$trt * ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0) +
  #  (1-n_dat$trt) * (1 - ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0))
  n_dat$follow[clean_dat$stage == 1] = 1
  n_dat = n_dat %>% group_by(ID) %>% mutate(follow_cum = cumprod(follow)) %>% as.data.frame
  
  ### define the estimation of t-year survival probability: see section 3.1 of the paper
  prob_surv <- function(i){
    sub1 = n_dat[n_dat$qal >a[i],]
    sub1 = sub1[!duplicated(sub1$ID),]
    numerator = sum(sub1$follow_cum * sub1$weight)
    
    if(numerator == 0){
      return(0)
    }
    else if(all(uncen_idx %in% unique(sub1$ID))){
      return(1)
    }
    else{
      uncen_supp = n_dat[(n_dat$ID %in% uncen_idx) & !(n_dat$ID %in% sub1$ID),]
      uncen_supp_last = uncen_supp[!duplicated(uncen_supp$ID, fromLast = T),]
      sub2 = rbind(sub1, uncen_supp_last)
      denominator = sum(sub2$follow_cum * sub2$weight)
      
      return(numerator/denominator)
    }
    
  }
  
  ### the discrete integration for RMST
  pp = unlist(lapply(1:(length(a) - 1), function(x) prob_surv(x)))
  
  return(pp)
}

### the opt case
sc_opt = surv_curve(temp1$par, index)
sc_no = surv_curve(c(-1, 0, 0), index)
sc_always = surv_curve(c(1, 0, 0), index)
sc_fvc = surv_curve(c(50, -1, 0), index)
sc_bmi = surv_curve(c(-10, 0, -1), index)

path <- paste0("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL")
pdf(paste0(path, "/real_survplot_double_surv.pdf"), width = 5, height  = 5)

par(mfrow = c(1,1))
plot(a[-length(a)], sc_always, type = "l", lty = 3, xlab = "Survival time (QAL, weeks)", ylab = "Surv.curve P(T > t)", lwd = 2)
points(a[-length(a)], sc_no, type = "l", lty = 2, lwd = 2)
points(a[-length(a)], sc_fvc, type = "l", lty = 2, col = 2, lwd = 2)
points(a[-length(a)], sc_bmi, type = "l", lty = 2, col = 4, lwd = 2)
points(a[-length(a)], sc_opt, type = "l", lty = 1, lwd = 2)
legend("bottomleft", c("IPW", "BMI", "FVC", "No PEG", "Always PEG"), lty = c(1,2,2,2,3), cex = 1, col = c(1,4,2,1,1), lwd = c(2,2,2,2,2))

dev.off()



#=====================================================================================================#
#================ generate the profile curve =========================================================#
#=====================================================================================================#
dat = clean_dat[clean_dat$ID == 23,]
score23 = par[1] + as.matrix(dat[,index]) %*% as.matrix(par[-1])

dat = clean_dat[clean_dat$ID == 27,]
score27 = par[1] + as.matrix(dat[,index]) %*% as.matrix(par[-1])

dat = clean_dat[clean_dat$ID == 63,]
score63 = par[1] + as.matrix(dat[,index]) %*% as.matrix(par[-1])
