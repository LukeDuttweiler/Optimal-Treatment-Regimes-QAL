### created at 05/09/2019 the real example for ALS data set

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
  mutate(QOL = alssqol_average_total/max(alssqol_average_total))

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
censor[last_week >= 120] = 0 # because the censoring occures after our upper bound, they are not censored before L

qol_dat = c()
for(i in 1:481){
  subdat = work_dat[work_dat$id == i, ]
  refined_week = censor[i] * round(last_week_l[i]/K)*K + (1 - censor[i]) * last_week_l[i]
  s = ceiling(refined_week/K) + censor[i]

  if(censor[i] == 0){
    filter_subdat  = subdat %>% filter(week %% K == 0 & week <= L)
    if(refined_week %% 4 == 0) filter_subdat = filter_subdat[-dim(filter_subdat)[1],]
    #if(refined_week %% 4 > 0)  filter_subdat = rbind(filter_subdat, filter_subdat[dim(filter_subdat)[1],])

    filter_subdat$stage = seq(0, (dim(filter_subdat)[1]-1))
    filter_subdat$cen = 0

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
  }

  else{
    filter_subdat  = subdat %>% filter(week %% K == 0 & week <= L)
    filter_subdat$stage = seq(0, (dim(filter_subdat)[1]-1))
    if(refined_week > last_week_l[i])  filter_subdat = rbind(filter_subdat, filter_subdat[dim(filter_subdat)[1],])
    filter_subdat$cen = c(rep(0, dim(filter_subdat)[1]-1),1)

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
qol_dat$time = qol_dat$stage + 1
qol_dat$stage2 = qol_dat$stage^2
qol_dat$stage3 = qol_dat$stage^3



#=========================#
#== step 2: fit ps model =#
#=========================#

### load the help functions
setwd("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL/code_github")
source("support_github.r")

### describe the data; check some basic features of the data
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

### specify the upper limit of RQAL
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

cen <- ipwtm_qal(exposure = cen, family = "binomial", link = "logit",
                 numerator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
                   Ceftriaxone + Riluzole + bmi0,
                 denominator = ~ stage + stage2 + stage3 + sex  + Age + race + siteofonset + months_symp_to_diag +
                   Ceftriaxone + Riluzole + bmi0 + als_total + als_bulbar + als_breathing + als_gross_motor + maxfvc +
                   bmi_pct_change_base + meangrip + alssqol_average_total,
                 id = id, tstart = stage, timevar = time, type = "first",
                 data = qol_dat[qol_dat$time > 1, ])

### construct the data for the main function.
# here, x1, x2, ... is user-specified. You can add/create the covariates you are interested in. These covariates can be
# considered in the optimal regimes.
# notice that x1 and x2 are the covariates used in the real example section of chapter 3.
clean_dat = qol_dat %>% select(ID = id, stage = time, x1 = maxfvc, x2 = bmi_pct_change_base,
                               x3 = als_bulbar, x1_square = maxfvc_square, qol = QOL, trt, outcome, censor = cen, start = stage)
### assign the weight function.
clean_dat$weight = 1
clean_dat$weight[clean_dat$stage > 1] = cen$ipw.weights * trt$ipw.weights
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
mod1 = loss_indicator2_als(n = 481, eta = temp1$par, dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))

temp2 <- genoud(fn=smooth_opt2_als, index = c(3, 4), nvars=3, default.domains=1, starting.values=rep(0, 3),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod2 = loss_indicator2_als(n = 481, eta = temp2$par, dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))

# these are the benchmar regimes
### no peg
loss_indicator2_als(n = 481, eta = c(-1, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))
### immediate peg
loss_indicator2_als(n = 481, eta = c(1, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))
### fvc less than 50
loss_indicator2_als(n = 481, eta = c(50, -1, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))
### bmi reduce more than 10%
loss_indicator2_als(n = 481, eta = c(-10, 0, -1), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4))


#=====================================================================================================#
#=========== the example when bmi, fvc and bulbar are considered =====================================#
#=========== this is the case suggested by emory =====================================================#
#=====================================================================================================#

temp3 <- genoud(fn=indicator_opt1_als, index = c(3, 4, 5), nvars=4, default.domains=1, starting.values=rep(0, 4),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod3 = loss_indicator2_als(n = 481, eta = temp1$par, dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))

temp4 <- genoud(fn=smooth_opt2_als, index = c(3, 4, 5), nvars=4, default.domains=1, starting.values=rep(0, 4),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod4 = loss_indicator2_als(n = 481, eta = temp2$par, dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))

# these are the benchmar regimes
### no peg
loss_indicator2_als(n = 481, eta = c(-1, 0, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))
### immediate peg
loss_indicator2_als(n = 481, eta = c(1, 0, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))
### fvc less than 50
loss_indicator2_als(n = 481, eta = c(50, -1, 0, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))
### bmi reduce more than 10%
loss_indicator2_als(n = 481, eta = c(-10, 0, -1, 0), dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5))


### save the optimation results
real_logit = list(temp1 = temp1, temp2 = temp2, temp3 = temp3, temp4 = temp4)
setwd("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL/code_github")
save(real_logit, file = "real_logit.rda")
