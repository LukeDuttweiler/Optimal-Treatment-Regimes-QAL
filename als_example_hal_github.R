### created at 05/09/2019 the real example for ALS data set using HAL

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


### highly adaptive lasso
library(hal9001)

tempdat = qol_dat[qol_dat$time > 1, ]
tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$id),function(x)if (!is.na(match(1, x)))
  return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))

tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
#tempdat_trt[,-c(1,4)] = apply(tempdat_trt[,-c(1,4)],2,function(x) (x - min(x))/(max(x) - min(x)))
# tempdat_trt$als_total = (tempdat_trt$als_total - min(tempdat_trt$als_total))/(max(tempdat_trt$als_total) - min(tempdat_trt$als_total))
# tempdat_trt$als_bulbar = (tempdat_trt$als_bulbar - min(tempdat_trt$als_bulbar))/(max(tempdat_trt$als_bulbar) - min(tempdat_trt$als_bulbar))
# tempdat_trt$stage = (tempdat_trt$stage - min(tempdat_trt$stage))/(max(tempdat_trt$stage) - min(tempdat_trt$stage))

### set model at each time point
tempdat_sub1 = tempdat_trt[tempdat_trt$id %in% c(1:30), ]
tempdat_sub2 = tempdat_trt[tempdat_trt$id %in% c(31:60), ]
tempdat_sub3 = tempdat_trt[tempdat_trt$id %in% c(61:90), ]
tempdat_sub4 = tempdat_trt[tempdat_trt$id %in% c(91:120), ]
tempdat_sub5 = tempdat_trt[tempdat_trt$id %in% c(121:150), ]
tempdat_sub6 = tempdat_trt[tempdat_trt$id %in% c(151:180), ]
tempdat_sub7 = tempdat_trt[tempdat_trt$id %in% c(181:210), ]
tempdat_sub8 = tempdat_trt[tempdat_trt$id %in% c(211:240), ]
tempdat_sub9 = tempdat_trt[tempdat_trt$id %in% c(241:270), ]
tempdat_sub10 = tempdat_trt[tempdat_trt$id %in% c(271:300), ]
tempdat_sub11 = tempdat_trt[tempdat_trt$id %in% c(301:330), ]
tempdat_sub12 = tempdat_trt[tempdat_trt$id %in% c(331:360), ]
tempdat_sub13 = tempdat_trt[tempdat_trt$id %in% c(361:390), ]
tempdat_sub14 = tempdat_trt[tempdat_trt$id %in% c(391:420), ]
tempdat_sub15 = tempdat_trt[tempdat_trt$id %in% c(421:450), ]
tempdat_sub16 = tempdat_trt[tempdat_trt$id > 450, ]
# tempdat_sub1 = tempdat_trt[tempdat_trt$stage == 1, ]
# tempdat_sub2 = tempdat_trt[tempdat_trt$stage == 2, ]
# tempdat_sub3 = tempdat_trt[tempdat_trt$stage == 3, ]
# tempdat_sub4 = tempdat_trt[tempdat_trt$stage == 4, ]
# tempdat_sub5 = tempdat_trt[tempdat_trt$stage == 5, ]
# tempdat_sub6 = tempdat_trt[tempdat_trt$stage == 6, ]
# tempdat_sub7 = tempdat_trt[tempdat_trt$stage == 7, ]
# tempdat_sub8 = tempdat_trt[tempdat_trt$stage == 8, ]
# tempdat_sub9 = tempdat_trt[tempdat_trt$stage %in% c(9, 10), ]
# tempdat_sub10 = tempdat_trt[tempdat_trt$stage %in% c(11,12), ]
# tempdat_sub11 = tempdat_trt[tempdat_trt$stage %in% c(13,14,15), ]
# tempdat_sub12 = tempdat_trt[tempdat_trt$stage %in% c(16,17,18), ]
# tempdat_sub13 = tempdat_trt[tempdat_trt$stage > 18, ]
#tempdat_sub1 = tempdat_trt[tempdat_trt$id <= 250,]
#tempdat_sub2 = tempdat_trt[tempdat_trt$id > 250, ]

# c(5,7,8,14,18)

#================================================================================================#
#===================== STEP 1: TRT MODEL ========================================================#
#================================================================================================#

idx1 = c(3,5,7,8,14,18)
### model for denominator
mod_trt1 = fit_hal(X = tempdat_sub1[,idx1],
                   Y = tempdat_sub1$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt2 = fit_hal(X = tempdat_sub2[,idx1],
                   Y = tempdat_sub2$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt3 = fit_hal(X = tempdat_sub3[,idx1],
                   Y = tempdat_sub3$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt4 = fit_hal(X = tempdat_sub4[,idx1],
                   Y = tempdat_sub4$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt5 = fit_hal(X = tempdat_sub5[,idx1],
                   Y = tempdat_sub5$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt6 = fit_hal(X = tempdat_sub6[,idx1],
                   Y = tempdat_sub6$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt7 = fit_hal(X = tempdat_sub7[,idx1],
                   Y = tempdat_sub7$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt8 = fit_hal(X = tempdat_sub8[,idx1],
                   Y = tempdat_sub8$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt9 = fit_hal(X = tempdat_sub9[,idx1],
                   Y = tempdat_sub9$trt,
                   family = "binomial",
                   yolo = FALSE)
mod_trt10 = fit_hal(X = tempdat_sub10[,idx1],
                    Y = tempdat_sub10$trt,
                    family = "binomial",
                    yolo = FALSE)
mod_trt11 = fit_hal(X = tempdat_sub11[,idx1],
                    Y = tempdat_sub11$trt,
                    family = "binomial",
                    yolo = FALSE)
mod_trt12 = fit_hal(X = tempdat_sub12[,idx1],
                    Y = tempdat_sub12$trt,
                    family = "binomial",
                    yolo = FALSE)
mod_trt13 = fit_hal(X = tempdat_sub13[,idx1],
                    Y = tempdat_sub13$trt,
                    family = "binomial",
                    yolo = FALSE)
mod_trt14 = fit_hal(X = tempdat_sub14[,idx1],
                    Y = tempdat_sub14$trt,
                    family = "binomial",
                    yolo = FALSE)
mod_trt15 = fit_hal(X = tempdat_sub15[,idx1],
                    Y = tempdat_sub15$trt,
                    family = "binomial",
                    yolo = FALSE)
mod_trt16 = fit_hal(X = tempdat_sub16[,idx1],
                    Y = tempdat_sub16$trt,
                    family = "binomial",
                    yolo = FALSE)
### models for numerator
mod_trt1_num = fit_hal(X = tempdat_sub1[,idx1[1:4]],
                       Y = tempdat_sub1$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt2_num = fit_hal(X = tempdat_sub2[,idx1[1:4]],
                       Y = tempdat_sub2$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt3_num = fit_hal(X = tempdat_sub3[,idx1[1:4]],
                       Y = tempdat_sub3$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt4_num = fit_hal(X = tempdat_sub4[,idx1[1:4]],
                       Y = tempdat_sub4$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt5_num = fit_hal(X = tempdat_sub5[,idx1[1:4]],
                       Y = tempdat_sub5$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt6_num = fit_hal(X = tempdat_sub6[,idx1[1:4]],
                       Y = tempdat_sub6$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt7_num = fit_hal(X = tempdat_sub7[,idx1[1:4]],
                       Y = tempdat_sub7$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt8_num = fit_hal(X = tempdat_sub8[,idx1[1:4]],
                       Y = tempdat_sub8$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt9_num = fit_hal(X = tempdat_sub9[,idx1[1:4]],
                       Y = tempdat_sub9$trt,
                       family = "binomial",
                       yolo = FALSE)
mod_trt10_num = fit_hal(X = tempdat_sub10[,idx1[1:4]],
                        Y = tempdat_sub10$trt,
                        family = "binomial",
                        yolo = FALSE)
mod_trt11_num = fit_hal(X = tempdat_sub11[,idx1[1:4]],
                        Y = tempdat_sub11$trt,
                        family = "binomial",
                        yolo = FALSE)
mod_trt12_num = fit_hal(X = tempdat_sub12[,idx1[1:4]],
                        Y = tempdat_sub12$trt,
                        family = "binomial",
                        yolo = FALSE)
mod_trt13_num = fit_hal(X = tempdat_sub13[,idx1[1:4]],
                        Y = tempdat_sub13$trt,
                        family = "binomial",
                        yolo = FALSE)
mod_trt14_num = fit_hal(X = tempdat_sub14[,idx1[1:4]],
                        Y = tempdat_sub14$trt,
                        family = "binomial",
                        yolo = FALSE)
mod_trt15_num = fit_hal(X = tempdat_sub15[,idx1[1:4]],
                        Y = tempdat_sub15$trt,
                        family = "binomial",
                        yolo = FALSE)
mod_trt16_num = fit_hal(X = tempdat_sub16[,idx1[1:4]],
                        Y = tempdat_sub16$trt,
                        family = "binomial",
                        yolo = FALSE)

preds_trt = preds_trt_num = c()
preds_trt[tempdat_trt$id %in% c(1:30)] = predict(mod_trt1, new_data = tempdat_sub1[,idx1])
preds_trt[tempdat_trt$id %in% c(31:60)] = predict(mod_trt2, new_data = tempdat_sub2[,idx1])
preds_trt[tempdat_trt$id %in% c(61:90)] = predict(mod_trt3, new_data = tempdat_sub3[,idx1])
preds_trt[tempdat_trt$id %in% c(91:120)] = predict(mod_trt4, new_data = tempdat_sub4[,idx1])
preds_trt[tempdat_trt$id %in% c(121:150)] = predict(mod_trt5, new_data = tempdat_sub5[,idx1])
preds_trt[tempdat_trt$id %in% c(151:180)] = predict(mod_trt6, new_data = tempdat_sub6[,idx1])
preds_trt[tempdat_trt$id %in% c(181:210)] = predict(mod_trt7, new_data = tempdat_sub7[,idx1])
preds_trt[tempdat_trt$id %in% c(211:240)] = predict(mod_trt8, new_data = tempdat_sub8[,idx1])
preds_trt[tempdat_trt$id %in% c(241:270)] = predict(mod_trt9, new_data = tempdat_sub9[,idx1])
preds_trt[tempdat_trt$id %in% c(271:300)] = predict(mod_trt10, new_data = tempdat_sub10[,idx1])
preds_trt[tempdat_trt$id %in% c(301:330)] = predict(mod_trt11, new_data = tempdat_sub11[,idx1])
preds_trt[tempdat_trt$id %in% c(331:360)] = predict(mod_trt12, new_data = tempdat_sub12[,idx1])
preds_trt[tempdat_trt$id %in% c(361:390)] = predict(mod_trt13, new_data = tempdat_sub13[,idx1])
preds_trt[tempdat_trt$id %in% c(391:420)] = predict(mod_trt14, new_data = tempdat_sub14[,idx1])
preds_trt[tempdat_trt$id %in% c(421:450)] = predict(mod_trt15, new_data = tempdat_sub15[,idx1])
preds_trt[tempdat_trt$id %in% c(451:481)] = predict(mod_trt16, new_data = tempdat_sub16[,idx1])

preds_trt_num[tempdat_trt$id %in% c(1:30)] = predict(mod_trt1_num, new_data = tempdat_sub1[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(31:60)] = predict(mod_trt2_num, new_data = tempdat_sub2[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(61:90)] = predict(mod_trt3_num, new_data = tempdat_sub3[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(91:120)] = predict(mod_trt4_num, new_data = tempdat_sub4[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(121:150)] = predict(mod_trt5_num, new_data = tempdat_sub5[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(151:180)] = predict(mod_trt6_num, new_data = tempdat_sub6[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(181:210)] = predict(mod_trt7_num, new_data = tempdat_sub7[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(211:240)] = predict(mod_trt8_num, new_data = tempdat_sub8[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(241:270)] = predict(mod_trt9_num, new_data = tempdat_sub9[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(271:300)] = predict(mod_trt10_num, new_data = tempdat_sub10[,idx1][1:4])
preds_trt_num[tempdat_trt$id %in% c(301:330)] = predict(mod_trt11_num, new_data = tempdat_sub11[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(331:360)] = predict(mod_trt12_num, new_data = tempdat_sub12[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(361:390)] = predict(mod_trt13_num, new_data = tempdat_sub13[,idx1[1:4]])
preds_trt_num[tempdat_trt$id %in% c(391:420)] = predict(mod_trt14_num, new_data = tempdat_sub14[,idx1])
preds_trt_num[tempdat_trt$id %in% c(421:450)] = predict(mod_trt15_num, new_data = tempdat_sub15[,idx1])
preds_trt_num[tempdat_trt$id %in% c(451:481)] = predict(mod_trt16_num, new_data = tempdat_sub16[,idx1])
# preds_trt[tempdat_trt$stage == 1] =
#   predict(mod_trt1, new_data = tempdat_sub1[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage == 2] =
#   predict(mod_trt2, new_data = tempdat_sub2[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage == 3] =
#   predict(mod_trt3, new_data = tempdat_sub3[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage == 4] =
#   predict(mod_trt4, new_data = tempdat_sub4[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage == 5] =
#   predict(mod_trt5, new_data = tempdat_sub5[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage == 6] =
#   predict(mod_trt6, new_data = tempdat_sub6[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage == 7] =
#   predict(mod_trt7, new_data = tempdat_sub7[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage == 8] =
#   predict(mod_trt8, new_data = tempdat_sub8[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage %in% c(9, 10)] =
#   predict(mod_trt9, new_data = tempdat_sub9[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage %in% c(11,12)] =
#   predict(mod_trt10, new_data = tempdat_sub10[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage %in% c(13,14,15)] =
#   predict(mod_trt11, new_data = tempdat_sub11[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage %in% c(16,17,18)] =
#   predict(mod_trt12, new_data = tempdat_sub12[,c(5,7,8,14,18)])
# preds_trt[tempdat_trt$stage > 18] =
#   predict(mod_trt13, new_data = tempdat_sub13[,c(5,7,8,14,18)])

#preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
#preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))

tempdat$p.numerator_trt <- vector("numeric", nrow(tempdat))
tempdat$p.numerator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt_num[tempdat$trt[tempdat$selvar_trt == 1] == 0]
tempdat$p.numerator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt_num[tempdat$trt[tempdat$selvar_trt == 1] == 1]
tempdat$p.numerator_trt[tempdat$selvar_trt == 0] <- 1
tempdat$w.numerator_trt <- unlist(lapply(split(tempdat$p.numerator_trt, tempdat$id), function(x)cumprod(x)))

tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$id), function(x)cumprod(x)))

tempdat$ipw.weights_trt <- tempdat$w.numerator_trt/tempdat$w.denominator_trt

#=======================================================================================#
#============================ STEP 2; CENSORING MODEL ==================================#
#=======================================================================================#

tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$cen, tempdat$id),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))

tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
tempdat_sub1 = tempdat_cen[tempdat_cen$id %in% c(1:30), ]
tempdat_sub2 = tempdat_cen[tempdat_cen$id %in% c(31:60), ]
tempdat_sub3 = tempdat_cen[tempdat_cen$id %in% c(61:90), ]
tempdat_sub4 = tempdat_cen[tempdat_cen$id %in% c(91:120), ]
tempdat_sub5 = tempdat_cen[tempdat_cen$id %in% c(121:150), ]
tempdat_sub6 = tempdat_cen[tempdat_cen$id %in% c(151:180), ]
tempdat_sub7 = tempdat_cen[tempdat_cen$id %in% c(181:210), ]
tempdat_sub8 = tempdat_cen[tempdat_cen$id %in% c(211:240), ]
tempdat_sub9 = tempdat_cen[tempdat_cen$id %in% c(241:270), ]
tempdat_sub10 = tempdat_cen[tempdat_cen$id %in% c(271:300), ]
tempdat_sub11 = tempdat_cen[tempdat_cen$id %in% c(301:330), ]
tempdat_sub12 = tempdat_cen[tempdat_cen$id %in% c(331:360), ]
tempdat_sub13 = tempdat_cen[tempdat_cen$id %in% c(361:390), ]
tempdat_sub14 = tempdat_cen[tempdat_cen$id %in% c(391:420), ]
tempdat_sub15 = tempdat_cen[tempdat_cen$id %in% c(421:450), ]
tempdat_sub16 = tempdat_cen[tempdat_cen$id > 450, ]
#tempdat_sub1 = tempdat_cen[tempdat_cen$stage > 1 & tempdat_cen$stage <= 9, ]
#tempdat_sub2 = tempdat_cen[tempdat_cen$stage > 9, ]
#tempdat_sub1 = tempdat_cen[tempdat_cen$id <= 250, ]
#tempdat_sub2 = tempdat_cen[tempdat_cen$id > 250, ]

### c(2,3,5:9,13:21,23) # c(5,7,8,13,23)

### denominator models

idx2 = c(5,7,8,13,14,18,23)
mod_cen1 = fit_hal(X = tempdat_sub1[,idx2],
                   Y = tempdat_sub1$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen2 = fit_hal(X = tempdat_sub2[,idx2],
                   Y = tempdat_sub2$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen3 = fit_hal(X = tempdat_sub3[,idx2],
                   Y = tempdat_sub3$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen4 = fit_hal(X = tempdat_sub4[,idx2],
                   Y = tempdat_sub4$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen5 = fit_hal(X = tempdat_sub5[,idx2],
                   Y = tempdat_sub5$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen6 = fit_hal(X = tempdat_sub6[,idx2],
                   Y = tempdat_sub6$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen7 = fit_hal(X = tempdat_sub7[,idx2],
                   Y = tempdat_sub7$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen8 = fit_hal(X = tempdat_sub8[,idx2],
                   Y = tempdat_sub8$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen9 = fit_hal(X = tempdat_sub9[,idx2],
                   Y = tempdat_sub9$cen,
                   family = "binomial",
                   yolo = FALSE)
mod_cen10 = fit_hal(X = tempdat_sub10[,idx2],
                    Y = tempdat_sub10$cen,
                    family = "binomial",
                    yolo = FALSE)
mod_cen11 = fit_hal(X = tempdat_sub11[,idx2],
                    Y = tempdat_sub11$cen,
                    family = "binomial",
                    yolo = FALSE)
mod_cen12 = fit_hal(X = tempdat_sub12[,idx2],
                    Y = tempdat_sub12$cen,
                    family = "binomial",
                    yolo = FALSE)
mod_cen13 = fit_hal(X = tempdat_sub13[,idx2],
                    Y = tempdat_sub13$cen,
                    family = "binomial",
                    yolo = FALSE)
mod_cen14 = fit_hal(X = tempdat_sub14[,idx2],
                    Y = tempdat_sub14$cen,
                    family = "binomial",
                    yolo = FALSE)
mod_cen15 = fit_hal(X = tempdat_sub15[,idx2],
                    Y = tempdat_sub15$cen,
                    family = "binomial",
                    yolo = FALSE)
mod_cen16 = fit_hal(X = tempdat_sub16[,idx2],
                    Y = tempdat_sub16$cen,
                    family = "binomial",
                    yolo = FALSE)

### numerator models
mod_cen1_num = fit_hal(X = tempdat_sub1[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub1$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen2_num = fit_hal(X = tempdat_sub2[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub2$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen3_num = fit_hal(X = tempdat_sub3[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub3$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen4_num = fit_hal(X = tempdat_sub4[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub4$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen5_num = fit_hal(X = tempdat_sub5[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub5$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen6_num = fit_hal(X = tempdat_sub6[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub6$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen7_num = fit_hal(X = tempdat_sub7[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub7$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen8_num = fit_hal(X = tempdat_sub8[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub8$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen9_num = fit_hal(X = tempdat_sub9[,idx2[c(1,2,3,7)]],
                       Y = tempdat_sub9$cen,
                       family = "binomial",
                       yolo = FALSE)
mod_cen10_num = fit_hal(X = tempdat_sub10[,idx2[c(1,2,3,7)]],
                        Y = tempdat_sub10$cen,
                        family = "binomial",
                        yolo = FALSE)
mod_cen11_num = fit_hal(X = tempdat_sub11[,idx2[c(1,2,3,7)]],
                        Y = tempdat_sub11$cen,
                        family = "binomial",
                        yolo = FALSE)
mod_cen12_num = fit_hal(X = tempdat_sub12[,idx2[c(1,2,3,7)]],
                        Y = tempdat_sub12$cen,
                        family = "binomial",
                        yolo = FALSE)
mod_cen13_num = fit_hal(X = tempdat_sub13[,idx2[c(1,2,3,7)]],
                        Y = tempdat_sub13$cen,
                        family = "binomial",
                        yolo = FALSE)
mod_cen14_num = fit_hal(X = tempdat_sub14[,idx2[c(1,2,3,7)]],
                        Y = tempdat_sub14$cen,
                        family = "binomial",
                        yolo = FALSE)
mod_cen15_num = fit_hal(X = tempdat_sub15[,idx2[c(1,2,3,7)]],
                        Y = tempdat_sub15$cen,
                        family = "binomial",
                        yolo = FALSE)
mod_cen16_num = fit_hal(X = tempdat_sub16[,idx2[c(1,2,3,7)]],
                        Y = tempdat_sub16$cen,
                        family = "binomial",
                        yolo = FALSE)

preds_cen = preds_cen_num = c()
#preds_cen[tempdat_cen$stage > 1 & tempdat_cen$stage <= 9] =
#  predict(mod_cen1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
#preds_cen[tempdat_cen$stage > 9] =
#  predict(mod_cen2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))

preds_cen[tempdat_cen$id %in% c(1:30)] = predict(mod_cen1, new_data = tempdat_sub1[,idx2])
preds_cen[tempdat_cen$id %in% c(31:60)] = predict(mod_cen2, new_data = tempdat_sub2[,idx2])
preds_cen[tempdat_cen$id %in% c(61:90)] = predict(mod_cen3, new_data = tempdat_sub3[,idx2])
preds_cen[tempdat_cen$id %in% c(91:120)] = predict(mod_cen4, new_data = tempdat_sub4[,idx2])
preds_cen[tempdat_cen$id %in% c(121:150)] = predict(mod_cen5, new_data = tempdat_sub5[,idx2])
preds_cen[tempdat_cen$id %in% c(151:180)] = predict(mod_cen6, new_data = tempdat_sub6[,idx2])
preds_cen[tempdat_cen$id %in% c(181:210)] = predict(mod_cen7, new_data = tempdat_sub7[,idx2])
preds_cen[tempdat_cen$id %in% c(211:240)] = predict(mod_cen8, new_data = tempdat_sub8[,idx2])
preds_cen[tempdat_cen$id %in% c(241:270)] = predict(mod_cen9, new_data = tempdat_sub9[,idx2])
preds_cen[tempdat_cen$id %in% c(271:300)] = predict(mod_cen10, new_data = tempdat_sub10[,idx2])
preds_cen[tempdat_cen$id %in% c(301:330)] = predict(mod_cen11, new_data = tempdat_sub11[,idx2])
preds_cen[tempdat_cen$id %in% c(331:360)] = predict(mod_cen12, new_data = tempdat_sub12[,idx2])
preds_cen[tempdat_cen$id %in% c(361:390)] = predict(mod_cen13, new_data = tempdat_sub13[,idx2])
preds_cen[tempdat_cen$id %in% c(391:420)] = predict(mod_cen14, new_data = tempdat_sub14[,idx2])
preds_cen[tempdat_cen$id %in% c(421:450)] = predict(mod_cen15, new_data = tempdat_sub15[,idx2])
preds_cen[tempdat_cen$id %in% c(451:481)] = predict(mod_cen16, new_data = tempdat_sub16[,idx2])

preds_cen_num[tempdat_cen$id %in% c(1:30)] = predict(mod_cen1_num, new_data = tempdat_sub1[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(31:60)] = predict(mod_cen2_num, new_data = tempdat_sub2[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(61:90)] = predict(mod_cen3_num, new_data = tempdat_sub3[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(91:120)] = predict(mod_cen4_num, new_data = tempdat_sub4[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(121:150)] = predict(mod_cen5_num, new_data = tempdat_sub5[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(151:180)] = predict(mod_cen6_num, new_data = tempdat_sub6[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(181:210)] = predict(mod_cen7_num, new_data = tempdat_sub7[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(211:240)] = predict(mod_cen8_num, new_data = tempdat_sub8[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(241:270)] = predict(mod_cen9_num, new_data = tempdat_sub9[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(271:300)] = predict(mod_cen10_num, new_data = tempdat_sub10[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(301:330)] = predict(mod_cen11_num, new_data = tempdat_sub11[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(331:360)] = predict(mod_cen12_num, new_data = tempdat_sub12[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(361:390)] = predict(mod_cen13_num, new_data = tempdat_sub13[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(391:420)] = predict(mod_cen14_num, new_data = tempdat_sub14[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(421:450)] = predict(mod_cen15_num, new_data = tempdat_sub15[,idx2[c(1,2,3,7)]])
preds_cen_num[tempdat_cen$id %in% c(451:481)] = predict(mod_cen16_num, new_data = tempdat_sub16[,idx2[c(1,2,3,7)]])


tempdat$p.numerator_cen <- vector("numeric", nrow(tempdat))
tempdat$p.numerator_cen[tempdat$cen == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen_num[tempdat$cen[tempdat$selvar_cen == 1] == 0]
tempdat$p.numerator_cen[tempdat$cen == 1 & tempdat$selvar_cen == 1] <- preds_cen_num[tempdat$cen[tempdat$selvar_cen == 1] == 1]
tempdat$p.numerator_cen[tempdat$selvar_cen == 0] <- 1
tempdat$w.numerator_cen <- unlist(lapply(split(tempdat$p.numerator_cen, tempdat$id), function(x)cumprod(x)))

tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
tempdat$p.denominator_cen[tempdat$cen == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$cen[tempdat$selvar_cen == 1] == 0]
tempdat$p.denominator_cen[tempdat$cen == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$cen[tempdat$selvar_cen == 1] == 1]
tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$id), function(x)cumprod(x)))
tempdat$ipw.weights_cen <- tempdat$w.numerator_cen/tempdat$w.denominator_cen

### assign the weight function
#clean_dat = qol_dat %>% select(ID = id, stage = time, x1 = maxfvc, x2 = bmi_pct_change_base, qol = QOL, trt, outcome, censor = cen, start = stage)
clean_dat = qol_dat %>% select(ID = id, stage = time, x1 = maxfvc, x2 = bmi_pct_change_base,
                               x3 = als_bulbar, x1_square = maxfvc_square, x1_cubic = maxfvc_cubic,
                               qol = QOL, trt, outcome, censor = cen, start = stage)
clean_dat$weight = 1
clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
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

#=====================================================================================================#
#=========== the example when bmi, fvc and bulbar are considered =====================================#
#=========== this is the case suggested by emory =====================================================#
#=====================================================================================================#

temp5 <- genoud(fn=indicator_opt1_als, index = c(3, 4, 5, 6, 7), nvars=6, default.domains=1, starting.values=rep(0, 6),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod5 = loss_indicator2_als(n = 481, eta = temp1$par, dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5, 6, 7))

temp6 <- genoud(fn=smooth_opt2_als, index = c(3, 4, 5, 6, 7), nvars=6, default.domains=1, starting.values=rep(0, 6),
                max=TRUE, print.level=1, BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130)
mod6 = loss_indicator2_als(n = 481, eta = temp2$par, dat = qol_dat, clean_dat, upp = upp, time, target = "rmst", est_var = T, index = c(3, 4, 5, 6, 7))


### save the optimation results
real_hal = list(temp1 = temp1, temp2 = temp2, temp3 = temp3, temp4 = temp4, temp5 = temp5, temp6 = temp6)
setwd("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL/code_github")
save(real_logit, file = "real_hal.rda")
