### clean the results for new mis-simulation
path <- paste0("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL")
setwd(path)
library(dplyr)
### the sd estimation
load("PO_60.rda")
load("PO_100.rda")

### calculate the RQAL under two cases
upp = 26
qal_time = PO_60 %>% group_by(ID) %>% summarise(outcome = sum(outcome * qol))
a = pmin(qal_time$outcome, upp)
base = c(0, sort(unique(a)))
pp=c()
for(i in 1:(length(base)-1)){
  pp[i]= mean(qal_time$outcome > base[i])
}
integ = diff(base)
rmst_60 = sum(integ * pp)

upp = 36
qal_time = PO_100 %>% group_by(ID) %>% summarise(outcome = sum(outcome * qol))
a = pmin(qal_time$outcome, upp)
base = c(0, sort(unique(a)))
pp=c()
for(i in 1:(length(base)-1)){
  pp[i]= mean(qal_time$outcome > base[i])
}
integ = diff(base)
rmst_100 = sum(integ * pp)

# ### check the data dist
# library(dplyr)
# #library(ipw)
# cen = dat %>% group_by(ID) %>% summarise(cen = sum(censor))
# cat("the censoring rate is", mean(cen$cen))
# obs_surv = dat %>% group_by(ID) %>% summarise(outcome = sum(outcome))
# cat("the mean observed survival time is", mean(obs_surv$outcome))
# quantile(obs_surv$outcome)
# hist(obs_surv$outcome)
# obs_qal = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
# cat("the mean observed qal is", mean(obs_qal$qal))
# a = dat %>% group_by(ID) %>% summarise(trt = sum(trt) > 0)
# cat("the mean observed treatment rate is", mean(a$trt))
# cat("the mean rate of following optimal treatment is", mean(dat$follow[dat$x1 > 0 & dat$stage > 1]))
# 
# dat0 = dat[dat$x1>0,]
# mr = dat0 %>% group_by(ID) %>% summarise(clas = length(ID) == sum(follow))
# mean(mr$clas)


#============================================================================================================#
# ======================================================LR ==================================================#
#============================================================================================================#

#==========================================================================================#
#======================== case 1:  K = 60, interval = 10, n = 250, indicator ==============#
#==========================================================================================#

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n250qaloptupp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_60_250 =  rqal2_60_250 = se_60_250 = c()
coef1_60_250 = coef1_norm_60_250 = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_250 = mr_po_60_250 = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_60_250[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_60_250[m] = result$mod$se
    
    coef1_60_250[m,] = result$mod$eta_opt
    coef1_norm_60_250[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_250[m] = result$rmst_po
    mr_po_60_250[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n500qaloptupp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 26
iter <- 500
effect <- c()

rqal1_60_500 =  rqal2_60_500 = se_60_500 = c()
coef1_60_500 = coef1_norm_60_500 = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_500 = mr_po_60_500 = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_60_500[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_60_500[m] = result$mod$se
    
    coef1_60_500[m,] = result$mod$eta_opt
    coef1_norm_60_500[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_500[m] = result$rmst_po
    mr_po_60_500[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n250qaloptupp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_250 =  rqal2_100_250= se_100_250 = c()
coef1_100_250 = coef1_norm_100_250 = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_250 = mr_po_100_250 = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_100_250[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_100_250[m] = result$mod$se
    
    coef1_100_250[m,] = result$mod$eta_opt
    coef1_norm_100_250[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_250[m] = result$rmst_po
    mr_po_100_250[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n500qaloptupp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_500 =  rqal2_100_500 = se_100_500 = c()
coef1_100_500 = coef1_norm_100_500 = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_500 = mr_po_100_500 = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_100_500[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_100_500[m] = result$mod$se
    
    coef1_100_500[m,] = result$mod$eta_opt
    coef1_norm_100_500[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_500[m] = result$rmst_po
    mr_po_100_500[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}


#==========================================================================================#
#======================== case 1:  K = 60, interval = 10, n = 250, indicator ==============#
#==========================================================================================#

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n250qalopt_bc_beta_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_60_250_s =  rqal2_60_250_s = se_60_250_s = c()
coef1_60_250_s = coef1_norm_60_250_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_250_s = mr_po_60_250_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_60_250_s[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_60_250_s[m] = result$mod$se
    
    coef1_60_250_s[m,] = result$mod$eta_opt
    coef1_norm_60_250_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_250_s[m] = result$rmst_po
    mr_po_60_250_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n500qalopt_bc_beta_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 26
iter <- 500
effect <- c()

rqal1_60_500_s =  rqal2_60_500_s = se_60_500_s = c()
coef1_60_500_s = coef1_norm_60_500_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_500_s = mr_po_60_500_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_60_500_s[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_60_500_s[m] = result$mod$se
    
    coef1_60_500_s[m,] = result$mod$eta_opt
    coef1_norm_60_500_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_500_s[m] = result$rmst_po
    mr_po_60_500_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n250qalopt_bc_beta_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_250_s =  rqal2_100_250_s= se_100_250_s = c()
coef1_100_250_s = coef1_norm_100_250_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_250_s = mr_po_100_250_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_100_250_s[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_100_250_s[m] = result$mod$se
    
    coef1_100_250_s[m,] = result$mod$eta_opt
    coef1_norm_100_250_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_250_s[m] = result$rmst_po
    mr_po_100_250_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/logit_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n500qalopt_bc_beta_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_500_s =  rqal2_100_500_s = se_100_500_s = c()
coef1_100_500_s = coef1_norm_100_500_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_500_s = mr_po_100_500_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250[m,] = result$rqal0
    rqal1_100_500_s[m] = result$mod$est
    #rqal2_60_250[m] = result$mod$est2
    se_100_500_s[m] = result$mod$se
    
    coef1_100_500_s[m,] = result$mod$eta_opt
    coef1_norm_100_500_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_500_s[m] = result$rmst_po
    mr_po_100_500_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}
# mean(rqal1_60_250, na.rm = T)
# sd(rqal1_60_250, na.rm = T)
# mean(se_60_250, na.rm = T)
# 
# mean(rqal2_60_250, na.rm = T)
# sd(rqal2_60_250, na.rm = T)
# mean(se_60_250, na.rm = T)
# 
# apply(coef1_norm_60_250,2,mean, na.rm  = T)
# apply(coef1_norm_60_250,2,sd, na.rm  = T)
# apply(coef1_60_250,2,mean, na.rm  = T)
# apply(coef1_60_250,2,sd, na.rm  = T)
# mean(coef1_60_250[,1]/coef1_60_250[,2], na.rm = T)
# mean(coef1_60_250[,1]/coef1_60_250[,3], na.rm = T)
# mean(coef1_60_250[,2]/coef1_60_250[,3], na.rm = T)
# sd(coef1_60_250[,1]/coef1_60_250[,2], na.rm = T)
# sd(coef1_60_250[,1]/coef1_60_250[,3], na.rm = T)
# sd(coef1_60_250[,2]/coef1_60_250[,3], na.rm = T)
# 
# cp = mean(rqal1_60_250 - 1.96 * se_60_250 <= rmst_60 & rqal1_60_250 + 1.96 * se_60_250 >= rmst_60, na.rm = T)
# mean(rmst_po_60_250, na.rm = T)
# sd(rmst_po_60_250, na.rm = T)
# mean(1-mr_po_60_250, na.rm = T)
# sd(1-mr_po_60_250, na.rm = T)
# 


#=====================================================================================================================#
#==================================================== HAL ============================================================#
#=====================================================================================================================#

#==========================================================================================#
#======================== case 1:  K = 60, interval = 10, n = 250, indicator ==============#
#==========================================================================================#

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n250qalopt_hal_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_60_250_hal =  rqal2_60_250_hal = se_60_250_hal = c()
coef1_60_250_hal = coef1_norm_60_250_hal = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_250_hal = mr_po_60_250_hal = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_60_250_hal[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_60_250_hal[m] = result$mod$se
    
    coef1_60_250_hal[m,] = result$mod$eta_opt
    coef1_norm_60_250_hal[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_250_hal[m] = result$rmst_po
    mr_po_60_250_hal[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n500qalopt_hal_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 26
iter <- 500
effect <- c()

rqal1_60_500_hal =  rqal2_60_500_hal = se_60_500_hal = c()
coef1_60_500_hal = coef1_norm_60_500_hal = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_500_hal = mr_po_60_500_hal = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_60_500_hal[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_60_500_hal[m] = result$mod$se
    
    coef1_60_500_hal[m,] = result$mod$eta_opt
    coef1_norm_60_500_hal[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_500_hal[m] = result$rmst_po
    mr_po_60_500_hal[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n250qalopt_hal_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_250_hal =  rqal2_100_250_hal= se_100_250_hal = c()
coef1_100_250_hal = coef1_norm_100_250_hal = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_250_hal = mr_po_100_250_hal = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_100_250_hal[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_100_250_hal[m] = result$mod$se
    
    coef1_100_250_hal[m,] = result$mod$eta_opt
    coef1_norm_100_250_hal[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_250_hal[m] = result$rmst_po
    mr_po_100_250_hal[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n500qalopt_hal_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_500_hal =  rqal2_100_500_hal = se_100_500_hal = c()
coef1_100_500_hal = coef1_norm_100_500_hal = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_500_hal = mr_po_100_500_hal = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_100_500_hal[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_100_500_hal[m] = result$mod$se
    
    coef1_100_500_hal[m,] = result$mod$eta_opt
    coef1_norm_100_500_hal[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_500_hal[m] = result$rmst_po
    mr_po_100_500_hal[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}


#==========================================================================================#
#======================== case 1:  K = 60, interval = 10, n = 250, indicator ==============#
#==========================================================================================#

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n250qalopt_bc_hal_beta_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_60_250_hal_s =  rqal2_60_250_hal_s = se_60_250_hal_s = c()
coef1_60_250_hal_s = coef1_norm_60_250_hal_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_250_hal_s = mr_po_60_250_hal_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_60_250_hal_s[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_60_250_hal_s[m] = result$mod$se
    
    coef1_60_250_hal_s[m,] = result$mod$eta_opt
    coef1_norm_60_250_hal_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_250_hal_s[m] = result$rmst_po
    mr_po_60_250_hal_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n500qalopt_bc_hal_beta_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 26
iter <- 500
effect <- c()

rqal1_60_500_hal_s =  rqal2_60_500_hal_s = se_60_500_hal_s = c()
coef1_60_500_hal_s = coef1_norm_60_500_hal_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_500_hal_s = mr_po_60_500_hal_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_60_500_hal_s[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_60_500_hal_s[m] = result$mod$se
    
    coef1_60_500_hal_s[m,] = result$mod$eta_opt
    coef1_norm_60_500_hal_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_500_hal_s[m] = result$rmst_po
    mr_po_60_500_hal_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n250qalopt_bc_hal_beta_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_250_hal_s =  rqal2_100_250_hal_s= se_100_250_hal_s = c()
coef1_100_250_hal_s = coef1_norm_100_250_hal_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_250_hal_s = mr_po_100_250_hal_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_100_250_hal_s[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_100_250_hal_s[m] = result$mod$se
    
    coef1_100_250_hal_s[m,] = result$mod$eta_opt
    coef1_norm_100_250_hal_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_250_hal_s[m] = result$rmst_po
    mr_po_100_250_hal_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/hal_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n500qalopt_bc_hal_beta_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_500_hal_s =  rqal2_100_500_hal_s = se_100_500_hal_s = c()
coef1_100_500_hal_s = coef1_norm_100_500_hal_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_500_hal_s = mr_po_100_500_hal_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_hal[m,] = result$rqal0
    rqal1_100_500_hal_s[m] = result$mod$est
    #rqal2_60_250_hal[m] = result$mod$est2
    se_100_500_hal_s[m] = result$mod$se
    
    coef1_100_500_hal_s[m,] = result$mod$eta_opt
    coef1_norm_100_500_hal_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_500_hal_s[m] = result$rmst_po
    mr_po_100_500_hal_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}
# mean(rqal1_60_250_hal, na.rm = T)




#=====================================================================================================================#
#==================================================== rf ============================================================#
#=====================================================================================================================#

#==========================================================================================#
#======================== case 1:  K = 60, interval = 10, n = 250, indicator ==============#
#==========================================================================================#

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n250qalopt_rf_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_60_250_rf =  rqal2_60_250_rf = se_60_250_rf = c()
coef1_60_250_rf = coef1_norm_60_250_rf = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_250_rf = mr_po_60_250_rf = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_60_250_rf[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_60_250_rf[m] = result$mod$se
    
    coef1_60_250_rf[m,] = result$mod$eta_opt
    coef1_norm_60_250_rf[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_250_rf[m] = result$rmst_po
    mr_po_60_250_rf[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n500qalopt_rf_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 26
iter <- 500
effect <- c()

rqal1_60_500_rf =  rqal2_60_500_rf = se_60_500_rf = c()
coef1_60_500_rf = coef1_norm_60_500_rf = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_500_rf = mr_po_60_500_rf = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_60_500_rf[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_60_500_rf[m] = result$mod$se
    
    coef1_60_500_rf[m,] = result$mod$eta_opt
    coef1_norm_60_500_rf[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_500_rf[m] = result$rmst_po
    mr_po_60_500_rf[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n250qalopt_rf_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_250_rf =  rqal2_100_250_rf= se_100_250_rf = c()
coef1_100_250_rf = coef1_norm_100_250_rf = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_250_rf = mr_po_100_250_rf = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_100_250_rf[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_100_250_rf[m] = result$mod$se
    
    coef1_100_250_rf[m,] = result$mod$eta_opt
    coef1_norm_100_250_rf[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_250_rf[m] = result$rmst_po
    mr_po_100_250_rf[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n500qalopt_rf_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_500_rf =  rqal2_100_500_rf = se_100_500_rf = c()
coef1_100_500_rf = coef1_norm_100_500_rf = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_500_rf = mr_po_100_500_rf = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_100_500_rf[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_100_500_rf[m] = result$mod$se
    
    coef1_100_500_rf[m,] = result$mod$eta_opt
    coef1_norm_100_500_rf[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_500_rf[m] = result$rmst_po
    mr_po_100_500_rf[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}


#==========================================================================================#
#======================== case 1:  K = 60, interval = 10, n = 250, indicator ==============#
#==========================================================================================#

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n250qalopt_bc_rf_beta_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_60_250_rf_s =  rqal2_60_250_rf_s = se_60_250_rf_s = c()
coef1_60_250_rf_s = coef1_norm_60_250_rf_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_250_rf_s = mr_po_60_250_rf_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_60_250_rf_s[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_60_250_rf_s[m] = result$mod$se
    
    coef1_60_250_rf_s[m,] = result$mod$eta_opt
    coef1_norm_60_250_rf_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_250_rf_s[m] = result$rmst_po
    mr_po_60_250_rf_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K60n500qalopt_bc_rf_beta_upp26"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 26
iter <- 500
effect <- c()

rqal1_60_500_rf_s =  rqal2_60_500_rf_s = se_60_500_rf_s = c()
coef1_60_500_rf_s = coef1_norm_60_500_rf_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_60_500_rf_s = mr_po_60_500_rf_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_60_500_rf_s[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_60_500_rf_s[m] = result$mod$se
    
    coef1_60_500_rf_s[m,] = result$mod$eta_opt
    coef1_norm_60_500_rf_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_60_500_rf_s[m] = result$rmst_po
    mr_po_60_500_rf_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n250qalopt_bc_rf_beta_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_250_rf_s =  rqal2_100_250_rf_s= se_100_250_rf_s = c()
coef1_100_250_rf_s = coef1_norm_100_250_rf_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_250_rf_s = mr_po_100_250_rf_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_100_250_rf_s[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_100_250_rf_s[m] = result$mod$se
    
    coef1_100_250_rf_s[m,] = result$mod$eta_opt
    coef1_norm_100_250_rf_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_250_rf_s[m] = result$rmst_po
    mr_po_100_250_rf_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}

################ the file to read simulated results from bluehive ######
############### in case 2, we include the potential outcome       ######
path <- paste0("F:/QAL/newsim2019/rf_s")
############ load the true potential outcome 
setwd(paste0(path, "/output_K100n500qalopt_bc_rf_beta_upp36"))

alloutputs <- paste0("QAL",1:500,".rda")

upp = 36
iter <- 500
effect <- c()

rqal1_100_500_rf_s =  rqal2_100_500_rf_s = se_100_500_rf_s = c()
coef1_100_500_rf_s = coef1_norm_100_500_rf_s = matrix(NA, ncol = 3, nrow = 500)

#qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)

rmst_po_100_500_rf_s = mr_po_100_500_rf_s = c()

for(m in 1:iter){
  if(file.exists(alloutputs[m])){
    effect[m] <- 1
    load(alloutputs[m])
    
    #rqal0_60_250_rf[m,] = result$rqal0
    rqal1_100_500_rf_s[m] = result$mod$est
    #rqal2_60_250_rf[m] = result$mod$est2
    se_100_500_rf_s[m] = result$mod$se
    
    coef1_100_500_rf_s[m,] = result$mod$eta_opt
    coef1_norm_100_500_rf_s[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
    
    rmst_po_100_500_rf_s[m] = result$rmst_po
    mr_po_100_500_rf_s[m] = mean(result$mr_po$clas)
  }
  else effect[m] <- 0
}


#==================================================================#
#=========== 01/05/2020 added plots for supp ======================#
#==================================================================#
path <- paste0("~/Desktop")
pdf(paste0(path, "/sim_s1_k6.pdf"), width = 15, height  = 15)
par(mfrow = c(3, 1),     # 2x2 layout
    oma = c(3, 5, 6, 0), # two rows of text at the outer left and bottom margin
    mar = c(1,1,1,1))
# potential outcome plots
po_result_60 = cbind(rmst_po_60_250, rmst_po_60_250_s, rmst_po_60_500, rmst_po_60_500_s,
                     rmst_po_60_250_hal, rmst_po_60_250_hal_s, rmst_po_60_500_hal, rmst_po_60_500_hal_s,
                     rmst_po_60_250_rf, rmst_po_60_250_rf_s, rmst_po_60_500_rf, rmst_po_60_500_rf_s)
              
colnames(po_result_60) = c("LR(250)", "BC-LR(250)", "LR(500)", "BC-LR(500)",
                           "HAL(250)", "BC-HAL(250)", "HAL(500)", "BC-HAL(500)",
                           "RF(250)", "BC-RF(250)", "RF(500)", "BC-RF(500)")
boxplot(po_result_60, ylab = "RQAL", cex.axis = 1.5, frame.plot = F)
box(bty="l")
abline(a = rmst_60, b = 0, col = 1, lwd = 2, lty = 2)

# est plots
est_result_60 = cbind(rqal1_60_250, rqal1_60_250_s, rqal1_60_500, rqal1_60_500_s,
                      rqal1_60_250_hal, rqal1_60_250_hal_s, rqal1_60_500_hal, rqal1_60_500_hal_s,
                      rqal1_60_250_rf, rqal1_60_250_rf_s, rqal1_60_500_rf, rqal1_60_500_rf_s)
colnames(est_result_60) = c("LR(250)", "BC-LR(250)", "LR(500)", "BC-LR(500)",
                            "HAL(250)", "BC-HAL(250)", "HAL(500)", "BC-HAL(500)",
                            "RF(250)", "BC-RF(250)", "RF(500)", "BC-RF(500)")
boxplot(est_result_60, ylab = "RQAL estimation", cex.axis = 1.5, frame.plot = F)
box(bty="l")
abline(a = rmst_60, b = 0, col = 1, lwd = 2, lty = 2)

### misclassification rate
mis_rate_60 = 1 - cbind(mr_po_60_250, mr_po_60_250_s, mr_po_60_500, mr_po_60_500_s,
                        mr_po_60_250_hal, mr_po_60_250_hal_s, mr_po_60_500_hal, mr_po_60_500_hal_s,
                        mr_po_60_250_rf, mr_po_60_250_rf_s, mr_po_60_500_rf, mr_po_60_500_rf_s)
colnames(mis_rate_60) = c("LR(250)", "BC-LR(250)", "LR(500)", "BC-LR(500)",
                          "HAL(250)", "BC-HAL(250)", "HAL(500)", "BC-HAL(500)",
                          "RF(250)", "BC-RF(250)", "RF(500)", "BC-RF(500)")
boxplot(mis_rate_60, ylab = "Mis-classification rate", cex.axis = 1.5, frame.plot = F)
box(bty="l")

#title(ylab = "Coverage Rate", outer = TRUE, line = 2, cex.lab = 1.4, adj = 0.07, font.lab = 2)
title(ylab = "Mis-classification Rate", outer = TRUE, line = 1.5, cex.lab = 2.5, adj = 0.07, font.lab = 1.5)
title(ylab = expression(hat(R)(hat(eta)[opt])), outer = TRUE, line = 1.2, cex.lab = 2.5, adj = 0.5,  font.lab = 2)
title(ylab = expression(R(hat(eta)[opt])), outer = TRUE, line = 1.2, cex.lab = 2.5, adj = 0.85,  font.lab = 2)
title(main = "K = 6", outer = TRUE, line = 1, cex.main = 3, adj = 0.5)

dev.off()


path <- paste0("~/Desktop")
pdf(paste0(path, "/sim_s1_k25.pdf"), width = 15, height  = 15)
par(mfrow = c(3, 1),     # 2x2 layout
    oma = c(3, 5, 6, 0), # two rows of text at the outer left and bottom margin
    mar = c(1,1,1,1))
# potential outcome plots
po_result_100 = cbind(rmst_po_100_250, rmst_po_100_250_s, rmst_po_100_500, rmst_po_100_500_s,
                      rmst_po_100_250_hal, rmst_po_100_250_hal_s, rmst_po_100_500_hal, rmst_po_100_500_hal_s,
                      rmst_po_100_250_rf, rmst_po_100_250_rf_s, rmst_po_100_500_rf, rmst_po_100_500_rf_s)

colnames(po_result_100) = c("LR(250)", "BC-LR(250)", "LR(500)", "BC-LR(500)",
                            "HAL(250)", "BC-HAL(250)", "HAL(500)", "BC-HAL(500)",
                            "RF(250)", "BC-RF(250)", "RF(500)", "BC-RF(500)")
boxplot(po_result_100, ylab = "RQAL", cex.axis = 1.5, frame.plot = F)
box(bty="l")
abline(a = rmst_100, b = 0, col = 1, lwd = 2, lty = 2)


# est plots
est_result_100 = cbind(rqal1_100_250, rqal1_100_250_s, rqal1_100_500, rqal1_100_500_s,
                       rqal1_100_250_hal, rqal1_100_250_hal_s, rqal1_100_500_hal, rqal1_100_500_hal_s,
                       rqal1_100_250_rf, rqal1_100_250_rf_s, rqal1_100_500_rf, rqal1_100_500_rf_s)
colnames(est_result_100) = c("LR(250)", "BC-LR(250)", "LR(500)", "BC-LR(500)",
                             "HAL(250)", "BC-HAL(250)", "HAL(500)", "BC-HAL(500)",
                             "RF(250)", "BC-RF(250)", "RF(500)", "BC-RF(500)")
boxplot(est_result_100, ylab = "RQAL estimation", cex.axis = 1.5, frame.plot = F)
box(bty="l")
abline(a = rmst_100, b = 0, col = 1, lwd = 2, lty = 2)

### misclassification rate
mis_rate_100 = 1 - cbind(mr_po_100_250, mr_po_100_250_s, mr_po_100_500, mr_po_100_500_s,
                         mr_po_100_250_hal, mr_po_100_250_hal_s, mr_po_100_500_hal, mr_po_100_500_hal_s,
                         mr_po_100_250_rf, mr_po_100_250_rf_s, mr_po_100_500_rf, mr_po_100_500_rf_s)
colnames(mis_rate_100) = c("LR(250)", "BC-LR(250)", "LR(500)", "BC-LR(500)",
                           "HAL(250)", "BC-HAL(250)", "HAL(500)", "BC-HAL(500)",
                           "RF(250)", "BC-RF(250)", "RF(500)", "BC-RF(500)")
boxplot(mis_rate_100, ylab = "Mis-classification rate", cex.axis = 1.5, frame.plot = F)
box(bty="l")

title(ylab = "Mis-classification Rate", outer = TRUE, line = 1.5, cex.lab = 2.5, adj = 0.07, font.lab = 1.5)
title(ylab = expression(hat(R)(hat(eta)[opt])), outer = TRUE, line = 1.2, cex.lab = 2.5, adj = 0.5,  font.lab = 2)
title(ylab = expression(R(hat(eta)[opt])), outer = TRUE, line = 1.2, cex.lab = 2.5, adj = 0.85,  font.lab = 2)
title(main = "K = 25", outer = TRUE, line = 1, cex.main = 3, adj = 0.5)

dev.off()
