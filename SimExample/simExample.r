#===========================================================================================#
#========================== simulation for our methods =====================================#
#===========================================================================================#
#================ this simulation is pretty computationally expensive ======================#
#================= and we use clusters to get the results ==================================#
#===========================================================================================#

#========================================================#
#========== some basic set up ===========================#
#========================================================#
source("~/Dropbox/Hsun/Research project/Causal Inference/ALS case study/QAL/code_github/SimExample/support_sim.r")

### load necessary packages
library(dplyr)
library(rgenoud)

### global parameter 
# sample size
n = 250
# study duration
K = 60
# gap between 2 landmarks
interval = 10
# number of stages
stage = K/interval
# the upper limit of the integral for RQAL
upp = 27 # 36
# the true optimal regime
eta = eta_opt = eta_opt = c(-0.6, -1, +1)
  
#c(1, -1, -1)
time = c(1:upp)


seed = 1
set.seed(100+seed)

### data generation and transformation
dat = data_gen_ch3(n,K, interval, eta_opt = eta_opt)


#=======================================================================#
#================ step 1: the eta_opt estimation =======================#
#=======================================================================#

mod_logit = optimal_sim1(n = 250, dat = dat, nuisance = "logit", est = "IPW")
mod_hal = optimal_sim1(n = 250, dat = dat, nuisance = "hal", est = "IPW")
mod_rf = optimal_sim1(n = 250, dat = dat, nuisance = "rf", est = "IPW")




#=======================================================================#
#=== step 2: caiculate the expected value under estimated regimes ======#
#=======================================================================#

# we approximate the true expected value by stochastic integration; this step would be pretty time consuming
po1 = data_gen_true_ch3(100000,K, interval, eta = mod_logit$eta_opt, eta_opt = eta_opt)
po1 = po1[po1$x1>0,]

# the mis-classification rate
mr1 = po1 %>% group_by(ID) %>% summarise(clas = length(ID) == sum(follow))
qal_time1 = po1 %>% group_by(ID) %>% summarise(outcome = sum(outcome * qol))
a = pmin(qal_time1$outcome, upp)
base = c(0, sort(unique(a)))
pp=c()
for(i in 1:(length(base)-1)){
  pp[i]= mean(qal_time1$outcome > base[i])
}
integ = diff(base)

# the expected RQAL 
rmst1 = sum(integ * pp)
qal_prob = c()
for(tt in 1:upp){
  qal_prob[tt] = mean(qal_time1$outcome > tt)
}
