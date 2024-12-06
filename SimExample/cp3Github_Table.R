#========================== table in chp 3 ===============================#

alloutputs <- paste0("QAL",5001:5500,".rda")

summ_tab <- function(stage=60, iter=500){
  upp = 30
  iter <- iter
  effect <- c()
  rqal0_60_250 = rqal1_60_250 =  matrix(NA, ncol = 1, nrow = iter)
  coef1_60_250 = coef1_norm_60_250 = matrix(NA, ncol = 3, nrow = iter)
  se_est = matrix(NA, ncol = 1, nrow = iter)
  #qal_surv = qal_surv2 = qal_surv_var = qal_surv_var2 = matrix(NA, nrow = 500, ncol = 5)
  
  rmst_po_60_250 = mr_po_60_250 = c()
  rr = matrix(NA, ncol = 4, nrow = iter)
  
  for(m in 1:iter){
    if(file.exists(alloutputs[m])){
      effect[m] <- 1
      load(alloutputs[m])
      rqal0_60_250[m,] = result$mod$est
      
      coef1_60_250[m,] = result$mod$eta_opt
      coef1_norm_60_250[m,] = result$mod$eta_opt/sqrt(sum(result$mod$eta_opt^2))
      se_est[m]  = result$mod$se
      
      rmst_po_60_250[m] = result$rmst_po
      mr_po_60_250[m] = mean(result$mr_po$clas)
      
      #rr[m,] = result$mod$count
    }
    else effect[m] <- 0
  }
  
  low = 0.02;high = 0.98
  ratio1 = -coef1_60_250[,1]/coef1_60_250[,2]
  q1 = quantile(ratio1, probs = c(low, high), na.rm = T)
  ratio1 = ratio1[(ratio1 >= q1[1])&(ratio1 <= q1[2])] # delete outliers
  ratio2 = coef1_60_250[,2]/coef1_60_250[,3]
  q2 = quantile(ratio2, probs = c(low, high), na.rm = T)
  ratio2 = ratio2[(ratio2 >= q2[1])&(ratio2 <= q2[2])] 
  ratio3 = -coef1_60_250[,3]/coef1_60_250[,1]
  q3 = quantile(ratio3, probs = c(low, high), na.rm = T)
  ratio3 = ratio3[(ratio3 >= q3[1])&(ratio3 <= q3[2])] 
  ratio4 = -coef1_60_250[,1]/coef1_60_250[,3]
  q4 = quantile(ratio4, probs = c(low, high), na.rm = T)
  ratio4 = ratio4[(ratio4 >= q4[1])&(ratio4 <= q4[2])] 
  
  mean_value = mean(rqal0_60_250, na.rm = T)
  sd_value = sd(rqal0_60_250, na.rm = T)
  mean_sd = mean(se_est, na.rm = T)
  mean_eta = apply(coef1_60_250,2,mean, na.rm = T)
  mean_eta_ratio = c(mean(ratio1, na.rm = T), mean(ratio2, na.rm = T), mean(ratio3, na.rm = T))
  sd_eta_ratio = c(sd(ratio1, na.rm = T), sd(ratio2, na.rm = T), sd(ratio3, na.rm = T))
  mean_norm_eta = apply(coef1_norm_60_250,2,mean, na.rm = T)
  sd_norm_eta = apply(coef1_60_250,2,sd, na.rm = T)
  true_val = ifelse(stage == 60, 21.04, 31.74)
  cp = mean((rqal0_60_250 - 1.96*se_est < true_val  & rqal0_60_250 + 1.96*se_est > true_val), na.rm = T)
  mean_po = mean(rmst_po_60_250, na.rm = T)
  sd_po = sd(rmst_po_60_250, na.rm = T)
  mean_mls = mean(1-mr_po_60_250, na.rm = T)
  sd_mls = sd(1-mr_po_60_250, na.rm = T)
  
  return(list(mean_value=mean_value, sd_value=sd_value, mean_sd=mean_sd, mean_eta_ratio=mean_eta_ratio, 
              sd_eta_ratio=sd_eta_ratio, true_val=true_val,cp=cp,mean_po=mean_po,
              sd_po=sd_po,mean_mls=mean_mls,sd_mls=sd_mls))
}


#========================= generating the tables with data =============================#
#================= table 1 ========================#
# ipw 50 250
setwd("F:/QAL/newsim2019/logit/output_K60n250qaloptupp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/hal/output_K60n250qalopt_hal_upp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/rf/output_K60n250qalopt_rf_upp26")
a=summ_tab(stage=60, iter=500)
a
# bc-ipw 50 250
setwd("F:/QAL/newsim2019/logit_s/output_K60n250qalopt_bc_beta_upp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/hal_s/output_K60n250qalopt_bc_hal_beta_upp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/rf_s/output_K60n250qalopt_bc_rf_beta_upp26")
a=summ_tab(stage=60, iter=500)
a

# ipw 50 500
setwd("F:/QAL/newsim2019/logit/output_K60n500qaloptupp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/hal/output_K60n500qalopt_hal_upp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/rf/output_K60n500qalopt_rf_upp26")
a=summ_tab(stage=60, iter=500)
a
# bc-ipw 50 500
setwd("F:/QAL/newsim2019/logit_s/output_K60n500qalopt_bc_beta_upp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/hal_s/output_K60n500qalopt_bc_hal_beta_upp26")
a=summ_tab(stage=60, iter=500)
a
setwd("F:/QAL/newsim2019/rf_s/output_K60n500qalopt_bc_rf_beta_upp26")
a=summ_tab(stage=60, iter=500)
a

# ipw 100 250
setwd("F:/QAL/newsim2019/logit/output_K100n250qaloptupp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/hal/output_K100n250qalopt_hal_upp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/rf/output_K100n250qalopt_rf_upp36")
a=summ_tab(stage=100, iter=500)
a
# bc-ipw 100 250
setwd("F:/QAL/newsim2019/logit_s/output_K100n250qalopt_bc_beta_upp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/hal_s/output_K100n250qalopt_bc_hal_beta_upp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/rf_s/output_K100n250qalopt_bc_rf_beta_upp36")
a=summ_tab(stage=100, iter=500)
a

# ipw 100 500
setwd("F:/QAL/newsim2019/logit/output_K100n500qaloptupp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/hal/output_K100n500qalopt_hal_upp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/rf/output_K100n500qalopt_rf_upp36")
a=summ_tab(stage=100, iter=500)
a
# bc-ipw 100 500
setwd("F:/QAL/newsim2019/logit_s/output_K100n500qalopt_bc_beta_upp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/hal_s/output_K100n500qalopt_bc_hal_beta_upp36")
a=summ_tab(stage=100, iter=500)
a
setwd("F:/QAL/newsim2019/rf_s/output_K100n500qalopt_bc_rf_beta_upp36")
a=summ_tab(stage=100, iter=500)
a
#================== end ==========================#
