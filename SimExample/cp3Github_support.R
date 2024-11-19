### Created by Hao Sun at 02/23/2020

#===========================================================================================#
#========================= main function ===================================================#
#===========================================================================================#

sim_gene <- function(n=250, K=60, interval=10, upp=26, eta_opt=c(1,-1,-1), scenario=1, nuisance='logit', est='IPW', 
                     seed=NULL){
  library(dplyr)
  library(rgenoud)
  
  ### global parameter 
  # sample size
  n = n
  # study duration
  K = K
  # gap between 2 landmarks
  interval = interval
  # number of stages
  stage = K/interval
  # the upper limit of the integral for RQAL
  upp = upp #26? or 36?
  # the true optimal regime
  eta = eta_opt
  time = c(1:upp)
  #scenario=scenario
  #nuisance=nuisance
  
  ### specify output folder
  today <- paste0("K",K, "n", n,"_scenario",scenario, "_nuisance",nuisance, "_est",est, "_upp",upp)
  output_folder <- paste0("output_",today)
  if (!file.exists(output_folder)){dir.create(output_folder)}
  
  
  ### simulation studies
  
  #Use SLURM_ID as seed if not specified
  if(is.null(seed)){
    ## SLURM_ARRAY_TASK_ID on bluehive
    SLURM_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    seed <- SLURM_ID + 5000
  }
  
  #Set seed
  set.seed(seed)
  
  
  ### data generation and transformation
  dat = data_gen_ch3_unif(n,K, interval, eta_opt = eta_opt, scenario = scenario, nuisance = nuisance)
  
  #=======================================================================#
  #================ step 1: the eta_opt estimation =======================#
  #=======================================================================#
  print('Parameter estimation...')
  if(TRUE){
    mod = optimal_sim1(n = n, dat = dat, nuisance = nuisance, est = est, upp = upp)
  }else{
    mod = optimal_sim2(n = n, dat = dat, nuisance = nuisance, est = est, upp = upp)
  } 
  
  
  #=======================================================================#
  #=== step 2: calculate the expected value under estimated regimes ======#
  #=======================================================================#
  if(!is.null(mod)){
    print('Counterfactural estimation...')
    # we approximate the true expected value by stochastic integration; this step would be pretty time consuming
    po1 = data_gen_true_ch3_unif(10000,K, interval, eta = mod$eta_opt, eta_opt = eta_opt)
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
    
    result <- list(mod = mod, rmst_po = rmst1, mr_po = mr1)
    
    save(result, file = paste0(output_folder,"/QAL",seed,".rda"))
  }
}

#===========================================================================================#
#================ support function file for the simulation study ===========================#
#===========================================================================================#
### the main function for estimation and ae estimation after we get the estimated optimal treatment regime
loss_indicator2 <- function(n, eta, dat, trt, cen, clean_dat, a, final_cen, upp = 26, time, target = "rmst", 
                            est_var = F, tole = 0, nuisance = "logit"){
  n_dat = clean_dat
  ### the un-smoothed surv-prob estimation
  pp = pp2 = c()
  
  ### the un-smoothed variance estimation
  var_est = var_est2 = c()
  
  ### construct the follow indicator
  n_dat$follow = n_dat$trt * ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0) +
    (1-n_dat$trt) * (1 - ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0))
  n_dat$follow[clean_dat$stage == 1] = 1
  n_dat = n_dat %>% group_by(ID) %>% mutate(follow_cum = cumprod(follow)) %>% as.data.frame
  
  ### the idx of subjects censored or not
  #final_cen = n_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
  
  if(est_var == T & nuisance == "logit"){
    ### the influence function of theta and beta model
    # step 2: calculate the hessian and influence function of treatment model
    eli_dat = n_dat[n_dat$stage > 1,]
    eli_dat$pa = trt$assign.prob
    eli_dat$score_a = eli_dat$trt - (eli_dat$pa * eli_dat$trt + (1 - eli_dat$pa) * (1 - eli_dat$trt))
    hessian_a =  t(cbind(1, eli_dat$x1, eli_dat$x2)[trt$selvar==1,]) %*% 
      diag((trt$assign.prob * (1 - trt$assign.prob))[trt$selvar==1]) %*% 
      cbind(1, eli_dat$x1, eli_dat$x2)[trt$selvar==1,]
    score_a = eli_dat %>% mutate(const = 1) %>% select(ID, const, x1, x2, score_a) %>% as.data.frame
    score_a[,2:4] = sweep(score_a[2:4], MARGIN=1, score_a$score_a, `*`)
    score_a_sub = lapply(sort(unique(score_a$ID)), function(i){
      apply(score_a[score_a$ID == i, 2:4],2,sum)
    })
    score_a_sub = do.call("rbind", score_a_sub)
    Score_a_sub = matrix(0, nrow = n, ncol = 3)
    # this is the influence function of theta model
    Score_a_sub[sort(unique(score_a$ID)),] = score_a_sub %*% (-solve(hessian_a))
    
    eli_dat$pc = cen$assign.prob
    eli_dat$score_c = eli_dat$cen - (eli_dat$pc * eli_dat$cen + (1 - eli_dat$pc) * (1 - eli_dat$cen))
    hessian_c =  t(cbind(1, eli_dat$x1, eli_dat$x2)[cen$selvar==1,]) %*% 
      diag((cen$assign.prob * (1 - cen$assign.prob))[cen$selvar==1]) %*% 
      cbind(1, eli_dat$x1, eli_dat$x2)[cen$selvar==1,]
    score_c = eli_dat %>% mutate(const = 1) %>% select(ID, const, x1, x2, score_c) %>% as.data.frame
    score_c[,2:4] = sweep(score_c[2:4], MARGIN=1, score_c$score_c, `*`)
    score_c_sub = lapply(sort(unique(score_c$ID)), function(i){
      apply(score_c[score_c$ID == i, 2:4],2,sum)
    })
    score_c_sub = do.call("rbind", score_c_sub)
    Score_c_sub = matrix(0, nrow = n, ncol = 3)
    # this is the influence function of theta model
    Score_c_sub[sort(unique(score_c$ID)),] = score_c_sub %*% (-solve(hessian_c)) 
  }
  
  ### the target is surv-prob: not implement
  if(target == "surv"){
    for(i in 1:length(time)){
      idx =  n_dat %>% group_by(ID) %>% summarise(idx = ifelse(sum(qal > a[i]) > 0, which(qal > a[i])[1], 0))
      idx = idx$idx
      
      ### the indicator of still being observable: it is I(C > t) * I(U > t)
      delta_o = idx > 0
      
      ### the indicator for denominator part
      id_d = which((1 - delta_o)*(1 - final_cen$final_cen) == 1)
      d_idx = unlist(lapply(id_d, function(x) dim(n_dat[n_dat$ID == x,])[1]))
      idx_d = idx
      idx_d[id_d] = d_idx
      delta_d = idx_d > 0
      
      final_weight = n_dat %>% group_by(ID) %>% summarise(f_weight = ifelse(delta_o[ID[1]] > 0, follow_cum[idx[ID[1]]]*weight[idx[ID[1]]],0))
      
      final_weight2 = n_dat %>% group_by(ID) %>% summarise(f_weight = ifelse(delta_d[ID[1]] > 0, follow_cum[idx_d[ID[1]]]*weight[idx_d[ID[1]]],0))
      
      #pp[i] = sum(final_weight$f_weight)/n
      pp2[i] = sum(final_weight$f_weight)/sum(final_weight2$f_weight)
      
      
      if(est_var == T & nuisance == "logit"){
        ### calculate the variance estimation: 12/02/2018
        # step 1: calculate the gradient function
        effect_id = which(final_weight$f_weight > 0)
        grad_a_theta = lapply(effect_id, function(i) {d = n_dat[n_dat$ID == i, ][1:idx[i], ];
        apply(sweep(cbind(1, d$x1, d$x2), MARGIN=1, d$grad_a, `*`),2,sum)})
        grad_a_theta = do.call("rbind", grad_a_theta)
        #grad_theta1 = apply(sweep(grad_a_theta, MARGIN=1, -final_weight$f_weight[effect_id], `*`),2,sum)/n
        grad_theta2 = apply(sweep(grad_a_theta, MARGIN=1, 
                                  final_weight$f_weight[effect_id]*(final_weight2$f_weight[effect_id] - sum(final_weight2$f_weight)), `*`),2,sum)/(sum(final_weight2$f_weight)^2)
        
        grad_c_beta = lapply(effect_id, function(i) {d = n_dat[n_dat$ID == i, ][1:idx[i], ];
        apply(sweep(cbind(1, d$x1, d$x2), MARGIN=1, d$grad_c, `*`),2,sum)})
        grad_c_beta = do.call("rbind", grad_c_beta)
        #grad_beta1 = apply(sweep(grad_c_beta, MARGIN=1, -final_weight$f_weight[effect_id], `*`),2,sum)/n
        grad_beta2 = apply(sweep(grad_c_beta, MARGIN=1, 
                                 final_weight$f_weight[effect_id]*(final_weight2$f_weight[effect_id] - sum(final_weight2$f_weight)), `*`),2,sum)/(sum(final_weight2$f_weight)^2)
        
        # step 3: empirical variance estimation
        #IC_eta = final_weight$f_weight - pp[i]
        #IC_theta = Score_a_sub %*% grad_theta1
        #IC_beta = Score_c_sub %*% grad_beta1
        
        IC_eta2 = (final_weight$f_weight - final_weight2$f_weight*pp2[i])/sum(final_weight2$f_weight)*n
        IC_theta2 = Score_a_sub %*% grad_theta2
        IC_beta2 = Score_c_sub %*% grad_beta2
        
        #IC = IC_eta + IC_theta + IC_beta
        IC2 = IC_eta2 + IC_theta2 + IC_beta2
        #var_est[i] = sum(IC^2)/n^2
        var_est2[i] = sum(IC2^2)/n^2
      }
      else if(est_var == T & nuisance == "hal"){
        IC_eta2 = (final_weight$f_weight - final_weight2$f_weight*pp2[i])/sum(final_weight2$f_weight)*n
        
        IC2 = IC_eta2 
        var_est2[i] = sum(IC2^2)/n^2
      }
    }
    return(list(pp2 = pp2, var_est2 = var_est2))
  }
  ### the target is rmst
  else if(target == "rmst"){
    IC_mat = IC_mat2 = matrix(NA, nrow = n, ncol = length(a)-1)
    ### here, we focus on rmst
    for(i in 1:(length(a) - 1)){
      ### find the row for each one when qal is a[i]
      idx =  n_dat %>% group_by(ID) %>% summarise(idx = ifelse(sum(qal > a[i]) > 0, which(qal > a[i])[1], 0))
      idx = idx$idx
      
      ### the indicator of still being observable: it is I(C > t) * I(U > t)
      delta_o = idx > 0
      
      ### the indicator for denominator part
      id_d = which((1 - delta_o)*(1 - final_cen$final_cen) == 1)
      d_idx = unlist(lapply(id_d, function(x) dim(n_dat[n_dat$ID == x,])[1]))
      idx_d = idx
      idx_d[id_d] = d_idx
      delta_d = idx_d > 0
      
      final_weight = n_dat %>% group_by(ID) %>% summarise(f_weight = ifelse(delta_o[ID[1]] > 0, follow_cum[idx[ID[1]]]*weight[idx[ID[1]]],0))
      
      final_weight2 = n_dat %>% group_by(ID) %>% summarise(f_weight = ifelse(delta_d[ID[1]] > 0, follow_cum[idx_d[ID[1]]]*weight[idx_d[ID[1]]],0))
      
      #pp[i] = sum(final_weight$f_weight)/n
      pp2[i] = sum(final_weight$f_weight)/sum(final_weight2$f_weight)
      
      if(est_var == T & nuisance == "logit"){
        ### calculate the variance estimation: 12/02/2018
        # step 1: calculate the gradient function
        effect_id = which(final_weight$f_weight > 0)
        
        if(length(effect_id) > 0){
          grad_a_theta = lapply(effect_id, function(i) {d = n_dat[n_dat$ID == i, ][1:idx[i], ];
          apply(sweep(cbind(1, d$x1, d$x2), MARGIN=1, d$grad_a, `*`),2,sum)})
          grad_a_theta = do.call("rbind", grad_a_theta)
          #grad_theta1 = apply(sweep(grad_a_theta, MARGIN=1, -final_weight$f_weight[effect_id], `*`),2,sum)/n
          grad_theta2 = apply(sweep(grad_a_theta, MARGIN=1, 
                                    final_weight$f_weight[effect_id]*(final_weight2$f_weight[effect_id] - sum(final_weight2$f_weight)), `*`),2,sum)/(sum(final_weight2$f_weight)^2)
          
          grad_c_beta = lapply(effect_id, function(i) {d = n_dat[n_dat$ID == i, ][1:idx[i], ];
          apply(sweep(cbind(1, d$x1, d$x2), MARGIN=1, d$grad_c, `*`),2,sum)})
          grad_c_beta = do.call("rbind", grad_c_beta)
          #grad_beta1 = apply(sweep(grad_c_beta, MARGIN=1, -final_weight$f_weight[effect_id], `*`),2,sum)/n
          grad_beta2 = apply(sweep(grad_c_beta, MARGIN=1, 
                                   final_weight$f_weight[effect_id]*(final_weight2$f_weight[effect_id] - sum(final_weight2$f_weight)), `*`),2,sum)/(sum(final_weight2$f_weight)^2)
          
          # step 3: empirical variance estimation
          #IC_eta = final_weight$f_weight - pp[i]
          #IC_theta = Score_a_sub %*% grad_theta1
          #IC_beta = Score_c_sub %*% grad_beta1
          
          IC_eta2 = (final_weight$f_weight - final_weight2$f_weight*pp2[i])/sum(final_weight2$f_weight)*n
          IC_theta2 = Score_a_sub %*% grad_theta2
          IC_beta2 = Score_c_sub %*% grad_beta2
          
          #IC = IC_eta + IC_theta + IC_beta
          IC2 = IC_eta2 + IC_theta2 + IC_beta2
          
          #IC_mat[,i] = IC
          IC_mat2[,i] = IC2 
        }
        else{
          IC_mat2[,i] = 0
        }
        # calculate the variance
        td = diff(a)
        #IC = apply(sweep(IC_mat, MARGIN=2, td, `*`),1,sum)
        IC2 = apply(sweep(IC_mat2, MARGIN=2, td, `*`),1,sum)
        var_2 = ifelse((upp == 36) & (sum(IC2^2)/n^2 < tole^2), NA, sum(IC2^2)/n^2)
      }
      else if(est_var == T & nuisance == "hal"){
        ### calculate the variance estimation: 12/02/2018
        # step 1: calculate the gradient function
        effect_id = which(final_weight$f_weight > 0)
        
        if(length(effect_id) > 0){
          IC_eta2 = (final_weight$f_weight - final_weight2$f_weight*pp2[i])/sum(final_weight2$f_weight)*n
          IC2 = IC_eta2
          IC_mat2[,i] = IC2 
        }
        else{
          IC_mat2[,i] = 0
        }
        # calculate the variance
        td = diff(a)
        #IC = apply(sweep(IC_mat, MARGIN=2, td, `*`),1,sum)
        IC2 = apply(sweep(IC_mat2, MARGIN=2, td, `*`),1,sum)
        var_2 = ifelse((upp == 36) & (sum(IC2^2)/n^2 < tole^2), NA, sum(IC2^2)/n^2)
      }
    }
    if(!is.na(var_2)){
      return(list(rqal2 = sum(td * pp2), var_est2 = sum(IC2^2)/n^2)) 
    }
    else{
      NA
    }
  }
}

### some wrapper function used for optimization prcedure
### loss function for IPW
indicator_opt1 <- function(eta, clean_dat, a, uncen_idx){
  n_dat = clean_dat
  
  pp = c()
  
  ### construct the follow indicator
  n_dat$follow = n_dat$trt * ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0) +
    (1-n_dat$trt) * (1 - ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0))
  n_dat$follow[clean_dat$stage == 1] = 1
  n_dat = n_dat %>% group_by(ID) %>% mutate(follow_cum = cumprod(follow)) %>% as.data.frame
  
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
  
  pp = unlist(lapply(1:(length(a) - 1), function(x) prob_surv(x)))
  #   }
  td = diff(a)
  return(sum(td * pp))
  
  
  # if(target == "rmst"){
  #   pp = unlist(lapply(1:(length(a) - 1), function(x) prob_surv(a[x])))
  #   #   }
  #   td = diff(a)
  #   return(sum(td * pp)) 
  # }
  # else{
  #   prob = c()
  #   for(i in 1:length(time)){
  #     prob[i] = prob_surv(time[i])
  #   }
  #   return(prob)
  # }
}
indicator_opt2 <- function(eta, clean_dat, a, uncen_idx){
  n_dat = clean_dat
  
  pp = c()
  
  ### construct the follow indicator
  n_dat$follow = n_dat$trt * ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0) +
    (1-n_dat$trt) * (1 - ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0))
  n_dat$follow[clean_dat$stage == 1] = 1
  n_dat = n_dat %>% group_by(ID) %>% mutate(follow_cum = cumprod(follow)) %>% as.data.frame
  
  prob_surv <- function(i){
    sub1 = n_dat[n_dat$qal >a[i],]
    sub1 = sub1[!duplicated(sub1$ID),]
    numerator = sum(sub1$follow_cum * sub1$weight)
    
    return(numerator/n)
  }
  
  pp = unlist(lapply(1:(length(a) - 1), function(x) prob_surv(x)))
  #   }
  td = diff(a)
  return(sum(td * pp))
}

### loss function for BC-IPW 
smooth_opt1 <- function(eta, clean_dat, a, uncen_idx,div, nSamp){
  n_dat = clean_dat
  
  pp = c()
  
  ### construct the trt_last record
  n_dat$trt_last = 0
  n_dat$trt_last[which(n_dat$stage > 1)] = n_dat$trt[which(n_dat$stage > 1)-1]
  
  score = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2)
  s = sd(score)
  s = ifelse(s > 0, s, eta[1]>=0)
  score_mono = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2)^(1 - n_dat$trt_last)
  bw = nSamp^(-1/3) * s/n_dat$stage*div
  
  n_dat$normal = pnorm(score_mono/bw)
  
  # the smooth indicator function w.r.t eta
  n_dat$smooth = n_dat$trt*n_dat$normal + (1-n_dat$trt) * (1 - n_dat$normal)
  ### because stage one is always non-treatment by design
  n_dat$smooth[clean_dat$stage == 1] = 1
  n_dat = n_dat %>% group_by(ID) %>% mutate(smooth_cum = cumprod(smooth)) %>% as.data.frame
  
  prob_surv <- function(i){
    sub1 = n_dat[n_dat$qal >a[i],]
    sub1 = sub1[!duplicated(sub1$ID),]
    numerator = sum(sub1$smooth_cum * sub1$weight)
    
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
      denominator = sum(sub2$smooth_cum * sub2$weight)
      
      return(numerator/denominator)
    }
    
  }
  
  pp = unlist(lapply(1:(length(a) - 1), function(x) prob_surv(x)))
  #   }
  td = diff(a)
  return(sum(td * pp))
}
smooth_opt2 <- function(eta, clean_dat, a, uncen_idx){
  n_dat = clean_dat
  
  pp = c()
  
  ### construct the trt_last record
  n_dat$trt_last = 0
  n_dat$trt_last[which(n_dat$stage > 1)] = n_dat$trt[which(n_dat$stage > 1)-1]
  
  score = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2)
  s = sd(score)
  s = ifelse(s > 0, s, eta[1]>=0)
  score_mono = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2)^(1 - n_dat$trt_last)
  bw = n^(-1/3) * s/stage*2
  
  n_dat$normal = pnorm(score_mono/bw)
  
  # the smooth indicator function w.r.t eta
  n_dat$smooth = n_dat$trt*n_dat$normal + (1-n_dat$trt) * (1 - n_dat$normal)
  ### because stage one is always non-treatment by design
  n_dat$smooth[clean_dat$stage == 1] = 1
  n_dat = n_dat %>% group_by(ID) %>% mutate(smooth_cum = cumprod(smooth)) %>% as.data.frame
  
  prob_surv <- function(i){
    sub1 = n_dat[n_dat$qal >a[i],]
    sub1 = sub1[!duplicated(sub1$ID),]
    numerator = sum(sub1$smooth_cum * sub1$weight)
    
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
      denominator = sum(sub2$smooth_cum * sub2$weight)
      
      return(numerator/n)
    }
    
  }
  
  pp = unlist(lapply(1:(length(a) - 1), function(x) prob_surv(x)))
  #   }
  td = diff(a)
  return(sum(td * pp))
}

### the function to generate similation data
data_gen_ch3_unif <- function(n,K, interval, eta_opt, scenario, nuisance){
  stage = K/interval
  
  y0 = -6
  y1 = -0.5
  y2 = -0.5
  y3 = -0.5
  
  a0 = 0.5 - stage * 0.1
  a1 = -0.5
  a2 = -0.5
  a3 = 1
  
  c0 = -2 + 0.5*a0
  c1 = -1
  c2 = -1
  c3 = 1
  
  dat = matrix(0, ncol = 12, nrow = K/interval * n)
  colnames(dat) = c("ID", "stage", "x1", "x2", "z1", "z2", "qol", "trt", "outcome","censor", "follow", "trt_opt")
  dat = as.data.frame(dat)
  
  dat[,1] <- rep(1:n, each = stage) #the ID variable
  dat[,2] <- rep(1:stage, time = n) #the time index
  
  for(i in 1:n){
    ### no treatment at the first stage
    x1 = runif(1, 0.6, 1)
    x2 = runif(1, 0.6, 1)
    #x3 = exp(4*x1)
    A = 0
    C = 0
    qol = runif(1, 0.6, 1) #adjusted quality of life
    
    d = dat[dat$ID == i, ]
    
    ### the indicator to record if the subject follow the optimal regime; at time 1, always follow
    Follow = 1
    Penalty = 0
    A_opt = 0
    
    ### specify baseline information
    next_stage = T
    t = 1 #stage indicator
    
    ### count the inner weeks of survival
    y = 0
    for(j in 1:interval){
      y = y + 1
      Y = (exp(y0 + t * (-y0*0.75)/stage + y1 * x1 + y2 * x2 + y3 * A * Follow + Penalty)/
             (1 + exp(y0 + y1 * x1 + y2 * x2 + y3 * A * Follow + Penalty))) > runif(1) #Indicator 1 if censored or treatment during a week
      
      if(Y == 1){
        next_stage = F
        break
      }
    }
    
    d[t,]$x1 = x1
    d[t,]$x2 = x2
    d[t,]$z1 = (x1+x2-1)^2
    d[t,]$z2 = (x1-0.5)^2
    #d[t,]$z1 = (x1+x2-1)^2
    #d[t,]$z2 = (x1-0.5)^2
    #d[t,]$z2 = exp(x1/2)
    #d[t,]$z1 = 10*(x1-0.5)^3
    #d[t,]$z1 = (x1-x2+0.5)^2
    #d[t,]$z1 = exp(x1/2)
    d[t,]$qol = qol
    d[t,]$outcome = y
    d[t,]$trt = A
    d[t,]$censor = C
    d[t,]$follow = Follow
    d[t,]$trt_opt = A_opt
    
    ### for the next interval
    while(next_stage == T & t < stage & C == 0){
      t = t + 1
      
      # update data
      # if Y = 0 (live), collect new data and do treatment assignment
      x1 = runif(1, 0.4 + 0.02 * stage, 1) * x1 # here, I only consider a simpliet case; scale the covariates
      x2 = runif(1, 0.4 + 0.02 * stage, 1) * x2
      #x3 = exp(4*x1)
      qol = runif(1, (0.4 + 0.02 * stage - Penalty*0.1), (1 - Penalty*0.1)) * qol
      # judge treatment assignment
      if(A != 1){
        A = (exp(a0 + a1 * x1 + a2 * x2 )/(1 + exp(a0 + a1 * x1 + a2 * x2 ))) > runif(1)
      }
      if(A_opt != 1){
        A_opt = I(eta_opt[1] + eta_opt[2] * x1 + eta_opt[3] * x2 > (scenario-1)*ifelse(nuisance=='logit',-.05, 0))
      }
      
      Follow = ifelse(Follow == 1, A == A_opt, 0)
      
      ### here, I just consider the simplist penalty, constant penalty for any one who don't follow the opt
      ### more complex penalty can depend on the distance to opt and direction(may consider it later)
      Penalty = ifelse(Follow == 1, 0, 1)
      
      # judge ceonsoring
      C =(exp(c0 + c1 * x1 + c2 * x2 )/(1 + exp(c0 + c1 * x1 + c2 * x2 ))) > runif(1) 
      
      ### coutnt the inner weeks of survival
      y = 0
      for(j in 1:interval){
        y = y + 1
        Y = (exp(y0 + t * (-y0*0.75)/stage + y1 * x1 + y2 * x2 + y3 * A * Follow + Penalty)/
               (1 + exp(y0 + y1 * x1 + y2 * x2 + y3 * A * Follow + Penalty))) > runif(1)
        if(Y == 1){
          next_stage = F
          break
        }
      }
      
      d[t,]$x1 = x1
      d[t,]$x2 = x2
      d[t,]$z1 = (x1+x2-1)^2
      d[t,]$z2 = (x1-0.5)^2
      #d[t,]$z1 = (x1+x2-1)^2
      #d[t,]$z2 = (x1-1)^2
      #d[t,]$z2 = exp(x1/2)
      #d[t,]$z1 = 10*(x1-0.5)^3
      #d[t,]$z1 = (x1-x2+0.5)^2
      #d[t,]$z1 = exp(x1/2)
      d[t,]$qol = qol
      d[t,]$outcome = y*(1-C)
      d[t,]$trt = A
      d[t,]$censor = C
      d[t,]$follow = Follow
      d[t,]$trt_opt = A_opt
    }
    
    dat[dat$ID == i,] = d
  }
  
  return(dat)
}

### the function to generate data under the true optimal regime; used to calculated the true valuse of some causal estimand
data_gen_true_ch3_unif <- function(n,K, interval, eta, eta_opt, seed = 12345){
  y0 = -6
  y1 = -0.5
  y2 = -0.5
  y3 = -0.5
  
  dat = matrix(0, ncol = 10, nrow = K/interval * n)
  colnames(dat) = c("ID", "stage", "x1", "x2", "qol", "trt", "outcome","censor", "follow", "trt_opt")
  dat = as.data.frame(dat)
  stage = K/interval
  
  dat[,1] <- rep(1:n, each = stage) #the ID variable
  dat[,2] <- rep(1:stage, time = n) #the time index
  
  for(i in 1:n){
    set.seed(seed+i)
    ### no treatment at the first stage
    x1 = runif(1, 0.6, 1)
    x2 = runif(1, 0.6, 1)
    A = 0
    C = 0
    
    ### the indicator to record if the subject follow the optimal regime; at time 1, always follow
    Follow = 1
    Penalty = 0
    A_opt = 0
    
    qol = runif(1, 0.6, 1)
    
    d = dat[dat$ID == i, ]
    
    ### specify baseline information
    next_stage = T
    t = 1
    
    ### count the inner weeks of survival
    y = 0
    for(j in 1:interval){
      y = y + 1
      Y = (exp(y0 + t * (-y0*0.75)/stage + y1 * x1 + y2 * x2 + y3 * A * Follow)/
             (1 + exp(y0 + y1 * x1 + y2 * x2 + y3 * A * Follow))) > runif(1)
      
      if(Y == 1){
        next_stage = F
        break
      }
    }
    
    d[t,]$x1 = x1
    d[t,]$x2 = x2
    d[t,]$qol = qol
    d[t,]$outcome = y
    d[t,]$trt = A
    d[t,]$censor = C
    d[t,]$follow = Follow
    d[t,]$trt_opt = A_opt
    ### for the next interval
    while(next_stage == T & t < stage & C == 0){
      t = t + 1
      
      # update data
      # if Y = 0 (live), collect new data and do treatment assignment
      x1 = runif(1, 0.4 + 0.02 * stage, 1) * x1 # here, I only consider a simpliet case; scale the covariates
      x2 = runif(1, 0.4 + 0.02 * stage, 1) * x2
      qol = runif(1, 0.4 + 0.02 * stage - Penalty*0.1, 1 - Penalty * 0.1) * qol
      # judge treatment assignment
      if(A != 1){
        A = I(eta[1] + eta[2] * x1 + eta[3] * x2 > 0)
      }
      if(A_opt != 1){
        A_opt = I(eta_opt[1] + eta_opt[2] * x1 + eta_opt[3] * x2 > 0)
      }
      
      Follow = ifelse(Follow == 1, A == A_opt, 0)
      
      ### here, I just consider the simplist penalty, constant penalty for any one who don't follow the opt
      ### more complex penalty can depend on the distance to opt and direction(may consider it later)
      Penalty = ifelse(Follow == 1, 0, 1)
      
      # judge ceonsoring
      #C =(exp(c0 + c1 * x1 + c2 * x2)/(1 + exp(c0 + c1 * x1 + c2 * x2))) > runif(1)
      
      ### coutnt the inner weeks of survival
      y = 0
      for(j in 1:interval){
        y = y + 1
        Y = (exp(y0 + t * (-y0*0.75)/stage + y1 * x1 + y2 * x2 + y3 * A * Follow + Penalty)/
               (1 + exp(y0 + y1 * x1 + y2 * x2 + y3 * A * Follow + Penalty))) > runif(1)
        if(Y == 1){
          next_stage = F
          break
        }
      }
      
      d[t,]$x1 = x1
      d[t,]$x2 = x2
      d[t,]$qol = qol
      d[t,]$outcome = y*(1-C)
      d[t,]$trt = A
      d[t,]$censor = C
      d[t,]$follow = Follow
      d[t,]$trt_opt = A_opt
    }
    
    dat[dat$ID == i,] = d
  }
  
  return(dat)
}

#===========================================================#
#====== some functions for simulation study ================#
#===========================================================#

### simulation function when logit model is correctly specified
### this function works for all simulation scenarios for logit and rf, and only n = 250 cases for HAL
optimal_sim1 <- function(n = 250, dat, nuisance = "logit", est = "IPW", upp){
  
  ### clean the data
  dat$start = dat$stage - 1
  clean_dat = dat[dat$x1 != 0,]
  KK = max(dat$stage)
  div = ifelse(KK > 10, 2, 1)
  
  if(nuisance == "logit"){
    ### fit logit model for treatment and censoring
    trt <- ipwtm_qal(exposure = trt, family = "binomial", link = "logit",
                     numerator = NULL, denominator = ~ x1+x2,
                     id = ID, tstart = start, timevar = stage, type = "first",
                     data = clean_dat[clean_dat$stage > 1, ])
    cen <- ipwtm_qal(exposure = censor, family = "binomial", link = "logit",
                     numerator = NULL, denominator = ~ x1+x2,
                     id = ID, tstart = start, timevar = stage, type = "first",
                     data = clean_dat[clean_dat$stage > 1, ]) 
    
    ### assign the weight function
    clean_dat$weight = 1
    clean_dat$weight[clean_dat$stage > 1] = cen$ipw.weights * trt$ipw.weights
    clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
    
    ### the influence function part for logistics regression, which will be used in the se estimation
    ### assign the gradient function of dP_a,i/d\theta and dP_c,i/d\beta:
    clean_dat$grad_a = 0
    clean_nohead = clean_dat[clean_dat$stage > 1,]
    clean_nohead$grad_a[trt$selvar==1 & clean_nohead$trt == 1] = 1 - trt$assign.prob[trt$selvar==1 & clean_nohead$trt == 1]
    clean_nohead$grad_a[trt$selvar==1 & clean_nohead$trt == 0] = trt$assign.prob[trt$selvar==1 & clean_nohead$trt == 0] - 1
    clean_dat$grad_a[clean_dat$stage > 1] = clean_nohead$grad_a
    
    clean_dat$grad_c = 0
    clean_nohead = clean_dat[clean_dat$stage > 1,]
    clean_nohead$grad_c[cen$selvar==1 & clean_nohead$censor == 0] = cen$assign.prob[cen$selvar==1 & clean_nohead$censor == 0] - 1
    clean_dat$grad_c[clean_dat$stage > 1] = clean_nohead$grad_c
    
    ### the test time points when computing RQAL
    a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
    a = c(0, sort(unique(pmin(a$qal, upp)))) 
    ### the idx of subjects censored or not
    final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
    uncen_idx = which(final_cen$final_cen == 0)
    
    if(est == "IPW"){
      temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx)
    } else{
      temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx, div = div, nSamp = n)
    }
    tole = ifelse((upp==36)&(est!='IPW'), 0.4, 0)
    mod1 = loss_indicator2(n, temp1$par, dat, trt=trt, cen=cen, a=a, final_cen=final_cen, clean_dat, upp = upp, time, 
                           target = "rmst", est_var = T, nuisance = "logit",tole=tole)
    if(all(all(!is.na(mod1)))){results = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))}
    
  }
  else if(nuisance == "hal"){
    rate = 0.9
    #==============================================================#
    #================== step 1: cv-hal results ====================#
    #==============================================================#
    if(TRUE){
      library(hal9001)
      tempdat = clean_dat[clean_dat$stage > 1, ]
      tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      mod_trt = fit_hal(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,] , Y = tempdat$trt[tempdat$selvar_trt==1], family = "binomial", yolo = FALSE, 
                        return_x_basis = T)
      
      preds_trt = predict(mod_trt, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,])
      
      tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
      tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
      tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
      tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      mod_cen = fit_hal(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,] , Y = tempdat$censor[tempdat$selvar_cen==1], family = "binomial", yolo = FALSE, 
                        return_x_basis = T)
      
      preds_cen = predict(mod_cen, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,])
      
      tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
      tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
      tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
      tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
      
      ### assign the weight function
      clean_dat$weight = 1
      clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
      clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
      
      ### the test time points when computing RMST
      a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
      a = c(0, sort(unique(pmin(a$qal, upp)))) 
      ### the idx of subjects censored or not
      final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
      uncen_idx = which(final_cen$final_cen == 0)
      
      
       
      if(est == "IPW"){
        temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx)
      } else{
        temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx, div = div, nSamp = n)
      }
      
      mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt, cen=mod_cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                             est_var = T, nuisance = "hal")
      
      cv_se = sqrt(mod1$var_est2)
      cv_lambda_trt = mod_trt$lambda_star
      coefs_trt = mod_trt$coefs 
      cv_lambda_cen = mod_cen$lambda_star
      coefs_cen = mod_cen$coefs 
      
      #==============================================================#
      #============= step 2: undersmooth-hal results ================#
      #==============================================================#
      
      
      #======================================#
      #== step 2.1: for trt model ===========#
      #======================================#
      
      #preds_trt = predict(mod_trt, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,])
      #basis_func_trt = as.matrix(cbind(1, mod_trt$x_basis)) ### n by p, the column 1 is for intercept term
      A = tempdat$trt[tempdat$selvar_trt==1]
      A_hat = preds_trt
      
      #NEW UNDERSMOOTHING (DEFINED AT BOTTOM OF FILE)
      underTrt <- newUndersmooth(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,],
                                 Y = tempdat$trt[tempdat$selvar_trt==1],
                                 A = A,
                                 cv_lambda = cv_lambda_trt)
      
      #OLD UNDERSMOOTHING
      #cutoff_trt = cv_se/log(n)/sqrt(n)##################################################################OLD UNDERSMOOTHING
      #cost_func_trt = sum(abs(coefs_trt) * abs(t(basis_func_trt) %*% as.matrix(A - A_hat)/n))############HERE
      #if(cost_func_trt <= cutoff_trt) {
      #  under_lambda_trt = cv_lambda_trt
      #  flag_trt = 1
      #} else{
      #  count = 1
      #  lambda_trt = mod_trt$lambda_star
      #  cost_func_trt_old = cost_func_trt + 1
      #  mod_trt_grid = "good"
      #  
      #  while (cost_func_trt > cutoff_trt & cost_func_trt < cost_func_trt_old) {
      #    lambda_trt = lambda_trt*rate
      #    
      #    tempdat = clean_dat[clean_dat$stage > 1, ]
      #    tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      #    
      #    mod_trt_grid = tryCatch(fit_hal(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,] , Y = tempdat$trt[tempdat$selvar_trt==1], family = "binomial", yolo = FALSE, return_x_basis = T,
      #                                    cv_select = F, lambda = lambda_trt),error=function(e) e, warning=function(w) w)
      #    
      #    if (is(mod_trt_grid,"warning")){
      #      break
      #    }
      #    
      #    preds_trt_grid = predict(mod_trt_grid, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,])
      #    basis_func_trt = as.matrix(cbind(1, mod_trt_grid$x_basis)) ### n by p, the row 1 is for intercept term
      #    A = tempdat$trt[tempdat$selvar_trt==1]
      #    A_hat = preds_trt_grid
      #    
      #    cost_func_trt_old = cost_func_trt
      #    cost_func_trt = sum(abs(coefs_trt) * abs(t(basis_func_trt) %*% as.matrix(A - A_hat)/n))
      #    count = count + 1
      #  }
      #  if(is(mod_trt_grid,"warning") | cost_func_trt >= cost_func_trt_old) { # in this case, cannot attain the boundary 
      #    flag_trt = 2
      #    under_lambda_trt = cv_lambda_trt
      #  } else{ # in this case, attain the boundary
      #    flag_trt = 0
      #    under_lambda_trt = lambda_trt
      #  }
      #}
      
      #======================================#
      #== step 2.2: for cen model ===========#
      #======================================#
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      preds_cen = predict(mod_cen, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,])
      basis_func_cen = as.matrix(cbind(1, mod_cen$x_basis)) ### n by p, the row 1 is for intercept term
      C = tempdat$censor[tempdat$selvar_cen==1]
      C_hat = preds_cen
      
      #NEW UNDERSMOOTHING (DEFINED AT BOTTOM OF FILE)
      underCen <- newUndersmooth(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,],
                                 Y = tempdat$censor[tempdat$selvar_cen==1],
                                 A = C,
                                 cv_lambda = cv_lambda_cen)
      
      #OLD UNDERSMOOTHING
      #cutoff_cen = cv_se/log(n)/sqrt(n)##################################################################OLD UNDERSMOOTHING
      #cost_func_cen = sum(abs(coefs_cen) * abs(t(basis_func_cen) %*% as.matrix(C - C_hat)/n))############HERE
      
      #if(cost_func_cen <= cutoff_cen) {
      #  under_lambda_cen = cv_lambda_cen
      #  flag_cen = 1
      #} else{
      #  count = 1
      #  lambda_cen = mod_cen$lambda_star
      #  cost_func_cen_old = cost_func_cen + 1
      #  mod_cen_grid = "good"
      #  
      #  while (cost_func_cen > cutoff_cen & cost_func_cen < cost_func_cen_old) {
      #    lambda_cen = lambda_cen*rate
      #    
      #    tempdat = clean_dat[clean_dat$stage > 1, ]
      #    tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      #    
      #    mod_cen_grid = tryCatch(fit_hal(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,] , Y = tempdat$censor[tempdat$selvar_cen==1], family = "binomial", yolo = FALSE, return_x_basis = T,
      #                                    cv_select = F, lambda = lambda_cen),error=function(e) e, warning=function(w) w)
      #    
      #    if (is(mod_cen_grid,"warning")){
      #      break
      #    }
      #    
      #    preds_cen_grid = predict(mod_cen_grid, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,])
      #    basis_func_cen = as.matrix(cbind(1, mod_cen_grid$x_basis)) ### n by p, the row 1 is for intercept term
      #    C = tempdat$censor[tempdat$selvar_cen==1]
      #    C_hat = preds_cen_grid
      #    
      #    cost_func_cen_old = cost_func_cen
      #    cost_func_cen = sum(abs(coefs_cen) * abs(t(basis_func_cen) %*% as.matrix(C - C_hat)/n))
      #    count = count + 1
      #  }
      #  if(is(mod_cen_grid,"warning") | cost_func_cen >= cost_func_cen_old) { # in this case, cannot attain the boundary 
      #    flag_cen = 2
      #    under_lambda_cen = cv_lambda_cen
      #  } else{
      #    flag_cen = 0
      #    under_lambda_cen = lambda_cen
      #  }
      #}
      
      
      #======================================#
      #== step 2.3: re-estimate the results =#
      #======================================#
      tole=ifelse(upp==36,0.4,0)
      if(FALSE){#flag_cen != 0 & flag_trt!=0){ ##BLOCKED THIS OUT AS I DON"T THINK WE NEED IT ANYMORE
        mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt, cen=mod_cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                               est_var = T, nuisance = "hal",tole = tole)
        if(all(!is.na(mod1))){results = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))}
      } else{
        tempdat = clean_dat[clean_dat$stage > 1, ]
        tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        #mod_trt = fit_hal(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,] , Y = tempdat$trt[tempdat$selvar_trt==1], family = "binomial", yolo = FALSE, 
        #                  cv_select = F, lambda = under_lambda_trt)
        #preds_trt = predict(mod_trt, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_trt==1,])
        mod_trt <- underTrt$halFit
        preds_trt <- underTrt$preds
        under_lambda_trt <- underTrt$underLam
        
        tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
        tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
        tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
        tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
        
        tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        #mod_cen = fit_hal(X = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,] , Y = tempdat$censor[tempdat$selvar_cen==1], family = "binomial", yolo = FALSE, 
        #                  cv_select = F, lambda = under_lambda_cen)
        #preds_cen = predict(mod_cen, new_data = cbind(1, tempdat$x1, tempdat$x2)[tempdat$selvar_cen==1,])
        mod_cen <- underCen$halFit
        preds_cen <- underCen$preds
        under_lambda_cen <- underCen$underLam
        
        tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
        tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
        tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
        tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
        
        ### assign the weight function
        clean_dat$weight = 1
        clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
        clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
        
        ### the test time points when computing RMST
        a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
        a = c(0, sort(unique(pmin(a$qal, upp)))) 
        ### the idx of subjects censored or not
        final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
        uncen_idx = which(final_cen$final_cen == 0)
        
        if(est == "IPW"){
          temp2 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx)
        } else{
          temp2 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx, div = div, nSamp = n)
        }
        mod2 = loss_indicator2(n, temp2$par, dat, trt=mod_trt, cen=mod_cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", est_var = T, 
                               nuisance = "hal",tole = tole)
        if(all(!is.na(mod2))){
          results = list(Status = 'SUCCESS', eta_opt = temp2$par, est = temp2$val, se = sqrt(mod2$var_est2),
                         trtCrit = underTrt$crits, cenCrit = underCen$crits, trtR = underTrt$r, cenR = underCen$r)
        }else{
            results = list(Status = 'ERROR')
          }
      }
    }
    else if(FALSE){ ### beyond the computation limit of hal, and we need to split the hal model into two parts
      
      #=========================================#
      #==== step 1: fit the cv hal model =======#
      #=========================================#
      
      library(hal9001)
      tempdat = clean_dat[clean_dat$stage > 1, ]
      tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
        return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
      tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
      tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250, ]
      mod_trt1 = fit_hal(X = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2) , 
                         Y = tempdat_sub1$trt, 
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      mod_trt2 = fit_hal(X = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2) , 
                         Y = tempdat_sub2$trt, 
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      preds_trt = c()
      preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
      preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
      
      tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
      tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
      tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
      tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
      tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
      tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
      mod_cen1 = fit_hal(X = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2) , 
                         Y = tempdat_sub1$censor, 
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      mod_cen2 = fit_hal(X = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2) , 
                         Y = tempdat_sub2$censor,
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      
      preds_cen = c()
      preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
      preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
      
      tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
      tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
      tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
      tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
      
      ### assign the weight function
      clean_dat$weight = 1
      clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
      clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
      
      #=========================================#
      #== step 2: find eta_opt for cv-hal ======#
      #=========================================#
      
      ### the test time points when computing RMST
      a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
      a = c(0, sort(unique(pmin(a$qal, upp)))) 
      ### the idx of subjects censored or not
      final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
      uncen_idx = which(final_cen$final_cen == 0)
      
      if(est == "IPW"){
        temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx)
      } else{
        temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx, div = div)
      }
      mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt1, cen=mod_cen1, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                             est_var = T, nuisance = "hal")
      results1 = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))
      #=========================================#
      #==== step 3: prepare the under_hal ======#
      #=========================================#
      
      cv_se = sqrt(mod1$var_est2)
      cv_lambda_trt1 = mod_trt1$lambda_star
      coefs_trt1 = mod_trt1$coefs
      cv_lambda_cen1 = mod_cen1$lambda_star
      coefs_cen1 = mod_cen1$coefs
      
      cv_lambda_trt2 = mod_trt2$lambda_star
      coefs_trt2 = mod_trt2$coefs
      cv_lambda_cen2 = mod_cen2$lambda_star
      coefs_cen2 = mod_cen2$coefs
      
      
      #======================================#
      #== step 3.1: for trt1 model ==========#
      #======================================#
      
      tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
      tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
      tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250, ]
      
      preds_trt = c()
      preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
      preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
      
      basis_func_trt1 = as.matrix(cbind(1, mod_trt1$x_basis)) ### n by p, the row 1 is for intercept term
      basis_func_trt2 = as.matrix(cbind(1, mod_trt2$x_basis)) ### n by p, the row 1 is for intercept term
      A1 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID <= 250]
      A2 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID > 250]
      A1_hat = preds_trt[tempdat_trt$ID <= 250]
      A2_hat = preds_trt[tempdat_trt$ID > 250]
      
      
      cutoff_trt1 = cv_se/log(n/2)/sqrt(n/2) ##################################################################OLD UNDERSMOOTHING
      cost_func_trt1 = sum(abs(coefs_trt1) * abs(t(basis_func_trt1) %*% as.matrix(A1 - A1_hat)/n/2))###########HERE
      if(cost_func_trt1 <= cutoff_trt1) {
        under_lambda_trt1 = cv_lambda_trt1
        flag_trt1 = 1
      } else{
        count = 1
        lambda_trt1 = mod_trt1$lambda_star
        cost_func_trt1_old = cost_func_trt1  + 1
        mod_trt1_grid = "good"
        
        while (cost_func_trt1 > cutoff_trt1 & cost_func_trt1 < cost_func_trt1_old) {
          lambda_trt1 = lambda_trt1*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
            return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
          tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
          
          mod_trt1_grid = tryCatch(fit_hal(X = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2), Y = tempdat_sub1$trt, family = "binomial",yolo = FALSE, return_x_basis = T,
                                           cv_select = F, lambda = lambda_trt1), error=function(e) e, warning=function(w) w)
          
          if (is(mod_trt1_grid,"warning")){
            break
          }
          
          preds_trt1_grid = predict(mod_trt1_grid, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
          basis_func_trt1 = as.matrix(cbind(1, mod_trt1_grid$x_basis)) ### n by p, the row 1 is for intercept term
          A1 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID <= 250]
          A1_hat = preds_trt1_grid
          
          cost_func_trt1_old = cost_func_trt1
          cost_func_trt1 = sum(abs(coefs_trt1) * abs(t(basis_func_trt1) %*% as.matrix(A1 - A1_hat)/n/2))
          count = count + 1
        }
        if(is(mod_trt1_grid,"warning") | cost_func_trt1 >= cost_func_trt1_old) { # in this case, cannot attain the boundary
          flag_trt1 = 2
          under_lambda_trt1 = cv_lambda_trt1
        } else{ # in this case, attain the boundary
          flag_trt1 = 0
          under_lambda_trt1 = lambda_trt1
        }
      }
      
      #======================================#
      #== step 3.2: for trt2 model ==========#
      #======================================#
      
      preds_trt = c()
      preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
      preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
      
      basis_func_trt1 = as.matrix(cbind(1, mod_trt1$x_basis)) ### n by p, the row 1 is for intercept term
      basis_func_trt2 = as.matrix(cbind(1, mod_trt2$x_basis)) ### n by p, the row 1 is for intercept term
      A1 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID <= 250]
      A2 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID > 250]
      A1_hat = preds_trt[tempdat_trt$ID <= 250]
      A2_hat = preds_trt[tempdat_trt$ID > 250]
      
      
      cutoff_trt2 = cv_se/log(n/2)/sqrt(n/2) ##################################################################OLD UNDERSMOOTHING
      cost_func_trt2 = sum(abs(coefs_trt2) * abs(t(basis_func_trt2) %*% as.matrix(A2 - A2_hat)/n/2))########### HERE
      if(cost_func_trt2 <= cutoff_trt2) {
        under_lambda_trt2 = cv_lambda_trt2
        flag_trt2 = 1
      } else{
        count = 1
        lambda_trt2 = mod_trt2$lambda_star
        cost_func_trt2_old = cost_func_trt2 + 1
        mod_trt2_grid = "good"
        
        while (cost_func_trt2 > cutoff_trt2 & cost_func_trt2 < cost_func_trt2_old) {
          lambda_trt2 = lambda_trt2*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
            return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
          tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250,]
          
          mod_trt2_grid = tryCatch(fit_hal(X = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2), Y = tempdat_sub2$trt, family = "binomial",yolo = FALSE, return_x_basis = T,
                                           cv_select = F, lambda = lambda_trt2), error=function(e) e, warning=function(w) w)
          
          if (is(mod_trt2_grid,"warning")){
            break
          }
          
          preds_trt2_grid = predict(mod_trt2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
          basis_func_trt2 = as.matrix(cbind(1, mod_trt2_grid$x_basis)) ### n by p, the row 1 is for intercept term
          A2 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID > 250]
          A2_hat = preds_trt2_grid
          
          cost_func_trt2_old = cost_func_trt2
          cost_func_trt2 = sum(abs(coefs_trt2) * abs(t(basis_func_trt2) %*% as.matrix(A2 - A2_hat)/n/2))
          count = count + 1
        }
        if(is(mod_trt2_grid,"warning") | cost_func_trt2 >= cost_func_trt2_old) { # in this case, cannot attain the boundary
          flag_trt2 = 2
          under_lambda_trt2 = cv_lambda_trt2
        } else{ # in this case, attain the boundary
          flag_trt2 = 0
          under_lambda_trt2 = lambda_trt2
        }
      }
      
      
      #======================================#
      #== step 3.3: for cen1 model ==========#
      #======================================#
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
      tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
      tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      preds_cen = c()
      preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
      preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
      
      basis_func_cen1 = as.matrix(cbind(1, mod_cen1$x_basis)) ### n by p, the row 1 is for intercept term
      C1 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID <= 250]
      C1_hat = preds_cen[tempdat_cen$ID <= 250]
      
      cutoff_cen1 = cv_se/log(n/2)/sqrt(n/2)##################################################################OLD UNDERSMOOTHING
      cost_func_cen1 = sum(abs(coefs_cen1) * abs(t(basis_func_cen1) %*% as.matrix(C1 - C1_hat)/n/2))##########HERE
      
      if(cost_func_cen1 <= cutoff_cen1) {
        under_lambda_cen1 = cv_lambda_cen1
        flag_cen1 = 1
      } else{
        count = 1
        lambda_cen1 = mod_cen1$lambda_star
        cost_func_cen1_old = cost_func_cen1  + 1
        mod_cen1_grid = "good"
        
        while (cost_func_cen1 > cutoff_cen1 & cost_func_cen1 < cost_func_cen1_old) {
          lambda_cen1 = lambda_cen1*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          mod_cen1_grid = tryCatch(fit_hal(X = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2), Y = tempdat_sub1$censor, family = "binomial", yolo = FALSE, return_x_basis = T,
                                           cv_select = F, lambda = lambda_cen1), error=function(e) e, warning=function(w) w)
          
          if (is(mod_cen1_grid,"warning")){
            break
          }
          
          preds_cen1_grid =  predict(mod_cen1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
          basis_func_cen1 = as.matrix(cbind(1, mod_cen1_grid$x_basis)) ### n by p, the row 1 is for intercept term
          C1 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID <= 250]
          C1_hat = preds_cen1_grid
          
          cost_func_cen1_old = cost_func_cen1
          cost_func_cen1 = sum(abs(coefs_cen1) * abs(t(basis_func_cen1) %*% as.matrix(C1 - C1_hat)/n/2))
          count = count + 1
        }
        if(is(mod_cen1_grid,"warning") | cost_func_cen1 >= cost_func_cen1_old) { # in this case, cannot attain the boundary
          flag_cen1 = 2
          under_lambda_cen1 = cv_lambda_cen1
        } else{
          flag_cen1 = 0
          under_lambda_cen1 = lambda_cen1
        }
      }
      
      
      #======================================#
      #== step 3.4: for cen2 model ==========#
      #======================================#
      tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
      tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
      tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
      
      preds_cen = c()
      preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
      preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
      
      basis_func_cen2 = as.matrix(cbind(1, mod_cen2$x_basis)) ### n by p, the row 1 is for intercept term
      C2 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID > 250]
      C2_hat = preds_cen[tempdat_cen$ID > 250]
      
      cutoff_cen2 = cv_se/log(n/2)/sqrt(n/2)##################################################################OLD UNDERSMOOTHING
      cost_func_cen2 = sum(abs(coefs_cen2) * abs(t(basis_func_cen2) %*% as.matrix(C2 - C2_hat)/n/2))##########HERE
      
      if(cost_func_cen2 <= cutoff_cen2) {
        under_lambda_cen2 = cv_lambda_cen2
        flag_cen2 = 1
      } else{
        count = 1
        lambda_cen2 = mod_cen2$lambda_star
        cost_func_cen2_old = cost_func_cen2 + 1
        mod_cen2_grid = "good"
        
        while (cost_func_cen2 > cutoff_cen2 & cost_func_cen2 < cost_func_cen2_old) {
          lambda_cen2 = lambda_cen2*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          mod_cen2_grid = tryCatch(fit_hal(X = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2), Y = tempdat_sub2$censor, family = "binomial", yolo = FALSE, return_x_basis = T, 
                                           cv_select = F, lambda = lambda_cen2),
                                   error=function(e) e, warning=function(w) w)
          
          if (is(mod_cen2_grid,"warning")){
            break
          }
          
          preds_cen2_grid =  predict(mod_cen2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
          basis_func_cen2 = as.matrix(cbind(1, mod_cen2_grid$x_basis)) ### n by p, the row 1 is for intercept term
          C2 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID > 250]
          C2_hat = preds_cen2_grid
          
          cost_func_cen2_old = cost_func_cen2
          cost_func_cen2 = sum(abs(coefs_cen2) * abs(t(basis_func_cen2) %*% as.matrix(C2 - C2_hat)/n/2))
          count = count + 1
        }
        if(is(mod_cen2_grid,"warning") | cost_func_cen2 >= cost_func_cen2_old) { # in this case, cannot attain the boundary
          flag_cen2 = 2
          under_lambda_cen2 = cv_lambda_cen2
        } else{
          flag_cen2 = 0
          under_lambda_cen2 = lambda_cen2
        }
      }
      
      
      #======================================#
      #=step 4: re-estimation with under_hal=#
      #======================================#
      if((upp==36)&(est=='IPW')) {tole=0.4
      } else if((upp==36)&(est!='IPW')) {tole=0.5
      } else {tole=0}
      
      if(flag_cen1 != 0 & flag_trt1!=0 & flag_cen2 != 0 & flag_trt2!=0){ # This MC data doesn't change between CV and UNDER HAL!
        mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt1, cen=mod_cen1, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                               est_var = T, nuisance = "hal", tole=tole)
        if(all(!is.na(mod1))){results = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))}
      }
      else{
        tempdat = clean_dat[clean_dat$stage > 1, ]
        tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
          return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
        tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
        tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250, ]
        mod_trt1 = fit_hal(X = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2) , 
                           Y = tempdat_sub1$trt, 
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_trt1)
        mod_trt2 = fit_hal(X = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2) , 
                           Y = tempdat_sub2$trt, 
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_trt2)
        preds_trt = c()
        preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
        preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
        
        tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
        tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
        tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
        tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
        
        tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
        tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
        tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
        mod_cen1 = fit_hal(X = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2) , 
                           Y = tempdat_sub1$censor, 
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_cen1)
        mod_cen2 = fit_hal(X = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2) , 
                           Y = tempdat_sub2$censor,
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_cen2)
        
        preds_cen = c()
        preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(1, tempdat_sub1$x1, tempdat_sub1$x2))
        preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(1, tempdat_sub2$x1, tempdat_sub2$x2))
        
        tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
        tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
        tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
        tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
        
        ### assign the weight function
        clean_dat$weight = 1
        clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
        clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
        
        ### the test time points when computing RMST
        a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
        a = c(0, sort(unique(pmin(a$qal, upp)))) 
        ### the idx of subjects censored or not
        final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
        uncen_idx = which(final_cen$final_cen == 0)
        
        if(est == "IPW"){
          temp2 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx)
        } else{
          temp2 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx, div = div)
        }
        mod2 = loss_indicator2(n, temp2$par, dat, trt=mod_trt1, cen=mod_cen1, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                               est_var = T, nuisance = "hal", tole=tole)
        if(all(!is.na(mod2))){results = list(eta_opt = temp2$par, est = temp2$val, se = sqrt(mod2$var_est2))}
        }
    }
  }
  else if(nuisance == "rf"){
    library(randomForest)
    tempdat = clean_dat[clean_dat$stage > 1, ]
    tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
    
    tempdat$trt <- as.factor(tempdat$trt)
    mod_trt <- randomForest(formula = trt~x1+x2, data=tempdat, subset = tempdat$selvar_trt == 1, ntree = 1000, nodesize = 10)
    preds_trt = predict(mod_trt, type="prob")
    # notice that the prediction of random forest has two columns
    preds_trt = predict(mod_trt, type = "prob")
    
    tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0, 1]
    tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1, 2]
    tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
    tempdat$p.denominator_trt[tempdat$selvar_trt == 1] <- pmin(tempdat$p.denominator_trt[tempdat$selvar_trt == 1], 0.99)
    tempdat$p.denominator_trt[tempdat$selvar_trt == 1] <- pmax(tempdat$p.denominator_trt[tempdat$selvar_trt == 1], 0.01)
    tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
    tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
    
    tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
    
    tempdat$cen <- as.factor(tempdat$cen)
    mod_cen <- randomForest(formula = cen~x1+x2, data=tempdat, subset = tempdat$selvar_cen == 1, ntree = 1000, nodesize = 10)
    preds_cen =  predict(mod_cen, type="prob")
    
    tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0, 1]
    tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1, 2]
    tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
    tempdat$p.denominator_cen[tempdat$selvar_cen == 1] <- pmin(tempdat$p.denominator_cen[tempdat$selvar_cen == 1], 0.99)
    tempdat$p.denominator_cen[tempdat$selvar_cen == 1] <- pmax(tempdat$p.denominator_cen[tempdat$selvar_cen == 1], 0.01)
    tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
    tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
    
    ### assign the weight function
    clean_dat$weight = 1
    clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
    clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
    
    ### the test time points when computing RMST
    a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
    a = c(0, sort(unique(pmin(a$qal, upp)))) 
    ### the idx of subjects censored or not
    final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
    uncen_idx = which(final_cen$final_cen == 0)
    
    if(est == "IPW"){
      temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx)
    } else{
      temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx, div = div)
    }
    mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt, cen=mod_cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", est_var = T, nuisance = "hal")
    
    results = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))
  }
  return(results)
}


### simulation function when logit model is mis-specified
### this function works for all simulation scenarios for logit and rf, and only n = 250 cases for HAL
optimal_sim2 <- function(n = 250, dat, nuisance = "logit", est = "IPW", upp){
  
  ### clean the data
  dat$start = dat$stage - 1
  clean_dat = dat[dat$x1 != 0,]
  KK = max(dat$stage)
  div = ifelse(KK > 10, 2, 1)
  
  if(nuisance == "logit"){
    ### fit logit model for treatment and censoring
    trt <- ipwtm_qal(exposure = trt, family = "binomial", link = "logit",
                     numerator = NULL, denominator = ~ z1+z2,
                     id = ID, tstart = start, timevar = stage, type = "first",
                     data = clean_dat[clean_dat$stage > 1, ])
    
    cen <- ipwtm_qal(exposure = censor, family = "binomial", link = "logit",
                     numerator = NULL, denominator = ~ z1+z2,
                     id = ID, tstart = start, timevar = stage, type = "first",
                     data = clean_dat[clean_dat$stage > 1, ]) 
    
    ### assign the weight function
    clean_dat$weight = 1
    clean_dat$weight[clean_dat$stage > 1] = cen$ipw.weights * trt$ipw.weights
    clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
    
    ### the influence function part for logistics regression, which will be used in the se estimation
    ### assign the gradient function of dP_a,i/d\theta and dP_c,i/d\beta:
    clean_dat$grad_a = 0
    clean_nohead = clean_dat[clean_dat$stage > 1,]
    clean_nohead$grad_a[trt$selvar==1 & clean_nohead$trt == 1] = 1 - trt$assign.prob[trt$selvar==1 & clean_nohead$trt == 1]
    clean_nohead$grad_a[trt$selvar==1 & clean_nohead$trt == 0] = trt$assign.prob[trt$selvar==1 & clean_nohead$trt == 0] - 1
    clean_dat$grad_a[clean_dat$stage > 1] = clean_nohead$grad_a
    
    clean_dat$grad_c = 0
    clean_nohead = clean_dat[clean_dat$stage > 1,]
    clean_nohead$grad_c[cen$selvar==1 & clean_nohead$censor == 0] = cen$assign.prob[cen$selvar==1 & clean_nohead$censor == 0] - 1
    clean_dat$grad_c[clean_dat$stage > 1] = clean_nohead$grad_c
    
    ### the test time points when computing RQAL
    a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
    a = c(0, sort(unique(pmin(a$qal, upp+1)))) 
    ### the idx of subjects censored or not
    final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
    uncen_idx = which(final_cen$final_cen == 0)
    
    if(est == "IPW"){
      temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx)
    } else{
      temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx, div = div)
    }
    tole = ifelse((upp==36)&(est!='IPW'), 0.4, 0)
    mod1 = loss_indicator2(n, temp1$par, dat, trt=trt, cen=cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", est_var = T, 
                           nuisance = "logit",tole=tole)
    
    results = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))
  }
  else if(nuisance == "hal"){
    rate = 0.95
    #==============================================================#
    #================== step 1: cv-hal results ====================#
    #==============================================================#
    if(n == 250){
      library(hal9001)
      
      tempdat = clean_dat[clean_dat$stage > 1, ]
      tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      mod_trt = fit_hal(X = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_trt==1,] , Y = tempdat$trt[tempdat$selvar_trt==1], family = "binomial", yolo = FALSE, return_x_basis = T)
      
      preds_trt = predict(mod_trt, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_trt==1,])
      
      tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
      tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
      tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
      tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      mod_cen = fit_hal(X = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_cen==1,] , Y = tempdat$censor[tempdat$selvar_cen==1], family = "binomial", yolo = FALSE, return_x_basis = T)
      preds_cen = predict(mod_cen, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_cen==1,])
      
      tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
      tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
      tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
      tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
      
      ### assign the weight function
      clean_dat$weight = 1
      clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
      clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
      
      ### the test time points when computing RMST
      a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
      a = c(0, sort(unique(pmin(a$qal, upp)))) 
      ### the idx of subjects censored or not
      final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
      uncen_idx = which(final_cen$final_cen == 0)
      
      
      if(est == "IPW"){
        temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx)
      } else{
        temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx, div = div)
      }
      mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt, cen=mod_cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", est_var = T, nuisance = "hal")
      cv_se = sqrt(mod1$var_est2)
      cv_lambda_trt = mod_trt$lambda_star
      coefs_trt = mod_trt$coefs 
      cv_lambda_cen = mod_cen$lambda_star
      coefs_cen = mod_cen$coefs 
      
      #==============================================================#
      #============= step 2: undersmooth-hal results ================#
      #==============================================================#
      
      
      #======================================#
      #== step 2.1: for trt model ===========#
      #======================================#
      
      preds_trt = predict(mod_trt, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_trt==1,])
      basis_func_trt = as.matrix(cbind(1, mod_trt$x_basis)) ### n by p, the row 1 is for intercept term
      A = tempdat$trt[tempdat$selvar_trt==1]
      A_hat = preds_trt
      
      cutoff_trt = cv_se/log(n)/sqrt(n)##################################################################OLD UNDERSMOOTHING
      cost_func_trt = sum(abs(coefs_trt) * abs(t(basis_func_trt) %*% as.matrix(A - A_hat)/n))############HERE
      if(cost_func_trt <= cutoff_trt) {
        under_lambda_trt = cv_lambda_trt
        flag_trt = 1
      } else{
        count = 1
        lambda_trt = mod_trt$lambda_star
        cost_func_trt_old = cost_func_trt + 1
        mod_trt_grid = "good"
        
        while (cost_func_trt > cutoff_trt & cost_func_trt < cost_func_trt_old) {
          lambda_trt = lambda_trt*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          mod_trt_grid = tryCatch(fit_hal(X = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_trt==1,] , Y = tempdat$trt[tempdat$selvar_trt==1], family = "binomial", yolo = FALSE, return_x_basis = T,
                                          cv_select = F, lambda = lambda_trt),error=function(e) e, warning=function(w) w)
          
          if (is(mod_trt_grid,"warning")){
            break
          }
          
          preds_trt_grid = predict(mod_trt_grid, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_trt==1,])
          basis_func_trt = as.matrix(cbind(1, mod_trt_grid$x_basis)) ### n by p, the row 1 is for intercept term
          A = tempdat$trt[tempdat$selvar_trt==1]
          A_hat = preds_trt_grid
          
          cost_func_trt_old = cost_func_trt
          cost_func_trt = sum(abs(coefs_trt) * abs(t(basis_func_trt) %*% as.matrix(A - A_hat)/n))
          count = count + 1
        }
        if(is(mod_trt_grid,"warning") | cost_func_trt >= cost_func_trt_old) { # in this case, cannot attain the boundary 
          flag_trt = 2
          under_lambda_trt = cv_lambda_trt
        } else{ # in this case, attain the boundary
          flag_trt = 0
          under_lambda_trt = lambda_trt
        }
      }
      
      #======================================#
      #== step 2.2: for cen model ===========#
      #======================================#
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      preds_cen = predict(mod_cen, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_cen==1,])
      basis_func_cen = as.matrix(cbind(1, mod_cen$x_basis)) ### n by p, the row 1 is for intercept term
      C = tempdat$censor[tempdat$selvar_cen==1]
      C_hat = preds_cen
      
      cutoff_cen = cv_se/log(n)/sqrt(n)##################################################################OLD UNDERSMOOTHING
      cost_func_cen = sum(abs(coefs_cen) * abs(t(basis_func_cen) %*% as.matrix(C - C_hat)/n))############HERE
      
      if(cost_func_cen <= cutoff_cen) {
        under_lambda_cen = cv_lambda_cen
        flag_cen = 1
      } else{
        count = 1
        lambda_cen = mod_cen$lambda_star
        cost_func_cen_old = cost_func_cen + 1
        mod_cen_grid = "good"
        
        while (cost_func_cen > cutoff_cen & cost_func_cen < cost_func_cen_old) {
          lambda_cen = lambda_cen*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          mod_cen_grid = tryCatch(fit_hal(X = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_cen==1,] , Y = tempdat$censor[tempdat$selvar_cen==1], family = "binomial", yolo = FALSE, return_x_basis = T,
                                          cv_select = F, lambda = lambda_cen),error=function(e) e, warning=function(w) w)
          
          if (is(mod_cen_grid,"warning")){
            break
          }
          
          preds_cen_grid = predict(mod_cen_grid, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_cen==1,])
          basis_func_cen = as.matrix(cbind(1, mod_cen_grid$x_basis)) ### n by p, the row 1 is for intercept term
          C = tempdat$censor[tempdat$selvar_cen==1]
          C_hat = preds_cen_grid
          
          cost_func_cen_old = cost_func_cen
          cost_func_cen = sum(abs(coefs_cen) * abs(t(basis_func_cen) %*% as.matrix(C - C_hat)/n))
          count = count + 1
        }
        if(is(mod_cen_grid,"warning") | cost_func_cen >= cost_func_cen_old) { # in this case, cannot attain the boundary 
          flag_cen = 2
          under_lambda_cen = cv_lambda_cen
        } else{
          flag_cen = 0
          under_lambda_cen = lambda_cen
        }
      }
      
      
      #======================================#
      #== step 2.3: re-estimate the results =#
      #======================================#
      if(flag_cen != 0 & flag_trt!=0){
        results = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))
      } else{
        tempdat = clean_dat[clean_dat$stage > 1, ]
        tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        mod_trt = fit_hal(X = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_trt==1,] , Y = tempdat$trt[tempdat$selvar_trt==1], family = "binomial", yolo = FALSE, cv_select = F, lambda = under_lambda_trt)
        preds_trt = predict(mod_trt, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_trt==1,])
        
        tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
        tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
        tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
        tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
        
        tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        mod_cen = fit_hal(X = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_cen==1,] , Y = tempdat$censor[tempdat$selvar_cen==1], family = "binomial", yolo = FALSE, cv_select = F, lambda = under_lambda_cen)
        preds_cen = predict(mod_cen, new_data = cbind(tempdat$z1, tempdat$z2)[tempdat$selvar_cen==1,])
        
        tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
        tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
        tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
        tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
        
        ### assign the weight function
        clean_dat$weight = 1
        clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
        clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
        
        ### the test time points when computing RMST
        a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
        a = c(0, sort(unique(pmin(a$qal, upp)))) 
        ### the idx of subjects censored or not
        final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
        uncen_idx = which(final_cen$final_cen == 0)
        
        if(est == "IPW"){
          temp2 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx)
        } else{
          temp2 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx, div = div)
        }
        
        mod2 = loss_indicator2(n, temp2$par, dat, trt=mod_trt, cen=mod_cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", est_var = T, nuisance = "hal")
        results = list(eta_opt = temp2$par, est = temp2$val, se = sqrt(mod2$var_est2))
      }
    }
    else if(n == 500){ ### beyond the computation limit of hal, and we need to split the hal model into two parts
      
      #=========================================#
      #==== step 1: fit the cv hal model =======#
      #=========================================#
      
      library(hal9001)
      tempdat = clean_dat[clean_dat$stage > 1, ]
      tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
        return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
      tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
      tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250, ]
      mod_trt1 = fit_hal(X = cbind(tempdat_sub1$z1, tempdat_sub1$z2) , 
                         Y = tempdat_sub1$trt, 
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      mod_trt2 = fit_hal(X = cbind(tempdat_sub2$z1, tempdat_sub2$z2) , 
                         Y = tempdat_sub2$trt, 
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      preds_trt = c()
      preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
      preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
      
      tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
      tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
      tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
      tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
      tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
      tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
      mod_cen1 = fit_hal(X = cbind(tempdat_sub1$z1, tempdat_sub1$z2) , 
                         Y = tempdat_sub1$censor, 
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      mod_cen2 = fit_hal(X = cbind(tempdat_sub2$z1, tempdat_sub2$z2) , 
                         Y = tempdat_sub2$censor,
                         family = "binomial",
                         yolo = FALSE, return_x_basis = T)
      
      preds_cen = c()
      preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
      preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
      
      tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
      tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
      tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
      tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
      tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
      tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
      
      ### assign the weight function
      clean_dat$weight = 1
      clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
      clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
      
      #=========================================#
      #== step 2: find eta_opt for cv-hal ======#
      #=========================================#
      
      ### the test time points when computing RMST
      a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
      a = c(0, sort(unique(pmin(a$qal, upp)))) 
      ### the idx of subjects censored or not
      final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
      uncen_idx = which(final_cen$final_cen == 0)
      
      if(est == "IPW"){
        temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx)
      } else{
        temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                        BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                        uncen_idx = uncen_idx, div = div)
      }
      mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt1, cen=mod_cen1, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                             est_var = T, nuisance = "hal")
      results1 = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))
      #=========================================#
      #==== step 3: prepare the under_hal ======#
      #=========================================#
      
      cv_se = sqrt(mod1$var_est2)
      cv_lambda_trt1 = mod_trt1$lambda_star
      coefs_trt1 = mod_trt1$coefs
      cv_lambda_cen1 = mod_cen1$lambda_star
      coefs_cen1 = mod_cen1$coefs
      
      cv_lambda_trt2 = mod_trt2$lambda_star
      coefs_trt2 = mod_trt2$coefs
      cv_lambda_cen2 = mod_cen2$lambda_star
      coefs_cen2 = mod_cen2$coefs
      
      
      #======================================#
      #== step 3.1: for trt1 model ==========#
      #======================================#
      
      tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
      tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
      tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250, ]
      
      preds_trt = c()
      preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
      preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
      
      basis_func_trt1 = as.matrix(cbind(1, mod_trt1$x_basis)) ### n by p, the row 1 is for intercept term
      basis_func_trt2 = as.matrix(cbind(1, mod_trt2$x_basis)) ### n by p, the row 1 is for intercept term
      A1 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID <= 250]
      A2 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID > 250]
      A1_hat = preds_trt[tempdat_trt$ID <= 250]
      A2_hat = preds_trt[tempdat_trt$ID > 250]
      
      
      cutoff_trt1 = cv_se/log(n/2)/sqrt(n/2) ##################################################################OLD UNDERSMOOTHING
      cost_func_trt1 = sum(abs(coefs_trt1) * abs(t(basis_func_trt1) %*% as.matrix(A1 - A1_hat)/n/2))###########HERE
      if(cost_func_trt1 <= cutoff_trt1) {
        under_lambda_trt1 = cv_lambda_trt1
        flag_trt1 = 1
      } else{
        count = 1
        lambda_trt1 = mod_trt1$lambda_star
        cost_func_trt1_old = cost_func_trt1  +1
        mod_trt1_grid = "good"
        
        while (cost_func_trt1 > cutoff_trt1 & cost_func_trt1 < cost_func_trt1_old) {
          lambda_trt1 = lambda_trt1*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
            return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
          tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
          
          mod_trt1_grid = tryCatch(fit_hal(X = cbind(tempdat_sub1$z1, tempdat_sub1$z2), Y = tempdat_sub1$trt, family = "binomial",yolo = FALSE, return_x_basis = T,
                                           cv_select = F, lambda = lambda_trt1), error=function(e) e, warning=function(w) w)
          
          if (is(mod_trt1_grid,"warning")){
            break
          }
          
          preds_trt1_grid = predict(mod_trt1_grid, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
          basis_func_trt1 = as.matrix(cbind(1, mod_trt1_grid$x_basis)) ### n by p, the row 1 is for intercept term
          A1 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID <= 250]
          A1_hat = preds_trt1_grid
          
          cost_func_trt1_old = cost_func_trt1
          cost_func_trt1 = sum(abs(coefs_trt1) * abs(t(basis_func_trt1) %*% as.matrix(A1 - A1_hat)/n/2))
          count = count + 1
        }
        if(is(mod_trt1_grid,"warning") | cost_func_trt1 >= cost_func_trt1_old) { # in this case, cannot attain the boundary
          flag_trt1 = 2
          under_lambda_trt1 = cv_lambda_trt1
        } else{ # in this case, attain the boundary
          flag_trt1 = 0
          under_lambda_trt1 = lambda_trt1
        }
      }
      
      #======================================#
      #== step 3.2: for trt2 model ==========#
      #======================================#
      
      preds_trt = c()
      preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
      preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
      
      basis_func_trt1 = as.matrix(cbind(1, mod_trt1$x_basis)) ### n by p, the row 1 is for intercept term
      basis_func_trt2 = as.matrix(cbind(1, mod_trt2$x_basis)) ### n by p, the row 1 is for intercept term
      A1 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID <= 250]
      A2 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID > 250]
      A1_hat = preds_trt[tempdat_trt$ID <= 250]
      A2_hat = preds_trt[tempdat_trt$ID > 250]
      
      
      cutoff_trt2 = cv_se/log(n/2)/sqrt(n/2)##################################################################OLD UNDERSMOOTHING
      cost_func_trt2 = sum(abs(coefs_trt2) * abs(t(basis_func_trt2) %*% as.matrix(A2 - A2_hat)/n/2))##########HERE
      if(cost_func_trt2 <= cutoff_trt2) {
        under_lambda_trt2 = cv_lambda_trt2
        flag_trt2 = 1
      } else{
        count = 1
        lambda_trt2 = mod_trt2$lambda_star
        cost_func_trt2_old = cost_func_trt2 + 1
        mod_trt2_grid = "good"
        
        while (cost_func_trt2 > cutoff_trt2 & cost_func_trt2 < cost_func_trt2_old) {
          lambda_trt2 = lambda_trt2*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
            return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
          tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250,]
          
          mod_trt2_grid = tryCatch(fit_hal(X = cbind(tempdat_sub2$z1, tempdat_sub2$z2), Y = tempdat_sub2$trt, family = "binomial",yolo = FALSE, return_x_basis = T,
                                           cv_select = F, lambda = lambda_trt2), error=function(e) e, warning=function(w) w)
          
          if (is(mod_trt2_grid,"warning")){
            break
          }
          
          preds_trt2_grid = predict(mod_trt2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
          basis_func_trt2 = as.matrix(cbind(1, mod_trt2_grid$x_basis)) ### n by p, the row 1 is for intercept term
          A2 = tempdat$trt[tempdat$selvar_trt==1][tempdat_trt$ID > 250]
          A2_hat = preds_trt2_grid
          
          cost_func_trt2_old = cost_func_trt2
          cost_func_trt2 = sum(abs(coefs_trt2) * abs(t(basis_func_trt2) %*% as.matrix(A2 - A2_hat)/n/2))
          count = count + 1
        }
        if(is(mod_trt2_grid,"warning") | cost_func_trt2 >= cost_func_trt2_old) { # in this case, cannot attain the boundary
          flag_trt2 = 2
          under_lambda_trt2 = cv_lambda_trt2
        } else{ # in this case, attain the boundary
          flag_trt2 = 0
          under_lambda_trt2 = lambda_trt2
        }
      }
      
      
      #======================================#
      #== step 3.3: for cen1 model ==========#
      #======================================#
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
      tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
      tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
      
      tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
      
      preds_cen = c()
      preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
      preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
      
      basis_func_cen1 = as.matrix(cbind(1, mod_cen1$x_basis)) ### n by p, the row 1 is for intercept term
      C1 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID <= 250]
      C1_hat = preds_cen[tempdat_cen$ID <= 250]
      
      cutoff_cen1 = cv_se/log(n/2)/sqrt(n/2)##################################################################OLD UNDERSMOOTHING
      cost_func_cen1 = sum(abs(coefs_cen1) * abs(t(basis_func_cen1) %*% as.matrix(C1 - C1_hat)/n/2))##########HERE
      
      if(cost_func_cen1 <= cutoff_cen1) {
        under_lambda_cen1 = cv_lambda_cen1
        flag_cen1 = 1
      } else{
        count = 1
        lambda_cen1 = mod_cen1$lambda_star
        cost_func_cen1_old = cost_func_cen1  + 1
        mod_cen1_grid = "good"
        
        while (cost_func_cen1 > cutoff_cen1 & cost_func_cen1 < cost_func_cen1_old) {
          lambda_cen1 = lambda_cen1*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          mod_cen1_grid = tryCatch(fit_hal(X = cbind(tempdat_sub1$z1, tempdat_sub1$z2), Y = tempdat_sub1$censor, family = "binomial", yolo = FALSE, return_x_basis = T,
                                           cv_select = F, lambda = lambda_cen1), error=function(e) e, warning=function(w) w)
          
          if (is(mod_cen1_grid,"warning")){
            break
          }
          
          preds_cen1_grid =  predict(mod_cen1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
          basis_func_cen1 = as.matrix(cbind(1, mod_cen1_grid$x_basis)) ### n by p, the row 1 is for intercept term
          C1 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID <= 250]
          C1_hat = preds_cen1_grid
          
          cost_func_cen1_old = cost_func_cen1
          cost_func_cen1 = sum(abs(coefs_cen1) * abs(t(basis_func_cen1) %*% as.matrix(C1 - C1_hat)/n/2))
          count = count + 1
        }
        if(is(mod_cen1_grid,"warning") | cost_func_cen1 >= cost_func_cen1_old) { # in this case, cannot attain the boundary
          flag_cen1 = 2
          under_lambda_cen1 = cv_lambda_cen1
        } else{
          flag_cen1 = 0
          under_lambda_cen1 = lambda_cen1
        }
      }
      
      
      #======================================#
      #== step 3.4: for cen2 model ==========#
      #======================================#
      tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
      tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
      tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
      
      preds_cen = c()
      preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
      preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
      
      basis_func_cen2 = as.matrix(cbind(1, mod_cen2$x_basis)) ### n by p, the row 1 is for intercept term
      C2 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID > 250]
      C2_hat = preds_cen[tempdat_cen$ID > 250]
      
      cutoff_cen2 = cv_se/log(n/2)/sqrt(n/2)##################################################################OLD UNDERSMOOTHING
      cost_func_cen2 = sum(abs(coefs_cen2) * abs(t(basis_func_cen2) %*% as.matrix(C2 - C2_hat)/n/2))##########HERE
      
      if(cost_func_cen2 <= cutoff_cen2) {
        under_lambda_cen2 = cv_lambda_cen2
        flag_cen2 = 1
      } else{
        count = 1
        lambda_cen2 = mod_cen2$lambda_star
        cost_func_cen2_old = cost_func_cen2 + 1
        mod_cen2_grid = "good"
        
        while (cost_func_cen2 > cutoff_cen2 & cost_func_cen2 < cost_func_cen2_old) {
          lambda_cen2 = lambda_cen2*rate
          
          tempdat = clean_dat[clean_dat$stage > 1, ]
          tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
          
          mod_cen2_grid = tryCatch(fit_hal(X = cbind(tempdat_sub2$z1, tempdat_sub2$z2), Y = tempdat_sub2$censor, family = "binomial", yolo = FALSE, return_x_basis = T, 
                                           cv_select = F, lambda = lambda_cen2),
                                   error=function(e) e, warning=function(w) w)
          
          if (is(mod_cen2_grid,"warning")){
            break
          }
          
          preds_cen2_grid =  predict(mod_cen2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
          basis_func_cen2 = as.matrix(cbind(1, mod_cen2_grid$x_basis)) ### n by p, the row 1 is for intercept term
          C2 = tempdat$censor[tempdat$selvar_cen==1][tempdat_cen$ID > 250]
          C2_hat = preds_cen2_grid
          
          cost_func_cen2_old = cost_func_cen2
          cost_func_cen2 = sum(abs(coefs_cen2) * abs(t(basis_func_cen2) %*% as.matrix(C2 - C2_hat)/n/2))
          count = count + 1
        }
        if(is(mod_cen2_grid,"warning") | cost_func_cen2 >= cost_func_cen2_old) { # in this case, cannot attain the boundary
          flag_cen2 = 2
          under_lambda_cen2 = cv_lambda_cen2
        } else{
          flag_cen2 = 0
          under_lambda_cen2 = lambda_cen2
        }
      }
      
      
      #======================================#
      #=step 4: re-estimation with under_hal=#
      #======================================#
      tole = ifelse((upp==36)&(est!='IPW'), 0.51, 0)
      if(flag_cen1 != 0 & flag_trt1!=0 & flag_cen2 != 0 & flag_trt2!=0){ # This MC data doesn't change between CV and UNDER HAL!
        mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt1, cen=mod_cen1, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                               est_var = T, nuisance = "hal", tole = tole)
        var_1 = ifelse(all(!is.na(mod1)), NA, sqrt(mod1$var_est2))
        results = list(eta_opt = temp1$par, est = temp1$val, se = var_1)
      }
      else{
        tempdat = clean_dat[clean_dat$stage > 1, ]
        tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) 
          return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        tempdat_trt = tempdat[tempdat$selvar_trt == 1, ]
        tempdat_sub1 = tempdat_trt[tempdat_trt$ID <= 250,]
        tempdat_sub2 = tempdat_trt[tempdat_trt$ID > 250, ]
        mod_trt1 = fit_hal(X = cbind(tempdat_sub1$z1, tempdat_sub1$z2) , 
                           Y = tempdat_sub1$trt, 
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_trt1)
        mod_trt2 = fit_hal(X = cbind(tempdat_sub2$z1, tempdat_sub2$z2) , 
                           Y = tempdat_sub2$trt, 
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_trt2)
        preds_trt = c()
        preds_trt[tempdat_trt$ID <= 250] = predict(mod_trt1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
        preds_trt[tempdat_trt$ID > 250] =  predict(mod_trt2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
        
        tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- 1 - preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0]
        tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1]
        tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
        tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
        
        tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
        
        tempdat_cen = tempdat[tempdat$selvar_cen == 1, ]
        tempdat_sub1 = tempdat_cen[tempdat_cen$ID <= 250, ]
        tempdat_sub2 = tempdat_cen[tempdat_cen$ID > 250, ]
        mod_cen1 = fit_hal(X = cbind(tempdat_sub1$z1, tempdat_sub1$z2) , 
                           Y = tempdat_sub1$censor, 
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_cen1)
        mod_cen2 = fit_hal(X = cbind(tempdat_sub2$z1, tempdat_sub2$z2) , 
                           Y = tempdat_sub2$censor,
                           family = "binomial",
                           yolo = FALSE, 
                           cv_select = F,
                           lambda = under_lambda_cen2)
        
        preds_cen = c()
        preds_cen[tempdat_cen$ID <= 250] = predict(mod_cen1, new_data = cbind(tempdat_sub1$z1, tempdat_sub1$z2))
        preds_cen[tempdat_cen$ID > 250] =  predict(mod_cen2, new_data = cbind(tempdat_sub2$z1, tempdat_sub2$z2))
        
        tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
        tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- 1 - preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0]
        tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1]
        tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
        tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
        tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
        
        ### assign the weight function
        clean_dat$weight = 1
        clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
        clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
        
        ### the test time points when computing RMST
        a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
        a = c(0, sort(unique(pmin(a$qal, upp)))) 
        ### the idx of subjects censored or not
        final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
        uncen_idx = which(final_cen$final_cen == 0)
        
        if(est == "IPW"){
          temp2 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx)
        } else{
          temp2 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                          BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                          uncen_idx = uncen_idx, div = div)
        }
        mod2 = loss_indicator2(n, temp2$par, dat, trt=mod_trt1, cen=mod_cen1, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", 
                               est_var = T, nuisance = "hal", tole=tole)
        var_2 = ifelse(all(!is.na(mod2)), NA, sqrt(mod2$var_est2))
        results = list(eta_opt = temp2$par, est = temp2$val, se = var_2)
        
      }
    }
  }
  else if(nuisance == "rf"){
    library(randomForest)
    tempdat = clean_dat[clean_dat$stage > 1, ]
    tempdat$selvar_trt <- do.call("c", lapply(split(tempdat$trt, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
    
    tempdat$trt <- as.factor(tempdat$trt)
    mod_trt <- randomForest(formula = trt~z1+z2, data=tempdat, subset = tempdat$selvar_trt == 1, ntree = 1000, nodesize = 10)
    preds_trt = predict(mod_trt, type="prob")
    # notice that the prediction of random forest has two columns
    preds_trt = predict(mod_trt, type = "prob")
    
    tempdat$p.denominator_trt <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator_trt[tempdat$trt == 0 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 0, 1]
    tempdat$p.denominator_trt[tempdat$trt == 1 & tempdat$selvar_trt == 1] <- preds_trt[tempdat$trt[tempdat$selvar_trt == 1] == 1, 2]
    tempdat$p.denominator_trt[tempdat$selvar_trt == 0] <- 1
    tempdat$p.denominator_trt[tempdat$selvar_trt == 1] <- pmin(tempdat$p.denominator_trt[tempdat$selvar_trt == 1], 0.99)
    tempdat$p.denominator_trt[tempdat$selvar_trt == 1] <- pmax(tempdat$p.denominator_trt[tempdat$selvar_trt == 1], 0.01)
    tempdat$w.denominator_trt <- unlist(lapply(split(tempdat$p.denominator_trt, tempdat$ID), function(x)cumprod(x)))
    tempdat$ipw.weights_trt <- 1/tempdat$w.denominator_trt
    
    tempdat$selvar_cen <- do.call("c", lapply(split(tempdat$censor, tempdat$ID),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))
    
    tempdat$cen <- as.factor(tempdat$cen)
    mod_cen <- randomForest(formula = cen~z1+z2, data=tempdat, subset = tempdat$selvar_cen == 1, ntree = 1000, nodesize = 10)
    preds_cen =  predict(mod_cen, type="prob")
    
    tempdat$p.denominator_cen <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator_cen[tempdat$censor == 0 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 0, 1]
    tempdat$p.denominator_cen[tempdat$censor == 1 & tempdat$selvar_cen == 1] <- preds_cen[tempdat$censor[tempdat$selvar_cen == 1] == 1, 2]
    tempdat$p.denominator_cen[tempdat$selvar_cen == 0] <- 1
    tempdat$p.denominator_cen[tempdat$selvar_cen == 1] <- pmin(tempdat$p.denominator_cen[tempdat$selvar_cen == 1], 0.99)
    tempdat$p.denominator_cen[tempdat$selvar_cen == 1] <- pmax(tempdat$p.denominator_cen[tempdat$selvar_cen == 1], 0.01)
    tempdat$w.denominator_cen <- unlist(lapply(split(tempdat$p.denominator_cen, tempdat$ID), function(x)cumprod(x)))
    tempdat$ipw.weights_cen <- 1/tempdat$w.denominator_cen
    
    ### assign the weight function
    clean_dat$weight = 1
    clean_dat$weight[clean_dat$stage > 1] = tempdat$ipw.weights_cen * tempdat$ipw.weights_trt
    clean_dat = clean_dat %>% group_by(ID) %>% mutate(qal = cumsum(qol * outcome)) %>% as.data.frame
    
    ### the test time points when computing RMST
    a = dat %>% group_by(ID) %>% summarise(qal = sum(qol * outcome))
    a = c(0, sort(unique(pmin(a$qal, upp)))) 
    ### the idx of subjects censored or not
    final_cen = clean_dat %>% group_by(ID) %>% summarise(final_cen = ifelse(sum(censor) > 0, 1, 0))
    uncen_idx = which(final_cen$final_cen == 0)
    
    if(est == "IPW"){
      temp1 <- genoud(fn=indicator_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx)
    } else{
      temp1 <- genoud(fn=smooth_opt1, nvars=3, default.domains=1, starting.values=rep(0, 3), max=TRUE, print.level=1, 
                      BFGS=FALSE, optim.method="Nelder-Mead", P9=0, unif.seed=1107, int.seed=0130, clean_dat = clean_dat, a = a,
                      uncen_idx = uncen_idx, div = div)
    }
    mod1 = loss_indicator2(n, temp1$par, dat, trt=mod_trt, cen=mod_cen, clean_dat, a=a, final_cen=final_cen, upp = upp, time, target = "rmst", est_var = T, nuisance = "hal")
    
    results = list(eta_opt = temp1$par, est = temp1$val, se = sqrt(mod1$var_est2))
  }
}


### the function used for logistic regression; modified from the ipwtm function from ipw package
ipwtm_qal <- function(
  exposure,
  family,
  link,
  numerator = NULL,
  denominator,
  id,
  tstart,
  timevar,
  type,
  data,
  corstr = "ar1",
  trunc = NULL,
  ...)
{
  #save input
  tempcall <- match.call()
  #some basic input checks
  if (!("exposure" %in% names(tempcall))) stop("No exposure variable specified")
  if (!("family" %in% names(tempcall)) | ("family" %in% names(tempcall) & !(tempcall$family %in% c("binomial", "survival", "multinomial", "ordinal", "gaussian")))) stop("No valid family specified (\"binomial\", \"survival\", \"multinomial\", \"ordinal\", \"gaussian\")")
  if (tempcall$family == "binomial") {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
  if (tempcall$family == "ordinal" ) {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}
  if (!("denominator" %in% names(tempcall))) stop("No denominator model specified")
  if (!is.null(tempcall$numerator) & !is(eval(tempcall$numerator), "formula")) stop("Invalid numerator formula specified")
  if (!is.null(tempcall$denominator) & !is(eval(tempcall$denominator), "formula")) stop("Invalid denominator formula specified")
  if (!("id" %in% names(tempcall))) stop("No patient id specified")
  if (tempcall$family == "survival" & !("tstart" %in% names(tempcall))) stop("No tstart specified, is necessary for family = \"survival\"")
  if (!("timevar" %in% names(tempcall))) stop("No timevar specified")
  if (!("type" %in% names(tempcall))) stop("No type specified (\"first\" or \"all\")")
  if (!(tempcall$type %in% c("first", "all"))) stop("No type specified (\"first\" or \"all\")")
  if (tempcall$family %in% c("survival", "multinomial", "ordinal") & tempcall$type == "all") stop(paste("Type \"all\" not yet implemented for family = ", deparse(tempcall$family, width.cutoff = 500), sep = ""))
  if (tempcall$family %in% c("gaussian") & tempcall$type == "first") stop(paste("Type \"first\" not implemented for family = ", deparse(tempcall$family, width.cutoff = 500), sep = ""))
  if (tempcall$family %in% c("gaussian") & !("numerator" %in% names(tempcall))) stop("Numerator necessary for family = \"gaussian\"")
  if (!("data" %in% names(tempcall))) stop("No data specified")
  if (!is.null(tempcall$trunc)) {if(tempcall$trunc < 0 | tempcall$trunc > 0.5) stop("Invalid truncation percentage specified (0-0.5)")}
  #record original order of dataframe so that the output can be returned in the same order
  order.orig <- 1:nrow(data)
  order.orig <- order.orig[order(
    eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
    eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
  )] #sort as below
  #sort dataframe on follow-up time within each individual, necessary for cumulative products below
  data <- data[order(
    eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
    eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
  ),]
  #make new dataframe for newly computed variables, to prevent variable name conflicts
  tempdat <- data.frame(
    id = data[,as.character(tempcall$id)],
    timevar = data[,as.character(tempcall$timevar)],
    exposure = data[,as.character(tempcall$exposure)]
  )
  #make selection variable, time points up to first switch from lowest value, or all time points
  if (type == "first" & (family == "binomial" | family == "survival"))
  {tempdat$selvar <- do.call("c", lapply(split(tempdat$exposure, tempdat$id),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))}
  if (type == "first" & (family == "multinomial" | family == "ordinal")){
    z <- unique(tempdat$exposure)[unique(tempdat$exposure) != sort(unique(tempdat$exposure))[1]]
    min2 <- function(x)ifelse(min(is.na(unique(x))) == 1, NA, min(x, na.rm = TRUE))
    tempdat$selvar <- do.call("c", lapply(split(tempdat$exposure, tempdat$id),function(x)if (!is.na(min2(match(z, x)))) return(c(rep(1,min2(match(z, x))),rep(0,length(x)-min2(match(z, x))))) else return(rep(1,length(x)))))
  }
  if (type == "all")
  {tempdat$selvar <- rep(1, nrow(tempdat))}
  #weights binomial, type "first"
  if (tempcall$family == "binomial" & tempcall$type == "first") {
    if(tempcall$link == "logit") lf <- binomial(link = logit)
    if(tempcall$link == "probit") lf  <- binomial(link = probit)
    if(tempcall$link == "cauchit") lf  <- binomial(link = cauchit)
    if(tempcall$link == "log") lf  <- binomial(link = log)
    if(tempcall$link == "cloglog") lf  <- binomial(link = cloglog)
    if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
    else {
      mod1 <- glm(
        formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        family = lf,
        data = data,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        ...)
      tempdat$p.numerator <- vector("numeric", nrow(tempdat))
      tempdat$p.numerator[tempdat$exposure == 0 & tempdat$selvar == 1] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 0]
      tempdat$p.numerator[tempdat$exposure == 1 & tempdat$selvar == 1] <- predict.glm(mod1, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      tempdat$w.numerator <- unlist(lapply(split(tempdat$p.numerator, tempdat$id), function(x)cumprod(x)))
      mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$family <- tempcall$link
      mod1$call$data <- tempcall$data
      mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    }
    mod2 <- glm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      family = lf,
      data = data,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      ...)
    tempdat$p.denominator <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator[tempdat$exposure == 0 & tempdat$selvar == 1] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 0]
    tempdat$p.denominator[tempdat$exposure == 1 & tempdat$selvar == 1] <- predict.glm(mod2, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    tempdat$w.denominator <- unlist(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)))
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$family <- tempcall$link
    mod2$call$data <- tempcall$data
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights binomial, type "all"
  if (tempcall$family == "binomial" & tempcall$type == "all") {
    if(tempcall$link == "logit") lf <- binomial(link = logit)
    if(tempcall$link == "probit") lf  <- binomial(link = probit)
    if(tempcall$link == "cauchit") lf  <- binomial(link = cauchit)
    if(tempcall$link == "log") lf  <- binomial(link = log)
    if(tempcall$link == "cloglog") lf  <- binomial(link = cloglog)
    if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
    else {
      mod1 <- glm(
        formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        family = lf,
        data = data,
        na.action = na.fail,
        ...)
      tempdat$p.numerator <- vector("numeric", nrow(tempdat))
      tempdat$p.numerator[tempdat$exposure == 0] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure == 0]
      tempdat$p.numerator[tempdat$exposure == 1] <- predict.glm(mod1, type = "response")[tempdat$exposure == 1]
      tempdat$w.numerator <- unlist(lapply(split(tempdat$p.numerator, tempdat$id), function(x)cumprod(x)))
      mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$family <- tempcall$link
      mod1$call$data <- tempcall$data
    }
    mod2 <- glm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      family = lf,
      data = data,
      na.action = na.fail,
      ...)
    tempdat$p.denominator <- vector("numeric", nrow(tempdat))
    tempdat$p.denominator[tempdat$exposure == 0] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure == 0]
    tempdat$p.denominator[tempdat$exposure == 1] <- predict.glm(mod2, type = "response")[tempdat$exposure == 1]
    tempdat$w.denominator <- unlist(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)))
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$family <- tempcall$link
    mod2$call$data <- tempcall$data
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights Cox
  if (tempcall$family == "survival") {
    if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
    else {
      mod1 <- coxph(
        formula = eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        data = data,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        method = "efron",
        ...)
      bh <- basehaz(mod1, centered = TRUE)
      temp <- data.frame(timevar = sort(unique(tempdat$timevar)))
      temp <- merge(temp, data.frame(timevar = bh$time, bashaz.cum.numerator = bh$hazard), by = "timevar", all.x = TRUE);rm(bh)
      if (is.na(temp$bashaz.cum.numerator[temp$timevar == min(unique(tempdat$timevar))])) temp$bashaz.cum.numerator[temp$timevar == min(unique(tempdat$timevar))] <- 0
      temp$bashaz.cum.numerator <- approx(x = temp$timevar, y = temp$bashaz.cum.numerator, xout = temp$timevar, method = "constant", rule = 2)$y
      temp$bashaz.numerator[1] <- temp$bashaz.cum.numerator[1]
      temp$bashaz.numerator[2:nrow(temp)] <- diff(temp$bashaz.cum.numerator, 1)
      temp$bashaz.cum.numerator <- NULL
      tempdat <- merge(tempdat, temp, by = "timevar", all.x = TRUE);rm(temp)
      tempdat <- tempdat[order(tempdat$id, tempdat$timevar),]
      tempdat$risk.numerator[tempdat$selvar == 1] <- predict(mod1, type="risk", centered = TRUE)
      tempdat$hazard.numerator[tempdat$selvar == 1] <- with(tempdat[tempdat$selvar == 1,], bashaz.numerator*risk.numerator)
      tempdat$p.numerator[with(tempdat, selvar == 1 & exposure == 0)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 0),], exp(-1*bashaz.numerator*risk.numerator))
      tempdat$p.numerator[with(tempdat, selvar == 1 & exposure == 1)] <- 1 - with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.numerator*risk.numerator))
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      tempdat$w.numerator <- unsplit(lapply(split(tempdat$p.numerator, tempdat$id), function(x) cumprod(x)), tempdat$id)
      mod1$call$formula <- eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
    }
    mod2 <- coxph(
      formula = eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      method = "efron",
      ...)
    bh <- basehaz(mod2, centered = TRUE)
    temp <- data.frame(timevar = sort(unique(tempdat$timevar)))
    temp <- merge(temp, data.frame(timevar = bh$time, bashaz.cum.denominator = bh$hazard), by = "timevar", all.x = TRUE);rm(bh)
    if (is.na(temp$bashaz.cum.denominator[temp$timevar == min(unique(tempdat$timevar))])) temp$bashaz.cum.denominator[temp$timevar == min(unique(tempdat$timevar))] <- 0
    temp$bashaz.cum.denominator <- approx(x = temp$timevar, y = temp$bashaz.cum.denominator, xout = temp$timevar, method = "constant", rule = 2)$y
    temp$bashaz.denominator[1] <- temp$bashaz.cum.denominator[1]
    temp$bashaz.denominator[2:nrow(temp)] <- diff(temp$bashaz.cum.denominator, 1)
    temp$bashaz.cum.denominator <- NULL
    tempdat <- merge(tempdat, temp, by = "timevar", all.x = TRUE);rm(temp)
    tempdat <- tempdat[order(tempdat$id, tempdat$timevar),]
    tempdat$risk.denominator[tempdat$selvar == 1] <- predict(mod2, type="risk", centered = TRUE)
    tempdat$hazard.denominator[tempdat$selvar == 1] <- with(tempdat[tempdat$selvar == 1,], bashaz.denominator*risk.denominator)
    tempdat$p.denominator[with(tempdat, selvar == 1 & exposure == 0)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 0),], exp(-1*bashaz.denominator*risk.denominator))
    tempdat$p.denominator[with(tempdat, selvar == 1 & exposure == 1)] <- 1 - with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.denominator*risk.denominator))
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    tempdat$w.denominator <- unsplit(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)), tempdat$id)
    mod2$call$formula <- eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
  }
  #weights multinomial
  if (tempcall$family == "multinomial") {
    if (is.null(tempcall$numerator)) tempdat$p.numerator <- 1
    else {
      mod1 <- multinom(
        formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        data = data,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        ...)
      pred1 <- as.data.frame(predict(mod1, type = "probs"))
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      for (i in 1:length(unique(tempdat$exposure)))tempdat$p.numerator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
      mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
      mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    }
    mod2 <- multinom(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      ...)
    pred2 <- as.data.frame(predict(mod2, type = "probs"))
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    for (i in 1:length(unique(tempdat$exposure)))tempdat$p.denominator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, p.numerator/p.denominator), tempdat$id), function(x)cumprod(x)), tempdat$id)
  }
  #weights ordinal
  if (tempcall$family == "ordinal") {
    if(tempcall$link == "logit") m <- "logistic"
    if(tempcall$link == "probit") m  <- "probit"
    if(tempcall$link == "cloglog") m  <- "cloglog"
    if(tempcall$link == "cauchit") m  <- "cauchit"
    if (is.null(tempcall$numerator)) tempdat$p.numerator <- 1
    else {
      mod1 <- polr(
        formula = eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
        data = data,
        method = m,
        subset = tempdat$selvar == 1,
        na.action = na.fail,
        ...)
      pred1 <- as.data.frame(predict(mod1, type = "probs"))
      tempdat$p.numerator[tempdat$selvar == 0] <- 1
      for (i in 1:length(unique(tempdat$exposure)))tempdat$p.numerator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
      mod1$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
      mod1$call$data <- tempcall$data
      mod1$call$method <- m
      mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    }
    mod2 <- polr(
      formula = eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      method = m,
      subset = tempdat$selvar == 1,
      na.action = na.fail,
      ...)
    pred2 <- as.data.frame(predict(mod2, type = "probs"))
    tempdat$p.denominator[tempdat$selvar == 0] <- 1
    for (i in 1:length(unique(tempdat$exposure)))tempdat$p.denominator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
    mod2$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$method <- m
    mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
    tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, p.numerator/p.denominator), tempdat$id), function(x)cumprod(x)), tempdat$id)
  }
  #weights gaussian
  if (tempcall$family == "gaussian") {
    mod1 <- geeglm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
      data = data,
      id = tempdat$id,
      corstr = tempcall$corstr,
      waves = tempdat$timevar,
      ...)
    mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
    mod1$call$data <- tempcall$data
    mod1$call$id <- tempcall$id
    mod1$call$corstr <- tempcall$corstr
    mod1$call$waves <- tempcall$waves
    mod2 <- geeglm(
      formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
      data = data,
      id = tempdat$id,
      corstr = tempcall$corstr,
      waves = tempdat$timevar,
      ...)
    mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
    mod2$call$data <- tempcall$data
    mod2$call$id <- tempcall$id
    mod2$call$corstr <- tempcall$corstr
    mod2$call$waves <- tempcall$waves
    tempdat$kdens1 <- dnorm(tempdat$exposure, predict(mod1), as.numeric(sqrt(summary(mod1)$dispersion)[1]))
    tempdat$kdens2 <- dnorm(tempdat$exposure, predict(mod2), as.numeric(sqrt(summary(mod2)$dispersion)[1]))
    tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, kdens1/kdens2), tempdat$id), function(x)cumprod(x)), tempdat$id)
  }
  #check for NA's in weights
  if (sum(is.na(tempdat$ipw.weights)) > 0) stop ("NA's in weights!")
  #truncate weights, when trunc value is specified (0-0.5)
  if (!(is.null(tempcall$trunc))){
    tempdat$weights.trunc <- tempdat$ipw.weights
    tempdat$weights.trunc[tempdat$ipw.weights <= quantile(tempdat$ipw.weights, 0+trunc)] <- quantile(tempdat$ipw.weights, 0+trunc)
    tempdat$weights.trunc[tempdat$ipw.weights >  quantile(tempdat$ipw.weights, 1-trunc)] <- quantile(tempdat$ipw.weights, 1-trunc)
  }
  #return results in the same order as the original input dataframe
  if (is.null(tempcall$trunc)){
    if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], assign.prob = tempdat$p.denominator, call = tempcall, selvar = tempdat$selvar[order(order.orig)], den.mod = mod2))
    else return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], assign.prob = tempdat$p.denominator, call = tempcall, selvar = tempdat$selvar[order(order.orig)], num.mod = mod1, den.mod = mod2))
  }
  else{
    if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], weights.trunc = tempdat$weights.trunc[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], den.mod = mod2))
    else return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], weights.trunc = tempdat$weights.trunc[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], num.mod = mod1, den.mod = mod2))
  }
}


#############################################################################################################################
#NEW UNDERSMOOTHING FUNCTION DEFINED
#############################################################################################################################
newUndersmooth <- function(X, Y, A, family = 'binomial', cv_lambda, lamGridLength = 40, rate = .9){
  nUnder <- nrow(X)
  
  #Create grid of lambdas
  lamGrid <- cv_lambda*cumprod(c(1, rep(rate, lamGridLength)))
  
  #Fit undersmoothed hals
  under_fit = fit_hal(X = X, Y = Y, family = family, 
                      yolo = FALSE, return_x_basis = T, lambda = lamGrid, fit_control = list(cv_select = FALSE))
  
  #Extract necessary pieces
  underFit_coefs <- under_fit$coefs
  underFit_basis <- cbind(1, as.matrix(under_fit$x_basis))
  underFit_preds <- predict(under_fit, new_data = X)
  
  #Loop over lambdas and return undersmoothing criterion
  criterion <- sapply(1:length(lamGrid), function(i){
    #Get values for this particular lambda
    lamCoefs <- underFit_coefs[,i]
    lamPreds <- underFit_preds[,i]
    
    #l1 norm of coefs
    l1Norm <- sum(abs(lamCoefs))
    
    #Undersmoothing fraction
    lamFrac <- (A - lamPreds)/(l1Norm*lamPreds)
    
    #Get lamFrac times each basis vector with a non-zero coefficent
    lamBasisFrac <- t(underFit_basis[,lamCoefs != 0]) %*% lamFrac
    
    #Take absolute value, sum and return
    return(sum(abs(lamBasisFrac)))
  })
  
  #Loop over lambdas and return number of non-zero features
  numFeat <- sapply(1:length(lamGrid), function(i){
    lamCoefs <- underFit_coefs[,i]
    return(sum(lamCoefs != 0))
  })
  
  #Select only criterion for which numFeatures <= n^{1/2}
  viableCrit <- criterion[numFeat <= sqrt(nUnder)]
  
  #Select chosen criterion
  lamSel <- which(viableCrit == min(viableCrit))
  
  #Select lambda
  underLam <- lamGrid[lamSel]
  
  return(list(halFit = under_fit, preds = underFit_preds[,lamSel], underLam = underLam, crits = criterion, r = lamSel))
}
