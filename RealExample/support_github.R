### created at 05/09/2019 clean and re-organize the code for chapter 3

#==========================================================================================#
#==================== the main functions ==================================================#
#==========================================================================================#

#******************************************************************************#
# indicator_opt1_als; smooth_opt2_als                                          #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  eta            a given coefficient vector for the dynamic regime            #
#  index          the column index of covariates considered in dynamic         #
#                 in clean_dat                                                 #
#******************************************************************************#
### these two inputs should  have the same length

### this is a simplified wrap of the main function: loss_indicator2_als()
# the reason I have it is because this function only include policy estimation and the optimation
# of this function is faster

### modify it for flexible covariates
# eta is a vector of coefficients, index is the index for covariates in optimal regime estimation in clean_dat
indicator_opt1_als <- function(eta, index){
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
  #   }
  td = diff(a)
  return(sum(td * pp))
}

### this is the function for bias-corrected estimation
smooth_opt2_als <- function(eta, index){
  n_dat = clean_dat
  stage = S
  pp = c()

  ### construct the design matrix for optiml regime (here, we include the intercept and interested time-varying covariates)
  regime_design = cbind(1, as.matrix(n_dat[, index]))

  ### construct the trt_last record
  n_dat$trt_last = 0
  n_dat$trt_last[which(n_dat$stage > 1)] = n_dat$trt[which(n_dat$stage > 1)-1]

  score = regime_design %*% as.matrix(eta)
  #score = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2)
  s = sd(score)
  s = ifelse(s > 0, s, eta[1]>=0)
  score_mono = score^(1 - n_dat$trt_last)
  #score_mono = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2)^(1 - n_dat$trt_last)
  bw = n^(-1/3) * s/stage

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


#******************************************************************************#
# loss_indicator2_als                                                          #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  n              the number of subjects                                       #
#  eta            a given coefficient vector for the dynamic regime            #
#  clean_dat      the input design matrix                                      #
#  upp            the upper limit of RQAL                                      #
#  time           the interested time when target = "surv                      #
#  target         "surv" for survival probability; "rmst" for RQAL             #
#  est_var        if do standard error estimation                              #
#  index          the column index of covariates considered in dynamic         #
#                 in clean_dat                                                 #
#******************************************************************************#

### this is the mian function. It can be used to do policy evaluatiion given regime estimation "eta"; and calculate the sd when
# eat_var = T
loss_indicator2_als <- function(n, eta, clean_dat, upp = 45, time, target = "rmst", est_var = F, index){
  n_dat = clean_dat
  ### the un-smoothed surv-prob estimation
  pp = pp2 = c()

  ### the un-smoothed variance estimation
  var_est = var_est2 = c()

  ### construct the design matrix for optiml regime (here, we include the intercept and interested time-varying covariates)
  regime_design = cbind(1, as.matrix(n_dat[, index]))

  ### construct the follow indicator
  id = clean_dat$ID

  trt_eta_raw  = (regime_design %*% as.matrix(eta)) > 0
  #trt_eta_raw  = (eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0
  trt_eta = unlist(lapply(1:length(unique(id)), function(x) pmin(cumsum(trt_eta_raw[id == x]), 1)))
  n_dat$follow = n_dat$trt * trt_eta_raw + (1-n_dat$trt) * (1 - trt_eta_raw)
  #n_dat$follow = n_dat$trt * ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0) +
  #  (1-n_dat$trt) * (1 - ((eta[1] + eta[2] * n_dat$x1 + eta[3] * n_dat$x2) > 0))
  n_dat$follow[clean_dat$stage == 1] = 1
  n_dat = n_dat %>% group_by(ID) %>% mutate(follow_cum = cumprod(follow)) %>% as.data.frame

  ### the target is surv-prob
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


      if(est_var == T){
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

      if(est_var == T){
        ### calculate the variance estimation: 12/02/2018
        # step 1: calculate the gradient function
        effect_id = which(final_weight$f_weight > 0)

        if(length(effect_id) > 0){
          # grad_a_theta = lapply(effect_id, function(i) {d = n_dat[n_dat$ID == i, ][1:idx[i], ];
          # apply(sweep(cbind(1, d$x1, d$x2), MARGIN=1, d$grad_a, `*`),2,sum)})
          # grad_a_theta = do.call("rbind", grad_a_theta)
          # #grad_theta1 = apply(sweep(grad_a_theta, MARGIN=1, -final_weight$f_weight[effect_id], `*`),2,sum)/n
          # grad_theta2 = apply(sweep(grad_a_theta, MARGIN=1,
          #                           final_weight$f_weight[effect_id]*(final_weight2$f_weight[effect_id] - sum(final_weight2$f_weight)), `*`),2,sum)/(sum(final_weight2$f_weight)^2)
          #
          # grad_c_beta = lapply(effect_id, function(i) {d = n_dat[n_dat$ID == i, ][1:idx[i], ];
          # apply(sweep(cbind(1, d$x1, d$x2), MARGIN=1, d$grad_c, `*`),2,sum)})
          # grad_c_beta = do.call("rbind", grad_c_beta)
          # #grad_beta1 = apply(sweep(grad_c_beta, MARGIN=1, -final_weight$f_weight[effect_id], `*`),2,sum)/n
          # grad_beta2 = apply(sweep(grad_c_beta, MARGIN=1,
          #                          final_weight$f_weight[effect_id]*(final_weight2$f_weight[effect_id] - sum(final_weight2$f_weight)), `*`),2,sum)/(sum(final_weight2$f_weight)^2)

          # step 3: empirical variance estimation
          #IC_eta = final_weight$f_weight - pp[i]
          #IC_theta = Score_a_sub %*% grad_theta1
          #IC_beta = Score_c_sub %*% grad_beta1

          IC_eta2 = (final_weight$f_weight - final_weight2$f_weight*pp2[i])/sum(final_weight2$f_weight)*n
          IC_theta2 = 0
          IC_beta2 = 0

          #IC = IC_eta + IC_theta + IC_beta
          IC2 = IC_eta2 + IC_theta2 + IC_beta2

          #IC_mat[,i] = IC
          IC_mat2[,i] = IC2
        }
        else{
          IC_mat2[,i] = 0
        }
      }
    }
    td = diff(a)
    #IC = apply(sweep(IC_mat, MARGIN=2, td, `*`),1,sum)
    IC2 = apply(sweep(IC_mat2, MARGIN=2, td, `*`),1,sum)
    return(list(rqal2 = sum(td * pp2), var_est2 = sum(IC2^2)/n^2))
  }
}

#==========================================================================================#
#================================ support functions =======================================#
#==========================================================================================#

### the function to calculate parametric longtidudinal weights: modified from the ipwtm function from ipw package
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
