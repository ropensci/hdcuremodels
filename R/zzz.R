c_stat_b <-
  function(b_u_hat=NULL, itct_hat, b_p_hat, testing_delta, testing_time, X_u=NULL, X_p){
    C_csw_num = 0
    C_csw_denom = 0
    testing_n = length(testing_time)
    if(all(b_p_hat==0)) {
      if(is.null(X_u)) xb = rep(itct_hat, testing_n) else xb = itct_hat+X_u %*% b_u_hat
    } else{
      if(is.null(X_u)) xb = itct_hat + X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0]
      else xb = itct_hat+X_u %*% b_u_hat + X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0]
    }
    for(i in 1:testing_n)
      for(j in 1:testing_n){
        if (j==i | !testing_delta[i] | testing_time[i]>testing_time[j]) next
        I_ij = testing_time[i]<testing_time[j] | (testing_time[i]==testing_time[j] & !testing_delta[j])
        if (!I_ij) next
        if (xb[i]>xb[j]) C_csw_num = C_csw_num + 1
        C_csw_denom = C_csw_denom + 1
      }
    return(C_csw_num / C_csw_denom)
  }

c_stat_beta <-
  function(beta_u_hat=NULL,  beta_p_hat, testing_delta, testing_time, W_u=NULL, W_p){
    C_csw_num = 0
    C_csw_denom = 0
    testing_n = length(testing_time)
    if(all(beta_p_hat==0)) {
      if(is.null(W_u)) W_beta = rep(0, testing_n) else W_beta = W_u %*% beta_u_hat
    } else{
      if(is.null(W_u)) W_beta = W_p[,beta_p_hat!=0,drop=FALSE] %*% beta_p_hat[beta_p_hat!=0]
      else W_beta = W_u %*% beta_u_hat + W_p[,beta_p_hat!=0,drop=FALSE] %*% beta_p_hat[beta_p_hat!=0]
    }
    for(i in 1:testing_n)
      for(j in 1:testing_n){
        if (j==i | !testing_delta[i] | testing_time[i]>testing_time[j]) next
        I_ij = testing_time[i]<testing_time[j] | (testing_time[i]==testing_time[j] & !testing_delta[j])
        if (!I_ij) next
        if (W_beta[i]>W_beta[j]) C_csw_num = C_csw_num + 1
        C_csw_denom = C_csw_denom + 1
      }
    return(C_csw_num / C_csw_denom)
  }



cox_l1 <-
  function(X_u=NULL, X_p=NULL, W_u=NULL, W_p=NULL, time, delta, mu_inc, mu_lat, inits=NULL,
           nIter=100, tol = 1e-4)
  {
    # mu: penalty parameter
    # tol: difference between log-likelihood
    N = length(time)
    J = ncol(X_p) # number of penalized incidence covariates
    M = ncol(W_p) # number of penalized latency covariates

    CAP = 10

    event_time = time[delta==1]
    uniq_event_time = unique(event_time)
    n0 = length(uniq_event_time)
    I0 = matrix(time, ncol = 1) %*% matrix(1, 1, n0) >=
      matrix(1, N, 1) %*% matrix(uniq_event_time, nrow = 1)    # N*n0
    d_j = as.numeric(table(event_time)[rank(unique(event_time))])   #  number of events at time T_j
    T_n0 = max(event_time)   # last event
    tail_ind = which(time>=T_n0 & delta==0)

    ######## initialization ########
    step = 1
    if (!is.null(X_p)) b_p <- rep(0,J) # KJA added to allow no penalized incidence variables
    if (!is.null(W_p)) beta_p <- rep(0,M) # KJA added added to allow no penalized latency variables
    if(is.null(inits)) inits = initialization_cox(X_u, W_u, time, delta)
    itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u; survprob = inits$survprob
    if(is.null(X_u)) b_x = rep(itct,N) else b_x = itct+X_u %*% b_u
    if(is.null(X_u)) { # KJA changed to allow no penalized incidence variables
      pen_fac_inc = rep(1, ncol(X_p))
    } else if (is.null(X_p)) {
      pen_fac_inc = rep(0, ncol(X_u))
    } else {
      pen_fac_inc = c(rep(0, ncol(X_u)), rep(1, ncol(X_p)))
    }
    if(is.null(W_u)) {
      pen_fac_lat = rep(1, ncol(W_p))
    } else if (is.null(W_p)) {
      pen_fac_lat = rep(0, ncol(W_u))
    } else {
      pen_fac_lat = c(rep(0, ncol(W_u)), rep(1, ncol(W_p)))
    }
    pir = rep(1, N)
    llp1_0 <- llp2_0 <- 0
    lik_inc <- lik_lat <- c()
    b_p_path <- beta_p_path <- b_u_path <- beta_u_path <- itct_path <- NULL
    conv1 <- conv2 <- FALSE

    ######## loop ########

    repeat{

      #### E-step
      pi_x = 1/(1+exp(-b_x))
      numerator = pi_x[delta==0] * survprob[delta==0]
      pir[delta==0] = numerator/(1-pi_x[delta==0] + numerator)

      #### M-step
      ### incidence
      if (!conv1){
        fit_inc = glmnet::glmnet(x = cbind(X_u, X_p), y = cbind(1-pir, pir), family = "binomial",
                                 penalty.factor = pen_fac_inc, lambda=mu_inc, standardize = FALSE)
        coef_inc = coef(fit_inc)
        itct = max(-CAP, min(CAP, coef_inc[1]))
        if(is.null(X_u)) {
          b_p = pmax(-CAP, pmin(CAP, coef_inc[-1]))
          b_x = itct + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
        } else if (is.null(X_p)) {
          b_u = pmax(-CAP, pmin(CAP, coef_inc[2:(ncol(X_u)+1)]))
          b_x = itct + X_u %*% b_u
        } else {
          b_u = pmax(-CAP, pmin(CAP, coef_inc[2:(ncol(X_u)+1)]))
          b_p = pmax(-CAP, pmin(CAP, coef_inc[-(1:(ncol(X_u)+1))]))
          b_x = itct + X_u %*% b_u + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
        }
        llp1_1 = sum(pir*b_x - log(1+exp(b_x))) - N * mu_inc *sum(abs(b_p))
      }
      lik_inc = c(lik_inc, llp1_1)

      ### latency
      if (!conv2){
        nz_pir = which(pir>0)
        fit_lat = glmnet::glmnet(x = cbind(W_u, W_p)[nz_pir,],
                                 y = Surv(time = time[nz_pir], event = delta[nz_pir]),
                                 family = "cox",
                                 offset = log(pir[nz_pir]),
                                 penalty.factor = pen_fac_lat,
                                 lambda = mu_lat, standardize = FALSE)
        coef_lat = coef(fit_lat)
        if (is.null(W_u)) {
          beta_p = pmax(-CAP, pmin(CAP, as.numeric(coef_lat)))
          beta_w = W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
        } else if (is.null(W_p)) {
          beta_u = pmax(-CAP, pmin(CAP, coef_lat[1:ncol(W_u)]))
          beta_w = W_u %*% beta_u
        } else {
          beta_u = pmax(-CAP, pmin(CAP, coef_lat[1:ncol(W_u)]))
          beta_p = pmax(-CAP, pmin(CAP, coef_lat[-(1:ncol(W_u))]))
          beta_w = W_u %*% beta_u + W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
        }

        ### survival
        exp_beta_w_p = exp(beta_w) * pir
        denom = matrix(exp_beta_w_p,1) %*% I0   # 1*n0
        hazard = rep(0, N)
        hazard[delta==1] = sapply(event_time, function(x) (d_j/denom)[uniq_event_time==x])
        accum_hazard = sapply(time, function(x) sum((d_j/denom)[uniq_event_time<=x]))
        surv_baseline = exp(-accum_hazard)
        if (length(tail_ind)>0){
          wtail_alpha = uniroot(function(a)
            sum(-exp_beta_w_p*max(accum_hazard)*log(time/T_n0)*(time/T_n0)^a + delta/a + delta*log(time/T_n0)),
            c(0.01,10), extendInt = "yes")$root
          wtail_lambda = (max(accum_hazard))^(1/wtail_alpha)/T_n0
          surv_baseline[tail_ind] = exp(-(wtail_lambda*time[tail_ind])^wtail_alpha)
        }
        survprob = surv_baseline^exp(beta_w)
        llp2_1 = sum(log(hazard[delta==1]) + beta_w[delta==1]) - sum(exp_beta_w_p*accum_hazard) -
          N * mu_lat *sum(abs(beta_p))
      }
      lik_lat = c(lik_lat, llp2_1)

      ## record updated parameters
      itct_path = c(itct_path, itct)
      if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
      if(!is.null(X_p)) b_p_path = rbind(b_p_path, b_p)
      if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)
      if(!is.null(W_p)) beta_p_path = rbind(beta_p_path, beta_p)
      if (!conv1 & abs(llp1_1- llp1_0)< tol) conv1 <- TRUE
      if (!conv2 & abs(llp2_1- llp2_0)< tol) conv2 <- TRUE
      if (step > 1 & (conv1 & conv2) | step >= nIter) { #
        break
      }
      llp1_0 <- llp1_1
      llp2_0 <- llp2_1
      step <- 1 + step
    }

    ######## output ########
    output <- list(b_p_path = b_p_path, beta_p_path = beta_p_path, b_u_path = b_u_path, itct_path = itct_path,
                   beta_u_path = beta_u_path,
                   lik_inc = lik_inc, lik_lat = lik_lat)
    output
  }

cox_mcp_scad <-
  function(X_u=NULL, X_p=NULL, W_u=NULL, W_p=NULL, time, delta, penalty = "MCP",
           mu_inc, mu_lat, gamma_inc = 3, gamma_lat = 3, inits = NULL,
           nIter=1000, tol = 1e-4)
  {
    # lambda: penalty parameter
    # tol: difference between log-likelihood
    N = length(time)
    J = ncol(X_p) # number of penalized incidence covariates
    M = ncol(W_p) # number of penalized latency covariates

    CAP = 20

    event_time = time[delta==1]
    uniq_event_time = unique(event_time)
    n0 = length(uniq_event_time)
    I0 = matrix(time, ncol = 1) %*% matrix(1, 1, n0) >=
      matrix(1, N, 1) %*% matrix(uniq_event_time, nrow = 1)    # N*n0
    d_j = as.numeric(table(event_time)[rank(uniq_event_time)])   #  number of events at time T_j
    Z = sapply(uniq_event_time, function(t) colSums(W_p[time==t & delta==1,,drop=FALSE])) # M*n0
    T_n0 = max(event_time)   # last event
    tail_ind = which(time>=T_n0 & delta==0)

    ######## initialization ########
    step = 1
    if (!is.null(X_p)) b_p <- rep(0,J) # KJA added to allow no penalized incidence variables
    if (!is.null(W_p)) beta_p <- rep(0,M) # KJA added added to allow no penalized latency variables
    if(is.null(inits)) inits = initialization_cox(X_u, W_u, time, delta)
    itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u; survprob = inits$survprob
    if(is.null(X_u)) b_x = rep(itct,N) else b_x = itct+X_u %*% b_u
    if(is.null(W_u)) beta_w = rep(0,N) else beta_w = W_u %*% beta_u
    pir = rep(1, N)
    exp_beta_w_p = exp(beta_w) * pir
    lik_inc <- lik_lat <- c()
    b_p_path <- beta_p_path <- b_u_path <- beta_u_path <- itct_path <- NULL
    llp1_0 <- llp2_0 <- 0
    conv1 <- conv2 <- FALSE
    ######## loop ########

    repeat{

      #### E-step
      pi_x = 1/(1+exp( pmin(CAP, -b_x) ))
      numerator = pi_x[delta==0] * survprob[delta==0]
      pir[delta==0] = numerator/pmax(1e-50,1-pi_x[delta==0] + numerator)

      #### M-step
      ### incidence
      if (!conv1){
        for (l in 1:J){
          bl = b_p[l]
          b_p[l] = max(-CAP, min(CAP, mcp_scad_cd_upd_inc(penalty=penalty, pir, pi_x, X_p[,l], bl,
                                                          gamma_inc, mu_inc)))
          b_x = b_x + X_p[,l,drop=FALSE] %*% (b_p[l] - bl)
          pi_x = 1/(1+exp( pmin(CAP, -b_x) ))
        }
        out_inc = optim(par = c(itct,b_u), fn = mcp_scad_negloglik_inc, gr = mcp_scad_gradient_inc,
                        b_p = b_p, X_u =X_u, X_p = X_p, pir = pir, CAP=CAP, method="BFGS") # incidence
        itct = max(-CAP, min(CAP, out_inc$par[1]))
        if(!is.null(X_u)) b_u = pmax(-CAP, pmin(CAP, out_inc$par[2:(ncol(X_u)+1)]))
        pen = ifelse(penalty=="MCP", mcp_penalty(b_p, gamma_inc, mu_inc),
                     scad_penalty(b_p, gamma_inc, mu_inc))
        llp1_1 = -out_inc$value - N*pen
        if(is.null(X_u)) b_x = itct + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
        else b_x = itct + X_u %*% b_u + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
      }
      lik_inc = c(lik_inc, llp1_1)

      ### latency
      if (!conv2){
        for (l in 1:M){
          betal = beta_p[l]
          beta_p[l] = max(-CAP, min(CAP, mcp_scad_cd_upd_lat(penalty, exp_beta_w_p, W_p[,l], betal,
                                                             Z[l,], I0, d_j, gamma_lat, mu_lat)))
          change = W_p[,l,drop=FALSE] %*% (beta_p[l] - betal)  # N*1
          beta_w = beta_w + change  # N*1
          exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir  # N*1
          l = l + 1
        }
        if(!is.null(W_u)){
          out_lat = optim(par = beta_u, fn = mcp_scad_negloglik_lat, gr = mcp_scad_gradient_lat,
                          beta_p=beta_p, W_u=W_u, W_p=W_p, delta=delta, pir=pir, I0=I0, d_j=d_j,
                          CAP=CAP, method="BFGS") # latency
          beta_u = pmax(-CAP, pmin(CAP, out_lat$par)) # unpenalized incidence
          beta_w = W_u %*% beta_u + W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
          exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir
        }

        ### survival
        denom = matrix(exp_beta_w_p,1) %*% I0   # 1*n0
        hazard = rep(0, N)
        hazard[delta==1] = sapply(event_time, function(x) (d_j/denom)[uniq_event_time==x])
        accum_hazard = sapply(time, function(x) sum((d_j/denom)[uniq_event_time<=x]))
        surv_baseline = exp(-accum_hazard)
        if (length(tail_ind)>0){
          tmp <- try(wtail_alpha <- uniroot(function(a)
            sum(-exp_beta_w_p*max(accum_hazard)*log(time/T_n0)*(time/T_n0)^a + delta/a + delta*log(time/T_n0)),
            c(0.01,10), extendInt = "yes")$root, silent = TRUE)
          if (!inherits(tmp, "try-error")){
            wtail_lambda = (max(accum_hazard))^(1/wtail_alpha)/T_n0
            surv_baseline[tail_ind] = exp(-(wtail_lambda*time[tail_ind])^wtail_alpha)
          }
        }
        survprob = surv_baseline^exp(pmin(CAP, beta_w))
        pen = ifelse(penalty=="MCP", mcp_penalty(beta_p, gamma_lat, mu_lat),
                     scad_penalty(beta_p, gamma_lat, mu_lat))
        llp2_1 = sum(log(hazard[delta==1]) + beta_w[delta==1])
        - sum(exp_beta_w_p*accum_hazard) - N*pen
      }
      lik_lat = c(lik_lat, llp2_1)

      ## record updated parameters
      itct_path = c(itct_path, itct)
      if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
      b_p_path = rbind(b_p_path, b_p)
      if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)
      beta_p_path = rbind(beta_p_path, beta_p)
      # cat("step=", step, "\n")
      if (!conv1 & abs(llp1_1- llp1_0)< tol) conv1 <- TRUE
      if (!conv2 & abs(llp2_1- llp2_0)< tol) conv2 <- TRUE
      if (step > 1 & (conv1 & conv2) | step >= nIter) { #
        break
      }
      llp1_0 <- llp1_1
      llp2_0 <- llp2_1
      step <- 1 + step
    }

    ######## output ########
    output <- list(b_p_path = b_p_path, b_u_path = b_u_path, itct_path = itct_path,
                   beta_p_path = beta_p_path, beta_u_path = beta_u_path,
                   lik_inc = lik_inc, lik_lat = lik_lat)
    output
  }

cure.em <-
  function(X_u=NULL, X_p, W_u=NULL, W_p, time, event, model="cox", penalty="lasso",
           penalty.factor.inc=rep(1,ncol(x.inc)), penalty.factor.lat=rep(1,ncol(x.lat)),
           thresh=1e-03, maxit=ifelse(penalty=="lasso", 100, 1000),
           lambda.inc=0.1, lambda.lat=0.1, gamma.inc=3, gamma.lat=3, inits=NULL){
    x.inc <- cbind(X_u, X_p)
    x.lat <- cbind(W_u, W_p)
    if (is.null(dim(x.inc)) | is.null(dim(x.lat)))
      stop("x.inc and x.lat should be matrices with 2 or more columns")
    if((ncol(x.inc) <= 1) | (ncol(x.lat) <= 1))
      stop("x.inc and x.lat should be matrices with 2 or more columns")
    if (nrow(x.inc) != nrow(x.lat) | nrow(x.lat) != length(time) | length(time)!= length(event))
      stop("Input dimension mismatch")
    if(class(x.inc)[1] == "data.frame" | class(x.lat)[1] == "data.frame"){
      x.inc = as.matrix(x.inc); x.lat = as.matrix(x.lat)
    }
    if(!model%in%c("cox","weibull","exponential"))
      stop("Only 'cox', 'weibull', 'exponential' available for 'model' parameter")
    if(!penalty%in%c("lasso","MCP","SCAD"))
      stop("Only 'lasso', 'MCP', 'SCAD' available for 'penalty' parameter")
    if(any(!c(penalty.factor.inc, penalty.factor.lat)%in%c(0,1)))
      stop("Penalty factors can only contain 0 or 1")
    if(any(c(lambda.inc, lambda.lat, gamma.inc, gamma.lat)<=0))
      stop("Penalty pamameters lambda and gamma should be positive")
    if(!is.null(inits))
      inits = inits_check(model, N=length(time), penalty.factor.inc, penalty.factor.lat, inits)
    if(model != "cox" & penalty != "lasso")    penalty = "lasso"
    if(model=="cox" & penalty=="lasso")
      fit = cox_l1(X_u, X_p, W_u, W_p, time, event, lambda.inc, lambda.lat, inits, maxit, thresh)
    else if(model=="cox" & penalty %in% c("MCP","SCAD"))
      fit = cox_mcp_scad(X_u, X_p, W_u, W_p, time, event, penalty, lambda.inc, lambda.lat,
                         gamma.inc, gamma.lat, inits, maxit, thresh)
    else if(model=="weibull")
      fit = weib_EM(X_u, X_p, W_u, W_p, time, event, lambda.inc, lambda.lat, inits, maxit, thresh)
    else if(model=="exponential")
      fit = exp_EM(X_u, X_p, W_u, W_p, time, event, lambda.inc, lambda.lat, inits, maxit, thresh)
    model_select = which.max(fit$lik_inc+fit$lik_lat)
    b = rep(NA, ncol(x.inc))
    beta = rep(NA, ncol(x.lat))
    b[penalty.factor.inc==0] = fit$b_u_path[model_select,]
    b[penalty.factor.inc==1] = fit$b_p_path[model_select,]
    beta[penalty.factor.lat==0] = fit$beta_u_path[model_select,]
    beta[penalty.factor.lat==1] = fit$beta_p_path[model_select,]
    output = list(b0=fit$itct_path[model_select], b=b, beta=beta)
    if(model %in% c("exponential","weibull"))
      output$rate = fit$lambda_path[model_select]
    if(model == "weibull")
      output$alpha = fit$alpha_path[model_select]
    output$logLik.inc = fit$lik_inc[model_select]
    output$logLik.lat = fit$lik_lat[model_select]
    return(output)
  }

cv.em.fdr <-
  function(X_u, X_p, W_u, W_p, time, delta, model=c("weibull","exponential"),
           penalty=c("lasso","MCP","SCAD"), fdr=0.2, thresh=1e-3,
           nIter=ifelse(penalty=="lasso", 100, 1000), penalty.factor.inc, penalty.factor.lat,
           grid.tuning = FALSE, lambda.inc.list=NULL, lambda.lat.list=NULL,
           nlambda.inc = ifelse(grid.tuning, 10, 50),
           nlambda.lat = ifelse(grid.tuning, 10, 50),
           lambda.min.ratio.inc = 0.1, lambda.min.ratio.lat = 0.1,
           gamma.inc=3, gamma.lat=3,
           inits=NULL,
           n_folds=5, measure.inc=c("c","auc"), one.se=FALSE, cure_cutoff=5,
           parallel=FALSE, seed=NULL, verbose=TRUE){
    if(!is.null(seed)) set.seed(seed)
    X_k = knockoff::create.second_order(X_p, method = "asdp", shrink = T)
    W_k = knockoff::create.second_order(W_p, method = "asdp", shrink = T)
    Xaug = cbind(X_p, X_k)
    Waug = cbind(W_p, W_k)
    fit = cv.em.nofdr(X_u, Xaug, W_u, Waug, time, delta, model, penalty, thresh, nIter,
                      grid.tuning, lambda.inc.list, lambda.lat.list, nlambda.inc, nlambda.lat,
                      lambda.min.ratio.inc, lambda.min.ratio.lat, gamma.inc, gamma.lat, inits,
                      n_folds, measure.inc, one.se, cure_cutoff, parallel, seed, verbose)
    if(is.null(X_u)) b_p = fit$b else {
      b_u = fit$b[1:ncol(X_u)]
      b_p = fit$b[(ncol(X_u)+1):(length(fit$b))]
    }
    if(is.null(W_u)) beta_p = fit$beta else {
      beta_u = fit$beta[1:ncol(W_u)]
      beta_p = fit$beta[(ncol(W_u)+1):(length(fit$beta))]
    }
    Z_b = abs(b_p)
    Z_beta = abs(beta_p)
    J = ncol(X_p)
    M = ncol(W_p)
    orig_b = 1:J
    orig_beta = 1:M
    W_b = Z_b[orig_b] - Z_b[orig_b+J]
    W_beta = Z_beta[orig_beta] - Z_beta[orig_beta+M]
    T_b = knockoff::knockoff.threshold(W_b, fdr = fdr, offset = 1)
    T_beta = knockoff::knockoff.threshold(W_beta, fdr = fdr, offset = 1)
    selected_b = (1:J)[W_b > T_b]
    selected_beta = (1:M)[W_beta > T_beta]
    # KJA 05-09 added b_u code above and this code to output b0, b, beta, rate, and alpha
    if (identical(unname(selected_b), integer(0))) {
      b_p[1:J] <- 0
    } else {
      b_p[-selected_b] <- 0
    }
    if (identical(unname(selected_beta), integer(0))) {
      beta_p[1:M] <- 0
    } else {
      beta_p[-selected_beta] <- 0
    }
    if (is.null(X_u)) {
      b = b_p[1:J]
    } else {
       b = c(b_u, b_p[1:J])
    }
    if (is.null(W_u)) {
      beta = beta_p[1:M]
    } else {
       beta = c(beta_u, beta_p[1:M])
    }
    return(list(b0 = fit$b0, b = b, beta = beta, rate = fit$rate, alpha = fit$alpha, selected_b=selected_b, selected_beta=selected_beta))
  }

cv.em.inner <-
  function(X_u, X_p, W_u, W_p, time, delta, folds_i, k,
           model=c("cox","weibull","exponential"), penalty=c("lasso","MCP","SCAD"),
           thresh=1e-5, nIter=1e4,
           grid.tuning = FALSE, tuning_sequence = NULL,
           gamma.inc=3, gamma.lat=3, inits=NULL,
           measure.inc=c("c","auc"), cure_cutoff=5){
    test_i = which(folds_i == k)
    test_delta = delta[test_i]
    test_time = time[test_i]
    test_X_p = X_p[test_i,,drop=FALSE]
    test_W_p = W_p[test_i,,drop=FALSE]
    if(is.null(X_u)) test_X_u = NULL else test_X_u = X_u[test_i,,drop=FALSE]
    if(is.null(W_u)) test_W_u = NULL else test_W_u = W_u[test_i,,drop=FALSE]
    if(is.null(X_u)) X_u_train = NULL else X_u_train = X_u[-test_i,,drop=FALSE]
    if(is.null(W_u)) W_u_train = NULL else W_u_train = W_u[-test_i,,drop=FALSE]
    if(is.null(X_u)) penalty.factor.inc = rep(1,ncol(X_p))
    else penalty.factor.inc = rep(0:1,c(ncol(X_u), ncol(X_p)))
    if(is.null(W_u)) penalty.factor.lat = rep(1,ncol(W_p))
    else penalty.factor.lat = rep(0:1,c(ncol(W_u), ncol(W_p)))
    initial_val = inits
    if(model=="cox") initial_val$survprob = initial_val$survprob[-test_i]
    cst = rep(NA, nrow(tuning_sequence))
    if(measure.inc=="auc") auc = rep(NA, nrow(tuning_sequence))
    for (i in 1:nrow(tuning_sequence)){
      lambda.inc = tuning_sequence[i,1]
      if(grid.tuning) lambda.lat = tuning_sequence[i,2] else lambda.lat = lambda.inc
      train_out = cure.em(X_u_train, X_p[-test_i,], W_u_train, W_p[-test_i,],
                          time[-test_i], delta[-test_i], model, penalty,
                          penalty.factor.inc, penalty.factor.lat, thresh, nIter,
                          lambda.inc, lambda.lat, gamma.inc, gamma.lat, initial_val)

      #train_out = cure.em(cbind(X_u_train, X_p[-test_i,]), cbind(W_u_train, W_p[-test_i,]),
      #                    time[-test_i], delta[-test_i], model, penalty,
      #                    penalty.factor.inc, penalty.factor.lat, thresh, nIter,
      #                    lambda.inc, lambda.lat, gamma.inc, gamma.lat, initial_val)
      if(is.null(X_u)) {b_u = NULL; b_p = train_out$b}
      else {b_u = train_out$b[1:ncol(X_u)]; b_p = train_out$b[(ncol(X_u)+1):(length(train_out$b))]}
      if(is.null(W_u)) {beta_u = NULL; beta_p = train_out$beta}
      else {beta_u = train_out$beta[1:ncol(W_u)]; beta_p = train_out$beta[(ncol(W_u)+1):(length(train_out$beta))]}
      cst[i] = C.stat(cure_cutoff, b_u, train_out$b0, b_p, beta_u, beta_p,
                      test_delta, test_time, test_X_u, test_X_p, test_W_u, test_W_p)
      if(measure.inc=="auc") auc[i] = AUC_msi(cure_cutoff = cure_cutoff, b_u, train_out$b0, b_p,
                                              test_delta, test_time, test_X_u, test_X_p)
    }
    if(measure.inc=="auc"){
      res = rbind(auc, cst)
      rownames(res) = c("AUC","Cstat")
    }
    else res = cst
    return(res)
  }

cv.em.nofdr <-
  function(X_u, X_p, W_u, W_p, time, delta, model=c("cox","weibull","exponential"),
           penalty=c("lasso","MCP","SCAD"),
           thresh=1e-03, nIter=ifelse(penalty=="lasso", 100, 1000),
           grid.tuning = FALSE, lambda.inc.list=NULL, lambda.lat.list=NULL,
           nlambda.inc = ifelse(grid.tuning, 10, 50),
           nlambda.lat = ifelse(grid.tuning, 10, 50),
           lambda.min.ratio.inc = 0.1, lambda.min.ratio.lat = 0.1,
           gamma.inc=3, gamma.lat=3,
           inits=NULL,
           n_folds=5, measure.inc=c("c","auc"), one.se=FALSE, cure_cutoff=5,
           parallel=FALSE, seed=NULL, verbose=TRUE){
    if(!is.null(seed)) set.seed(seed)
    if(!grid.tuning){
      if((!is.null(lambda.inc.list) & !is.null(lambda.lat.list) &
          !identical(lambda.inc.list, lambda.lat.list)) | (nlambda.inc!=nlambda.lat) |
         (lambda.min.ratio.inc!=lambda.min.ratio.lat))
        warning("Grid tuning is off. Same lambda sequence for incidence and latency was used.")
      if(length(lambda.inc.list)>length(lambda.lat.list)) lambda.lat.list=lambda.inc.list
      else lambda.inc.list=lambda.lat.list
      nlambda.inc <- nlambda.lat <- max(nlambda.inc, nlambda.lat)
      lambda.min.ratio.inc <- lambda.min.ratio.lat <- min(lambda.min.ratio.inc, lambda.min.ratio.lat)
    }

    # decide lambda grid or sequence for tuning
    N = length(time)
    if(is.null(inits)){
      if(model=="cox") inits = initialization_cox(X_u, W_u, time, delta)
      else inits = initialization_parametric(X_u, W_u, time, delta, model)
    }
    if(grid.tuning){
      if(is.null(lambda.inc.list) | is.null(lambda.lat.list)){
        pir = rep(1, N)
        if(is.null(X_u)) pi_x = rep(mean(delta), N)
        else pi_x = 1/(1+exp(-inits$itct-X_u %*% inits$b_u))
        if(is.null(W_u)) beta_w = rep(0,N) else beta_w = W_u %*% inits$beta_u
        if(model=="cox") survprob = inits$survprob
        else if(model=="weibull")
          survprob = exp(-(inits$lambda * time)^(inits$alpha) * exp(beta_w))
        else if(model=="exponential") survprob = exp(-inits$lambda * time * exp(beta_w))
        numerator = pi_x[delta==0] * survprob[delta==0]
        pir[delta==0] = numerator/(1-pi_x[delta==0] + numerator)
        if(is.null(lambda.inc.list)){
          lambda.max.inc = max(abs(matrix(pir-0.5,1,N)%*%X_p))/N
          lambda.min.inc = lambda.max.inc * lambda.min.ratio.inc
          k1 = (0:(nlambda.inc-1)) / nlambda.inc
          lambda.inc.list = lambda.max.inc * lambda.min.ratio.inc^k1
        }
        if(is.null(lambda.lat.list)){
          nz_pir = which(pir>0)
          lambda.max.lat = get_cox_lambda_max(x = cbind(W_u, W_p)[nz_pir,],
                                              time = time[nz_pir], event = delta[nz_pir],
                                              offset = log(pir[nz_pir]))
          lambda.min.lat = lambda.max.lat * lambda.min.ratio.lat
          k2 = (0:(nlambda.lat-1)) / nlambda.lat
          lambda.lat.list = lambda.max.lat * lambda.min.ratio.lat^k2
        }
      }
      lambda.grid = expand.grid(lambda.inc.list, lambda.lat.list)
    }else{ # use same sequence to tune lambda.inc and lambda.lat
      if(is.null(lambda.inc.list)){
        lambda_max = max(colSums(abs(X_p)), colSums(abs(W_p)))/(2*N)
        lambda_min = lambda_max * lambda.min.ratio.inc
        k1 = (0:(nlambda.inc-1)) / nlambda.inc
        lambda.list = lambda_max * lambda.min.ratio.inc^k1
      }else lambda.list = lambda.inc.list
    }
    if(grid.tuning) tuning_sequence = lambda.grid
    else tuning_sequence = matrix(lambda.list, ncol = 1)

    # cross-validation
    folds_i = sample(rep(1:n_folds, length.out = N))
    Cstat = matrix(NA, nrow(tuning_sequence), n_folds)
    if(measure.inc=="auc") AUC = matrix(NA, nrow(tuning_sequence), n_folds)
    if(parallel){ # parallel computing
      ncores = n_folds
      registerDoMC(ncores)
      res_list = foreach(k = 1:n_folds) %dopar% {
        res = cv.em.inner(X_u, X_p, W_u, W_p, time, delta, folds_i, k,
                          model,penalty,thresh, nIter, grid.tuning, tuning_sequence,
                          gamma.inc, gamma.lat, inits,
                          measure.inc, cure_cutoff)
        return(res)
      }
      if(measure.inc=="auc"){
        for (k in 1:n_folds){
          AUC[,k] = res_list[[k]][1,]
          Cstat[,k] = res_list[[k]][2,]
        }
      }else{
        for (k in 1:n_folds)  Cstat[,k] = res_list[[k]]
      }
    }else{ # no parallel computing
      for (k in 1:n_folds) {
        if(verbose) cat("Fold", k, "out of", n_folds, "training...\n")
        res = cv.em.inner(X_u, X_p, W_u, W_p, time, delta, folds_i, k,
                          model,penalty,thresh, nIter, grid.tuning, tuning_sequence,
                          gamma.inc, gamma.lat, inits,
                          measure.inc, cure_cutoff)
        if(measure.inc=="auc"){
          AUC[,k] = res[1,]
          Cstat[,k] = res[2,]
        }else Cstat[,k] = res
      }
    }

    if(measure.inc=="auc"){ # use AUC to tune incidence and C-statistic to tune latency
      auc = rowMeans(AUC, na.rm = T)
      c_stat = rowMeans(Cstat, na.rm = T)
      if(one.se){
        auc_sd = sqrt(apply(AUC, 1, var, na.rm = T)*(n_folds-1)/n_folds/(nrow(X_p)-1))
        c_stat_sd = sqrt(apply(Cstat, 1, var, na.rm = T)*(n_folds-1)/n_folds/(nrow(X_p)-1))
        opt_step_b = which.max(auc)
        opt_step_beta = which.max(c_stat)
        model_select_b = which(auc >= (auc[opt_step_b] - auc_sd[opt_step_b]))[1]
        model_select_beta = which(c_stat >= (c_stat[opt_step_beta] - c_stat_sd[opt_step_beta]))[1]
      }else{
        model_select_b = which.max(auc)
        model_select_beta = which.max(c_stat)
      }
    }
    else{ # use C-statistic to tune both components
      c_stat = rowMeans(Cstat, na.rm = T)
      if(one.se){
        c_stat_sd = sqrt(apply(Cstat, 1, var, na.rm = T)*(n_folds-1)/n_folds/(nrow(X_p)-1))
        opt_step = which.max(c_stat)
        model_select <- model_select_b <- model_select_beta <-
          which(c_stat >= (c_stat[opt_step] - c_stat_sd[opt_step]))[1]
      }else model_select <- model_select_b <- model_select_beta <- which.max(c_stat)
    }
    optimal_lambda_b = tuning_sequence[model_select_b,1]
    if(grid.tuning) optimal_lambda_beta = tuning_sequence[model_select_beta,2]
    else optimal_lambda_beta = tuning_sequence[model_select_beta,1]
    if(verbose){
      cat("Selected lambda for incidence:",round(optimal_lambda_b,3),"\n")
      cat("Selected lambda for latency:", round(optimal_lambda_beta,3),"\n")
      cat("Maximum C-statistic:",max(c_stat),"\n")
    }
    if(measure.inc=="auc") cat("Maximum AUC:",max(auc),"\n")
    if(verbose) cat("Fitting a final model...\n")
    if(is.null(X_u)) penalty.factor.inc = rep(1,ncol(X_p))
    else penalty.factor.inc = rep(0:1,c(ncol(X_u), ncol(X_p)))
    if(is.null(W_u)) penalty.factor.lat = rep(1,ncol(W_p))
    else penalty.factor.lat = rep(0:1,c(ncol(W_u), ncol(W_p)))
    output = cure.em(X_u, X_p, W_u, W_p,
                     time, delta, model, penalty,
                     penalty.factor.inc, penalty.factor.lat, thresh, nIter,
                     optimal_lambda_b, optimal_lambda_beta, gamma.inc, gamma.lat, inits)
    #output = cure.em(cbind(X_u, X_p), cbind(W_u, W_p),
    #                 time, delta, model, penalty,
    #                 penalty.factor.inc, penalty.factor.lat, thresh, nIter,
    #                 optimal_lambda_b, optimal_lambda_beta, gamma.inc, gamma.lat, inits)
    output$selected.lambda.inc=optimal_lambda_b
    output$selected.lambda.lat=optimal_lambda_beta
    output$max.c = max(c_stat)
    if(measure.inc=="auc") output$max.auc = max(auc)
    return(output)
  }

cv.gmifs.fdr <-
  function(X_u, X_p, W_u, W_p, time, delta, model=c("weibull","exponential"),
           fdr=0.2, thresh=1e-5, nIter=1e4, epsilon=0.001, inits=NULL,
           n_folds=5, measure.inc=c("c","auc"), one.se=FALSE, cure_cutoff=5,
           parallel=FALSE, seed=NULL, verbose=TRUE){
    if(!is.null(seed)) set.seed(seed)
    X_k = knockoff::create.second_order(X_p, method = "asdp", shrink = T)
    W_k = knockoff::create.second_order(W_p, method = "asdp", shrink = T)
    Xaug = cbind(X_p, X_k)
    Waug = cbind(W_p, W_k)
    fit = cv.gmifs.nofdr(X_u, Xaug, W_u, Waug, time, delta, model, thresh, nIter, epsilon, inits,
                         n_folds, measure.inc, one.se, cure_cutoff,
                         parallel, seed, verbose)
    b_p <- fit$b_p
    b_u <- fit$b_u
    beta_p <- fit$beta_p
    beta_u <- fit$beta_u
    Z_b = abs(fit$b_p)
    Z_beta = abs(fit$beta_p)
    J = ncol(X_p)
    M = ncol(W_p)
    orig_b = 1:J
    orig_beta = 1:M
    W_b = Z_b[orig_b] - Z_b[orig_b+J]
    W_beta = Z_beta[orig_beta] - Z_beta[orig_beta+M]
    T_b = knockoff::knockoff.threshold(W_b, fdr = fdr, offset = 1)
    T_beta = knockoff::knockoff.threshold(W_beta, fdr = fdr, offset = 1)
    selected_b = (1:J)[W_b > T_b]
    selected_beta = (1:M)[W_beta > T_beta]
    # KJA 05-09 added b_u code above and this code to output b0, b, beta, rate, and alpha
    b_p <- b_p[1:J]
    beta_p <- beta_p[1:M]
    if (identical(unname(selected_b), integer(0))) {
      b_p[1:J] <- 0
    } else {
      b_p[-selected_b] <- 0
    }
    if (identical(unname(selected_beta), integer(0))) {
      beta_p[1:M] <- 0
    } else {
      beta_p[-selected_beta] <- 0
    }
    return(list(b0 = fit$b0, b_u = b_u, b_p = b_p, beta_u = beta_u, beta_p = beta_p, rate = fit$rate, alpha = fit$alpha, selected_b=selected_b, selected_beta=selected_beta))
  }

cv.gmifs.inner <-
  function(X_u, X_p, W_u, W_p, time, delta, folds_i, k,
           model=c("weibull","exponential"),thresh=1e-5, nIter=1e4, epsilon=0.001,
           inits=NULL, measure.inc=c("c","auc"), cure_cutoff=5, verbose=TRUE){
    test_i = which(folds_i == k)
    test_delta = delta[test_i]
    test_time = time[test_i]
    test_X_p = X_p[test_i,,drop=F]
    test_W_p = W_p[test_i,,drop=F]
    if(is.null(X_u)) test_X_u = NULL else test_X_u = X_u[test_i,,drop=F]
    if(is.null(W_u)) test_W_u = NULL else test_W_u = W_u[test_i,,drop=F]
    if(is.null(X_u)) X_u_train = NULL else X_u_train = X_u[-test_i,,drop=F]
    if(is.null(W_u)) W_u_train = NULL else W_u_train = W_u[-test_i,,drop=F]
    if(model=="exponential")
      train_out = exp_cure(X_u_train, X_p[-test_i,],
                           W_u_train, W_p[-test_i,],
                           time[-test_i], delta[-test_i],
                           epsilon, thresh, nIter, inits, verbose)
    else if(model=="weibull")
      train_out = weibull.cure(X_u_train, X_p[-test_i,],
                               W_u_train, W_p[-test_i,],
                               time[-test_i], delta[-test_i],
                               epsilon, thresh, nIter, inits, verbose)
    cst = sapply(1:length(train_out$itct_path),
                 function(x) {
                   if(is.null(X_u)) b_u = NULL else b_u = train_out$b_u_path[x,]
                   if(is.null(W_u)) beta_u = NULL else beta_u = train_out$beta_u_path[x,]
                   cs = C.stat(cure_cutoff=cure_cutoff,
                               b_u, train_out$itct_path[x], train_out$b_p_path[x,],
                               beta_u, train_out$beta_p_path[x,],
                               test_delta, test_time, test_X_u, test_X_p, test_W_u, test_W_p)
                   return(cs)
                 })
    if(measure.inc=="auc"){
      auc = sapply(1:length(train_out$itct_path),
                   function(x) {
                     if(is.null(X_u)) b_u = NULL else b_u = train_out$b_u_path[x,]
                     if(is.null(W_u)) beta_u = NULL else beta_u = train_out$beta_u_path[x,]
                     a = AUC_msi(cure_cutoff = cure_cutoff, b_u, train_out$itct_path[x],
                                 train_out$b_p_path[x,], test_delta, test_time, test_X_u, test_X_p)
                     return(a)
                   })
      res = rbind(c(auc, rep(NA, nIter-length(auc))), c(cst, rep(NA, nIter-length(cst))))
      rownames(res) = c("AUC","Cstat")
    }
    else res = c(cst, rep(NA, nIter-length(cst)))
    return(res)
  }

cv.gmifs.nofdr <-
  function(X_u, X_p, W_u, W_p, time, delta, model=c("weibull","exponential"),
           thresh=1e-5, nIter=1e4, epsilon=0.001, inits=NULL,
           n_folds=5, measure.inc=c("c","auc"), one.se=FALSE, cure_cutoff=5,
           parallel=FALSE, seed=NULL, verbose=TRUE){
    if(!is.null(seed)) set.seed(seed)
    folds_i = sample(rep(1:n_folds, length.out = length(time)))
    Cstat = matrix(NA, nIter, n_folds)
    if(measure.inc=="auc") AUC = matrix(NA, nIter, n_folds)
    if(parallel){ # parallel computing
      ncores = n_folds
      registerDoMC(ncores)
      res_list = foreach(k = 1:n_folds) %dopar% {
        res = cv.gmifs.inner(X_u, X_p, W_u, W_p, time, delta, folds_i, k,
                             model,thresh, nIter, epsilon, inits, measure.inc, cure_cutoff, verbose)
        return(res)
      }
      if(measure.inc=="auc"){
        for (k in 1:n_folds){
          AUC[,k] = res_list[[k]][1,]
          Cstat[,k] = res_list[[k]][2,]
        }
      }else{
        for (k in 1:n_folds)  Cstat[,k] = res_list[[k]]
      }
    }else{ # no parallel computing
      for (k in 1:n_folds) {
        if(verbose) cat("Fold", k, "out of", n_folds, "training...\n")
        res = cv.gmifs.inner(X_u, X_p, W_u, W_p, time, delta, folds_i, k,
                             model,thresh, nIter, epsilon, inits, measure.inc, cure_cutoff, verbose)
        if(measure.inc=="auc"){
          AUC[,k] = res[1,]
          Cstat[,k] = res[2,]
        }else Cstat[,k] = res
      }
    }
    if(measure.inc=="auc"){ # use AUC to tune incidence and C-statistic to tune latency
      auc = rowMeans(AUC, na.rm = T)
      c_stat = rowMeans(Cstat, na.rm = T)
      if(one.se){
        auc_sd = sqrt(apply(AUC, 1, var, na.rm = T)*(n_folds-1)/n_folds/(nrow(X_p)-1))
        c_stat_sd = sqrt(apply(Cstat, 1, var, na.rm = T)*(n_folds-1)/n_folds/(nrow(X_p)-1))
        opt_step_b = which.max(auc)
        opt_step_beta = which.max(c_stat)
        model_select_b = which(auc >= (auc[opt_step_b] - auc_sd[opt_step_b]))[1]
        model_select_beta = which(c_stat >= (c_stat[opt_step_beta] - c_stat_sd[opt_step_beta]))[1]
      }else{
        model_select_b = which.max(auc)
        model_select_beta = which.max(c_stat)
      }
      if(verbose){
        cat("Selected step for incidence:",model_select_b,"\n")
        cat("Selected step for latency:",model_select_beta,"\n")
        cat("Maximum AUC:",max(auc, na.rm = T),"\n")
        cat("Maximum C-statistic:",max(c_stat, na.rm = T),"\n")
      }
      if(verbose) cat("Fitting a final model...\n")
      if(model=="exponential")
        fit = exp_cure(X_u, X_p, W_u, W_p, time, delta, epsilon, thresh,
                       nIter=max(model_select_b, model_select_beta), inits, verbose)
      else if(model=="weibull")
        fit = weibull.cure(X_u, X_p, W_u, W_p, time, delta, epsilon, thresh,
                           nIter=max(model_select_b, model_select_beta), inits, verbose)
    }
    else{ # use C-statistic to tune both components
      c_stat = rowMeans(Cstat, na.rm = T)
      if(one.se){
        c_stat_sd = sqrt(apply(Cstat, 1, var, na.rm = T)*(n_folds-1)/n_folds/(nrow(X_p)-1))
        opt_step = which.max(c_stat)
        model_select <- model_select_b <- model_select_beta <-
          which(c_stat >= (c_stat[opt_step] - c_stat_sd[opt_step]))[1]
      }else model_select <- model_select_b <- model_select_beta <- which.max(c_stat)
      if(verbose){
        cat("Selected step:",model_select,"\n")
        cat("Maximum C-statistic:",max(c_stat, na.rm = T),"\n")
        cat("Fitting a final model...\n")
      }
      if(model=="exponential")
        fit = exp_cure(X_u, X_p, W_u, W_p, time, delta, epsilon, thresh, nIter=model_select, inits,
                       verbose)
      else if(model=="weibull")
        fit = weibull.cure(X_u, X_p, W_u, W_p, time, delta, epsilon, thresh, nIter=model_select,
                           inits, verbose)
    }
    nstep = length(fit$logLikelihood)
    b_p = fit$b_p_path[nstep,]
    b_u = fit$b_u_path[nstep,]
    b0 = fit$itct_path[nstep]
    beta_p = fit$beta_p_path[nstep,]
    beta_u = fit$beta_u_path[nstep,]
    rate = fit$lambda_path[nstep]
    if(model=="weibull") alpha = fit$alpha_path[nstep]
    else alpha = 1
    output = list(model.select.inc=model_select_b, model.select.lat=model_select_beta,
                  b0=b0, b_u=b_u, b_p=b_p, beta_u=beta_u, beta_p=beta_p, alpha=alpha, rate=rate,
                  logLik = fit$logLikelihood[nstep], max.c = max(c_stat, na.rm = T))
    if(measure.inc=="auc") output$max.auc = max(auc, na.rm = T)
    return(output)
  }

exp_cure <-
  function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, epsilon = 0.001, tol = 1e-05,
           nIter=1e4, inits=NULL, verbose=TRUE)
  {

    # X_u: N by nonp, non-penalized covariate matrix associated with incidence
    # X_p: N by J, penalized covariate matrix associated with incidence
    # W_u: N by nonp, non-penalized covariate matrix associated with latency
    # W_p:  N by M, penalized covariate matrix associated with incidence
    # time: vector of length N, observed survival time
    # delta: vector of length N, censoring status (not censored = 0)
    # epsilon: incremental size
    # tol: difference between log-likelihood

    N = length(time)
    J = ncol(X_p) # number of penalized incidence covariates
    M = ncol(W_p) # number of penalized latency covariates
    X_p = cbind(X_p, -X_p) # N by 2J
    W_p = cbind(W_p, -W_p) # N by 2M

    ######## initialization ########
    step = 1
    b_p = rep(0, 2*J) # penalized incidence
    beta_p = rep(0, 2*M) # penalized latency
    if(is.null(inits)) inits = initialization_parametric(X_u, W_u, time, delta, model="weibull")
    itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u; lambda = inits$lambda
    log_lambda = log(max(lambda, 1e-15))
    LL0 = 0

    b_p_path <- beta_p_path <- lambda_path <- b_u_path <- beta_u_path <- itct_path <- NULL
    logLikelihood <- numeric()

    ######## loop ########

    repeat{

      #### update penalized parameters
      upd = exp_update(lambda, b_p, beta_p, b_u, itct, beta_u, X_u , X_p, W_u, W_p, time, delta, epsilon)
      b_p = upd$b_p
      beta_p = upd$beta_p
      b_p_path = rbind(b_p_path, b_p)
      beta_p_path = rbind(beta_p_path, beta_p)

      #### update other parameters
      out = optim(par = c(log_lambda, b_u, itct, beta_u), fn = exp_negloglik, gr = exp_gradient,
                  b_p = b_p, beta_p = beta_p, X_u = X_u, X_p = X_p, W_u = W_u, W_p = W_p,
                  time = time, delta = delta, method="BFGS")
      log_lambda = out$par[1]
      lambda = exp(log_lambda)
      if(!is.null(X_u)) {
        b_u = out$par[2:(ncol(X_u)+1)] # unpenalized incidence
        itct = out$par[ncol(X_u)+2]
        if(!is.null(W_u)) beta_u = out$par[(ncol(X_u)+3):(ncol(X_u)+ncol(W_u)+2)] # unpenalized latency
      } else{
        itct = out$par[2]
        if(!is.null(W_u)) beta_u = out$par[3:(ncol(W_u)+2)] # unpenalized latency
      }

      lambda_path = c(lambda_path, lambda)
      if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
      itct_path = c(itct_path, itct)
      if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)

      LL1 = -out$value
      logLikelihood = c(logLikelihood, LL1)

      if(verbose & step%%1000==0) cat("step=", step, "\n")
      if (step > 1 && (abs(LL1 - LL0) < tol | step >= nIter)) { #
        break
      }
      LL0 <- LL1
      step <- 1 + step
    }

    ######## output ########
    b_p_path = b_p_path[,1:J] - b_p_path[,(J+1):(2*J)]
    beta_p_path = beta_p_path[,1:M] - beta_p_path[,(M+1):(2*M)]
    output <- list(b_p_path = b_p_path, beta_p_path = beta_p_path,
                   lambda_path = lambda_path, b_u_path = b_u_path, itct_path = itct_path,
                   beta_u_path = beta_u_path, logLikelihood = logLikelihood)
    output
  }


exp_EM <-
  function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, mu_inc, mu_lat,
           inits=NULL, nIter = 100, tol = 1e-4){
    # mu: penalty parameter
    # tol: difference between log-likelihood
    N = length(time)
    J = ncol(X_p) # number of penalized incidence covariates
    M = ncol(W_p) # number of penalized latency covariates

    ######## initialization ########
    step = 1
    b_p <- rep(0,J)
    beta_p <- rep(0,M)
    b_p_ext <- rep(0,2*J)
    beta_p_ext <- rep(0,2*M)
    if(is.null(inits)) inits = initialization_parametric(X_u, W_u, time, delta, model="weibull")
    itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u; lambda = inits$lambda
    log_lambda = log(max(lambda, 1e-15))
    pir = rep(1, N)
    llp1_0 <- llp2_0 <-0
    lik_inc <- lik_lat <- c()
    b_p_path <- beta_p_path <- b_u_path <- beta_u_path <- itct_path <- lambda_path <- NULL
    conv1 <- conv2 <- FALSE

    ######## loop ########

    repeat{

      #### E-step
      if(is.null(X_u)) b_x = itct + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
      else b_x = itct + X_u %*% b_u + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
      if(is.null(W_u)) beta_w = W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
      else beta_w = W_u %*% beta_u + W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
      pir[delta==0] = 1/(1+exp(-b_x[delta==0]+lambda * time[delta==0] * exp(beta_w[delta==0])))
      #### update penalized parameters
      if (!conv1){
        b_p_ext = optim(par = b_p_ext, fn = l1_negloglik_inc_b, gr = l1_grad_b,
                        itct = itct, b_u = b_u, X_u = X_u, X_p = X_p,
                        pir = pir, mu = mu_inc, method = 'L-BFGS-B', lower = rep(0, 2*J))$par
        b_p = b_p_ext[1:J]-b_p_ext[(J+1):(2*J)]
      }

      if (!conv2){
        beta_p_ext = optim(par = beta_p_ext, fn = exp_negloglik_lat_beta, gr = exp_grad_beta,
                           lambda = lambda,
                           beta_u = beta_u, W_u = W_u, W_p = W_p, time = time, delta = delta, pir = pir, mu = mu_lat,
                           method = 'L-BFGS-B', lower = rep(0, 2*M))$par
        beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
      }

      #### update nonpenalized parameters
      if (!conv1){
        out_inc = optim(par = c(itct,b_u), fn = l1_negloglik_inc, gr = l1_gradient_inc,
                        b_p = b_p, X_u =X_u, X_p = X_p, pir = pir, mu = mu_inc, method="BFGS") # incidence
        itct = out_inc$par[1]
        if(!is.null(X_u)) b_u = out_inc$par[2:(ncol(X_u)+1)]
        llp1_1 = -out_inc$value
      }
      lik_inc = c(lik_inc, llp1_1)

      if (!conv2){
        out_lat = optim(par = c(log_lambda, beta_u), fn = exp_negloglik_lat, gr = exp_gradient_lat,
                        beta_p=beta_p, W_u=W_u, W_p=W_p, time=time, delta=delta, pir=pir, mu=mu_lat,
                        method="BFGS") # latency
        log_lambda = out_lat$par[1]
        lambda = exp(log_lambda)
        if(!is.null(W_u)) beta_u = out_lat$par[2:(ncol(W_u)+1)] # unpenalized incidence
        llp2_1 = -out_lat$value
      }
      lik_lat = c(lik_lat, llp2_1)

      ## record updated parameters
      itct_path = c(itct_path, itct)
      if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
      b_p_path = rbind(b_p_path, b_p)
      if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)
      beta_p_path = rbind(beta_p_path, beta_p)
      lambda_path = c(lambda_path, lambda)
      if (!conv1 & abs(llp1_1- llp1_0)< tol) conv1 <- TRUE
      if (!conv2 & abs(llp2_1- llp2_0)< tol) conv2 <- TRUE
      if (step > 1 & (conv1 & conv2) | step >= nIter) { #
        break
      }
      llp1_0 <- llp1_1
      llp2_0 <- llp2_1
      step <- 1 + step
    }

    ######## output ########
    output <- list(b_p_path = b_p_path, b_u_path = b_u_path, itct_path = itct_path,
                   beta_p_path = beta_p_path, beta_u_path = beta_u_path, lambda_path = lambda_path,
                   lik_inc = lik_inc, lik_lat = lik_lat)
    output
  }

exp_grad_beta <-
  function(theta, lambda, beta_u, W_u, W_p, time, delta, pir, mu){
    N = nrow(W_p)
    M = ncol(W_p)
    beta_p_ext = theta
    beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    temp1 = pir*lambda*time*exp(betaw)
    grad = matrix(delta-temp1,1)%*% W_p
    return(c(-grad+ N*mu, grad+ N*mu))
  }

exp_gradient_lat <-
  function(theta, beta_p, W_u, W_p, time, delta, pir, mu){ # latency
    log_lambda = theta[1]
    lambda = exp(log_lambda)
    if(!is.null(W_u)) beta_u = theta[2:(ncol(W_u)+1)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    temp1 = pir*lambda*time*exp(betaw)
    temp2 = delta-temp1
    grad2 = sum(temp2)
    if(!is.null(W_u)) grad3 = matrix(temp2,1) %*% W_u else grad3=NULL
    return(-c(grad2, grad3))
  }

exp_gradient <-
  function(theta, b_p, beta_p, X_u , X_p, W_u, W_p, time, delta)
  {
    N = length(time)
    lambda = exp(theta[1])
    if(!is.null(X_u)) {
      b_u = theta[2:(ncol(X_u)+1)] # unpenalized incidence
      itct = theta[ncol(X_u)+2]
      if(!is.null(W_u)) beta_u = theta[(ncol(X_u)+3):(ncol(X_u)+ncol(W_u)+2)] # unpenalized latency
    } else{
      itct = theta[2]
      if(!is.null(W_u)) beta_u = theta[3:(ncol(W_u)+2)] # unpenalized latency
    }
    b_nonzero = which(b_p!=0)
    beta_nonzero = which(beta_p!=0)
    if(!is.null(X_u)) C_b = exp(itct + X_u %*% b_u + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    else C_b = exp(itct + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    if(!is.null(W_u)) C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    else C_beta = exp(W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    C_lbeta = exp(-lambda * time * C_beta)
    temp1 = log (pmax(lambda*time,rep(1e-100,N)))
    temp2 = lambda * time * C_beta
    temp3 = 1 / (1 + C_b*C_lbeta)
    temp4 = (1-delta) * C_b * C_lbeta * temp2 * temp3
    temp5 = delta * (1-temp2)
    temp6 = 1/(1+C_b)
    temp7 = temp6*(delta - (1-delta) * C_b * (1-C_lbeta) * temp3)
    grad2 = sum(temp5 - temp4)
    if(!is.null(X_u)) grad3 = matrix(temp7,1) %*% X_u else grad3=NULL
    grad4 = sum(temp7) # intercept
    if(!is.null(W_u)) grad5 = matrix(temp5 - temp4, 1) %*% W_u else grad5=NULL
    return(-c(grad2, grad3, grad4, grad5))
  }

exp_negloglik_lat_beta <-
  function(theta, lambda, beta_u, W_u, W_p, time, delta, pir, mu){ # latency
    N = nrow(W_p)
    M = ncol(W_p)
    beta_p_ext = theta
    beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    ll1 = delta*( log(lambda)+betaw )
    ll2 = -pir*lambda*time*exp(betaw)
    return(-sum(ll1+ll2)+N*mu*sum(beta_p_ext))
  }

exp_negloglik_lat <-
  function(theta, beta_p, W_u, W_p, time, delta, pir, mu){ # latency
    N = nrow(W_p)
    log_lambda = theta[1]
    lambda = exp(log_lambda)
    if(!is.null(W_u)) beta_u = theta[2:(ncol(W_u)+1)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    ll1 = delta*( log_lambda +betaw )
    ll2 = -pir*lambda*time*exp(betaw)
    return(-sum(ll1+ll2)+N*mu*sum(abs(beta_p)))
  }

exp_negloglik <-
  function(theta, b_p, beta_p, X_u , X_p, W_u, W_p, time, delta)
  {
    N = length(time)
    lambda = exp(theta[1])
    if(!is.null(X_u)) {
      b_u = theta[2:(ncol(X_u)+1)] # unpenalized incidence
      itct = theta[ncol(X_u)+2]
      if(!is.null(W_u)) beta_u = theta[(ncol(X_u)+3):(ncol(X_u)+ncol(W_u)+2)] # unpenalized latency
    } else{
      itct = theta[2]
      if(!is.null(W_u)) beta_u = theta[3:(ncol(W_u)+2)] # unpenalized latency
    }
    b_nonzero = which(b_p!=0)
    beta_nonzero = which(beta_p!=0)
    if(!is.null(X_u)) logC_b = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else logC_b = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    if(!is.null(W_u)) logC_beta = W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero]
    else logC_beta = W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero]
    temp1 = 1/(1+exp(-logC_b)) # exp(logC_b) / (1+exp(logC_b)) #Mar17
    logC_lbeta = - lambda * time * exp(logC_beta)
    ll1 = delta * (log(temp1) + theta[1] + logC_beta + logC_lbeta)
    ll2 = (1-delta) * log(pmax(1 - temp1 * (1-exp(logC_lbeta)),rep(1e-100,N)))
    return(-sum(ll1+ll2))
  }

exp_update <-
  function(lambda, b_p, beta_p, b_u, itct, beta_u, X_u , X_p, W_u, W_p, time, delta, epsilon)
  {
    b_nonzero = which(b_p!=0)
    beta_nonzero = which(beta_p!=0)
    if(!is.null(X_u)) C_b = exp(itct + X_u %*% b_u + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    else C_b = exp(itct + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    if(!is.null(W_u)) C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    else C_beta = exp(W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    C_lbeta = exp(-lambda * time * C_beta)
    temp2 = lambda * time * C_beta
    temp3 = 1 / (1 + C_b*C_lbeta)
    temp4 = (1-delta) * C_b * C_lbeta * temp2 * temp3
    temp5 = delta * (1-temp2)
    temp6 = 1/(1+C_b)
    ### update b_p
    grad_b = matrix(temp6*(delta - (1-delta) * C_b * (1-C_lbeta) * temp3),1) %*% X_p # length J
    j_b = which.max(grad_b)
    b_p[j_b] = b_p[j_b] + epsilon
    ### update beta_p
    grad_beta = matrix(temp5 - temp4,1) %*% W_p # length M
    j_beta = which.max(grad_beta)
    beta_p[j_beta] = beta_p[j_beta] + epsilon

    return(list(b_p = b_p, beta_p = beta_p))
  }

initialization_cox <-
  function(X_u=NULL, W_u=NULL, time, delta){
    data_init = as.data.frame(cbind(time, delta, X_u, W_u))
    if(is.null(X_u) & is.null(W_u)){
      fit = glm(delta~1, family="binomial")
      surv = survival::survfit(Surv(time, delta)~1)
      survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
      inits = list(itct=coef(fit), b_u=NULL, beta_u=NULL, survprob=survprob)
    } else if(is.null(X_u)){
      colnames(data_init) <- c("time","delta",paste0("lat",1:ncol(W_u)))
      formula_lat = as.formula(paste("Surv(time, delta) ~", paste(paste0("lat",1:ncol(W_u)),collapse=" + ")))
      fit1 = survival::coxph(formula_lat, data = data_init)
      surv = survival::survfit(fit1)
      survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
      fit2 = glm(delta~1, family="binomial")
      inits = list(itct=coef(fit2), b_u=NULL, beta_u=coef(fit1), survprob=survprob)
    } else if(is.null(W_u)){
      colnames(data_init) <- c("time","delta",paste0("inc",1:ncol(X_u)))
      formula_inc = as.formula(paste("delta ~", paste(paste0("inc",1:ncol(X_u)),collapse=" + ")))
      fit = glm(formula_inc, family="binomial", data = data_init)
      surv = survival::survfit(Surv(time, delta)~1)
      survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
      itct = coef(fit)[1]
      b_u = coef(fit)[-1]
      inits = list(itct=itct, b_u=b_u, beta_u=NULL, survprob=survprob)
    } else{
      colnames(data_init) <- c("time","delta",paste0("inc",1:ncol(X_u)), paste0("lat",1:ncol(W_u)))
      formula_lat = as.formula(paste("Surv(time, delta) ~", paste(paste0("lat",1:ncol(W_u)),collapse=" + ")))
      formula_inc = as.formula(paste("delta ~", paste(paste0("inc",1:ncol(X_u)),collapse=" + ")))
      fit1 = glm(formula_inc, family="binomial", data = data_init)
      itct = coef(fit1)[1]
      b_u = coef(fit1)[-1]
      fit2 = survival::coxph(formula_lat, data = data_init)
      beta_u=coef(fit2)
      surv = survival::survfit(fit2)
      survprob = sapply(time, function(x) summary(surv, times = x, extend=T)$surv)
      inits = list(itct=itct, b_u=b_u, beta_u=beta_u, survprob=survprob)
    }
    return(inits)
  }

get_cox_lambda_max <- function (x, time, event, offset = rep(0, nrow(x)))
{
  N <- nrow(x)
  beta_w <- offset
  beta_w <- beta_w - mean(beta_w)

  event_time = time[event==1]
  uniq_event_time = unique(event_time)
  n0 = length(uniq_event_time)
  I0 = matrix(time, ncol = 1) %*% matrix(1, 1, n0) >=
    matrix(1, N, 1) %*% matrix(uniq_event_time, nrow = 1)    # N*n0
  d_j = as.numeric(table(event_time)[rank(uniq_event_time)])

  exp_beta_w = exp(beta_w)
  sum0 = pmax(1e-50, matrix(exp_beta_w,1) %*% I0)   # 1*n0
  sum1 = t(x) %*% diag(as.numeric(exp_beta_w)) %*% I0   # ncol(x) * n0
  null_grad = colSums(x[event==1,,drop=FALSE])- rowSums(sum1 %*% diag(as.numeric(d_j/sum0))) # ncol(x) * 1
  g <- abs(null_grad)
  lambda_max <- max(g)/N
  return(lambda_max)
}

initialization_parametric <-
  function(X_u=NULL, W_u=NULL, time, delta, model="weibull"){
    if(!is.null(X_u)){
      fit1 = glm(as.formula(paste("delta ~", paste(colnames(as.data.frame(X_u)),collapse=" + "))),
                 data =as.data.frame(X_u), family = binomial(link="logit"))
      b_u = fit1$coefficients[-1]
      itct = fit1$coefficients[1]
    } else{
      itct = log(mean(delta)/(1-mean(delta)))
      b_u = NULL
    }
    if(!is.null(W_u)) {
      fit2 = survival::coxph(as.formula(paste("Surv(time, delta) ~", paste(colnames(as.data.frame(W_u)),collapse=" + "))),
                             data =as.data.frame(W_u),ties = "breslow")
      beta_u <- coef(fit2)
    } else beta_u = NULL
    if(model=="weibull"){
      alpha = uniroot(function(a) # MOM estimates
        log(gamma(1+2/a))-2*log(gamma(1+1/a))-log(var(time)+(mean(time))^2)+2*log(mean(time)),
        c(0.01,10))$root
      lambda = gamma(1+1/alpha)/mean(time)
    }else if(model=="exponential"){
      lambda = 1/mean(time)
    }
    inits = list(itct=itct, b_u=b_u, beta_u=beta_u, lambda=lambda)
    if(model=="weibull") inits$alpha = alpha
    return(inits)
  }

inits_check <-
  function(model, N, penalty.factor.inc, penalty.factor.lat, inits){
    if(is.null(inits$itct)){
      warning("Initial value of intercept is missing. Initial values were generated by the program.")
      return(NULL)
    }
    if(model=="cox"){
      if(is.null(inits$survprob)){
        warning("Initial value of survprob is missing. Initial values were generated by the program.")
        return(NULL)
      }else if(length(inits$survprob)!=N){
        warning("Initial value of survprob has incorrect dimension. Initial values were generated by the program.")
        return(NULL)
      }
    }
    if(model %in% c("weibull","exponential")){
      if(is.null(inits$lambda)){
        warning("Initial value of lambda is missing. Initial values were generated by the program.")
        return(NULL)
      }else if(inits$lambda<=0){
        warning("Initial value of lambda should be positive. Initial values were generated by the program.")
        return(NULL)
      }
    }
    if(model=="weibull"){
      if(is.null(inits$alpha)){
        warning("Initial value of alpha is missing. Initial values were generated by the program.")
        return(NULL)
      }else if(inits$alpha<=0){
        warning("Initial value of alpha should be positive. Initial values were generated by the program.")
        return(NULL)
      }
    }
    if(any(penalty.factor.inc==0))
      if(is.null(inits$b_u)){
        warning("Initial value of b_u is missing. Initial values were generated by the program.")
        return(NULL)
      } else if(length(inits$b_u)!=sum(penalty.factor.inc==0)){
        warning("Initial value of b_u has incorrect dimension. Initial values were generated by the program.")
        return(NULL)
      }
    if(any(penalty.factor.lat==0))
      if(is.null(inits$beta_u)){
        warning("Initial value of beta_u is missing. Initial values were generated by the program.")
        return(NULL)
      }else if(length(inits$beta_u)!=sum(penalty.factor.lat==0)){
        warning("Initial value of beta_u has incorrect dimension. Initial values were generated by the program.")
        return(NULL)
      }
    return(inits)
  }

l1_grad_b <-
  function(theta, itct, b_u, X_u , X_p, pir, mu){
    N = nrow(X_p)
    J = ncol(X_p)
    b_p_ext = theta
    b_p = b_p_ext[1:J]-b_p_ext[(J+1):(2*J)]
    b_nonzero = which(b_p!=0)
    if(!is.null(X_u)) bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    pix = 1/(1+exp(-bx))
    grad = matrix(pir-pix,1) %*% X_p
    return(c(-grad+N*mu, grad+N*mu))
  }

l1_gradient_inc <-
  function(theta, b_p, X_u , X_p, pir, mu) # incidence
  {
    itct = theta[1]
    if(!is.null(X_u)) b_u = theta[2:(ncol(X_u)+1)]
    b_nonzero = which(b_p!=0)
    if(!is.null(X_u)) bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    pix = 1/(1+exp(-bx))
    grad1 = sum(pir-pix)
    if(!is.null(X_u)) grad2 = matrix(pir-pix,1) %*% X_u else grad2=NULL
    return(-c(grad1, grad2))
  }

l1_negloglik_inc_b <-
  function(theta, itct, b_u, X_u , X_p, pir, mu) # incidence
  {
    N = nrow(X_p)
    J = ncol(X_p)
    b_p_ext = theta
    b_p = b_p_ext[1:J]-b_p_ext[(J+1):(2*J)]
    b_nonzero = which(b_p!=0)
    if(!is.null(X_u)) bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    ll = sum(pir*bx - log(1+exp(bx)))-N*mu*sum(b_p_ext)
    return(-ll)
  }

l1_negloglik_inc <-
  function(theta, b_p, X_u , X_p, pir, mu) # incidence
  {
    itct = theta[1]
    if(!is.null(X_u)) b_u = theta[2:(ncol(X_u)+1)]
    N = nrow(X_p)
    b_nonzero = which(b_p!=0)
    if(!is.null(X_u)) bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    ll = sum(pir*bx - log(1+exp(bx)))-N*mu*sum(abs(b_p))
    return(-ll)
  }

nonparametric_time_generator <-
  function (N, beta_w, maxT = 20, knots = 8)
  {
    time <- 0:maxT
    k <- c(0, sort(sample(time[2:maxT], size = knots, replace = FALSE)), maxT)
    heights <- c(0, sort(1-pmin(rexp(knots, rate = 10),1-1e-15)), 1)
    # heights <- c(0, sort(runif(knots)), 1)
    tk <- merge(data.frame(time), data.frame(time = k, heights),
                by = "time", all = FALSE)
    MonotonicSpline <- stats::splinefun(x = tk$time, y = tk$heights,
                                        method = "hyman")
    t = rep(NA, N)
    for(i in 1:N){
      u = runif(1)
      t[i] = uniroot(function(a) (1-MonotonicSpline(a))^exp(beta_w[i])-u, c(0,maxT))$root
    }
    return(t)
  }

mcp_penalty <-
  function(b_p, gamma, lambda){
    idx = which(b_p !=0 & abs(b_p) <= gamma*lambda)
    s1 = lambda*sum(abs(b_p[idx]))-sum((b_p[idx])^2)/(2*gamma)
    s2 = gamma*lambda^2/2 * sum(abs(b_p) > gamma*lambda)
    return(s1 + s2)
  }

mcp_scad_cd_upd_inc <-
  function(penalty, pir, pi_x, x_l, bl, gamma, lambda){
    d1 = mean((pir - pi_x)*x_l)
    vl = mean(pi_x*(1-pi_x)*x_l^2)
    zl = d1 + vl * bl
    if(penalty == "MCP"){
      if(abs(bl) <= gamma*lambda) return(soft(zl, lambda)/(vl-1/gamma))
      else return(ifelse(vl==0,0,zl/vl))
    }  else if(penalty == "SCAD"){
      if(abs(bl)<=lambda) return(ifelse(vl==0,0,soft(zl, lambda)/vl))
      else if(abs(bl)<=gamma*lambda) return(soft(zl, gamma*lambda/(gamma-1))/(vl-1/(gamma-1)))
      else return(ifelse(vl==0,0,zl/vl))
    }
  }

mcp_scad_cd_upd_lat <-
  function(penalty, exp_beta_w_p, w_l, betal, Zl, I0, d_j, gamma, lambda){
    N = length(exp_beta_w_p)
    sum0 = pmax(1e-50, matrix(exp_beta_w_p,1) %*% I0)   # 1*n0
    sum1 = matrix(exp_beta_w_p * w_l,1) %*% I0   # 1*n0
    sum2 = matrix(exp_beta_w_p * w_l^2,1) %*% I0   # 1*n0
    d1 = sum(Zl - d_j * sum1/sum0)/N
    vl = -sum(d_j * (sum1^2/sum0^2 - sum2/sum0))/N
    zl = d1 + vl * betal
    if(penalty == "MCP"){
      if(abs(betal) <= gamma*lambda) return(soft(zl, lambda)/(vl-1/gamma))
      else return(ifelse(vl==0,0,zl/vl))
    }  else if(penalty == "SCAD"){
      if(abs(betal)<=lambda) return(ifelse(vl==0,0,soft(zl, lambda)/vl))
      else if(abs(betal)<=gamma*lambda) return(soft(zl, gamma*lambda/(gamma-1))/(vl-1/(gamma-1)))
      else return(ifelse(vl==0,0,zl/vl))
    }
  }

mcp_scad_gradient_inc <-
  function(theta, b_p, X_u , X_p, pir, CAP) # incidence
  {
    itct = theta[1]
    if(!is.null(X_u)) b_u = theta[2:(ncol(X_u)+1)]
    b_nonzero = which(b_p!=0)
    if(is.null(X_u)) bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    pix = 1/(1+exp( pmin(CAP, -bx) ))
    grad1 = sum(pir-pix)
    if(!is.null(X_u)) grad2 = matrix(pir-pix,1) %*% X_u else grad2 = NULL
    return(-c(grad1, grad2))
  }

mcp_scad_gradient_lat <-
  function(theta, beta_p, W_u, W_p, delta, pir, I0, d_j, CAP){ # latency
    beta_u = theta
    N = nrow(W_p)
    beta_nonzero = which(beta_p!=0)
    beta_w = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir
    sum0 = pmax(1e-50, matrix(exp_beta_w_p,1) %*% I0)   # 1*n0
    sum1 = t(W_u) %*% diag(as.numeric(exp_beta_w_p)) %*% I0   # ncol(W_u) * n0
    du = colSums(W_u[delta==1,,drop=FALSE])- rowSums(sum1 %*% diag(as.numeric(d_j/sum0))) # ncol(W_u) * 1
    return(-du)
  }

mcp_scad_negloglik_inc <-
  function(theta, b_p, X_u , X_p, pir, CAP) # incidence
  {
    itct = theta[1]
    if(!is.null(X_u)) b_u = theta[2:(ncol(X_u)+1)]
    N = nrow(X_p)
    b_nonzero = which(b_p!=0)
    if(is.null(X_u)) bx = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else bx = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    ll = sum(pir*bx - log(1+exp(pmin(CAP, bx))))
    return(-ll)
  }

mcp_scad_negloglik_lat <-
  function(theta, beta_p, W_u, W_p, delta, pir, I0, d_j, CAP){ # latency
    beta_u = theta
    N = nrow(W_p)
    beta_nonzero = which(beta_p!=0)
    beta_w = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    exp_beta_w_p = exp( pmin(CAP, beta_w) ) * pir
    sum0 = pmax(1e-50, matrix(exp_beta_w_p,1) %*% I0)  # 1*n0
    ll = sum(beta_w[delta==1]) - sum(d_j * log(sum0))
    return(-ll)
  }

scad_penalty <-
  function(b_p, gamma, lambda){
    idx1 = which(b_p !=0 & abs(b_p) <= lambda)
    idx2 = which(abs(b_p) > lambda & abs(b_p)<=gamma*lambda)
    s1 = lambda*sum(abs(b_p[idx1]))
    s2 = (gamma*lambda*sum(abs(b_p[idx2])) - 0.5*(sum((b_p[idx2])^2 + lambda^2)))/(gamma-1)
    s3 = lambda^2*(gamma^2-1)/(2*(gamma-1)) * sum(abs(b_p) > gamma*lambda)
    return(s1 + s2 + s3)
  }

self_scale <-
  function(X, scale){
    if(is.null(X)) return(NULL)
    if(ncol(X)==0) return(NULL)
    if (scale) {
    n = nrow(X)
    Xs = apply(X, 2, function(x) if(sd(x)==0) return(rep(0,n)) else return(scale(x)))
    } else {
      Xs = X
    }
    return(Xs)
  }

soft <-
  function(z, gamma){
    if(gamma >= abs(z)) return(0)
    else return(ifelse(z>0, z-gamma, z+gamma))
  }

weib_EM <-
  function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, mu_inc, mu_lat,
           inits=NULL, nIter = 100, tol = 1e-4)
  {
    # mu: penalty parameter
    # tol: difference between log-likelihood
    N = length(time)
    J = ncol(X_p) # number of penalized incidence covariates
    M = ncol(W_p) # number of penalized latency covariates

    ######## initialization ########
    step = 1
    b_p <- rep(0,J)
    beta_p <- rep(0,M)
    b_p_ext <- rep(0,2*J)
    beta_p_ext <- rep(0,2*M)
    if(is.null(inits)) inits = initialization_parametric(X_u, W_u, time, delta, model="weibull")
    itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u
    lambda = inits$lambda; alpha = inits$alpha
    log_alpha = log(max(alpha, 1e-15))
    log_lambda = log(max(lambda, 1e-15))
    pir = rep(1, N)
    llp1_0 <- llp2_0 <-0
    lik_inc <- lik_lat <- c()
    b_p_path <- beta_p_path <- b_u_path <- beta_u_path <- itct_path <- alpha_path<-lambda_path<-NULL
    conv1 <- conv2 <- FALSE

    ######## loop ########

    repeat{

      #### E-step
      if(is.null(X_u)) b_x = itct + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
      else b_x = itct + X_u %*% b_u + X_p[,b_p!=0,drop=FALSE] %*% b_p[b_p!=0]
      if(is.null(W_u)) beta_w = W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
      else beta_w = W_u %*% beta_u + W_p[,beta_p!=0,drop=FALSE] %*% beta_p[beta_p!=0]
      pir[delta==0] = 1/(1+exp(-b_x[delta==0]+
                                 (lambda * time[delta==0])^alpha * exp(beta_w[delta==0])))
      #### update penalized parameters
      if (!conv1){
        b_p_ext = optim(par = b_p_ext, fn = l1_negloglik_inc_b, gr = l1_grad_b,
                        itct = itct, b_u = b_u, X_u = X_u, X_p = X_p,
                        pir = pir, mu = mu_inc, method = 'L-BFGS-B', lower = rep(0, 2*J))$par
        b_p = b_p_ext[1:J]-b_p_ext[(J+1):(2*J)]
      }

      if (!conv2){
        beta_p_ext = optim(par = beta_p_ext, fn = weib_negloglik_lat_beta, gr = weib_grad_beta,
                           alpha = alpha, lambda = lambda,
                           beta_u = beta_u, W_u = W_u, W_p = W_p, time = time, delta = delta, pir = pir, mu = mu_lat,
                           method = 'L-BFGS-B', lower = rep(0, 2*M))$par
        beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
      }

      #### update nonpenalized parameters
      if (!conv1){
        out_inc = optim(par = c(itct,b_u), fn = l1_negloglik_inc, gr = l1_gradient_inc,
                        b_p = b_p, X_u =X_u, X_p = X_p, pir = pir, mu = mu_inc, method="BFGS") # incidence
        itct = out_inc$par[1]
        if(!is.null(X_u)) b_u = out_inc$par[2:(ncol(X_u)+1)]
        llp1_1 = -out_inc$value
      }
      lik_inc = c(lik_inc, llp1_1)

      if (!conv2){
        out_lat = optim(par = c(log_alpha, log_lambda, beta_u),
                        fn = weib_negloglik_lat, gr = weib_gradient_lat,
                        beta_p=beta_p, W_u=W_u, W_p=W_p, time=time, delta=delta, pir=pir, mu=mu_lat,
                        method="BFGS") # latency
        log_alpha = out_lat$par[1]
        alpha = exp(log_alpha)
        log_lambda = out_lat$par[2]
        lambda = exp(log_lambda)
        if(!is.null(W_u)) beta_u = out_lat$par[3:(ncol(W_u)+2)] # unpenalized incidence
        llp2_1 = -out_lat$value
      }
      lik_lat = c(lik_lat, llp2_1)

      ## record updated parameters
      itct_path = c(itct_path, itct)
      if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
      b_p_path = rbind(b_p_path, b_p)
      if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)
      beta_p_path = rbind(beta_p_path, beta_p)
      alpha_path = c(alpha_path, alpha)
      lambda_path = c(lambda_path, lambda)
      if (!conv1 & abs(llp1_1- llp1_0)< tol) conv1 <- TRUE
      if (!conv2 & abs(llp2_1- llp2_0)< tol) conv2 <- TRUE
      if (step > 1 & (conv1 & conv2) | step >= nIter) { #
        break
      }
      llp1_0 <- llp1_1
      llp2_0 <- llp2_1
      step <- 1 + step
    }

    ######## output ########
    output <- list(b_p_path = b_p_path, b_u_path = b_u_path, itct_path = itct_path,
                   beta_p_path = beta_p_path, beta_u_path = beta_u_path,
                   alpha_path = alpha_path, lambda_path = lambda_path,
                   lik_inc = lik_inc, lik_lat = lik_lat)
    output
  }

weib_grad_beta <-
  function(theta, alpha, lambda, beta_u, W_u, W_p, time, delta, pir, mu){
    N = nrow(W_p)
    M = ncol(W_p)
    beta_p_ext = theta
    beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    temp1 = pir*(lambda*time)^alpha*exp(betaw)
    grad = matrix(delta-temp1,1)%*% W_p
    return(c(-grad+ N*mu, grad+ N*mu))
  }

weib_gradient_lat <-
  function(theta, beta_p, W_u, W_p, time, delta, pir, mu){ # latency
    log_alpha = theta[1]
    log_lambda = theta[2]
    alpha = exp(log_alpha)
    lambda = exp(log_lambda)
    if(!is.null(W_u)) beta_u = theta[3:(ncol(W_u)+2)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    temp1 = pir*(lambda*time)^alpha*exp(betaw)
    temp2 = delta-temp1
    grad1 = sum(delta+log(pmax(lambda * time,1e-100))*temp2*alpha)
    grad2 = sum(alpha*temp2)
    if(!is.null(W_u)) grad3 = matrix(temp2,1) %*% W_u else grad3=NULL
    return(-c(grad1, grad2, grad3))
  }

weib.cure.update <-
  function(alpha, lambda, b_p, beta_p, b_u, itct, beta_u, X_u , X_p, W_u, W_p, time, delta, epsilon)
  {
    b_nonzero = which(b_p!=0)
    beta_nonzero = which(beta_p!=0)
    if(!is.null(X_u)) C_b = exp(itct + X_u %*% b_u + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    else C_b = exp(itct + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    if(!is.null(W_u)) C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    else C_beta = exp(W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    C_lbeta = exp(-(lambda * time)^alpha * C_beta)
    temp2 = (lambda * time)^alpha * C_beta
    temp3 = 1 / (1 + C_b*C_lbeta)
    temp4 = (1-delta) * C_b * C_lbeta * temp2 * temp3
    temp5 = delta * (1-temp2)
    temp6 = 1/(1+C_b)
    ### update b_p
    grad_b = matrix(temp6*(delta - (1-delta) * C_b * (1-C_lbeta) * temp3),1) %*% X_p # length J
    j_b = which.max(grad_b)
    b_p[j_b] = b_p[j_b] + epsilon
    ### update beta_p
    grad_beta = matrix(temp5 - temp4,1) %*% W_p # length M
    j_beta = which.max(grad_beta)
    beta_p[j_beta] = beta_p[j_beta] + epsilon

    return(list(b_p = b_p, beta_p = beta_p))
  }

weibull.cure <-
  function(X_u=NULL, X_p, W_u=NULL, W_p, time, delta, epsilon = 0.001, tol = 1e-05,
           nIter=1e4, inits=NULL, verbose=TRUE)
  {

    # X_u: N by nonp, non-penalized covariate matrix associated with incidence
    # X_p: N by J, penalized covariate matrix associated with incidence
    # W_u: N by nonp, non-penalized covariate matrix associated with latency
    # W_p:  N by M, penalized covariate matrix associated with incidence
    # time: vector of length N, observed survival time
    # delta: vector of length N, censoring status (not censored = 0)
    # epsilon: incremental size
    # tol: difference between log-likelihood

    N = length(time)
    J = ncol(X_p) # number of penalized incidence covariates
    M = ncol(W_p) # number of penalized latency covariates
    X_p = cbind(X_p, -X_p) # N by 2J
    W_p = cbind(W_p, -W_p) # N by 2M

    ######## initialization ########
    step = 1
    b_p = rep(0, 2*J) # penalized incidence
    beta_p = rep(0, 2*M) # penalized latency
    if(is.null(inits)) inits = initialization_parametric(X_u, W_u, time, delta, model="weibull")
    itct = inits$itct; b_u = inits$b_u; beta_u = inits$beta_u
    lambda = inits$lambda; alpha = inits$alpha
    log_alpha = log(max(alpha, 1e-15))
    log_lambda = log(max(lambda, 1e-15))
    LL0 = 0

    b_p_path <- beta_p_path <- alpha_path <- lambda_path <- b_u_path <- beta_u_path <- itct_path <- NULL
    logLikelihood <- numeric()

    ######## loop ########

    repeat{

      #### update penalized parameters
      upd = weib.cure.update(alpha, lambda, b_p, beta_p, b_u, itct, beta_u, X_u , X_p, W_u, W_p, time, delta, epsilon)
      b_p = upd$b_p
      beta_p = upd$beta_p
      b_p_path = rbind(b_p_path, b_p)
      beta_p_path = rbind(beta_p_path, beta_p)

      #### update other parameters
      out = optim(par = c(log_alpha, log_lambda, b_u, itct, beta_u), fn = weib.cure.negloglik, gr = weib.cure.gradient,
                  b_p = b_p, beta_p = beta_p, X_u = X_u, X_p = X_p, W_u = W_u, W_p = W_p,
                  time = time, delta = delta, method="BFGS")
      log_alpha = out$par[1]
      alpha = exp(log_alpha)
      log_lambda = out$par[2]
      lambda = exp(log_lambda)
      if(!is.null(X_u)) {
        b_u = out$par[3:(ncol(X_u)+2)] # unpenalized incidence
        itct = out$par[ncol(X_u)+3]
        if(!is.null(W_u)) beta_u = out$par[(ncol(X_u)+4):(ncol(X_u)+ncol(W_u)+3)] # unpenalized latency
      } else{
        itct = out$par[3]
        if(!is.null(W_u)) beta_u = out$par[4:(ncol(W_u)+3)] # unpenalized latency
      }

      alpha_path = c(alpha_path, alpha)
      lambda_path = c(lambda_path, lambda)
      if(!is.null(X_u)) b_u_path = rbind(b_u_path, b_u)
      itct_path = c(itct_path, itct)
      if(!is.null(W_u)) beta_u_path = rbind(beta_u_path, beta_u)

      LL1 = -out$value
      logLikelihood = c(logLikelihood, LL1)

      if(verbose & step%%1000==0) cat("step=", step, "\n")
      if (step > 1 && (abs(LL1 - LL0) < tol | step >= nIter)) { #
        break
      }
      LL0 <- LL1
      step <- 1 + step
    }

    ######## output ########
    b_p_path = b_p_path[,1:J] - b_p_path[,(J+1):(2*J)]
    beta_p_path = beta_p_path[,1:M] - beta_p_path[,(M+1):(2*M)]
    output <- list(b_p_path = b_p_path, beta_p_path = beta_p_path, alpha_path = alpha_path,
                   lambda_path = lambda_path, b_u_path = b_u_path, itct_path = itct_path,
                   beta_u_path = beta_u_path, logLikelihood = logLikelihood)
    output
  }

weib.cure.gradient <-
  function(theta, b_p, beta_p, X_u , X_p, W_u, W_p, time, delta)
  {
    N = length(time)
    alpha = exp(theta[1])
    lambda = exp(theta[2])
    if(!is.null(X_u)) {
      b_u = theta[3:(ncol(X_u)+2)] # unpenalized incidence
      itct = theta[ncol(X_u)+3]
      if(!is.null(W_u)) beta_u = theta[(ncol(X_u)+4):(ncol(X_u)+ncol(W_u)+3)] # unpenalized latency
    } else{
      itct = theta[3]
      if(!is.null(W_u)) beta_u = theta[4:(ncol(W_u)+3)] # unpenalized latency
    }
    b_nonzero = which(b_p!=0)
    beta_nonzero = which(beta_p!=0)
    if(!is.null(X_u)) C_b = exp(itct + X_u %*% b_u + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    else C_b = exp(itct + X_p[,b_nonzero,drop=FALSE] %*% b_p[b_nonzero])
    if(!is.null(W_u)) C_beta = exp(W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    else C_beta = exp(W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero])
    C_lbeta = exp(-(lambda * time)^alpha * C_beta)
    temp1 = log (pmax(lambda*time,rep(1e-100,N)))
    temp2 = (lambda * time)^alpha * C_beta
    temp3 = 1 / (1 + C_b*C_lbeta)
    temp4 = (1-delta) * C_b * C_lbeta * temp2 * temp3
    temp5 = delta * (1-temp2)
    temp6 = 1/(1+C_b)
    temp7 = temp6*(delta - (1-delta) * C_b * (1-C_lbeta) * temp3)
    grad1 = sum(delta * (1/alpha + temp1 - temp2 * temp1) - temp4 * temp1) * alpha
    grad2 = alpha * sum(temp5 - temp4)
    if(!is.null(X_u)) grad3 = matrix(temp7,1) %*% X_u else grad3=NULL
    grad4 = sum(temp7) # intercept
    if(!is.null(W_u)) grad5 = matrix(temp5 - temp4, 1) %*% W_u else grad5=NULL
    return(-c(grad1, grad2, grad3, grad4, grad5))
  }

weib.cure.negloglik <-
  function(theta, b_p, beta_p, X_u , X_p, W_u, W_p, time, delta)
  {
    N = length(time)
    alpha = exp(theta[1])
    lambda = exp(theta[2])
    if(!is.null(X_u)) {
      b_u = theta[3:(ncol(X_u)+2)] # unpenalized incidence
      itct = theta[ncol(X_u)+3]
      if(!is.null(W_u)) beta_u = theta[(ncol(X_u)+4):(ncol(X_u)+ncol(W_u)+3)] # unpenalized latency
    } else{
      itct = theta[3]
      if(!is.null(W_u)) beta_u = theta[4:(ncol(W_u)+3)] # unpenalized latency
    }
    b_nonzero = which(b_p!=0)
    beta_nonzero = which(beta_p!=0)
    if(!is.null(X_u)) logC_b = itct + X_u %*% b_u + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    else logC_b = itct + X_p[,b_nonzero, drop=FALSE] %*% b_p[b_nonzero]
    if(!is.null(W_u)) logC_beta = W_u %*% beta_u + W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero]
    else logC_beta = W_p[,beta_nonzero,drop=FALSE] %*% beta_p[beta_nonzero]
    temp1 = 1/(1+exp(-logC_b)) # exp(logC_b) / (1+exp(logC_b)) #Mar17
    logC_lbeta = -(lambda * time)^alpha * exp(logC_beta)
    ll1 = delta * (log(temp1) + theta[2] + theta[1] +
                     (alpha-1)*log(pmax(lambda * time,rep(1e-100,N))) +
                     logC_beta + logC_lbeta)
    ll2 = (1-delta) * log(pmax(1 - temp1 * (1-exp(logC_lbeta)),rep(1e-100,N)))
    return(-sum(ll1+ll2))
  }

weib_negloglik_lat_beta <-
  function(theta, alpha, lambda, beta_u, W_u, W_p, time, delta, pir, mu){ # latency
    N = nrow(W_p)
    M = ncol(W_p)
    beta_p_ext = theta
    beta_p = beta_p_ext[1:M]-beta_p_ext[(M+1):(2*M)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    ll1 = delta*( log(alpha)+log(lambda) +(alpha-1)*log(pmax(lambda * time,1e-100))+betaw )
    ll2 = -pir*(lambda*time)^alpha*exp(betaw)
    return(-sum(ll1+ll2)+N*mu*sum(beta_p_ext))
  }

weib_negloglik_lat <-
  function(theta, beta_p, W_u, W_p, time, delta, pir, mu){ # latency
    N = nrow(W_p)
    log_alpha = theta[1]
    log_lambda = theta[2]
    alpha = exp(log_alpha)
    lambda = exp(log_lambda)
    if(!is.null(W_u)) beta_u = theta[3:(ncol(W_u)+2)]
    beta_nonzero = which(beta_p!=0)
    if(!is.null(W_u)) betaw = W_u %*% beta_u + W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    else betaw = W_p[,beta_nonzero, drop=FALSE] %*% beta_p[beta_nonzero]
    ll1 = delta*( log_alpha+log_lambda +(alpha-1)*log(pmax(lambda * time,1e-100))+betaw )
    ll2 = -pir*(lambda*time)^alpha*exp(betaw)
    return(-sum(ll1+ll2)+N*mu*sum(abs(beta_p)))
  }

# If B = NULL use exponential; otherwise use uniform with B=upper limit
sim_cure <-
  function(n, mu=1, censor.mu=3, B=NULL, Reps=10000, seed=NULL) {
    susceptible.prop<-numeric(length=Reps)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    for (i in 1:Reps) {
      surv<-rexp(n, 1/mu)
      if (is.null(B)) {
        u<-rexp(n, 1/censor.mu) ## rate is 1/mu
      } else {
        u<-runif(n,0,B)
      }
      time<-ifelse(surv<=u,surv,u)
      censor<-ifelse(surv<=u,1,0)
      fit<-survival::survfit(survival::Surv(time,censor)~1)
      susceptible.prop[i]<-1-min(fit$surv)
    }
    susceptible.prop
  }

AUC_msi <- function(cure_cutoff = 5, b_u_hat=NULL, itct_hat, b_p_hat, testing_delta, testing_time,
                    X_u=NULL, X_p){
  testing_n = length(testing_time)
  v = rep(0, testing_n)
  y = rep(999, testing_n)
  y[testing_time>cure_cutoff] = 0
  y[testing_time<=cure_cutoff & testing_delta==1] = 1
  v[y<2] = 1
  if(all(b_p_hat==0)) {
    if(is.null(X_u)) xb = rep(itct_hat, testing_n) else xb = itct_hat+X_u %*% b_u_hat
  } else{
    if(is.null(X_u)) xb = itct_hat + X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0]
    else xb = itct_hat+X_u %*% b_u_hat + X_p[,b_p_hat!=0,drop=FALSE] %*% b_p_hat[b_p_hat!=0]
  }
  p_hat = 1/(1+exp(-xb))
  temp = v*y + (1-v)*p_hat
  temp1 = temp[order(p_hat, decreasing = T)]
  temp_f = v*(1-y) + (1-v) * (1-p_hat)
  temp1f = temp_f[order(p_hat, decreasing = T)]
  TPR = c(0, cumsum(temp1)/cumsum(temp1)[testing_n])
  FPR = c(0, cumsum(temp1f)/cumsum(temp1f)[testing_n])
  height = (TPR[-1]+TPR[-length(TPR)])/2
  width = diff(FPR)
  auc_msi = sum(height*width)
  return(auc_msi)
}


C.stat <-
function (cure_cutoff = 5, b_u_hat = NULL, itct_hat, b_p_hat,
          beta_u_hat = NULL, beta_p_hat, testing_delta, testing_time,
          X_u = NULL, X_p, W_u = NULL, W_p)
{
  C_csw_num = 0
  C_csw_denom = 0
  testing_n = length(testing_time)
  v = rep(0, testing_n)
  y = rep(999, testing_n)
  y[testing_time > cure_cutoff] = 0
  y[testing_time <= cure_cutoff & testing_delta == 1] = 1
  v[y < 2] = 1
  if (all(b_p_hat == 0)) {
    if (is.null(X_u))
      p_hat = 1/(1 + exp(-itct_hat))
    else p_hat = 1/(1 + exp(-itct_hat - X_u %*% b_u_hat))
  }
  else {
    if (is.null(X_u))
      p_hat = 1/(1 + exp(-itct_hat - X_p[, b_p_hat != 0,
                                         drop = FALSE] %*% b_p_hat[b_p_hat != 0]))
    else p_hat = 1/(1 + exp(-itct_hat - X_u %*% b_u_hat -
                              X_p[, b_p_hat != 0, drop = FALSE] %*% b_p_hat[b_p_hat !=
                                                                              0]))
  }
  temp = v * y + (1 - v) * as.vector(p_hat)
  if (all(beta_p_hat == 0)) {
    if (is.null(W_u))
      W_beta = rep(0, testing_n)
    else W_beta = W_u %*% beta_u_hat
  }
  else {
    if (is.null(W_u))
      W_beta = W_p[, beta_p_hat != 0, drop = FALSE] %*%
        beta_p_hat[beta_p_hat != 0]
    else W_beta = W_u %*% beta_u_hat + W_p[, beta_p_hat !=
                                             0, drop = FALSE] %*% beta_p_hat[beta_p_hat != 0]
  }
  for (i in 1:testing_n) for (j in 1:testing_n) {
    if (j == i | !testing_delta[i] | testing_time[i] > testing_time[j])
      next
    I_ij = testing_time[i] < testing_time[j] | (testing_time[i] ==
                                                  testing_time[j] & !testing_delta[j])
    if (!I_ij)
      next
    if (W_beta[i] > W_beta[j])
      C_csw_num = C_csw_num + temp[j]
    C_csw_denom = C_csw_denom + temp[j]
  }
  return(C_csw_num/C_csw_denom)
}
