.clip <- function(val, bound) {
  pmin(pmax(val, bound[1]), bound[2])
}

#------------------
# tmleLogLike -  calculate the quasi-loglikelihood loss
#		for the TMLE updating step
# Here are using a 2-dimensional covariate for fluctuation bc we want to calculate
#	the marginal risk difference, risk ratio and odds ratio
# input: eps (flucutation parameter), outcome Y, initial estimates of the Qbar(A,W)
#		the clever covariate under no exposure H.0W, the clever covariate under exposure H.1W, weights
# returns: negative quasi-loglikelihood loss
#-----------------
tmleLogLike<- function(eps, Y, Q.AW, H.0W, H.1W, wt=1){

  pi<- plogis( qlogis(Q.AW) + eps[1]*H.0W + eps[2]*H.1W)
  pi[pi==0] <- .Machine$double.neg.eps		# prevent from taking the log of 0 or 1
  pi[pi==1]<- 1-.Machine$double.neg.eps
  logLike<- sum( wt*( Y*log(pi)+ (1-Y)*log(1-pi)) )

  if(is.na(logLike) | is.nan(logLike) | is.infinite(logLike)){
    eps<- c(0,0)	# set epsilon1 and epsilon2 to zero
    logLike<- 99
  }
  return(-logLike)		# log likelihood loss function.

}

#---------------------------
# tmlegrad - corresponding gradient function.
# 	other optimization routines do not require a gradient
#-----------------------------
tmlegrad<- function(eps, Y, Q.AW, H.0W, H.1W, wt=1){

  pi<- plogis( qlogis(Q.AW) + eps[1]*H.0W + eps[2]*H.1W)
  pi[pi==0] <- .Machine$double.neg.eps
  pi[pi==1]<- 1-.Machine$double.neg.eps
  resid<- wt*(Y-pi) # weighted residuals
  gr<- crossprod( cbind(H.0W, H.1W), resid)

  if( sum(is.na(gr))>0 | sum(is.nan(gr))>0 | sum(is.infinite(gr))>0 ){
    gr<- c(99, 99)
  }
  return(-gr)

}

create_custom_glm <- function(var_name) {
  force(var_name)
  function(Y, X, newX, ...) {
    SL.glm(Y = Y, X = X[, var_name, drop = FALSE],
           newX = newX[, var_name, drop = FALSE], ...)
  }
}

rtmle <- function(W,
                  A,
                  Y,
                  a,
                  b,
                  wt) {

  #browser()

  #-----------------
  # Step 0 -  transform outcome
  Y.tilde <- (Y-a)/(b-a)

  #---------------
  # Step 1 -  initial estimation of conditional mean outcome
  data_A1 <- cbind(A = 1, W)
  data_A0 <- cbind(A = 0, W)
  W_names <- names(W)
  for (var in W_names) {
    assign(paste0("SL.glm.Q.", var), create_custom_glm(var))
    assign(paste0("SL.glm.Q.", var, ".A"), create_custom_glm(c("A", var)))
  }
  Q.SL.lib <- c(paste0("SL.glm.Q.", W_names), paste0("SL.glm.Q.", W_names, ".A"))
  Q.SL.lib <- c("SL.mean", Q.SL.lib)
  Q.SL <- SuperLearner(Y = Y.tilde, X = cbind(A, W),
                       family = gaussian(),
                       SL.library = Q.SL.lib)
  QbarAW <- .clip(as.numeric(Q.SL$SL.predict), c(0.01, 0.99))
  Qbar1W <- .clip(as.numeric(predict(Q.SL, newdata = data_A1)$pred), c(0.01, 0.99))
  Qbar0W <- .clip(as.numeric(predict(Q.SL, newdata = data_A0)$pred), c(0.01, 0.99))
  Qbar <- data.frame(QbarAW, Qbar1W, Qbar0W)

  # g SL library
  for (var in W_names) {
    assign(paste0("SL.glm.g.", var), create_custom_glm(var))
  }
  g.SL.lib <- paste0("SL.glm.g.", W_names)
  g.SL.lib <- c("SL.mean", g.SL.lib)
  g.SL <- SuperLearner(Y = A, X = W,
                       family = binomial(),
                       SL.library = g.SL.lib)
  g1W <- .clip(as.numeric(g.SL$SL.predict), c(0.01, 0.99))

  #----------------------------
  # Step 2 - targeting
  g0W <- 1-g1W

  # calculate clever covariate -> going to do a 2-dim update here
  H.1W <- A/g1W
  H.0W <- (1-A)/g0W  # note not negative here.
  # it's easier to code a 2-dim update when going after Risk Diff, Risk Ratio, Odds Ratio

  # initial parameter estimates are 0
  opt.out <- optim(par=c(0,0), fn=tmleLogLike, gr=tmlegrad, Y=Y.tilde,
                   Q.AW=Qbar$QbarAW, H.0W=H.0W, H.1W=H.1W, wt=wt, method="BFGS",
                   hessian=F, control=list(maxit=500, trace=F))

  eps <- opt.out$par
  #print(paste('two-dim epsilon', eps, sep=' '))

  # update. still on the transformed scale.
  QbarAW <- plogis(qlogis(Qbar$QbarAW)+eps[1]*H.0W + eps[2]*H.1W)
  Qbar0W <- plogis(qlogis(Qbar$Qbar0W)+eps[1]/g0W)
  Qbar1W <- plogis(qlogis(Qbar$Qbar1W)+eps[2]/g1W)

  Q <- data.frame(cbind(QbarAW, Qbar1W, Qbar0W))
  colnames(Q)<- c("QbarAW", "Qbar1W", "Qbar0W")

  # transform back
  Qstar <- Q*(b-a)+a

  #------------------
  # point estimate
  # for CCDesigns, need to average over case-control weighted dist of cov
  # see van der Laan 2008
  Risk1 <- weighted.mean(Qstar$Qbar1W, w=wt)
  Risk0 <- weighted.mean(Qstar$Qbar0W, w=wt)
  RiskDiff <- Risk1 - Risk0
  RiskRatio <- Risk1/Risk0
  OddsRatio <- (Risk1/(1-Risk1))/(Risk0/(1-Risk0))

  # inference
  eic <- Qstar$Qbar1W-Qstar$Qbar0W-mean(Qstar$Qbar1W-Qstar$Qbar0W)+(A/g1W-(1-A)/g0W)*(Y-Qstar$QbarAW)
  se <- sqrt(var(eic)/length(eic))
  lower <- RiskDiff - 1.96*se
  upper <- RiskDiff + 1.96*se

  # AIPW
  aipw_nuisance <- AIPW_nuis$new(Y = Y,
                                 A = A,
                                 mu0 = Qbar0W*(b-a)+a,
                                 mu1 = Qbar1W*(b-a)+a,
                                 raw_p_score = g1W*A+(1-A)*g0W,
                                 verbose = FALSE)$summary()

  data.frame(Risk1,
             Risk0,
             RiskDiff,
             lower,
             upper,
             se,
             AIPW_Risk1 = aipw_nuisance$result["Risk of exposure", "Estimate"],
             AIPW_Risk0 = aipw_nuisance$result["Risk of control", "Estimate"],
             AIPW_RiskDiff = aipw_nuisance$result["Risk Difference", "Estimate"],
             AIPW_lower = aipw_nuisance$result["Risk Difference", "95% LCL"],
             AIPW_upper = aipw_nuisance$result["Risk Difference", "95% UCL"],
             AIPW_se = aipw_nuisance$result["Risk Difference", "SE"])
}
