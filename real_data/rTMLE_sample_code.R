################# 
# Sample R code for estimation of the marginal risk difference, 
#	risk ratio and odds ratio using MLE, TMLE and rare outcomes rTMLE
#
# For further details about the rare outcome TMLE procedure, 
#		see http://biostats.bepress.com/ucbbiostat/paper310/
#
# The code includes weights to correct for case-control sampling (if applicable). 
# Weights are included in the ltmle package (version 0.9-6)
#	http://cran.r-project.org/web/packages/ltmle/index.html
#
# Programmer: Laura Balzer (lbbalzer@hsph.harvard.edu)
#	Please email with questions or additional code regarding SuperLearner and/or
#	 inference via the estimated influence curves
#
# Acknowledgements: Thank you to Dr. Mark van der Laan and Sam Lendle for their 
#	help 
#
# Last update: 08.03.15
###################


#-------------------------
# mle.wt: function to run a simple substitution estimator 
# (a.k.a. parametric MLE, parametric G-comp)
# for case-control designs, this will be weighted logistic regression
# 	input: data, Qform (formula for the outcome regression), wt (weights)
#	output: estimates of RiskDiff, Risk Ratio, OddsRatio
#-------------------------
mle.wt<- function(data, Qform, wt){
		
	# parametric estimation for the conditional risk P(Y=1|A,W)=E(Y|A,W)=Qbar(A,W)
	glm.out<- suppressWarnings(glm(Qform, family='binomial', data=data, weights= wt))
	# ignore warning message

	# To obtain the predicted outcome given (A=1, W) and given (A=0, W)
	# Create two new data frames X1 and X0; reassign the exposure
	X1 = X0 = data
	X1$A<- 1; X0$A<- 0
	
	# get predicted probabilities of outcome, under the exposure and control
	Qbar1W<- predict(glm.out, newdata=X1, type='response')
	Qbar0W<- predict(glm.out, newdata=X0, type='response')

	# for CCDesigns, need to average over case-control weighted dist of cov
	# see van der Laan 2008 - Estimation based on Case-Control designs with known
	#	prevalence probability. IJB
	Risk1<- weighted.mean(Qbar1W, w=wt)
	Risk0<- weighted.mean(Qbar0W, w=wt)
	RiskDiff<- Risk1 - Risk0
	RiskRatio <- Risk1/Risk0
	OddsRatio<- (Risk1/(1-Risk1))/(Risk0/(1-Risk0))		
	
	data.frame(RiskDiff, RiskRatio, OddsRatio)
}
		
#-------------------------
# tmle.wt: function to run tmle with case-control weights
# 	input:  data, Qform (formula for outcome regression), 
#		gform (formula for exposure mechanism), 	wt (weights)
#	output: estimates of RiskDiff, Risk Ratio, OddsRatio
#--------------------------

# please note, the ltmle package now includes weights 
#	http://cran.r-project.org/web/packages/ltmle/index.html
# 	so it's not necesary to program this function yourself.

tmle.wt<- function(data, Qform, gform, wt){
	
	#-------------------
	# Step 1: initial estimation of the conditional risk P(Y=1|A,W)= E(Y|A,W)=Qbar(A,W)
	glm.out<-  suppressWarnings(glm(Qform, family='binomial', data=data, weights= wt))
	# data-adaptive estimation (e.g. SuperLearner) is recommended
	# we are showing parametric estimation for pedagogic purposes

	# To obtain the predicted outcome given the observed exposure & baseline covariates,
	# 	given A=1 and W, given A=0 and W
	# Create two new data frames X1 and X0; reassign the exposure
	X1 = X0 = data
	X1$A<- 1; X0$A<- 0
	
	# get predicted probabilities of outcome, under the observed exposure,
	#	the intervention, and control
	QbarAW<- predict(glm.out, newdata=data, type='response')
	Qbar1W<- predict(glm.out, newdata=X1, type='response')
	Qbar0W<- predict(glm.out, newdata=X0, type='response')

	#----------------------------
	# Step 2: targeting
	
	#parametric estimation of the propensity score P(A=1|W)=g(A=1|W) 
	#	again we recommend using more data-adaptive methods in practice
	glm.out<-  suppressWarnings(glm(gform, family='binomial', data=data, weights=wt))
	g1W<-  predict( glm.out, type='response' ) 
	g0W<- 1-g1W
	
	# calculate the clever covariates (two-dim)
	# this targets E(Y1) and E(Y0) separately, but maps into locally efficient estimators of RD, RR, OR
	H.1W<- data$A/g1W
	H.0W<- (1-data$A)/g0W
	logitUpdate<-  suppressWarnings(glm(data$Y ~ -1 +offset(qlogis(QbarAW)) 
		+ H.0W+ H.1W, family=binomial, weights=wt))
	eps<-logitUpdate$coef
	#print(eps)
	
	#-----------------------
	# Update
	Qbar0W.star<-plogis( qlogis(Qbar0W)+eps[1]/g0W)
	Qbar1W.star<-plogis( qlogis(Qbar1W)+eps[2]/g1W)
	
	#------------------
	# point estimate
	# for CCDesigns, need to average over case-control weighted dist of cov
	# see van der Laan 2008
	Risk1<- weighted.mean(Qbar1W.star, w=wt)
	Risk0<- weighted.mean(Qbar0W.star, w=wt)
	RiskDiff<- Risk1 - Risk0
	RiskRatio <- Risk1/Risk0
	OddsRatio<- (Risk1/(1-Risk1))/(Risk0/(1-Risk0))
	
	data.frame(RiskDiff, RiskRatio, OddsRatio)

}

#=========================================
# For further details about the rare outcome TMLE procedure, 
#		see http://biostats.bepress.com/ucbbiostat/paper310/
#
# The following functions needed for rare outcomes TMLE.
# More succint code is available from Sam Lendle at 
#		https://github.com/lendle/ltmle/blob/rare_outcomes/R/boundedlogistic.R
#=====================================

#------------------------------------
# rtmle.wt function to do rTMLE with  weights
#		input: covariates W, exp A, outcome Y, a (lower bound), b (upper bound), 
#			Qform (initial formula for the outcome regression),
#			gform (formula for the exposure mechanism), weights wt
#		output: estimates of RiskDiff, Risk Ratio, OddsRatio
#----------------------------------

rtmle.wt<- function(W, A, Y, a, b, Qform, gform,  wt){
	
	#-----------------
	# Step 0 -  transform outcome
	Y.tilde <- (Y-a)/(b-a)
	
	#---------------
	# Step 1 -  initial estimation of conditional mean outcome
	# as before, we recommend using more flexible methods (e.g SuperLearner)
	X.model <- model.matrix(Qform, data=data.frame(W,A,Y) )

	optim.out <- optim(par=rep(0, ncol(X.model)), fn=logLike, gr=grad, Y=Y.tilde, X=X.model, 
			wt= wt, method="BFGS", hessian=F,  control = list(maxit = 500, trace=F)) 
			
	beta<- optim.out$par
	X1= X0= X.model; X1[,"A"]=1; X0[,"A"]=0
	QbarAW<- plogis(X.model%*% beta); Qbar1W<- plogis(X1%*%beta); Qbar0W<- plogis(X0%*%beta) 
	Q.opt<- data.frame(QbarAW, Qbar1W, Qbar0W); 			
		
	#----------------------------
	# Step 2 - targeting
	
	#parametric estimation of P(A|W)=g(A|W)
	glm.out<-  suppressWarnings(glm(gform, family='binomial', data=data.frame(W,A), weights=wt))
	g1W<-  predict( glm.out, type='response' )
	g0W<- 1-g1W
	
	# calculate clever covariate -> going to do a 2-dim update here
	H.1W<- A/g1W
	H.0W<- (1-A)/g0W  # note not negative here.
	# it's easier to code a 2-dim update when going after Risk Diff, Risk Ratio, Odds Ratio

	# initial parameter estimates are 0 
	opt.out<- optim(par=c(0,0), fn=tmleLogLike, gr=tmlegrad, Y=Y.tilde, Q.AW=Q.opt$QbarAW, 
		H.0W=H.0W, H.1W=H.1W, wt=wt,  method="BFGS", hessian=F,  control = list(maxit = 500, trace=F)) 
	
	eps<-opt.out$par
	#print(paste('two-dim epsilon', eps, sep=' ')) 
	
	# update. still on the transformed scale.
	QbarAW<-plogis(qlogis(Q.opt$QbarAW)+eps[1]*H.0W + eps[2]*H.1W)
	Qbar0W<-plogis(qlogis(Q.opt$Qbar0W)+eps[1]/g0W)
	Qbar1W<-plogis(qlogis(Q.opt$Qbar1W)+eps[2]/g1W)
	
	Q<- data.frame(cbind(QbarAW, Qbar1W, Qbar0W))
	colnames(Q)<- c("QbarAW", "Qbar1W", "Qbar0W")

	# transform back
	Qstar<- Q*(b-a) +a

	#------------------
	# point estimate
	# for CCDesigns, need to average over case-control weighted dist of cov
	# see van der Laan 2008
	Risk1<- weighted.mean(Qstar$Qbar1W, w=wt)
	Risk0<- weighted.mean(Qstar$Qbar0W, w=wt)
	RiskDiff<- Risk1 - Risk0
	RiskRatio <- Risk1/Risk0
	OddsRatio<- (Risk1/(1-Risk1))/(Risk0/(1-Risk0))

	data.frame(RiskDiff, RiskRatio, OddsRatio)

}


#-------------- 
#logLike -  calculate quasi-loglikelihood loss for a bounded Y with values possibly
#	outside of [0,1]
# this is based on a parametric model for the conditional mean Qbar(A,W)
# input: beta (initial values for beta), outcome Y, design matrix X, weights
# returns: negative quasi-loglikelihood loss
#-------------
logLike<- function(beta, Y, X, wt=1){

	pi<- plogis( X%*%beta ) #Qbar(A,W)= expit(beta0 +beta1*X1 + beta2*X2... )
	pi[pi==0] <- .Machine$double.neg.eps	# prevent taking a log of 0 or 1
	pi[pi==1]<- 1-.Machine$double.neg.eps
	logLike<- sum( wt*(Y*log(pi)  + (1-Y)*log(1-pi))  )
	return(-logLike)		# negative quasi-loglikelihood loss function
}

#-------------------
# grad - corresponding function to calculate the gradient
#	other optimization methods do not require the use of the gradient
#	see ?optim for more help
#--------------------
grad<- function(beta, Y, X, wt=1){
	
	pi<- plogis( X%*%beta)  #Qbar(A,W)= expit(beta0 +beta1*X1 + beta2*X2... )
	pi[pi==0] <- .Machine$double.neg.eps # prevent from taking the log of 0 or 1
	pi[pi==1]<- 1-.Machine$double.neg.eps
	resid<- wt*(Y-pi) # weighted residuals
	gr<- crossprod(X, resid)
	return(-gr)
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

#------------------------------------------
###########################
# SIMULATION 1
##########################
#--------------------------------------------

# Example: Data generating experiment for Simulation 1
set.seed(123)
n=2500
W1<- rnorm(n, 0, .25)
W2<- runif(n, 0, 1)
W3<- rbinom(n, size=1, 0.5)
A<- rbinom(n, size=1, prob= plogis(-.5+ W1+W2+W3) )
pi<- plogis(-3+ 2*A + 1*W1 + 2*W2 - 4*W3 + .5*A*W1)/15
Y<- rbinom(n, size=1, prob= pi)
sum(Y) # total number of outcomes in the data set

#----------------------------------
# Demonstration of the minimizing the quasi log-likelihood 
#---------------------------
# 
# Qbounds (l,u)= (0,0.075)
l=0; u=0.075
#create the design matrix 
X <- model.matrix(as.formula(Y~W1+W2+W3+A*W1))
# transform Y to Y.tilde in between (l,u)
Y.tilde<- (Y - l)/(u-l)
summary(Y.tilde)
 
# call to the optim function. 
# par: initial parameter estimates; f:function to minimize; gr: gradient
# arguments to LogLikelihood() & grad() are Y and X
optim.out <- optim(par=rep(0, ncol(X)), fn= logLike, gr=grad, 
	Y=Y.tilde, X=X, method="BFGS")
# see optim help files for more details and other optimization routines

# get parameter estimates
beta<- optim.out$par
# get predicted values and transform to proper scale
pred.prob.optim <- plogis(X%*%beta)*(u-l) + l
# compare with standard logistic regression
pred.prob.glm<- predict(glm(Y~W1+W2+W3+A*W1, family="binomial"), type="response")
predictions<- data.frame(optim=pred.prob.optim,  glm=pred.prob.glm)
summary(predictions*100)

#--------------
# running the full rTMLE algorithm
#-------------
W<- data.frame(W1,W2,W3)
Qform<-as.formula(Y~W1+W2+W3+A*W1) 
gform<- as.formula(A~ W1+W2+W3)
wt<- rep(1, n)
rtmle.out<- rtmle.wt(W=W, A=A, Y=Y, 
	a=l, b=u, Qform=Qform, gform=gform,  wt=wt) 
 
# Please contact the authors for code to get inference based on the 
# estimated influence curve.

#-----------------------------------------
###################################
# SIMULATION 2 - CASE-CONTROL DATA
##################################
#-----------------------------------------

#-------------
# generate data for the underlying cohort
#-------------
set.seed(123)
nCohort<- 500000

# generate three baseline covariates
W1<- rnorm(nCohort, 0, .25)
W2<- runif(nCohort, 0, 1)
W3<- rbinom(nCohort, 1, .5)

# generate the exposure (binary) 
A<- rbinom(nCohort, size=1, prob= plogis(-.5+W1+W2+W3))

# generate the outcome (binary)
Y<- rbinom(nCohort, size=1, prob= plogis(-3+2*A+W1+2*W2-4*W3+.5*A*W1)/15 )

# true marginal probability of the outcome
q<- mean(Y)
q

Cohort<- data.frame(W1, W2, W3, A, Y)
rm(W1, W2, W3, A, Y)

########################
# Draw the case-control sample
########################
# specify number of cases and number of controls
nCa<- 500
nCo<- 650
n<- nCa + nCo
J<- nCo/nCa # ratio of number controls to cases

# draw the case-control sample
cases<- sample( which(Cohort$Y==1), size=nCa)
controls<- sample( which(Cohort$Y==0), size=nCo)

data<- Cohort[c(cases, controls),]
sum(data$Y==1) #should have 500 cases

# specify the weights
wt<- rep(NA, n)
wt[data$Y==1]<- q # for cases
wt[data$Y==0]<- (1-q)/J# for controls

# specify regression models for the conditional mean outcome E(Y|A,W)
#		and the exposure mechanism g(A|W)
Qform<- as.formula(Y~A+W1+W2+W3) 
gform<- as.formula(A~W1+W2+W3)

#####
MLE.est<- mle.wt(data=data, Qform=Qform, wt=wt)
TMLE.est<- tmle.wt(data=data, Qform=Qform, gform=gform, wt=wt)
# rTMLE with bounds on P(Y=1|A,W) as a=0 and b=7*q
rTMLE.est<- rtmle.wt(W=data[,c('W1','W2','W3')], A=data$A, Y=data$Y, 
	a=0, b=7*q, Qform=Qform, gform=gform, wt=wt)

#---------------------
MLE.est
TMLE.est
rTMLE.est

#---------------
# Bootstrap based on 100 samples
#	setting the number of bootstrapped samples low for illustration
#-----------
numReps<- 100

# matrices to hold estimates
RD<- data.frame(matrix(NA, nrow=numReps, ncol=3))
colnames(RD)<- c('MLE', 'TMLE', 'rTMLE')
RR<- OR<- RD

for(i in 1:numReps){
	
	cases<- sample(which(data$Y==1), size=nCa, replace=T)
	controls<- sample(which(data$Y==0), size=nCo, replace=T )
	
	data.b<- data[c(cases, controls), ]
	wt.b<- wt[c(cases,controls)]	
		
	MLE.b<- mle.wt(data=data.b, Qform=Qform, wt=wt.b)
	TMLE.b<- tmle.wt(data=data.b, Qform=Qform, gform=gform, wt=wt.b)
	rTMLE.b<- rtmle.wt(W=data.b[,c('W1','W2','W3')], A=data.b$A, Y=data.b$Y, 
		a=0, b=7*q, Qform=Qform, gform=gform, wt=wt.b)

	RD[i, ]<- c(MLE.b$RiskDiff, TMLE.b$RiskDiff, rTMLE.b$RiskDiff)
	RR[i, ]<- c(MLE.b$RiskRatio, TMLE.b$RiskRatio, rTMLE.b$RiskRatio)
	OR[i, ]<- c(MLE.b$OddsRatio, TMLE.b$OddsRatio, rTMLE.b$OddsRatio)
	
}

colMeans(RD)
sqrt(apply(RD,2, var))
# After suffficient replications, can then use the variance of the bootstrap estimates
#	as an estimator of the variance... or directly get the 2.5% and 97.5% quantiles
