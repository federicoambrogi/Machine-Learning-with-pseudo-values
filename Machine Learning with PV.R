## .-------------------------------------------------------------------------------------------.
## | This file contains the R-code underlying the analysis in                                  |
## | Ambrogi & Andersen (2016), Predicting survival probabilities using pseudo-observations    |
## | Federico Ambrogi <federico.ambrogi@unimi.it>                                              |
## `-------------------------------------------------------------------------------------------'

rm(list=ls())
library(prodlim)
library(pec) 
library(mgcv)
library(caret)
library(nnet)
library(survival)
library(party)

# function to compute quantiles
quantili = function(x){
  c(1:x)/(x+1)
}

# ---------------------------------------------------------------------------------- #
# function to prepare data with pseudo observations									 #
# time: follow-up time (right censored).                                             #
# status: the status indicator, normally 1:death; 0: alive.                          #
# data: a data.frame in with the variables time, status and the covariates.          #
# prob: quantiles of the time distribution at which to compute the pseudo-values,    #
#       also called pseudo-times.                                                    #
# ---------------------------------------------------------------------------------- #

preparePV <- function(time, status, data, prob=c(.04, seq(.15,.85,.1), .96)){
  eventi <- sort(unique(time[status==1]))
  #    cutoffs = quantile(eventi, probs = prob)
  cutoffs = quantile(time, probs = prob)
  NPSEUDO = length(cutoffs)
  pseudo = pseudo::pseudosurv(time, status, tmax=cutoffs )
  pka <- NULL
  for(it in 1:length(pseudo$time)){
    pka <- rbind(pka,cbind(data, pseudo = pseudo$pseudo[,it],
                           tpseudo = pseudo$time[it], id=1:nrow(data)))
  }
  data <- pka[order(pka$id),]
}

# ---------------------------------------------------------------------------------- #
# function to fit a GAM to survival data using pseudo-observations.					 #
# in this version a pseudo-time is placed each 10 event times.
# formula: a formula object, with the response (pseudo) on the left of a ~ operator, #
#          and the terms on the right. The first term is tpseudo.                    #
# data: a data.frame with the variables time, status and the covariates listed in    #
#       the formula. 																 #
# index: a two element vector. The first element gives the position of the column    #
#       with time, the second element gives the position of the column with status   #
# gamma: see gam function															 #
# ---------------------------------------------------------------------------------- #

gamPV <- function(formula, data, index, gamma=1,...){
  time <- data[,index[1]]
  status <- data[,index[2]]
  
  var.n <- all.vars(formula)
  X <- data[,var.n[-index]]
  
  conta <- table(status)
  numero <- round(conta[2]/10)
  data <- preparePV(time, status, X, prob = quantili(numero))
  #    data <- preparePV(time, status, X)
  
  fitted = gam(formula, data=data, 
               family=quasi(link="cloglog"),  scale=1, etastart= rep(.1, nrow(data)), 
               control=gam.control(keepData=TRUE), gamma=gamma)
  fitted    
}

# ---------------------------------------------------------------------------------- #
# function to get survival probabilities from models estimated with gamPV.			 #
# For details see:																	 #
# Mogensen UB, Ishwaran H, Gerds TA. Evaluating random forests for survival analysis #
# using prediction error curves. Journal of Statistical Software 2012; 50(11):1–23.  #
# ---------------------------------------------------------------------------------- #
 
predictSurvProb.gam <- function(object,newdata,times,...){
  dfnew <- data.frame(newdata[rep(1:nrow(newdata), each=length(times)),], tpseudo=rep(times, nrow(newdata)) )    
  gampred.object <- predict.gam(object, newdata = dfnew, type = "response", outer.ok = TRUE, se.fit=FALSE)    
  p <- t(matrix(gampred.object, length(times), nrow(newdata)))
}

# ---------------------------------------------------------------------------------- #
# function to fit a neural network to survival data using pseudo-observations		 #
# formula: a formula object, with the response (pseudo) on the left of a ~ operator, #
#          and the terms on the right. The first term is tpseudo.                    #
# data: a data.frame with the variables time, status and the covariates listed in    #
#       the formula. 																 #
# index: a two element vector. The first element gives the position of the column    #
#       with time, the second element gives the position of the column with status   #
# ---------------------------------------------------------------------------------- #

nnetPV <- function(formula,data, index,...){
    time <- data[,index[1]]
    status <- data[,index[2]]
    
    var.n <- all.vars(formula)
    X <- data[,var.n[-index]]
    
    conta <- table(status)
    numero <- round(conta[2]/10)
    data2 <- preparePV(time, status, X, prob = quantili(numero))
    #    data <- preparePV(time, status, X)
    #    data2 <- data2[1:200,]
    control <- trainControl(method="repeatedcv", number=10, repeats=10)
    my.grid <- expand.grid(.decay = seq(.01, .001, by=-.003), .size = c(4,6))
    model <- train(formula, data=data2, 
                   method = "nnet", tuneGrid=my.grid, trControl=control, maxit=1000, 
                   preProcess=c("center", "scale"))
    model
}

# ---------------------------------------------------------------------------------- #
# function to get survival probabilities from models estimated with nnetPV.			 #
# For details see:																	 #
# Mogensen UB, Ishwaran H, Gerds TA. Evaluating random forests for survival analysis #
# using prediction error curves. Journal of Statistical Software 2012; 50(11):1–23.  #
# ---------------------------------------------------------------------------------- #

predictSurvProb.train <- function(object,newdata,times,...){
    dfnew <- data.frame(newdata[rep(1:nrow(newdata), each=length(times)),], 
    tpseudo=rep(times, nrow(newdata)) )    
    nnetpred.object <- predict(object, newdata = dfnew)    
    p <- t(matrix(nnetpred.object, length(times), nrow(newdata)))
}

# ---------------------------------------------------------------------------------- #
# function to fit a conditional random forest to survival data using 				 #
# pseudo-observations																 #
# formula: a formula object, with the response (pseudo) on the left of a ~ operator, #
#          and the terms on the right. The first term is tpseudo.                    #
# data: a data.frame with the variables time, status and the covariates listed in    #
#       the formula. 																 #
# index: a two element vector. The first element gives the position of the column    #
#       with time, the second element gives the position of the column with status   #
# ---------------------------------------------------------------------------------- #

cfPV <- function(formula,data, index,...){
  time <- data[,index[1]]
  status <- data[,index[2]]
  
  var.n <- all.vars(formula)
  X <- data[,var.n[-index]]
  
  conta <- table(status)
  numero <- round(conta[2]/10)
  data2 <- preparePV(time, status, X, prob = quantili(numero))
  #    data <- preparePV(time, status, X)
  #    data2 <- data2[1:200,]
  model <- list(forest = cforest(formula, data=data2, 
  controls=cforest_classical(ntree=500)))
  class(model) <- "cfPV"
  model$call <- match.call()
  model
}

# ---------------------------------------------------------------------------------- #
# function to get survival probabilities from models estimated with cfPV.			 #
# For details see:																	 #
# Mogensen UB, Ishwaran H, Gerds TA. Evaluating random forests for survival analysis #
# using prediction error curves. Journal of Statistical Software 2012; 50(11):1–23.  #
# ---------------------------------------------------------------------------------- #

predictSurvProb.cfPV <- function(object,newdata,times,...){
  times <- as.numeric(times)
  dfnew <- data.frame(newdata[rep(1:nrow(newdata), each=length(times)),], 
  tpseudo=rep(times, nrow(newdata)) )    
  cfpred.object <- predict(object$forest, newdata = dfnew)    
  cfpred <- pmax(0, pmin(1, cfpred.object))  
  p <- t(matrix(cfpred.object, length(times), nrow(newdata)))
}

# loading the data from the GBSG2 study. Available in package pec
data(GBSG2)
GBSG2$status <- GBSG2$cens
GBSG2 <- GBSG2[order(GBSG2$time,-GBSG2$status),]
GBSG2$grade.bin <- as.factor(as.numeric(GBSG2$tgrade!="I"))
levels(GBSG2$grade.bin) <- c("I","II/III")
GBSG2 <- GBSG2[, c(9,11, 1:8, 12)]

# ---------------------------------------------------------------------------------- #
# Covariates in GBSG2 data															 #
# AGE: age																			 #
# Tumor size: tsize																	 #
# No. positive lymph nodes: pnodes													 #
# Progesterone receptor: progrec													 #
# Estrogen receptor: estrec															 #
# Tumor grade: tgrade																 #
# ---------------------------------------------------------------------------------- #

#GBSG2$tn <- exp(-.12*GBSG2$pnodes)

# ---------------------------------------------------------
# lets get a random sample of observations 
# from this dataset to speed the computational time
# ---------------------------------------------------------
set.seed(200)
GBSG2.s <- GBSG2[sample(nrow(GBSG2), 100), ]
set.seed(2006)

# ---------------------------------------------------------------------------------- #
# cforest: estimate a conditional forest model for survival data                     #
# (see Mogensen UB, Ishwaran H, Gerds TA. 2012); 									 #
# Cox: multiple cox regression model.  												 #
# gam: GAM estimated with pseudo-values without interactions covariates by time. 	 #
#      The time effect is estimated with a P-spline on pseudo times.				 #
# gam2: GAM estimated with pseudo-values without interactions covariates by time. 	 #
#      The time effect is estimated with a P-spline on pseudo times. 				 #
#      Covariate effects are estimated using thin plane splines.					 #
# gam3: GAM estimated with pseudo-values with interactions covariates by time. 		 #
#      tensor product splines are used to model interactions.						 #
# ann: artificial neural network estimated with pseudo-values. 						 #
# cf: conditional forest model estimated with pseudo-values.						 #
# ---------------------------------------------------------------------------------- #

gbsg2.models <- list("cforest"=pecCforest(Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,
                      data=GBSG2.s,controls=cforest_classical(ntree=500)),
                     "Cox"=coxph(Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2.s),
                     "gam"  = gamPV(pseudo ~ s(tpseudo, bs="ps", m=c(2,2), k=4) + age+tsize+grade.bin+pnodes+progrec+estrec, 
                      data=GBSG2.s, index=c(1,2)),
                     "gam2" = gamPV(pseudo ~ s(tpseudo, bs="ps", m=c(2,2), k=5) + s(age, k=5) + s(tsize, k=5) + 
                      grade.bin + s(pnodes, k=5) + s(progrec, k=5) + s(estrec, k=5), data=GBSG2.s, index=c(1,2)),
                     "gam3" = gamPV(pseudo ~ s(tpseudo, bs="ps", m=c(2,2), k=5, by=grade.bin) + 
                      te(tpseudo, age, bs="ps", m=c(2,2), k=5) + 
                      te(tpseudo, tsize, bs="ps", m=c(2,2), k=5) + 
                      te(tpseudo, pnodes, bs="ps", m=c(2,2), k=5) + 
                      te(tpseudo, progrec, bs="ps", m=c(2,2), k=5) + 
                      te(tpseudo, estrec, bs="ps", m=c(2,2), k=5), data=GBSG2.s, index=c(1,2)),
                     "ann"=nnetPV(pseudo ~ tpseudo+age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2.s, index=c(1,2)),
                     "cf" = cfPV(pseudo ~ tpseudo+age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2.s, index=c(1,2))
)

# ---------------------------------------------------------------------------------- #
# here we obtain the prediction for a speficic covariate pattern from the different	 #
# methods in study - this is for illustrative purposes								 #
# ---------------------------------------------------------------------------------- #

newdata <- data.frame(age = as.integer(round(mean(GBSG2.s$age))),
                      tsize = as.integer(round(mean(GBSG2.s$tsize))),
                      grade.bin = GBSG2.s$grade.bin[1],
                      pnodes = as.integer(round(mean(GBSG2.s$pnodes))),
                      progrec = as.integer(round(mean(GBSG2.s$progrec))),
                      estrec = as.integer(round(mean(GBSG2.s$estrec))))

ann.pred <- predictSurvProb.train(gbsg2.models$ann,newdata,times=1:2000)
cf.pred <- predictSurvProb.cfPV(gbsg2.models$cf,newdata,times=1:2000)
cforest.pred <- predictSurvProb(gbsg2.models$cforest,newdata,times=1:2000)
gam3.pred <- predictSurvProb.gam(gbsg2.models$gam3,newdata,times=1:2000)
gam2.pred <- predictSurvProb.gam(gbsg2.models$gam2,newdata,times=1:2000)

plot(1:2000, cforest.pred, type="s", ylim=c(0,1), lty=1)
lines(1:2000, cf.pred, lty=2)
lines(1:2000, gam3.pred, lty=3)
lines(1:2000, gam2.pred, lty=4)
lines(1:2000, ann.pred, lty=5)
legend(0, .5, legend=c("cforest", "cf PV", "GAM TD", "GAM noTD", "annPV"), lty=1:5)

# --------------------------------------------------------------------- #
# Computing the apparent error and No Information Error					#
# Function pec computes the No Information Error only when 				#
# splitMethod ="Boot632plus"											#
# We therefore use B=1, just for convenience but the error estimates 	#
# are not used															#
# --------------------------------------------------------------------- #

GBSG2.noinf <- pec(object=gbsg2.models,
                      formula=Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,
                      data=GBSG2.s,
                      cens.model="cox",
                      splitMethod ="Boot632plus",
                      verbose = TRUE,
                      B=1,
                      maxtime=2000,
                      exact=FALSE)

# --------------------------------------------------------------------- #
# The following code was extrapolated from the function "pec"			#
# to perform the resampling for the models under study and compute		#
# the prediction error on the left out samples.							#
# The function pec uses an update of the model object 					#
# and does not work properly with models estimated with PV				#
# --------------------------------------------------------------------- #

# ------------------------ WARNING ------------------------------------ #
# This piece of code takes some time to run								#
# use load("GBSG_sample_res.RData") with all computations already done	#
# --------------------------------------------------------------------- #

set.seed(2006)
N=nrow(GBSG2.s)
B=100
ResampleIndex <- do.call("cbind", lapply(1:B, function(b) {
    sort(sample(1:N, size = N, replace = TRUE))
}))
colnames(ResampleIndex) <- paste("Train", 1:B, sep = ".")

ptm <- proc.time()
boot.r <- do.call(cbind, lapply(1:B, function(x){
    # this bootstrap version selects the patients and recalculates the PV
    print(x)
    train.b <- GBSG2.s[ResampleIndex[,x],]
    vindex.b.2 <- (1:N)[!(1:N) %in% unique(ResampleIndex[,x])]
    val.b <- GBSG2.s[vindex.b.2,,drop=FALSE]
    
    gbsg2.models.b <- try(list("cforest"=pecCforest(Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,data=train.b,controls=cforest_classical(ntree=1000)),
                               "Cox"=coxph(Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,data=train.b,x=TRUE,y=TRUE),
                               "gam"  = gamPV(pseudo ~ s(tpseudo, bs="ps", m=c(2,2), k=5) + s(age, k=3) + s(tsize, k=3) + 
                                grade.bin + s(pnodes, k=3) + s(progrec, k=3) + s(estrec, k=3), data=train.b, index=c(1,2)),
                               "gam2" = gamPV(pseudo ~ s(tpseudo, bs="ps", m=c(2,2), k=5) + s(age, k=5) + s(tsize, k=5) + 
                                grade.bin + s(pnodes, k=5) + s(progrec, k=5) + s(estrec, k=5), data=train.b, index=c(1,2)),
                               "gam3" = gamPV(pseudo ~ s(tpseudo, bs="ps", m=c(2,2), k=5) + 
                               te(tpseudo, age, bs="ps", m=c(2,2), k=5) + 
                               te(tpseudo, tsize, bs="ps", m=c(2,2), k=5) + 
                               grade.bin + 
                               te(tpseudo, pnodes, bs="ps", m=c(2,2), k=5) + 
                               te(tpseudo, progrec, bs="ps", m=c(2,2), k=5) + 
                               te(tpseudo, estrec, bs="ps", m=c(2,2), k=5), data=train.b, index=c(1,2)),
                               "ann"=nnetPV(pseudo ~ tpseudo+age+tsize+grade.bin+pnodes+progrec+estrec,data=train.b, index=c(1,2)),
                               "cf" = cfPV(pseudo ~ tpseudo+age+tsize+grade.bin+pnodes+progrec+estrec,data=train.b, index=c(1,2))

    ))
    if (inherits(gbsg2.models.b, 'try-error')) {boot632 <- NULL}
    else{
        boot632 <- pec(object=gbsg2.models.b,
                       formula=Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,
                       data=val.b,
                       cens.model="cox",
                       splitMethod ="none",
                       maxtime=2000,
                       exact=FALSE)
        if (inherits(boot632, 'try-error')) {boot632 <- NULL}
    }
    return(boot632$AppErr)
    
}
))

# save(AP, boot.r, file="GBSG_sample_res.RData")
load("GBSG_sample_res.RData")
res <- boot.r
NSIM=100

# --------------------------------------------------------------------- #
# The following code was taken from the function "pec"	        		#
# to perform the .632 plus prediction error calculations				#
# --------------------------------------------------------------------- #

BE.ref  <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$Reference))
BCVE <- apply(BE.ref, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$Reference)
overfit <- (Err1 - AP$AppErr$Reference) / (AP$NoInfErr$Reference - AP$AppErr$Reference)
overfit[!(Err1>AP$AppErr$Reference)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.Reference <- (1-w) * AP$AppErr$Reference  + w * Err1
B632Err.Reference <- 0.368 * AP$AppErr$Reference + 0.632 * BCVE

BE.cforest  <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$cforest))
BCVE <- apply(BE.cforest, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$cforest)
overfit <- (Err1 - AP$AppErr$cforest) / (AP$NoInfErr$cforest - AP$AppErr$cforest)
overfit[!(Err1>AP$AppErr$cforest)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.cforest <- (1-w) * AP$AppErr$cforest  + w * Err1
B632Err.cforest <- 0.368 * AP$AppErr$cforest + 0.632 * BCVE

BE.cf  <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$cf))
BCVE <- apply(BE.cf, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$cf)
overfit <- (Err1 - AP$AppErr$cf) / (AP$NoInfErr$cf - AP$AppErr$cf)
overfit[!(Err1>AP$AppErr$cf)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.cf <- (1-w) * AP$AppErr$cf  + w * Err1
B632Err.cf <- 0.368 * AP$AppErr$cf + 0.632 * BCVE

BE.Cox  <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$Cox))
BCVE <- apply(BE.Cox, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$Cox)
overfit <- (Err1 - AP$AppErr$Cox) / (AP$NoInfErr$Cox - AP$AppErr$Cox)
overfit[!(Err1>AP$AppErr$Cox)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.Cox <- (1-w) * AP$AppErr$Cox  + w * Err1
B632Err.Cox <- 0.368 * AP$AppErr$Cox + 0.632 * BCVE

BE.gam <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$gam))
BCVE <- apply(BE.gam, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$gam)
overfit <- (Err1 - AP$AppErr$gam) / (AP$NoInfErr$gam - AP$AppErr$gam)
overfit[!(Err1>AP$AppErr$gam)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.gam <- (1-w) * AP$AppErr$gam  + w * Err1
B632Err.gam <- 0.368 * AP$AppErr$gam + 0.632 * BCVE

BE.gam2 <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$gam2))
BCVE <- apply(BE.gam2, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$gam2)
overfit <- (Err1 - AP$AppErr$gam2) / (AP$NoInfErr$gam2 - AP$AppErr$gam2)
overfit[!(Err1>AP$AppErr$gam2)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.gam2 <- (1-w) * AP$AppErr$gam2  + w * Err1
B632Err.gam2 <- 0.368 * AP$AppErr$gam2 + 0.632 * BCVE

BE.gam3 <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$gam3))
BCVE <- apply(BE.gam3, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$gam3)
overfit <- (Err1 - AP$AppErr$gam3) / (AP$NoInfErr$gam3 - AP$AppErr$gam3)
overfit[!(Err1>AP$AppErr$gam3)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.gam3 <- (1-w) * AP$AppErr$gam3  + w * Err1
B632Err.gam3 <- 0.368 * AP$AppErr$gam3 + 0.632 * BCVE

BE.ann <- do.call(cbind, lapply(1:NSIM, function(x) res[,x]$ann))
BCVE <- apply(BE.ann, 1, mean, na.rm=TRUE)
Err1 <- pmin(BCVE,AP$NoInfErr$ann)
overfit <- (Err1 - AP$AppErr$ann) / (AP$NoInfErr$ann - AP$AppErr$ann)
overfit[!(Err1>AP$AppErr$ann)] <- 0
w <- .632 / (1 - .368 * overfit)
B632plusErr.ann <- (1-w) * AP$AppErr$ann  + w * Err1
B632Err.ann <- 0.368 * AP$AppErr$ann + 0.632 * BCVE

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

plot(AP$time, AP$AppErr$Reference, type="l", lty=1, col="gray", 
     xlab="Time", ylab="Apparent Error", lwd=1, ylim=c(0, .3))
lines(AP$time, AP$AppErr$Cox, lty=2, lwd=1)
#lines(AP$time, AP$AppErr$gam, lty=3, lwd=2)
lines(AP$time, AP$AppErr$gam2, lty=3, lwd=1)
lines(AP$time, AP$AppErr$gam3, lty=4, lwd=1)
lines(AP$time, AP$AppErr$cforest, lty=5, lwd=1)
lines(AP$time, AP$AppErr$cf, lty=6, lwd=2)
lines(AP$time, AP$AppErr$ann, lty=7, lwd=2)
legend(0,.3, c("Reference", "Cox", "GAM no TD", "GAM TD", "cforest", "cforest PV", "nnet PN"), lty=c(1:7),
       bty="n", col=c("gray", rep("black", 6)), lwd=2, cex=.5)
grid()

plot(AP$time, AP$NoInfErr$Reference, type="l", lty=1, col="gray",
     xlab="Time", ylab="No Inf Error", lwd=1, ylim=c(0, .4))
lines(AP$time, AP$NoInfErr$Cox, lty=2, lwd=1)
#lines(AP$time, AP$NoInfErr$gam, lty=3, lwd=2)
lines(AP$time, AP$NoInfErr$gam2, lty=3, lwd=1)
lines(AP$time, AP$NoInfErr$gam3, lty=4, lwd=1)
lines(AP$time, AP$NoInfErr$cforest, lty=5, lwd=1)
lines(AP$time, AP$NoInfErr$cf, lty=6, lwd=2)
lines(AP$time, AP$NoInfErr$ann, lty=7, lwd=2)
grid()

plot(AP$time, B632plusErr.Reference, type="l", lty=1, col="gray", xlab="Time", ylab=".632+ Bootstrap Error", lwd=1, ylim=c(0, .4))
lines(AP$time, B632plusErr.Cox, lty=2, lwd=1)
lines(AP$time, B632plusErr.gam2, lty=3, lwd=1)
lines(AP$time, B632plusErr.gam3, lty=4, lwd=1)
lines(AP$time, B632plusErr.cforest, lty=5, lwd=1)
lines(AP$time, B632plusErr.cf, lty=6, lwd=2)
lines(AP$time, B632plusErr.ann, lty=7, lwd=2)
legend(0,.25, c("Reference", "Cox", "GAM no TD", "GAM TD", "cforest", "cforest PV", "nnet PV"), lty=c(1:7),
       bty="n", col=c("gray", rep("black", 6)), lwd=c(1,1,1,1,1,2,2), cex=.5)
grid()

plot(AP$time, B632Err.Reference, type="l", lty=1, col="gray", xlab="Time", ylab=".632 Bootstrap Error", lwd=1, ylim=c(0, .4))
lines(AP$time, B632Err.Cox, lty=2, lwd=1)
lines(AP$time, B632Err.gam2, lty=3, lwd=1)
lines(AP$time, B632Err.gam3, lty=4, lwd=1)
lines(AP$time, B632Err.cforest, lty=5, lwd=1)
lines(AP$time, B632Err.cf, lty=6, lwd=2)
lines(AP$time, B632Err.ann, lty=7, lwd=2)
#legend(0,.25, c("Reference", "Cox", "cforest", "GAM no TD", "GAM TD", "ann"), lty=c(1:6),
#       bty="n", col=c("gray", rep("black", 5)), lwd=2)
grid()