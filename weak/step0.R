### R code from vignette source '/localhome/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step0'
### Encoding: UTF-8

###################################################
### code chunk number 1: Stangle (eval = FALSE)
###################################################
## #setwd("/home/ann-sophie/Documents/wip/lasso/interactions/Weak/")
## setwd("/run/user/551438008/gvfs/smb-share:domain=SUND,server=sund.root.ku.dk,share=users,user=zpt949/LINUX/ifsv/zpt949/Documents/wip/lasso/interactions/Weak")
## #setwd("/localhome/ifsv/zpt949/Documents/wip/lasso/interactions/Weak")
## Stangle("/localhome/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step0")
## #Stangle("/run/user/551438008/gvfs/smb-share:domain=SUND,server=sund.root.ku.dk,share=users,user=zpt949/LINUX/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step0")
## #Stangle("/home/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step0")
## #Stangle("/home/ann-sophie/Documents/wip/lasso/interactions/Weak/step0")
## #Stangle("step0")


###################################################
### code chunk number 2: Preamble
###################################################
library(glmnet)
library(Matrix)
library(grpreg)
library(MASS)
library(xtable)
library(plyr)
library(parallel)

colo <- "black"

red <- green <- blue <- 152
R <- c(0,red,42,153,156,255,144,54,144,70,225,144,120)
G <- c(0,green,103,150,85,95,26,90,180,95,166,125,26,120)
B <- c(0,blue,236,50,141,24,30,165,255,24,60,154,30,120)

cbPalette <- rgb(R,G,B, max = 255)
cbPaletteT <- rgb(R,G,B, alpha=51, maxColorValue=255)


a4plotWidth <- 8.3-0.9-1.8

clean <- function(string, nums = 1:10, alph = LETTERS, others = NULL) {
  st1 <- string[string != "(Intercept)"]
  for(l in others){
    st1 <- sub(l, "", st1)
  }
  st <- substr(st1,nchar(st1),nchar(st1))
  for(l in c(nums, alph)){
    st <- sub(l, "", st)
    ST <- paste0(substr(st1,1,nchar(st1)-1),st)
  }
  return(sort(unique(ST)))
}

# Screening: The global strong rule
gsr <- Vectorize(function(x, y, lambda=0) {
  res <- abs(as.vector(crossprod(x, y)))
  lambdamax <- max(abs(res))
  discard <- res<(2*lambda-lambdamax)
  seq.int(1, length(res))[!discard]
}, "lambda")
# Screening: The sequential DPP
sdpp <- Vectorize(function(x, y, lambda=0, lambdaprim=n) {
  fit <- glmnet(x, y, lambda=lambdaprim) 
  B <- coef(fit)[-1,]
  res <- abs(as.vector(crossprod(x, y-x%*%B)))
  discard <- res<2*lambda-lambdaprim
  R <- seq.int(1, length(res))[!discard]
  r <- length(R)
  list(R=R, r=r)
}, "lambda")


###################################################
### code chunk number 3: Parameters
###################################################
p <- 600
n <- 200
# Working variables
nz <- arrange(expand.grid(a=1:p,b=1:p),a)
# Features
set.seed(1)
X <- matrix(sample(0:2, n*p, replace=TRUE), nrow=n, ncol=p)
colnames(X) <- paste0("X[, ", 1:p, "]")
# Interactions
ind <- t(combn(ncol(X),2))
out <- apply(ind, 1, function(x) X[,x[1]] * X[,x[2]])
colnames(out) <- apply(ind, 1, function(x) paste0("X[, ",x[1],"]:X[, ",x[2],"]"))
# Design
design <- cbind(X, out)



###################################################
### code chunk number 4: function
###################################################

step1 <- function(truth, strat, xaxis, scree, X, design, i,
                     W=NULL, N, p, n, bet1, bet2=NULL, bet3=NULL, bet12=NULL, rh=NULL, mu=NULL){
  a <- Sys.time()
  
  # Parameters
  MEX <- ifelse(strat=="none", 1, length(eval(parse(text=strat))))
  K <- length(eval(parse(text=xaxis)))
  M <- choose(p, 2)
  beta <- rep(0, p)
  beta12 <- rep(0,M)
  
  # Call lapply
  IFDR <- lapply(1:MEX, function(EX){
    # Stratification settings
    if(strat == "bet2"){
      beta[2] <- bet2[EX]
      beta12[1] <- bet12
    } else if(strat == "bet12"){
      beta12[1] <- bet12[EX]
    } else if(strat == "bet3"){
      beta[3] <- bet3[EX]
    } else if(strat == "rh"){
      if(xaxis == "bet1"){
        beta12[1] <- bet12
      } else if(xaxis == "bet12"){
        beta[1] <- bet1
        beta[2] <- bet1
      }
      rho <- rh[EX] 
      Sigma <- matrix(rho, nrow=p, ncol=p) + diag(p)*(1-rho)
      # Features / design
      set.seed(1)
      X <- mvrnorm(n=n, mu=mu, Sigma=Sigma)
      colnames(X) <- paste0("X[, ", 1:p, "]")
      ind <- t(combn(ncol(X),2))
      out <- apply(ind, 1, function(x) X[,x[1]] * X[,x[2]])
      colnames(out) <- apply(ind, 1, function(x) paste0("X[, ",x[1],"]:X[, ",x[2],"]"))
      design <- cbind(X, out)
    } 
    ifdr <- c(); isem <- c(); ip <- c(); iSelectB1 <- list(); 
    RESAM <- list()
    coefs <- list()
    # Loop over x-axis values
    for(j in 1:K){
      if(truth == "Strong"){
        if(xaxis == "bet1") {
          beta[1] <- bet1[j]
          if(strat == "bet12") beta[2] <- bet1[j]
        } else if(xaxis == "bet12"){
          beta12[1] <- bet12[j]
        } else if(xaxis == "bet3"){
          beta[3] <- bet3[j]
        }
      } else if(truth == "Weak"){
        beta[1] <- bet1[j]
      } else if(truth == "Anti"){
        beta[3] <- bet1[j]
      } else if(truth == "All"){
        beta[1] <- bet1[j]
        beta12[1] <- bet12[j]
        beta12[2] <- bet12[j]
        beta12[1198] <- bet12[j]
      } else if(truth=="PI"){
        beta12[1] <- bet12[j]
      } else if(truth=="PM"){
        beta[1] <- bet1[j]
      }
      
        # Gaussian noise
        set.seed(20000+i)
        eps <- rnorm(n)
        # Outcomes
        y <- design %*% c(beta,beta12) + eps
        # Lasso regression
        if(W==1){
          set.seed(2000+i) #
          fit <- try(cv.glmnet(X, y, family="gaussian", standardize = TRUE), silent=TRUE)
          if(inherits(fit, "try-error")) next
          # Zero and non-zero coefficients
          c1 <- coef(fit, s = "lambda.1se") #
        } else {
          MSEs <- matrix(ncol=W, nrow=100)
          LAMBDAs <- matrix(ncol=W, nrow=100)
          for(l in 1:W){
            set.seed(i*2000+l) 
            fit <- cv.glmnet(X, y, family="gaussian", standardize = TRUE)
            MSEs[1:length(fit$cvm),l] <- fit$cvm
            LAMBDAs[1:length(fit$lambda),l] <- fit$lambda
          }
          lambda.min <- rowMeans(LAMBDAs, na.rm=TRUE)[which.min(rowMeans(MSEs, na.rm=TRUE))]
          # Zero and non-zero coefficients
          c1 <- coef(fit, s = lambda.min) 
        }
        # Non-zero coefficients 
        nZ1 <- names(as.matrix(c1)[as.matrix(c1) != 0,])
        nZ1 <- nZ1[nZ1 != "(Intercept)"]
        # Selected main effects
        iSelectA1 <- colnames(X)[colnames(X) %in% clean(nZ1)]
        # Keep track
        cat(".")
        # Return
        RESAM[[j]] <- list(iSelectA1=iSelectA1, y=y)
      # Resampling done
    }
    # Return
    return(list(RESAM=RESAM))
  })
  
  # Running time
  tim <- difftime(Sys.time(),a,units="mins")
  cat(tim, "mins")
  #Return
  list(IFDR=IFDR, ETIM=tim, n=n, p=p, N=N, W=W)
  
}


###################################################
### code chunk number 5: step0.Rnw:210-390
###################################################
step2 <- function(hierarchy, truth, strat, xaxis, scree, X, design, i,
                     W=NULL, N, p, n, bet1, bet2=NULL, bet3=NULL, bet12=NULL, rh=NULL, mu=NULL){
  a <- Sys.time()
  
  # Parameters
  MEX <- ifelse(strat=="none", 1, length(eval(parse(text=strat))))
  K <- length(eval(parse(text=xaxis)))
  M <- choose(p, 2)
  beta <- rep(0, p)
  beta12 <- rep(0,M)
  # Working variables
  nz <- arrange(expand.grid(a=1:p,b=1:p),a)
  
  # Call lapply
  IFDR <- lapply(1:MEX, function(EX){
    ifdr <- c(); Vi <- c(); Si <- c(); 
    # Loop over x-axis values
    for(j in 1:K){
         
        ###################################################
        # Feature selection step 2: Interaction selection #
        ###################################################
        iS1 <- STEP1$IFDR[[EX]]$RESAM[[j]]$iSelectA1
        y <- STEP1$IFDR[[EX]]$RESAM[[j]]$y
        if(hierarchy == "Weak"){
          # Features: Include only interactions between at least 
          # one selected main effect and another main effect and 
          # all main effects with no penalty on selected main effects. 
          if(length(iS1) == 0){
            t4 <- 1:p
          }else{
            t0 <- gsub("\\[", "\\\\[", iS1)
            t1 <- gsub("\\]", "\\\\]", t0)
            t2 <- list()
            for(k in 1:length(t1)){
              t2[[k]] <- grep(t1[k], colnames(design), value=FALSE)  
            }
            t3 <- unlist(t2)
            t4 <- unique(c(1:p,t3))
          }
          XX <- design[,t4]
        } else if(hierarchy == "Strong"){
          # Features: Include only interactions between selected main effects 
          # and all main effects with no penalty on selected main effects. 
          if(length(iS1) %in% 0:1){
            t4 <- 1:p
            XX <- design[,t4]
          }else{
            t0 <- arrange(expand.grid(a=iS1,b=iS1),a)
            t1 <- t0[t0[1]!=t0[2],]  
            t2 <- apply(t1, 1, function(x) paste0(x[1], ":", x[2]))
            t3 <- gsub("\\[", "\\\\[", t2)
            t4 <- gsub("\\]", "\\\\]", t3)
            t5 <- list()
            for(k in 1:length(t3)){
              t5[[k]] <- grep(t4[k], colnames(design), value=FALSE)  
            }
            t6 <- unlist(t5)
            t7 <- unique(c(1:p,t6))
            XX <- design[,t7]
          }
        } else if(hierarchy == "Anti") {
          # Features: Include only interactions between non-selected main effects 
          if(length(iS1) != 0){
            t0 <- gsub("\\[", "\\\\[", iS1)
            t1 <- gsub("\\]", "\\\\]", t0)
            t2 <- list()
            for(k in 1:length(t1)){
              t2[[k]] <- grep(t1[k], colnames(design), value=FALSE)  
            }
            t3 <- unlist(t2)
            t4 <- unique(c(1:p,(1:ncol(design))[!(1:ncol(design) %in% t3)]))
            xx <- design[,t4]
          }else{
            xx <- design
          }
          lambda.length <- 20
          Xs <- gsr(x=xx, y=y, lambda=seq(0.1,0.9,length.out = lambda.length)*n)
          if(length(Xs) != lambda.length){
            XX <- design
          } else {
            zz <- lapply(1:length(Xs), function(x) length(Xs[[x]]))
            sel <- max(which(unlist(zz) > p))
            XX <- design[,Xs[[sel]]]
          }
        } else if(hierarchy == "No"){
          # Features: Include all main effects and interactions 
          # with no shrinkage of selected main effects
          if(scree == "GSR"){
            # Screening: The global strong rule
            lambda.length <- 20
            Xs <- gsr(x=design, y=y, lambda=seq(0.1,0.9,length.out = lambda.length)*n)
            if(length(Xs) != lambda.length){
              XX <- design
            } else {
              zz <- lapply(1:length(Xs), function(x) length(Xs[[x]]))
              sel <- max(which(unlist(zz) > p))
              XX <- design[,Xs[[sel]]]
            }
          } else if(scree == "sDPP"){
            # Screening: The sequential DPP
            for(l in seq(0.2,1,length.out = 10)){
              XX <- design[,sdpp(x=design, y=y, lambda=l*n, lambdaprim=n)]
              if(ncol(XX) < ncol(design)) break
            }
          } else {
            XX <- design
          }
        } else if(hierarchy == "Mains") {
          # Number of true interaction discoveries
          Si <- NA
          # False discovery rate (interactions)
          ifdr1 <- NA
        } else if(hierarchy == "Interactions") {
          XX <- design[,-(1:p)]
        }
        
        # Penalty factors that multiplies lambda to allow 
        # no shrinkage of selected main effects
        if(hierarchy != "No"){
          pf <- ifelse(colnames(XX) %in% iS1, 0, 1)
        } else {
          pf <- rep(1, ncol(XX))
        }
        
        # Lasso regression 
        if(W==1){
          set.seed(2000+i) #
          fit2c <- try(cv.glmnet(XX, y, family="gaussian", standardize = TRUE, 
                                 penalty.factor = pf), silent=TRUE)
          if(inherits(fit2c, "try-error")) next
          # Zero and non-zero coefficients
          c2c <- coef(fit2c, s = "lambda.1se") #
        } else {
          MSEs <- matrix(ncol=W, nrow=100)
          LAMBDAs <- matrix(ncol=W, nrow=100)
          for(l in 1:W){
            set.seed(i*20000+l) 
            fit2c <- cv.glmnet(XX, y, family="gaussian", standardize = TRUE, 
                               penalty.factor = pf)
            #fit2c <- try(cv.glmnet(XX, y, family="gaussian", standardize = TRUE, 
            #                       nlambda = 100,
            #                       penalty.factor = pf), silent=TRUE)
            #if(inherits(fit2c, "try-error")) next
            # Hvordan bestemmes hvor lang lambda sekvensen skal være??????????????
            MSEs[1:length(fit2c$cvm),l] <- fit2c$cvm
            LAMBDAs[1:length(fit2c$lambda),l] <- fit2c$lambda
            #ncols <- ncols+1
          }
          lambda.min <- rowMeans(LAMBDAs, na.rm=TRUE)[which.min(rowMeans(MSEs, na.rm=TRUE))]
          # Zero and non-zero coefficients
          c2c <- coef(fit2c, s = lambda.min) #coef(fit2c, s = "lambda.1se") #
        }
        # Non-zero coefficients
        nZ <- names(as.matrix(c2c)[as.matrix(c2c) != 0,])
        nZ <- nZ[nZ != "(Intercept)"]
        # Non-zero interactions
        nZi <- nZ[nZ %in% paste0("X[, ", nz$a, "]:X[, ", nz$b, "]")]
        # Selected features (mains and interactions)
        iSelect <- colnames(XX)[colnames(XX) %in% clean(nZ)]
        # Number of false interaction discoveries
        Vi[j] <- sum(!(nZi %in% c("X[, 1]:X[, 2]")))
        # Number of true interaction discoveries
        Si[j] <- sum(nZi %in% c("X[, 1]:X[, 2]"))
        # False discovery rate (interactions)
        ifdr[j] <- Vi[j]/(Vi[j]+Si[j])
        # Keep track
        cat(".")
    }
    # Return
    return(list(ifdr=ifdr, Vi=Vi, Si=Si))
  })
    
  # Running time
  tim <- difftime(Sys.time(),a,units="mins")
  cat(tim, "mins")
  #Return
  list(IFDR=IFDR, ETIM=tim, n=n, p=p, N=N, W=W)
  
 }


###################################################
### code chunk number 6: step0.Rnw:393-394
###################################################
W <- 10


###################################################
### code chunk number 7: step0.Rnw:397-415
###################################################
#############################
# Truth is strong hierarchy #
#############################
truth <- "Weak"
######################
# Stratify on beta12 #
######################
# Function of beta1 #
#####################
bet1 <- seq(0.1,1.1,by=0.1)
bet2 <- NULL
bet12 <- seq(0.2,2.2,by=0.4)

strat <- "bet12"
xaxis <- "bet1"

scree <- "No"



