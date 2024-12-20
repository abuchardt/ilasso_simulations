### R code from vignette source '/localhome/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step1'
### Encoding: UTF-8

###################################################
### code chunk number 1: Stangle (eval = FALSE)
###################################################
## #setwd("/home/ann-sophie/Documents/wip/lasso/interactions/Weak/")
## setwd("/run/user/551438008/gvfs/smb-share:domain=SUND,server=sund.root.ku.dk,share=users,user=zpt949/LINUX/ifsv/zpt949/Documents/wip/lasso/interactions/Weak")
## #setwd("/localhome/ifsv/zpt949/Documents/wip/lasso/interactions/Weak")
## Stangle("/localhome/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step1")
## #Stangle("/run/user/551438008/gvfs/smb-share:domain=SUND,server=sund.root.ku.dk,share=users,user=zpt949/LINUX/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step1")
## #Stangle("/home/ifsv/zpt949/Documents/wip/lasso/interactions/Weak/step1")
## #Stangle("/home/ann-sophie/Documents/wip/lasso/interactions/Weak/step1")
## #Stangle("step1")


###################################################
### code chunk number 2: step1.Rnw:19-28
###################################################
rm(list=ls())
gc()
source(file="step0.R")

seed.in <- 1234

i <- as.numeric(Sys.getenv("SGE_TASK_ID")) # Index of current process
N <- as.numeric(Sys.getenv("SGE_TASK_LAST")) # Total number of processes
name <- Sys.getenv("JOB_NAME")


###################################################
### code chunk number 3: step1.Rnw:32-41
###################################################
#load(file=paste0("catStep1",strat,xaxis,truth,W,N,p,n,".RData"))

STEP1 <- step1(truth=truth, strat=strat, xaxis=xaxis, 
                 scree=scree, X=X, design=design, i=i, W=W, N=N, p=p, n=n, 
                 bet1=bet1, bet2=bet2, bet3=NULL, bet12=bet12, rh=NULL, mu=NULL)   
#save(STEP1, file=paste0("catStep1",strat,xaxis,truth,W,N,p,n,".RData")) 
name2 <- paste0("catStep1",strat,xaxis,truth,W,N,p,n)
filename <- paste('output/',name2,'-',i,'.rda',sep="")
save(STEP1,file=filename)


