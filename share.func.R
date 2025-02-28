library(parallel)
library(nleqslv)
library(GLDEX)
library(dplyr)
library(pbmcapply)
library(sampling)
library(partitions)
library(data.table)
library(msm)
library(cubature)
library(rje)
library(xtable)
library(MASS)
library(caret)
library(rpart)
library(randomForest)
library(survey)
library(WeightIt)
# Get parameter setting-----------------------------------------------------------

getstats.f <- function(nsample, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina) { # p.z=P(Z(A)=last state; X)
  
  f_gamma0 <- function(gamma0, gamma1, p.x1.1, p.x2.1) {
    expit(gamma0+gamma1*1)*p.x2.1+
      expit(gamma0+gamma1*0)*(1-p.x2.1)-p.x1.1
  }
  
  gamma0 <- uniroot(f_gamma0, lower = -1000, upper = 1000, 
                    gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1)$root
  
  
  f_lam1 <- function(lam1, r, beta1, beta2, p.z, p.x2.1, gamma0, gamma1, inA){
    lam1 <- lam1
    lam2 <- (r^1)*lam1
    Q <- matrix(0, nrow=3, ncol=3)
    Q[1,2] <- lam1
    Q[2,3] <- lam2
    diag(Q) <- (-1)*apply(Q, 1, sum)
    P11 <- MatrixExp(Q*exp(beta1*1+beta2*1), t=inA, method="series")
    P10 <- MatrixExp(Q*exp(beta1*1+beta2*0), t=inA, method="series")
    P01 <- MatrixExp(Q*exp(beta1*0+beta2*1), t=inA, method="series")
    P00 <- MatrixExp(Q*exp(beta1*0+beta2*0), t=inA, method="series")
    out <- P11[1,3]*expit(gamma0+gamma1*1)*p.x2.1 +
      P10[1,3]*expit(gamma0+gamma1*0)*(1-p.x2.1) +
      P01[1,3]*(1-expit(gamma0+gamma1*1))*p.x2.1 +
      P00[1,3]*(1-expit(gamma0+gamma1*0))*(1-p.x2.1) - p.z
    return(out)
  }
  
  lam1 <- uniroot(f_lam1, lower = 0, upper = 100,
                  r=r, beta1=beta1, beta2=beta2, p.z=p.z, p.x2.1=p.x2.1, gamma0=gamma0, gamma1=gamma1, inA=inA)$root
  lam2 <- (r^1)*lam1
  Q <- matrix(0, nrow=3, ncol=3)
  Q[1,2] <- lam1
  Q[2,3] <- lam2
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  outstats <- NULL
  outstats$nsample <- nsample
  outstats$lam <- data.frame(lam1, lam2)
  outstats$Q <- Q
  outstats$A <- inA
  outstats$beta1 <- beta1
  outstats$beta2 <- beta2
  outstats$gamma0 <- gamma0
  outstats$gamma1 <- gamma1
  outstats$p.x1.1 <- p.x1.1
  outstats$p.x2.1 <- p.x2.1
  outstats$lam.ina <- lam.ina
  outstats$p.z <- p.z
  outstats$r <- r
  return(outstats)
}



# Generate a data set-----------------------------------------------------------


generatedata.f <- function(nsample, inQ, inA, lam.ina, gamma0, gamma1, p.x2.1, beta1, beta2) {
  
  # a <- rexp(100, lam.ina)
  # cumsum(a)
  # ina <- c(cumsum(a)[cumsum(a)<1], inA)
  
  x2 <- rbinom(nsample, 1, p.x2.1)
  x1 <- rbinom(nsample, 1, expit(gamma0+gamma1*x2))
  
  maxstate <- ncol(inQ)
  
  rawdata <- NULL
  for (i in 1:nsample) {
    inQx <- inQ*exp(beta1*x1[i]+beta2*x2[i])
    mat <- sim.msm(inQx, maxtime=inA, start=1, mintime=0)
    mat <- data.frame(times=mat$times, states=mat$states)
    
    ina <- 0
    while (ina == 0) {
      a <- rexp(1000, lam.ina)
      inaa <- cumsum(a)[cumsum(a)<=inA]
      if (length(inaa) >= 1) { break; }
    }
    
    ina <- c(0, inaa)
    
    states <- 1
    times  <- ina[1]
    for (j in 2:length(ina)) {
      inwhichstate <- mat$states[sum(mat$times <= ina[j])]
      
      states <- c(states, inwhichstate) 
      times  <- c(times, ina[j])
      if ( inwhichstate == maxstate ) { break }
    }
    
    datai <- data.frame(id=i, states, times, x1=x1[i], x2=x2[i])
    rawdata <- rbind(rawdata, datai)
  } 
  rawdata <- rawdata %>% unique()
  
  event <- rawdata
  event$estop  <- event$times
  event$estart <- c(0, event$times[-nrow(event)])
  event$to    <- event$states
  event$from  <- c(0, event$states[-nrow(event)])
  event <- event[event$estop > 0, c("id","estart","estop","from","to", "x1", "x2")]
  
  rawdata_strata <- rawdata %>% group_by(id) %>% mutate(Z.am = max(states)) %>% dplyr::select(id, x2, Z.am) %>% slice(1) %>% 
    mutate(strata = case_when(Z.am == 1 & x2 == 0 ~ 1,
                              Z.am == 1 & x2 == 1 ~ 2,
                              Z.am == 2 & x2 == 0 ~ 3,
                              Z.am == 2 & x2 == 1 ~ 4,
                              Z.am == 3 & x2 == 0 ~ 5,
                              Z.am == 3 & x2 == 1 ~ 6)) %>% dplyr::select(id, strata)
  
  event <- left_join(event, rawdata_strata, by="id")
  rawdata <- left_join(rawdata, rawdata_strata, by="id")
  
  outdata <- NULL
  outdata$event <- event
  outdata$rawdata <- rawdata
  return(outdata)
}


# Select phase 2 samples
p2data.SRS.f <- function(indata, nsample2) {
  event <- indata$event
  rawdata <- indata$rawdata
  nsample <- length(unique(event$id))
  
  # Simple random sampling
  id.random <- sample(1:nsample, nsample2, replace = FALSE)
  event$SRS <- ifelse(event$id %in% id.random, 1, 0)
  rawdata$SRS <- ifelse(rawdata$id %in% id.random, 1, 0)
  
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}



p2data.OPT.f <- function(indata, nsample2, inQ, gamma0, gamma1) {
  event <- indata$event
  rawdata <- indata$rawdata
  maxstate <- ncol(inQ)
  nsample <- length(unique(event$id))
  
  # optimal sampling
  fitm <- msm(states ~ times, subject=id, data=rawdata, 
              covariates = ~ x2,
              constraint = list(x2 = c(1,1)),
              qmatrix=inQ, obstype=1, opt.method="fisher",
              gen.inits=TRUE, center = FALSE)
  # fitm$estimates.t
  # qmatrix.msm(fitm,covariates=list (x2 = 0), ci="none" )
  # pmatrix.msm(fitm, t=0.2, covariates=list (x2 = 0))
  
  SUBJ  <- unique(as.character(event$id))
  
  Mmu.f <- function(whichsubj, locdata, maxstate) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x2 <- locdata$x2[1]
    nlen <- length(locdata$w)
    
    Qxmat <- qmatrix.msm(fitm, covariates=list (x2 = x2), ci="none" )
    
    lamx <- NULL
    for (j in 1:maxstate) {
      for (k in 1:maxstate) {
        if ( (j != k) && (inQ[j,k] > 0) ) {
          lamx <- c(lamx, Qxmat[j,k])
        }
      }
    } 
    
    lamx <- c(lamx, 0)
    Mmu <- 0
    for (i in 1:nlen) {
      Pxmat <- pmatrix.msm(fitm, t=locdata$w[i], covariates=list (x2 = x2))
      
      # dP_klr(X)/dmu at j
      dPxmuj.f <- function(k,l,j){
        Ckjl <- prod( lamx[k:(l-1)] )/prod( fun.zero.omit(lamx[k:l]-lamx[j]) ) 
        Ckjl <- ifelse(k==l & k==j, 1, Ckjl)
        
        dPxmuj <- Ckjl*exp(-lamx[j]*locdata$w[i])*(-lamx[j]*locdata$w[i])
        return(dPxmuj)
      }
      
      # dPxmu <- 0
      # for (g in locdata$from[i]:locdata$to[i]) {
      #   dPxmu <- dPxmu + dPxmuj.f( k=locdata$from[i], l=locdata$to[i], j=g)
      # }
      
      dPxmu <- sum(apply(X=as.matrix(c(locdata$from[i]:locdata$to[i])), 1, FUN=dPxmuj.f, k=locdata$from[i], l=locdata$to[i]))
      #print(1/Pxmat[locdata$from[i],locdata$to[i]]*dPxmu)
      Mmu <- Mmu + 1/Pxmat[locdata$from[i],locdata$to[i]]*dPxmu
    }
    
    return(Mmu)
  }
  
  Mmu <- lapply(SUBJ, Mmu.f, locdata=event, maxstate=maxstate)
  # Mmu <- mclapply(SUBJ, Mmu.f, locdata=event, maxstate=maxstate, mc.cores=cores.num)
  Mmu <- as.numeric(Mmu)
  
  # 1/(lamx[1]*(exp(-lamx[2]*0.2) - exp(-lamx[1]*0.2) )/(lamx[1]-lamx[2]) ) * lamx[1]/(lamx[1]-lamx[2])*(-lamx[2]*0.2*exp(-lamx[2]*0.2) + lamx[1]*0.2*exp(-lamx[1]*0.2))+
  # 1/exp(-lamx[2]*0.2)*(-lamx[2]*0.2*exp(-lamx[2]*0.2))*2+
  # 1/(1-exp(-lamx[2]*0.2))*(lamx[2]*0.2*exp(-lamx[2]*0.2))
  
  x2_dat <- rawdata %>% group_by(id) %>% slice(1) 
  x2 <- x2_dat$x2
  
  Mmu[which(x2==0)] <- Mmu[which(x2==0)]*sqrt(expit(gamma0+gamma1*0)*(1-expit(gamma0+gamma1*0)))
  Mmu[which(x2==1)] <- Mmu[which(x2==1)]*sqrt(expit(gamma0+gamma1*1)*(1-expit(gamma0+gamma1*1)))
  
  dat <- data.frame(id=1:nsample, Mmu=Mmu, x2=x2) %>% arrange(Mmu)
  dat.x2.1 <- dat %>% filter(x2==1) %>% arrange(Mmu)
  dat.x2.0 <- dat %>% filter(x2==0) %>% arrange(Mmu)
  
  id.optimal <- c(dat.x2.1$id[1:(nsample2/4)], dat.x2.1$id[(nrow(dat.x2.1)-(nsample2/4)+1):nrow(dat.x2.1)],
                  dat.x2.0$id[1:(nsample2/4)], dat.x2.0$id[(nrow(dat.x2.0)-(nsample2/4)+1):nrow(dat.x2.0)])
  
  event$OPT <- ifelse(event$id %in% id.optimal, 1, 0)
  rawdata$OPT <- ifelse(rawdata$id %in% id.optimal, 1, 0)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}


p2data.OPT.new.f <- function(indata, nsample2, inQ, gamma0, gamma1) {
  event <- indata$event
  rawdata <- indata$rawdata
  maxstate <- ncol(inQ)
  nsample <- length(unique(event$id))
  
  # optimal sampling
  fitm <- msm(states ~ times, subject=id, data=rawdata, 
              covariates = ~ x2,
              constraint = list(x2 = c(1,1)),
              qmatrix=inQ, obstype=1, opt.method="fisher",
              gen.inits=TRUE, center = FALSE)
  # fitm$estimates.t
  # qmatrix.msm(fitm,covariates=list (x2 = 0), ci="none" )
  # pmatrix.msm(fitm, t=0.2, covariates=list (x2 = 0))
  
  SUBJ  <- unique(as.character(event$id))
  
  Mmu.f <- function(whichsubj, locdata, maxstate) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x2 <- locdata$x2[1]
    nlen <- length(locdata$w)
    
    Qxmat <- qmatrix.msm(fitm, covariates=list (x2 = x2), ci="none" )
    
    lamx <- NULL
    for (j in 1:maxstate) {
      for (k in 1:maxstate) {
        if ( (j != k) && (inQ[j,k] > 0) ) {
          lamx <- c(lamx, Qxmat[j,k])
        }
      }
    } 
    
    lamx <- c(lamx, 0)
    Mmu <- 0
    for (i in 1:nlen) {
      Pxmat <- pmatrix.msm(fitm, t=locdata$w[i], covariates=list (x2 = x2))
      
      # dP_klr(X)/dmu at j
      dPxmuj.f <- function(k,l,j){
        Ckjl <- prod( lamx[k:(l-1)] )/prod( fun.zero.omit(lamx[k:l]-lamx[j]) ) 
        Ckjl <- ifelse(k==l & k==j, 1, Ckjl)
        
        dPxmuj <- Ckjl*exp(-lamx[j]*locdata$w[i])*(-lamx[j]*locdata$w[i])
        return(dPxmuj)
      }
      
      # dPxmu <- 0
      # for (g in locdata$from[i]:locdata$to[i]) {
      #   dPxmu <- dPxmu + dPxmuj.f( k=locdata$from[i], l=locdata$to[i], j=g)
      # }
      
      dPxmu <- sum(apply(X=as.matrix(c(locdata$from[i]:locdata$to[i])), 1, FUN=dPxmuj.f, k=locdata$from[i], l=locdata$to[i]))
      #print(1/Pxmat[locdata$from[i],locdata$to[i]]*dPxmu)
      Mmu <- Mmu + 1/Pxmat[locdata$from[i],locdata$to[i]]*dPxmu
    }
    
    return(Mmu)
  }
  
  Mmu <- lapply(SUBJ, Mmu.f, locdata=event, maxstate=maxstate)
  # Mmu <- mclapply(SUBJ, Mmu.f, locdata=event, maxstate=maxstate, mc.cores=cores.num)
  Mmu <- as.numeric(Mmu)
  
  # 1/(lamx[1]*(exp(-lamx[2]*0.2) - exp(-lamx[1]*0.2) )/(lamx[1]-lamx[2]) ) * lamx[1]/(lamx[1]-lamx[2])*(-lamx[2]*0.2*exp(-lamx[2]*0.2) + lamx[1]*0.2*exp(-lamx[1]*0.2))+
  # 1/exp(-lamx[2]*0.2)*(-lamx[2]*0.2*exp(-lamx[2]*0.2))*2+
  # 1/(1-exp(-lamx[2]*0.2))*(lamx[2]*0.2*exp(-lamx[2]*0.2))
  
  x2_dat <- rawdata %>% group_by(id) %>% slice(1) 
  x2 <- x2_dat$x2
  
  Mmu[which(x2==0)] <- Mmu[which(x2==0)]*sqrt(expit(gamma0+gamma1*0)*(1-expit(gamma0+gamma1*0)))
  Mmu[which(x2==1)] <- Mmu[which(x2==1)]*sqrt(expit(gamma0+gamma1*1)*(1-expit(gamma0+gamma1*1)))
  
  dat <- data.frame(id=1:nsample, Mmu=Mmu, x2=x2) %>% arrange(Mmu)
  dat.x2.1 <- dat %>% filter(x2==1) %>% arrange(Mmu)
  dat.x2.0 <- dat %>% filter(x2==0) %>% arrange(Mmu)
  
  id.optimal <- c(dat.x2.1$id[1:(nsample2/4)], dat.x2.1$id[(nrow(dat.x2.1)-(nsample2/4)+1):nrow(dat.x2.1)],
                  dat.x2.0$id[1:(nsample2/4)], dat.x2.0$id[(nrow(dat.x2.0)-(nsample2/4)+1):nrow(dat.x2.0)])
  
  event$OPT <- ifelse(event$id %in% id.optimal, 1, 0)
  rawdata$OPT <- ifelse(rawdata$id %in% id.optimal, 1, 0)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}



p2data.PROP.f <- function(indata, nsample2) {
  event <- indata$event
  rawdata <- indata$rawdata
  maxstate <- ncol(inQ)
  nsample <- length(unique(event$id))
  
  # proportional sampling
  rawdata_strata <- rawdata %>% group_by(id) %>% mutate(Z.am = max(states)) %>% dplyr::select(id, x2, Z.am) %>% slice(1) %>% 
    mutate(strata = case_when(Z.am == 1 & x2 == 0 ~ 1,
                              Z.am == 1 & x2 == 1 ~ 2,
                              Z.am == 2 & x2 == 0 ~ 3,
                              Z.am == 2 & x2 == 1 ~ 4,
                              Z.am == 3 & x2 == 0 ~ 5,
                              Z.am == 3 & x2 == 1 ~ 6))
  strata.prop <- round(summary(as.factor(rawdata_strata$strata))/nsample*nsample2)
  num.prop <- c(strata.prop[1:5], nsample2-sum(strata.prop[1:5]))
  id.prop <- c(sample(rawdata_strata$id[rawdata_strata$strata == 1], strata.prop[1], replace = FALSE),
               sample(rawdata_strata$id[rawdata_strata$strata == 2], strata.prop[2], replace = FALSE),
               sample(rawdata_strata$id[rawdata_strata$strata == 3], strata.prop[3], replace = FALSE),
               sample(rawdata_strata$id[rawdata_strata$strata == 4], strata.prop[4], replace = FALSE),
               sample(rawdata_strata$id[rawdata_strata$strata == 5], strata.prop[5], replace = FALSE),
               sample(rawdata_strata$id[rawdata_strata$strata == 6], strata.prop[6], replace = FALSE))
  event$PROP <- ifelse(event$id %in% id.prop, 1, 0)
  rawdata$PROP <- ifelse(rawdata$id %in% id.prop, 1, 0)
  
  # # balanced sampling
  # num.bal <- c(rep(round(nsample2/6), 5), nsample2-sum(rep(round(nsample2/6), 5)))
  # 
  # id.bal <- c(sample(rawdata_strata$id[rawdata_strata$strata == 1], num.bal[1], replace = FALSE),
  #             sample(rawdata_strata$id[rawdata_strata$strata == 2], num.bal[2], replace = FALSE),
  #             sample(rawdata_strata$id[rawdata_strata$strata == 3], num.bal[3], replace = FALSE),
  #             sample(rawdata_strata$id[rawdata_strata$strata == 4], num.bal[4], replace = FALSE),
  #             sample(rawdata_strata$id[rawdata_strata$strata == 5], num.bal[5], replace = FALSE),
  #             sample(rawdata_strata$id[rawdata_strata$strata == 6], num.bal[6], replace = FALSE))
  # 
  # event$BAL <- ifelse(event$id %in% id.bal, 1, 0)
  # rawdata$BAL <- ifelse(rawdata$id %in% id.bal, 1, 0)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}


p2data.BAL.f <- function(indata, nsample2) {
  num.bal <- nsample2/6
  
  rawdata <- indata$rawdata
  event <- indata$event
  rawdata_strata <- rawdata %>% group_by(id)  %>% dplyr::select(id, strata) %>% slice(1) 
  nstrata <- summary(as.factor(rawdata_strata$strata))
  
  id.bal1 <- NULL
  for (i in 1:6) {
    if (nstrata[i] <= num.bal) {
      id.select <- rawdata_strata$id[rawdata_strata$strata == i]
    } else {
      id.select <- sample(rawdata_strata$id[rawdata_strata$strata == i], num.bal, replace = FALSE)
    }
    id.bal1 <- c(id.bal1,  id.select)
  }
  
  # sample additional items from larger groups with a surplus of items to achieve the desired total sample size
  
  remaining_sample_size <- nsample2 - length(id.bal1)
  if (remaining_sample_size == 0){
    id.bal <- id.bal1
  } else {
    remaining_sample <- rawdata_strata %>% filter(!id %in% id.bal1)
    remaining_nstrata <- summary(as.factor(remaining_sample$strata))
    remaining_strata <- sort(unique(remaining_sample$strata))
    
    # find the allocation which meets the criteria that the remaining group sizes are sufficient to choose from and has the lower variance
    cmp <- compositions(n = remaining_sample_size, m = length(remaining_nstrata))
    std <- apply(cmp, 2, sd)
    cmp2 <- cmp[ , order(std)]
    
    valid_indexes <- apply(cmp2, 2, function(x) all(x <= remaining_nstrata))
    sampling_sizes <- cmp2[, which(valid_indexes)[1]]
    
    rm(cmp, cmp2, std)
    
    for (j in 1:length(remaining_nstrata)) {
      strata.j <- remaining_strata[j]
      id.select <- sample(remaining_sample$id[remaining_sample$strata == strata.j], sampling_sizes[j], replace = FALSE)
      id.bal1 <- c(id.bal1,  id.select)
    }
    id.bal <- id.bal1
  }
  
  
  event$BAL <- ifelse(event$id %in% id.bal, 1, 0)
  rawdata$BAL <- ifelse(rawdata$id %in% id.bal, 1, 0)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}


# optimal design for IPW method
p2data.IPW.f <- function(indata, nsample2, inQ, nlambda, nbeta, gamma0, gamma1) {
  event <- indata$event
  rawdata <- indata$rawdata
  event <- event %>% mutate(strata2 = case_when(strata == 1 | strata == 2 ~ 1,
                                                strata == 3 | strata == 4 ~ 2,
                                                strata == 5 | strata == 6 ~ 3))
  rawdata <- rawdata %>% mutate(strata2 = case_when(strata == 1 | strata == 2 ~ 1,
                                                strata == 3 | strata == 4 ~ 2,
                                                strata == 5 | strata == 6 ~ 3))
  maxstate <- ncol(inQ)
  nsample <- length(unique(event$id))
  ntheta <- nlambda + nbeta
  
  n2a <- round(nsample2*1/3/3)
  
  rawdata_strata <- rawdata %>% group_by(id)  %>% dplyr::select(id, strata, strata2) %>% slice(1) 
  nstrata <- summary(as.factor(rawdata_strata$strata2))
  
  id.2a <- NULL
  for (i in 1:3) {
    id.2a.select <- sample(rawdata_strata$id[rawdata_strata$strata2 == i], n2a, replace = FALSE)
    id.2a <- c(id.2a,  id.2a.select)
  }
  
  rawdata.2a <- rawdata %>% filter(id %in% id.2a)
  
  fitm <- msm(states ~ times, subject=id, data=rawdata.2a, 
              covariates = ~ x1+x2,
              constraint = list(x1 = c(1,1), x2 = c(1,1)),
              qmatrix=inQ, obstype=1, opt.method="fisher",
              gen.inits=TRUE, center = FALSE)
  # fitm$estimates.t
  # qmatrix.msm(fitm,covariates=list (x2 = 0), ci="none" )
  # pmatrix.msm(fitm, t=0.2, covariates=list (x2 = 0))
  
  beta <- fitm$estimates.t[c(3,5)]
  lam <- fitm$estimates.t[1:2]
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ[j, k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- lam[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  SUBJ  <- unique(as.character(event$id))
  
  MAT2.f <- function(whichsubj, locdata, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
  
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    MMi <- EMMi <- matrix(0, nrow=ntheta, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              EMMi[u,j] <- EMMi[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
            }
          }
        }
        
      }
    }
    
    outMAT <- rbind(SSi1, MMi, EMMi)
    return(outMAT)
  }
  
  MAT2 <- lapply(SUBJ, MAT2.f, locdata=event, nlambda=nlambda, nbeta=nbeta, 
                 maxstate=maxstate, Q=Q, inbeta=beta)
  nSUBJ <- length(unique(as.character(event$id)))
  
  Ai <- lapply(MAT2, function(x){ x[c(6:9), c(1:(ntheta))]   })
  EA <- Reduce('+', Ai)/nSUBJ
  hi <- lapply(MAT2, function(x){ (ginv(EA) %*% x[1, c(1:(ntheta))])[3]   })
  hi <- unlist(hi)
  
  
  
  rawdata_strata$hi <- hi
  rawdata_strata <- rawdata_strata %>% group_by(strata2) %>% mutate(sd.hi = sd(hi))
  sd.hi <- rawdata_strata %>% group_by(strata2) %>% slice(1)
  # test <- rawdata_strata %>% filter(strata == 5)
  
  nstrata.tmp <- ceiling(nstrata*sd.hi$sd.hi/sum(nstrata*sd.hi$sd.hi)*(nsample2)-n2a)
  nstrata.tmp[nstrata.tmp < 0] <- 0
  
  remaining_deduction <- -(nsample2-n2a*3-sum(nstrata.tmp))
  remaining_addition <- nsample2-n2a*3-sum(nstrata.tmp)
  
  # Loop until the total deduction reaches `remaining_deduction`
  while (remaining_deduction > 0) {
    # Find the index of the minimum value in X that is greater than 0
    min_index <- which.min(nstrata.tmp[nstrata.tmp > 0])
    actual_index <- which(nstrata.tmp > 0)[min_index]
    
    # Determine the amount we can deduct from this element
    deduct_amount <- min(remaining_deduction, nstrata.tmp[actual_index])
    
    # Deduct and update remaining deduction
    nstrata.tmp[actual_index] <- nstrata.tmp[actual_index] - deduct_amount
    remaining_deduction <- remaining_deduction - deduct_amount
  }
  
  
  # Loop until the total addition reaches `remaining_addition`
  while (remaining_addition > 0) {
    # Find the index of the maximum value in X
    max_index <- which.max(nstrata.tmp)
    
    # Add as much as possible to this element
    nstrata.tmp[max_index] <- nstrata.tmp[max_index] + remaining_addition
    
    # Set remaining addition to 0, as we've added all to the largest element
    remaining_addition <- 0
  }
  
  
  rawdata_remain <- rawdata_strata %>% filter(!id %in% id.2a)
  
  id.2b <- NULL
  for (i in 1:3) {
    id.2b.select <- sample(rawdata_remain$id[rawdata_remain$strata2 == i], nstrata.tmp[i], replace = FALSE)
    id.2b <- c(id.2b,  id.2b.select)
  }
  
  id.optimal <- c(id.2a, id.2b)

  event$IPW.OPT <- ifelse(event$id %in% id.optimal, 1, 0)
  rawdata$IPW.OPT <- ifelse(rawdata$id %in% id.optimal, 1, 0)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}




p2data.IPW.true.f <- function(indata, stats, nsample2, inQ, nlambda, nbeta, gamma0, gamma1) {
  event <- indata$event
  rawdata <- indata$rawdata
  event <- event %>% mutate(strata2 = case_when(strata == 1 | strata == 2 ~ 1,
                                                strata == 3 | strata == 4 ~ 2,
                                                strata == 5 | strata == 6 ~ 3))
  rawdata <- rawdata %>% mutate(strata2 = case_when(strata == 1 | strata == 2 ~ 1,
                                                    strata == 3 | strata == 4 ~ 2,
                                                    strata == 5 | strata == 6 ~ 3))
  maxstate <- ncol(inQ)
  nsample <- length(unique(event$id))
  ntheta <- nlambda + nbeta
  
  
  rawdata_strata <- rawdata %>% group_by(id)  %>% dplyr::select(id, strata, strata2) %>% slice(1) 
  nstrata <- summary(as.factor(rawdata_strata$strata2))
  
  beta <- c(stats$beta1, stats$beta2) 
  lam <- as.numeric(stats$lam)
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ[j, k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- lam[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  SUBJ  <- unique(as.character(event$id))
  
  MAT2.f <- function(whichsubj, locdata, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    MMi <- EMMi <- matrix(0, nrow=ntheta, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              EMMi[u,j] <- EMMi[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
            }
          }
        }
        
      }
    }
    
    outMAT <- rbind(SSi1, MMi, EMMi)
    return(outMAT)
  }
  
  MAT2 <- lapply(SUBJ, MAT2.f, locdata=event, nlambda=nlambda, nbeta=nbeta, 
                 maxstate=maxstate, Q=Q, inbeta=beta)
  nSUBJ <- length(unique(as.character(event$id)))
  
  Ai <- lapply(MAT2, function(x){ x[c(6:9), c(1:(ntheta))]   })
  EA <- Reduce('+', Ai)/nSUBJ

  
  hi <- lapply(MAT2, function(x){ (ginv(EA) %*% x[1, c(1:(ntheta))])[3]   })
  hi <- unlist(hi)
  
  rawdata_strata$hi <- hi
  rawdata_strata <- rawdata_strata %>% group_by(strata2) %>% mutate(sd.hi = sd(hi))
  sd.hi <- rawdata_strata %>% group_by(strata2) %>% slice(1)
  # test <- rawdata_strata %>% filter(strata == 5)
  
  nstrata.tmp <- ceiling(nstrata*sd.hi$sd.hi/sum(nstrata*sd.hi$sd.hi)*(nsample2))

  remaining_deduction <- -(nsample2-sum(nstrata.tmp))
  remaining_addition <- nsample2-sum(nstrata.tmp)
  
  # Loop until the total deduction reaches `remaining_deduction`
  while (remaining_deduction > 0) {
    # Find the index of the minimum value in X that is greater than 0
    min_index <- which.min(nstrata.tmp[nstrata.tmp > 0])
    actual_index <- which(nstrata.tmp > 0)[min_index]
    
    # Determine the amount we can deduct from this element
    deduct_amount <- min(remaining_deduction, nstrata.tmp[actual_index])
    
    # Deduct and update remaining deduction
    nstrata.tmp[actual_index] <- nstrata.tmp[actual_index] - deduct_amount
    remaining_deduction <- remaining_deduction - deduct_amount
  }
  
  
  # Loop until the total addition reaches `remaining_addition`
  while (remaining_addition > 0) {
    # Find the index of the maximum value in X
    max_index <- which.max(nstrata.tmp)
    
    # Add as much as possible to this element
    nstrata.tmp[max_index] <- nstrata.tmp[max_index] + remaining_addition
    
    # Set remaining addition to 0, as we've added all to the largest element
    remaining_addition <- 0
  }
  
  
  id.optimal <- NULL
  for (i in 1:3) {
    id.optimal.select <- sample(rawdata_strata$id[rawdata_strata$strata2 == i], nstrata.tmp[i], replace = FALSE)
    id.optimal <- c(id.optimal,  id.optimal.select)
  }
  
  event$IPW.OPT.TRUE <- ifelse(event$id %in% id.optimal, 1, 0)
  rawdata$IPW.OPT.TRUE <- ifelse(rawdata$id %in% id.optimal, 1, 0)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}



p2data.IPW.true6.f <- function(indata, stats, nsample2, inQ, nlambda, nbeta, gamma0, gamma1) {
  event <- indata$event
  rawdata <- indata$rawdata

  maxstate <- ncol(inQ)
  nsample <- length(unique(event$id))
  ntheta <- nlambda + nbeta
  
  
  rawdata_strata <- rawdata %>% group_by(id)  %>% dplyr::select(id, strata) %>% slice(1) 
  nstrata <- summary(as.factor(rawdata_strata$strata))
  
  beta <- c(stats$beta1, stats$beta2) 
  lam <- as.numeric(stats$lam)
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ[j, k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- lam[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  SUBJ  <- unique(as.character(event$id))
  
  MAT2.f <- function(whichsubj, locdata, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    MMi <- EMMi <- matrix(0, nrow=ntheta, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              EMMi[u,j] <- EMMi[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
            }
          }
        }
        
      }
    }
    
    outMAT <- rbind(SSi1, MMi, EMMi)
    return(outMAT)
  }
  
  MAT2 <- lapply(SUBJ, MAT2.f, locdata=event, nlambda=nlambda, nbeta=nbeta, 
                 maxstate=maxstate, Q=Q, inbeta=beta)
  nSUBJ <- length(unique(as.character(event$id)))
  
  Ai <- lapply(MAT2, function(x){ x[c(6:9), c(1:(ntheta))]   })
  EA <- Reduce('+', Ai)/nSUBJ
  
  
  hi <- lapply(MAT2, function(x){ (ginv(EA) %*% x[1, c(1:(ntheta))])[3]   })
  hi <- unlist(hi)
  
  rawdata_strata$hi <- hi
  rawdata_strata <- rawdata_strata %>% group_by(strata) %>% mutate(sd.hi = sd(hi))
  sd.hi <- rawdata_strata %>% group_by(strata) %>% slice(1)
  # test <- rawdata_strata %>% filter(strata == 5)
  
  nstrata.tmp <- ceiling(nstrata*sd.hi$sd.hi/sum(nstrata*sd.hi$sd.hi)*(nsample2))
  
  remaining_deduction <- -(nsample2-sum(nstrata.tmp))
  remaining_addition <- nsample2-sum(nstrata.tmp)
  
  # Loop until the total deduction reaches `remaining_deduction`
  while (remaining_deduction > 0) {
    # Find the index of the minimum value in X that is greater than 0
    min_index <- which.min(nstrata.tmp[nstrata.tmp > 0])
    actual_index <- which(nstrata.tmp > 0)[min_index]
    
    # Determine the amount we can deduct from this element
    deduct_amount <- min(remaining_deduction, nstrata.tmp[actual_index])
    
    # Deduct and update remaining deduction
    nstrata.tmp[actual_index] <- nstrata.tmp[actual_index] - deduct_amount
    remaining_deduction <- remaining_deduction - deduct_amount
  }
  
  
  # Loop until the total addition reaches `remaining_addition`
  while (remaining_addition > 0) {
    # Find the index of the maximum value in X
    max_index <- which.max(nstrata.tmp)
    
    # Add as much as possible to this element
    nstrata.tmp[max_index] <- nstrata.tmp[max_index] + remaining_addition
    
    # Set remaining addition to 0, as we've added all to the largest element
    remaining_addition <- 0
  }
  
  
  id.optimal <- NULL
  for (i in 1:6) {
    id.optimal.select <- sample(rawdata_strata$id[rawdata_strata$strata == i], nstrata.tmp[i], replace = FALSE)
    id.optimal <- c(id.optimal,  id.optimal.select)
  }
  
  event$IPW.OPT.TRUE6 <- ifelse(event$id %in% id.optimal, 1, 0)
  rawdata$IPW.OPT.TRUE6 <- ifelse(rawdata$id %in% id.optimal, 1, 0)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
  
  return(outdata)
}
# Kalbfleisch-Lawless Approach------------------------------------------------------

dQx.f <- function(states, maxstate, nlambda, nbeta, x1, x2, inQ, inbeta1, inbeta2) {
  ntheta <- nlambda+nbeta
  
  outQx <- inQ*exp(inbeta1*x1+inbeta2*x2)
  
  outdQx <- array(data=matrix(0, nrow=maxstate, ncol=maxstate),
                  dim=c(maxstate, maxstate, ntheta))
  
  for (k in 1:nrow(states)) {
    outdQx[states$from[k], states$from[k], k] <- (-1)*exp(inbeta1*x1+inbeta2*x2)
    outdQx[states$from[k], states$to[k], k]   <- 1*exp(inbeta1*x1+inbeta2*x2)
  }  
  
  outdQx[,, nrow(states)+1] <- outQx*x1
  outdQx[,, ntheta] <- outQx*x2
  
  outdata <- NULL
  outdata$outdQx <- outdQx
  outdata$outQx   <- outQx
  
  return(outdata)
}


dPxt.f <- function(u, w, inQx, indQx, maxstate) { # uth parameter
  
  decompQx <- eigen(inQx)
  
  dd    <- decompQx$values
  AA    <- decompQx$vectors
  invAA <- ginv(AA, tol=1e-16)
  
  GG <- invAA %*% indQx[,,u] %*% AA
  
  VV <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (i in 1:maxstate) {
    VV[i,1:maxstate] <- GG[i,1:maxstate]*( exp(dd[i]*w) - exp(dd*w) )/( dd[i] - dd )
    VV[i,1:maxstate] <- ifelse((dd[i] - dd) == 0, 0, VV[i,1:maxstate])
  }
  diag(VV) <- diag(GG)*w*exp(dd*w)
  
  outdPx <- AA %*% VV %*% invAA 
  outdPx <- round(outdPx, 16)
  return(outdPx) 
}


# Likelihood-EM -----------------------------------------------------------

scoring.f <- function(indata, selectR, nlambda, nbeta, ngamma, inQ0, inbeta, ingamma) {
  maxstate <- nrow(inQ0)
  ntheta <- nlambda + nbeta
  
  cur.lam <- NULL
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        cur.lam <- c(cur.lam, inQ0[j,k])
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  cur.beta <- inbeta
  cur.gamma <- ingamma
  
  SUBJ  <- unique(as.character(indata$id))
  nSUBJ <- length(SUBJ)
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, ngamma, maxstate, 
                    cur.Q, pre.Q, cur.beta, pre.beta, cur.gamma, pre.gamma) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    
    cur.outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                      x1=x1, x2=x2, inQ=cur.Q, inbeta1=cur.beta[1], inbeta2=cur.beta[2])
    cur.Qx <- cur.outQ$outQx
    cur.dQx <- cur.outQ$outdQx
    
    
    pre.Qx.x1.1 <- pre.Q*exp(pre.beta[1]*1+pre.beta[2]*x2)
    
    pre.Qx.x1.0 <- pre.Q*exp(pre.beta[1]*0+pre.beta[2]*x2)
    
    # pseudo
    cur.outQ.x1.1 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                           x1=1, x2=x2, inQ=cur.Q, inbeta1=cur.beta[1], inbeta2=cur.beta[2])
    cur.Qx.x1.1 <- cur.outQ.x1.1$outQx
    cur.dQx.x1.1 <- cur.outQ.x1.1$outdQx
    
    cur.outQ.x1.0 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                           x1=0, x2=x2, inQ=cur.Q, inbeta1=cur.beta[1], inbeta2=cur.beta[2])
    cur.Qx.x1.0 <- cur.outQ.x1.0$outQx
    cur.dQx.x1.0 <- cur.outQ.x1.0$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    SSi1.x1.1 <- matrix(0, nrow=1, ncol=ntheta)
    SSi1.x1.0 <- matrix(0, nrow=1, ncol=ntheta)
    
    PD.x1.1 <- 1 # P(D | X1=1, X2)
    PD.x1.0 <- 1 # P(D | X1=0, X2)
    
    SSi2 <- matrix(0, nrow=1, ncol=ngamma)
    SSi2.x1.1 <- matrix(0, nrow=1, ncol=ngamma)
    SSi2.x1.0 <- matrix(0, nrow=1, ncol=ngamma)
    
    for (i in 1:nlen) {
      cur.Pxmat <- round(MatrixExp(cur.Qx, t=locdata$w[i], method="series"), 16)
      
      cur.Pxmat.x1.1 <- round(MatrixExp(cur.Qx.x1.1, t=locdata$w[i], method="series"), 16)
      cur.Pxmat.x1.0 <- round(MatrixExp(cur.Qx.x1.0, t=locdata$w[i], method="series"), 16)
      
      pre.Pxmat.x1.1 <- round(MatrixExp(pre.Qx.x1.1, t=locdata$w[i], method="series"), 16)
      pre.Pxmat.x1.0 <- round(MatrixExp(pre.Qx.x1.0, t=locdata$w[i], method="series"), 16)
      
      PD.x1.1 <- PD.x1.1*pre.Pxmat.x1.1[locdata$from[i], locdata$to[i]]
      PD.x1.0 <- PD.x1.0*pre.Pxmat.x1.0[locdata$from[i], locdata$to[i]]
      
      for (u in 1:ntheta) {
        cur.dPxmatu <- dPxt.f(u, locdata$w[i], cur.Qx, cur.dQx, maxstate)
        cur.dPxmatu.x1.1 <- dPxt.f(u, locdata$w[i], cur.Qx.x1.1, cur.dQx.x1.1, maxstate)
        cur.dPxmatu.x1.0 <- dPxt.f(u, locdata$w[i], cur.Qx.x1.0, cur.dQx.x1.0, maxstate)
        
        SSi1[1,u] <- SSi1[1,u] + cur.dPxmatu[locdata$from[i],locdata$to[i]] / cur.Pxmat[locdata$from[i],locdata$to[i]]
        SSi1.x1.1[1,u] <- SSi1.x1.1[1,u] + cur.dPxmatu.x1.1[locdata$from[i],locdata$to[i]] / cur.Pxmat.x1.1[locdata$from[i],locdata$to[i]]
        SSi1.x1.0[1,u] <- SSi1.x1.0[1,u] + cur.dPxmatu.x1.0[locdata$from[i],locdata$to[i]] / cur.Pxmat.x1.0[locdata$from[i],locdata$to[i]]
      }
    }
    
    
    SSi2[1,1] <- x1-expit(cur.gamma[1]+cur.gamma[2]*x2)
    SSi2[1,2] <- (x1-expit(cur.gamma[1]+cur.gamma[2]*x2))*x2
    
    SSi2.x1.1[1,1] <- 1-expit(cur.gamma[1]+cur.gamma[2]*x2)
    SSi2.x1.1[1,2] <- (1-expit(cur.gamma[1]+cur.gamma[2]*x2))*x2
    
    SSi2.x1.0[1,1] <- 0-expit(cur.gamma[1]+cur.gamma[2]*x2)
    SSi2.x1.0[1,2] <- (0-expit(cur.gamma[1]+cur.gamma[2]*x2))*x2
    
    
    pre.w <- ( PD.x1.1*expit(pre.gamma[1]+pre.gamma[2]*x2) )/( PD.x1.1*expit(pre.gamma[1]+pre.gamma[2]*x2) + PD.x1.0*(1-expit(pre.gamma[1]+pre.gamma[2]*x2)) )
    
    Si1 <- R*SSi1+(1-R)*( SSi1.x1.1*pre.w + SSi1.x1.0*(1-pre.w) )
    Si2 <- R*SSi2+(1-R)*( SSi2.x1.1*pre.w + SSi2.x1.0*(1-pre.w) )
    
    Si <- cbind(Si1, Si2)
    
    return(Si)
  }
  
  iter <- 0
  tol  <- 9999
  while (tol > 1e-05) {
    iter <- iter + 1
    # print(iter)
    
    pre.lam <- cur.lam
    pre.beta <- cur.beta
    pre.gamma <- cur.gamma
    
    pre.all <- c(pre.lam, pre.beta, pre.gamma)
    
    pre.Q <- matrix(0, nrow=maxstate, ncol=maxstate)
    for (k in 1:nlambda) {
      pre.Q[states$from[k], states$to[k]] <- pre.lam[k]
    }
    diag(pre.Q) <- (-1)*apply(pre.Q, 1, sum)
    
    solve.cur.f <- function(x) {
      cur.Q <- matrix(0, nrow=maxstate, ncol=maxstate)
      for (k in 1:nlambda) {
        cur.Q[states$from[k], states$to[k]] <- x[k]
      }
      diag(cur.Q) <- (-1)*apply(cur.Q, 1, sum)
      
      cur.beta <- c(x[nlambda+1], x[nlambda+2])
      cur.gamma <- c(x[ntheta+1], x[ntheta+2])
      
      MAT <- lapply(SUBJ, MAT.f, locdata=indata, selectR=selectR, nlambda=nlambda, nbeta=nbeta, ngamma= ngamma, maxstate=maxstate, 
                    cur.Q=cur.Q, pre.Q=pre.Q, cur.beta=cur.beta, pre.beta=pre.beta, cur.gamma=cur.gamma, pre.gamma=pre.gamma)
      
      S <- Reduce('+', MAT)
      return(S)
    }
    
    cur.out <- nleqslv(c(pre.lam, pre.beta, pre.gamma), solve.cur.f, control = list(allowSingular=TRUE))
    cur.all <- cur.out$x
    
    cur.lam <- cur.all[1: nlambda]
    cur.beta <- cur.all[(nlambda+1): (nlambda+2)]
    cur.gamma <- cur.all[(ntheta+1): (ntheta+2)]
    
    tol <- max(abs( (cur.all - pre.all)/pre.all )[is.finite(abs( (cur.all - pre.all)/pre.all ))])
    
    # print(tol)
    if ( iter >= 500 ) { break }
  }
  
  # # covairance
  # cur.Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  # for (k in 1:nlambda) {
  #   cur.Q[states$from[k], states$to[k]] <- cur.lam[k]
  # }
  # diag(cur.Q) <- (-1)*apply(cur.Q, 1, sum)
  # 
  # MAT <- mclapply(SUBJ, MAT.f, locdata=indata, selectR=selectR, nlambda=nlambda, nbeta=nbeta, ngamma= ngamma, maxstate=maxstate, 
  # cur.Q=cur.Q, pre.Q=cur.Q, cur.beta=cur.beta, pre.beta=cur.beta, cur.gamma=cur.gamma, pre.gamma=cur.gamma, mc.cores=100)
  # 
  # M <- matrix(0, ntheta+ngamma, ntheta+ngamma)
  # for (i in 1:nSUBJ) {
  #   M <- M + t(MAT[[i]]) %*% MAT[[i]]
  # }
  # 
  # obs.I <- 1/M
  
  # cat("method = ", selectR, "\n")
  # cat("iter = ", iter, ": tol = ", tol, "\n")
  # cat("lam = ", round(cur.lam, 6), "\n")
  # cat("beta = ", round(cur.beta, 6), "\n")
  # cat("gamma = ", round(cur.gamma, 6), "\n")
  # cat("var.beta1 = ", round(obs.I[(nlambda+1), (nlambda+1)], 6), "\n")
  
  scoring <- NULL
  scoring$iter     <- iter
  scoring$tol      <- tol
  scoring$lam      <- as.vector(cur.lam)
  scoring$beta     <- as.vector(cur.beta)
  scoring$gamma    <- as.vector(cur.gamma)
  # scoring$cov      <- obs.I
  
  return(scoring)
  
  # cat(iter, tol, as.vector(cur.lam), as.vector(cur.beta), as.vector(cur.gamma),
  #     sep=" ", "\n", append=TRUE, file=file.out)
  
}


sim.likeli.f <- function(selectR, sim.num, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
  file.out <- paste("results/", "LIK_",selectR, "_NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                    "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "iter", "tol", "lam1"," lam2", "beta1", "beta2","gamma0", "gamma1",
      sep=" ", "\n", append=TRUE, file=file.out)
  
  
  parallel.sim.f <- function(selectR, nsim, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
    set.seed(nsim)
    
    stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
    stats$nsample2 <- nsample2
    
    simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                                  gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
    if (selectR == "SRS"){
      simdata <- p2data.SRS.f(indata=simdata_raw, nsample2=stats$nsample2)
    } else if(selectR == "OPT") {
      simdata <- p2data.OPT.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, gamma0=stats$gamma0, gamma1=stats$gamma1)
    } else if(selectR == "PROP"){
      simdata <- p2data.PROP.f(indata=simdata_raw, nsample2=stats$nsample2)
    }
    
    fit_S <- scoring.f(indata=simdata$event, selectR=selectR, nlambda=length(stats$lam), nbeta=2, ngamma=2,
                       inQ0=stats$Q, inbeta=c(stats$beta1, stats$beta2), 
                       ingamma= c(stats$gamma0, stats$gamma1) )
    cat(nsim, fit_S$iter, fit_S$tol, fit_S$lam, fit_S$beta, fit_S$gamma,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
  }
  
  mclapply(1:sim.num, parallel.sim.f, selectR=selectR, nsample=nsample, nsample2=nsample2, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, 
           p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina, mc.cores=cores.num)
  
  return()
}


# Likelihood-Direct Max -----------------------------------------------------------


evalLi.f <- function(indata, selectR, qmat, beta1, beta2, gamma0, gamma1) {
  indata$w <- indata$estop - indata$estart
  nlen <- length(indata$w)
  
  R <- indata[, selectR][1]
  x1 <- indata$x1[1]
  x2 <- indata$x2[1]
  
  qmat1 <- qmat*exp(beta1*x1+beta2*x2)
  qmat2 <- qmat*exp(beta1*1+beta2*x2)
  qmat3 <- qmat*exp(beta1*0+beta2*x2)
  
  # R=1
  PE <- 1
  for (j in 1:nlen) {
    Pmat <- MatrixExp(qmat1, t=indata$w[j], method="series")
    PE <- PE*Pmat[indata$from[j],indata$to[j]]
  }
  PX1 <- ifelse(x1==1, expit(gamma0+gamma1*x2), 1-expit(gamma0+gamma1*x2))
  PP.R1 <- PE*PX1
  
  # R=0
  PE1 <- PE2 <- 1
  for (j in 1:nlen) {
    Pmat1 <- MatrixExp(qmat2, t=indata$w[j], method="series")
    PE1 <- PE1*Pmat1[indata$from[j],indata$to[j]]
    
    Pmat2 <- MatrixExp(qmat3, t=indata$w[j], method="series")
    PE2 <- PE2*Pmat2[indata$from[j],indata$to[j]]
  }
  PX1.1 <- expit(gamma0+gamma1*x2)
  PX1.0 <- 1-expit(gamma0+gamma1*x2)
  PP.R0 <- PE1*PX1.1+PE2*PX1.0
  
  PP <- ifelse(R==1, PP.R1, PP.R0)
  return(PP)
}


logL.f <- function(p, maxstate, states, indata, selectR) {
  lam <- exp(p[1:2])
  beta1 <- p[3]
  beta2 <- p[4]
  gamma0 <- p[5]
  gamma1 <- p[6]
  
  qmat <- matrix(0, nrow=maxstate, ncol=maxstate)
  
  pos <- 0
  for (k in 1:nrow(states)) {
    pos <- pos + 1
    qmat[states$from[k], states$to[k]] <- lam[pos]
  }
  diag(qmat) <- (-1)*apply(qmat, 1, sum)
  
  L  <- lapply(sort(unique(as.character(indata$id))),
               function(id, indata, selectR, qmat, beta1, beta2, gamma0, gamma1) {
                 Li <- evalLi.f(indata=indata[as.character(indata$id) == as.character(id),], 
                                selectR=selectR,
                                qmat=qmat, 
                                beta1=beta1, beta2=beta2,
                                gamma0=gamma0, gamma1=gamma1)
                 return(Li) 
               }, indata=indata, selectR=selectR, qmat=qmat, 
               beta1=beta1, beta2=beta2,
               gamma0=gamma0, gamma1=gamma1)
  logL <- sum( log(unlist(L)) )
  return( (-1)*logL )
}


logL.msm.f <- function(qmat, indata, selectR, beta1, beta2, gamma0, gamma1) {
  maxstate <- nrow(qmat)
  
  p0 <- NULL
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (qmat[j,k] > 0) ) {
        p0 <- c(p0, qmat[j,k])
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  fit <- nlm(logL.f, p=c(log(p0), beta1, beta2, gamma0, gamma1), hessian=TRUE,
             maxstate=maxstate, states=states, indata=indata, 
             selectR=selectR,
             steptol=1e-06, stepmax= 5)
  
  fit$cov <- ginv(fit$hessian)
  return(fit)
}



asymp.likeli.f <- function(fit, indata, inQ0, selectR) {
  event <- indata
  
  maxstate <- nrow(inQ0)
  
  lam <- exp(fit$estimate[1:2])
  beta <- fit$estimate[3:4]
  gamma <- fit$estimate[5:6]
  
  nlambda <- length(lam)
  nbeta <- length(beta)
  ngamma <- length(gamma)
  
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- lam[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  ntheta <- nlambda + nbeta
  
  SUBJ  <- unique(as.character(event$id))
  nSUBJ <- length(SUBJ)
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, ngamma, ntheta, maxstate, inQ0, beta, gamma) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                  x1=x1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    # pseudo
    outQ.x1.1 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                       x1=1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx.x1.1 <- outQ.x1.1$outQx
    dQx.x1.1 <- outQ.x1.1$outdQx
    
    outQ.x1.0 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                       x1=0, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx.x1.0 <- outQ.x1.0$outQx
    dQx.x1.0 <- outQ.x1.0$outdQx
    
    nlen <- length(locdata$w)
    PD.x1.1 <- 1 # P(D | X1=1, X2)
    PD.x1.0 <- 1 # P(D | X1=0, X2)
    
    SSi1 <- SSi1.x1.1 <- SSi1.x1.0 <- matrix(0, nrow=1, ncol=ntheta)
    SSi2 <- SSi2.x1.1 <- SSi2.x1.0 <- matrix(0, nrow=1, ncol=ngamma)
    
    IIi1 <- IIi1.x1.1 <- IIi1.x1.0 <- matrix(0, nrow=ntheta, ncol=ntheta)
    IIi2 <- matrix(0, nrow=ngamma, ncol=ngamma)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      Pxmat.x1.1 <- round(MatrixExp(Qx.x1.1, t=locdata$w[i], method="series"), 16)
      Pxmat.x1.0 <- round(MatrixExp(Qx.x1.0, t=locdata$w[i], method="series"), 16)
      
      PD.x1.1 <- PD.x1.1*Pxmat.x1.1[locdata$from[i], locdata$to[i]]
      PD.x1.0 <- PD.x1.0*Pxmat.x1.0[locdata$from[i], locdata$to[i]]
      
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        dPxmatu.x1.1 <- dPxt.f(u, locdata$w[i], Qx.x1.1, dQx.x1.1, maxstate)
        dPxmatu.x1.0 <- dPxt.f(u, locdata$w[i], Qx.x1.0, dQx.x1.0, maxstate)
        
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        SSi1.x1.1[1,u] <- SSi1.x1.1[1,u] + dPxmatu.x1.1[locdata$from[i],locdata$to[i]] / Pxmat.x1.1[locdata$from[i],locdata$to[i]]
        SSi1.x1.0[1,u] <- SSi1.x1.0[1,u] + dPxmatu.x1.0[locdata$from[i],locdata$to[i]] / Pxmat.x1.0[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          dPxmatj.x1.1 <- dPxt.f(j, locdata$w[i], Qx.x1.1, dQx.x1.1, maxstate)
          dPxmatj.x1.0 <- dPxt.f(j, locdata$w[i], Qx.x1.0, dQx.x1.0, maxstate)
          # MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              IIi1[u,j] <- IIi1[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
              IIi1.x1.1[u,j] <- IIi1.x1.1[u,j] + ( (1/Pxmat.x1.1[locdata$from[i],j2])*dPxmatu.x1.1[locdata$from[i],j2]*dPxmatj.x1.1[locdata$from[i],j2] )
              IIi1.x1.0[u,j] <- IIi1.x1.0[u,j] + ( (1/Pxmat.x1.0[locdata$from[i],j2])*dPxmatu.x1.0[locdata$from[i],j2]*dPxmatj.x1.0[locdata$from[i],j2] )
              
            }
          }
          
        }
        
      }
    }
    
    
    SSi2[1,1] <- x1-expit(gamma[1]+gamma[2]*x2)
    SSi2[1,2] <- (x1-expit(gamma[1]+gamma[2]*x2))*x2
    SSi2.x1.1[1,1] <- 1-expit(gamma[1]+gamma[2]*x2)
    SSi2.x1.1[1,2] <- (1-expit(gamma[1]+gamma[2]*x2))*x2
    SSi2.x1.0[1,1] <- 0-expit(gamma[1]+gamma[2]*x2)
    SSi2.x1.0[1,2] <- (0-expit(gamma[1]+gamma[2]*x2))*x2
    
    IIi2[1,1] <- expit(gamma[1]+gamma[2]*x2)*( 1/(1+exp(gamma[1]+gamma[2]*x2)) )
    IIi2[1,2] <- IIi2[2,1] <- expit(gamma[1]+gamma[2]*x2)*( 1/(1+exp(gamma[1]+gamma[2]*x2)) )*x2
    IIi2[2,2] <- expit(gamma[1]+gamma[2]*x2)*( 1/(1+exp(gamma[1]+gamma[2]*x2)) )*x2^2
    
    eta <- ( PD.x1.1*expit(gamma[1]+gamma[2]*x2) )/( PD.x1.1*expit(gamma[1]+gamma[2]*x2) + PD.x1.0*(1-expit(gamma[1]+gamma[2]*x2)) )
    
    Si1 <- R*SSi1 + (1-R)*( SSi1.x1.1*eta + SSi1.x1.0*(1-eta) )
    Si2 <- R*SSi2 + (1-R)*( SSi2.x1.1*eta + SSi2.x1.0*(1-eta) )
    Si <- cbind(Si1, Si2)
    
    Ii1 <- R*IIi1 + (1-R)*( IIi1.x1.1*eta +  IIi1.x1.0*(1-eta) - t(SSi1.x1.1)%*%SSi1.x1.1*eta - t(SSi1.x1.0)%*%SSi1.x1.0*(1-eta) +
                              t((SSi1.x1.1*eta + SSi1.x1.0*(1-eta)))%*%(SSi1.x1.1*eta + SSi1.x1.0*(1-eta)) )
    Ii2 <- R*IIi2 + (1-R)*( IIi2 -  t(SSi2.x1.1)%*%SSi2.x1.1*eta - t(SSi2.x1.0)%*%SSi2.x1.0*(1-eta) +
                              t(( SSi2.x1.1*eta + SSi2.x1.0*(1-eta) ))%*%( SSi2.x1.1*eta + SSi2.x1.0*(1-eta) ) )
    
    outI <- matrix(0, nrow = ntheta+ngamma, ncol = ntheta+ngamma)
    outI[1:ntheta, 1:ntheta] <- Ii1
    outI[(ntheta+1):(ntheta+ngamma), (ntheta+1):(ntheta+ngamma)] <- Ii2
    
    out <- rbind(Si, outI)
    
    return(out)
    #return(Si)
  }
  
  MAT <- lapply(SUBJ, MAT.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, ngamma=ngamma, ntheta=ntheta, maxstate=maxstate,
                inQ0=Q, beta=beta, gamma=gamma)
  
  SxSi <- lapply(MAT, function(x){ x[1, c(1:(ntheta+ngamma))]%*% t(x[1, c(1:(ntheta+ngamma))]) })
  SE1 <- sqrt((ginv(Reduce("+", SxSi)/nSUBJ)/nSUBJ)[3,3]) # score*score
  
  Ii <- lapply(MAT, function(x){ x[c((1:(ntheta+ngamma))+1), c(1:(ntheta+ngamma))] })
  SE2 <- sqrt(ginv(Reduce("+", Ii))[3,3]) # Louise method
  
  return(c(SE1, SE2))
}


sim.likeli.dir.f <- function(selectR, sim.num, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
  file.out <- paste("results/LIK/", "LIK_DIR_",selectR, "_NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                    "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "lam1"," lam2", "beta1", "beta2","gamma0", "gamma1", "asymp.se1", "asymp.se2",
      sep=" ", "\n", append=TRUE, file=file.out)
  
  
  parallel.sim.f <- function(nsim, selectR, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
    set.seed(nsim)
    
    stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
    stats$nsample2 <- nsample2
    
    simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                                  gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
    if (selectR == "SRS"){
      simdata <- p2data.SRS.f(indata=simdata_raw, nsample2=stats$nsample2)
    } else if(selectR == "OPT") {
      simdata <- p2data.OPT.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, gamma0=stats$gamma0, gamma1=stats$gamma1)
    } else if(selectR == "PROP"){
      simdata <- p2data.PROP.f(indata=simdata_raw, nsample2=stats$nsample2)
    } else if(selectR == "BAL") {
      simdata <- p2data.BAL.f(indata=simdata_raw, nsample2=stats$nsample2)
    }
    
    
    fitL <- logL.msm.f(indata=simdata$event, qmat=stats$Q, selectR=selectR, beta1=stats$beta1, beta2=stats$beta2, 
                       gamma0=stats$gamma0, gamma1=stats$gamma1)
    
    asymp.se <- asymp.likeli.f(fit=fitL, indata=simdata$event, inQ0=stats$Q, selectR=selectR)
    cat(nsim, exp(fitL$estimate[1:2]), fitL$estimate[3:6], asymp.se,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
  }
  
  mclapply(1:sim.num, parallel.sim.f, selectR=selectR, nsample=nsample, nsample2=nsample2, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, 
           p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina, mc.cores=cores.num)
  
  return()
}

# Con-Likelihood-Direct  -----------------------------------------------------------

pi.f <- function(indata, selectR, inQ0) { # indata=simdata, inQ0=stats$Q
  rawdata <- indata$rawdata
  event <- indata$event
  
  maxstate <- nrow(inQ0)
  
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  dat <- rawdata %>% group_by(id) %>% mutate(Z.am = max(states),
                                             am = max(times)) %>% dplyr::select(id, x2, am, Z.am, selectR) %>% slice(1)
  dat <- as.data.frame(dat)
  dat[, selectR] <- as.factor(dat[, selectR])
  dat[, "Z.am"] <- as.factor(dat[, "Z.am"])
  dat[, "x2"] <- as.factor(dat[, "x2"])
  
  fit.pi <- glm(dat[, selectR] ~ x2*Z.am, data = dat, family = "binomial")
  # fit.pi$coefficients
  # summary(fit.pi)
  dat$pi <- predict(fit.pi, newdata = dat, type = "response")
  dat <- dat %>% dplyr::select(id, am, Z.am, pi)
  event <- left_join(event, dat, by="id")
  # event$x2 <- as.factor(event$x2)
  
  out <- NULL
  out$event <- event
  out$fit.pi <- fit.pi
  
  return(out)
}

evalConLi.f <- function(indata, fit.pi, selectR, qmat, beta1, beta2) { # indata=event for i
  indata$w <- indata$estop - indata$estart
  nlen <- length(indata$w)
  
  R <- indata[, selectR][1]
  x1 <- indata$x1[1]
  x2 <- indata$x2[1]
  pi <- indata$pi[1]
  
  qmat1 <- qmat*exp(beta1*x1+beta2*x2)
  
  # R=1
  PE <- 1
  for (j in 1:nlen) {
    Pmat <- MatrixExp(qmat1, t=indata$w[j], method="series")
    PE <- PE*Pmat[indata$from[j],indata$to[j]]
  }
  nume <- PE*pi
  
  pi.state1 <- predict(fit.pi, newdata = data.frame(Z.am=as.factor(1), x2=as.factor(x2) ), type = "response")
  pi.state2 <- predict(fit.pi, newdata = data.frame(Z.am=as.factor(2), x2=as.factor(x2) ), type = "response")
  pi.state3 <- predict(fit.pi, newdata = data.frame(Z.am=as.factor(3), x2=as.factor(x2) ), type = "response")
  
  Pmat.am <- MatrixExp(qmat1, t=indata$am[1], method="series")
  
  deno <- pi.state1*Pmat.am[1,1] + pi.state2*Pmat.am[1,2] + pi.state3*Pmat.am[1,3]
  
  PP <- nume/deno
  return(PP)
}


logConL.f <- function(p, maxstate, states, indata, selectR, fit.pi) {
  lam <- exp(p[1:2])
  beta1 <- p[3]
  beta2 <- p[4]
  
  qmat <- matrix(0, nrow=maxstate, ncol=maxstate)
  
  pos <- 0
  for (k in 1:nrow(states)) {
    pos <- pos + 1
    qmat[states$from[k], states$to[k]] <- lam[pos]
  }
  diag(qmat) <- (-1)*apply(qmat, 1, sum)
  
  SUBJ  <- unique(as.character(indata$id[indata[, selectR]==1]))
  
  L  <- lapply(sort(SUBJ),
               function(id, indata, fit.pi, selectR, qmat, beta1, beta2) {
                 Li <- evalConLi.f(indata=indata[as.character(indata$id) == as.character(id),], 
                                   fit.pi=fit.pi,
                                   selectR=selectR,
                                   qmat=qmat, 
                                   beta1=beta1, beta2=beta2)
                 return(Li) 
               }, indata=indata, fit.pi=fit.pi, selectR=selectR, qmat=qmat, 
               beta1=beta1, beta2=beta2)
  logL <- sum( log(unlist(L)) )
  return( (-1)*logL )
}


logConL.msm.f <- function(qmat, indata, selectR, beta1, beta2) { # indata=simdata, inQ0=stats$Q
  out <- pi.f(indata=indata, selectR=selectR, inQ0=qmat)
  event <- out$event
  fit.pi <- out$fit.pi
  
  maxstate <- nrow(qmat)
  
  p0 <- NULL
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (qmat[j,k] > 0) ) {
        p0 <- c(p0, qmat[j,k])
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  fit <- nlm(logConL.f, p=c(log(p0), beta1, beta2), hessian=TRUE,
             maxstate=maxstate, states=states, indata=event, 
             selectR=selectR, fit.pi=fit.pi,
             steptol=1e-06, stepmax= 5)
  
  # fit$cov <- ginv(fit$hessian)
  return(fit)
}


# function(indata, fit.pi, selectR, qmat, beta1, beta2)

asymp.conlikeli.f <- function(fit, indata, inQ0, selectR) {
  out <- pi.f(indata=indata, selectR=selectR, inQ0=inQ0)
  event <- out$event
  fit.pi <- out$fit.pi
  
  maxstate <- nrow(inQ0)
  
  lam <- exp(fit$estimate[1:2])
  beta <- fit$estimate[3:4]
  
  nlambda <- length(lam)
  nbeta <- length(beta)
  
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- lam[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  ntheta <- nlambda + nbeta
  
  SUBJ  <- unique(as.character(event$id[event[, selectR]==1])) 
  nSUBJ <- length(SUBJ)
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, ntheta, maxstate, inQ0, beta, fit.pi) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    pi <- locdata$pi[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                  x1=x1, x2=x2, inQ=Q, inbeta1=beta[1], inbeta2=beta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
      }
    }
    
    
    
    pi.state1 <- predict(fit.pi, newdata = data.frame(Z.am=as.factor(1), x2=as.factor(x2) ), type = "response")
    pi.state2 <- predict(fit.pi, newdata = data.frame(Z.am=as.factor(2), x2=as.factor(x2) ), type = "response")
    pi.state3 <- predict(fit.pi, newdata = data.frame(Z.am=as.factor(3), x2=as.factor(x2) ), type = "response")
    
    Pmat.am <- MatrixExp(Qx, t=locdata$am[1], method="series")
    deno <- pi.state1*Pmat.am[1,1] + pi.state2*Pmat.am[1,2] + pi.state3*Pmat.am[1,3]
    
    SSi2 <- matrix(0, nrow=1, ncol=ntheta)
    for (u in 1:ntheta) {
      dPxmatu <- dPxt.f(u, locdata$am[1], Qx, dQx, maxstate)
      SSi2[1, u] <- SSi2[1, u] + pi.state1*dPxmatu[1, 1] + pi.state2*dPxmatu[1, 2] + pi.state3*dPxmatu[1, 3]
    }
    
    Sc <- R*(SSi1 - SSi2/deno)
    
    
    return(Sc)
  }
  
  MAT <- lapply(SUBJ, MAT.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta,  ntheta=ntheta, maxstate=maxstate,
                inQ0=Q, beta=beta, fit.pi=fit.pi)
  
  SxSi <- lapply(MAT, function(x){ x[1, 1:4]%*% t(x[1, 1:4]) })
  SE1 <- sqrt((ginv(Reduce("+", SxSi)/nSUBJ)/nSUBJ)[3,3]) # score*score
  
  
  return(SE1)
}



sim.conlikeli.dir.f <- function(selectR, sim.num, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
  file.out <- paste("results/CONLIK/", "CONLIK_DIR_",selectR, "_NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                    "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "lam1"," lam2", "beta1", "beta2", "se",
      sep=" ", "\n", append=TRUE, file=file.out)
  
  
  parallel.sim.f <- function(nsim, selectR, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
    set.seed(nsim)
    
    stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
    stats$nsample2 <- nsample2
    
    simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                                  gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
    if (selectR == "SRS"){
      simdata <- p2data.SRS.f(indata=simdata_raw, nsample2=stats$nsample2)
    } else if(selectR == "OPT") {
      simdata <- p2data.OPT.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, gamma0=stats$gamma0, gamma1=stats$gamma1)
    } else if(selectR == "PROP"){
      simdata <- p2data.PROP.f(indata=simdata_raw, nsample2=stats$nsample2)
    } else if(selectR == "BAL") {
      simdata <- p2data.BAL.f(indata=simdata_raw, nsample2=stats$nsample2)
    } else if(selectR == "IPW.OPT.TRUE6") {
      simdata <- p2data.IPW.true6.f(indata=simdata_raw, stats=stats, nsample2=stats$nsample2, inQ=stats$Q, nlambda=2, nbeta=2, gamma0=stats$gamma0, gamma1=stats$gamma1)
    }
    
    fitL <- logConL.msm.f(indata=simdata, qmat=stats$Q, selectR=selectR, beta1=stats$beta1, beta2=stats$beta2)
    asymp.se <- asymp.conlikeli.f(fit=fitL, indata=simdata, inQ0=stats$Q, selectR=selectR)
    
    cat(nsim, exp(fitL$estimate[1:2]), fitL$estimate[3:4], asymp.se,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
  }
  
  mclapply(1:sim.num, parallel.sim.f, selectR=selectR, nsample=nsample, nsample2=nsample2, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, 
           p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina, mc.cores=cores.num)
  
  return()
}




# IPW ---------------------------------------------------------------------
ipw.f <- function(indata, selectR, nlambda, nbeta, inQ0, inbeta) {
  rawdata <- indata$rawdata
  event <- indata$event
  
  maxstate <- nrow(inQ0)
  ntheta <- nlambda + nbeta
  
  lam <- NULL
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        lam <- c(lam, inQ0[j,k])
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }  
  
  dat <- rawdata %>% group_by(id) %>% mutate(Z.am = max(states),
                                             am = max(times)) %>% dplyr::select(id, x2, am, Z.am, selectR) %>% slice(1)
  dat <- as.data.frame(dat)
  dat[, selectR] <- as.factor(dat[, selectR])
  dat[, "Z.am"] <- as.factor(dat[, "Z.am"])
  dat[, "x2"] <- as.factor(dat[, "x2"])
  
  fit.pi <- glm(dat[, selectR] ~ x2*Z.am, data = dat, family = "binomial")
  # summary(fit.pi)
  dat$pi <- predict(fit.pi, newdata = dat, type = "response")
  dat <- dat %>% dplyr::select(id, pi)
  event <- left_join(event, dat, by="id")
  
  SUBJ  <- unique(as.character(event$id[event[, selectR]==1])) 
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    pi <- locdata$pi[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
      }
    }
    
    Si <- R/pi*SSi1
    
    return(Si)
  }
  
  
  solve.f <- function(x) {
    Q <- matrix(0, nrow=maxstate, ncol=maxstate)
    for (k in 1:nlambda) {
      Q[states$from[k], states$to[k]] <- x[k]
    }
    diag(Q) <- (-1)*apply(Q, 1, sum)
    
    beta <- c(x[nlambda+1], x[nlambda+2])
    
    MAT <- lapply(SUBJ, MAT.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, maxstate=maxstate, 
                  Q=Q, inbeta=beta)
    
    S <- Reduce('+', MAT)
    return(S)
  }
  
  out.theta <- nleqslv(c(lam, inbeta), solve.f, control = list(allowSingular=TRUE))
  
  # asymptotic variance
  MAT2.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    pi <- locdata$pi[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    MMi <- EMMi <- matrix(0, nrow=ntheta, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              EMMi[u,j] <- EMMi[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
            }
          }
        }
        
      }
    }
    
    Si <- R/pi*SSi1
    Mi <- R/pi*MMi
    EMi <- R/pi*EMMi
    outMAT <- rbind(Si, Mi, EMi)
    return(outMAT)
  }
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- out.theta$x[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  beta <- c(out.theta$x[nlambda+1], out.theta$x[nlambda+2])
  
  MAT2 <- lapply(SUBJ, MAT2.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, maxstate=maxstate, 
                 Q=Q, inbeta=beta)
  nSUBJ <- length(unique(as.character(event$id)))
  
  
  SxSi <- lapply(MAT2, function(x){ x[1, c(1:(ntheta))] %*% t(x[1, c(1:(ntheta))]) })
  B <- Reduce("+", SxSi)/nSUBJ
  
  Mi <- lapply(MAT2, function(x){ x[c(2:5), c(1:(ntheta))] })
  EMi <- lapply(MAT2, function(x){ x[c(6:9), c(1:(ntheta))] })
  
  A.M <- Reduce("+", Mi)/nSUBJ
  A.EM <- Reduce("+", EMi)/nSUBJ
  
  SE.M <- sqrt( (ginv(A.M) %*% B %*% t(ginv(A.M))/nSUBJ)[3,3] )
  SE.EM <- sqrt( (ginv(A.EM) %*% B %*% t(ginv(A.EM))/nSUBJ)[3,3] )
  
  out.lam <- out.theta$x[1: nlambda]
  out.beta <- out.theta$x[(nlambda+1): (nlambda+2)]
  
  scoring <- NULL
  scoring$conv      <- out.theta$termcd
  scoring$fvec      <- out.theta$fvec
  scoring$lam      <- as.vector(out.lam)
  scoring$beta     <- as.vector(out.beta)
  scoring$se1      <- SE.EM
  scoring$se2      <- SE.M
  return(scoring)
  
}



sim.ipw.f <- function(selectR, sim.num, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
  file.out <- paste("results/IPW/", "IPW_",selectR, "_NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                    "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "conv", "fvec.lam1", "fvec.lam2", "fvec.beta1", "fvec.beta2", "lam1"," lam2", "beta1", "beta2", "se1", "se2",
      sep=" ", "\n", append=TRUE, file=file.out)
  
  
  parallel.sim.f <- function(selectR, nsim, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
    set.seed(nsim)
    
    stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
    stats$nsample2 <- nsample2
    
    simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                                  gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
    if (selectR == "SRS"){
      simdata <- p2data.SRS.f(indata=simdata_raw, nsample2=stats$nsample2)
    } else if(selectR == "IPW.OPT.TRUE") {
      simdata <- p2data.IPW.true.f(indata=simdata_raw, stats=stats, nsample2=stats$nsample2, inQ=stats$Q, nlambda=2, nbeta=2, gamma0=stats$gamma0, gamma1=stats$gamma1)
    } else if(selectR == "IPW.OPT.TRUE6") {
      simdata <- p2data.IPW.true6.f(indata=simdata_raw, stats=stats, nsample2=stats$nsample2, inQ=stats$Q, nlambda=2, nbeta=2, gamma0=stats$gamma0, gamma1=stats$gamma1)
    } else if(selectR == "IPW.OPT"){
      simdata <- p2data.IPW.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, nlambda=2, nbeta=2, gamma0=stats$gamma0, gamma1=stats$gamma1)
    } else if(selectR == "BAL") {
      simdata <- p2data.BAL.f(indata=simdata_raw, nsample2=stats$nsample2)
    }

    fit_S <- ipw.f(indata=simdata, selectR=selectR, nlambda=length(stats$lam), nbeta=2, 
                   inQ0=stats$Q, inbeta=c(stats$beta1, stats$beta2))
    
    cat(nsim, fit_S$conv, fit_S$fvec, fit_S$lam, fit_S$beta, fit_S$se1, fit_S$se2,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
  }
  
  mclapply(1:sim.num, parallel.sim.f, selectR=selectR, nsample=nsample, nsample2=nsample2, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, 
           p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina, mc.cores=cores.num)
  
  return()
}


# GR ---------------------------------------------------------------------
check.pred.f <- function(indata, selectR, inQ0, nsim, inbeta, sim.data) {
  beta <- c(sim.data$beta1[sim.data$nsim == nsim], sim.data$beta2[sim.data$nsim == nsim])
  gamma <- c(sim.data$gamma0[sim.data$nsim == nsim], sim.data$gamma1[sim.data$nsim == nsim])
  lam <- c(sim.data$lam1[sim.data$nsim == nsim], sim.data$lam2[sim.data$nsim == nsim])
  
  maxstate <- nrow(inQ0)
  
  event <- indata$event
  rawdata <- indata$rawdata
  
  event.non.phase2 <- event %>% filter(!!sym(selectR) == 0)
  
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }
  
  
  nlambda <- length(lam)
  nbeta <- length(beta)
  ngamma <- length(gamma)
  ntheta <- nlambda + nbeta
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- lam[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  
  SUBJ  <- unique(as.character(event.non.phase2$id))
  nSUBJ <- length(SUBJ)
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, ngamma, maxstate, inQ0, beta, gamma) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                  x1=x1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    # pseudo
    outQ.x1.1 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                       x1=1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx.x1.1 <- outQ.x1.1$outQx
    dQx.x1.1 <- outQ.x1.1$outdQx
    
    outQ.x1.0 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                       x1=0, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx.x1.0 <- outQ.x1.0$outQx
    dQx.x1.0 <- outQ.x1.0$outdQx
    
    nlen <- length(locdata$w)
    PD.x1.1 <- 1 # P(D | X1=1, X2)
    PD.x1.0 <- 1 # P(D | X1=0, X2)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      Pxmat.x1.1 <- round(MatrixExp(Qx.x1.1, t=locdata$w[i], method="series"), 16)
      Pxmat.x1.0 <- round(MatrixExp(Qx.x1.0, t=locdata$w[i], method="series"), 16)
      
      PD.x1.1 <- PD.x1.1*Pxmat.x1.1[locdata$from[i], locdata$to[i]]
      PD.x1.0 <- PD.x1.0*Pxmat.x1.0[locdata$from[i], locdata$to[i]]
      
    }
    
    
    eta <- ( PD.x1.1*expit(gamma[1]+gamma[2]*x2) )/( PD.x1.1*expit(gamma[1]+gamma[2]*x2) + PD.x1.0*(1-expit(gamma[1]+gamma[2]*x2)) )
    
    
    
    return(eta)
  }
  
  eta <- lapply(SUBJ, MAT.f, locdata=event.non.phase2, selectR=selectR, nlambda=nlambda, nbeta=nbeta, ngamma=ngamma, maxstate=maxstate,
                inQ0=Q, beta=beta, gamma=gamma)
  
  eta <- unlist(eta)
  # summary(eta)
  
  dat.temp <- data.frame(id=SUBJ, eta=eta, x1.pred=rbinom(length(eta), 1, eta))
  
  
  out <- pi.f(indata=indata, selectR=selectR, inQ0=inQ0)
  event <- out$event
  fit.pi <- out$fit.pi
  rawdata <- indata$rawdata
  
  nsample <- length(unique(event$id))
  
  dat <- event %>% dplyr::select(id, x1, x2, selectR, am, Z.am, pi) %>% group_by(id) %>% slice(1)
  dat <- as.data.frame(dat)
  dat[, selectR] <- as.factor(dat[, selectR])
  dat[, "Z.am"] <- as.factor(dat[, "Z.am"])
  dat[, "x1"] <- as.factor(dat[, "x1"])
  dat[, "x2"] <- as.factor(dat[, "x2"])
  
  dat.phase2 <- dat %>% filter(!!sym(selectR) == 1) %>% mutate(x1.pred=x1)
  dat.non.phase2 <- dat %>% filter(!!sym(selectR) == 0)
  
  dat.non.phase2$x1.pred <- factor(dat.temp$x1.pred, levels = c(0, 1))
  dat.pred <- rbind(dat.phase2, dat.non.phase2)
  dat.pred <- dat.pred[order(dat.pred$id), ]
  
  dat.pred <- dat.pred %>% dplyr::select(id, x1.pred)
  
  rawdata <- left_join(rawdata, dat.pred, by = "id")
  event <- left_join(event, dat.pred, by = "id")
  
  rawdata[, "x1.pred"] <- as.numeric(as.character(rawdata[, "x1.pred"]))
  event[, "x1.pred"] <- as.numeric(as.character(event[, "x1.pred"]))
  
  fitm <- msm(states ~ times, subject=id, data=rawdata,
              covariates = ~ x1.pred + x2,
              constraint = list(x1.pred = c(1,1), x2 = c(1,1)),
              qmatrix=inQ0, obstype=1, opt.method="fisher",
              gen.inits=TRUE, center = FALSE)
  # fitm$estimates.t
  # I.inverse <- fitm$covmat[3,3]
  
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1.pred[1] # x1 <- locdata$x1.pred[1]
    x2 <- locdata$x2[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    ntheta <- nlambda + nbeta
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    MMi <- EMMi <- matrix(0, nrow=ntheta, ncol=ntheta)
    
    # SSi2 <- matrix(0, nrow=1, ncol=ngamma)
    # IIi2 <- matrix(0, nrow=ngamma, ncol=ngamma)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              EMMi[u,j] <- EMMi[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
            }
          }
        }
      }
    }
    
    
    outMAT <- rbind(SSi1, EMMi, MMi)
    return(outMAT)
    # return(SSi1[3])
  }
  
  SUBJ  <- unique(as.character(event$id))
  
  lam <- NULL
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        lam <- c(lam, inQ0[j,k])
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }
  
  Q.impute <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q.impute[states$from[k], states$to[k]] <- fitm$estimates.t[k]
  }
  diag(Q.impute) <- (-1)*apply(Q.impute, 1, sum)
  
  beta.impute <- fitm$estimates.t[c(3,5)]
  
  # gamma <- fitm2$coefficients
  
  MAT <- lapply(SUBJ, MAT.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, maxstate=maxstate,
                Q=Q.impute, inbeta=beta.impute)
  
  efficient.info <- lapply(MAT, function(x){ (ginv(x[c(2:5), c(1:4)]) %*% x[1,])[3] })
  # efficient.info <- lapply(MAT, function(x){ (ginv(x[c(6:9), c(1:4)]) %*% x[1,])[3] })
  # efficient.info <- lapply(MAT, function(x){ x[1,3] })
  
  efficient.info <- unlist(efficient.info)
  # efficient.info <- trim(efficient.info, at = 0.995)
  # summary(efficient.info)
  # efficient.info <- Reduce('+', MAT)
  
  
  dat <- event %>% dplyr::select(id, x1, x2, x1.pred, selectR, am, Z.am, pi) %>% group_by(id) %>% slice(1)
  dat$efficient.info <- efficient.info
  
  Xs <- dat$efficient.info[dat[, selectR]==1]       # sample matrix of calibration variables
  piks <- dat$pi[dat[,selectR]==1]              # sample inclusion probabilities
  # vector of population totals of the calibration variables
  total <- sum(efficient.info)
  # summary(efficient.info)
  # computation of the g-weights
  g2 <- calib(Xs, d=1/piks, total, method="logit", bounds=c(0.5, 1.5), max_iter=500)
  summary(g2)
  
  check <- length(g2) == length(unique(event$id))-nSUBJ
  
  out <- NULL
  out$check <- check
  out$g2 <- g2
  return(out)
}



calib.f <- function(indata, selectR, inQ0, nsim, inbeta, sim.data, g2) {
  
  beta <- c(sim.data$beta1[sim.data$nsim == nsim], sim.data$beta2[sim.data$nsim == nsim])
  gamma <- c(sim.data$gamma0[sim.data$nsim == nsim], sim.data$gamma1[sim.data$nsim == nsim])
  lam <- c(sim.data$lam1[sim.data$nsim == nsim], sim.data$lam2[sim.data$nsim == nsim])
  
  maxstate <- nrow(inQ0)
  
  event <- indata$event
  rawdata <- indata$rawdata
  
  event.non.phase2 <- event %>% filter(!!sym(selectR) == 0)
  
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }
  
  
  nlambda <- length(lam)
  nbeta <- length(beta)
  ngamma <- length(gamma)
  ntheta <- nlambda + nbeta
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- lam[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  
  
  SUBJ  <- unique(as.character(event.non.phase2$id))
  nSUBJ <- length(SUBJ)
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, ngamma, maxstate, inQ0, beta, gamma) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                  x1=x1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    # pseudo
    outQ.x1.1 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                       x1=1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx.x1.1 <- outQ.x1.1$outQx
    dQx.x1.1 <- outQ.x1.1$outdQx
    
    outQ.x1.0 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                       x1=0, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
    Qx.x1.0 <- outQ.x1.0$outQx
    dQx.x1.0 <- outQ.x1.0$outdQx
    
    nlen <- length(locdata$w)
    PD.x1.1 <- 1 # P(D | X1=1, X2)
    PD.x1.0 <- 1 # P(D | X1=0, X2)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      Pxmat.x1.1 <- round(MatrixExp(Qx.x1.1, t=locdata$w[i], method="series"), 16)
      Pxmat.x1.0 <- round(MatrixExp(Qx.x1.0, t=locdata$w[i], method="series"), 16)
      
      PD.x1.1 <- PD.x1.1*Pxmat.x1.1[locdata$from[i], locdata$to[i]]
      PD.x1.0 <- PD.x1.0*Pxmat.x1.0[locdata$from[i], locdata$to[i]]
      
    }
    
    
    eta <- ( PD.x1.1*expit(gamma[1]+gamma[2]*x2) )/( PD.x1.1*expit(gamma[1]+gamma[2]*x2) + PD.x1.0*(1-expit(gamma[1]+gamma[2]*x2)) )
    
    
    
    return(eta)
  }
  
  eta <- lapply(SUBJ, MAT.f, locdata=event.non.phase2, selectR=selectR, nlambda=nlambda, nbeta=nbeta, ngamma=ngamma, maxstate=maxstate,
                inQ0=Q, beta=beta, gamma=gamma)
  
  eta <- unlist(eta)
  # summary(eta)
  
  dat.temp <- data.frame(id=SUBJ, eta=eta, x1.pred=rbinom(length(eta), 1, eta))
  
  
  out <- pi.f(indata=indata, selectR=selectR, inQ0=inQ0)
  event <- out$event
  fit.pi <- out$fit.pi
  rawdata <- indata$rawdata
  
  nsample <- length(unique(event$id))
  
  dat <- event %>% dplyr::select(id, x1, x2, selectR, am, Z.am, pi) %>% group_by(id) %>% slice(1)
  dat <- as.data.frame(dat)
  dat[, selectR] <- as.factor(dat[, selectR])
  dat[, "Z.am"] <- as.factor(dat[, "Z.am"])
  dat[, "x1"] <- as.factor(dat[, "x1"])
  dat[, "x2"] <- as.factor(dat[, "x2"])
  
  dat.phase2 <- dat %>% filter(!!sym(selectR) == 1) %>% mutate(x1.pred=x1)
  dat.non.phase2 <- dat %>% filter(!!sym(selectR) == 0)
  
  dat.non.phase2$x1.pred <- factor(dat.temp$x1.pred, levels = c(0, 1))
  dat.pred <- rbind(dat.phase2, dat.non.phase2)
  dat.pred <- dat.pred[order(dat.pred$id), ]
  
  dat.pred <- dat.pred %>% dplyr::select(id, x1.pred)
  
  rawdata <- left_join(rawdata, dat.pred, by = "id")
  event <- left_join(event, dat.pred, by = "id")
  
  rawdata[, "x1.pred"] <- as.numeric(as.character(rawdata[, "x1.pred"]))
  event[, "x1.pred"] <- as.numeric(as.character(event[, "x1.pred"]))
  
  fitm <- msm(states ~ times, subject=id, data=rawdata,
              covariates = ~ x1.pred + x2,
              constraint = list(x1.pred = c(1,1), x2 = c(1,1)),
              qmatrix=inQ0, obstype=1, opt.method="fisher",
              gen.inits=TRUE, center = FALSE)
  # fitm$estimates.t
  # I.inverse <- fitm$covmat[3,3]
  
  
  MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1.pred[1] # x1 <- locdata$x1.pred[1]
    x2 <- locdata$x2[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    ntheta <- nlambda + nbeta
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    MMi <- EMMi <- matrix(0, nrow=ntheta, ncol=ntheta)
    
    # SSi2 <- matrix(0, nrow=1, ncol=ngamma)
    # IIi2 <- matrix(0, nrow=ngamma, ncol=ngamma)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              EMMi[u,j] <- EMMi[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
            }
          }
        }
      }
    }
    
    
    outMAT <- rbind(SSi1, EMMi, MMi)
    return(outMAT)
    # return(SSi1[3])
  }
  
  SUBJ  <- unique(as.character(event$id))
  
  lam <- NULL
  states <- NULL
  for (j in 1:maxstate) {
    for (k in 1:maxstate) {
      if ( (j != k) && (inQ0[j,k] > 0) ) {
        lam <- c(lam, inQ0[j,k])
        states <- rbind(states, data.frame(from=j, to=k))
      }
    }
  }
  
  Q.impute <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q.impute[states$from[k], states$to[k]] <- fitm$estimates.t[k]
  }
  diag(Q.impute) <- (-1)*apply(Q.impute, 1, sum)
  
  beta.impute <- fitm$estimates.t[c(3,5)]
  
  # gamma <- fitm2$coefficients
  
  MAT <- lapply(SUBJ, MAT.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, maxstate=maxstate,
                Q=Q.impute, inbeta=beta.impute)
  
  efficient.info <- lapply(MAT, function(x){ (ginv(x[c(2:5), c(1:4)]) %*% x[1,])[3] })
  # efficient.info <- lapply(MAT, function(x){ (ginv(x[c(6:9), c(1:4)]) %*% x[1,])[3] })
  # efficient.info <- lapply(MAT, function(x){ x[1,3] })
  
  efficient.info <- unlist(efficient.info)
  # efficient.info <- trim(efficient.info, at = 0.995)
  # summary(efficient.info)
  # efficient.info <- Reduce('+', MAT)
  
  
  dat <- event %>% dplyr::select(id, x1, x2, x1.pred, selectR, am, Z.am, pi) %>% group_by(id) %>% slice(1)
  dat$efficient.info <- efficient.info
  
  Xs <- dat$efficient.info[dat[, selectR]==1]       # sample matrix of calibration variables
  piks <- dat$pi[dat[,selectR]==1]              # sample inclusion probabilities
  # vector of population totals of the calibration variables
  total <- sum(efficient.info)
  # summary(efficient.info)
  # computation of the g-weights
  # g2 <- calib(Xs, d=1/piks, total, method="logit", bounds=c(0.5, 1.5), max_iter=500)
  
  # g2 <- calib(Xs, d=1/piks, total, method="raking", description = TRUE)
  # g2 <- trim(g2, at = 2, lower = TRUE)
  # summary(g2)
  # summary(piks)
  calib.w <-  g2*1/piks
  # summary(calib.w)
  # dat.phase2$g2 <- g2
  # dat.phase2$efficient.info <- Xs
  
  SUBJ  <- unique(as.character(event$id[event[, selectR]==1]))
  nSUBJ <- length(SUBJ)
  
  dat.s <- data.frame(id=as.numeric(SUBJ), calib.w=calib.w)
  event <- left_join(event, dat.s, by="id")
  
  MAT.f2 <- function(whichsubj, locdata, selectR, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    calib.w <- locdata$calib.w[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta,
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
      }
    }
    
    Si <- R*calib.w*SSi1
    
    return(Si)
  }
  
  
  solve.f <- function(x) {
    Q <- matrix(0, nrow=maxstate, ncol=maxstate)
    for (k in 1:nlambda) {
      Q[states$from[k], states$to[k]] <- x[k]
    }
    diag(Q) <- (-1)*apply(Q, 1, sum)
    
    beta <- c(x[nlambda+1], x[nlambda+2])
    MAT <- lapply(SUBJ, MAT.f2, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, maxstate=maxstate,
                  Q=Q, inbeta=beta)
    
    S <- Reduce('+', MAT)
    return(S)
  }
  
  out.theta <- nleqslv(c(lam, inbeta), solve.f, control = list(allowSingular=TRUE))
  
  
  # asymptotic variance
  MAT2.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, maxstate, Q, inbeta) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x1 <- locdata$x1[1]
    x2 <- locdata$x2[1]
    R <- locdata[, selectR][1]
    calib.w <- locdata$calib.w[1]
    
    outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
                  x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
    Qx <- outQ$outQx
    dQx <- outQ$outdQx
    
    nlen <- length(locdata$w)
    
    SSi1 <- matrix(0, nrow=1, ncol=ntheta)
    MMi <- EMMi <- matrix(0, nrow=ntheta, ncol=ntheta)
    
    for (i in 1:nlen) {
      Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
      for (u in 1:ntheta) {
        dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
        SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
        
        for (j in 1:ntheta) {
          dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
          MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
          
          for (j2 in 1:maxstate) {
            if ( Pxmat[locdata$from[i], j2] != 0 ) {
              EMMi[u,j] <- EMMi[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
            }
          }
        }
        
      }
    }
    
    Si <- R*calib.w*SSi1
    Mi <- R*calib.w*MMi
    EMi <- R*calib.w*EMMi
    outMAT <- rbind(Si, Mi, EMi)
    return(outMAT)
  }
  
  Q <- matrix(0, nrow=maxstate, ncol=maxstate)
  for (k in 1:nlambda) {
    Q[states$from[k], states$to[k]] <- out.theta$x[k]
  }
  diag(Q) <- (-1)*apply(Q, 1, sum)
  beta <- c(out.theta$x[nlambda+1], out.theta$x[nlambda+2])
  
  MAT2 <- lapply(SUBJ, MAT2.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, maxstate=maxstate, 
                 Q=Q, inbeta=beta)
  nSUBJ <- length(unique(as.character(event$id)))
  
  
  SxSi <- lapply(MAT2, function(x){ x[1, c(1:(ntheta))] %*% t(x[1, c(1:(ntheta))]) })
  B <- Reduce("+", SxSi)/nSUBJ
  
  Mi <- lapply(MAT2, function(x){ x[c(2:5), c(1:(ntheta))] })
  EMi <- lapply(MAT2, function(x){ x[c(6:9), c(1:(ntheta))] })
  
  A.M <- Reduce("+", Mi)/nSUBJ
  A.EM <- Reduce("+", EMi)/nSUBJ
  
  SE.M <- sqrt( (ginv(A.M) %*% B %*% t(ginv(A.M))/nSUBJ)[3,3] )
  SE.EM <- sqrt( (ginv(A.EM) %*% B %*% t(ginv(A.EM))/nSUBJ)[3,3] )
  
  out.lam <- out.theta$x[1: nlambda]
  out.beta <- out.theta$x[(nlambda+1): (nlambda+2)]
  
  scoring <- NULL
  scoring$conv      <- out.theta$termcd
  scoring$fvec      <- out.theta$fvec
  scoring$lam      <- as.vector(out.lam)
  scoring$beta     <- as.vector(out.beta)
  scoring$se1      <- SE.EM
  scoring$se2      <- SE.M
  return(scoring)
  
}


sim.calib.f <- function(selectR, sim.num, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
  file.out <- paste("results/CALIB/", "CALIB_", selectR, "_NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                    "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "conv", "fvec.lam1", "fvec.lam2", "fvec.beta1", "fvec.beta2", "lam1"," lam2", "beta1", "beta2", "se1", "se2",
      sep=" ", "\n", append=TRUE, file=file.out)
  
  parallel.sim.f <- function(selectR, nsim, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
    
    stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
    stats$nsample2 <- nsample2
    
    if (selectR == "SRS"){
      sim.data <- read.table(paste("results/LIK/", "LIK_DIR_SRS_", "NN", stats$nsample, "_n",  stats$nsample2, "_pz", stats$p.z, "_r", stats$r , "_beta1_", round(stats$beta1, 2), "_beta2_", round(stats$beta2, 2),
                                   "_gamma1_", round(stats$gamma1, 2), "_PX11_", stats$p.x1.1, "_PX21_", stats$p.x2.1,  ".dat", sep=""), header=TRUE)
    } else if (selectR == "BAL"){
      sim.data <- read.table(paste("results/LIK/", "LIK_DIR_BAL_", "NN", stats$nsample, "_n",  stats$nsample2, "_pz", stats$p.z, "_r", stats$r , "_beta1_", round(stats$beta1, 2), "_beta2_", round(stats$beta2, 2),
                                   "_gamma1_", round(stats$gamma1, 2), "_PX11_", stats$p.x1.1, "_PX21_", stats$p.x2.1,  ".dat", sep=""), header=TRUE)
    }
    
    i=0
    check <- FALSE
    while (check==FALSE) {
      
      print(i)
      set.seed(nsim+i*1000)
      
      simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                                    gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
      
      if (selectR == "SRS"){
        simdata <- p2data.SRS.f(indata=simdata_raw, nsample2=stats$nsample2)
      } else if(selectR == "OPT") {
        simdata <- p2data.OPT.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, gamma0=stats$gamma0, gamma1=stats$gamma1)
      } else if(selectR == "PROP"){
        simdata <- p2data.PROP.f(indata=simdata_raw, nsample2=stats$nsample2)
      } else if(selectR == "BAL") {
        simdata <- p2data.BAL.f(indata=simdata_raw, nsample2=stats$nsample2)
      }
      
      out <- check.pred.f(indata=simdata, selectR=selectR, inQ0=stats$Q, nsim=nsim, inbeta=c(stats$beta1, stats$beta2), sim.data=sim.data)
      check <- out$check
      i <- i+1
    }
    
    
    fit_S <- calib.f(indata=simdata, selectR=selectR, inQ0=stats$Q, nsim=nsim, inbeta=c(stats$beta1, stats$beta2), sim.data=sim.data, g2=out$g2)
    
    warnings()
    cat(nsim, fit_S$conv, fit_S$fvec, fit_S$lam, fit_S$beta, fit_S$se1, fit_S$se2,
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
  }
  
  mclapply(1:sim.num, parallel.sim.f, selectR=selectR, nsample=nsample, nsample2=nsample2, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, 
           p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina, mc.cores=cores.num)
  
  return()
}

# check strata dist -------------------------------------------------------

sim.strata.f <- function(sim.num, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
  file.out <- paste("results/strata/", "STRATA", "_NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                    "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
  
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }
  
  cat("nsim", "all1"," all2", "all3", "all4", "all5", "all6", 
      "SRS1"," SRS2", "SRS3", "SRS4", "SRS5", "SRS6", 
      "OPT1"," OPT2", "OPT3", "OPT4", "OPT5", "OPT6", 
      "BAL1"," BAL2", "BAL3", "BAL4", "BAL5", "BAL6", 
      sep=" ", "\n", append=TRUE, file=file.out)
  
  
  parallel.sim.f <- function(nsim, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
    set.seed(nsim)
    
    stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
    stats$nsample2 <- nsample2
    
    simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                                  gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
    simdata_all <- simdata_raw$rawdata %>% dplyr::select(id, strata) %>% group_by(id) %>% slice(1)
    
    simdata_SRS <- p2data.SRS.f(indata=simdata_raw, nsample2=stats$nsample2)
    simdata_SRS <- simdata_SRS$rawdata  %>% filter(SRS==1) %>% dplyr::select(id, strata)%>% group_by(id) %>% slice(1)
    
    simdata_OPT <- p2data.OPT.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, gamma0=stats$gamma0, gamma1=stats$gamma1)
    simdata_OPT <- simdata_OPT$rawdata  %>% filter(OPT==1) %>% dplyr::select(id, strata)%>% group_by(id) %>% slice(1)
    
    simdata_BAL <- p2data.BAL.f(indata=simdata_raw, nsample2=stats$nsample2)
    simdata_BAL <- simdata_BAL$rawdata  %>% filter(BAL==1) %>% dplyr::select(id, strata)%>% group_by(id) %>% slice(1)
    
    cat(nsim, summary(factor(simdata_all$strata, levels = c("1", "2", "3", "4", "5", "6"))),  summary(factor(simdata_SRS$strata, levels = c("1", "2", "3", "4", "5", "6"))), 
        summary(factor(simdata_OPT$strata, levels = c("1", "2", "3", "4", "5", "6"))), summary(factor(simdata_BAL$strata, levels = c("1", "2", "3", "4", "5", "6"))),
        sep=" ", "\n", append=TRUE, file=file.out)
    
    return()
  }
  
  mclapply(1:sim.num, parallel.sim.f, nsample=nsample, nsample2=nsample2, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, 
           p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina, mc.cores=cores.num)
  
  return()
}



# analytical se -----------------------------------------------------------


# summary -----------------------------------------------------------------
summary.ipw.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  print(paste("results/IPW/", "IPW_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.SRS <- read.table(paste("results/IPW/", "IPW_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  print(paste("results/IPW/", "IPW_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.BAL <- read.table(paste("results/IPW/", "IPW_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  
  # stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(sd(sim.SRS$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL$beta1), 3), nsmall=3),   
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}

summary.lik.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  print(paste("results/LIK/", "LIK_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.SRS <- read.table(paste("results/LIK/", "LIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  print(paste("results/LIK/", "LIK_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.BAL <- read.table(paste("results/LIK/", "LIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  
  print(paste("results/LIK/", "LIK_OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.OPT <- read.table(paste("results/LIK/", "LIK_DIR_OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(sd(sim.SRS$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL$beta1), 3), nsmall=3),  " & ", format(round(sd(sim.OPT$beta1), 3), nsmall=3),  
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}


summary.conlik.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  print(paste("results/CONLIK/", "CONLIK_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.SRS <- read.table(paste("results/CONLIK/", "CONLIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  print(paste("results/CONLIK/", "LIK_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.BAL <- read.table(paste("results/CONLIK/", "CONLIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  print(paste("results/CONLIK/", "LIK_OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""))
  sim.OPT <- read.table(paste("results/CONLIK/", "CONLIK_DIR_OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  
  # stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(sd(sim.SRS$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL$beta1), 3), nsmall=3),
      # " & ", format(round(sd(sim.OPT$beta1), 3), nsmall=3), 
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}

summary.calib.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  sim.SRS <- read.table(paste("results/CALIB/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  sim.BAL <- read.table(paste("results/CALIB/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                              "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(sd(sim.SRS$beta1), 3), nsmall=3), " & ", format(round(mean(sim.SRS$se1), 3), nsmall=3), " & ", 
      format(round(sd(sim.BAL$beta1), 3), nsmall=3), " & ", format(round(mean(sim.BAL$se1), 3), nsmall=3), " & ",
      #" & ", format(round(mean(sim.SRS$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL$beta1)-beta1, 3), nsmall=3),
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}

summary.calib.all.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  # true
  sim.SRS.true <- read.table(paste("results/CALIB/true-para/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                   "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL.true <- read.table(paste("results/CALIB/true-para/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                   "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # impute
  sim.SRS.impute <- read.table(paste("results/CALIB/impute/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                     "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL.impute <- read.table(paste("results/CALIB/impute/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                     "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # # impute Zam/am
  # sim.SRS.impute.Zam <- read.table(paste("results/CALIB/impute-Z-am/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
  #                                          "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  # 
  # sim.BAL.impute.Zam <- read.table(paste("results/CALIB/impute-Z-am/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
  #                                          "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # EM-EM-est
  sim.SRS.EM.est <- read.table(paste("results/CALIB/EM-EM-est/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                     "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL.EM.est <- read.table(paste("results/CALIB/EM-EM-est/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                     "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # EM-re-est
  sim.SRS.EM.reest <- read.table(paste("results/CALIB/EM-re-est/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                       "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL.EM.reest <- read.table(paste("results/CALIB/EM-re-est/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                       "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(sd(sim.SRS.true$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL.true$beta1), 3), nsmall=3), " && ",
      format(round(sd(sim.SRS.impute$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL.impute$beta1), 3), nsmall=3), " && ",
      format(round(sd(sim.SRS.EM.est$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL.EM.est$beta1), 3), nsmall=3), " && ",
      format(round(sd(sim.SRS.EM.reest$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL.EM.reest$beta1), 3), nsmall=3), 
      #" & ", format(round(mean(sim.SRS$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL$beta1)-beta1, 3), nsmall=3),
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}

summary.all.ESE.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  # ML
  sim.SRS_ML <- read.table(paste("results/LIK/", "LIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_ML <- read.table(paste("results/LIK/", "LIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.OPT_ML <- read.table(paste("results/LIK/", "LIK_DIR_OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  # CML
  sim.SRS_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.OPT_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_IPW.OPT.TRUE6_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # IPW
  sim.SRS_IPW <- read.table(paste("results/IPW/", "IPW_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_IPW <- read.table(paste("results/IPW/", "IPW_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.OPT_IPW <- read.table(paste("results/IPW/", "IPW_IPW.OPT.TRUE6_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(sd(sim.SRS_ML$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL_ML$beta1), 3), nsmall=3), " & ", format(round(sd(sim.OPT_ML$beta1), 3), nsmall=3), " && ",
      format(round(sd(sim.SRS_CML$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL_CML$beta1), 3), nsmall=3), " & ", format(round(sd(sim.OPT_CML$beta1), 3), nsmall=3), " && ",
      format(round(sd(sim.SRS_IPW$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL_IPW$beta1), 3), nsmall=3), " & ", format(round(sd(sim.OPT_IPW$beta1), 3), nsmall=3), 
      # format(round(sd(sim.SRS_GR$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL_GR$beta1), 3), nsmall=3), 
      #" & ", format(round(mean(sim.SRS$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL$beta1)-beta1, 3), nsmall=3),
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}

summary.all.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  # ML
  sim.SRS_ML <- read.table(paste("results/LIK/", "LIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_ML <- read.table(paste("results/LIK/", "LIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.OPT_ML <- read.table(paste("results/LIK/", "LIK_DIR_OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  # CML
  sim.SRS_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.OPT_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_IPW.OPT.TRUE6_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # IPW
  sim.SRS_IPW <- read.table(paste("results/IPW/", "IPW_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_IPW <- read.table(paste("results/IPW/", "IPW_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.OPT_IPW <- read.table(paste("results/IPW/", "IPW_IPW.OPT.TRUE6_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  # CALIB
  # sim.SRS_GR <- read.table(paste("results/CALIB/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
  #                                "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  # sim.BAL_GR <- read.table(paste("results/CALIB/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
  #                                "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(sd(sim.SRS_ML$beta1), 3), nsmall=3), " (", format(round(mean(sim.SRS_ML$asymp.se1),3), nsmall=3), ")", " & ", format(round(sd(sim.BAL_ML$beta1), 3), nsmall=3), " (", format(round(mean(sim.BAL_ML$asymp.se1),3), nsmall=3), ")", " & " ,format(round(sd(sim.OPT_ML$beta1), 3), nsmall=3), " (", format(round(mean(sim.OPT_ML$asymp.se1),3), nsmall=3),")", " && ",
      format(round(sd(sim.SRS_CML$beta1), 3), nsmall=3), " (", format(round(mean(sim.SRS_CML$se), 3), nsmall=3), ")", " & ",  format(round(sd(sim.BAL_CML$beta1), 3), nsmall=3), " (", format(round(mean(sim.BAL_CML$se), 3), nsmall=3), ")", " & ",  format(round(sd(sim.OPT_CML$beta1), 3), nsmall=3), " (", format(round(mean(sim.OPT_CML$se), 3), nsmall=3), ")", " && ",
      format(round(sd(sim.SRS_IPW$beta1), 3), nsmall=3), " (", format(round(mean(sim.SRS_IPW$se1), 3), nsmall=3), ")", " & ", format(round(sd(sim.BAL_IPW$beta1), 3), nsmall=3), " (", format(round(mean(sim.BAL_IPW$se1), 3), nsmall=3), ")", " & ", format(round(sd(sim.OPT_IPW$beta1), 3), nsmall=3), " (", format(round(mean(sim.OPT_IPW$se1), 3), nsmall=3), ")",
      # format(round(sd(sim.SRS_GR$beta1), 3), nsmall=3), " (", format(round(mean(sim.SRS_GR$se1), 3), nsmall=3), ")", " & ", format(round(sd(sim.BAL_GR$beta1), 3), nsmall=3), " (", format(round(mean(sim.BAL_GR$se1), 3), nsmall=3), ")", 
      # format(round(sd(sim.SRS_GR$beta1), 3), nsmall=3), " & ", format(round(sd(sim.BAL_GR$beta1), 3), nsmall=3), 
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}





summary.allmean.f <- function(nsample, nsample2,  r=1.1, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4){
  # ML
  sim.SRS_ML <- read.table(paste("results/LIK/", "LIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_ML <- read.table(paste("results/LIK/", "LIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.OPT_ML <- read.table(paste("results/LIK/", "LIK_DIR_OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  # CML
  sim.SRS_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_CML <- read.table(paste("results/CONLIK/", "CONLIK_DIR_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  # IPW
  sim.SRS_IPW <- read.table(paste("results/IPW/", "IPW_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_IPW <- read.table(paste("results/IPW/", "IPW_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                  "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  # CALIB
  sim.SRS_GR <- read.table(paste("results/CALIB/", "CALIB_SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  sim.BAL_GR <- read.table(paste("results/CALIB/", "CALIB_BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z, "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
                                 "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep=""), header=TRUE)
  
  cat(nsample2," & ", p.z, " & ", p.x1.1," & ", round(beta1,2), " & ", 
      format(round(mean(sim.SRS_ML$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL_ML$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.OPT_ML$beta1)-beta1, 3), nsmall=3), " && ",
      format(round(mean(sim.SRS_CML$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL_CML$beta1)-beta1, 3), nsmall=3), " && ", 
      format(round(mean(sim.SRS_IPW$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL_IPW$beta1)-beta1, 3), nsmall=3), " && ", 
      format(round(mean(sim.SRS_GR$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL_GR$beta1)-beta1, 3), nsmall=3), 
      #" & ", format(round(mean(sim.SRS$beta1)-beta1, 3), nsmall=3), " & ", format(round(mean(sim.BAL$beta1)-beta1, 3), nsmall=3),
      "\\\\ ", "\n", sep="", append=TRUE, file=file.out)
  
  return()
}

# illustration plot -------------------------------------------------------

plot.f <- function(nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
  set.seed(nsim)
  
  stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
  stats$nsample2 <- nsample2
  
  simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                                gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
  # simdata <- phase2data.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, gamma0=stats$gamma0, gamma1=stats$gamma1)
  event <- simdata_raw$event
  rawdata <- simdata_raw$rawdata
  maxstate <- ncol(inQ)
  nsample <- length(unique(event$id))
  
  # Simple random sampling
  id.random <- sample(1:nsample, nsample2, replace = FALSE)
  event$R.random <- ifelse(event$id %in% id.random, 1, 0)
  rawdata$R.random <- ifelse(rawdata$id %in% id.random, 1, 0)
  
  # optimal sampling
  fitm <- msm(states ~ times, subject=id, data=rawdata, 
              covariates = ~ x2,
              constraint = list(x2 = c(1,1)),
              qmatrix=inQ, obstype=1, opt.method="fisher",
              gen.inits=TRUE, center = FALSE)
  
  SUBJ  <- unique(as.character(event$id))
  
  Mmu.f <- function(whichsubj, locdata, maxstate) {
    locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
    locdata$w <- locdata$estop - locdata$estart
    
    x2 <- locdata$x2[1]
    nlen <- length(locdata$w)
    
    Qxmat <- qmatrix.msm(fitm, covariates=list (x2 = x2), ci="none" )
    
    lamx <- NULL
    for (j in 1:maxstate) {
      for (k in 1:maxstate) {
        if ( (j != k) && (inQ[j,k] > 0) ) {
          lamx <- c(lamx, Qxmat[j,k])
        }
      }
    } 
    
    lamx <- c(lamx, 0)
    Mmu <- 0
    for (i in 1:nlen) {
      Pxmat <- pmatrix.msm(fitm, t=locdata$w[i], covariates=list (x2 = x2))
      
      # dP_klr(X)/dmu at j
      dPxmuj.f <- function(k,l,j){
        Ckjl <- prod( lamx[k:(l-1)] )/prod( fun.zero.omit(lamx[k:l]-lamx[j]) ) 
        Ckjl <- ifelse(k==l & k==j, 1, Ckjl)
        
        dPxmuj <- Ckjl*exp(-lamx[j]*locdata$w[i])*(-lamx[j]*locdata$w[i])
        return(dPxmuj)
      }
      
      
      dPxmu <- sum(apply(X=as.matrix(c(locdata$from[i]:locdata$to[i])), 1, FUN=dPxmuj.f, k=locdata$from[i], l=locdata$to[i]))
      Mmu <- Mmu + 1/Pxmat[locdata$from[i],locdata$to[i]]*dPxmu
    }
    return(Mmu)
  }
  
  Mmu <- lapply(SUBJ, Mmu.f, locdata=event, maxstate=maxstate)
  Mmu <- as.numeric(Mmu)
  
  x2_dat <- rawdata %>% group_by(id) %>% slice(1) 
  x2 <- x2_dat$x2
  
  Mmu[which(x2==0)] <- Mmu[which(x2==0)]*sqrt(expit(gamma0+gamma1*0)*(1-expit(gamma0+gamma1*0)))
  Mmu[which(x2==1)] <- Mmu[which(x2==1)]*sqrt(expit(gamma0+gamma1*1)*(1-expit(gamma0+gamma1*1)))
  
  dat <- data.frame(id=1:nsample, Mmu=Mmu, x2=x2) %>% arrange(Mmu)
  dat.x2.1 <- dat %>% filter(x2==1) %>% arrange(Mmu)
  dat.x2.0 <- dat %>% filter(x2==0) %>% arrange(Mmu)
  
  id.optimal <- c(dat.x2.1$id[1:(nsample2/4)], dat.x2.1$id[(nrow(dat.x2.1)-(nsample2/4)+1):nrow(dat.x2.1)],
                  dat.x2.0$id[1:(nsample2/4)], dat.x2.0$id[(nrow(dat.x2.0)-(nsample2/4)+1):nrow(dat.x2.0)])
  
  
  event$R.optimal <- ifelse(event$id %in% id.optimal, 1, 0)
  rawdata$R.optimal <- ifelse(rawdata$id %in% id.optimal, 1, 0)
  
  #test <- rawdata %>% filter(R.optimal == 1)
  
  outdata <- NULL
  outdata$rawdata <- rawdata
  outdata$event   <- event
}


# sim.likeli.f <- function(sim.num, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
#   # file.out0 <- paste("results/SRS_NAIVE_", "NN", nsample, "_n",  nsample2, "_pz", p.z ,"_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
#   #                    "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
#   file.out1 <- paste("results/SRS_", "NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
#                      "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
#   file.out2 <- paste("results/OPT_", "NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
#                      "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
#   file.out3 <- paste("results/PROP_", "NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
#                      "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
#   file.out4 <- paste("results/BAL_", "NN", nsample, "_n",  nsample2, "_pz", p.z , "_r", r , "_beta1_", round(beta1, 2), "_beta2_", round(beta2, 2),
#                      "_gamma1_", round(gamma1, 2), "_PX11_", p.x1.1, "_PX21_", p.x2.1,  ".dat", sep="")
#   
#   
#   # if ( file.exists(as.character(file.out0)) ) { unlink(as.character(file.out0)) }
#   if ( file.exists(as.character(file.out1)) ) { unlink(as.character(file.out1)) }
#   if ( file.exists(as.character(file.out2)) ) { unlink(as.character(file.out2)) }
#   if ( file.exists(as.character(file.out3)) ) { unlink(as.character(file.out3)) }
#   if ( file.exists(as.character(file.out4)) ) { unlink(as.character(file.out4)) }
#   
#   # cat("lam1", "lam2", "lam1_beta1", "lam2_beta1", "lam1_beta2", "lam2_beta2",
#   #     sep=" ", "\n", append=TRUE, file=file.out0)
#   cat("nsim", "iter", "tol", "lam1"," lam2", "beta1", "beta2","gamma1", "gamma2",
#       sep=" ", "\n", append=TRUE, file=file.out1)
#   cat("nsim", "iter", "tol", "lam1"," lam2", "beta1", "beta2","gamma1", "gamma2",
#       sep=" ", "\n", append=TRUE, file=file.out2)
#   cat("nsim", "iter", "tol", "lam1"," lam2", "beta1", "beta2","gamma1", "gamma2",
#       sep=" ", "\n", append=TRUE, file=file.out3)
#   cat("nsim", "iter", "tol", "lam1"," lam2", "beta1", "beta2","gamma1", "gamma2",
#       sep=" ", "\n", append=TRUE, file=file.out4)
#   
#   parallel.sim.f <- function(nsim, nsample, nsample2, r, p.z, beta1, beta2, gamma1, p.x1.1, p.x2.1, inA, lam.ina){
#     set.seed(nsim)
#     
#     stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
#     #stats$a <- seq(0,stats$A, by=(1/5))
#     stats$nsample2 <- nsample2
#     
#     simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
#                                   gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)
#     simdata <- phase2data.f(indata=simdata_raw, nsample2=stats$nsample2, inQ=stats$Q, gamma0=stats$gamma0, gamma1=stats$gamma1)
#     
#     
#     # fit_SRS <- msm(states ~ times, subject=id, data=simdata$rawdata[simdata$rawdata$R.random==1, ],
#     #                qmatrix=stats$Q, covariates = ~ x1 + x2, constraint = list(x1 = c(1,1), x2 = c(1,1)), 
#     #                opt.method="fisher", gen.inits=TRUE, center=FALSE)
#     # cat(fit_SRS$estimates.t, sep=" ", "\n", append=TRUE, file=file.out0)
#     
#     # SRS
#     fit_SRS_S <- scoring.f(indata=simdata$event, selectR="R.random", nlambda=length(stats$lam), nbeta=2, ngamma=2,
#                            inQ0=stats$Q, inbeta=c(stats$beta1, stats$beta2), 
#                            ingamma= c(stats$gamma0, stats$gamma1) )
#     cat(nsim, fit_SRS_S$iter, fit_SRS_S$tol, fit_SRS_S$lam, fit_SRS_S$beta, fit_SRS_S$gamma,
#         sep=" ", "\n", append=TRUE, file=file.out1)
#     
#     # OPT
#     fit_OPT_S <- scoring.f(indata=simdata$event, selectR="R.optimal", nlambda=length(stats$lam), nbeta=2, ngamma=2,
#                            inQ0=stats$Q, inbeta=c(stats$beta1, stats$beta2), 
#                            ingamma= c(stats$gamma0, stats$gamma1) )
#     cat(nsim, fit_OPT_S$iter, fit_OPT_S$tol, fit_OPT_S$lam, fit_OPT_S$beta, fit_OPT_S$gamma,
#         sep=" ", "\n", append=TRUE, file=file.out2)
#     
#     # PROP
#     fit_PROP_S <- scoring.f(indata=simdata$event, selectR="R.prop", nlambda=length(stats$lam), nbeta=2, ngamma=2,
#                             inQ0=stats$Q, inbeta=c(stats$beta1, stats$beta2), 
#                             ingamma= c(stats$gamma0, stats$gamma1) )
#     cat(nsim, fit_PROP_S$iter, fit_PROP_S$tol, fit_PROP_S$lam, fit_PROP_S$beta, fit_PROP_S$gamma,
#         sep=" ", "\n", append=TRUE, file=file.out3)
#     
#     # BAL
#     fit_BAL_S <- scoring.f(indata=simdata$event, selectR="R.bal", nlambda=length(stats$lam), nbeta=2, ngamma=2,
#                            inQ0=stats$Q, inbeta=c(stats$beta1, stats$beta2), 
#                            ingamma= c(stats$gamma0, stats$gamma1) )
#     cat(nsim, fit_BAL_S$iter, fit_BAL_S$tol, fit_BAL_S$lam, fit_BAL_S$beta, fit_BAL_S$gamma,
#         sep=" ", "\n", append=TRUE, file=file.out4)
#     
#     return()
#   }
#   
#   mclapply(1:sim.num, parallel.sim.f, nsample=nsample, nsample2=nsample2, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, 
#            p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina, mc.cores=cores.num)
#   
#   return()
# }


# auxiliary.f <- function(indata, selectR, inQ0, inbeta, ingamma) {
#   maxstate <- nrow(inQ0)
#   
#   event <- indata$event
#   rawdata <- indata$rawdata
#   
#   lam <- NULL
#   states <- NULL
#   for (j in 1:maxstate) {
#     for (k in 1:maxstate) {
#       if ( (j != k) && (inQ0[j,k] > 0) ) {
#         lam <- c(lam, inQ0[j,k])
#         states <- rbind(states, data.frame(from=j, to=k))
#       }
#     }
#   }  
#   
#   beta <- inbeta
#   gamma <- ingamma
#   
#   nlambda <- length(lam)
#   nbeta <- length(beta)
#   ngamma <- length(gamma)
#   ntheta <- nlambda + nbeta
#   
#   SUBJ  <- unique(as.character(event$id))
#   nSUBJ <- length(SUBJ)
#   
#   MAT.f <- function(whichsubj, locdata, selectR, nlambda, nbeta, ngamma, maxstate, inQ0, beta, gamma) {
#     locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
#     locdata$w <- locdata$estop - locdata$estart
#     
#     x1 <- locdata$x1[1]
#     x2 <- locdata$x2[1]
#     R <- locdata[, selectR][1]
#     
#     outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
#                   x1=x1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
#     Qx <- outQ$outQx
#     dQx <- outQ$outdQx
#     
#     # pseudo
#     outQ.x1.1 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
#                        x1=1, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
#     Qx.x1.1 <- outQ.x1.1$outQx
#     dQx.x1.1 <- outQ.x1.1$outdQx
#     
#     outQ.x1.0 <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
#                        x1=0, x2=x2, inQ=inQ0, inbeta1=beta[1], inbeta2=beta[2])
#     Qx.x1.0 <- outQ.x1.0$outQx
#     dQx.x1.0 <- outQ.x1.0$outdQx
#     
#     nlen <- length(locdata$w)
#     PD.x1.1 <- 1 # P(D | X1=1, X2)
#     PD.x1.0 <- 1 # P(D | X1=0, X2)
#     
#     SSi1 <- SSi1.x1.1 <- SSi1.x1.0 <- matrix(0, nrow=1, ncol=ntheta)
#     SSi2 <- SSi2.x1.1 <- SSi2.x1.0 <- matrix(0, nrow=1, ncol=ngamma)
#     
#     IIi1 <- IIi1.x1.1 <- IIi1.x1.0 <- matrix(0, nrow=ntheta, ncol=ntheta)
#     IIi2 <- matrix(0, nrow=ngamma, ncol=ngamma)
#     
#     for (i in 1:nlen) {
#       Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
#       Pxmat.x1.1 <- round(MatrixExp(Qx.x1.1, t=locdata$w[i], method="series"), 16)
#       Pxmat.x1.0 <- round(MatrixExp(Qx.x1.0, t=locdata$w[i], method="series"), 16)
#       
#       PD.x1.1 <- PD.x1.1*Pxmat.x1.1[locdata$from[i], locdata$to[i]]
#       PD.x1.0 <- PD.x1.0*Pxmat.x1.0[locdata$from[i], locdata$to[i]]
#       
#       for (u in 1:ntheta) {
#         dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
#         dPxmatu.x1.1 <- dPxt.f(u, locdata$w[i], Qx.x1.1, dQx.x1.1, maxstate)
#         dPxmatu.x1.0 <- dPxt.f(u, locdata$w[i], Qx.x1.0, dQx.x1.0, maxstate)
#         
#         SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
#         SSi1.x1.1[1,u] <- SSi1.x1.1[1,u] + dPxmatu.x1.1[locdata$from[i],locdata$to[i]] / Pxmat.x1.1[locdata$from[i],locdata$to[i]]
#         SSi1.x1.0[1,u] <- SSi1.x1.0[1,u] + dPxmatu.x1.0[locdata$from[i],locdata$to[i]] / Pxmat.x1.0[locdata$from[i],locdata$to[i]]
#         
#         for (j in 1:ntheta) {
#           dPxmatj <- dPxt.f(j, locdata$w[i], Qx, dQx, maxstate)
#           dPxmatj.x1.1 <- dPxt.f(j, locdata$w[i], Qx.x1.1, dQx.x1.1, maxstate)
#           dPxmatj.x1.0 <- dPxt.f(j, locdata$w[i], Qx.x1.0, dQx.x1.0, maxstate)
#           # MMi[u,j] <- MMi[u,j] + ( (1/(Pxmat[locdata$from[i],locdata$to[i]]^2))*dPxmatu[locdata$from[i],locdata$to[i]]*dPxmatj[locdata$from[i],locdata$to[i]] )
#           
#           for (j2 in 1:maxstate) {
#             if ( Pxmat[locdata$from[i], j2] != 0 ) {
#               IIi1[u,j] <- IIi1[u,j] + ( (1/Pxmat[locdata$from[i],j2])*dPxmatu[locdata$from[i],j2]*dPxmatj[locdata$from[i],j2] )
#               IIi1.x1.1[u,j] <- IIi1.x1.1[u,j] + ( (1/Pxmat.x1.1[locdata$from[i],j2])*dPxmatu.x1.1[locdata$from[i],j2]*dPxmatj.x1.1[locdata$from[i],j2] )
#               IIi1.x1.0[u,j] <- IIi1.x1.0[u,j] + ( (1/Pxmat.x1.0[locdata$from[i],j2])*dPxmatu.x1.0[locdata$from[i],j2]*dPxmatj.x1.0[locdata$from[i],j2] )
#               
#             }
#           }
#           
#         }
#         
#       }
#     }
#     
#     
#     SSi2[1,1] <- x1-expit(gamma[1]+gamma[2]*x2)
#     SSi2[1,2] <- (x1-expit(gamma[1]+gamma[2]*x2))*x2
#     SSi2.x1.1[1,1] <- 1-expit(gamma[1]+gamma[2]*x2)
#     SSi2.x1.1[1,2] <- (1-expit(gamma[1]+gamma[2]*x2))*x2
#     SSi2.x1.0[1,1] <- 0-expit(gamma[1]+gamma[2]*x2)
#     SSi2.x1.0[1,2] <- (0-expit(gamma[1]+gamma[2]*x2))*x2
#     
#     IIi2[1,1] <- expit(gamma[1]+gamma[2]*x2)*( 1/(1+exp(gamma[1]+gamma[2]*x2)) )
#     IIi2[1,2] <- IIi2[2,1] <- expit(gamma[1]+gamma[2]*x2)*( 1/(1+exp(gamma[1]+gamma[2]*x2)) )*x2
#     IIi2[2,2] <- expit(gamma[1]+gamma[2]*x2)*( 1/(1+exp(gamma[1]+gamma[2]*x2)) )*x2^2
#     
#     eta <- ( PD.x1.1*expit(gamma[1]+gamma[2]*x2) )/( PD.x1.1*expit(gamma[1]+gamma[2]*x2) + PD.x1.0*(1-expit(gamma[1]+gamma[2]*x2)) )
#     
#     Si1 <- R*SSi1 + (1-R)*( SSi1.x1.1*eta + SSi1.x1.0*(1-eta) )
#     Si2 <- R*SSi2 + (1-R)*( SSi2.x1.1*eta + SSi2.x1.0*(1-eta) )
#     Si <- cbind(Si1, Si2)
#     
#     Ii1 <- R*IIi1 + (1-R)*( IIi1.x1.1*eta +  IIi1.x1.0*(1-eta) - t(SSi1.x1.1)%*%SSi1.x1.1*eta - t(SSi1.x1.0)%*%SSi1.x1.0*(1-eta) + 
#                               t((SSi1.x1.1*eta + SSi1.x1.0*(1-eta)))%*%(SSi1.x1.1*eta + SSi1.x1.0*(1-eta)) )
#     Ii2 <- R*IIi2 + (1-R)*( IIi2 -  t(SSi2.x1.1)%*%SSi2.x1.1*eta - t(SSi2.x1.0)%*%SSi2.x1.0*(1-eta) + 
#                               t(( SSi2.x1.1*eta + SSi2.x1.0*(1-eta) ))%*%( SSi2.x1.1*eta + SSi2.x1.0*(1-eta) ) )
#     
#     outI <- matrix(0, nrow = ntheta+ngamma, ncol = ntheta+ngamma)
#     outI[1:ntheta, 1:ntheta] <- Ii1
#     outI[(ntheta+1):(ntheta+ngamma), (ntheta+1):(ntheta+ngamma)] <- Ii2
#     
#     # solve(outI+diag(rep(1e-6,6)))%*% t(Si)
#     auxi <-  ginv(outI) %*% t(Si)
#     
#     return(auxi[3])
#   }
#   
#   auxi <- lapply(SUBJ, MAT.f, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, ngamma=ngamma, maxstate=maxstate, 
#                  inQ0=inQ0, beta=beta, gamma=gamma)
#   
#   auxi <- unlist(auxi)
#   
#   # summary(auxi)
#   return(auxi)
# }  
# 
# 
# calib.f <- function(indata, selectR, inQ0,  inbeta, ingamma) { # nlambda, nbeta,
#   
#   out <- pi.f(indata=indata, selectR=selectR, inQ0=inQ0)
#   event <- out$event
#   fit.pi <- out$fit.pi
#   # summary(fit.pi)
#   rawdata <- indata$rawdata
#   rawdata.phase2 <- rawdata %>% filter(!!sym(selectR) == 1)
#   
#   maxstate <- nrow(inQ0)
#   nsample <- length(unique(event$id))
#   
#   lam <- NULL
#   states <- NULL
#   for (j in 1:maxstate) {
#     for (k in 1:maxstate) {
#       if ( (j != k) && (inQ0[j,k] > 0) ) {
#         lam <- c(lam, inQ0[j,k])
#         states <- rbind(states, data.frame(from=j, to=k))
#       }
#     }
#   }  
#   
#   # fitm <- msm(states ~ times, subject=id, data=rawdata.phase2,
#   #             covariates = ~ x1 + x2,
#   #             constraint = list(x1 = c(1,1), x2 = c(1,1)),
#   #             qmatrix=inQ0, obstype=1, opt.method="fisher",
#   #             gen.inits=TRUE, center = FALSE)
#   # 
#   # beta <- fitm$estimates.t[c(3,5)]
#   # lam <- fitm$estimates.t[c(1,2)]
#   # 
#   # fitx1 <- glm(x1 ~ x2, data = rawdata %>% group_by(id) %>% slice(1), family = "binomial")
#   # gamma <- fitx1$coefficients
#   
#   
#   
#   # Q <- matrix(0, nrow=maxstate, ncol=maxstate)
#   # for (k in 1:length(lam)) {
#   #   Q[states$from[k], states$to[k]] <- lam[k]
#   # }
#   # diag(Q) <- (-1)*apply(Q, 1, sum)
#   
#   
#   
#   nlambda <- 2
#   nbeta <- length(inbeta)
#   ngamma <- length(ingamma)
#   ntheta <- nlambda + nbeta
#   
#   
#   auxi <- auxiliary.f(indata=indata, selectR=selectR, inQ0=inQ0, inbeta=inbeta, ingamma=ingamma)
#   auxi.trim <- trim(auxi, at = 0.995, lower = TRUE)
#   # sum(auxi) sum(auxi.trim)
#   # summary(auxi)
#   
#   dat <- event %>% dplyr::select(id, x1, x2, selectR, am, Z.am, pi) %>% group_by(id) %>% slice(1)
#   dat$auxi <- auxi
#   dat$auxi.trim <- auxi.trim
#   
#   
#   Xs <- dat$auxi.trim[dat[, selectR]==1]       # sample matrix of calibration variables
#   piks <- dat$pi[dat[, selectR]==1]              # sample inclusion probabilities
#   # vector of population totals of the calibration variables
#   total <- sum(auxi.trim)
#   # summary(auxi)
#   # computation of the g-weights
#   g2 <- calib(Xs, d=1/piks, total, method="raking")
#   # g2 <- calib(Xs, d=1/piks, total, method="raking", description=TRUE)
#   # summary(g2)
#   # summary(piks)
#   calib.w <-  g2*1/piks
#   # summary(calib.w) summary(1/piks)
#   
#   SUBJ  <- unique(as.character(event$id[event[, selectR]==1]))
#   nSUBJ <- length(SUBJ)
#   
#   dat.s <- data.frame(id=as.numeric(SUBJ), calib.w=calib.w)
#   event <- left_join(event, dat.s, by="id")
#   
#   MAT.f2 <- function(whichsubj, locdata, selectR, nlambda, nbeta, maxstate, Q, inbeta) {
#     locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
#     locdata$w <- locdata$estop - locdata$estart
#     
#     x1 <- locdata$x1[1]
#     x2 <- locdata$x2[1]
#     R <- locdata[, selectR][1]
#     calib.w <- locdata$calib.w[1]
#     
#     outQ <- dQx.f(states=states, maxstate=maxstate, nlambda=nlambda, nbeta=nbeta, 
#                   x1=x1, x2=x2, inQ=Q, inbeta1=inbeta[1], inbeta2=inbeta[2])
#     Qx <- outQ$outQx
#     dQx <- outQ$outdQx
#     
#     nlen <- length(locdata$w)
#     
#     SSi1 <- matrix(0, nrow=1, ncol=ntheta)
#     
#     for (i in 1:nlen) {
#       Pxmat <- round(MatrixExp(Qx, t=locdata$w[i], method="series"), 16)
#       
#       for (u in 1:ntheta) {
#         dPxmatu <- dPxt.f(u, locdata$w[i], Qx, dQx, maxstate)
#         SSi1[1,u] <- SSi1[1,u] + dPxmatu[locdata$from[i],locdata$to[i]] / Pxmat[locdata$from[i],locdata$to[i]]
#       }
#     }
#     
#     Si <- R*calib.w*SSi1
#     
#     return(Si)
#   }
#   
#   
#   solve.f <- function(x) {
#     Q <- matrix(0, nrow=maxstate, ncol=maxstate)
#     for (k in 1:nlambda) {
#       Q[states$from[k], states$to[k]] <- x[k]
#     }
#     diag(Q) <- (-1)*apply(Q, 1, sum)
#     
#     beta <- c(x[nlambda+1], x[nlambda+2])
#     MAT <- lapply(SUBJ, MAT.f2, locdata=event, selectR=selectR, nlambda=nlambda, nbeta=nbeta, maxstate=maxstate, 
#                   Q=Q, inbeta=beta)
#     
#     S <- Reduce('+', MAT)
#     return(S)
#   }
#   
#   out.theta <- nleqslv(c(lam, inbeta), solve.f, control = list(allowSingular=TRUE))
#   
#   out.lam <- out.theta$x[1: nlambda]
#   out.beta <- out.theta$x[(nlambda+1): (nlambda+2)]
#   
#   scoring <- NULL
#   scoring$conv      <- out.theta$termcd
#   scoring$fvec      <- out.theta$fvec
#   scoring$lam      <- as.vector(out.lam)
#   scoring$beta     <- as.vector(out.beta)
#   
#   return(scoring)
#   
# }