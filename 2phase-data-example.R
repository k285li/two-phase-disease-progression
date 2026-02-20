# Example of two-phase design with 
# pseudo-score residual-dependent subsampling using
# maximum likelihood estimation

# Transition intensities are modeled as
# lam_k(Xi; Theta) = lam_{k0}*exp(beta1*Xi1+ beta2*Xi2), k = 0, 1

# load functions
source("share.func.R")

# example of parameter setting
nsample2 <- 300
p.z <- 0.4
p.x1.1 <- 0.2
gamma1 <- log(2)
beta1 <- log(2)
nsample <- 2000
r <- 1.1
beta2 <- log(2)
p.x2.1 <- 0.5
inA <- 1
lam.ina <- 4


# simulated data generation
stats <- getstats.f(nsample=nsample, r=r, p.z=p.z, beta1=beta1, beta2=beta2, gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=p.x2.1, inA=inA, lam.ina=lam.ina)
stats$nsample2 <- nsample2
simdata_raw <- generatedata.f(nsample=stats$nsample, inQ=stats$Q, inA=stats$A, lam.ina=stats$lam.ina, 
                              gamma0=stats$gamma0, gamma1=stats$gamma1, p.x2.1=stats$p.x2.1, beta1=stats$beta1, beta2=stats$beta2)

# simulated data with selected subjects using pseudo-score residual-dependent subsampling
simdata <- p2data.SRS.f(indata=simdata_raw, nsample2=stats$nsample2)


# estimate the biomarker effect using maximum likelihood estimation
fitL <- logL.msm.f(indata=simdata$event, qmat=stats$Q, selectR="SRS", beta1=stats$beta1, beta2=stats$beta2, 
                   gamma0=stats$gamma0, gamma1=stats$gamma1)

# beta1 is the biomarker effect of interest
cat("lam1"," lam2", "beta1", "beta2", sep=" ", "\n")
cat(exp(fitL$estimate[1:2]), fitL$estimate[3:4])
