library(msm)
library(cubature)
library(rje)
library(xtable)
library(MASS)

cores.num <- 1
source("share.func.R")


for (nsample2 in c(300, 600)) {
  for (p.z in c(0.4, 0.8)) {
    for (p.x1.1 in c(0.2, 0.5)) {
      for (gamma1 in c(log(2))) { # log(2), log(4)
        for (beta1 in c(0, log(1.25), log(2))) {
          sim.likeli.dir.f(selectR="SRS", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.likeli.dir.f(selectR="BAL", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.likeli.dir.f(selectR="OPT", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
        }
      }
    }
  }
}

for (nsample2 in c(300, 600)) { 
  for (p.z in c(0.4, 0.8)) {
    for (p.x1.1 in c(0.2, 0.5)) {
      for (gamma1 in c(log(2), )) { # log(2), log(4)
        for (beta1 in c(0, log(1.25), log(2))) {
          sim.ipw.f(selectR="SRS", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.ipw.f(selectR="BAL", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
        }
      }
    }
  }
}

for (nsample2 in c(300, 600)) { 
  for (p.z in c(0.4, 0.8)) {
    for (p.x1.1 in c(0.2, 0.5)) {
      for (gamma1 in c(log(2))) { # log(2), log(4)
        for (beta1 in c(0, log(1.25), log(2))) {
          sim.conlikeli.dir.f(selectR="SRS", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.conlikeli.dir.f(selectR="BAL", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
        }
      }
    }
  }
}



for (nsample2 in c(300, 600)) { 
  for (p.z in c(0.4)) {
    for (p.x1.1 in c(0.2)) {
      for (gamma1 in c(log(2))) { # log(2), log(4)
        for (beta1 in c(0)) {
          sim.calib.f(selectR="SRS", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.calib.f(selectR="BAL", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          
        }
      }
    }
  }
}

# mycov <- fitL$cov[upper.tri(fitL$cov)]
# myvar <- diag(fitL$cov)
# 
mat<- matrix(1, 4, 4)
diag(mat) <- as.numeric(sim.SRS[1, 6:9])
mat[upper.tri(mat)] <- as.numeric(sim.SRS[1, 10:15])
mat <- t(mat)
mat[upper.tri(mat)] <- as.numeric(sim.SRS[1,10:15])
mat
diag(as.numeric(sim.SRS[1, 2:5])) %*% mat %*% t(diag(as.numeric(sim.SRS[1, 2:5])))

sim.ipw.f(selectR="SRS", sim.num=2, nsample=2000, nsample2=400, r=1.1, p.z=0.4, beta1=0, beta2=log(2), gamma1=log(2), p.x1.1=0.2, p.x2.1=0.5, inA=1, lam.ina=4)
sim.calib.f(selectR="SRS", sim.num=2, nsample=2000, nsample2=400, r=1.1, p.z=0.4, beta1=0, beta2=log(2), gamma1=log(2), p.x1.1=0.2, p.x2.1=0.5, inA=1, lam.ina=4)



source("share.func.R")
selectR="IPW.OPT.TRUE6"
nsample2=300
p.z=0.4
p.x1.1=0.2
gamma1=log(2)
beta1=0
nsample=2000
r=1.1
beta2=log(2)
p.x2.1=0.5
inA=1
lam.ina=4
nsim=79
