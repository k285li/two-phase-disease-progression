
cores.num <- 100
source("share.func.R")


# Maximum likelihood 
for (nsample2 in c(300, 600)) {
  for (p.z in c(0.4, 0.8)) {
    for (p.x1.1 in c(0.2, 0.5)) {
      for (gamma1 in c(log(2))) { 
        for (beta1 in c(0, log(1.25), log(2))) {
          sim.likeli.dir.f(selectR="SRS", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.likeli.dir.f(selectR="BAL", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.likeli.dir.f(selectR="OPT", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
        }
      }
    }
  }
}


# IPW
for (nsample2 in c(300, 600)) {
  for (p.z in c(0.4, 0.8)) {
    for (p.x1.1 in c(0.2, 0.5)) {
      for (gamma1 in c(log(2))) { 
        for (beta1 in c(0, log(1.25), log(2))) {
          sim.ipw.f(selectR="SRS", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.ipw.f(selectR="BAL", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.ipw.f(selectR="IPW.OPT.TRUE6", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          
        }
      }
    }
  }
}



# Conditional maximum likelihood 
for (nsample2 in c(300, 600)) { 
  for (p.z in c(0.4, 0.8)) {
    for (p.x1.1 in c(0.2, 0.5)) {
      for (gamma1 in c(log(2))) { 
        for (beta1 in c(0, log(1.25), log(2))) {
          sim.conlikeli.dir.f(selectR="SRS", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.conlikeli.dir.f(selectR="BAL", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
          sim.conlikeli.dir.f(selectR="IPW.OPT.TRUE6", sim.num=1000, nsample=2000, nsample2=nsample2, r=1.1, p.z=p.z, beta1=beta1, beta2=log(2), gamma1=gamma1, p.x1.1=p.x1.1, p.x2.1=0.5, inA=1, lam.ina=4)
        }
      }
    }
  }
}

