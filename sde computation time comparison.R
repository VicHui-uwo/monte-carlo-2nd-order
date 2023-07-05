library(sde)
library(ggplot2)
library(dplyr)

seed <- 20230331
seed <- 10

set.seed(seed)
### Changing mu
T = 1
N = 252 * T
delta = T/N
sigma = 0.25
r = 0.25
s=expression(0.25 * x)
d=expression(0.05 * x)

# Black Scholes True solution function
sde.bs.true <- function(X0=100, theta=c(0.05, 0.25), delta.t=0.003968254, T=1, seed=0){
  set.seed(seed)
  t <- seq(0, T, delta.t)
  dW <- c(0)
  dW <- append(dW, rnorm(T/delta.t, 0, sqrt(delta.t)))
  B <- cumsum(dW)
  return (X0*exp((theta[1]-(theta[2]**2)/2)*t + theta[2]*B))
}

# Runge Kutta function
sde.rk <- function(X0=100, theta=c(0.05, 0.25), delta.t=0.003968254, T=1, seed=0){
  set.seed(seed)
  mu <- theta[1]
  sigma <- theta[2]
  dW <- c(0)
  dW <- append(dW, rnorm(T / delta.t, 0, sqrt(delta.t)))
  t <- seq(0, T, delta.t)
  sde.runge.kutta <- numeric(length(t))
  sde.runge.kutta[1] <- X0
  for (i in 2:length(t)){
    upsilon <- sde.runge.kutta[i-1] + mu*sde.runge.kutta[i-1]*delta.t +
      sigma*sde.runge.kutta[i-1]*sqrt(delta.t)
    
    sde.runge.kutta[i] <- sde.runge.kutta[i-1] + mu*sde.runge.kutta[i-1]*delta.t +
      sigma*sde.runge.kutta[i-1]*dW[i] + 
      0.5*sigma*(upsilon-sde.runge.kutta[i-1])*(dW[i]**2-delta.t)*(delta.t**(-0.5))
  }
  return (sde.runge.kutta)
}

time.true <- system.time(replicate(5000, sde.bs.true(delta.t=delta, T=T)))
time.rk <- system.time(replicate(5000, sde.rk(delta.t=delta, T=T)))
defaultW <- getOption("warn") 
options(warn = -1) 
time.mil <- system.time(replicate(5000, sde.sim(X0=100, drift=d, sigma=s, sigma.x=sigma, method="milstein", N=N, T=T)))
time.mil2 <- system.time(replicate(5000, sde.sim(X0=100, drift=d, sigma=s, sigma.x=sigma, drift.x = r, sigma.xx=0, drift.xx=0, method="milstein2", N=N, T=T)))
time.euler <- system.time(replicate(5000, sde.sim(X0=100, drift=d, sigma=s, sigma.x=sigma, method="euler", N=N, T=T)))
options(warn = defaultW)

par.fun.mil <- function(x) sde.sim(X0=100, drift=expression(0.05*x), sigma = expression(0.25 * x), sigma.x = expression(0.25), method="euler", N=252, T=1)


cl <- makeCluster(detectCores())
clusterSetRNGStream(cl, iseed=seed)
time.par.true <- system.time(parSapply(cl,rep(100,5000),sde.bs.true))
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterSetRNGStream(cl, iseed=seed)
time.par.rk <- system.time(parSapply(cl,rep(100,5000),sde.rk))
stopCluster(cl)

cl <- makeCluster(detectCores())
clusterSetRNGStream(cl, iseed=seed)
clusterExport(cl, sde.sim)
time.par.mil <- system.time(parSapply(cl,rep(100,5000),par.fun.mil))
stopCluster(cl)




stopCluster(cl)