library(sde)
library(ggplot2)
library(dplyr)

seed <- 20230331
set.seed(seed)
d <- expression(0.05 * x)
s <- expression(0.25 * x)
s.x <- expression(0.25)
d.x <- expression(0.05)


seed.list <- sample(1:10000, 1000, replace=TRUE)
delta.ts <- c(0.001, 0.01, 0.05, 0.1)
mus <- c(0.05, 0)
sigmas <- c(0.5, 0.25, 0.1, 0.05)

# Black Scholes True solution function
sde.bs.true <- function(X0=100, theta=c(0.05, 0.25), delta.t=0.004, T=1, seed=0){
  set.seed(seed)
  t <- seq(0, T, delta.t)
  dW <- c(0)
  dW <- append(dW, rnorm(T/delta.t, 0, sqrt(delta.t)))
  B <- cumsum(dW)
  return (X0*exp((theta[1]-(theta[2]**2)/2)*t + theta[2]*B))
}

# Runge Kutta function
sde.rk <- function(X0=100, theta=c(0.05, 0.25), delta.t=0.004, T=1, seed=0){
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


# Simple Checking - No parameter change
T = 0.75
N = 252 * T
delta = T/N
sim.true <- sde.bs.true(T=T, delta.t = T/N, seed=seed)
sim.rk <- sde.rk(T=T, delta.t = T/N, seed=seed)
set.seed(seed)
sim.mil <- sde.sim(X0=100, drift=d, sigma=s, 
                   method="milstein", N=N, T=T, delta = T/N)
set.seed(seed)
sim.mil2 <- sde.sim(X0=100, drift=d, sigma=s, 
                    method="milstein2", N=N, T=T, delta = T/N)
set.seed(seed)
sim.euler <- sde.sim(X0=100, drift=d, sigma=s,
                     method="euler", N=N, T=T, delta = T/N)

df.mil = data.frame(t=c(seq(0,T,T/N)),method=c(rep("mil",N+1)),err=c(abs(sim.mil-sim.true)))
df.mil2 = data.frame(t=c(seq(0,T,T/N)),method=c(rep("mil2",N+1)),err=c(abs(sim.mil2-sim.true)))
df.euler = data.frame(t=c(seq(0,T,T/N)),method=c(rep("euler",N+1)),err=c(abs(sim.euler-sim.true)))
df.rk = data.frame(t=c(seq(0,T,T/N)),method=c(rep("rk",N+1)),err=c(abs(sim.rk-sim.true)))

df.err = rbind(df.mil, df.mil2, df.euler, df.rk)

p <-df.err %>%
  mutate( bin=cut_width(t, width=0.05, boundary=0) ) %>%
  ggplot( aes(x=bin, y=err, fill=method) ) +
    geom_boxplot() +
    xlab("Time") +
    ylab("Absolute Error") +
    stat_summary(fun = median,
               geom = "line",
               aes(group = method, color = method)) +
    scale_fill_discrete(labels=c('Euler Approximation', 'Milstein Scheme',
                                 'Second Milstein Scheme', 'Runge-Kutta Method'))
p
