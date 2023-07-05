library(sde)
library(ggplot2)
library(dplyr)

seed <- 20230331
set.seed(seed)

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

### Changing mu
T = 1
N = 252 * T
delta = T/N
mu = 0.03

sigma.lower = 0.5
sigma.upper = 3
sigma.delta = 0.25

sigma.list <- seq(sigma.lower, sigma.upper, sigma.delta)

df.err = data.frame(mu=numeric(0), method=character(0), err=numeric(0))

for (sigma in sigma.list) {
  d <- parse(text = paste(mu,"* x"))
  s <- parse(text = paste(sigma,"* x"))

  sim.true <- sde.bs.true(T=T, delta.t=T/N, seed=seed, theta=c(mu,sigma))
  sim.rk <- sde.rk(T=T, delta.t = T/N, seed=seed, theta=c(mu,sigma))
  set.seed(seed)
  sim.mil <- as.numeric(sde.sim(X0=100, drift=d, sigma=s,
                     method="milstein", N=N, T=T, delta = T/N))
  set.seed(seed)
  sim.mil2 <- as.numeric(sde.sim(X0=100, drift=d, sigma=s, 
                      method="milstein2", N=N, T=T, delta = T/N))
  set.seed(seed)
  sim.euler <- as.numeric(sde.sim(X0=100, drift=d, sigma=s,
                       method="euler", N=N, T=T, delta = T/N))
  
  df.mil = data.frame(sigma=rep(sigma,N+1),method=rep("mil",N+1),err=abs(sim.mil-sim.true))
  df.mil2 = data.frame(sigma=rep(sigma,N+1),method=rep("mil2",N+1),err=abs(sim.mil2-sim.true))
  df.euler = data.frame(sigma=rep(sigma,N+1),method=rep("euler",N+1),err=abs(sim.euler-sim.true))
  df.rk = data.frame(sigma=rep(sigma,N+1),method=rep("rk",N+1),err=abs(sim.rk-sim.true))
  
  df.err = rbind(df.err, df.mil, df.mil2, df.euler, df.rk)
  
}

p <-df.err %>%
  mutate( bin=cut_width(sigma, width=sigma.delta, boundary=0) ) %>%
  ggplot( aes(x=bin, y=err, fill=method) ) +
    geom_boxplot() +
    xlab("Volatility") +
    ylab("Absolute Error") +
    stat_summary(fun = median,
               geom = "line",
               aes(group = method, color = method)) +
    scale_fill_discrete(labels=c('Euler Approximation', 'Milstein Scheme',
                                 'Second Milstein Scheme', 'Runge-Kutta Method')) +
    scale_color_discrete(labels=c('Euler Approximation', 'Milstein Scheme',
                                'Second Milstein Scheme', 'Runge-Kutta Method')) +
    guides(fill=guide_legend(title="Method"), color=guide_legend(title="Method")) +
    scale_x_discrete(labels = sigma.list) +
    ggtitle("Varying Volatility (50% & Above)")
p
