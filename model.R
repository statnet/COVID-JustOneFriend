
## COVID-19 household low-degree connections simulator
## Authors: Steven M. Goodreau, Samuel M. Jenness, Martina Morris
## Utilizing EpiModel-Gallery SEIR Model by 
##            Samuel M. Jenness, Venkata R. Duvvuri
## Date: March 2020

## Load EpiModel
suppressMessages(library(EpiModel))
suppressMessages(library(ndtv))

rm(list = ls())
#eval(parse(text = print(commandArgs(TRUE)[1])))
source("module-fx.R")


nsims <- 1
ncores <- 1
nsteps <- 100  # Make sure divisible by 5
n <- 100


# Network model estimation ------------------------------------------------

# Set up features common across all models
nw <- network.initialize(n, directed = FALSE)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 100)
nw %v% 'status' <- rep(c('s','e','c','i','r'), time=n/5)  #  Just for estimation, not simulation
  
# Set up specific model features
formation.1 = ~edges
meandeg.1 <- 4
target.stats.1 <- c(meandeg.1*n/2)

formation.2 = ~edges + nodefactor('status', levels=c('e','c','i','r'))
meandeg.by.status.2 <- c(4, 4, 4, 2, 4)
target.stats.2 <- c(mean(meandeg.by.status.2)*n/2, meandeg.by.status.2[2:5]*n/2/5)

est.1 <- netest(nw, formation.1, target.stats.1, coef.diss)
est.2 <- netest(nw, formation.2, target.stats.2, coef.diss)

# Model diagnostics
dx.1 <- netdx(est.1, nsims = 2, ncores = 2, nsteps = 500)
print(dx.1)
plot(dx.1, plots.joined = FALSE, qnts.alpha = 0.8)

dx.2 <- netdx(est.2, nsims = 2, ncores = 2, nsteps = 500)
print(dx.2)
plot(dx.2, plots.joined = FALSE, qnts.alpha = 0.8)


# Epidemic model simulation -----------------------------------------------

# Model parameters, initial conditions and control variables
param <- param.net(inf.prob = 0.1, act.rate = 1,
                   ec.rate = 0.5, ci.rate = 0.1, ir.rate = 0.1)
init <- init.net(i.num = 5)
control <- control.net(nsteps = nsteps,
                       nsims = nsims,
                       ncores = ncores,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       recovery.FUN = NULL)

# Run the network model simulation with netsim
sim.1 <- netsim(est.1, param, init, control)
print(sim.1)

sim.2 <- netsim(est.2, param, init, control)
print(sim.2)




# Plot outcomes
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim.1,
     mean.col = 1:5, mean.lwd = 1, mean.smooth = FALSE,
     qnts = 1, qnts.col = 1:5, qnts.alpha = 0.25, qnts.smooth = FALSE,
     legend = TRUE)

plot(sim.1, y = c("se.flow", "ec.flow", "ci.flow", "ir.flow"),
     mean.col = 1:4, mean.lwd = 1, mean.smooth = FALSE,
     qnts.col = 1:4, qnts.alpha = 0.25, qnts.smooth = TRUE,
     ylim = c(0, 3), legend = TRUE)








































####################################### Make movies

nw.1 <- get_network(sim.1)

nw.1 <- color_tea_SECIR(nw.1, verbose = FALSE)

slice.par <- list(start = 1, end = 10, interval = 1, 
                  aggregate.dur = 1, rule = "any")
render.par <- list(tween.frames = 10, show.time = FALSE)
plot.par <- list(mar = c(0, 0, 0, 0))

nwanim.1 <- compute.animation(nw.1, slice.par = slice.par, verbose = TRUE)

render.d3movie(
  nw.1,
  render.par = render.par,
  plot.par = plot.par,
  vertex.col = "ndtvcol",
  edge.col = "darkgrey",
  vertex.border = "lightgrey",
  displaylabels = FALSE,
  filename = paste0(getwd(), "/movie.html"))
