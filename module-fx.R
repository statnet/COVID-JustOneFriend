
##
## COVID-19 household low-degree connections simulator
## Authors: Steven M. Goodreau, Samuel M. Jenness, Martina Morris
## Incorporating code from EpiModel-Gallery SEIR Model by Samuel M. Jenness, Venkata R. Duvvuri
## Date: March 2020


##
## SEIR Model: Adding an Exposed State to an SIR
## EpiModel Gallery (https://github.com/statnet/EpiModel-Gallery)
##
## Authors: Samuel M. Jenness, Venkata R. Duvvuri
## Date: August 2018
##

# SECIR
# S = susceptible
# E = exposed: infected, not infectious, not symptomatic
# C = incubating: infectious, not symptomatic
# I = infectious and symptomatic
# R = recovered


# Replacement infection/transmission module -------------------------------

infect <- function(dat, at) {
  
  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  
  ## Parameters ##
  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate
  
  ## Find infected nodes ##
  idsInf <- which(active == 1 & status %in% c("c", "i"))
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  ## Initialize default incidence at 0 ##
  nInf <- 0
  
  ## If any infected nodes, proceed with transmission ##
  if (nElig > 0 && nElig < nActive) {
    
    ## Look up discordant edgelist ##
    del <- discord_edgelist_SECIR(dat, at)
    
    ## If any discordant pairs, proceed ##
    if (!(is.null(del))) {
      
      # Set parameters on discordant edgelist data frame
      del$transProb <- inf.prob
      del$actRate <- act.rate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      
      # Stochastic transmission process
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      
      # Keep rows where transmission occurred
      del <- del[which(transmit == 1), ]
      
      # Look up new ids if any transmissions occurred
      idsNewInf <- unique(del$sus)
      nInf <- length(idsNewInf)
      
      # Set new attributes for those newly infected
      if (nInf > 0) {
        dat$attr$status[idsNewInf] <- "e"
        dat$attr$infTime[idsNewInf] <- at
      }
    }
  }
  
  ## Save summary statistic for S->E flow
  dat$epi$se.flow[at] <- nInf
  
  return(dat)
}

# New disease progression module ------------------------------------------
# (Replaces the recovery module)

progress <- function(dat, at) {

  ## Attributes ##
  active <- dat$attr$active
  status <- dat$attr$status
  
  ## Parameters ##
  ec.rate <- dat$param$ec.rate
  ci.rate <- dat$param$ci.rate
  ir.rate <- dat$param$ir.rate
  
  ## E to C progression process ##
  nInc <- 0
  idsEligInc <- which(active == 1 & status == "e")
  nEligInc <- length(idsEligInc)
  
  if (nEligInc > 0) {
    vecInc <- which(rbinom(nEligInc, 1, ec.rate) == 1)
    if (length(vecInc) > 0) {
      idsInc <- idsEligInc[vecInc]
      nInc <- length(idsInc)
      status[idsInc] <- "c"
    }
  }
  
  ## C to I progression process ##
  nInf <- 0
  idsEligInf <- which(active == 1 & status == "c")
  nEligInf <- length(idsEligInf)
  
  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ci.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }

  ## I to R progression process ##
  nRec <- 0
  idsEligRec <- which(active == 1 & status == "i")
  nEligRec <- length(idsEligRec)
  
  if (nEligRec > 0) {
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1)
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
    }
  }
  
  ## Write out updated status attribute ##
  dat$attr$status <- status
  
  ## Save summary statistics ##
  dat$epi$ec.flow[at] <- nInc
  dat$epi$ci.flow[at] <- nInf
  dat$epi$ir.flow[at] <- nRec
  #dat$epi$s.num[at] <- sum(active == 1 & status == "s")
  dat$epi$e.num[at] <- sum(active == 1 & status == "e")
  dat$epi$c.num[at] <- sum(active == 1 & status == "c")
  #dat$epi$i.num[at] <- sum(active == 1 & status == "i")
  dat$epi$r.num[at] <- sum(active == 1 & status == "r")
  
  return(dat)
}


# Create Discordant edgelist for the SECIR model -------------------------
# (Replaces the discord_esgelist function)


discord_edgelist_SECIR <- function (dat, at) {
  status <- dat$attr$status
  el <- get.dyads.active(dat$nw, at = at)
  del <- NULL
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
    stat <- matrix(status[el], ncol = 2)
    isInf <- matrix(stat %in% c("c","i"), ncol = 2)
    isSus <- matrix(stat %in% "s", ncol = 2)
    SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
    ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
    pairs <- rbind(SIpairs, ISpairs[, 2:1])
    if (nrow(pairs) > 0) {
      sus <- pairs[, 1]
      inf <- pairs[, 2]
      del <- data.frame(at, sus, inf)
    }
  }
  return(del)
}


# Re-color nodes for network movie -------------------------
# (Replaces the color_tea function)

color_tea_SECIR <- function (nd, old.var = "testatus", 
                             old.sus = "s", old.exp = "e", old.inc = "c", old.inf = "i", old.rec = "r", 
                             new.var = "ndtvcol", 
                             new.sus, new.exp, new.inc, new.inf, new.rec, 
                             verbose = TRUE) 
{
  if (missing(new.sus)) {
    new.sus <- transco("steelblue", 0.75)
  }
  if (missing(new.exp)) {
    new.exp <- transco("black", 0.75)
  }
  if (missing(new.inc)) {
    new.inc <- transco("brown", 0.75)
  }
  if (missing(new.inf)) {
    new.inf <- transco("firebrick", 0.75)
  }
  if (missing(new.rec)) {
    new.rec <- transco("seagreen", 0.75)
  }
  
  times <- 1:max(get.change.times(nd))
  for (at in times) {
    stat <- get.vertex.attribute.active(nd, old.var, at = at)
    susceptible <- which(stat == old.sus)
    exposed <- which(stat == old.exp)
    incubating <- which(stat == old.inc)
    infectious.symptomatic <- which(stat == old.inf)
    recovered <- which(stat == old.rec)
    nd <- activate.vertex.attribute(nd, prefix = new.var, 
                                    value = new.sus, onset = at, terminus = Inf, v = susceptible)
    nd <- activate.vertex.attribute(nd, prefix = new.var, 
                                    value = new.exp, onset = at, terminus = Inf, v = exposed)
    nd <- activate.vertex.attribute(nd, prefix = new.var, 
                                    value = new.inc, onset = at, terminus = Inf, v = incubating)
    nd <- activate.vertex.attribute(nd, prefix = new.var, 
                                    value = new.inf, onset = at, terminus = Inf, v = infectious.symptomatic)
    nd <- activate.vertex.attribute(nd, prefix = new.var, 
                                    value = new.rec, onset = at, terminus = Inf, v = recovered)
    if (verbose == TRUE) {
      cat("\n", at, "/", max(times), "\t", 
          sep = "")
    }
  }
  return(nd)
}

