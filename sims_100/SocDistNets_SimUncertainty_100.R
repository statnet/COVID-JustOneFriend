
library(statnet)
library(sna)

# Number of sims ----------------------------------
nsim <- 100

# Matrix for results ------------------------------
nets.to.loop <- c('Pre-Covid', 
               'Essential Only', 'Just One Friend', 
               'Just One Friend Per HH')
nets.to.loop.ids <- c('precov', 'essl', 'comb.1', 'comb.2')
empty.results.mat <- replicate(length(nets.to.loop), rep(NA,nsim))
colnames(empty.results.mat) <- nets.to.loop.ids
largest <- geo.3 <- geo.6 <- empty.results.mat

# Model Fits --------------------------------------
set.seed(0)
n <- 200
# Pre-Covid
emptynet <- network.initialize(n, directed=FALSE)
meandeg.precov <- 15
fit.precov <- ergm(emptynet~edges, target.stats=meandeg.precov*n/2)
# Essl
meandeg.essl <- 4
prop.essl <- 0.1
essl <- rep(0, n)
essl[sample(1:n, round(prop.essl*n,0), replace=FALSE)] <- 1
set.vertex.attribute(emptynet, 'essl', essl)
fit.essl <- ergm(emptynet~edges+ nodematch("essl", diff=TRUE, levels=1),
              target.stats= c(
                meandeg.essl*prop.essl*n, 0),
              control=control.ergm(MCMC.burnin = 1e6))
# non.essl.1
meandeg.non.essl.1 <- 2
fit.non.essl.1 <- ergm(emptynet~edges, target.stats= meandeg.non.essl.1*n/2)
# non.essl.2
meandeg.non.essl.2 <- 0.999
fit.non.essl.2 <- ergm(emptynet~edges+ concurrent, target.stats= c(meandeg.non.essl.2*n/2, 0))

# Chunk 3 - define geo.sep
geo.sep <- function(net, netname, distances=c(3,6)) {
  # Geodesic matrix - if not reachable, use NA
  geo.mat <- sna:::geodist(net, inf.replace=NA)$gdist
  # How many households are separated by distance x or less?
  # Do a node-based, not dyad-based, mean
  means <- sapply(distances, function(x) {
    geo.count <- (geo.mat<=x)
    # geo.count[lower.tri(geo.count, diag=TRUE)] <- NA
    diag(geo.count) <- NA
    return(mean(colSums(geo.count, na.rm=TRUE)))
  })
    
    return(means)
}

# Code Loop ---------------------------------------
# Have copied and pasted code from the aggregated chunks, 
# removing the plotting and kpath calculations
start <- Sys.time()
for (i in 1:nsim) {
  
   # Chunk 2 - Pre-COVID
  net.precov <- simulate(fit.precov)
  largcomp.precov <- sum(component.largest(net.precov))
  
  
  # Chunk 4 - compute geo.precov
  geo.precov <-  geo.sep(net.precov)
  
  # Chunk 5 - plot emptynet
  
  # Chunk 6 - essl.net 
  net.essl <- simulate(fit.essl, control = control.simulate(MCMC.burnin = 1e6))
  net.essl.el <- as.edgelist(net.essl)
  largcomp.essl <- sum(component.largest(net.essl))
  geo.essl <-  geo.sep(net.essl)

  # Chunk 7 - non.essl.net.1
  net.non.essl.1 <- simulate(fit.non.essl.1)
  net.comb.1 <- net.non.essl.1
  net.comb.1 <- add.edges(net.comb.1, tail=net.essl.el[,1], head=net.essl.el[,2])
  largcomp.comb.1 <- sum(component.largest(net.comb.1))
  geo.comb.1 <-  geo.sep(net.comb.1)
  
  # Chunk 8 - non.essl.net.2
  net.non.essl.2 <- simulate(fit.non.essl.2, control=control.simulate(MCMC.burnin=1e8))
  net.comb.2 <- net.non.essl.2
  net.comb.2 <- add.edges(net.comb.2, tail=net.essl.el[,1], head=net.essl.el[,2])
  largcomp.comb.2 <- sum(component.largest(net.comb.2))
  geo.comb.2 <-  geo.sep(net.comb.2) 
  
  # Fill results matrices: largest
  largest[i,'precov'] <- largcomp.precov
  largest[i,'essl'] <- largcomp.essl
  largest[i,'comb.1'] <- largcomp.comb.1
  largest[i,'comb.2'] <- largcomp.comb.2
    
  # Fill results matrix: geo.3
  geo.3[i,'precov'] <- geo.precov[1]
  geo.3[i,'essl'] <- geo.essl[1]
  geo.3[i,'comb.1'] <- geo.comb.1[1]
  geo.3[i,'comb.2'] <- geo.comb.2[1]
    
  # Fill results matrix: geo.6
  geo.6[i,'precov'] <- geo.precov[2]
  geo.6[i,'essl'] <- geo.essl[2]
  geo.6[i,'comb.1'] <- geo.comb.1[2]
  geo.6[i,'comb.2'] <- geo.comb.2[2]
  
  if (i==1) append.logical <- FALSE else append.logical <- TRUE
  cat('\n',i,'....', as.character(Sys.time()), 
      file='SocDistNets_SimUncertainty_Progress.txt', append=append.logical)
}

end <- Sys.time()

end-start

write.csv(largest, file=paste0('SocDistNets_SimUncertainty_',
                                nsim, '_largest.csv', row.names=FALSE))
write.csv(geo.3, file=paste0('SocDistNets_SimUncertainty_',
                                nsim, '_geo.3.csv', row.names=FALSE))
write.csv(geo.6, file=paste0('SocDistNets_SimUncertainty_',
                              nsim, '_geo.6.csv', row.names=FALSE))
