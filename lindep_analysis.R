
################################
################################
#### Code written by Maarten J. Bijlsma
#### Max Planck Institute for Demographic Research
#### Maarten.Bijlsma@Gmail.com
#### Last updated: 2020/Dec/14
#### For the Project: Linear Dependency Issues
#### with Fixed Effects in Family and Sibling models
#### Free to use / share but only with due reference to the Author
#### If intended to be used in scientific work: please check with the Author
#### to settle intellectual property rights
################################
################################

## Analysis for the paper on linear dependency issues
# that emerge once you adjust for family fixed effects

# source the functions that were written for this
source('lindep_functions.R')

# load library for analysis with fixed effect intercepts
library(plm)

# load library for clustered resampling
library(cfdecomp)

####### generate data and analyze it many times ########
####### i.e. do a simulation study              #######

# generate some initial data
# generate data
# DGP1: effect.fambc = 0; effect.bc = 0; effect.bc.sex = 0
DGP.1 <- DGP(n=30000,sibsize=2,sibspacing=2,
             bcstart=1940,bcend=1980,
             famaycenter=1960,famaysd=2,effect.fambc=0,
             famfemean=1,famfesd=1,
             effect.aar=1,effect.bc=0,
             effect.sex=1,effect.bc.sex=0,
             effect.famfe=1,
             yint=0,ysd=1,
             seed=1234)

# DGP2: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = 0
DGP.2 <- DGP(n=30000,sibsize=2,sibspacing=2,
             bcstart=1940,bcend=1980,
             famaycenter=1955,famaysd=2,effect.fambc=0.5,
             famfemean=1,famfesd=1,
             effect.aar=1,effect.bc=1,
             effect.sex=1,effect.bc.sex=0,
             effect.famfe=1,
             yint=0,ysd=1,
             seed=1234)

# DGP3: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = nonzero
DGP.3 <- DGP(n=30000,sibsize=2,sibspacing=2,
             bcstart=1940,bcend=1980,
             famaycenter=1955,famaysd=2,effect.fambc=0.5,
             famfemean=1,famfesd=1,
             effect.aar=1,effect.bc=1,
             effect.sex=1,effect.bc.sex=1,
             effect.famfe=1,
             yint=0,ysd=1,
             seed=1234)

# how many simulations?
it.size <- 999

# set birth cohort groups
min(DGP.1$bc,DGP.2$bc,DGP.3$bc);max(DGP.1$bc,DGP.2$bc,DGP.3$bc)
bc.groups <- seq(1940,1986,5) # I choose 5-year categories

# set siblingsize here already
sibsize <- 2
clustersizemod <- 1/1:sibsize

# make arrays in which to save the results
results.array <- rep(NA,3*5*it.size)
dim(results.array) <- c(3,5,it.size)
results.a.array <- rep(NA,3*5*it.size)
dim(results.a.array) <- c(3,5,it.size)
results.b.array <- rep(NA,3*5*it.size)
dim(results.b.array) <- c(3,5,it.size)

results.a.Ghom.array <- rep(NA,3*5*it.size)
dim(results.a.Ghom.array) <- c(3,5,it.size)
results.a.Gmix.array <- rep(NA,3*5*it.size)
dim(results.a.Gmix.array) <- c(3,5,it.size)

t1 <- Sys.time()
for(i in 1:it.size) {
  print(i)
  
  ## generate data
  # DGP1: effect.fambc = 0; effect.bc = 0; effect.bc.sex = 0
  DGP.1 <- DGP(n=30000,sibsize=sibsize,sibspacing=2,
               bcstart=1940,bcend=1980,
               famaycenter=1960,famaysd=2,effect.fambc=0,
               famfemean=1,famfesd=1,
               effect.aar=1,effect.bc=0,
               effect.sex=1,effect.bc.sex=0,
               effect.famfe=1,
               yint=0,ysd=1,
               seed=i*(1234-i))
  
  # DGP2: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = 0
  DGP.2 <- DGP(n=30000,sibsize=sibsize,sibspacing=2,
               bcstart=1940,bcend=1980,
               famaycenter=1960,famaysd=2,effect.fambc=0.5,
               famfemean=1,famfesd=1,
               effect.aar=1,effect.bc=1,
               effect.sex=1,effect.bc.sex=0,
               effect.famfe=1,
               yint=0,ysd=1,
               seed=i*(1234-i))
  
  # DGP3: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = nonzero
  DGP.3 <- DGP(n=30000,sibsize=sibsize,sibspacing=2,
               bcstart=1940,bcend=1980,
               famaycenter=1960,famaysd=2,effect.fambc=0.5,
               famfemean=1,famfesd=1,
               effect.aar=1,effect.bc=1,
               effect.sex=1,effect.bc.sex=1,
               effect.famfe=1,
               yint=0,ysd=1,
               seed=i*(1234-i))
  
  DGP.1a <- subset(DGP.1, aar<= 18) # removes all those with an age at arrival greater than 18
  DGP.1b <- subset(DGP.1, aar<= 18 & aar>0) # removes those with a zero (those concordant on the exposure)
  DGP.2a <- subset(DGP.2, aar<= 18) # removes all those with an age at arrival greater than 18
  DGP.2b <- subset(DGP.2, aar<= 18 & aar>0) # removes those with a zero (those concordant on the exposure)
  DGP.3a <- subset(DGP.3, aar<= 18) # removes all those with an age at arrival greater than 18
  DGP.3b <- subset(DGP.3, aar<= 18 & aar>0) # removes those with a zero (those concordant on the exposure)
  
  # further subsample DGP.1a, to determine which family type drives identification
  X <- gensplit(DGP.1a); DGP.1a.Gmix <- X[[1]]; DGP.1a.Ghom <- X[[2]]; rm(X)
  X <- gensplit(DGP.2a); DGP.2a.Gmix <- X[[1]]; DGP.2a.Ghom <- X[[2]]; rm(X)
  X <- gensplit(DGP.3a); DGP.3a.Gmix <- X[[1]]; DGP.3a.Ghom <- X[[2]]; rm(X)
  
  # sample size reduction is done so that all samples have the same size
  # due to MC error reduction this will have no bearing on the findings
  # but co-authors feel that this is necessary
  DGP.1 <- cluster.reduce.sample(DGP.1,10000,"famid",1/sibsize)
  DGP.2 <- cluster.reduce.sample(DGP.2,10000,"famid",1/sibsize)
  DGP.3 <- cluster.reduce.sample(DGP.3,10000,"famid",1/sibsize) 
  
  DGP.1a <- cluster.reduce.sample(DGP.1a,10000,"famid",clustersizemod)
  DGP.2a <- cluster.reduce.sample(DGP.2a,10000,"famid",clustersizemod)
  DGP.3a <- cluster.reduce.sample(DGP.3a,10000,"famid",clustersizemod)
  
  DGP.1b <- cluster.reduce.sample(DGP.1b,10000,"famid",clustersizemod)
  DGP.2b <- cluster.reduce.sample(DGP.2b,10000,"famid",clustersizemod)
  DGP.3b <- cluster.reduce.sample(DGP.3b,10000,"famid",clustersizemod)
  
  # also for the hom group in a.G (the mix group is too small to start with)
  DGP.1a.Ghom <- cluster.reduce.sample(DGP.1a.Ghom,10000,"famid",clustersizemod)
  DGP.2a.Ghom <- cluster.reduce.sample(DGP.2a.Ghom,10000,"famid",clustersizemod)
  DGP.3a.Ghom <- cluster.reduce.sample(DGP.3a.Ghom,10000,"famid",clustersizemod)
  
  # analyse the DGPs and save results in an array
  results.array[,,i] <-   analyse.DGP(DGP.1 ,DGP.2 ,DGP.3 ,bc.groups)[[1]]
  results.a.array[,,i] <- analyse.DGP(DGP.1a,DGP.2a,DGP.3a,bc.groups)[[1]]
  results.b.array[,,i] <- analyse.DGP(DGP.1b,DGP.2b,DGP.3b,bc.groups)[[1]]
  
  results.a.Ghom.array[,,i] <- analyse.DGP(DGP.1a.Ghom,DGP.2a.Ghom,DGP.3a.Ghom,bc.groups)[[1]]
  results.a.Gmix.array[,,i] <- analyse.DGP(DGP.1a.Gmix,DGP.2a.Gmix,DGP.3a.Gmix,bc.groups)[[1]]
}
t2 <- Sys.time()
t2 - t1

# results averaged over the simulations:
DGP123 <- apply(results.array,c(1,2),mean)
DGP123a <- apply(results.a.array,c(1,2),mean)
DGP123b <- apply(results.b.array,c(1,2),mean)

DGP123a.Ghom <- apply(results.a.Ghom.array,c(1,2),mean)
DGP123a.Gmix <- apply(results.a.Gmix.array,c(1,2),mean)

# test the effect of different ratios of discordant on bc but concordant on aar
# generate a very large pool of data
# sibsize has to be 2, otherwise the filtering function (below) fails
DGP.2.pool <- DGP(n=1000000,sibsize=2,sibspacing=2,
                  bcstart=1940,bcend=1980,
                  famaycenter=1955,famaysd=2,effect.fambc=0.5,
                  famfemean=1,famfesd=1,
                  effect.aar=1,effect.bc=1,
                  effect.sex=1,effect.bc.sex=0,
                  effect.famfe=1,
                  yint=0,ysd=1,
                  seed=1234)

# let's take as percentages of discon: 0.05, 0.10, 0.25, 0.5, 0.85
min(DGP.2.pool$bc);max(DGP.2.pool$bc)
bc.groups <- seq(1940,1986,5) # I choose 5-year categories
it.size <- 999
results.array <- rep(NA,it.size*5*5) # it.size, 5 analytical model results, 5 percentages of discon
dim(results.array) <- c(it.size,5,5) # it.size, 5 analytical model results, 5 percentages of discon
for(d in 1:5) {
  
  print(d)
  
  p <- c(0.05,0.10,0.25,0.5,0.85)[d]
  
  for(i in 1:it.size) {
    
    print(i)
    
    e <- 1
    while(e <= 1) {
      DGP.2 <- disconfilter(DGP.2.pool,n=5000,perc.discon=p)
      # 5000 is a sample size of 10.000 because sibsize here is kept at 2
      e <- sum(table(table(DGP.2$aar)))
    }
    # the above while loop only exists in case it mysteriously generates only aar's of 0
    
    results.array[i,,d] <- analyse.DGP.2(DGP.2,bc.groups)
    
  }
}

# results averaged over the re-sampled data
apply(results.array,c(2,3),mean)


