
################################
################################
#### Code written by Maarten J. Bijlsma
#### Max Planck Institute for Demographic Research
#### Maarten.Bijlsma@Gmail.com
#### Last updated: 2020/Dec/14
#### For the Project: Linear Dependency Issues in
#### with Fixed Effects in Family and Sibling models
#### Free to use / share but only with due reference to the Author
#### If intended to be used in scientific work: please check with the Author
#### to settle intellectual property rights
################################
################################

Sys.setenv(LANG = "en")

## Data Generating Processes (DGPs)
# DGP1 (no birth cohort):  	          Age at arrival + sibling FE + sex
# DGP2 (with birth cohort):           Age at arrival + birth cohort + sibling FE + sex
# DGP3 (differential birth cohort):   Age at arrival + birth cohort + birth cohort * sex + sibling FE + sex

# DGP 4 was planned to be done at a later time (more complicated)
# DGP4:   Age at arrival + birth cohort + age of respondent  + period of measurement t + sibling FE + sex
# was ultimately deemed unnecessary for the goals of our analysis

# I created a function that generates data. E.g. a DGP (data generating process).
# another function will be written that analyses the DGP output following some
# functional form

# parameters:
# n = number of families
# sibsize = number of siblings per family
# sibspacing = mean difference in age between sibling s and s+1, poisson distributed
# bcstart = birth cohort start year for sibling 1
# bcend = birth cohort end year for sibling 1
# famaycenter = in what year on average do families arrive?
# famaysd = what is the standard deviation around the famaycenter?
# effect.fambc = what is the effect of birth cohort on arrival year? default = 0
#       set this to 1 would mean that more recent birth cohorts have a higher chance to arrive more recently
#       which means (age at arrival = year of arrival - birth cohort) that they have a higher age at arrival
#       if kept at 0, it means that there is no causal effect of birth cohort on age at arrival (DGP1)
# famfemean = family fixed effect mean, default is 0
# famfesd = family fixed effect standard deviation
#       the family fixed effects are assumed to be normally distributed
# effect.aar = age at arrival effect on outcome y (single number)
# effect.bc = birth cohort effect on outcome y (single number)
# effect.sex = sex effect on outcome y (single number)
# effect.bc.sex = effect of interaction between birth cohort and sex on outcome y, set to non-zero for DGP3, default is 0
# effect.famfe = family fixef effect on outcome y, default is 1
# yint = intercept of the outcome y, default is 0
# ysd = standard deviation of the residuals of y

# seed = set.seed, i.e. if we want the DGP to generate under the same conditions 
# for scenario 1, 2 and 3, set to the same value
# in a loop, seed should be randomized

# * age of respondent is not included, but in DGPs 1-3 this is e.g. age 30 (education at age 30)
# * period of measurement is not included, can be dependent on other variables

DGP <- function(n,
         sibsize,
         sibspacing,
         bcstart,
         bcend,
         famaycenter,
         famaysd,
         effect.fambc,
         famfemean,
         famfesd,
         effect.aar,
         effect.bc,
         effect.sex,
         effect.bc.sex,
         effect.famfe,
         yint,
         ysd,
         seed) {
  set.seed(seed)
  # generate id
  id <- 1:(n*sibsize) # overall ID
  famid <- rep(1:n,each=sibsize) # generate family identifier
  birthorder <- rep(1:sibsize,n) # generate birth order within families
  # start generate birth cohorts
  bc <- rep(NA,n*sibsize)
  bc[seq(1,n*sibsize,by=sibsize)] <- floor(runif(n,bcstart,bcend-0.0000001)) # generate birth cohort of sibling 1
  for(s in 2:sibsize) {
    bc[seq(s,n*sibsize,by=sibsize)] <- bc[seq(s-1,n*sibsize,by=sibsize)]+rpois(n,sibspacing)
  } # generates birth cohort of remaining siblings by adding years to each successive siblings birth cohort
  # end generate birth cohort
  # start generate age at arrival
  # for birth cohort to affect age at arrival (DGP2), ideally we would say something like:
  # "more recent birth cohorts have a lower age at arrival", but since birth cohorts are generated
  # semi-independently of one another, this would mean that different members of the family could have
  # a different arrival year (since arrival year = age at arrival + birth cohort)
  # if we want families to enter the country at the same time, then we can only create this dependency
  # by making arrival year dependent on mean family birth cohort
  # calculate mean birth cohort within families:
  fambc <- rep(tapply(bc,famid,mean),each=sibsize)
  # if family arrival year (and hence age at arrival) is dependent on birth cohort, then:
  famay <- round(rnorm(n,famaycenter,famaysd) + (fambc[seq(1,n*sibsize,by=sibsize)]-mean(fambc))*effect.fambc)
  famay <- rep(famay,each=sibsize) # give each  individual their family arrival year
  # ! it seems mathematically and conceptually silly to have an effect (causal effect?)
  # of birth cohort on age at arrival, when families arrive at the same time
  # we should probably discuss that in the paper and then perhaps only have this situation
  # when family members can arrive at different years
  # we should also discuss if we don't want this to be normally distributed, but e.g. a uniform distribution
  aar <- bc - famay # age at arrival variable
  aar <- ifelse(aar < 0, 0, aar) # if age at arrival is negative (arrival year is before birth), set it to 0 instead
  # end generate age at arrival
  # start generate family fixed effect
  famfe <- rep(rnorm(n,famfemean,famfesd),each=sibsize)
  # end generate family fixed effct
  # start generate sex
  sex <- rbinom(n*sibsize,1,0.5) # 50% male and 50% female
  # end generate sex
  # start generate outcome Y
  y <- aar*effect.aar + bc*effect.bc + sex*effect.sex + bc*sex*effect.bc.sex + famfe*effect.famfe + rnorm(n,yint,ysd)
  
  output <- data.frame(cbind(id,famid,birthorder,bc,fambc,famay,aar,famfe,sex,y))
  
  return(output)
}

# longitudinal sampling function, useful for when sampling families with replacement
long.sample <- function(originaldata, originaldataid, size) {
  # select a bunch of IDs
  IDs <- unique(originaldataid)
  y <- sample(IDs,size,replace=T)
  z <- table(table(y))
  
  # from there, select a group once
  selectID <- sample(IDs,size=z[1],replace=F)
  newdata <- originaldata[which(originaldataid %in% selectID),]
  
  if(length(z) > 1) {
    
    for(i in 2:length(z)) {
      
      # select a new group of IDs that was not yet selected
      IDs2 <- setdiff(IDs,selectID)
      
      # from there, randomly select a group of people of the right size
      selectID2 <- sample(IDs2,size=z[i],replace=F)
      selectID <- c(selectID,selectID2) # so we don't re-select the newly
      # selected people either
      
      for(j in 1:i) {
        
        # copy the new dataset i number of times
        newdata <- rbind(newdata,originaldata[which(originaldataid %in% selectID2),])
        
      }
      
    }
    
  }
  return(newdata)
}

# DGP.pool is data produced by the DGP function, ideally of a very large size
# n is the number of families you want in the analysis
# perc.discon is the percentage of families that are discordant on birth cohort but concordant on aar
# ! note, this function  only works for DGPs with sibsize of 2
# ! note, this function requires the long.sample function to be loaded

disconfilter <- function(DGP.pool,n,perc.discon) {
  # calculate differences between sibling birth cohorts
  fambcdiff <- tapply(DGP.pool$bc,DGP.pool$famid,diff) # this part is especially tricky with sibsize > 2
  DGP.pool$fambcdiff <- rep(fambcdiff,each=2) # each is based on sibsize
  # throw out those families where birth cohorts are equal (fambcdiff == 0)
  # (these are 'twins')
  DGP.pool <- DGP.pool[DGP.pool$fambcdiff!=0,]
  # identify G1 and G2 individuals
  # individuals are G2 if aar==0, and G1 otherwise
  DGP.pool$G <- ifelse(DGP.pool$aar==0,2,1)
  # identify families that are 2xG2 or not. Since we threw out twins,
  # families that are 2xG2 will be concordant on aar, but discordant on birth cohort
  # and other families will be discordant on both
  # We just check if they have equal AAR
  # if they do, we know they are concordant on AAR but they are discordant on BC (since we threw out twins)
  # which means they are both G2
  # if they do not have equal AAR, then they are discordant on both AAR and BC
  famaardiff <- tapply(DGP.pool$aar,DGP.pool$famid,diff)
  DGP.pool$famaardiff <- rep(famaardiff,each=2)
  # now we can select families based on their discordant/concordant status
  n.discon <- n*perc.discon
  n.disdis <- n*(1-perc.discon)
  DGP.discon <- long.sample(DGP.pool[DGP.pool$famaardiff==0,],DGP.pool$famid[DGP.pool$famaardiff==0],size=n.discon)
  DGP.disdis <- long.sample(DGP.pool[DGP.pool$famaardiff!=0,],DGP.pool$famid[DGP.pool$famaardiff!=0],size=n.disdis)
  DGP.new <- rbind(DGP.discon,DGP.disdis)
  return(DGP.new)
}

# generate data
# DGP1: effect.fambc = 0; effect.bc = 0; effect.bc.sex = 0
DGP.1 <- DGP(n=20000,sibsize=3,sibspacing=2,
             bcstart=1940,bcend=1980,
             famaycenter=1955,famaysd=2,effect.fambc=0,
             famfemean=1,famfesd=1,
             effect.aar=1,effect.bc=0,
             effect.sex=1,effect.bc.sex=0,
             effect.famfe=1,
             yint=0,ysd=1,
             seed=1234)

# DGP2: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = 0
DGP.2 <- DGP(n=20000,sibsize=3,sibspacing=2,
             bcstart=1940,bcend=1980,
             famaycenter=1955,famaysd=2,effect.fambc=0.5,
             famfemean=1,famfesd=1,
             effect.aar=1,effect.bc=1,
             effect.sex=1,effect.bc.sex=0,
             effect.famfe=1,
             yint=0,ysd=1,
             seed=1234)

# DGP3: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = nonzero
DGP.3 <- DGP(n=20000,sibsize=3,sibspacing=2,
             bcstart=1940,bcend=1980,
             famaycenter=1955,famaysd=2,effect.fambc=0.5,
             famfemean=1,famfesd=1,
             effect.aar=1,effect.bc=1,
             effect.sex=1,effect.bc.sex=1,
             effect.famfe=1,
             yint=0,ysd=1,
             seed=1234)


# code to create subsets of the data

DGP.1a <- subset(DGP.1, aar<= 18) # removes all those with an age at arrival greater than 18
DGP.1b <- subset(DGP.1, aar<= 18 & aar>0) # removes those with a zero (those concordant on the exposure)

table1 <- table(DGP.1$aar)
table1
table1a <- table(DGP.1a$aar)
table1a
table1b <- table(DGP.1b$aar)
table1b



DGP.2a <- subset(DGP.2, aar<= 18) # removes all those with an age at arrival greater than 18
DGP.2b <- subset(DGP.2, aar<= 18 & aar>0) # removes those with a zero (those concordant on the exposure)

table2 <- table(DGP.2$aar)
table2
table2a <- table(DGP.2a$aar)
table2a
table2b <- table(DGP.2b$aar)
table2b

DGP.3a <- subset(DGP.3, aar<= 18) # removes all those with an age at arrival greater than 18
DGP.3b <- subset(DGP.3, aar<= 18 & aar>0) # removes those with a zero (those concordant on the exposure)

table3 <- table(DGP.3$aar)
table3
table3a <- table(DGP.3a$aar)
table3a
table3b <- table(DGP.3b$aar)
table3b

# we want all three groups to have the same size
# not super important (since we don't investigate standard errors)
# but it might be neater - note, however, that this selecting may also
# impact results. To make sure that the effect of selection is minimized
# I check the empirical distribution of families within each subsample
# and then make sure that this empirical distribution remains
# I will reduce the sample size down to the same sample size for all
# by doing stratified clustered sampling within DGPs. The relevant strata here are
# famsize/sibsize. The clusters are families.
# I will make a function for this

cluster.reduce.sample <- function(DGP,n.reduced,clustername,clustersizemod) {
  # DGP is the dataset, could be produced using the DGP function
  # n.reduced is the new size of the dataset we want
  # clustername is the name that we want to do stratified resampling over
  # clustersizemod is whether we want the size of each group 
  # (now sampled proportionally) to instead be modified somehow
  
  tablecluster <- table(DGP[,clustername])
  proptabletablecluster <- prop.table(table(tablecluster))
  cluster.sizes <- n.reduced * proptabletablecluster
  DGP.new <- NULL
  for(i in 1:length(names(proptabletablecluster))) {
    
    k <- names(proptabletablecluster)[i]
    
    cluster.id <- names(tablecluster)[tablecluster==k]
    # sample with appropriate size (and without replacement):
    reduced.cluster.id <- sample(cluster.id,size=round(cluster.sizes[names(proptabletablecluster)==k]*clustersizemod[i]))
    DGPcluster <- DGP[which(DGP[,clustername] %in% reduced.cluster.id),]
    
    # we don't need to sample with replacement, but I'll keep the code below:
    # long.sample(DGPcluster,DGPcluster[,clustername],size=cluster.sizes[names(proptabletablecluster)==k])
    
    DGP.new <- rbind(DGP.new,DGPcluster)
  }
  return(DGP.new)
}

clustersizemod <- c(1,1/2,1/3)
test <- cluster.reduce.sample(DGP.1a,10000,"famid",clustersizemod)
# works! Since test has a size of 10.000
rm(test)

# now let's analyze these datasets following different models:

library(plm)

## Analysis
# A1: No controls
# A2: Standardized outcome (now with z-score rather than quantiles)
# A3: linear controls (continuous variables)
# A4: factor controls full (discrete variables)
# A5: factor controls (grouped)

# note that A4 is automatically an APC model and 'adding a reference'
# could add parameter to set the reference manually

# parameters
# DGP.1 = the DGP.1 from the DGP function
# DGP.2 = the DGP.2 from the DGP function
# DGP.3 = the DGP.3 from the DGP function
# bc.groups = the cutpoints for the birth cohort groupings we want
#             do we also want to group age at arrival? That complicates things
#             but is at least not strictly needed to break linear dependency, as long as BC is grouped

analyse.DGP <- function(DGP.1,
                        DGP.2,
                        DGP.3,
                        bc.groups,
                        plot.resid=FALSE) {
  
  # make dataframe in which to put results by DGP (3) and analytical model (5)
  results <- rep(NA,3*5)
  dim(results) <- c(3,5)
  
  # results for sd of residuals
  results.sd <- rep(NA,3*5)
  dim(results.sd) <- c(3,5)
  
  # ! not sure if best solution: some DGPs might produce bcs outside of this range? test
  # ! likely a minor issue though
  DGP.1$bc.grouped <- as.character(cut(DGP.1$bc,breaks=bc.groups))
  DGP.2$bc.grouped <- as.character(cut(DGP.2$bc,breaks=bc.groups))
  DGP.3$bc.grouped <- as.character(cut(DGP.3$bc,breaks=bc.groups))
  
  # determine Z-score by birth cohort (A2):
  for(bc in unique(DGP.1$bc)) {
    mu <- mean(DGP.1$y[DGP.1$bc==bc])
    std <- sd(DGP.1$y[DGP.1$bc==bc])
    DGP.1$y.q[DGP.1$bc==bc] <- (DGP.1$y[DGP.1$bc==bc]-mu)/std
  }
  for(bc in unique(DGP.2$bc)) {
    mu <- mean(DGP.2$y[DGP.2$bc==bc])
    std <- sd(DGP.2$y[DGP.2$bc==bc])
    DGP.2$y.q[DGP.2$bc==bc] <- (DGP.2$y[DGP.2$bc==bc]-mu)/std
  }
  for(bc in unique(DGP.3$bc)) {
    mu <- mean(DGP.3$y[DGP.3$bc==bc])
    std <- sd(DGP.3$y[DGP.3$bc==bc])
    DGP.3$y.q[DGP.3$bc==bc] <- (DGP.3$y[DGP.3$bc==bc]-mu)/std
  }
  
  # I store only the aar parameter, which is what we're substantively interested in
  # analyse DGP.1 data (first line is analytical model A1, second is A2, etc.)
  fit1 <- plm(y ~ aar, effect='individual',model='within',index='famid',family='gaussian',data=DGP.1)
  fit2 <- plm(y.q ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.1)
  fit3 <- plm(y ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.1)
  fit4 <- plm(y ~ aar + as.factor(bc) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.1)
  fit5 <- plm(y ~ aar + as.factor(bc.grouped) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.1)
  
  sdy <- sd(DGP.1$y)
  sdy.q <- sd(DGP.1$y.q)
  
  if(plot.resid==TRUE) {
  plot(ksmooth(DGP.1$aar,fit1$residuals/sdy,kernel="normal",bandwidth=1),type='l',ylab='residuals',xlab='age at arrival',ylim=c(-0.1,0.1),lwd=2,
       main='DGP.1')
  abline(a=0,b=0)
  lines(ksmooth(DGP.1$aar,fit2$residuals/sdy.q,kernel="normal",bandwidth=1),col='blue',lwd=2)
  lines(ksmooth(DGP.1$aar,fit3$residuals/sdy,kernel="normal",bandwidth=1),col='red',lwd=2)
  lines(ksmooth(DGP.1$aar,fit4$residuals/sdy,kernel="normal",bandwidth=1),col='green',lwd=2)
  lines(ksmooth(DGP.1$aar,fit5$residuals/sdy,kernel="normal",bandwidth=1),col='purple',lwd=2)
  PLOT.aar.DGP.1 <- recordPlot()
  
  plot(ksmooth(DGP.1$bc,fit1$residuals/sdy,kernel="normal",bandwidth=1),type='l',ylab='residuals',xlab='birth cohort',ylim=c(-0.1,0.1),lwd=2,
       main='DGP.1')
  abline(a=0,b=0)
  lines(ksmooth(DGP.1$bc,fit2$residuals/sdy.q,kernel="normal",bandwidth=1),col='blue',lwd=2)
  lines(ksmooth(DGP.1$bc,fit3$residuals/sdy,kernel="normal",bandwidth=1),col='red',lwd=2)
  lines(ksmooth(DGP.1$bc,fit4$residuals/sdy,kernel="normal",bandwidth=1),col='green',lwd=2)
  lines(ksmooth(DGP.1$bc,fit5$residuals/sdy,kernel="normal",bandwidth=1),col='purple',lwd=2)
  PLOT.bc.DGP.1 <- recordPlot()}
  
  # ! checking this around the third dimension would also be good
  # ! not that the residuals from fit2 check against y.q rather than against actual y
  # ! this is possibly bad: y.q are on a different scale than actual y
  # ! therefore, to compare the standard deviations of the residuals of the 5 models
  # ! I  normalize using the sd of the simulated y and the simulated y.q

  results[1,1] <- fit1$coef[1]
  results[1,2] <- fit2$coef[1]
  results[1,3] <- fit3$coef[1]
  results[1,4] <- fit4$coef[1]
  results[1,5] <- fit5$coef[1]

  results.sd[1,1] <- sd(fit1$residuals)/sdy
  results.sd[1,2] <- sd(fit2$residuals)/sdy.q
  results.sd[1,3] <- sd(fit3$residuals)/sdy
  results.sd[1,4] <- sd(fit4$residuals)/sdy
  results.sd[1,5] <- sd(fit5$residuals)/sdy
  
  rm(fit1,fit2,fit3,fit4,fit5)
  
  # analyse DGP.2 data (first line is analytical model A1, second is A2, etc.)
  fit1 <- plm(y ~ aar, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)
  fit2 <- plm(y.q ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)
  fit3 <- plm(y ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)
  fit4 <- plm(y ~ aar + as.factor(bc) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)
  fit5 <- plm(y ~ aar + as.factor(bc.grouped) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)
  
  sdy <- sd(DGP.1$y)
  sdy.q <- sd(DGP.1$y.q)
  
  if(plot.resid==TRUE) {
  plot(ksmooth(DGP.2$aar,fit1$residuals/sdy,kernel="normal",bandwidth=1),type='l',ylab='residuals',xlab='age at arrival',ylim=c(-0.1,0.1),lwd=2,
       main='DGP.2')
  abline(a=0,b=0)
  lines(ksmooth(DGP.2$aar,fit2$residuals/sdy.q,kernel="normal",bandwidth=1),col='blue',lwd=2)
  lines(ksmooth(DGP.2$aar,fit3$residuals/sdy,kernel="normal",bandwidth=1),col='red',lwd=2)
  lines(ksmooth(DGP.2$aar,fit4$residuals/sdy,kernel="normal",bandwidth=1),col='green',lwd=2)
  lines(ksmooth(DGP.2$aar,fit5$residuals/sdy,kernel="normal",bandwidth=1),col='purple',lwd=2)
  PLOT.aar.DGP.2 <- recordPlot()
  
  plot(ksmooth(DGP.2$bc,fit1$residuals/sdy,kernel="normal",bandwidth=1),type='l',ylab='residuals',xlab='birth cohort',ylim=c(-0.1,0.1),lwd=2,
       main='DGP.2')
  abline(a=0,b=0)
  lines(ksmooth(DGP.2$bc,fit2$residuals/sdy.q,kernel="normal",bandwidth=1),col='blue',lwd=2)
  lines(ksmooth(DGP.2$bc,fit3$residuals/sdy,kernel="normal",bandwidth=1),col='red',lwd=2)
  lines(ksmooth(DGP.2$bc,fit4$residuals/sdy,kernel="normal",bandwidth=1),col='green',lwd=2)
  lines(ksmooth(DGP.2$bc,fit5$residuals/sdy,kernel="normal",bandwidth=1),col='purple',lwd=2)
  PLOT.bc.DGP.2 <- recordPlot()}
  
  # ! checking this around the third dimension would also be good
  # ! not that the residuals from fit2 check against y.q rather than against actual y
  # ! this is possibly bad: y.q are on a different scale than actual y
  # ! therefore, to compare the standard deviations of the residuals of the 5 models
  # ! I  normalize using the sd of the simulated y and the simulated y.q
  
  results[2,1] <- fit1$coef[1]
  results[2,2] <- fit2$coef[1]
  results[2,3] <- fit3$coef[1]
  results[2,4] <- fit4$coef[1]
  results[2,5] <- fit5$coef[1]

  results.sd[2,1] <- sd(fit1$residuals)/sdy
  results.sd[2,2] <- sd(fit2$residuals)/sdy.q
  results.sd[2,3] <- sd(fit3$residuals)/sdy
  results.sd[2,4] <- sd(fit4$residuals)/sdy
  results.sd[2,5] <- sd(fit5$residuals)/sdy
  
  rm(fit1,fit2,fit3,fit4,fit5)
  
  # analyse DGP.3 data (first line is analytical model A1, second is A2, etc.)
  fit1 <- plm(y ~ aar, effect='individual',model='within',index='famid',family='gaussian',data=DGP.3)
  fit2 <- plm(y.q ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.3)
  fit3 <- plm(y ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.3)
  fit4 <- plm(y ~ aar + as.factor(bc) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.3)
  fit5 <- plm(y ~ aar + as.factor(bc.grouped) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.3)
  
  sdy <- sd(DGP.1$y)
  sdy.q <- sd(DGP.1$y.q)
  
  if(plot.resid==TRUE) {
  plot(ksmooth(DGP.3$aar,fit1$residuals/sdy,kernel="normal",bandwidth=1),type='l',ylab='residuals',xlab='age at arrival',ylim=c(-0.1,0.1),lwd=2,
       main='DGP.3')
  abline(a=0,b=0)
  lines(ksmooth(DGP.3$aar,fit2$residuals/sdy.q,kernel="normal",bandwidth=1),col='blue',lwd=2)
  lines(ksmooth(DGP.3$aar,fit3$residuals/sdy,kernel="normal",bandwidth=1),col='red',lwd=2)
  lines(ksmooth(DGP.3$aar,fit4$residuals/sdy,kernel="normal",bandwidth=1),col='green',lwd=2)
  lines(ksmooth(DGP.3$aar,fit5$residuals/sdy,kernel="normal",bandwidth=1),col='purple',lwd=2)
  PLOT.aar.DGP.3 <- recordPlot()
  
  plot(ksmooth(DGP.3$bc,fit1$residuals/sdy,kernel="normal",bandwidth=1),type='l',ylab='residuals',xlab='birth cohort',ylim=c(-0.1,0.1),lwd=2,
       main='DGP.3')
  abline(a=0,b=0)
  lines(ksmooth(DGP.3$bc,fit2$residuals/sdy.q,kernel="normal",bandwidth=1),col='blue',lwd=2)
  lines(ksmooth(DGP.3$bc,fit3$residuals/sdy,kernel="normal",bandwidth=1),col='red',lwd=2)
  lines(ksmooth(DGP.3$bc,fit4$residuals/sdy,kernel="normal",bandwidth=1),col='green',lwd=2)
  lines(ksmooth(DGP.3$bc,fit5$residuals/sdy,kernel="normal",bandwidth=1),col='purple',lwd=2)
  PLOT.bc.DGP.3 <- recordPlot()}
  
  # ! checking this around the third dimension would also be good
  # ! not that the residuals from fit2 check against y.q rather than against actual y
  # ! this is possibly bad: y.q are on a different scale than actual y
  # ! therefore, to compare the standard deviations of the residuals of the 5 models
  # ! I  normalize using the sd of the simulated y and the simulated y.q
  
  results[3,1] <- fit1$coef[1]
  results[3,2] <- fit2$coef[1]
  results[3,3] <- fit3$coef[1]
  results[3,4] <- fit4$coef[1]
  results[3,5] <- fit5$coef[1]

  results.sd[3,1] <- sd(fit1$residuals)/sdy
  results.sd[3,2] <- sd(fit2$residuals)/sdy.q
  results.sd[3,3] <- sd(fit3$residuals)/sdy
  results.sd[3,4] <- sd(fit4$residuals)/sdy
  results.sd[3,5] <- sd(fit5$residuals)/sdy
  
  rm(fit1,fit2,fit3,fit4,fit5)
  
  if(plot.resid==FALSE) {
    return(list(results,results.sd))
  } else if (plot.resid==TRUE) {
    return(list(results=results,
                results.sd=results.sd,
                PLOT.aar.DGP.1=PLOT.aar.DGP.1,
                PLOT.bc.DGP.1=PLOT.bc.DGP.1,
                PLOT.aar.DGP.2=PLOT.aar.DGP.2,
                PLOT.bc.DGP.2=PLOT.bc.DGP.2,
                PLOT.aar.DGP.3=PLOT.aar.DGP.3,
                PLOT.bc.DGP.3=PLOT.bc.DGP.3))
  }
  
}

# before function, set bc grouping categories
min(DGP.1$bc,DGP.2$bc,DGP.3$bc);max(DGP.1$bc,DGP.2$bc,DGP.3$bc)
bc.groups <- seq(1940,1986,5) # I choose 5-year categories

# run the function
results <-   analyse.DGP(DGP.1 ,DGP.2 ,DGP.3 ,bc.groups) 
results.a <- analyse.DGP(DGP.1a,DGP.2a,DGP.3a,bc.groups)
results.b <- analyse.DGP(DGP.1b,DGP.2b,DGP.3b,bc.groups)
# two sets of results since the second set of results is about SD
# there we see NA for the second analytical strategy (column 2)
# I will ignore that here since we are not interested in SD anymore

# rows are the DGPs, columns are the models used to analyze the data
# the original results tables (gives coefficients)
results[1]
results.a[1]
results.b[1]
# ! the additional results table that gives normalized sd of the residuals
results[2]
results.a[2]
results.b[2]

# ! what about making plots the residuals?
# ! then we need to run the code with the plots turned on
# ! this is turned off by default because (currently) it makes no senso doing th
# ! this by Monte Carlo (below). Plotting should be done just once
# ! we COULD do it more often by e.g. saving residuals for each Monte Carlo iteration
# ! and then averaging (to get more stable results)
results <- analyse.DGP(DGP.1 ,DGP.2 ,DGP.3 ,bc.groups,plot.resid = TRUE)
# ! plots should appear naturally on the right (use arrows to move through them)
# note that:
# black = no control
# blue = standardized outcome (z-scores)
# red = linear controls (continuous variables)
# green = factor controls full (discrete variables)
# purple = factor controls (grouped)
results[2]







####### now we can generate data and analyze it many times ########
####### i.e. do a simulation study #######

# how many simulations?
it.size <- 999
#it.size <- 250

# make an array in which to save the results
results.array <- rep(NA,3*5*it.size)
dim(results.array) <- c(3,5,it.size)
# make an array in which to save the results
results.a.array <- rep(NA,3*5*it.size)
dim(results.a.array) <- c(3,5,it.size)
# make an array in which to save the results
results.b.array <- rep(NA,3*5*it.size)
dim(results.b.array) <- c(3,5,it.size)

clustersizemod <- c(1,1/2,1/3)
t1 <- Sys.time()
for(i in 1:it.size) {
  print(i)
  
  ## generate data
  # DGP1: effect.fambc = 0; effect.bc = 0; effect.bc.sex = 0
  DGP.1 <- DGP(n=20000,sibsize=3,sibspacing=2,
               bcstart=1940,bcend=1980,
               famaycenter=1955,famaysd=2,effect.fambc=0,
               famfemean=1,famfesd=1,
               effect.aar=1,effect.bc=0,
               effect.sex=1,effect.bc.sex=0,
               effect.famfe=1,
               yint=0,ysd=1,
               seed=i*(1234-i))
  
  # DGP2: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = 0
  DGP.2 <- DGP(n=20000,sibsize=3,sibspacing=2,
               bcstart=1940,bcend=1980,
               famaycenter=1955,famaysd=2,effect.fambc=0.5,
               famfemean=1,famfesd=1,
               effect.aar=1,effect.bc=1,
               effect.sex=1,effect.bc.sex=0,
               effect.famfe=1,
               yint=0,ysd=1,
               seed=i*(1234-i))
  
  # DGP3: effect.fambc = nonzero; effect.bc = nonzero; effect.bc.sex = nonzero
  DGP.3 <- DGP(n=20000,sibsize=3,sibspacing=2,
               bcstart=1940,bcend=1980,
               famaycenter=1955,famaysd=2,effect.fambc=0.5,
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
  
  # sample size reduction is done so that all samples have the same size
  # due to MC error reduction this will have no bearing on the findings
  # but co-authors feel that this is necessary
  DGP.1 <- cluster.reduce.sample(DGP.1,10000,"famid",1/3) # everyone has 3 siblings
  DGP.2 <- cluster.reduce.sample(DGP.2,10000,"famid",1/3) # so 10.000 * 3 / 3 = 10.000
  DGP.3 <- cluster.reduce.sample(DGP.3,10000,"famid",1/3) 
  
  DGP.1a <- cluster.reduce.sample(DGP.1a,10000,"famid",clustersizemod)
  DGP.2a <- cluster.reduce.sample(DGP.2a,10000,"famid",clustersizemod)
  DGP.3a <- cluster.reduce.sample(DGP.3a,10000,"famid",clustersizemod)
  
  DGP.1b <- cluster.reduce.sample(DGP.1b,10000,"famid",clustersizemod)
  DGP.2b <- cluster.reduce.sample(DGP.2b,10000,"famid",clustersizemod)
  DGP.3b <- cluster.reduce.sample(DGP.3b,10000,"famid",clustersizemod)
    
  # analyse the DGPs and save results in an array
  results.array[,,i] <-   analyse.DGP(DGP.1 ,DGP.2 ,DGP.3 ,bc.groups)[[1]]
  results.a.array[,,i] <- analyse.DGP(DGP.1a,DGP.2a,DGP.3a,bc.groups)[[1]]
  results.b.array[,,i] <- analyse.DGP(DGP.1b,DGP.2b,DGP.3b,bc.groups)[[1]]

}
t2 <- Sys.time()
t2 - t1


# these are the results by simulation:
results.array
results.a.array
results.b.array

# results averaged over the simulations:
DGP123 <- apply(results.array,c(1,2),mean)
DGP123a <- apply(results.a.array,c(1,2),mean)
DGP123b <- apply(results.b.array,c(1,2),mean)

write.csv(DGP123,"DGP123 999it.csv")
write.csv(DGP123a,"DGP123a 999it.csv")
write.csv(DGP123b,"DGP123b 999it.csv")

DGP123 <- read.csv("DGP123 999it.csv")

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

# make a function just to analyse DGP2, varying the proportion discon
# this is to test the effect of different ratios of discordant on bc but concordant on aar
analyse.DGP.2 <- function(DGP.2,
                        bc.groups) {
  
  # make dataframe in which to put results by DGP (1) and analytical model (5)
  results <- rep(NA,1*5)
  # dim(results) <- c(1,5)
  
  # ! not sure if best solution: some DGPs might produce bcs outside of this range? test
  # ! likely a minor issue though
  DGP.2$bc.grouped <- as.character(cut(DGP.2$bc,breaks=bc.groups))

  # determine quantiles by birth cohort (A2):
  for(bc in unique(DGP.2$bc)) {
    DGP.2$y.q[DGP.2$bc==bc] <- rank(DGP.2$y[DGP.2$bc==bc])/length(DGP.2$y[DGP.2$bc==bc])
  }

  # analyse DGP.2 data (first line is analytical model A1, second is A2, etc.)
  results[1] <- plm(y ~ aar, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)$coef[1]
  results[2] <- plm(y.q ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)$coef[1]
  results[3] <- plm(y ~ aar + bc + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)$coef[1]
  results[4] <- plm(y ~ aar + as.factor(bc) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)$coef[1]
  results[5] <- plm(y ~ aar + as.factor(bc.grouped) + sex, effect='individual',model='within',index='famid',family='gaussian',data=DGP.2)$coef[1]

  return(results)
}

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

apply(results.mat,c(2,3),mean)

conv.mean <- function(x) {cumsum(x)/(1:length(x))}

# is analysis method 5 highly unstable?
mean(results.array[,5,1])
plot(conv.mean(results.array[,5,1]))
plot(conv.mean(results.array[,5,5]))
# the estimate is unstable, so it is the method that is faulty!
