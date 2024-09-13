graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
setwd("~/Documents/MATH2269_Bayesian/2024/presentations/Module 6")
source("DBDA2E-utilities.R")      

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

smryMCMC = function(  codaSamples , 
                      compValBeta0=NULL , ropeBeta0=NULL , 
                      compValBeta1=NULL , ropeBeta1=NULL , 
                      compValVar=NULL , ropeVar=NULL , 
                      saveName=NULL ) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = rbind( summaryInfo , 
                       "beta0" = summarizePost( mcmcMat[,"beta0"] , 
                                                compVal=compValBeta0 , 
                                                ROPE=ropeBeta0 ) )
  summaryInfo = rbind( summaryInfo , 
                       "beta1" = summarizePost( mcmcMat[,"beta1"] , 
                                                compVal=compValBeta1 , 
                                                ROPE=ropeBeta1 ) )
  summaryInfo = rbind( summaryInfo , 
                       "Var" = summarizePost( mcmcMat[,"Var"] , 
                                              compVal=compValVar , 
                                              ROPE=ropeVar ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , data , xName="x" , yName="y" ,
                     compValBeta0=NULL , ropeBeta0=NULL , 
                     compValBeta1=NULL , ropeBeta1=NULL , 
                     compValVar=NULL , ropeVar=NULL , 
                     showCurve=FALSE ,  pairsPlot=FALSE ,
                     saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = data[,xName]
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta1 = mcmcMat[,"zbeta1"]
  zVar = mcmcMat[,"zVar"]
  beta0 = mcmcMat[,"beta0"]
  beta1 = mcmcMat[,"beta1"]
  Var = mcmcMat[,"Var"]
  
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta1 , Var )[plotIdx,] ,
           labels=c( expression(beta[0]) , expression(beta[1]) , 
                     expression(Var) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  # Set up window and layout:
  nPtToPlot = 1000
  plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
  openGraph(width=8,height=5)
  layout( matrix( 1:6 , nrow=2, byrow=TRUE ) )
  par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta0 , ROPE=ropeBeta0 ,
                       xlab=bquote(beta[0]) , main=paste("Intercept") )
  histInfo = plotPost( beta1 , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValBeta1 , ROPE=ropeBeta1 ,
                       xlab=bquote(beta[1]) , main=paste("Slope") )
  plot( beta1[plotIdx] , beta0[plotIdx] , 
        xlab=bquote(beta[1]) , ylab=bquote(beta[0]) ,
        col="skyblue" , cex.lab = 1.75 )
  histInfo = plotPost( Var , cex.lab = 1.75 , showCurve=showCurve ,
                       compVal=compValVar , ROPE=ropeVar ,
                       xlab=bquote(Var) , main=paste("Scale") )
  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"PostMarg",sep=""), type=saveType)
  }
  #-----------------------------------------------------------------------------
}

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

myData <- read.csv("activities.csv")
myData$Type <- as.factor(myData$Type)
head(myData)

# Scatter plots
p1 <- ggplot(myData, aes(x=Distance.km., y=Elapsed.Time.sec.,  color=Type)) +
  geom_point()

p2 <- ggplot(myData, aes(x=ElevationGain.m., y=Elapsed.Time.sec.,  color=Type)) +
  geom_point()

p3 <- ggplot(myData, aes(x=AverageGrade.percent., y=Elapsed.Time.sec.,  color=Type)) +
  geom_point()

p4 <- ggplot(myData, aes(x=AverageWatts, y=Elapsed.Time.sec.,  color=Type)) +
  geom_point()
figure <- ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
figure

ggplot(myData[which(myData$Type == "Road"),], aes(x=Distance.km., y=Elapsed.Time.sec.)) +
  geom_point()

# Kernel density estimation
plot(kde(myData$Elapsed.Time.sec.), xlab = "Elapsed.Time.sec.") # with default settings


# THE DATA.
y = myData[,"Elapsed.Time.sec."]
x = myData[,"Distance.km."]

# Specify the data in a list, for later shipment to JAGS:
dataList <- list(
  x = x ,
  y = y ,
  Ntotal = length(y)
  
)

# First run without initials!
initsList <- list(
  beta0 = 2000,
  beta1 = 250,
  Var = 12000000
)

# WE WILL RUN THE MODEL WITHOUT SCALING FIRST and THEN WITH SCALING!

# THE MODEL without scaling
modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dgamma( (mu[i]^2)/Var , mu[i]/Var ) 
    mu[i] <- beta0 + beta1 * x[i]
  }
  # Priors vague on standardized scale:
  beta0 ~ dnorm( 0 , 1/1000 ) # When the prior variances are reduced, the impact of an informative prior can be seen!
  beta1 ~ dnorm( 0 , 1/1000 )
  Var ~ dgamma( 0.1 , 0.1 )
}
" # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )

parameters = c( "beta0" ,  "beta1" ,  "Var")
# First run parameters:
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 4 
thinSteps = 7
numSavedSteps = 1000

# Second run parameters:
adaptSteps = 1500  # A smaller value shows the impact of adaptSteps
burnInSteps = 5000 # A smaller value shows the impact of burnInSteps
nChains = 3
thinSteps = 11
numSavedSteps = 1000

# Slow run
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
# Create, initialize, and adapt the model:
# First run: do not provide inits!
startTime = proc.time()
jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nIter , thin=thinSteps )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)

# Parallel run
startTime = proc.time()
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=c("beta0" ,  "beta1" ,  "Var")  ,
                        data=dataList ,
                        inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)

diagMCMC( codaSamples , parName="beta0" )
diagMCMC( codaSamples , parName="beta1" )
diagMCMC( codaSamples , parName="Var" )
graphics.off()

# THE MODEL with scaling

dataList <- list(
  x = x ,
  y = y 
)
# 
# xsd <- sd(x)
# ysd <- sd(y)
# zx <- array(NA, Ntotal)
# zy <- array(NA, Ntotal)
# for ( i in 1:length(y) ) {
#   zx[i] <- x[i] / xsd
#   zy[i] <- y[i] / ysd
# }
# 
# summary(zy)
# summary(zx)
# var(zy)/var(zx)

initsList <- list(
  zbeta0 = 0.1,
  zbeta1 = 0.5,
  zVar = 1
)


modelString <- "
  # Scale the data:
  data {
    Ntotal <- length(y)
    xsd <- sd(x)
    ysd <- sd(y)
    for ( i in 1:length(y) ) {
      zx[i] <- x[i] / xsd
      zy[i] <- y[i] / ysd
    }
  }
  model {
    for ( i in 1:Ntotal ) {
      zy[i] ~ dgamma( (mu[i]^2)/zVar , mu[i]/zVar ) 
      mu[i] <- zbeta0 + zbeta1 * zx[i] 
    }
    # Priors vague on scaled scale:
    zbeta0 ~ dnorm( 0 , 1/(10)^2 ) # 1/ variance for normal distribution
    zbeta1 ~ dnorm( 0 , 1/(10)^2 ) # 1/ variance for normal distribution
    zVar ~ dgamma( 1 , 1/(10)^2 ) # 1/ variance for normal distribution
    
    # Transform to original scale:
    beta1 <- zbeta1 * ysd / xsd
    beta0 <- zbeta0 * ysd  
    Var <- zVar * ysd^2
}"
writeLines( modelString , con="TEMPmodel.txt" )

parameters = c( "zbeta0" ,  "zbeta1" , "beta0" ,  "beta1" ,  "Var", "zVar")
# First run parameters:
adaptSteps = 1500  # Number of steps to "tune" the samplers
burnInSteps = 5000
nChains = 2 
thinSteps = 11
numSavedSteps = 1500
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

# First run: do not provide inits!
startTime = proc.time()
jagsModel = jags.model( "TEMPmodel.txt" , data=dataList ,  inits=initsList ,
                        n.chains=nChains , n.adapt=adaptSteps )
# Burn-in:
cat( "Burning in the MCMC chain...\n" )
update( jagsModel , n.iter=burnInSteps )
# The saved MCMC chain:
cat( "Sampling final MCMC chain...\n" )
codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nIter , thin=thinSteps )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)

# Parallel run
startTime = proc.time()
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=c( "zbeta0" ,  "zbeta1" , "beta0" ,  "beta1" ,  "Var", "zVar")  ,
                        data=dataList ,
                        inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)


diagMCMC( codaSamples , parName="beta0" )
diagMCMC( codaSamples , parName="beta1" )
diagMCMC( codaSamples , parName="Var" )
diagMCMC( codaSamples , parName="zbeta0" )
diagMCMC( codaSamples , parName="zbeta1" )
graphics.off()

smryMCMC( codaSamples , compValBeta0=0 , compValBeta1=0.5  )

plotMCMC( codaSamples , data = myData, xName="Distance.km." , yName="Elapsed.Time.sec.", 
          compValBeta1=250)

# Predictive check

y2 <- 3720 +192 *x

plot(y,y2)

MAE <- mean(abs(y2-y)); MAE
RMSE <- sqrt(mean((y2-y)^2)); RMSE

