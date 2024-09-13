graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
library(rjags)
library(runjags)
library(coda)
setwd("~/Documents/MATH2269_Bayesian/2024/presentations/Module 5")
source("DBDA2E-utilities.R")   

# data <- read.csv("~/Documents/MATH2269_Bayesian/2020/presentations/Module 4/COVID19Cases.csv") # Analysed in Week 4
data <- read.csv("COVID19Cases20Aug.csv") # Full dataset
# data <- read.csv("COVID19Cases20Aug2ndWave.csv") # Only the 2nd wave
data <- data$y
summary(data)
min(data) # There are 0 values in data. These will create problem with Gamma prior.
data[which(data == 0)] = 0.0001
hist(data, main = "Histogram of COVID case numbers")

# Load data
Ntotal = length(data)
dataList = list(
  y = data ,
  Ntotal = Ntotal 
)

# Specify model
modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dgamma( mu^2/tau, mu/tau )
  }
  mu ~ dnorm(250 , 1/25) # Change the variance for the degree of informativeness
  tau ~ dgamma(123.45, 0.2469)
}
"
writeLines( modelString , con="TEMPmodel2.txt" )

# Initialize chains

initsList = list(
  mu = 1, 
  tau = 10
)

# Generate chains
parameters = c( "mu" ,  "tau")

adaptSteps = 500
burnInSteps = 1000
nChains = 3
thinSteps = 11
numSavedSteps = 2000

runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel2.txt" ,
                        monitor=parameters  ,
                        data=dataList ,
                        inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )

# Examine chains
diagMCMC( codaSamples , parName="mu" )
diagMCMC( codaSamples , parName="tau" )

# Interpret the posteriors
# First prepare the functions to get summary statistics and plots of posterior distributions.
smryMCMC = function(  codaSamples, compVal=NULL , rope=NULL ) {
# This function is adapted from the scripts given by 'J. Kruschke, Doing Bayesian Data Analysis, 2014, Elsevier.'
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,pName] , 
                                        compVal=unlist(compVal[pName]) , ROPE=rope ) )
  }
  rownames(summaryInfo) = paramName
  return( summaryInfo)
}

decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                            nRow=2 , nCol=3 ) {
# This function is adapted from the scripts given by 'J. Kruschke, Doing Bayesian Data Analysis, 2014, Elsevier.'
  # If finishing a set:
  if ( finished==TRUE ) {
    if ( !is.null(saveName) ) {
      saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                 type=saveType)
    }
    panelCount = 1 # re-set panelCount
    return(panelCount)
  } else {
    # If this is first panel of a graph:
    if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
      # If previous graph was open, save previous one:
      if ( panelCount>1 & !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                   type=saveType)
      }
      # Open new graph
      openGraph(width=nCol*7.0/3,height=nRow*2.0)
      layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
      par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
    }
    # Increment and return panel count:
    panelCount = panelCount+1
    return(panelCount)
  }
}


plotMCMC = function( codaSamples , data , compVal=NULL , rope=NULL , 
                     showCurve=FALSE , saveName = NULL) {
# This function is adapted from the scripts given by 'J. Kruschke, Doing Bayesian Data Analysis, 2014, Elsevier.'
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  mu = mcmcMat[,"mu"]
  tau = mcmcMat[,"tau"]

  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( mu , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(mu) , main="mu" ,compVal = 50 )
  histInfo = plotPost( tau , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(tau) , main="tau" )
  
}

# Then use the functions to get the results to talk on
smryMCMC( codaSamples , compVal=list(mu =3, tau = 49)  )

plotMCMC( codaSamples , data=myData )



