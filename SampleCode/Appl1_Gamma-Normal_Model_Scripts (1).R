graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
library(rjags)
library(runjags)
library(coda)
setwd("~/Documents/MATH2269_Bayesian/2024/presentations/Module 5")
source("DBDA2E-utilities.R")   
path ="/Users/haydardemirhan/Documents/MATH2269_Bayesian/2024/presentations/Module 5/"

myData = read.csv("rides.csv")
y = myData$y
type = as.numeric(as.factor(myData$Type)) 
# converts character to consecutive integer levels

# Load data

Ntotal = length(y)
Nt = length(unique(type))
dataList = list(
  y = y ,
  type = type ,
  Ntotal = Ntotal ,
  Nt = Nt
)

# Specify model
modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dnorm( mu[type[i]], 1/sigma.sq[type[i]] )
    # notice nested indexing and the variance of normal distribution
  }
  mu[1] ~ dnorm(7500 , 1/5) # MTB is coded 1
  mu[2] ~ dnorm(6000 , 1/5) # MTB is coded 2
  sigma.sq[1] ~ dgamma(25, 22.5)
  sigma.sq[2] ~ dgamma(25, 22.5)
}
" # close quote for modelString
writeLines( modelString , con="modelText.txt" )

# Initialize chains
initsList = list(
  mu = c(-10,0),
  sigma.sq = c(-10000,1)
)

# Generate chains
nChains = 3
nThinSteps = 2
nUseSteps = 1000
startTime = proc.time()
jagsModel = jags.model(
  file="modelText.txt" ,
  data = dataList ,
  # inits=initsList , #Let JAGS generate the inits
  n.chains = nChains ,
  n.adapt = 500 
)
update(jagsModel , n.iter=1000 )
codaSamples = coda.samples( jagsModel ,
                            variable.names=c("mu[1]","mu[2]","sigma.sq[1]","sigma.sq[2]") ,
                            n.iter=1000,#ceiling(nUseSteps*nThinSteps/nChains) ,
                            #initsList = initsList,
                            thin=nThinSteps )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)

# Alternative parallel run
library(runjags)
startTime = proc.time()
runJagsOut <- run.jags( method="parallel" ,
                        model="modelText.txt" ,
                        monitor=c("mu[1]","mu[2]","sigma.sq[1]",
                                  "sigma.sq[2]")  ,
                        data=dataList ,
                        # inits=initsList ,
                        n.chains=nChains ,
                        adapt=500 ,
                        burnin=1000 ,
                        sample=1000,#ceiling(nUseSteps/nChains) ,
                        thin=nThinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )
stopTime = proc.time()

elapsedTime = stopTime - startTime
show(elapsedTime)
# Alternative parallel run

# Examine chains
diagMCMC( codaSamples , parName="mu[1]" )
diagMCMC( codaSamples , parName="mu[2]" )
diagMCMC( codaSamples , parName="sigma.sq[1]" )
diagMCMC( codaSamples , parName="sigma.sq[2]" )
graphics.off()
          
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
                                        compVal=compVal , ROPE=rope ) )
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
  mu1 = mcmcMat[,"mu[1]"]
  mu2 = mcmcMat[,"mu[2]"]
  sigma.sq1 = mcmcMat[,"sigma.sq[1]"]
  sigma.sq2 = mcmcMat[,"sigma.sq[2]"]
  
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( mu1 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(mu[1]) , main="mu1" )
  histInfo = plotPost( mu2 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(mu[2]) , main="mu2" )
  histInfo = plotPost( sigma.sq1 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(sigma[1]) , main="sigma.sq1" )
  histInfo = plotPost( sigma.sq2 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(sigma[2]) , main="sigma.sq2" )
  
}

# Then use the functions to get the results to talk on
smryMCMC( codaSamples , compVal=6500  )

plotMCMC( codaSamples , data=myData , compVal=NULL )
