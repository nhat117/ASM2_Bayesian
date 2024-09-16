# Clear R's environment and close all graphics
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

# Load required libraries
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
# Load the necessary library
library(dplyr)


# Set the working directory#setwd("~/Documents/MATH2269_Bayesian/2024/presentations/Module 6/Application2")
source("DBDA2E-utilities.R")
property_data <- read.csv("./Assignment2PropertyPrices.csv")
set.seed(42)  # Setting a seed for reproducibility
property_data <- sample_n(property_data, 1000)
#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================
smryMCMC_HD = function(  codaSamples , compVal = NULL,  saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] , 
                                                          compVal = as.numeric(compVal[pName]) ))
      }
      else {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
      }
    } else {
      summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
    }
  }
  rownames(summaryInfo) = paramName
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================
# Include new prediction data

# Define the new prediction data (from the image you provided):
xPred <- array(NA, dim = c(5, 5))  # 5 properties with 5 predictors

# Populate the data (Area, Bedrooms, Bathrooms, CarParks, PropertyType)
xPred[1,] <- c(600, 2, 2, 1, 1)  # Property 1 (Unit)
xPred[2,] <- c(800, 3, 1, 2, 0)  # Property 2 (House)
xPred[3,] <- c(1500, 2, 1, 1, 0)  # Property 3 (House)
xPred[4,] <- c(2500, 5, 4, 4, 0)  # Property 4 (House)
xPred[5,] <- c(250, 3, 2, 1, 1)   # Property 5 (Unit)

# Prepare the data for JAGS
dataList <- list(
  x = cbind(property_data$Area, property_data$Bedrooms, property_data$Bathrooms, property_data$CarParks, property_data$PropertyType),
  y = property_data$SalePrice,  # SalePrice in 100k AUD
  Nx = 5,  # Number of predictors
  Ntotal = nrow(property_data),  # Number of observations
  xPred = xPred  # Prediction data for the 5 new properties
)

# Initial values based on expert knowledge
initsList <- list(
  zbeta0 = 2000,  # Intercept
  zbeta = c(0.9, 1, 0, 1.2, -1.5),  # Expert knowledge for Area, Bedrooms, Bathrooms, CarParks, PropertyType
  Var = 12000000  # Variance
)

# JAGS model string with priors scaled for 100k AUD units
modelString = "
data {
  ysd <- sd(y)
  for (i in 1:Ntotal) {
    zy[i] <- y[i] / ysd
  }
  for (j in 1:Nx) {
    xsd[j] <- sd(x[,j])
    for (i in 1:Ntotal) {
      zx[i,j] <- x[i,j] / xsd[j]
    }
  }
}
model {
  for (i in 1:Ntotal) {
    zy[i] ~ dgamma((mu[i]^2) / zVar, mu[i] / zVar)
    mu[i] <- zbeta0 + sum(zbeta[1:Nx] * zx[i,1:Nx])
  }
  
  # Priors for intercept and predictors, scaled to 100k AUD units
  zbeta0 ~ dnorm(0, 1/2^2)
  zbeta[1] ~ dnorm(0.9 / xsd[1], 1/(4/xsd[1]^2))  # Area: 90,000 AUD -> 0.9 in 100k units
  zbeta[2] ~ dnorm(1 / xsd[2], 1/(4/xsd[2]^2))    # Bedrooms: 100,000 AUD -> 1 in 100k units
  zbeta[3] ~ dnorm(0, 1/4)                        # Bathrooms: no expert knowledge
  zbeta[4] ~ dnorm(1.2 / xsd[4], 1/(4/xsd[4]^2))  # CarParks: 120,000 AUD -> 1.2 in 100k units
  zbeta[5] ~ dnorm(-1.5 / xsd[5], 1/(4/xsd[5]^2)) # PropertyType: -150,000 AUD -> -1.5 in 100k units
  
  zVar ~ dgamma(0.01, 0.01)
  
  # Back-transform to the original scale
  beta[1:Nx] <- (zbeta[1:Nx] / xsd[1:Nx]) * ysd
  beta0 <- zbeta0 * ysd
  tau <- zVar * (ysd)^2
  
  # Predictions for new properties
  for (i in 1:5) {  # 5 new properties based on the table you provided
    pred[i] <- beta0 + beta[1] * xPred[i,1] + beta[2] * xPred[i,2] + beta[3] * xPred[i,3] + beta[4] * xPred[i,4] + beta[5] * xPred[i,5]
  }
}
"

# Write the model to a file
writeLines(modelString, con="TEMPmodel.txt")

# Define the parameters to monitor
parameters <- c("zbeta0", "zbeta", "beta0", "beta", "tau", "zVar", "pred")

# MCMC settings
adaptSteps = 500       # Number of adaptation steps
burnInSteps = 1000     # Number of burn-in steps
nChains = 2            # Number of chains
thinSteps = 3          # Thinning parameter
numSavedSteps = 10000  # Number of saved MCMC steps
nIter = ceiling((numSavedSteps * thinSteps) / nChains)  # Total number of iterations per chain

# Run JAGS using the run.jags function
runJagsOut <- run.jags(method = "parallel",
                       model = "TEMPmodel.txt",
                       monitor = parameters,
                       data = dataList,
                       inits = initsList,
                       n.chains = nChains,
                       adapt = adaptSteps,
                       burnin = burnInSteps,
                       sample = numSavedSteps,
                       thin = thinSteps,
                       summarise = FALSE,
                       plots = FALSE)

# Convert the output to coda samples for further analysis
codaSamples = as.mcmc.list(runJagsOut)

#================PREDICTIONS===================
# Extract and summarize predictions for the new properties
summary(codaSamples)
predictions <- as.matrix(codaSamples)[,grep("pred", colnames(as.matrix(codaSamples)))]

# Present the posterior means for the predicted sale prices
apply(predictions, 2, mean)  # These are the predicted sale prices

#=============== Plot MCMC HD =================
plotMCMC_HD = function( codaSamples , data , xName="x" , yName="y" ,
                        showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL,
                        saveName=NULL , saveType="jpg" ) {
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples, chains=TRUE)
  chainLength = NROW(mcmcMat)
  
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[", colnames(mcmcMat))]
  if (ncol(x)==1) { zbeta = matrix(zbeta, ncol=1) }
  
  zVar = mcmcMat[,"zVar"]
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[", colnames(mcmcMat))]
  if (ncol(x)==1) { beta = matrix(beta, ncol=1) }
  
  tau = mcmcMat[,"tau"]
  pred1 = mcmcMat[,"pred[1]"]
  pred2 = mcmcMat[,"pred[2]"]
  
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
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
  # Plot marginal histograms:
  panelCount = 1
  panelCount = decideOpenGraph(panelCount, saveName=paste0(saveName,"PostMarg"))
  histInfo = plotPost(beta0, cex.lab = 1.75, showCurve=showCurve, xlab=bquote(beta[0]), main="Intercept")
  
  for (bIdx in 1:ncol(beta)) {
    panelCount = decideOpenGraph(panelCount, saveName=paste0(saveName,"PostMarg"))
    histInfo = plotPost(beta[,bIdx], cex.lab = 1.75, showCurve=showCurve, xlab=bquote(beta[.(bIdx)]), main=xName[bIdx])
  }
  
  panelCount = decideOpenGraph(panelCount, saveName=paste0(saveName,"PostMarg"))
  histInfo = plotPost(tau, cex.lab = 1.75, showCurve=showCurve, xlab=bquote(tau), main="Scale")
  
  # Plot predictions:
  for (i in 1:2) {
    panelCount = decideOpenGraph(panelCount, saveName=paste0(saveName,"PostMarg"))
    histInfo = plotPost(mcmcMat[,paste0("pred[", i, "]")], cex.lab = 1.75, showCurve=showCurve, xlab=paste("pred", i), main=paste("Prediction", i))
  }
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
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
}}


diagMCMC( codaSamples , parName="beta0" )
diagMCMC( codaSamples , parName="beta[1]" )
diagMCMC( codaSamples , parName="beta[2]" )
diagMCMC( codaSamples , parName="beta[3]" )
diagMCMC( codaSamples , parName="beta[4]" )
diagMCMC( codaSamples , parName="tau" )
diagMCMC( codaSamples , parName="pred[1]" )
diagMCMC( codaSamples , parName="pred[2]" )
diagMCMC( codaSamples , parName="zbeta0" )
diagMCMC( codaSamples , parName="zbeta[1]" )
diagMCMC( codaSamples , parName="zbeta[2]" )
diagMCMC( codaSamples , parName="zbeta[3]" )
diagMCMC( codaSamples , parName="zbeta[4]" )

# Call the plotMCMC_HD function with the coda samples
plotMCMC_HD(codaSamples = codaSamples, data = property_data, xName=c("Area","Bedrooms","Bathrooms","CarParks","PropertyType"), yName="SalePrice")


