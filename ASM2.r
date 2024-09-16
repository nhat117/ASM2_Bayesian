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
# Prepare the data for JAGS
# Assuming you have the property_data dataframe with 'SalePrice' as the target variable:
# Assuming 'SalePrice' is the name of the column with sale prices
# Load your dataset

property_data <- read.csv("./Assignment2PropertyPrices.csv")
property_data <- sample_n(property_data, 500)
y = property_data$SalePrice * 100000  # Sale price in AUD
x =  as.matrix(property_data[, c("Area", "Bedrooms", "Bathrooms", "CarParks", "PropertyType")])
# Prepare data for JAGS
dataList <- list(
  x = x,
  y = y,
  xPred = xPred,
  Nx = dim(x)[2],
  Ntotal = dim(x)[1],
  Npred = dim(xPred)[1]
)

# JAGS model string with priors scaled for 100k AUD units
model_string <- "
  # Standardize the data:
  data {
    # Define Mean and Variance of y
    ym <- mean(y)
    ysd <- sd(y)
    
    # Standardize Y
    for (i in 1:Ntotal) {
      zy[i] <- (y[i] - ym) / ysd
    }

    for (j in 1:Nx) {
      xm[j] <- mean(x[,j])
      xsd[j] <- sd(x[,j])
      for (i in 1:Ntotal) {
        zx[i,j] <- (x[i,j] - xm[j]) / xsd[j]
      }
    }

    # Prior locations to reflect the expert information
    mu[1] <- 90        # Area
    mu[2] <- 100000    # Bedrooms
    mu[3] <- 100000    # Bathrooms
    mu[4] <- 120000    # CarParks
    mu[5] <- -150000   # PropertyType
    
    # Prior variances to reflect the expert information
    Var[1] <- 10^2     # Area
    Var[2] <- 10^5     # Bedrooms
    Var[3] <- 10^7     # Bathrooms
    Var[4] <- 999      # CarParks
    Var[5] <- 10       # PropertyType
    
    # Compute corresponding prior means and variances for the standardized parameters
    muZ[1:Nx] <- mu[1:Nx] * xsd[1:Nx] / ysd
  }

  # Model block for standardized data
  model {
    for (i in 1:Ntotal) {
      zy[i] ~ dt(zbeta0 + sum(zbeta[1:Nx] * zx[i,1:Nx]), 1/zsigma^2, nu)
    }

    # Priors on standardized scale
    zbeta0 ~ dnorm(0, 1/10^6)  

    for (j in 1:Nx) {
      zbeta[j] ~ dnorm(muZ[j], 1/Var[j])
    }

    zsigma ~ dunif(1.0E-5, 1.0E+5)
    nu ~ dexp(1/30.0)

    # Transform to original scale
    beta[1:Nx] <- (zbeta[1:Nx] / xsd[1:Nx]) * ysd
    beta0 <- zbeta0 * ysd + ym - sum(zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx]) * ysd
    sigma <- zsigma * ysd
    
    # Compute predictions at every step of the MCMC
    for ( i in 1:Npred){
      pred[i] <- beta0 + beta[1] * xPred[i,1] + beta[2] * xPred[i,2] + beta[3] * xPred[i,3] + beta[4] * xPred[i,4] + beta[5] * xPred[i,5]
    }
  }
"
# Write the model to a file
writeLines(model_string, con="TEMPmodel.txt")

# Define the parameters to monitor
parameters <- c("beta0" ,  "beta" , "sigma", "zbeta0" ,"zbeta", "zsigma", "nu","pred")

# MCMC settings
adaptSteps = 500       # Number of adaptation steps
burnInSteps = 1000     # Number of burn-in steps
nChains = 2            # Number of chains
thinSteps = 4          # Thinning parameter
numSavedSteps = 10000  # Number of saved MCMC steps
nIter = ceiling((numSavedSteps * thinSteps) / nChains)  # Total number of iterations per chain run

# Measure the runtime of the MCMC process using system.time()
runtime <- system.time({
  runJagsOut <- run.jags(method="parallel",
                         model = "TEMPmodel.txt", 
                         data = dataList, 
                         monitor = parameters, 
                         n.chains=nChains ,
                         adapt=adaptSteps ,
                         burnin=burnInSteps ,
                         sample=numSavedSteps ,
                         thin=thinSteps , summarise=FALSE , plots=FALSE )
})
# Print the runtime
print(paste("Runtime:", runtime["elapsed"], "seconds"))

codaSamples = as.mcmc.list( runJagsOut)
# Do further thinning from the codaSamples
furtherThin <- 7
thiningSequence <- seq(1,nrow(codaSamples[[1]]), furtherThin)
newCodaSamples <- mcmc.list()
for ( i in 1:nChains){
  newCodaSamples[[i]] <- as.mcmc(codaSamples[[i]][thiningSequence,])
}
#================PREDICTIONS===================
# Extract and summarize predictions for the new properties
summary(newCodaSamples)
predictions <- as.matrix(newCodaSamples)[,grep("pred", colnames(as.matrix(newCodaSamples)))]

# Present the posterior means for the predicted sale prices
apply(predictions, 2, mean)  # These are the predicted sale prices

#=============== Plot MCMC HD =================
plotMCMC_HD = function(codaSamples , data , xName="x" , yName="y" ,
                        showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL,
                        saveName=NULL , saveType="jpg" ) {
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples, chains=TRUE)
  chainLength = NROW(mcmcMat)
  
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[", colnames(mcmcMat))]
  if (ncol(x)==1) { zbeta = matrix(zbeta, ncol=1) }
  
  zVar = mcmcMat[,"nu"]
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[", colnames(mcmcMat))]
  if (ncol(x)==1) { beta = matrix(beta, ncol=1) }
  
  sigma = mcmcMat[,"sigma"]
  pred1 = mcmcMat[,"pred[1]"]
  pred2 = mcmcMat[,"pred[2]"]
  pred3 = mcmcMat[,"pred[3]"]
  pred4 = mcmcMat[,"pred[4]"]
  pred5 = mcmcMat[,"pred[5]"]
  
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
  histInfo = plotPost(sigma, cex.lab = 1.75, showCurve=showCurve, xlab=bquote(sigma), main="Precision")
  
  # Plot predictions:
  for (i in 1:5) {
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

# Define the list of parameters to check
param_names <- c("beta0","beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","sigma" ,"zsigma", "nu")

# Loop through the parameter names and call diagMCMC for each one
for (parName in param_names) {
  diagMCMC(newCodaSamples, parName = parName)
}
# Call the plotMCMC_HD function with the coda samples
plotMCMC_HD(codaSamples = newCodaSamples, data = property_data, xName=c("Area","Bedrooms","Bathrooms","CarParks","PropertyType"), yName="SalePrice")


