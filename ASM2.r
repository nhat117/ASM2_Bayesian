# Load necessary libraries
graphics.off()  # Close all open graphics windows
rm(list=ls())   # Clear all objects from memory
library(ggplot2)
library(rjags)
library(runjags)
library(coda)
library(ggpubr)
source("DBDA2E-utilities.R")


# Define the smryMCMC_HD function
smryMCMC_HD = function(codaSamples, compVal = NULL, saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for (pName in paramName) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind(summaryInfo, summarizePost(paramSampleVec = mcmcMat[,pName], compVal = as.numeric(compVal[pName])))
      } else {
        summaryInfo = rbind(summaryInfo, summarizePost(paramSampleVec = mcmcMat[,pName]))
      }
    } else {
      summaryInfo = rbind(summaryInfo, summarizePost(paramSampleVec = mcmcMat[,pName]))
    }
  }
  rownames(summaryInfo) = paramName
  if (!is.null(saveName)) {
    write.csv(summaryInfo, file=paste0(saveName, "SummaryInfo.csv"))
  }
  return(summaryInfo)
}

# Define the plotMCMC_HD function
plotMCMC_HD = function(codaSamples, data, xName="x", yName="y", showCurve=FALSE, pairsPlot=FALSE, compVal = NULL, saveName=NULL, saveType="jpg") {
  y = data[, yName]
  x = as.matrix(data[, xName])
  mcmcMat = as.matrix(codaSamples, chains=TRUE)
  chainLength = NROW(mcmcMat)
  beta0 = mcmcMat[,"beta0"]
  beta = mcmcMat[,grep("^beta$|^beta\\[", colnames(mcmcMat))]
  if (ncol(x) == 1) { beta = matrix(beta, ncol=1) }
  tau = mcmcMat[,"tau"]
  
  # Plot the parameters
  panelCount = 1
  decideOpenGraph = function(panelCount, saveName, finished=FALSE, nRow=2, nCol=3) {
    if (finished == TRUE) {
      if (!is.null(saveName)) {
        saveGraph(file=paste0(saveName, ceiling((panelCount-1)/(nRow*nCol))), type=saveType)
      }
      panelCount = 1
      return(panelCount)
    } else {
      if ((panelCount %% (nRow*nCol)) == 1) {
        if (panelCount > 1 & !is.null(saveName)) {
          saveGraph(file=paste0(saveName, (panelCount %/% (nRow * nCol))), type=saveType)
        }
        openGraph(width=nCol*7.0/3, height=nRow*2.0)
        layout(matrix(1:(nRow*nCol), nrow=nRow, byrow=TRUE))
        par(mar=c(4,4,2.5,0.5), mgp=c(2.5,0.7,0))
      }
      panelCount = panelCount + 1
      return(panelCount)
    }
  }
  
  # Plot for beta0
  panelCount = decideOpenGraph(panelCount, saveName=paste0(saveName, "PostMarg"))
  histInfo = plotPost(beta0, cex.lab = 1.75, showCurve=showCurve, xlab="beta[0]", main="Intercept")
  
  # Plot for each beta
  for (bIdx in 1:ncol(beta)) {
    panelCount = decideOpenGraph(panelCount, saveName=paste0(saveName, "PostMarg"))
    histInfo = plotPost(beta[,bIdx], cex.lab = 1.75, showCurve=showCurve, xlab=paste0("beta[", bIdx, "]"), main=xName[bIdx])
  }
  
  # Plot for tau
  panelCount = decideOpenGraph(panelCount, saveName=paste0(saveName, "PostMarg"))
  histInfo = plotPost(tau, cex.lab = 1.75, showCurve=showCurve, xlab="tau", main="Scale")
  
  panelCount = decideOpenGraph(panelCount, finished=TRUE, saveName=paste0(saveName, "PostMarg"))
}

# Load the dataset
# Load the CSV file
data <- read.csv("Assignment2PropertyPrices.csv")

# Randomly select 500 rows
#data <- data[sample(1:nrow(data), 500, replace = FALSE), ]
# Prepare standardized data for JAGS
x <- as.matrix(data[, c("Area", "Bedrooms", "Bathrooms", "CarParks", "PropertyType")])
y <- data$SalePrice

dataList <- list(
  SalePrice = y,
  x = x,
  Nx = ncol(x),
  Ntotal = nrow(data)
)

# THE MODEL WITH STANDARDIZATION
modelString = "
  data {
    ysd <- sd(SalePrice)
    for (i in 1:Ntotal) {
      zy[i] <- SalePrice[i] / ysd
    }

    for (j in 1:Nx) {
      xsd[j] <- sd(x[, j])
      for (i in 1:Ntotal) {
        zx[i, j] <- x[i, j] / xsd[j]
      }
    }
  }

  model {
    for (i in 1:Ntotal) {
      zy[i] ~ dnorm(mu[i], tau)
      mu[i] <- zbeta0 + sum(zbeta[1:Nx] * zx[i, 1:Nx])
    }

    zbeta0 ~ dnorm(0, 1/10^2)
    zbeta[1] ~ dnorm(90/xsd[1], 1/10^2)
    zbeta[2] ~ dnorm(100000/xsd[2], 1/50000^2)
    zbeta[3] ~ dnorm(0, 1/10^6)
    zbeta[4] ~ dnorm(120000/xsd[4], 1/10000^2)
    zbeta[5] ~ dnorm(-150000/xsd[5], 1/5000^2)

    zVar ~ dgamma(0.01, 0.01)
    tau <- 1 / zVar

    beta0 <- zbeta0 * ysd
    for (j in 1:Nx) {
      beta[j] <- (zbeta[j] * ysd) / xsd[j]
    }
  }
"

model_string = modelString = "
  data {
    # Standardizing SalePrice
    ysd <- sd(SalePrice)
    for (i in 1:Ntotal) {
      zy[i] <- SalePrice[i] / ysd
    }

    # Standardizing predictors (Area, Bedrooms, etc.)
    for (j in 1:Nx) {
      xsd[j] <- sd(x[, j])
      for (i in 1:Ntotal) {
        zx[i, j] <- x[i, j] / xsd[j]
      }
    }
  }

  model {
    # Likelihood for standardized data
    for (i in 1:Ntotal) {
      zy[i] ~ dnorm(mu[i], tau)
      mu[i] <- zbeta0 + sum(zbeta[1:Nx] * zx[i, 1:Nx])
    }

    # Priors for coefficients
    zbeta0 ~ dnorm(0, 1/10^2)  # Weak prior for intercept

    # Area: Strong belief, centered at 90 AUD/m²
    zbeta[1] ~ dnorm(90 / xsd[1], 1 / 10^2)  # Small variance for strong belief

    # Bedrooms: Weak belief, centered at 100,000 AUD/bedroom
    zbeta[2] ~ dnorm(100000 / xsd[2], 1 / 50000^2)  # Larger variance for weak belief

    # Bathrooms: No expert knowledge
    zbeta[3] ~ dnorm(0, 1 / 10^6)  # Non-informative prior

    # CarParks: Strong belief, centered at 120,000 AUD per car park
    zbeta[4] ~ dnorm(120000 / xsd[4], 1 / 10000^2)  # Small variance for strong belief

    # PropertyType: Strong belief, 150,000 AUD less for units
    zbeta[5] ~ dnorm(-150000 / xsd[5], 1 / 5000^2)  # Small variance for strong belief

    # Residual variance (inverse gamma)
    zVar ~ dgamma(0.01, 0.01)
    tau <- 1 / zVar

    # Transforming back to original scale
    beta0 <- zbeta0 * ysd
    for (j in 1:Nx) {
      beta[j] <- (zbeta[j] * ysd) / xsd[j]
    }
  }
"

# Write out the model to a text file
writeLines(modelString, con="BayesianModel.txt")
# Track start time
start_time <- Sys.time()

# Define the model string (your Bayesian model definition)
writeLines(modelString, con="TEMPmodel.txt")

# MCMC Parameters
#adaptSteps <- 500       # Number of steps to "tune" the samplers
#burnInSteps <- 1000     # Burn-in steps
#nChains <- 2            # Number of chains
#thinSteps <- 3         # Thinning interval
#numSavedSteps <- 10000 # Number of saved steps

#adaptSteps <- 2000       # Number of steps to "tune" the samplers (lower for simpler models)
#burnInSteps <- 5000      # Burn-in steps (check convergence diagnostics)
#nChains <- 3             # Number of chains (increase to 3 for better diagnostics)
#thinSteps <- 5           # Thinning interval (reduce if autocorrelation is low)
#numSavedSteps <- 5000    # Number of saved steps (adjust based on required precision)

#adaptSteps <- 5000        # Increase for more complex models
#burnInSteps <- 5000       # Adjust based on diagnostics
#nChains <- 4              # Increase chains for better convergence
#thinSteps <- 1            # Adjust based on autocorrelation (lower for low autocorrelation)
#numSavedSteps <- 10000    # Increase to get more precise estimates

adaptSteps <- 10000        # Increased adaptation for complex models
burnInSteps <- 10000       # Increased burn-in to ensure convergence
nChains <- 4               # Four chains should be sufficient, increase if needed
thinSteps <- 1             # Keeping thinning at 1 based on low autocorrelation
numSavedSteps <- 15000     # More saved steps for higher precision
nIter <- ceiling((numSavedSteps * thinSteps) / nChains)  # Total iterations

# Define initial values for the chains (optional)
initsList <- list(
  list(zbeta0 = 0, zbeta = rep(0, 5), zVar = 1),
  list(zbeta0 = 1, zbeta = rep(1, 5), zVar = 2)
)

# Data for the model (assuming dataList has been defined earlier)
dataList <- list(
  SalePrice = y,
  x = x,
  Nx = ncol(x),
  Ntotal = nrow(data)
)

# Run the model using run.jags
runJagsOut <- run.jags(method = "parallel",   # Use parallel processing
                       model = "TEMPmodel.txt",  # JAGS model file
                       monitor = c("zbeta0", "zbeta", "beta0", "beta", "tau", "zVar"),  # Parameters to monitor
                       data = dataList,         # Data for the model
                       inits = initsList,       # Initial values
                       n.chains = nChains,      # Number of chains
                       adapt = adaptSteps,      # Adaptation steps
                       burnin = burnInSteps,    # Burn-in steps
                       sample = numSavedSteps,  # Number of saved samples
                       thin = thinSteps,        # Thinning interval
                       summarise = FALSE,       # Disable automatic summary
                       plots = FALSE)           # Disable automatic plotting

# Convert to coda object for further analysis
codaSamples <- as.mcmc.list(runJagsOut)

# Display summary statistics of the sampled parameters
summary(codaSamples)

# Further thinning from the codaSamples
furtherThin <- 7  # Additional thinning interval
thinningSequence <- seq(1, nrow(codaSamples[[1]]), furtherThin)

# Apply further thinning
newCodaSamples <- mcmc.list()
for (i in 1:nChains) {
  newCodaSamples[[i]] <- as.mcmc(codaSamples[[i]][thinningSequence, ])
}

# Display summary statistics after further thinning
summary(newCodaSamples)

# Diagnostics for MCMC chains
param_names <- c("beta0", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "tau")
for (param in param_names) {
  diagMCMC(newCodaSamples, parName = param)
}

# Save the workspace for future use
save.image(file = "MCMCSamplingResults.RData")

# Track end time
end_time <- Sys.time()


# Print start and end times
print(paste("MCMC sampling started at:", start_time))
print(paste("MCMC sampling ended at:", end_time))

# Calculate and print total duration
duration <- end_time - start_time
print(paste("Total duration of MCMC sampling:", duration))

# Summarize the MCMC Results
summaryInfo <- smryMCMC_HD(codaSamples = newCodaSamples)
print(summaryInfo)


# Plot MCMC results using the custom plot function
plotMCMC_HD(codaSamples = newCodaSamples, 
            data = data, 
            xName = c("Area", "Bedrooms", "Bathrooms", "CarParks", "PropertyType"), 
            yName = "SalePrice", 
            saveName = "BayesianModel", saveType = "jpg")




# Conduct a predictive check
coefficients <- summaryInfo[7:11,3]  # Get the model coefficients from posterior
Variance <- summaryInfo[12,3]        # Get the variance

# Generate random data from the posterior distribution
meanGamma <- as.matrix(cbind(rep(1, nrow(x)), x)) %*% as.vector(coefficients)
randomData <- rgamma(n = nrow(x), shape = meanGamma^2 / Variance, rate = meanGamma / Variance)

# Display the density plot of observed data and posterior predictive distribution
predicted <- data.frame(SalePrice = randomData)
observed <- data.frame(SalePrice = y)
predicted$type <- "Predicted"
observed$type <- "Observed"
dataPred <- rbind(predicted, observed)

# Plot observed vs predicted densities
ggplot(dataPred, aes(SalePrice, fill = type)) + geom_density(alpha = 0.2)