# Install necessary libraries if not already installed
# install.packages("rjags")
# install.packages("coda")
#install.packages("psych")  # For descriptive statistics

# Load necessary libraries
library(rjags)
library(coda)
library(psych)  # For descriptive statistics

# Load your dataset (you can modify the path as needed)
file_path <- './Assignment2PropertyPrices.csv'
data <- read.csv(file_path)

# Randomly sample 500 rows from the dataset
set.seed(123)  # Set seed for reproducibility
data <- data[sample(nrow(data), 500), ]

# Rename columns for easier reference
colnames(data) <- c('SalePrice', 'Area', 'Bedrooms', 'Bathrooms', 'CarParks', 'PropertyType')

# Convert SalePrice from 100K units to full AUD
data$SalePrice <- data$SalePrice * 100000

# Descriptive statistics
print(describe(data))  # Displays descriptive statistics of the dataset

#===============================================================================
# Function to generate MCMC samples using JAGS
genMCMC <- function(data, numSavedSteps = 1500, adaptSteps = 1500, burnInSteps = 5000, thinSteps = 11, nChains = 2, saveName = NULL) {
  require(rjags)
  
  # Prepare the data for the JAGS model
  jags_data <- list(
    N = nrow(data),  
    y = data$SalePrice,  
    Area = data$Area,  
    Bedrooms = data$Bedrooms,  
    Bathrooms = data$Bathrooms,  
    CarParks = data$CarParks,  
    PropertyType = data$PropertyType  # PropertyType is binary (0 or 1)
  )
  
  # Define the JAGS model as a string
  jags_model_code <- "
  model {
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)  # Normal distribution for SalePrice
      
      mu[i] <- beta0 + beta1 * Area[i] + beta2 * Bedrooms[i] + beta3 * Bathrooms[i] +
               beta4 * CarParks[i] + beta5 * PropEffect[i]
      
      # PropertyType is a Bernoulli variable (categorical with 0 and 1)
      PropEffect[i] <- beta5 * PropertyType[i]
    }

    # Priors for regression coefficients
    beta0 ~ dnorm(0, 1.0E-6)
    beta1 ~ dnorm(90, 0.01)
    beta2 ~ dnorm(100000, 1.0E-4)
    beta3 ~ dnorm(0, 1.0E-6)
    beta4 ~ dnorm(120000, 0.01)
    beta5 ~ dnorm(-150000, 0.001)  # PropertyType effect modeled with normal distribution

    tau <- 1 / sigma2  # Precision (inverse of variance)
    sigma2 ~ dgamma(1, 1)  # Gamma distribution for variance
  }
  "
  
  # Write the model to a file
  writeLines(jags_model_code, con = "model.bug")
  
  # Compile the JAGS model
  model <- jags.model("model.bug", data = jags_data, n.chains = nChains, n.adapt = adaptSteps)
  
  # Burn-in phase
  cat("Burning in the MCMC chain...\n")
  update(model, n.iter = burnInSteps)
  
  # Sample MCMC
  cat("Sampling final MCMC chain...\n")
  codaSamples <- coda.samples(model, variable.names = c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "sigma2"),
                              n.iter = numSavedSteps * thinSteps / nChains, thin = thinSteps)
  
  if (!is.null(saveName)) {
    save(codaSamples, file = paste0(saveName, "_Mcmc.Rdata"))
  }
  
  return(codaSamples)
}

#===============================================================================
# Function to summarize MCMC samples
smryMCMC <- function(codaSamples, saveName = NULL) {
  summaryInfo <- summary(codaSamples)
  
  if (!is.null(saveName)) {
    write.csv(summaryInfo$statistics, file = paste0(saveName, "_Summary.csv"))
  }
  
  return(summaryInfo)
}

#===============================================================================
# Function to plot MCMC diagnostics for each parameter
plotMCMC <- function(codaSamples, data, saveName = NULL) {
  mcmcMat <- as.matrix(codaSamples, chains = TRUE)
  parNames <- colnames(mcmcMat)
  
  # Create diagnostic plots for each parameter
  for (parName in parNames) {
    diagMCMC(codaSamples, parName = parName)
  }
  
  if (!is.null(saveName)) {
    saveGraph(file = paste0(saveName, "_Diagnostics"), type = "jpg")
  }
}

#===============================================================================
# Function to generate diagnostic plots (trace, density, autocorrelation)
diagMCMC <- function(codaSamples, parName) {
  parValues <- as.matrix(codaSamples)[, parName]
  
  # Trace plot
  plot(parValues, type = "l", main = paste("Trace plot of", parName), ylab = parName)
  
  # Density plot
  plot(density(parValues), main = paste("Density plot of", parName), xlab = parName)
  
  # Autocorrelation plot
  acf(parValues, main = paste("Autocorrelation of", parName))
}

#===============================================================================
# Generate MCMC samples
codaSamples <- genMCMC(data, numSavedSteps = 1500, adaptSteps = 1500, burnInSteps = 5000, thinSteps = 11, nChains = 2)

# Summarize MCMC samples
summaryInfo <- smryMCMC(codaSamples)

# Diagnostic plots for beta0
diagMCMC(codaSamples, parName="beta0")

# Diagnostic plots for beta1
diagMCMC(codaSamples, parName="beta1")

# Diagnostic plots for beta2
diagMCMC(codaSamples, parName="beta2")

# Diagnostic plots for beta3
diagMCMC(codaSamples, parName="beta3")

# Diagnostic plots for beta4
diagMCMC(codaSamples, parName="beta4")

# Diagnostic plots for beta5
diagMCMC(codaSamples, parName="beta5")

# Diagnostic plots for sigma2
diagMCMC(codaSamples, parName="sigma2")
