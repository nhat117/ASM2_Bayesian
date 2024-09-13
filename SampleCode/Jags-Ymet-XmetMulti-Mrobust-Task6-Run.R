graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)

setwd("~/Documents/MATH2269_Bayesian/tasks/Module6")
myData = read.csv("~/Documents/MATH2269_Bayesian/tasks/Module6/concrete.csv")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # This one does the same job
myData = read.csv("concrete.csv")

yName = "CCS" ; xName = c("Cement",	"Blast.Furnace.Slag",	"Fly.Ash", 	"Water",  	"Superplasticizer", 	
                          "Coarse.Aggregate", 	"Fine.Aggregate" )
fileNameRoot = "Task6"

library(ggplot2)
head(myData)

hist(myData$CCS)
plot(kde(myData$CCS))

# Scatter plots
p1 <- ggplot(myData, aes(x=Cement, y=CCS)) +
  geom_point()

p2 <- ggplot(myData, aes(x=Blast.Furnace.Slag, y=CCS)) +
  geom_point()

p3 <- ggplot(myData, aes(x=Fly.Ash, y=CCS)) +
  geom_point()

p4 <- ggplot(myData, aes(x=Water, y=CCS)) +
  geom_point()

p5 <- ggplot(myData, aes(x=Superplasticizer, y=CCS)) +
  geom_point()

p6 <- ggplot(myData, aes(x=Coarse.Aggregate, y=CCS)) +
  geom_point()

p7 <- ggplot(myData, aes(x=Fine.Aggregate, y=CCS)) +
  geom_point()

figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, nrow = 4, ncol = 2)
figure

summary(myData)

numSavedSteps = 1000 
thinSteps = 7
nChains = 2

graphFileType = "eps" 


#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-XmetMulti-Mrobust-Task6.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
startTime = proc.time()
xPred = c(168 ,	42.1 ,	163.8 ,	121.8 ,	5.7 ,	1058.7 ,	780.1  )

mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , 
                    numSavedSteps=numSavedSteps , thinSteps=thinSteps , 
                    nChains = nChains , xPred = xPred )
stopTime = proc.time()
duration = stopTime - startTime
show(duration)

# save.image(file='Task6Chains.RData')
# numSavedSteps <- 10000; thinSteps <- 27; nChains <- 2; #Run time: 175.974 sec
# adaptSteps = 5000; burnInSteps = 5000 # These are controlled in "Jags-Ymet-XmetMulti-Mrobust-Task6.R" 
# load('Task6Chains.RData') # Longer run with better diagnostics


#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
graphics.off()
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:

summaryInfo = smryMCMC( mcmcCoda , 
                        saveName=fileNameRoot  )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 

summaryInfo[2:9,3]
summaryInfo[10,3]

# ============ Predictive check ============
library(LaplacesDemon) # To generate data from three parameter t-distribution

coefficients <- summaryInfo[2:9,3] # Get the model coefficients out
Variance <- summaryInfo[10,3]^2 # Get the variance out
nu <- summaryInfo[20,3] # Get the degrees fof freedom out

# Since we imposed the regression model on the mean of the student-t likelihood,
# we use the model (X*beta) to generate the mean of student-t population for each 
# observed x vector. 
x <- myData[,xName]
y <- myData[,yName]
mean <- as.matrix(cbind(rep(1,nrow(x)),  x)) %*% as.vector(coefficients)
# Generate random data from the posterior distribution. Here I take the 
# reparameterisation back to alpha and beta.
randomData <- rstp(n = nrow(myData), mu=mean, tau=1/Variance, nu=nu) # The second parameter is precision! See: https://search.r-project.org/CRAN/refmans/LaplacesDemon/html/dist.Student.t.Precision.html
  
# Display the density plot of observed data and posterior distribution:
predicted <- data.frame(CSS = randomData)
observed <- data.frame(CSS = y)
predicted$type <- "Predicted"
observed$type <- "Observed"
dataPred <- rbind(predicted, observed)

ggplot(dataPred, aes(CSS, fill = type)) + geom_density(alpha = 0.2)

