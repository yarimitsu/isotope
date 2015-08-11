####IsotopeR code and directions updated by Jake Ferguson and Jack Hopkins: 9/30/11####
##IsotopeR developed from code provided by Semmens et al. (2009)
##3rd source is not estimated -- it has a fixed mean & SD.
##This code does not use the IsotopeR GUI. See R vignette for GUI instructions and link for example data. 
##Vignette can be viewed in R: 

vignette("IsotopeR")

####Directions for running R code (for example in text)####
##To run the IsotopeR code for the example in text, extract the folder from the zip file; then, do the following: 
##(1) Install the JAGS software from http://sourceforge.net/projects/mcmc-jags/ 
##(2) Load the following R packages: coda, lattice, R2WinBUGS, rjags
##You can install these packages with the following line:

install.packages("R2jags", dep=T)

##You may need to change the current working directory in R to the folder that you unzipped using:

setwd('path') ##Where 'path' is the path to the directory containing the Isotoper code. You can run the code by pasting the entire text from the IsotopeR.r file into the R console, or sourcing it by:
 
source(“IsotopeR.r”) 

##Once the code is running a counter will run on the bottom indicating the progress. Note: If convergence did not occur, you can increase the size of the mcmc chains by changing from mcmc.chainLength (line 55) <- 5e3 to 5e5.

library(R2jags) 

#Read in needed data -- this section of code will read in files that have been saved as .csv files in openoffice, libreoffice, or in ms office.
X 		      <- try(as.matrix(read.table("Data/ModelPaperBearsInd.csv",sep='\t',header=T)), silent=T) #data file
if(dim(X)[2] == 1) {
    X           <- as.matrix(read.table("Data/ModelPaperBearsInd.csv", sep=',', header=TRUE)) #mixture data file
}
Z 		      <- try(as.matrix(read.table("Data/Jack_mixerror.csv", sep='\t', header=T)), silent=T) #file for process error
if(dim(Z)[2] == 1) {
    Z           <- as.matrix(read.table("Data/Jack_mixerror.csv", sep=',', header=TRUE)) #mixture data file
}
sources 	      <- try(as.matrix(read.table("Data/ModelPaperSources.csv", sep='\t', header=T)), silent=T)
if(dim(sources)[2] == 1) {
    sources           <- as.matrix(read.table("Data/ModelPaperSources.csv", sep=',', header=TRUE)) #mixture data file
}
D 		      <- try(as.matrix(read.table("Data/ModelPaperSourceDigestCD.csv", sep='\t', header=T)), silent=T)
if(dim(D)[2] == 1) {
    D           <- as.matrix(read.table("Data/ModelPaperSourceDigestCD.csv", sep=',', header=TRUE)) #mixture data file
}
frac.mean	      <- try(as.matrix(read.csv("Data/fracMax.csv",sep='\t',header=T)), silent=T)#file with mean of fractionation in sources
if(dim(frac.mean)[2] == 1) {
    frac.mean           <- as.matrix(read.table("Data/fracMax.csv", sep=',', header=TRUE)) #mixture data file
}
frac.sd               <- try(as.matrix(read.csv("Data/Jack_datafracsd1.csv", sep='\t', header=F), colClasses="numeric"), silent=T)#file with standard deviation of fractionation in sources
if(dim(frac.sd)[2] == 1) {
    frac.sd           <- as.matrix(read.table("Data/Jack_datafracsd1.csv", sep=',', header=TRUE)) #mixture data file
}
image.name  <- "Data/IsotopeData.Rdata" #this saves the chains and all the model output into this image file

jags.params <- c("D.mu","D.tau","cov.source","corr.source","mu.source","p.pop","p","sigmaZ", "mix.mu", "res.err","pop.var","mu","CNratio")
model.loc   <- "IsotopeR.txt" #name and location of the model file passed to JAGS.

# These are the parameters that need to be set by the user
mcmc.chains		<- 3 #number of markov chains (each chain has a different starting value)
mcmc.burn 		<- 1e3 #length of chain discarded at the beginning of the run
mcmc.chainLength 	<- 5e3+mcmc.burn #total number of iterations per chain (includes burnin)
mcmc.thin 		<- 100 #thinning rate. reduces the sample size to every nth iteration. use if there is autocorrelation in the sample.


#digestibility data
source.mean <- c(52.83,6.88)
source.var <- c(2.45^2, 1.10^2)


#extract useful info from the read in data to pass to JAGS
N		<- dim(X)[1] #number of individuals in the sample
num.iso		<- dim(X)[2]-2 #number of isotopes in the sample
num.sources <- nlevels(as.factor(sources[,num.iso+1]))

#center Z (observation error) around 0.
Z[,1] <- Z[,1] - mean(Z[,1])
Z[,2] <- Z[,2] - mean(Z[,2])
Nz <- dim(Z)[1]

#prior parameters for sources        
mu.prior.mu <- c(-1,1) #source mean prior mean values
mu.prior.cov <- solve(diag(num.iso)*1000)#source mean prior covariance matrix

#prior parameters for concentrations
dmu.prior.tau <- diag(num.iso)*1/1000
dmu.prior.mu <- matrix(c(50,1),num.iso,num.sources)
dmu.prior.mu[2,2] <- 12

#population prior
alpha <- c(1,1,1)/3#prior proportions
alpha.clr <- log(alpha/prod(alpha)^(1/length(alpha))) #transform onto CLR scale

#priors for concentration covariance matrix
D.R <- diag(num.iso)*1/1000

#########Below here is data formatting for JAGS. Do not mess with this unless you understand what you are doing###########
#now put sources into an array makes it easier to use in JAGS. this is just processing stuff and not very important
source.ss <- vector('numeric', num.sources)
for(i in 1:num.sources) {
  source.ss[i] <- length(which(sources[,num.iso+1] == i))
}
source.array <- array(NA, c(max(source.ss),num.iso,num.sources))
for(i in 1:num.sources) {
  source.array[1:source.ss[i],,i] <- sources[which(sources[,num.iso+1] == i),1:num.iso]
}

#build array of istope concentrations as well.
cd.ss <- vector('numeric', num.sources)
for(i in 1:num.sources) {
  cd.ss[i] <- length(which(D[,num.iso+1] == i))
}
cd.array <- array(NA, c(max(cd.ss),num.iso,num.sources))
for(i in 1:num.sources) {
  cd.array[1:cd.ss[i],,i] <- D[which(D[,num.iso+1] == i),1:num.iso]
}

#builds fractionation variation matrices. These are added to estimators of the standard deviation in order to account for fractionation variation
fracvar.mat <- array(0,c(num.iso,num.iso,num.sources))
for(i in 1:num.sources) {
  fracvar.mat[,,i] <- diag(frac.sd[i,])
}

#get number of individual observations
num.inds <- nlevels(as.factor(X[,4]))
ind.counts <- vector('numeric',num.inds)
for(i in 1:num.inds) {
  ind.counts[i] <- length(which(i == X[,4]))
}

#individual observation matrix
ind.array <- array(NA,c(num.iso,num.inds,max(ind.counts)))
for(i in 1:num.inds) {
  ind.array[1:num.iso,i,] <- X[which(X[,4]==i),1:num.iso]
}
N <- num.inds

jags.data = list("ind.counts", "ind.array", "N", "num.sources", "num.iso", "X", "Z",  "dmu.prior.mu", "Nz", "source.array", "cd.array", "mu.prior.mu", "mu.prior.cov", "frac.sd", "dmu.prior.tau", "alpha.clr", "source.ss", "cd.ss", "source.mean", "source.var", "fracvar.mat")  #data passed into jags

#function used to initialize parameters
jags.inits <- function(){
    list(dmu.prior.mu=dmu.prior.mu, mu.prior.mu=mu.prior.mu, p.transform=runif(num.sources), region.sig=0.5, ind.sig=0.5, p.ind = matrix(runif(N*num.sources), N, num.sources))
}

#run jags
jags.1 = jags(jags.data, inits=jags.inits, parameters.to.save=jags.params, model.file=model.loc, n.chains=mcmc.chains, n.burnin=mcmc.burn, n.thin=mcmc.thin, n.iter=mcmc.chainLength, DIC=F)#, n.adapt = mcmc.adapt) 
print(jags.1) ##print out parameter estimates
save.image(image.name) #save image file with all the info from the MCMC 
