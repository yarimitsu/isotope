# Clean up R
rm(list=ls())

# R libraries
library(LaplacesDemon)
library(MCMCpack)
library(compositions)

# R code
source("Simulator.R")
source("Model.R")

# =====================================================================#
# DATA
# =====================================================================#
dat <- read.csv("data/krill.csv", header = TRUE)

dH_obs <- as.vector(dat$d2H)
dC_obs <- as.vector(dat$d13C)
dN_obs <- as.vector(dat$d15N)
dH_w <- dat$d2H_w
tau <- dat$TL
gr <- dat$gr # grouping variable - Genus
N <- nrow(dat) 

# =====================================================================#
# SIMULATION
# =====================================================================#
Initial.Values <- c(-23.7, -24.9, -19.1, -7.1, -113, -15, 1.626544, -14.256495, 7.616156, 2.044173, 1.728355, -8.700835, -6.765336, 7.443521, -6.809056, -6.449024, -5.089111, -1.374354, 0.875469, 5.106157)
Sim <- Model_Sim(Initial.Values, N, seed = 42)
dH_sim <- dH_obs <- Sim$dH_sim; dH_sim
dC_sim <- dC_obs <- Sim$dC_sim; dC_sim


# =====================================================================#
# MCMC
# =====================================================================#
parm.names <- c("dC_M", "dC_F", "dC_C", "dH_M", "dH_F", "dH_C", "phi_modM","phi_modF", "phi_modC","logvar_CM", "logvar_CF",  "logvar_CC",  "logvar_HM", "logvar_HF","logvar_HC", "logvar_DeltaC", "logvar_omega", "logomega", "logDelta_C", "logdH_p")
mon.names <- c("LogPrior","LogLike","LogPosterior", "phi_M", "phi_F", "phi_C", "omega", "dH_p", "Delta_C", "u_C", "u_H", "sigma_C", "sigma_H")
Initial.Values <- c(-23.7, -24.9, -19.1, -7.1, -113, -15, 1.626544, -14.256495, 7.616156, 2.044173, 1.728355, -8.700835, -6.765336, 7.443521, -6.809056, -6.449024, -5.089111, -1.374354, 0.875469, 5.106157)
#Initial.Values <- c ( -23.7, -24.9, -19.1, -7.1,   -113,   -15,        3,       1,            7,      log(1),       log(1),        log(1),      log(1),         log(1),     log(1),       log(1),          log(1),      log(.2),    log(0.5),     log(163))
Data <- list(N = N, mon.names = mon.names, parm.names = parm.names, dH_obs = dH_obs, dC_obs = dC_obs, dH_w = dH_w, tau = tau)

set.seed(42)
Fit <- LaplacesDemon(Model, Data = Data,
                     Initial.Values = Initial.Values, Covar = NULL,
                     Iterations = 100000, Status = 10000, Thinning = 100,
                     Algorithm = "HARM",
                     Specs = list(alpha.star = .234))
                     #Algorithm = "HARM", Specs = list(alpha.star = .234, B = NULL))
Initial.Values <- as.initial.values(Fit)
plot(Fit, BurnIn = 100, Data, PDF = FALSE, Parms = NULL, ask = FALSE)
Consort(Fit)


#########mean SD posterior predictive plot
labN<-expression(paste(delta^{15}, "N (\u2030)"))
labC<-expression(paste(delta^{13}, "C (\u2030)"))
labH<-expression(paste(delta^{2}, "H (\u2030)"))

Burn <- 200

Posterior_pred <- matrix(NA, nrow = N, ncol = 2)
colnames(Posterior_pred) <- c("dC", "dH")
Posterior_pred[,"dH"] <- rnorm(N, Fit$Monitor[(Burn + 1):nrow(Fit$Monitor),"u_H"], Fit$Monitor[(Burn + 1):nrow(Fit$Monitor),"sigma_H"])
Posterior_pred[,"dC"] <- rnorm(N, Fit$Monitor[(Burn + 1):nrow(Fit$Monitor),"u_C"], Fit$Monitor[(Burn + 1):nrow(Fit$Monitor),"sigma_C"])
dim(Fit$Monitor)

par(mar = c(4, 5, 4, 3), ask = F)
plot(dC_obs, dH_obs, xlim = range(Posterior_pred[,"dC"], dC_obs, Fit$Monitor[(Burn + 1):nrow(Fit$Monitor),"u_C"]), ylim = range(Posterior_pred[,"dH"], dH_obs, Fit$Monitor[(Burn + 1):nrow(Fit$Monitor),"u_H"]), pch = 16, col = "grey", ylab = labH, xlab = labC, main = "Data vs. Mean(SD) Posterior Predictive Distribution")
points(mean(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]), pch = 16, cex = 1.2)
arrows(mean(Posterior_pred[,"dC"]) - sd(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]), mean(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]), angle = 90, code = 1, length = 0.05)
arrows(mean(Posterior_pred[,"dC"]) + sd(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]), mean(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]), angle = 90, code = 1, length = 0.05)
arrows(mean(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]) - sd(Posterior_pred[,"dH"]), mean(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]), angle = 90, code = 1, length = 0.05)
arrows(mean(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]) + sd(Posterior_pred[,"dH"]), mean(Posterior_pred[,"dC"]), mean(Posterior_pred[,"dH"]), angle = 90, code = 1, length = 0.05)

