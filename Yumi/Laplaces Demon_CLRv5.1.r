rm(list=ls())

library(LaplacesDemon)
library(MCMCpack)
library(compositions)
source("Simulator.v5.1.r")

dat <- read.csv("krill.csv", header = T)

dH_obs <- as.vector(dat$d2H)
dC_obs <- as.vector(dat$d13C)
dN_obs <- as.vector(dat$d15N)
dH_w <- dat$d2H_w
tau <- dat$TL
gr <- dat$gr #grouping variable - Genus
N <- nrow(dat) 

# Simulation model
parm.names <-     c("dC_M", "dC_F", "dC_C", "dH_M", "dH_F", "dH_C", "phi_modM","phi_modF", "phi_modC","logvar_CM", "logvar_CF",  "logvar_CC",  "logvar_HM", "logvar_HF","logvar_HC", "logvar_DeltaC", "logvar_omega", "logomega", "logDelta_C", "logdH_p")
Initial.Values <- c(-23.7, -24.9, -19.1, -7.1, -113, -15, 1.626544, -14.256495, 7.616156, 2.044173, 1.728355, -8.700835, -6.765336, 7.443521, -6.809056, -6.449024, -5.089111, -1.374354, 0.875469, 5.106157)
Sim <- Model_Sim(Initial.Values, N, seed = 42)
dH_sim <- dH_obs <- Sim$dH_sim; dH_sim
dC_sim <- dC_obs <- Sim$dC_sim; dC_sim

mon.names <- c("LogPrior","LogLike","LogPosterior", "phi_M", "phi_F", "phi_C", "omega", "dH_p", "Delta_C", "u_C", "u_H", "sigma_C", "sigma_H")
parm.names <-     c("dC_M", "dC_F", "dC_C", "dH_M", "dH_F", "dH_C", "phi_modM","phi_modF", "phi_modC","logvar_CM", "logvar_CF",  "logvar_CC",  "logvar_HM", "logvar_HF","logvar_HC", "logvar_DeltaC", "logvar_omega", "logomega", "logDelta_C", "logdH_p")
Initial.Values <- c(-23.7, -24.9, -19.1, -7.1, -113, -15, 1.626544, -14.256495, 7.616156, 2.044173, 1.728355, -8.700835, -6.765336, 7.443521, -6.809056, -6.449024, -5.089111, -1.374354, 0.875469, 5.106157)
#Initial.Values <- c ( -23.7, -24.9, -19.1, -7.1,   -113,   -15,        3,       1,            7,      log(1),       log(1),        log(1),      log(1),         log(1),     log(1),       log(1),          log(1),      log(.2),    log(0.5),     log(163))
Data <- list(N = N, mon.names = mon.names, parm.names = parm.names, dH_obs = dH_obs, dC_obs = dC_obs, dH_w = dH_w, tau = tau)


Model <- function(parm, Data)
{
  # DATA
  dH_obs <- Data$dH_obs # i = 1,..,n
  dC_obs <- Data$dC_obs # i = 1,..,n
  tau  <- Data$tau  # trophic level for each individual i= 1,...,n
  dH_w <- Data$dH_w # hydrogen isotope measurement of environmental water for each individual
  # PARAMETERS
  dC_M <- parm[1]
  dC_F <- parm[2]
  dC_C <- parm[3]
  dH_M <- parm[4]
  dH_F <- parm[5]
  dH_C <- parm[6]
  phi <- parm[7:9] #centered-log ratio transformed source proportions, parameters of interest
  var   <- exp(parm[10:17])
  omega <- exp(parm[18]) #proportion of hydrogen composition due to environmental water dH_w
  Delta_C <- exp(parm[19]) #change in carbon due to trophic level
  dH_p  <- exp(parm[20]) #change in hydrogen between water and trophic level 1
  # TRANSFORMATIONS
  phi_T <- clr(phi[1:3]) #centered-log transformation of phi
  phis <- exp(phi_T)/sum(exp(phi_T)) #back tranform to mixture proportion
  phi_M <- phis[1] # marine source proportion
  phi_F <- phis[2] # freshwater source proportion
  phi_C <- phis[3] # coastal source proportion
  # PRIORS
  dC_M_prior <- dnorm(dC_M, 23.7, 0.7) #mean SD of marine carbon (Kline 2010)
  dC_F_prior <- dnorm(dC_F, 24.9, 2.1) #mean of POM at tau ~ 1 in FW
  dC_C_prior <- dnorm(dC_C, 19.1, 1.17) # mean of POM at tau ~ 1 in coastal
  dH_M_prior <- dnorm(dH_M, 7.4, 1) #mean of dH_w of 3 marine water
  dH_F_prior <- dnorm(dH_F, 113.1, 10.9) #mean of dH_w of 4 glacier streams
  dH_C_prior <- dnorm(dH_C, 15.4, 3.1)  #mean of 10 coastal water samples
  source_prior <- dC_M_prior + dC_F_prior + dC_C_prior + dH_M_prior + dH_F_prior + dH_C_prior 
  phi_modM_prior  <- dunif(phi_T[1], -10, 10, log = TRUE)
  phi_modF_prior  <- dunif(phi_T[2], -10, 10, log = TRUE)
  phi_modC_prior  <- dunif(phi_T[3], -10, 10, log = TRUE)
  phi_prior <- phi_modM_prior + phi_modF_prior + phi_modC_prior
  var_CM_prior <- dinvgamma(var[1], 0.001, 0.001) #variance of carbon marine source
  var_CF_prior <- dinvgamma(var[2], 0.001, 0.001) #variance of carbon freshwater source
  var_CC_prior <- dinvgamma(var[3], 0.001, 0.001) #variance of carbon coastal source
  var_HM_prior <- dinvgamma(var[4], 0.001, 0.001) #variance of hydrogen marine source
  var_HF_prior <- dinvgamma(var[5], 0.001, 0.001) #variance of hydrogen freshwater source
  var_HC_prior <- dinvgamma(var[6], 0.001, 0.001) #variance of hydrogen coastal source
  var_Delta_C_prior <- dinvgamma(var[7], 0.001, 0.001) #var carbon trophic fractionation
  var_dH_p_prior <- dinvgamma(var[8], 0.001, 0.001) #var hydrogen water-phyto fractionation
  var_prior <- var_CM_prior + var_CF_prior + var_CC_prior + var_HM_prior + var_HF_prior + var_HC_prior + var_Delta_C_prior + var_dH_p_prior
  omega_prior <- dnorm(omega, 0.23, 0.03, log = TRUE) 
  Delta_C_prior <- dnorm(Delta_C, 2.4, 0.5, log = TRUE)
  dH_p_prior <- dnorm(dH_p, 163.7, 27)
  # MODEL
#   C_tot <- Delta_C * (tau - 2)
#   u_C <- mean(phi_M*(dC_M + C_tot) + phi_F*(dC_F + C_tot) + phi_C*(dC_C + C_tot))
  u_C <- mean(phi_M*(dC_M + Delta_C) + phi_F*(dC_F + Delta_C) + phi_C*(dC_C + Delta_C))
  u_H <- mean((omega * dH_w) + (1 - omega)*(phi_M*(dH_M - dH_p) + phi_F*(dH_F - dH_p) + phi_C*(dH_C - dH_p)))
  sigma_C <- sqrt(phi_M*phi_M*(var[1] + var[7]) + phi_F*phi_F*(var[2] + var[7]) + phi_C*phi_C*(var[3] + var[7]))
  sigma_H <- sqrt(phi_M*phi_M*(var[4] + var[8]) + phi_F*phi_F*(var[5] + var[8]) + phi_C*phi_C*(var[6] + var[8]))
  epsilon_H <- dH_obs - u_H
  epsilon_C <- dC_obs - u_C
  # LIKELIHOOD
  C_like <- dnorm(epsilon_C, 0, sigma_C, log = TRUE)
  H_like <- dnorm(epsilon_H, 0, sigma_H, log = TRUE)
  LogLike <- sum(H_like + C_like)
  LogPrior <- source_prior + phi_prior + var_prior + omega_prior + Delta_C_prior + dH_p_prior
  LogPosterior <- LogLike + LogPrior
  # MODEL RETURN
  Modelout <- list(LP = LogPosterior, 
                   Dev = -2*LogLike,
                   Monitor = c(LogPrior, LogLike, LogPosterior, phi_M, phi_F, phi_C, omega, dH_p, Delta_C, u_C, u_H,  sigma_C, sigma_H),
                   yhat = NULL, 
                   parm = parm)
  return(Modelout)
}

set.seed(42)
Fit <- LaplacesDemon(Model, Data = Data, Initial.Values, Covar = NULL,
                     Iterations=100000, Status=10000, Thinning=100,
                     Algorithm = "HARM", Specs = list(alpha.star = .234, B = NULL))
Initial.Values <- as.initial.values(Fit)
plot(Fit, BurnIn=100, Data, PDF=F, Parms=NULL, ask = F)
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

