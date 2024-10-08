## Scenario:Broiler chickens following POMW 50mg/kg for 10 doses wihthin 5 days

## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(tidyr)       # R-package for tidy messy data
library(tidyverse)   # R-package for tidy messy data
library(truncnorm)   # R package for Truncated normal distribution
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(ggpubr)      # R package for plotting the data
library(truncnorm)   # R package for Truncated normal distribution
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(ggpubr)      # R package for plotting the data
library(ggprism)
library(ggplot2)
library(tibble)
library(rlang)
library(patchwork)

setwd("C:/Users/zhichengzhang/Desktop/Tyloasin_08_18/Figures/Figure 7_Sensitivity analysis/Broiler chicken/POMW")
pars <- log(c(
  Kst = 1.47,
  Kint = 2.10,
  KurineC = 1.38 , 
  KmC = 0.004 ,
  Ka = 0.48,
  
  QCC = 10.17,
  QLC = 0.35,
  QKC = 0.25,
  QFC = 0.015,
  QMC = 0.35,
  
  BW = 2.2,
  VLC = 0.0204,
  VKC = 0.0057,
  VFC = 0.05,
  VMC = 0.4015,
  VbloodC = 0.0483,
  Htc = 0.307,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  PF       = 2.73,
  PRest    = 6.20 
))


OUT<- pred(pars, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)

length(pars)


## Define the sensitivity function for AUC_V
NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    R.G            <- pred(Gpars, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_V
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_V
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(pars, pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "BroilersensV.csv")

## Define the sensitivity function for AUC_L
pars <- log(c(
  Kst = 1.47,
  Kint = 2.10,
  KurineC = 1.38 , 
  KmC = 0.004 ,
  Ka = 0.48,
  
  QCC = 10.17,
  QLC = 0.35,
  QKC = 0.25,
  QFC = 0.015,
  QMC = 0.35,
  
  BW = 2.2,
  VLC = 0.0204,
  VKC = 0.0057,
  VFC = 0.05,
  VMC = 0.4015,
  VbloodC = 0.0483,
  Htc = 0.307,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  PF       = 2.73,
  PRest    = 6.20 
))

NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    R.G            <- pred(Gpars, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_L
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_L
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}


NSC_GCA<- NSC_func(pars, pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "BroilersensL.csv")

## Define the sensitivity function for AUC_K
pars <- log(c(
  Kst = 1.47,
  Kint = 2.10,
  KurineC = 1.38 , 
  KmC = 0.004 ,
  Ka = 0.48,
  
  QCC = 10.17,
  QLC = 0.35,
  QKC = 0.25,
  QFC = 0.015,
  QMC = 0.35,
  
  BW = 2.2,
  VLC = 0.0204,
  VKC = 0.0057,
  VFC = 0.05,
  VMC = 0.4015,
  VbloodC = 0.0483,
  Htc = 0.307,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  PF       = 2.73,
  PRest    = 6.20 
))


NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    R.G            <- pred(Gpars, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_K
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_K
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(pars, pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "BroilersensK.csv")

## Define the sensitivity function for AUC_M
pars <- log(c(
  Kst = 1.47,
  Kint = 2.10,
  KurineC = 1.38 , 
  KmC = 0.004 ,
  Ka = 0.48,
  
  QCC = 10.17,
  QLC = 0.35,
  QKC = 0.25,
  QFC = 0.015,
  QMC = 0.35,
  
  BW = 2.2,
  VLC = 0.0204,
  VKC = 0.0057,
  VFC = 0.05,
  VMC = 0.4015,
  VbloodC = 0.0483,
  Htc = 0.307,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  PF       = 2.73,
  PRest    = 6.20 
))

NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    R.G            <- pred(Gpars, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_M
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_M
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}


NSC_GCA<- NSC_func(pars, pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "BroilersensM.csv")

## Define the sensitivity function for AUC_F
pars <- log(c(
  Kst = 1.47,
  Kint = 2.10,
  KurineC = 1.38 , 
  KmC = 0.004 ,
  Ka = 0.48,
  
  QCC = 10.17,
  QLC = 0.35,
  QKC = 0.25,
  QFC = 0.015,
  QMC = 0.35,
  
  BW = 2.2,
  VLC = 0.0204,
  VKC = 0.0057,
  VFC = 0.05,
  VMC = 0.4015,
  VbloodC = 0.0483,
  Htc = 0.307,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  PF       = 2.73,
  PRest    = 6.20 
))

NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    R.G            <- pred(Gpars, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_F
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_F
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(pars, pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "BroilersensF.csv")


library(viridis)
library(RColorBrewer)
## load package "gplots" containing "heatmap.2"
library(gplots)
library(ggplot2)
library(reshape2)

## loading data saved as ".CSV" format
IT_Sensitivity_Analysis <- read.csv("Broilersens.csv")

## name the loaded data and covert it to a data table by "data.frame"
ttt <- IT_Sensitivity_Analysis

# Assuming your data is in a dataframe called 'ttt'
# Convert your data from wide to long format
ttt_long <- melt(ttt, id.vars = "X")

# Now, plot the heatmap with a better color scale
heatmap3 <- ggplot(ttt_long, aes(x = variable, y = X, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = brewer.pal(11, "RdYlBu")[1], 
                       high = brewer.pal(11, "RdYlBu")[11], 
                       mid = "darkslategray3", # Shiny middle color
                       midpoint = 0,   # Value at which the mid color should be placed
  ) + # Set your limits accordingly
  labs(title = "(A) Normalized Sensitivity Coefficients for broiler chicken (POMW)", x = "Target tissues", y = "Parameters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "right")
heatmap3

