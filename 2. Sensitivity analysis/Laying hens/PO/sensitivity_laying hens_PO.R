## Scenario:laying hens following oral gavage administration 50mg/kg for 10 consecutive days wihthin 5 days

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

setwd("C:/Users/zhichengzhang/Desktop/Tyloasin_08_18/Figures/Figure 7_Sensitivity analysis/Laying hens/PO")

pars <- c(
  Kst = 2.11,
  Kint = 1.61,
  KurineC = 2.10, 
  KmC = 0.007,
  Ka = 0.24,
  
  QCC = 10.17,
  QLC = 0.2526,
  QKC = 0.2012,
  QFC = 0.015,
  QMC = 0.20,
  QovaryC = 0.1646,
  QoviductC = 0.1381,
  
  BW = 2.2,
  VLC = 0.0249,
  VKC = 0.0076,
  VFC = 0.05,
  VMC = 0.40,
  VbloodC = 0.0483,
  Htc = 0.307,
  VovaryC = 0.0191,
  VoviductC = 0.0258,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  
  Povary = 0.71 ,
  Poviduct = 1.40 ,
  PF       = 2.73,
  PRest    = 6.20,    
  
  Ay = 25,
  Ky = 0.02 ,                        
  Kw = 0.09  ,                         
  tlag = 1*24,                      
  tsig = 2*24,                  
  s    = 1/24,                    
  talbumen = 10,        
  tlay1  = 12,               
  Rwhitefor = (34*0.001)/10
)


OUT<- pred(log(pars),tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)

length(pars)

## Define the sensitivity function for AUC_V
NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    R.G            <- pred(Gpars, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_V
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_V
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}


NSC_GCA<- NSC_func(log(pars), pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "poultrysensV.csv")


## Define the sensitivity function for AUC_L
pars <- c(
  Kst = 2.11,
  Kint = 1.61,
  KurineC = 2.10, 
  KmC = 0.007,
  Ka = 0.24,
  
  QCC = 10.17,
  QLC = 0.2526,
  QKC = 0.2012,
  QFC = 0.015,
  QMC = 0.20,
  QovaryC = 0.1646,
  QoviductC = 0.1381,
  
  BW = 2.2,
  VLC = 0.0249,
  VKC = 0.0076,
  VFC = 0.05,
  VMC = 0.40,
  VbloodC = 0.0483,
  Htc = 0.307,
  VovaryC = 0.0191,
  VoviductC = 0.0258,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  
  Povary = 0.71 ,
  Poviduct = 1.40 ,
  PF       = 2.73,
  PRest    = 6.20,    
  
  Ay = 25,
  Ky = 0.02 ,                        
  Kw = 0.09  ,                         
  tlag = 1*24,                      
  tsig = 2*24,                  
  s    = 1/24,                    
  talbumen = 10,        
  tlay1  = 12,               
  Rwhitefor = (34*0.001)/10
)
NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    R.G            <- pred(Gpars, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_L
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_L
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}
  
NSC_GCA<- NSC_func(log(pars), pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "poultrysensL.csv")

## Define the sensitivity function for AUC_K
pars <- c(
  Kst = 2.11,
  Kint = 1.61,
  KurineC = 2.10, 
  KmC = 0.007,
  Ka = 0.24,
  
  QCC = 10.17,
  QLC = 0.2526,
  QKC = 0.2012,
  QFC = 0.015,
  QMC = 0.20,
  QovaryC = 0.1646,
  QoviductC = 0.1381,
  
  BW = 2.2,
  VLC = 0.0249,
  VKC = 0.0076,
  VFC = 0.05,
  VMC = 0.40,
  VbloodC = 0.0483,
  Htc = 0.307,
  VovaryC = 0.0191,
  VoviductC = 0.0258,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  
  Povary = 0.71 ,
  Poviduct = 1.40 ,
  PF       = 2.73,
  PRest    = 6.20,    
  
  Ay = 25,
  Ky = 0.02 ,                        
  Kw = 0.09  ,                         
  tlag = 1*24,                      
  tsig = 2*24,                  
  s    = 1/24,                    
  talbumen = 10,        
  tlay1  = 12,               
  Rwhitefor = (34*0.001)/10
)

NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    R.G            <- pred(Gpars, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_K
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_K
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(log(pars), pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "poultrysensK.csv")

## Define the sensitivity function for AUC_M
pars <- c(
  Kst = 2.11,
  Kint = 1.61,
  KurineC = 2.10, 
  KmC = 0.007,
  Ka = 0.24,
  
  QCC = 10.17,
  QLC = 0.2526,
  QKC = 0.2012,
  QFC = 0.015,
  QMC = 0.20,
  QovaryC = 0.1646,
  QoviductC = 0.1381,
  
  BW = 2.2,
  VLC = 0.0249,
  VKC = 0.0076,
  VFC = 0.05,
  VMC = 0.40,
  VbloodC = 0.0483,
  Htc = 0.307,
  VovaryC = 0.0191,
  VoviductC = 0.0258,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  
  Povary = 0.71 ,
  Poviduct = 1.40 ,
  PF       = 2.73,
  PRest    = 6.20,    
  
  Ay = 25,
  Ky = 0.02 ,                        
  Kw = 0.09  ,                         
  tlag = 1*24,                      
  tsig = 2*24,                  
  s    = 1/24,                    
  talbumen = 10,        
  tlay1  = 12,               
  Rwhitefor = (34*0.001)/10
)
NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    R.G            <- pred(Gpars, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_M
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_M
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(log(pars), pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "poultrysensM.csv")

## Define the sensitivity function for AUC_F
pars <- c(
  Kst = 2.11,
  Kint = 1.61,
  KurineC = 2.10, 
  KmC = 0.007,
  Ka = 0.24,
  
  QCC = 10.17,
  QLC = 0.2526,
  QKC = 0.2012,
  QFC = 0.015,
  QMC = 0.20,
  QovaryC = 0.1646,
  QoviductC = 0.1381,
  
  BW = 2.2,
  VLC = 0.0249,
  VKC = 0.0076,
  VFC = 0.05,
  VMC = 0.40,
  VbloodC = 0.0483,
  Htc = 0.307,
  VovaryC = 0.0191,
  VoviductC = 0.0258,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  
  Povary = 0.71 ,
  Poviduct = 1.40 ,
  PF       = 2.73,
  PRest    = 6.20,    
  
  Ay = 25,
  Ky = 0.02 ,                        
  Kw = 0.09  ,                         
  tlag = 1*24,                      
  tsig = 2*24,                  
  s    = 1/24,                    
  talbumen = 10,        
  tlay1  = 12,               
  Rwhitefor = (34*0.001)/10
)

NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    R.G            <- pred(Gpars, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_F
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_F
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(log(pars), pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "poultrysensF.csv")

## Define the sensitivity function for AUC_Ovary
pars <- c(
  Kst = 2.11,
  Kint = 1.61,
  KurineC = 2.10, 
  KmC = 0.007,
  Ka = 0.24,
  
  QCC = 10.17,
  QLC = 0.2526,
  QKC = 0.2012,
  QFC = 0.015,
  QMC = 0.20,
  QovaryC = 0.1646,
  QoviductC = 0.1381,
  
  BW = 2.2,
  VLC = 0.0249,
  VKC = 0.0076,
  VFC = 0.05,
  VMC = 0.40,
  VbloodC = 0.0483,
  Htc = 0.307,
  VovaryC = 0.0191,
  VoviductC = 0.0258,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  
  Povary = 0.71 ,
  Poviduct = 1.40 ,
  PF       = 2.73,
  PRest    = 6.20,    
  
  Ay = 25,
  Ky = 0.02 ,                        
  Kw = 0.09  ,                         
  tlag = 1*24,                      
  tsig = 2*24,                  
  s    = 1/24,                    
  talbumen = 10,        
  tlay1  = 12,               
  Rwhitefor = (34*0.001)/10
)

NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    R.G            <- pred(Gpars, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_Ovary
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_Ovary
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(log(pars), pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "poultrysensOvary.csv")

## Define the sensitivity function for AUC_Oviduct
pars <- c(
  Kst = 2.11,
  Kint = 1.61,
  KurineC = 2.10, 
  KmC = 0.007,
  Ka = 0.24,
  
  QCC = 10.17,
  QLC = 0.2526,
  QKC = 0.2012,
  QFC = 0.015,
  QMC = 0.20,
  QovaryC = 0.1646,
  QoviductC = 0.1381,
  
  BW = 2.2,
  VLC = 0.0249,
  VKC = 0.0076,
  VFC = 0.05,
  VMC = 0.40,
  VbloodC = 0.0483,
  Htc = 0.307,
  VovaryC = 0.0191,
  VoviductC = 0.0258,
  
  PB = 0.2753, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  
  Povary = 0.71 ,
  Poviduct = 1.40 ,
  PF       = 2.73,
  PRest    = 6.20,    
  
  Ay = 25,
  Ky = 0.02 ,                        
  Kw = 0.09  ,                         
  tlag = 1*24,                      
  tsig = 2*24,                  
  s    = 1/24,                    
  talbumen = 10,        
  tlay1  = 12,               
  Rwhitefor = (34*0.001)/10
)

NSC_func <- function (Gpars, pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 1)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Rnew.G         <- pred(Gpars.new, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    R.G            <- pred(Gpars, tinterval = 12, TDOSE = 10, BW = 2.2, route = 'oral gavage', Dose =50, tlay1 = 12)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GCA.new       =  Rnew.G[Rnew.G$Time==120,]$AUC_Oviduct
    AUC.GCA.ori       =  R.G[R.G$Time==120,]$AUC_Oviduct
    
    delta.AUC.GCA     =  AUC.GCA.new - AUC.GCA.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GCA/AUC.GCA.ori) * delta.Gpars)
    
  }
  
  return (NSC_GCA = NSC_GCA)
}

NSC_GCA<- NSC_func(log(pars), pred)

pars[]=NSC_GCA[,1]

pars_01=as.data.frame(pars)

NSC_Table<- abs(cbind.data.frame(pars_01))

write.csv(pars_01, file = "poultrysensOviduct.csv")


library(viridis)
library(RColorBrewer)
## load package "gplots" containing "heatmap.2"
library(gplots)
## loading data saved as ".CSV" format
IT_Sensitivity_Analysis <- read.csv("poultrysens.csv")
## name the loaded data and covert it to a data table by "data.frame"
ttt <- IT_Sensitivity_Analysis

library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Assuming your data is in a dataframe called 'ttt'
# Convert your data from wide to long format
ttt_long <- melt(ttt, id.vars = "X")

# Now, plot the heatmap with a better color scale
heatmap4 <- ggplot(ttt_long, aes(x = variable, y = X, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = brewer.pal(11, "RdYlBu")[1], 
                       high = brewer.pal(11, "RdYlBu")[11], 
                       mid = "darkslategray3", # Shiny middle color
                       midpoint = 0,   # Value at which the mid color should be placed
                       ) + # Set your limits accordingly
  labs(title = "(B) Normalized Sensitivity Coefficients for laying hens (PO)", x = "Target tissues", y = "Parameters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "right")
heatmap4
