## Species: Broiler chicken
## Tylosin, a hydrosoluble compound, with a molecular weight of 916 g/mol, has a low oral bioavailability and a low apparent volume of distribution (https://doi.org/10.4315/0362-028X.JFP-13-440)

# Clear current environment -----------------------------------------------
rm(list=ls(all=TRUE))

getwd()
setwd("C:/Users/zhichengzhang/Desktop/Tyloasin_08_18/1. Pharmacokinetics data/Broiler chicken")

{
  library(readxl)
  library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' package
  library(magrittr)    # Package for the pipe, %>% 
  library(dplyr)       # Package for data manipulation and transformation using a grammar of data manipulation
  library(ggplot2)     # Package for plots
  library(FME)         # Package for MCMC simulation and model fitting
  library(minpack.lm)  # Package for model fitting
  library(reshape)     # Package for melt function to reshape the table
  library(truncnorm)   # Package for the truncated normal distribution function   
  library(EnvStats)    # Package for Environmental Statistics, including the US EPA Guidance
  library(invgamma)    # Package for inverse gamma distribution function
  library(foreach)     # Package for executing loop iterations in parallel, facilitating efficient computation
  library(doParallel)  # Package for parallel computing
  library(bayesplot)   # Package for MCMC traceplot
  library(tidyr)       # R-package for tidy messy data
  library(tidyverse)   # R-package for tidy messy data
  library(truncnorm)   # R package for truncated normal distribution
  library(EnvStats)    # Package for Environmental Statistics, including the US EPA Guidance
  library(ggpubr)      # R package for plotting the data
  library(readxl)      # Package for reading Excel files, supporting both .xls and .xlsx file formats
}


# Build mrgsolve-based PBPK Model -----------------------------------------
SolvePBPK <- '
$PARAM @annotated
 
// Oral absorption rate constants
Kst    : 2        :                   /h, gastric emptying rate constant, Lin et al. 2015
Kint   : 0.108    :                   /h, intestinal transit rate constant, Table 10, Wang et al. 2021

// IV infusion rate constants
Timeiv  : 0.01     :                  /h, IV injection/infusion time

// Urinary elimination rate constant converted from elimination half life
KurineC : 0.1199   :                  L/h/kg

// Metabolic rate constant adjusted by body weight
KmC     : 0.025        :             L/h/kg, metabolic rate constant , Li et al. 2017

// Dosing, multiple oral gavage
Ka      : 0.23      :                 /h, intestinal absorption rate constant, Ronaghinia et al. 2021

// Physiological parameters
// Blood flow rates
QCC     : 10.17      :                Cardiac output (L/h/kg), Table 13, Wang et al. 2021
QLC     : 0.35    :                 Fraction of blood flow to the liver, Table 16, Wang et al. 2021
QKC     : 0.25    :                 Fraction of blood flow to the kidneys, Table 16, Wang et al. 2021
QFC     : 0.015     :                 Fraction of blood flow to the fat, Zeng et al. 2019
QMC     : 0.35      :                 Fraction of blood flow to the muscle, Table 16, Wang et al. 2021

// Tissue  volumes
BW      : 2.2       :                 Body weight (kg)
VLC     : 0.0204    :                 Fractional liver weight, Table 1, Wang et al. 2021
VKC     : 0.0057    :                 Fractional kidney weight, Table 1, Wang et al. 2021
VFC     : 0.05      :                 Fractional fat weight, Table 1, Wang et al. 2021
VMC     : 0.4015    :                 Fractional muscle weight, Table 1, Wang et al. 20211
VbloodC : 0.0483    :                 Fractional blood volume, Table 24, Wang et al. 2021
Htc     : 0.307     :                 Hematocrit, Table 22, Wang et al. 2021

// Mass Transfer Parameters (Chemical-specific parameters)
// Partition coefficients (PC, tissue:plasma)
PB      : 0.324     :                Percentage of drug bound to plasma proteins
PL      : 10.60     :                Liver plasma PC, calculated using AUCtissue:AUCplasma method 
PK      : 13.92     :                Kidney plasma PC, calculated using AUCtissue:AUCplasma method 
PM      : 2.90      :                Muscle plasma PC, breast muscle, calculated using AUCtissue:AUCplasma method 
PF      : 2.90      :                Fat plasma PC, assumed to be the same as muscle
PRest   : 7.58      :                Rest of the body tissues:plasma PC, Yuan et al. 2022

$MAIN

// Cardiac output and blood flows to tissues/organs (L/h)
double QC = QCC * BW;                  // Cardiac output
double QL = QLC * QC;                  // Liver
double QK = QKC * QC;                  // Kidney
double QF = QFC * QC;                  // Fat
double QM = QMC * QC;                  // Muscle
double QRest = QC - QK - QL - QF - QM; // Rest of the tissues


//  Tissue volumes (Kg)
double VL = VLC * BW;                  // Liver
double VK = VKC * BW;                  // Kidney
double VF = VFC * BW;                  // Fat
double VM = VMC * BW;                  // Muscle
double Vblood = VbloodC*BW;            // Blood
double Vplasma = VplasmaC * BW;        // Plasma
double VRest = 1 * BW - VL - VK - VF - VM - Vblood; // Rest of the tissues
double VplasmaC = (1-Htc) * VbloodC;   // Plasma volume, fraction of BW, Table 22, Wang et al. 2021
double Kmet = KmC * BW;	               // L/h Metabolic rate constant adjusted by body weight
double Kurine = KurineC * BW;          // L/h Urinary elimination rate
double Free = 1-PB;                    //Free drug percentage

      
$CMT AIV ADOSE AST AI Acolon AAO Aplas_free AL AK Aurine Amet AM AUCCV AUCCL AUCCK AUCCM AUCCF AF ARest

$ODE

dxdt_ADOSE=0;

 // The changing rate of the amount of dose via oral
 
 // Concentration of the chemical in vein and tissue compartments
      double CVL = AL/(VL * PL);      // Concentration in liver / PC of liver
      double CVK = AK/(VK * PK);      // Concentration in Kidney / PC of kidney
      double CVF = AF/(VF * PF);   // Concentration in fat / PC of fat
      double CVRest = ARest/(VRest * PRest); // Concentration in Rest  / PC of rest of tissues
      double CVM = AM/(VM * PM);      // Concentration in muscle / PC of muscle 
      double CV = ((QL * CVL  + QK  * CVK  + QF * CVF  + QM * CVM  + QRest * CVRest)/QC);
      double Cplas_free = Aplas_free/Vplasma;  // Free drug concentration in plasma = amount of free drug in plasma / volume of plasma
      double Cplasma = Cplas_free/Free ; // drug concentration in plasma 
      double CL = AL/VL;              // Concentration in liver
      double CK = AK/VK;              // Concentration in kidney
      double CM = AM/VM;              // Concentration in muscle
      double CF = AF/VF;              // Concentration in fat 
      double CRest = ARest/VRest;     // Concentration in rest of body

  // Rate of concentration changes in different compartments    
    double RAST = - Kst * AST;        // oral to the gastric 
    double RAI = Kst * AST - Kint * AI - Ka * AI; // Absorption in the intestine 
    double Rcolon = Kint * AI ;       // Intestinal transit to Colon 
    double RAO = Ka * AI; // Drug absorbed 
    double Rplas_free = QC * (CV - Cplasma) * Free; // rate of change in amount of the drug in plasma of tissues
    double RL = QL * (Cplasma - CVL)* Free + RAO - Kmet * AL ;// rate of change in amount of the drug in liver
    double Rmet = Kmet * AL;          // Rmet the metabolic rate in liver (mg/h)
    double Rurine = Kurine * CVK;   // Rurine the elimination rate in liver (mg/h)
    double RK = QK * (Cplasma - CVK) * Free - Rurine; // rate of change in amount of the drug in kidney
    double RM = QM * (Cplasma - CVM) * Free;  // rate of change in amount of the drug in muscle
    double RF = QF * (Cplasma - CVF) * Free;  // rate of change in amount of the drug in fat
    double RRest = QRest * (Cplasma - CVRest) * Free; // rate of change in amount of the drug in rest of the body
   
    dxdt_AST = RAST;    // amount of drug in gastric
    dxdt_AI = RAI; // amount of drug in intestinal
    dxdt_Acolon = Rcolon;    // amount of drug in colon
    dxdt_AAO = RAO; // amount of drug absorbed
    dxdt_Aplas_free = Rplas_free; //  amount of free drug in plasma
    dxdt_AL = RL;   // amount of drug in liver
    dxdt_Aurine = Rurine; // amount of drug in urine
    dxdt_AK = RK; // amount of drug in kidney
    dxdt_AM = RM; // amount of drug in muscle
    dxdt_AF = RF; // amount of drug in fat
    dxdt_ARest = RRest; // amount of drug in rest of the body
    dxdt_Amet  = Rmet;  // amount of drug in metabolism
      
     // Equation for the AUC of chemical in tissue compartment 
    dxdt_AUCCV = CV;                // AUC of chemical in venous blood
    dxdt_AUCCL = CL;                // AUC of chemical in liver compartment 
    dxdt_AUCCK = CK;                // AUC of chemical in kidney compartment  
    dxdt_AUCCM = CM;                // AUC of chemical in muscle compartment
    dxdt_AUCCF = CF;                // AUC of chemical in fat compartment
 
    double Qbal = QC - QL - QK - QM - QF - QRest; 
    double Tmass = Aplas_free + AL + AK + Aurine + AM + AF + ARest + Amet + Acolon;
    double bal   = ADOSE - Tmass;

$TABLE
capture Plasma = CV;
capture Liver  = CL;
capture Kidney = CK;
capture Muscle = CM;
capture Fat  = CF;
capture Rest = CRest;
capture Bal  = bal;
capture tmass = Tmass;
'


## Load Model
mod <- mcode_cache("pbpk", SolvePBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline

pred<- function (pars, tinterval, BW, TDOSE, route, Dose){
  
  DOSE = Dose * BW
  ## Get out of log domain
  pars <- exp(pars)   
  
  if (route == "oral gavage") {
    
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    
    ex <- ev_1+ev_2
    
  }
  
  if (route == "iv") {
    
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "Aplas_free", replicate = FALSE)
    
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    
    ex <- ev_1+ev_2
    
    
  }
  
  ## set up the exposure time
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*20, 0.01)
  ## Simulation
  out <- mod %>%param(pars)%>%
    
    mrgsim_d (data = ex, tgrid = tsamp)
  
  output <- cbind.data.frame(Time = out$time,
                             CP  = out$Plasma,
                             CL  = out$Liver,
                             CK  = out$Kidney,
                             CM  = out$Muscle,
                             CF  = out$Fat,
                             Crest = out$Rest,
                             AUC_V = out$AUCCV,
                             AUC_L = out$AUCCL,
                             AUC_K = out$AUCCK,
                             AUC_M = out$AUCCM,
                             AUC_F = out$AUCCF,
                             Bal = out$Bal)
  Final_time <<- tinterval*1
  return(output)
}

# Input the dataset for calibration and evaluation ------------
{
  ## Input the dataset for iv infusion 
  data1 <- read_excel("data1_10 mgkg (IV single)_plasma_tartrate.xls", col_types = c("numeric", "numeric", "numeric"))
  data1.1<- cbind.data.frame(Time = data1$Time, CP  = data1$CP) #tartrate BW: 2.3 – 2.6 kg
  
  data2 <- read_excel("data2_10 mgkg (IV single)_plasma_phosphate.xlsx", col_types = c("numeric", "numeric"))
  data2<-  cbind.data.frame(Time = data2$Time, CP  = data2$CP) #phosphate BW: 1.4 – 1.6 kg
  
  data3 <- read_excel("data2_10 mgkg (IV single)_plasma_tartrate.xls", col_types = c("numeric", "numeric"))
  data3<-  cbind.data.frame(Time = data3$Time, CP  = data3$CP) #tartrate  BW: 1.4 – 1.6 kg
  
  data5 <- read_excel("data6_50 mgkg (IV single)_plasma_tartrate.xls", col_types = c("numeric", "numeric"))
  data5<-  cbind.data.frame(Time = data5$Time, CP  = data5$CP) #tartrate  BW: 1.6 - 1.8 kg
  
  ## Input the dataset for oral gavage
  data6 <- read_excel("data2_10 mgkg (Oral single)_plasma_phosphate.xls", col_types = c("numeric", "numeric"))
  data6<-  cbind.data.frame(Time = data6$Time, CP  = data6$CP) #phosphate  BW: 1.6 - 1.8kg
  
  data7 <- read_excel("data2_10 mgkg (Oral single)_plasma_tartrate.xls", col_types = c("numeric", "numeric"))
  data7<-  cbind.data.frame(Time = data7$Time, CP  = data7$CP) #tartrate   BW: 1.4 – 1.6 kg
  
  data8 <- read_excel("data3_24 mgkg (Oral single)_plasma_tartrate.xls", col_types = c("numeric", "numeric"))
  data8<-  cbind.data.frame(Time = data8$Time, CP  = data8$CP) #tartrate   BW: NA
  
  data9 <- read_excel("data4_20 mgkg (Oral single)_plasma_phosphate.xls", col_types = c("numeric", "numeric", "numeric"))
  data9.1<- cbind.data.frame(Time = data9$Time, CP  = data9$CP) #phosphate  BW: NA
  
  data10 <- read_excel("data6_50 mgkg (Oral single)_plasma.xlsx", col_types = c("numeric", "numeric"))
  data10<-  cbind.data.frame(Time = data10$Time, CP  = data10$CP) #tartrate  BW: 1.6 - 1.8 kg
  
  data12 <- read_excel("data5_50 mgkg(PO 5days)_Kidney.xls", col_types = c("numeric", "numeric", "numeric"))
  data12<-  cbind.data.frame(Time = data12$Time, CK  = data12$Concentration) #tartrate  BW: 1.8 - 2.2 kg
  
  data13 <- read_excel("data5_50 mgkg(PO 5days)_Liver.xls", col_types = c("numeric", "numeric", "numeric"))
  data13<-  cbind.data.frame(Time = data13$Time, CL  = data13$Concentration)  #tartrate  BW: 1.8 - 2.2 kg
  
  data14 <- read_excel("data5_50 mgkg(PO 5days)_Muscle.xls", col_types = c("numeric", "numeric", "numeric"))
  data14<-  cbind.data.frame(Time = data14$Time, CM  = data14$Concentration)  #tartrate  BW: 1.8 - 2.2 kg
}

Cost <- function(pars,w) {
  out_A <-  pred (pars,  tinterval = 12, BW = 1.5, TDOSE = 1, route = 'iv', Dose =10)
  out_C <-  pred (pars,  tinterval = 12, BW = 1.5, TDOSE = 1, route = 'oral gavage', Dose =10)
  
  
  cost<- modCost(model = out_A, obs = data3, x ="Time", weight = w)
  cost<- modCost(model = out_C, obs = data7, x ="Time", weight = w, cost = cost)
  
  return(cost)
}

# Input parameters --------------------------------------------------------
pars <- c(
  Kst = 2,
  Kint = 0.108,
  KurineC = 0.1199, 
  KmC = 0.025,
  Ka = 0.23,
  PB = 0.324, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  PF       = 2.90,
  PRest    = 7.58
)

Cost(log(pars), w="mean")

Sns <- sensFun(func = Cost,w = "mean",
               parms = log(pars), varscale = 1)

Sen_1 <- summary(Sns)
Sen_1
plot(Sen_1)

# Select sensitive parameter ----------------------------------------------
theta <- pars[abs(Sen_1$Mean) > mean(abs(Sen_1$Mean))]
theta

Sen.rank_hum <- Sen_1 %>% mutate(rank.L1 = rank(abs(L1)),
                                 rank.L2 = rank(abs(L2)),
                                 rank.Mean = rank(abs(Mean)),
                                 rank.Min = rank(abs(Min)),
                                 rank.Max = rank(abs(Max)))
pars <- c(
  Kst = 2,
  Kint = 0.108,
  KurineC = 0.1199, 
  KmC = 0.025,
  Ka = 0.23
)

Fit <- modFit(f = Cost, p = log(pars),
              method = "Port",  w="mean",
              
              #lower = c(log(0.1199 * 0.1),log(0.025 * 0.1),log(0.23 * 0.1),log(2 * 0.1)),
              #upper = c(log(0.1199 * 10),log(0.025 * 10),log(0.23 * 10),log(2 * 10)),
              
              control = list(maxit = 10000, nprint = 1))

summary(Fit)                           ## Summary of fit
exp(Fit$par)                           ## Get the arithmetic value out of the log domain

# Model -------------------------------------------------------
#IV infusion
out_Aini <- pred (log(pars),  tinterval = 12, TDOSE = 1,  BW = 2.45, route = 'iv', Dose =10)
out_Afin <-  pred (Fit$par,  tinterval = 12, TDOSE = 1,  BW = 2.45, route = 'iv', Dose =10)
plot(out_Afin$Time, out_Afin$CP, type = 'l',xlim=c(0,15),
     xlab = "Time(h)", ylab = "Concentration(mg/L)",main = "Plasma, 10mg/kg IV single dose")
lines(out_Aini$Time, out_Aini$CP,lty = 2)
points(data1$Time, data1$CP) #Evaluation


out_A1ini <- pred (log(pars),  tinterval = 12, TDOSE = 1,  BW = 1.5, route = 'iv', Dose =10)
out_A1fin <-  pred (Fit$par,  tinterval = 12, TDOSE = 1,  BW = 1.5, route = 'iv', Dose =10)
plot(out_A1fin$Time, out_A1fin$CP, type = 'l',xlim=c(0,15),
     xlab = "Time(h)", ylab = "Concentration(mg/L)",main = "Plasma, 10mg/kg IV single dose")
lines(out_A1ini$Time, out_A1ini$CP,lty = 2)
points(data2$Time, data2$CP) #Evaluation #phosphate
points(data3$Time, data3$CP) #Calibration

out_Bfin <-  pred (Fit$par,  tinterval = 12, TDOSE = 1, BW = 2.0, route = 'iv', Dose =50)
plot(out_Bfin$Time, out_Bfin$CP, type = 'l',xlim=c(0,10),
     xlab = "Time(h)", ylab = "Concentration(mg/L)",main = "Plasma, 50mg/kg iv single dose")
points(data5$Time, data5$CP) #Evaluation

#(PO) --------------------------------------------------------
out_Cfin <-  pred (Fit$par,  tinterval = 12, TDOSE = 1, BW = 1.7, route = 'oral gavage', Dose =10)
out_Cini <-  pred (log(pars),  tinterval = 12, TDOSE = 1, BW = 1.7, route = 'oral gavage', Dose =10)
plot(out_Cfin$Time, out_Cfin$CP, type = 'l',xlim=c(0,10),
     xlab = "Time(h)", ylab = "Concentration(mg/L)",main = "Plasma, 10mg/kg PO single dose")
lines(out_Cini$Time, out_Cini$CP,lty = 2)
points(data6$Time, data6$CP) #Evaluation #phosphate
points(data7$Time, data7$CP) #Calibration


out_Dfin <-  pred (Fit$par,  tinterval = 12, TDOSE = 1, BW = 5.0, route = 'oral gavage', Dose =24)
plot(out_Dfin$Time, out_Dfin$CP, type = 'l',xlim=c(0,10),
     xlab = "Time(h)", ylab = "Concentration(mg/L)",main = "Plasma, 24mg/kg PO single dose")
points(data8$Time, data8$CP) #Evaluation


out_Efin <-  pred (Fit$par,  tinterval = 12, TDOSE = 1, BW = 2.5, route = 'oral gavage', Dose =20)
plot(out_Efin$Time, out_Efin$CP, type = 'l',xlim=c(0,10),
     xlab = "Time(h)", ylab = "Concentration(mg/L)",main = "Plasma, 20mg/kg PO single dose")
points(data9$Time, data9$CP) #Evaluation #phosphate


Cost <- function(pars,w) {
  
  output <- pred(pars, tinterval = 24, TDOSE = 5, BW = 2.0, route = 'oral gavage', Dose =50)
  
  out2 <- output%>%select(Time = Time, CL = CL)
  out3 <- output%>%select(Time = Time, CK = CK)
  out4 <- output%>%select(Time = Time, CM = CM)
  
  cost<- modCost(model = out2, obs = data13, x ="Time", weight = w)
  cost<- modCost(model = out3, obs = data12, x ="Time", weight = w, cost = cost)
  cost<- modCost(model = out4, obs = data14, x ="Time", weight = w, cost = cost)
  
  return(cost)
}

# Input parameters --------------------------------------------------------
pars <- c(
  Kst = 1.465123777  ,
  Kint = 1.075947416  ,
  KurineC = 1.381790038  , 
  KmC = 0.004342775 ,
  Ka = 0.389387102 ,
  PB = 0.324, 
  PL  = 10.60,
  PK  = 13.92,
  PM  = 2.90,
  PF       = 2.90,
  PRest    = 7.58
)


Cost(log(pars), w="mean")

Sns <- sensFun(func = Cost,w = "mean",
               parms = log(pars), varscale = 1)

Sen_1 <- summary(Sns)
Sen_1
plot(Sen_1)

# Select sensitive parameter ----------------------------------------------
theta <- pars[abs(Sen_1$Mean) > mean(abs(Sen_1$Mean))]
theta

Sen.rank_hum <- Sen_1 %>% mutate(rank.L1 = rank(abs(L1)),
                                 rank.L2 = rank(abs(L2)),
                                 rank.Mean = rank(abs(Mean)),
                                 rank.Min = rank(abs(Min)),
                                 rank.Max = rank(abs(Max)))

# Model fitting
pars <- c(
  PF       = 2.90,
  PRest    = 7.58
)

Fit <- modFit(f = Cost, p = log(pars),
              method = "Port",  w="mean",
              
              lower = c(log(2.90  * 0.1),log(7.58 * 0.1)),
              upper = c(log(2.90 * 10),log(7.58 * 10)),
              
              control = list(maxit = 10000, nprint = 1))

summary(Fit)                           ## Summary of fit
exp(Fit$par)                           ## Get the arithmetic value out of the log domain

# Oral/9016821_50 mg kg (Oral)  tissue concentration
# Liver
times <- c(24+120,48+120,72+120,96+120,120+120)
means <- c(42.0,26.0,13.0,4.0,1.0)
errors <- c(4.11,3.18,1.11,0.5,0.04)
out_Ffin <- pred (Fit$par,  tinterval = 24, TDOSE = 5,  BW = 2.0, route = 'oral gavage', Dose =50)
out_F <- pred (log(pars),  tinterval = 24, TDOSE = 5,  BW = 2.0, route = 'oral gavage', Dose =50)
plot(out_Ffin$Time, out_Ffin$CL, type = 'l', main = '50 mg/kg (Oral) 5 days',xlim=c(0,300),
     xlab = "Time(h)", ylab = "Liver Concentration(mg/L)")
lines(out_F$Time, out_F$CP,lty = 2)
points(times, means, pch=16, cex=1.5, col="red")  # adjust 'cex' for size, 'col' for color
arrows(times, means - errors, times, means + errors, angle=90, code=3, length=0.05)

# kidney (mg/L)
times <- c(24+120,48+120,72+120,96+120,120+120)
means <- c(55.0,33.0,16.0,7.0,2.0)
errors <- c(5.93,4.48,2.59,0.46,0.06)
out_Gfin <- pred (Fit$par,  tinterval = 24, TDOSE = 5,  BW = 2.0, route = 'oral gavage', Dose =50)
plot(out_Gfin$Time, out_Gfin$CK, type = 'l', main = '50 mg/kg (Oral) 5 days',
     xlab = "Time(h)", ylab = "Kidney Concentration(mg/L)")
points(times, means, pch=16, cex=1.5, col="red")  # adjust 'cex' for size, 'col' for color
arrows(times, means - errors, times, means + errors, angle=90, code=3, length=0.05)

# Muscle (mg/L)
times <- c(24+120,48+120,72+120,96+120)
means <- c(15.0,6.0,2.0,0.5)
errors <- c(1.13,1.25,0.45,0.03)
out_Hfin <- pred (Fit$par,  tinterval = 24, TDOSE = 5,  BW = 2.0, route = 'oral gavage', Dose =50)
plot(out_Hfin$Time, out_Hfin$CM, type = 'l', main = '50 mg/kg (Oral) 5 days',
     xlab = "Time(h)", ylab = "Muscle Concentration(mg/L)")
points(times, means, pch=16, cex=1.5, col="red")  # adjust 'cex' for size, 'col' for color
arrows(times, means - errors, times, means + errors, angle=90, code=3, length=0.05)

calibration_IV_10_1 <- out_A1fin[, c("Time","CP")]
write.csv(calibration_IV_10_1, "calibration_IV_10_1.csv")

calibration_PO_10 <- out_Cfin[, c("Time","CP")]
write.csv(calibration_PO_10, "calibration_PO_10.csv")

calibration_liver <- out_Ffin[, c("Time","CL")]
write.csv(calibration_liver, "calibration_liver.csv")

calibration_kidney <- out_Gfin[, c("Time","CK")]
write.csv(calibration_kidney, "calibration_kidney.csv")

calibration_muscle <- out_Hfin[, c("Time","CM")]
write.csv(calibration_muscle, "calibration_muscle.csv")

evaluation_IV_10_1 <- out_A1fin[, c("Time","CP")]
write.csv(evaluation_IV_10_1, "evaluation_IV_10_1.csv")

evaluation_IV_10 <- out_Afin[, c("Time","CP")]
write.csv(evaluation_IV_10, "evaluation_IV_10.csv")

evaluation_IV_50 <- out_Bfin[, c("Time","CP")]
write.csv(evaluation_IV_50, "evaluation_IV_50.csv")

evaluation_PO_10 <- out_Cfin[, c("Time","CP")]
write.csv(evaluation_PO_10, "evaluation_PO_10.csv")

evaluation_PO_24 <- out_Dfin[, c("Time","CP")]
write.csv(evaluation_PO_24, "evaluation_PO_24.csv")

evaluation_PO_20 <- out_Efin[, c("Time","CP")]
write.csv(evaluation_PO_20, "evaluation_PO_20.csv")

