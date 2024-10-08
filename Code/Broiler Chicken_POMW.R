## Species: Broiler chicken
## Tylosin, a hydrosoluble compound, with a molecular weight of 916 g/mol, has a low oral bioavailability and a low apparent volume of distribution (https://doi.org/10.4315/0362-028X.JFP-13-440)

# Clear current environment -----------------------------------------------
rm(list=ls(all=TRUE))

getwd()
setwd("C:/Users/zhichengzhang/Desktop/Tyloasin_08_18/1. Pharmacokinetics data/Broiler chicken")

library(readxl)

{
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
Kst    : 1.465123777        :                   /h, gastric emptying rate constant
Kint   : 1.075947416        :                   /h, intestinal transit rate constant

// IV infusion rate constants
Timeiv  : 0.01     :                  /h, IV injection/infusion time

// Urinary elimination rate constant
KurineC : 1.381790038        :         L/h/kg

// Metabolic rate constant adjusted by body weight
KmC     : 0.004342775        :         L/h/kg, metabolic rate constant

// Dosing, multiple oral gavage
Ka      : 0.38938710        :         /h, intestinal absorption rate constant

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
PF      : 2.7317    :                Fat plasma PC, assumed to be the same as muscle
PRest   : 6.1966    :                Rest of the body tissues:plasma PC

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

pred<- function (pars, BW, TDOSE, route, Dose){
  
  DOSE = Dose * BW/15 # Calculate the hourly dose during the 15 hours of daytime
  ## Get out of log domain
  pars <- exp(pars)   
  ## Create events for dosing from 07:00 to 24:00, restarting each day at 07:00
  evs <- list()  # Initialize list to hold dosing events
  total_hours = TDOSE * 24  # Total hours to simulate based on the number of days
  # Create dosing events for each hour between 07:00 and 24:00 each day
  for (day in 0:(TDOSE - 1)) {
    for (hour in 7:22) {  # From 7 AM to 11 PM
      time_of_dose = 24 * day + hour
  
  if (route == "oral gavage") {
    
    evs[[length(evs) + 1]] <- ev(ID = 1, amt = DOSE, time = time_of_dose, cmt = "AST")
    
    }
  
  else if  (route == "iv") {
    
    evs[[length(evs) + 1]] <- ev(ID = 1, amt = DOSE, time = time_of_dose, tinf = 0.01, cmt = "Aplas_free")
    
      }
    }
  }
  
  # Combine all events into a single event object
  ex <- Reduce(`+`, evs)
  ## set up the exposure time
  tsamp  = tgrid(0, total_hours, by = 1)
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
  return(output)
}


# Input the dataset for calibration and evaluation (POMW) -----------------
{
  data18 <- read_excel("data7_51 mgkg (POMW 5days)_Plasma.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data18.1 <- cbind.data.frame(Time = data18$Time, CP  = data18$CP) #tartrate BW: 754 +- 4.6 g
  
  data19 <- read_excel("data7_102 mgkg (POMW 5days)_Plasma.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data19.1 <- cbind.data.frame(Time = data19$Time, CP  = data19$CP) #tartrate BW: 754 +- 4.6 g
}

# Create the Cost function for POMW ---------------------------------------
Cost <- function(pars,w) {
  out_D <-  pred (pars, TDOSE = 5, BW= 0.75, route = 'oral gavage', Dose =102)
  
  # Calculate costs
  cost <- modCost(model = out_D, obs = data19.1, x = "Time")
  
  return(cost)
}

pars <- c(
  Kint = 1.075947416 ,
  Ka = 0.38938710)

Fit <- modFit(f = Cost, p = log(pars),
              method = "Port",  w="mean",
              
              lower = c(log(1.075947416 * 0.5), log(0.38938710 *0.5)),
              upper = c(log(1.075947416 * 2), log(0.38938710 *2)),

              control = list(maxit = 10000, nprint = 1))
summary(Fit)                           ## Summary of fit
exp(Fit$par)                           ## Get the arithmetic value out of the log domain

# POMW - Calibration -------------------------------------------------------
#51 mg/kg (water 200mgL) 5 days Plasma
{
  out_Lfin <- pred (Fit$par, TDOSE = 5, BW= 0.75, route = 'oral gavage', Dose =51)
  
  p13 <- ggplot(data=data18, aes(x=Time, y=CP), scale = "Free")+
    geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
    geom_errorbar(aes(ymin=CP-SD, ymax=CP+SD), colour="black") + # add this line for error bars
    geom_line(data=out_Lfin, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("51 mg/kg oral 5 days (POMW)") +
    ylab("Plasma Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold"), axis.text    
                                                                                 = element_text(size = 15,  color = "black", face="bold"),
                                                                                 axis.title   = element_text(size = 18, face = "bold", color = "black"),
                                                                                 strip.text.x = element_text(size=1, face="bold", color="black"),
                                                                                 #strip.text.y = element_text(size=20, face="bold", color="black"),
                                                                                 legend.position = "none")+ 
    scale_x_continuous(limits=(c(60,120)))+
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line())
  p13
}

# POMW - Evaluation -------------------------------------------------------
#102 mg/kg (water 200mgL) 5 days Plasma
{
  out_Mfin <- pred (Fit$par, TDOSE = 5, BW= 0.75, route = 'oral gavage', Dose =102)
  
  p14 <- ggplot(data=data19, aes(x=Time, y=CP), scale = "Free")+
    geom_point(shape = 21, colour = "black", fill = "red", size = 1, stroke =2)+
    geom_errorbar(aes(ymin=CP-SD, ymax=CP+SD), colour="black") + # add this line for error bars
    geom_line(data=out_Mfin, aes(x=Time, y=CP), size=0.8, colour = "black")+ ggtitle("102 mg/kg oral 5 days (POMW)") +
    ylab("Plasma Concentration (mg/L)")+ xlab("Time (hours)") +theme_bw() +theme(plot.title = element_text(size=15, face="bold"), axis.text    
                                                                                 = element_text(size = 15,  color = "black", face="bold"),
                                                                                 axis.title   = element_text(size = 18, face = "bold", color = "black"),
                                                                                 strip.text.x = element_text(size=1, face="bold", color="black"),
                                                                                 #strip.text.y = element_text(size=20, face="bold", color="black"),
                                                                                 legend.position = "none")+ 
    scale_x_continuous(limits=(c(60,120)))+
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line())
  p14
}

calibration_POMW_51 <- out_Lfin[, c("Time","CP")]
write.csv(calibration_POMW_51, "calibration_POMW_51.csv")

evaluation_POMW_102 <- out_Mfin[, c("Time","CP")]
write.csv(evaluation_POMW_102, "evaluation_POMW_102.csv")
