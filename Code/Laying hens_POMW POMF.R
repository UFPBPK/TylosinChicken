## Project title: Development of a physiologically based pharmacokinetic model for tylosin in laying hens following oral exposure
## Objective: To create a validated PBPK model to accurately estimate chemical residues and extralabel withdrawal intervals of tylosin in eggs from laying hens.
## Tylosin, a hydrosoluble compound, with a molecular weight of 916 g/mol, has a low oral bioavailability and a low apparent volume of distribution (https://doi.org/10.4315/0362-028X.JFP-13-440)

## Four administration methods: iv; oral gavage; POMF; POMW

# Clear current environment -----------------------------------------------
rm(list=ls(all=TRUE))

# Set routine -------------------------------------------------------------
getwd()
setwd("C:/Users/zhichengzhang/Desktop/Tyloasin_08_18/1. Pharmacokinetics data/Laying hens")

# Set packages ------------------------------------------------------------
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


#Build mrgsolve-based PBPK Model --------------------------------------
SolvePBPK <- '
$PARAM @annotated

// Oral absorption rate constants
Kst    : 1.465123777:               /h, gastric emptying rate constant
Kint   : 1.075947416     :          /h, intestinal transit rate constant

// IV infusion rate constants
Timeiv  : 0.01     :                       /h, IV injection/infusion time

// Urinary elimination rate constant converted from elimination half life (5.78h)
KurineC : 1.381790038     :                L/h/kg, From broiler chicken model

// Metabolic rate constant adjusted by body weight
KmC     : 0.004342775      :                     L/h/kg, From broiler chicken model

// Dosing, multiple oral gavage
Ka      : 0.38938710          :               /h, intestinal absorption rate constant


// Physiological parameters
// Blood flow rates
QCC     : 10.17      :                Cardiac output (L/h/kg), Table 13, Wang et al. 2021
QLC     : 0.2526    :                 Fraction of blood flow to the liver, Table 16, Wang et al. 2021
QKC     : 0.2012    :                 Fraction of blood flow to the kidneys, Table 16, Wang et al. 2021
QFC     : 0.015     :                 Fraction of blood flow to the fat, Zeng et al. 2019
QMC     : 0.20      :                 Fraction of blood flow to the muscle, Table 16, Wang et al. 2021
QovaryC : 0.1646    :                 Fraction of blood flow to the ovary, Table 20, Wang et al. 2021
QoviductC : 0.1381  :                 Fraction of blood flow to the oviduct (infundlbulum/magnum/isthmus/uterus/Shell Gland/vagina), Table 20, Wang et al. 2021

// Tissue  volumes
BW      : 2.2       :                 Body weight (kg)
VLC     : 0.0249    :                 Fractional liver weight, Table 3, Wang et al. 2021
VKC     : 0.0076    :                 Fractional kidney weight, Table 3, Wang et al. 2021
VFC     : 0.05      :                 Fractional fat weight, Zeng et al. 2019
VMC     : 0.40      :                 Fractional muscle weight, Table 2, Wang et al. 2021
VovaryC : 0.0191    :                 Fractional ovary weight, Table 3, Wang et al. 2021
VbloodC : 0.0483    :                 Fractional blood volume, Table 22, Wang et al. 2021
Htc     : 0.307     :                 Hematocrit, Table 22, Wang et al. 2021
VoviductC : 0.0258  :                 Fractional the oviduct (infundlbulum/magnum/isthmus/uterus/vagina), Table 16, Wang et al. 2021

// Mass Transfer Parameters (Chemical-specific parameters)
// Partition coefficients (PC, tissue:plasma)
PB      : 0.324       :                Percentage of drug bound to plasma proteins
PL      : 10.60     :                Liver plasma PC, from broiler chicken 
PK      : 13.92     :                Kidney plasma PC, from broiler chicken 
PM      : 2.90      :                Muscle plasma PC, breast muscle, calculated using AUCtissue:AUCplasma method 
PF      : 2.7317     :                Fat plasma PC, assumed to be the same as muscle

Povary  : 0.19190526       :           
Poviduct: 2.34127347       :      
PRest   : 6.1966           :                Adapted from broiler chicken

// Egg yolk and white parameters
Ay   : 25                             : apparent maximum follicle weight g
Ky   : 0.05560858                          : transport constant into yolk
Kw   : 0.10268711                          : transport constant into albumen
tlag : 1*24                           : h
tsig : 2*24                           : h
s    : 1/24                           : /h
talbumen : 10                         : h
tlay1  : 12                           : tlay1 default
Rwhitefor : (34*0.001)/10             : rate of albumen formation Kg/h 

$MAIN

// Cardiac output and blood flows to tissues/organs (L/h)
double QC = QCC * BW;                  // Cardiac output
double QL = QLC * QC;                  // Liver
double QK = QKC * QC;                  // Kidney
double QF = QFC * QC;                  // Fat
double QM = QMC * QC;                  // Muscle
double Qovary = QovaryC * QC;          // Ovary
double Qoviduct = QoviductC * QC;      // Oviduct
double QRest = QC - QK - QL - QF - QM - Qovary - Qoviduct; // Rest of the tissues


//  Tissue volumes (L)
double VL = VLC * BW;                  // Liver
double VK = VKC * BW;                  // Kidney
double VF = VFC * BW;                  // Fat
double VM = VMC * BW;                  // Muscle
double Vblood = VbloodC*BW;            // Blood
double Vplasma = VplasmaC * BW;        // Plasma
double Vovary = VovaryC * BW;          // Ovary
double Voviduct = VoviductC * BW;      // Oviduct
double VRest = 1 * BW - VL - VK - VF - VM - Vovary - Vblood - Voviduct; // Rest of the tissues
double VplasmaC = (1-Htc) * VbloodC;   // Plasma volume, fraction of BW, Table 22, Wang et al. 2021
double Kmet = KmC * BW;	               // L/h Metabolic rate constant adjusted by body weight
double Kurine = KurineC * BW;          // L/h Urinary elimination rate
double Free = 1-PB;                    // Free drug percentage


// Egg yolk and white functions
double tlay   =  tlay1*24;             // h
double Wy     = (25/(1 + exp(-s*(TIME - (tlay - tsig - tlag)))))*0.001;  //  Kg weight of the follicle
double Ww = (TIME >= tlay-tlag && TIME <= tlay - tlag + talbumen) ?  // Kg weight of the egg white
       Rwhitefor * (TIME - tlay + tlag) : 
       (TIME > tlay - tlag + talbumen) ? 
       Rwhitefor * talbumen : 
       0;

// Only transport drug into yolk during the yolk formation period
double Ryolk = (TIME >= 0 && TIME <= tlay - tlag) ? Covary * Ky * Wy : 0; // rate of drug deposition into yolk mg/h     

// Only transport drug into white during the white formation period
double Rwhite = (TIME >= tlay-tlag && TIME <= tlay - tlag + talbumen) ? Coviduct * Kw * Ww : 0;   // rate of drug deposition into white mg/h     

double Cyolk   = Ayolk/Wy;    // Tylosin concentration in Yolk
double Cwhite  = Awhite/Ww;   // Tylosin concentration in white
double Cegg    = (Ayolk + Awhite) / (Ww + Wy);   // Tylosin concentration in egg


$CMT ADOSE AST AI Acolon AAO Aplas_free AL AK Aurine Amet AM Aovary Ayolk Awhite AWw AWy AUCCV AUCCL AUCCK AUCCM AUCCF AUCCovary AUCCoviduct AF ARest Aoviduct Aegg_ex
             
$ODE

dxdt_ADOSE=0;

 // The changing rate of the amount of dose via oral

 // Concentration of the chemical in vein and tissue compartments
      double CVL = AL/(VL * PL);      // Concentration in liver / PC of liver
      double CVK = AK/(VK * PK);      // Concentration in Kidney / PC of kidney
      double CVF = AF/(VF * PF);   // Concentration in fat / PC of fat
      double CVRest = ARest/(VRest * PRest); // Concentration in Rest  / PC of rest of tissues
      double CVM = AM/(VM * PM);      // Concentration in muscle / PC of muscle 
      double CVovary = Aovary/(Vovary * Povary);  // Concentration in ovary/PC of ovary
      double CVoviduct = Aoviduct/(Voviduct * Poviduct);  // Concentration in oviduct vein
      double CV = ((QL * CVL  + QK  * CVK  + QF * CVF  + QM * CVM  + Qovary * CVovary + QRest * CVRest + Qoviduct * CVoviduct)/QC);
      double Cplas_free = Aplas_free/Vplasma;  // Free drug concentration in plasma = amount of free drug in plasma / volume of plasma
      double Cplasma = Cplas_free/Free ; // drug concentration in plasma 
      double CL = AL/VL;              // Concentration in liver
      double CK = AK/VK;              // Concentration in kidney
      double CM = AM/VM;              // Concentration in muscle
      double CF = AF/VF;              // Concentration in fat 
      double Covary = Aovary/Vovary;  // Concentration in ovary 
      double Coviduct = Aoviduct/Voviduct;  // Concentration in oviduct
      double CRest = ARest/VRest;     // Concentration in rest of body

  // Rate of concentration changes in different compartments    
    double RAST = - Kst * AST;        // oral to the stomach 
    double RAI = Kst * AST - Kint * AI - Ka * AI; // in Intestine 
    double Rcolon = Kint * AI ;       // in Colon 
    double RAO = Ka * AI;
    double Rplas_free = QC * (CV - Cplasma) * Free; // rate of change in amount of chem in plasma compartment of tissues
    double RL = QL * (Cplasma - CVL)* Free + RAO - Kmet * AL ;// rate of change in amount of the drug in liver
    double Rmet = Kmet * AL;          // Rmet the metabolic rate in liver (mg/h)
    double Rurine = Kurine * CVK;   // Rurine the elimination rate in liver (mg/h)
    double RK = QK * (Cplasma - CVK) * Free - Rurine; // rate of change in amount of the drug in kidney
    double RM = QM * (Cplasma - CVM) * Free;  // rate of change in amount of the drug in muscle
    double RF = QF * (Cplasma - CVF) * Free;  // rate of change in amount of the drug in fat
    double RRest = QRest * (Cplasma - CVRest) * Free; // rate of change in amount of the drug in rest of the body
    double Rovary = Qovary * (Cplasma - CVovary) * Free - Ryolk; // rate of change in amount of the drug in ovary 
    double Roviduct = Qoviduct * (Cplasma - CVoviduct) * Free - Rwhite; // rate of change in amount of the drug in oviduct 
    double Regg_ex = Ryolk + Rwhite; // mg/h


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
    dxdt_Aovary = Rovary; // amount of drug in ovary
    dxdt_ARest = RRest; // amount of drug in rest of the body
    dxdt_Amet  = Rmet;  // amount of drug in metabolism
    dxdt_Ayolk = Ryolk; // amount of drug in yolk
    dxdt_AWw = Ww; // the weight of white
    dxdt_AWy = Wy; // the weight of yolk
    dxdt_Awhite = Rwhite; // amount of drug in white
    dxdt_Aovary = Rovary; // amount of drug in ovary
    dxdt_Aoviduct = Roviduct; // amount of drug in oviduct
    dxdt_Aegg_ex = Regg_ex;

      
     // Equation for the AUC of chemical in tissue compartment 
    dxdt_AUCCV = CV;                // AUC of chemical in venous blood
    dxdt_AUCCL = CL;                // AUC of chemical in liver compartment 
    dxdt_AUCCK = CK;                // AUC of chemical in kidney compartment  
    dxdt_AUCCM = CM;                // AUC of chemical in muscle compartment
    dxdt_AUCCF = CF;                // AUC of chemical in fat compartment
    dxdt_AUCCovary = Covary;        // AUC of chemical in ovary compartment
    dxdt_AUCCoviduct = Coviduct;    // AUC of chemical in oviduct compartment
      
    double Qbal = QC - QL - QK - QM - QF - Qovary - Qoviduct - QRest; 
    double Tmass = Aplas_free + AL + AK + Aurine + AM + AF + Aovary + Aoviduct + ARest + Amet + Acolon + AI + AST;
    double bal   = ADOSE - Tmass;

$TABLE
capture massE = Ww;
capture massF = Wy;
capture massW = Awhite;
capture massY = Ayolk;
capture Plasma = CV;
capture Liver  = CL;
capture Kidney = CK;
capture Muscle = CM;
capture Fat  = CF;
capture Ovary  = Covary;
capture Oviduct  = Coviduct;
capture Rest = CRest;
capture Bal  = bal;
capture tmass = Tmass;
capture yolk = Cyolk;
capture white = Cwhite;
capture BloodBalance = Qbal;
capture egg = Cegg;
'


# Load Model --------------------------------------------------------------
{
  mod <- mcode_cache("pbpk", SolvePBPK)
  
  pred<- function (pars, BW, TDOSE, route, Dose, tlay1 ){
    
    DOSE = Dose * BW  / 15  # Calculate hourly dose during the 15 hours of daytime
    
    pars <- exp(pars)   
    
    evs <- list()
    for (day in 0:(TDOSE - 1)) {
      for (hour in 7:22) {  # From 7 AM to 11 PM
        time_of_dose = 24 * day + hour
        if (route == "oral gavage") {
          evs[[length(evs) + 1]] <- ev(ID = 1, amt = DOSE, time = time_of_dose, cmt = "AST", replicate = FALSE)
        } else if (route == "iv") {
          evs[[length(evs) + 1]] <- ev(ID = 1, amt = DOSE, time = time_of_dose, tinf = 0.01, cmt = "Aplas_free", replicate = FALSE)
        }
      }
    }
    
    ex <- Reduce(`+`, evs)
    
    
    # Set up the exposure time ------------------------------------------------
    
    tsamp  = tgrid(0, 24 * TDOSE + 240, 0.1)
    
    # Simulation --------------------------------------------------------------
    
    out <- mod %>% param(c(tlay1=tlay1))%>%param(pars)%>%
      
      mrgsim_d (data = ex, tgrid = tsamp)
    
    output <- cbind.data.frame(Time = out$time,
                               Wy = out$massF,
                               Ww = out$massE,
                               Awhite = out$massW,
                               Ayolk = out$massY,
                               CP  = out$Plasma,
                               CL  = out$Liver,
                               CK  = out$Kidney,
                               CM  = out$Muscle,
                               CF  = out$Fat,
                               Crest = out$Rest,
                               Cyolk = out$yolk, 
                               Cwhite = out$white, 
                               Covary = out$Ovary,
                               Coviduct = out$Oviduct,
                               AUC_V = out$AUCCV,
                               AUC_L = out$AUCCL,
                               AUC_K = out$AUCCK,
                               AUC_M = out$AUCCM,
                               AUC_F = out$AUCCF,
                               AUC_Ovary = out$AUCCovary,
                               AUC_Oviduct = out$AUCCoviduct,
                               Bal = out$Bal,
                               BloodBalance = out$BloodBalance,
                               Cegg = out$egg)
    return(output)
    
  }
}

# Administration route: POMW ----------------------------------------------
# Input the dataset for POMW ----------------------------------------------
{
  Water_100_5_Yolk <- read_excel("data 13_100 mgkg (POMW 5days)_Yolk.xlsx", col_types = c("numeric", "numeric", "numeric", "numeric"))
  data28 <- cbind.data.frame(Time = Water_100_5_Yolk$Time, Cyolk  = Water_100_5_Yolk$Yolk) #BW:1.7 – 2.4 kg
  
  Water_100_5_Albumin <- read_excel("data 13_100 mgkg (POMW 5days)_Albumi.xlsx", col_types = c("numeric", "numeric", "numeric", "numeric"))
  data29 <- cbind.data.frame(Time = Water_100_5_Albumin$Time, Cwhite  = Water_100_5_Albumin$Albumin) #BW:1.7 – 2.4 kg
  
  Water_100_5_Egg <- read_excel("data 13_100 mgkg (POMW 5days)_Egg.xlsx", col_types = c("numeric", "numeric", "numeric", "numeric"))
  data30 <- cbind.data.frame(Time = Water_100_5_Egg$Time, Cegg  = Water_100_5_Egg$Egg) #BW:1.7 – 2.4 kg
  
  Water_77.25_5_Egg_W4 <- read_excel("data12_77.25 mgkg (POMW 5days)_Egg_W4.xlsx",col_types = c("numeric", "numeric"))
  data32 <- cbind.data.frame(Time = Water_77.25_5_Egg_W4$Time, Cegg  = Water_77.25_5_Egg_W4$Concentration) #BW:1.7 – 2.4 kg
  
  Water_77.25_5_Egg_W5 <- read_excel("data12_77.25 mgkg (POMW 5days)_Egg_W5.xlsx",col_types = c("numeric", "numeric"))
  data33 <- cbind.data.frame(Time = Water_77.25_5_Egg_W5$Time, Cegg  = Water_77.25_5_Egg_W5$Concentration) #BW:1.7 – 2.4 kg
}

Cost <- function(pars,w) {
  t1  <- data28$Time
  tout1 <- list()
  for (i in 1:length(t1)) {
    output1 <- pred(pars, TDOSE = 5, route = 'oral gavage', Dose =100, BW = 2.05, tlay1 = t1[i]/24)
    tout1[[i]] <-output1 %>% filter(Time == t1[i])
  }
  
  yout1 <- do.call(rbind,tout1)%>%select(Time = Time, Cwhite = Cwhite)
  yout2 <- do.call(rbind,tout1)%>%select(Time = Time, Cyolk = Cyolk)
  yout3 <- do.call(rbind,tout1)%>%select(Time = Time, Cegg = Cegg)
  
  cost<- modCost(model = yout1, obs = data29, x ="Time", weight = w)
  cost<- modCost(model = yout2, obs = data28, x ="Time", weight = w, cost = cost)
  cost<- modCost(model = yout3, obs = data30, x ="Time", weight = w, cost = cost)
  
  return(cost)
}


pars <- c(
  kst = 1.465123777,
  Ka = 0.38938710,
  Kint = 1.075947416  ,
  Povary = 1.08,
  Poviduct = 1.08,
  KurineC = 1.381790038 ,
  KmC = 0.004342775 ,
  
  Kw = 0.0640       ,
  Ky = 0.0682     
)

Fit<- modFit(f=Cost, p= log(pars), method ="Nelder-Mead",  w="mean",
             #lower = c( log(0.00108), log(0.023), log(0.108), log(0.108), log(0.00640), log(0.00682)),
             #upper = c( log(1.8),    log(2.3),    log(10.8), log(10.8), log(6.40), log(6.82)),
             
             control = nls.lm.control(nprint=1))

summary(Fit)                                 ## Summary of fit 
exp(Fit$par)                                 ## Get the arithmetic value out of the log domain

Cost(Fit$par, w="mean")

# POMW - Calibration --------------------------------------------------------
#100 mg/kg (water) 5days Egg
t <- c((1:20)*24)
tout1 <- list()
for (i in 1:length(t)) {
  out_G <- pred(Fit$par,  TDOSE = 5, route = 'oral gavage', Dose =100, BW = 2.05, tlay1 = t[i]/24)
  tout1[[i]] <-out_G%>% filter(Time == t[i])
}

yout3 <- do.call(rbind,tout1)%>%select(Time = Time, Cwhite = Cwhite, Cyolk = Cyolk, Cegg = Cegg, Coviduct = Coviduct, Covary = Covary)
plot(yout3$Time, yout3$Cegg, xlab = "time", ylab = "Cegg" , type = 'l')
points(data30$Time,data30$Cegg)

plot(yout3$Time, yout3$Cyolk, xlab = "time", ylab = "Cyolk" , type = 'l')
points(data28$Time,data28$Cyolk)

plot(yout3$Time, yout3$Cwhite, xlab = "time", ylab = "Cwhite" , type = 'l')
points(data29$Time,data29$Cwhite)

# POMW - Evaluation --------------------------------------------------------
#77.25 mgkg (POMW 5days) Egg
t <- c((1:20)*24)
tout1 <- list()
for (i in 1:length(t)) {
  out_G <- pred(Fit$par, TDOSE = 5, route = 'oral gavage', Dose =77.25, BW = 2, tlay1 = t[i]/24)
  tout1[[i]] <-out_G%>% filter(Time == t[i])
}

yout4 <- do.call(rbind,tout1)%>%select(Time = Time, Cwhite = Cwhite, Cyolk = Cyolk, Cegg = Cegg, Coviduct = Coviduct, Covary = Covary)
plot(yout4$Time, yout4$Cegg, xlab = "time", ylab = "Cegg" , type = 'l')

points(data32$Time,data32$Cegg)
points(data33$Time,data33$Cegg)



# Input the dataset for POMF ----------------------------------------------
{
  Feed_30.625_Egg <- read_excel("data8_30.625 mgkg (POMF 14days)_Egg.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data21 <- cbind.data.frame(Time = Feed_30.625_Egg$Time-17, Cegg  = Feed_30.625_Egg$Egg, SD = Feed_30.625_Egg$SD)
  data21.1 <- data21 %>% select(Time, Cegg) #BW : 1.6
  
  Feed_96_Egg <- read_excel("data12_96 mgkg (POMF 5days)_Egg_F5.xlsx",col_types = c("numeric", "numeric"))
  data22 <- cbind.data.frame(Time = Feed_96_Egg$Time-17, Cegg  = Feed_96_Egg$Concentration) #BW : 1.7-2.4
  
  Feed_96_Egg1 <- read_excel("data12_96 mgkg (POMF 5days)_Egg_F9.xlsx",col_types = c("numeric", "numeric"))
  data23 <- cbind.data.frame(Time = Feed_96_Egg1$Time-17, Cegg  = Feed_96_Egg1$Concentration) #BW : 1.7-2.4
}
# Feed diet/9010766_30.625 mg kg (Feed) 14days_Egg
t <- c((1:40)*24)
tout1 <- list()
for (i in 1:length(t)) {
  output1 <- pred(Fit$par, TDOSE = 14, route = 'oral gavage', BW = 1.6, Dose =30.625, tlay1 = t[i]/24)
  tout1[[i]] <-output1 %>% filter(Time == t[i])
}
yout1 <- do.call(rbind,tout1)%>%select(Time = Time, Cyolk = Cyolk, Cwhite = Cwhite, Cegg = Cegg)
plot(yout1$Time, yout1$Cegg, xlab = "Time", ylab = "Cegg" , type = 'l')
points(data21$Time, data21$Cegg)

# 9011717_96 mg kg (Feed) 5days_Egg (270mg) phosphate
t <- c((1:40)*24)
tout1 <- list()
for (i in 1:length(t)) {
  output1 <- pred(Fit$par, TDOSE = 5, route = 'oral gavage', BW = 2.05, Dose =96, tlay1 = t[i]/24)
  tout1[[i]] <-output1 %>% filter(Time == t[i])
}
yout2 <- do.call(rbind,tout1)%>%select(Time = Time, Cyolk = Cyolk, Cwhite = Cwhite, Cegg = Cegg)
plot(yout2$Time, yout2$Cegg, xlab = "Time", ylab = "Cegg" , type = 'l')
points(data22$Time, data22$Cegg)
points(data23$Time, data23$Cegg)


evaluation_laying_POMF_30.625 <- yout1[, c("Time","Cegg")]
write.csv(evaluation_laying_POMF_30.625, "evaluation_laying_POMF_30.625.csv")

evaluation_laying_POMF_96 <- yout2[, c("Time","Cegg")]
write.csv(evaluation_laying_POMF_96, "evaluation_laying_POMF_96.csv")



calibration_laying_POMW <- yout3[, c("Time","Cwhite","Cyolk","Cegg")]
write.csv(calibration_laying_POMW, "calibration_laying_POMW.csv")

evaluation_laying_POMW <- out_G[, c("Time","Cegg")]
write.csv(evaluation_laying_POMW, "evaluation_laying_POMW_77.25.csv")

