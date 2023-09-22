##Intrahepatic GCDCA concentrations
## Coupled drug and bile acid PBK model
## Generic drug PBK model is based on Punt, Ans, et al. "Predictive performance of next generation human physiologically based kinetic 
## (PBK) models based on in vitro and in silico input data." Altex 39.2 (2022): 221-234.
## Bile acid PBK model is based on De Bruijn, Véronique MP, et al. "Intestinal in vitro transport assay combined with physiologically 
## based kinetic modeling as a tool to predict bile acid levels in vivo." ALTEX-Alternatives to animal experimentation (2023).


## load packages
library(tidyverse)
library(cowplot)

## set working directory
setwd()
source("./functions/20230713_cmp_liver_cBA_work.R")

## sens scaling factor (1=default)
sens <- 1

####### simulations for maximal prescribed daily dose
model_results <- read.csv("./input/inputDataFile_final.csv",  fileEncoding="UTF-8-BOM") %>% 
  filter(reference=="PDD.max") %>% #filter for maximal prescribed daily dose
  filter(method.BSEP=="SHH") %>% #use suspension human hepatocytes for inhibitory potency of hepatic bile acid efflux inhibition
  group_by(compound, severity, Cmax.plasma, dose,reference, 
           method.KaFa,method.phys.chem.Fup,  method.Fup, method.phys.chem, 
           method.partition, method.Clint, scaling.for.liver, 
           MW, logP, pKa1, pKa2, frac.unionized, frac.ionized.acid, frac.ionized.base, Fup, Fuinc, Ka, Fa,
           Kpad, Kpbo, Kpbr, Kpgu, 
           Kphe, Kpki, Kpli,Kplu,Kpmu,Kpsk,Kpsp, 
           Clint, method.reference.Clint, 
           original.reference.Clint, BP, Ki, method.BSEP) %>%
  drop_na(Ki) %>% 
  do(model(pars = c(
    #physiological parameters
    BW = 70,
    dose =.$dose*1e-3/.$MW*1e6*70*.$Fa,
    FVad = 0.21, FVbo = 0.085629, FVbr = 0.02, FVgu = 0.0171,	
    FVhe = 0.0047, FVki = 0.0044, FVli = 0.021, FVliw=0.017535, FVlew=0.003465, FVlu = 0.0076,		
    FVmu = 0.4, FVsk = 0.0371, FVsp = 0.0026, FVte = 0.01,                	 
    FVve = 0.0514, FVpv=0.0064, FVar = 0.0257, 
    FVre =  0.0997711, 
    Vil=5.9, #volume ileum, L
    QC = 389.988, 
    FQad = 0.05, FQbo = 0.05, FQbr = 0.12, FQgu = 0.146462,     
    FQhe = 0.04, FQki = 0.19, FQh = 0.215385, FQlu = 1,	         	
    FQmu = 0.17, FQsk = 0.05, FQsp = 0.017231, FQte = 0.01076, 
    FQre = 0.103855,        
    
    #Compound specific parameters
    Kpad = .$Kpad, Kpbo = .$Kpbo, Kpbr = .$Kpbr, Kpgu = .$Kpgu,
    Kphe = .$Kphe,  Kpki = .$Kpki,  Kpli = .$Kpli,	Kplu = .$Kplu,	
    Kpmu = .$Kpmu, Kpsk = .$Kpsk, Kpsp = .$Kpsp,  Kpte = .$Kpmu, 
    Kpre = .$Kpmu, 
    fup= .$Fup, 
    fuhep = .$Fuinc,
    CLint = .$Clint,
    scaling.for.liver=.$scaling.for.liver,
    CLrenal = 6.7,
    SF= .$scaling.factor,
    MW = .$MW, 
    BP = .$BP,
    Ka = .$Ka,
    Fuinh=.$Fuinh,
    Ki=.$Ki, 
    
    #BA specific parameters
    Gdose=3020,
    Kpad_cBA=0.05, Kpgu_cBA=0.18, 
    Kpslp_cBA=0.19, Kpra_cBA=0.136,  	
    Kplew_cBA=0.31,  BP_BA=0.55,#  
    QGb=0.5, 
    ka_cBA=0.09, #cBA absorption rate constant from upper intestine and colon
    ktj_fed=2.155, #upper intestinal transit rate in fed state, /hr
    ktj_fasted=2.44, #upper intestinal transit rate in fasting state, /hr
    kti_fed=1.2, #lower intestinal transit rate in fed state, /hr
    kti_fasted=2.76, #lower intestinal transit rate in fasting state, /hr
    ge=1.2, #gall bladder emptying 1/hr
    fec=0.05, #fraction of ileal BA that is excreted via feces, equal to de novo synthesis
    aBSEP=0.839, # BSEP protein abundance, pmol/10^6 hepatocytes
    MWBSEP=140000, #molecular weight of BSEP,g/mol
    Hep=99, #hepatocellularity, 10^6 hepatocytes/g liver
    WL=1400, #weight of liver
    VmaxBSEPc=5.848, #maximum rate of BSEP-mediated efflux, umole/min/mg BSEP
    KmBSEP=4.3, #umol/L
    VmaxNTCPc=510, #maximal NTCP-mediated uptake, pmol/10^6 hepatocytes/min
    KmNTCP=5.284, #Km of NTCP-mediated uptake, umol/L
    SF_OATP=1.25, #scaling factor to correct NTCP-mediated uptake also for OATP
    VmaxASBTc=203.2, #maximal ASBT-mediated uptake, pmol/cm2/min
    KmASBT=0.662, #Km of ASBT-mediated uptake, umol/L
    CBfs_cBA=0.45 #fasting concentration in plasma
  ), 
  tout  = seq(0, 24, by = 0.1),
  state = c(Aad=0, Abo=0, Abr=0, Agu=0, 
            Ahe=0, Aki=0, Ali=0, Alu=0, 
            Amu=0, Ask=0, Asp=0, Ate=0, 
            Ave=0, Aar=0, Are=0, D= .$dose*1e-3/.$MW*1e6*70*.$Fa,
            AliClearance = 0, AkiClearance = 0,
            Aad_cBA=0, Agu_cBA=0,
            Aliw_cBA=0, Alew_cBA=0, Apv_cBA=0, Aslp_cBA=0, Ara_cBA=0, 
            AGdose=3020*sens,  Alo_cBA=0, Aup_cBA=0, Aco_cBA=0,
            Abl_cBA=0) # dose of BA produced in liver, umol #AGdose=Gdose, Ilu_BA=0
  ))

model_results
