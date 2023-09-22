## Sensitivity analysis
.libPaths("C:/ProgramData/R/win-library/4.1")

library(deSolve)
library(tidyverse)
library(FME)

## generate input parameters
input.sens <- read.csv("./input/inputDataFile_final.csv",  fileEncoding="UTF-8-BOM") %>%
  group_by(method.partition, method.Clint, method.Fup, compound) %>% slice_sample(n=1) %>%
  drop_na(Kpad) %>% 
  mutate(dose=1)

names(input.sens)[8] <- "fup"
names(input.sens)[40] <- "CLint"
names(input.sens)[60] <-  "fuhep"
names(input.sens)[48] <- "SF"

pars <- c(
  #physiological parameters
  BW = 70,
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
  FQre = 0.103855, CLrenal=6.7,       
  
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
  SF_OATP=1.25, #scaling factor to correct NTCP-mediated uptake als for OATP
  VmaxASBTc=203.2, #maximal ASBT-mediated uptake, pmol/cm2/min
  KmASBT=0.662, #Km of ASBT-mediated uptake, umol/L
  CBfs_cBA=0.45)

parms <- as_tibble(as.list(pars)) %>%
  cbind(., input.sens)

nums <- unlist(lapply(parms, is.numeric), use.names = FALSE) 
parms <- parms[, nums] 

#step 4: differential equations
sens.model <- function (parms){#input required for the solver
  derivs <- function(t, state, parms) { 
    # returns rate of change
    with(as.list(c(state, parms)), {
      #Volume fractions
      Vad = BW*FVad	#adipose tissue 
      Vbo = BW*FVbo #bone       	 
      Vbr = BW*FVbr #brain		
      Vgu = BW*FVgu #gut          	 
      Vhe = BW*FVhe #heart        	
      Vki = BW*FVki #kidney
      Vli = BW*FVli #liver 
      Vliw =BW*FVliw #liver intracellular water
      Vlew =BW*FVlew #liver extracellular water
      Vlu = BW*FVlu #lung
      Vmu = BW*FVmu #muscle
      Vsk = BW*FVsk #skin  
      Vsp = BW*FVsp #spleen
      Vte = BW*FVte #testes
      Vve = BW*FVve #venous blood
      Var = BW*FVar #arterial blood
      Vbl = BW*(FVve+FVar-FVpv) #volume of blood excluding portal vein
      Vpv = BW*FVpv #portal vein 
      Vre = BW*FVre #rest
      Vslp = BW*(FVbo+FVsk+FVre) #slowly perfused tissue volume
      Vra=BW*(1-FVad-FVgu-FVli-FVve-FVar-FVbo-FVsk-FVre) #rapidly perfused tissue volume
      
      #Blood flows in L/hr
      Qad = QC*FQad   #adipose tissue         
      Qbo = QC*FQbo   #bone             
      Qbr = QC*FQbr   #brain             
      Qgu = QC*FQgu   #gut         
      Qhe = QC*FQhe   #heart          
      Qki = QC*FQki   #kidney           
      Qh = QC*FQh     #total blood flow to liver
      Qpv = QC*FQgu+QC*FQsp #portal vein
      Qha = QC*FQh - QC*FQgu - QC*FQsp #hepatic artery
      Qlu = QC*FQlu		#lung 
      Qmu = QC*FQmu   #muscle 
      Qsk = QC*FQsk   #skin 
      Qsp = QC*FQsp   #spleen 
      Qte = QC*FQte   #testes 
      Qre = QC*FQre		#rest
      Qslp = QC*(FQbo+FQsk+FQre)   #slowly perfused 
      Qrp = QC*(1-FQad-FQgu-FQh-FQbo-FQsk-FQre)   #rapidly perfused
      
      #Concentrations for generic model
      Cadipose = Aad/Vad		
      Cbone = Abo/Vbo		 
      Cbrain = Abr/Vbr	
      Cgut = Agu/Vgu		
      Cheart = Ahe/Vhe	 
      Ckidney = Aki/Vki		
      Cliver = Ali/Vli	 
      Clung = Alu/Vlu			 
      Cmuscle = Amu/Vmu		
      Cskin = Ask/Vsk		 
      Cspleen = Asp/Vsp		 
      Ctestes = Ate/Vte		 
      Cvenous = Ave/Vve
      Carterial = Aar/Var	
      Crest = Are/Vre 		 
      Cplasmavenous = Cvenous/BP #convert to plasma concentration, umol/L
      CLmet=(CLint/fuhep)*SF*(Vli*scaling.for.liver)*60/1000
      
      Cliverfree =(Cliver/Kpli*BP)*fup
      Ckidneyfree = Ckidney*fup
      Cplasmavenousfree = Cplasmavenous*fup
      
      dAad = Qad*(Carterial - Cadipose/Kpad*BP)										
      dAbo = Qbo*(Carterial - Cbone/Kpbo*BP)										
      dAbr = Qbr*(Carterial - Cbrain/Kpbr*BP)										
      dAgu = Ka*D + Qgu*(Carterial - Cgut/Kpgu*BP) #			
      dAhe = Qhe*(Carterial - Cheart/Kphe*BP)										
      dAki = Qki*(Carterial - Ckidney/Kpki*BP) - CLrenal*Cplasmavenousfree							
      dAli=  Qha*Carterial + Qgu*(Cgut/Kpgu*BP) + Qsp*(Cspleen/Kpsp*BP) - Qh*(Cliver/Kpli*BP) - Cliverfree*CLmet	 #sumUptake+
      dAlu = (Qad+Qbo+Qbr+Qhe+Qki+Qha+Qgu+Qsp+Qmu+Qsk+Qte+Qre)*Cvenous - (Qad+Qbo+Qbr+Qhe+Qki+Qha+Qgu+Qsp+Qmu+Qsk+Qte+Qre)*(Clung/Kplu*BP)
      dAmu = Qmu*(Carterial - Cmuscle/Kpmu*BP)
      dAsk = Qsk*(Carterial - Cskin/Kpsk*BP)	
      dAsp = Qsp*(Carterial - Cspleen/Kpsp*BP)
      dAte = Qte*(Carterial - Ctestes/Kpmu*BP)
      dAve = Qad*(Cadipose/Kpad*BP) + Qbo*(Cbone/Kpbo*BP) + Qbr*(Cbrain/Kpbr*BP) +  Qhe*(Cheart/Kphe*BP) + Qki*(Ckidney/Kpki*BP) + Qh*(Cliver/Kpli*BP)  + Qmu*(Cmuscle/Kpmu*BP) + Qsk*(Cskin/Kpsk*BP) + Qte*(Ctestes/Kpmu*BP) + Qre*(Crest/Kpmu*BP) - (Qad+Qbo+Qbr+Qhe+Qki+Qha+Qgu+Qsp+Qmu+Qsk+Qte+Qre)*Cvenous 						
      dAar = (Qad+Qbo+Qbr+Qhe+Qki+Qha+Qgu+Qsp+Qmu+Qsk+Qte+Qre)*(Clung/Kplu*BP) - (Qad+Qbo+Qbr+Qhe+Qki+Qha+Qgu+Qsp+Qmu+Qsk+Qte+Qre)*Carterial
      dAre = Qre*(Carterial - Crest/Kpmu*BP)		
      dD = - Ka*D
      
      #{Defining amount metabolized and cleared by the kidney for the mass balance equation}
      dAliClearance =   Cliverfree*CLmet
      dAkiClearance =   CLrenal*Cplasmavenousfree
      
      #parameters BA model
      #Ki=(IC50_BSEP/Fuinh)/(1+(60/KmBSEP)) #competitive inhibition, Cheng-Prusoff equation
      SF_BSEP=aBSEP*MWBSEP*Hep*WL*10^-9 #scaling factor, mg BSEP/liver
      SF_NTCP=Hep*WL*10^-6 #scaling factor, 10^6 hepatocytes/liver
      SF_ASBT=4712*2.8*10^-6 #scaling factor, cm2/ileum
      
      VmaxBSEP=VmaxBSEPc*SF_BSEP*60 #umol/h/entire liver
      VmaxBSEPapp=VmaxBSEP#maximal BSEP-mediated efflux, umole/liver/hr
      KmBSEP=KmBSEP*(1+(Cliver/Kpli*BP)/Ki)
      VmaxNTCP=VmaxNTCPc*SF_NTCP*SF_OATP*60 #maximal NTCP-mediated uptake rate, umole/liver/hr
      VmaxASBT=VmaxASBTc*SF_ASBT*60 #maximal ASBT-mediated uptake rate, umole/ileum/hr
      
      QIb=1-QGb #fraction of bile acids going directly from gallbladder to intestinal lumen
      
      ##definining fast or fed state
      kti=(t<1.5)*kti_fed+(t>=1.5)*kti_fasted 
      ktj=(t<1.5)*ktj_fed+(t>=1.5)*ktj_fasted
      
      # Concentrations cBA model
      Cadipose_cBA = Aad_cBA /Vad
      Cgut_cBA  = Agu_cBA /Vgu
      Cileum_cBA=Alo_cBA/Vil
      Cliverew_cBA= Alew_cBA / Vlew
      Cliveriw_cBA=Aliw_cBA/Vliw
      Cslowlyper_cBA=Aslp_cBA/Vslp
      Crapidlyper_cBA=Ara_cBA/Vra
      Cportal_cBA=Apv_cBA/Vpv
      Cblood_cBA  = Abl_cBA /Vbl
      Cplasmavenous_cBA  = Cblood_cBA/BP_BA + CBfs_cBA
      
      #differential equations cBA
      dAad_cBA=Qad*(Cblood_cBA - Cadipose_cBA/Kpad_cBA*BP_BA)
      dAgu_cBA = Qgu*(Cblood_cBA - Cgut_cBA/Kpgu_cBA*BP_BA)
      dAliw_cBA=VmaxNTCP*Cliverew_cBA/Kplew_cBA*BP_BA/(KmNTCP+Cliverew_cBA/Kplew_cBA*BP_BA) +fec*Aco_cBA-VmaxBSEPapp*Cliveriw_cBA/(KmBSEP+ Cliveriw_cBA)
      dAlew_cBA=Qha*(Cblood_cBA-Cliverew_cBA/Kplew_cBA*BP_BA) + Qpv*(Cportal_cBA-Cliverew_cBA/Kplew_cBA*BP_BA) -VmaxNTCP*Cliverew_cBA/Kplew_cBA*BP_BA/(KmNTCP+Cliverew_cBA/Kplew_cBA*BP_BA) 
      dApv_cBA= ka_cBA*Aup_cBA + VmaxASBT *Cileum_cBA /(KmASBT + Cileum_cBA)+ka_cBA*Aco_cBA - Qpv*Cportal_cBA
      dAslp_cBA= Qslp*(Cblood_cBA - Cslowlyper_cBA/Kpslp_cBA*BP_BA)
      dAra_cBA= Qrp*(Cblood_cBA - Crapidlyper_cBA/Kpra_cBA*BP_BA)
      dAGdose=-ge*(t<1.5)*AGdose+VmaxBSEPapp*Cliveriw_cBA/(KmBSEP+ Cliveriw_cBA)*QGb
      dAup_cBA=ge*(t<1.5)*AGdose + VmaxBSEPapp*Cliveriw_cBA/(KmBSEP+ Cliveriw_cBA)*QIb-ka_cBA*Aup_cBA-ktj*Aup_cBA
      dAlo_cBA= ktj*Aup_cBA - kti*Alo_cBA -VmaxASBT *Cileum_cBA /(KmASBT + Cileum_cBA)
      dAco_cBA= kti*Alo_cBA-ka_cBA*Aco_cBA-fec*Aco_cBA
      dAbl_cBA= Qad*Cadipose_cBA/Kpad_cBA*BP_BA+Qgu*Cgut_cBA/Kpgu_cBA*BP_BA+Qh*Cliverew_cBA/Kplew_cBA*BP_BA+Qslp*Cslowlyper_cBA/Kpslp_cBA*BP_BA+Qrp*Crapidlyper_cBA/Kpra_cBA*BP_BA-(Qad+Qha+Qslp+Qrp+Qgu)*Cblood_cBA
      
      #{Mass Balance}
      #Generic
      Total = dose
      Calculated = Aad+Abo+Abr+Agu+Ahe+Aki+Ali+Alu+
        Amu+Ask+Asp+Ate+Ave+Aar+Are+
        D+AliClearance+AkiClearance
      ERROR=((Total-Calculated)/Total+1E-30)*100
      MASSBBAL=Total-Calculated + 1 
      
      #Bile acids
      Total_BA=Gdose
      Calculated_cBA=Aad_cBA+Agu_cBA+Alew_cBA+Aliw_cBA+Ara_cBA+Aslp_cBA+AGdose+Aup_cBA+Alo_cBA+Aco_cBA+Abl_cBA+Apv_cBA
      ERROR_BA=((Total_BA-Calculated_cBA)/Total_BA+1E-30)*100
      MASSBBAL_BA=Total_BA-Calculated_cBA + 1 
      
      res<-c(dAad, dAbo, dAbr, dAgu, 
             dAhe, dAki, dAli, dAlu, 
             dAmu, dAsk, dAsp, dAte, 
             dAve, dAar, dAre, dD, 
             dAliClearance, dAkiClearance, dAad_cBA, dAgu_cBA,
             dAliw_cBA, dAlew_cBA, dApv_cBA, dAslp_cBA, dAra_cBA, 
             dAGdose,  dAlo_cBA, dAup_cBA, dAco_cBA,
             dAbl_cBA)
      
      return(list(res, ERROR = ERROR, 
                  MASSBBAL = MASSBBAL,
                  Calculated = Calculated,
                  Cplasmavenous = Cplasmavenous,
                  Cadipose = Cadipose, 
                  Cbone = Cbone,	 
                  Cbrain = Cbrain,	
                  Cgut = Cgut,		
                  Cheart = Cheart,	 
                  Ckidney = Ckidney,		
                  Cliver = Cliver,	 
                  Clung = Clung,			 
                  Cmuscle = Cmuscle,		
                  Cskin = Cskin,	 
                  Cspleen = Cspleen,		 
                  Cgonads = Ctestes,		 
                  Cvenous = Cvenous,		
                  Carterial = Carterial, 
                  Cplasmavenous_cBA=Cplasmavenous_cBA,
                  Cadipose_cBA = Cadipose_cBA,
                  Cgut_cBA = Cgut_cBA,
                  Cliverew_cBA= Cliverew_cBA,
                  Cliveriw_cBA=Cliverew_cBA,
                  Cileum_cBA=Cileum_cBA, 
                  Cblood_cBA = Cblood_cBA,
                  Total_BA=Total_BA,
                  Calculated_cBA=Calculated_cBA,
                  ERROR_BA=ERROR_BA,
                  MASSBBAL_BA=MASSBBAL_BA))
      
    }) #end with
  }
  
  state <- c(Aad=0, Abo=0, Abr=0, Agu=0, 
             Ahe=0, Aki=0, Ali=0, Alu=0, 
             Amu=0, Ask=0, Asp=0, Ate=0, 
             Ave=0, Aar=0, Are=0, D= parms$dose/parms$MW*parms$BW*parms$Fa*1e3,
             AliClearance = 0, AkiClearance = 0,
             Aad_cBA=0, Agu_cBA=0,
             Aliw_cBA=0, Alew_cBA=0, Apv_cBA=0, Aslp_cBA=0, Ara_cBA=0, 
             AGdose=parms$Gdose,  Alo_cBA=0, Aup_cBA=0, Aco_cBA=0,
             Abl_cBA=0)
  tout  <-  seq(0, 16, by = 0.1)
  solvePBPKCmax<- as.data.frame(ode(y = state, times = tout, parms = parms,
                                    atol=1e-50,
                                    func = derivs)) #%>%
    summarize(auc=sum(Aliw_cBA*0.1)) ##to speed up calculations, remove if interested in another endpoint
  return(solvePBPKCmax)
}

out <- sens.model(parms[1,])
plot(out$MASSBBAL_BA)
plot(out$Cliveriw_cBA)


##check whether the function works for one compound
senspar <- c("Ka", "Fa", "fup", "Fut", "Kpad", "Kpbo", "Kpbr", "Kpgu", "Kphe", "Kpki", "Kpli", "Kplu", "Kpsk", "Kpsp", "BP", "CLint", "SF", "Ki", "fuhep", "dose", names(pars))
lsa <- sensFun(sens.model, parms[1,], sensvar="auc", senspar=senspar, tiny=0.05)

## continue if lsa returns data

## run the sensitivity analysis for all compounds. Might take several hours, depending on your input
senspar <- c("Ka", "Fa", "fup", "Fut", "Kpad", "Kpbo", "Kpbr", "Kpgu", "Kphe", "Kpki", "Kpli", "Kplu", "Kpsk", "Kpsp", "BP", "CLint", "SF", "Ki", "fuhep", "dose", 
             "BW", "FVad", "FVbo", "FVbr", "FVgu", "FVhe", "FVki", "FVli", "FVliw", "FVlew", "FVlu", "FVmu", "FVsk", "FVsp", "FVte", "FVve", "FVpv", "FVar", "FVre",
             "QC", "FQad", "FQbo", "FQbr", "FQgu", "FQhe", "FQki", "FQh", "FQlu", "FQmu", "FQsk", "FQsp", "FQte", "FQre", "CLrenal")
df <- data.frame(matrix(NA, nrow = nrow(parms), ncol = length(senspar) + 2))
colnames(df) <- c("x", "var", senspar)
result_list <- vector("list", nrow(parms))
pb <- txtProgressBar(min = 0, max = nrow(parms), style = 3)

# Record the start time
start_time <- Sys.time()

for(i in 1:nrow(parms)) {
  lsa <- sensFun(sens.model, parms[i, ], sensvar=c("Cplasmavenous"), senspar=senspar, tiny=0.05) 
  result_list[[i]] <- lsa
  setTxtProgressBar(pb, i)
  
  # Calculate time taken for this iteration
  elapsed_time <- difftime(Sys.time(), start_time, units = "secs")
  time_per_iteration <- elapsed_time / i
  
  # Calculate remaining time
  iterations_left <- nrow(parms) - i
  time_left <- time_per_iteration * iterations_left
  
  # Calculate remaining time in hours
  iterations_left <- nrow(parms) - i
  time_left <- time_per_iteration * iterations_left
  time_left_hours <- as.numeric(time_left, units = "hours")
  
  cat("Estimated time left:", round(time_left_hours, 2), "hours\n")
}

#Close progress bar
close(pb)

df <- do.call(rbind, result_list)
write.csv(df, "results/sensitivityanalysis_Cmax.csv")