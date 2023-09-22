.libPaths("C:/ProgramData/R/win-library/4.1")
library(deSolve)

##Bile acid model for GCDCA
##Organs: adipose,  gut (compartmentalized), liver (compartmentalized), slowly and rapidly perfused
##Compartmentalized organs: liver (intracellular and extracellular water) and intestine (upper and lower intestine + colon)

# Differential equations
model <- function (pars, tout, state, dosing){#input required for the solver
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
      
      #{Defining amount metabolitezed and cleared by the kidney for the mass balance equation}
      dAliClearance =   Cliverfree*CLmet
      dAkiClearance =   CLrenal*Cplasmavenousfree
      
      #parameters BA model
      SF_BSEP=aBSEP*MWBSEP*Hep*WL*10^-9 #scaling factor, mg BSEP/liver
      SF_NTCP=Hep*WL*10^-6 #scaling factor, 10^6 hepatocytes/liver
      SF_ASBT=4712*2.8*10^-6 #scaling factor, cm2/ileum
      
      VmaxBSEP=VmaxBSEPc*SF_BSEP*60 #umol/h/entire liver
      KmBSEPapp=KmBSEP*(1+Cliverfree/Ki) #competitive inhibition, Km is modified
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
      dAliw_cBA=VmaxNTCP*Cliverew_cBA/Kplew_cBA*BP_BA/(KmNTCP+Cliverew_cBA/Kplew_cBA*BP_BA) +fec*Aco_cBA-VmaxBSEP*Cliveriw_cBA/(KmBSEPapp+ Cliveriw_cBA)
      dAlew_cBA=Qha*(Cblood_cBA-Cliverew_cBA/Kplew_cBA*BP_BA) + Qpv*(Cportal_cBA-Cliverew_cBA/Kplew_cBA*BP_BA) -VmaxNTCP*Cliverew_cBA/Kplew_cBA*BP_BA/(KmNTCP+Cliverew_cBA/Kplew_cBA*BP_BA) 
      dApv_cBA= ka_cBA*Aup_cBA + VmaxASBT *Cileum_cBA /(KmASBT + Cileum_cBA)+ka_cBA*Aco_cBA - Qpv*Cportal_cBA
      dAslp_cBA= Qslp*(Cblood_cBA - Cslowlyper_cBA/Kpslp_cBA*BP_BA)
      dAra_cBA= Qrp*(Cblood_cBA - Crapidlyper_cBA/Kpra_cBA*BP_BA)
      dAGdose=-ge*(t<1.5)*AGdose+VmaxBSEP*Cliveriw_cBA/(KmBSEPapp+ Cliveriw_cBA)*QGb
      dAup_cBA=ge*(t<1.5)*AGdose + VmaxBSEP*Cliveriw_cBA/(KmBSEPapp+ Cliveriw_cBA)*QIb-ka_cBA*Aup_cBA-ktj*Aup_cBA
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
  } #end derivs
  solvePBPKCmax<- as.data.frame(ode(y = state, times = tout, parms = pars,
                                    atol=1e-50,
                                    func = derivs, method="lsoda")) #%>%
  #filter(Cplasmavenous == max(Cplasmavenous)) 
  return(solvePBPKCmax)
}

model



