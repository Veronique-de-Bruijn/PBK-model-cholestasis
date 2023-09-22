library(tidyverse)
library(readxl)

#f.input.data<-function(){
######Intestinal uptake######
R<-1 #radius small intestine
Tsi<-3.32 #small intestine transit time 
kt<-1/(Tsi/7)
#background:
#https://pubmed.ncbi.nlm.nih.gov/10486429/
#https://pub.iapchem.org/ojs/index.php/admet/article/view/638
insilico.KaFa<-read.csv("./input/physchem_input_Ka_calculations.csv",  fileEncoding="UTF-8-BOM") %>%
  mutate(Papp.cm.power.minus.six.per.s = 10^(-4.36 - 0.01*T_PSA)*1000000)%>%
  mutate(Peff.cm.power.minus.four.per.s = (10^(0.4926*log10(Papp.cm.power.minus.six.per.s)- 0.1454)))%>%
  mutate(Peff.cm.per.hr = Peff.cm.power.minus.four.per.s/10000*3600)%>%
  mutate(Ka= (Peff.cm.per.hr*2)/R) %>%
  mutate(Fa = 1-(1+Ka/kt)**-7) %>%
  select(compound, CAS, Ka, Fa)%>%
  mutate(method.KaFa = "in silico") 

### do not use in vitro papp
# invitro.KaFa<-read.csv("./input/invitro_Papp_input.csv") %>%
#   mutate(Peff.cm.power.minus.four.per.s = (10^(0.4926*log10(PappAB)- 0.1454)))%>%
#   mutate(Peff.cm.per.hr = Peff.cm.power.minus.four.per.s/10000*3600)%>%
#   mutate(Ka= (Peff.cm.per.hr*2)/R) %>%
#   mutate(Fa = 1-(1+Ka/kt)**-7) %>%
#   select(compound, CAS, Ka, Fa)%>%
#   mutate(method.KaFa = "in vitro")


#KaFa<-rbind(insilico.KaFa,invitro.KaFa)

KaFa <- insilico.KaFa

######Hepatic clearance######
clearance<-read.csv("./input/invitro_insilico_clearance_input.csv", fileEncoding="UTF-8-BOM")  %>%
  mutate(scaling.for.liver=ifelse(reference.Clint=="pkCSM", yes=0, no=1)) 

######BSEP inhibition########
BSEPIC50 <- read.csv("input/IC50_effluxinhibition_input.csv", fileEncoding="UTF-8-BOM") %>%
  mutate(Ki=IC50_BSEP/(1+substrate.uM/4.3)) %>% 
  group_by(compound,  CAS, severity, DILI.annotation, conc.enzymes.inh, unit.conc.enzymes.inh, method.BSEP) %>%
  summarize(Ki=mean(Ki))

######Fup and partition coefficients######
source("./functions/ionization_calculations.R", local =TRUE)
source("./functions/Fup.R", local = TRUE)
source("./functions/Fuinc.R", local = TRUE)
source("./functions/partition_RodgersRowland.R", local = TRUE)
source("./functions/partition_Berezhkovskiy.R", local=TRUE)
tissue_comp<-read.csv("./input/unified_tissue_comp.csv")

phys.chem<-read.csv("./input/physchem_input.csv", fileEncoding="UTF-8-BOM") %>%
  #calculate the fractions ionized (acid/base) from these data as input to further calculations
  rowwise() %>%
  mutate(frac.unionized =  f.frac.unionized(ionization, pH = 7.4, pKa1, pKa2))%>%
  mutate(frac.ionized.acid = f.frac.ionized.acid(ionization, pH = 7.4, pKa1, pKa2)) %>%
  mutate(frac.ionized.base = f.frac.ionized.base(ionization, pH = 7.4, pKa1, pKa2)) %>%
  mutate(type = ifelse(ionization == "neutral", 1,
                       ifelse(ionization == "monoproticAcid", 2, 
                              ifelse(ionization == "monoproticBase", 3,
                                     ifelse(ionization == "diproticAcid", 4,
                                            ifelse(ionization == "diproticBase", 5, 6)))))) 

invitro.Fup <- read.csv("./input/invitro_Fup.csv", fileEncoding="UTF-8-BOM") %>% 
  drop_na(Fup) %>% 
  mutate(Fut = 1/(1+(((1-Fup)/Fup)*0.5))) %>%
  add_column(method.phys.chem.Fup = NA)

######PBK-model input calculations######
df.BP <- phys.chem %>% mutate(BP=f.BP(ionization, frac.ionized.acid, frac.ionized.base)) %>%
  mutate(method.BP="assumption")

Fup<-phys.chem %>%
  #calculate the fraction unbound in plasma based on Lobell, M.; Sivarajah, V. 
  mutate(Fup = f.Fup(ionization , frac.unionized, frac.ionized.acid, frac.ionized.base, logP, quat.nitrogen)) %>%
  mutate(Fut = 1/(1+(((1-Fup)/Fup)*0.5)))%>%
  mutate(reference.Fup = "LobellSivarajah")%>%
  mutate(original.reference.Fup = "LobellSivarajah")%>%
  mutate(method.Fup = "in silico") %>%
  select(compound, CAS, "method.phys.chem.Fup" = method.phys.chem, Fup,Fut, method.Fup, reference.Fup, original.reference.Fup) %>%
  bind_rows(.,invitro.Fup) %>%
  inner_join(., df.BP, by=c("compound", "CAS"))

Berezhkovskiy.partition.coefficient<-
  full_join(Fup, phys.chem)  %>%
  rowwise() %>%
  mutate(calcKp_Berez_Result= calcKp_Berez(logP, pKa1, pKa2, pKa3=NA,Fup, BP, type, dat = tissue_comp)) %>%
  unnest_wider(calcKp_Berez_Result) %>%
  mutate(method.partition="Berezhkovskiy") %>%
  select(compound, CAS, method.phys.chem.Fup, Fup, Fut, method.Fup,
         reference.Fup, original.reference.Fup, ionization,
         pKa1, pKa2, logP, logD, method.phys.chem, quat.nitrogen,
         frac.unionized, frac.ionized.acid, frac.ionized.base, type,
         Kpad, Kpbo, Kpbr, Kpgu, Kphe, Kpki,
         Kpli, Kplu, Kpmu, Kpsk, Kpsp, method.partition, MW, BP, method.BP)

RodgersRowland.partition.coefficient<-
  full_join(Fup, phys.chem)  %>%
  rowwise() %>%
  mutate(calcKp_RR_Result= calcKp_RR(logP, pKa1, pKa2, pKa3=NA,Fup, BP, type, dat = tissue_comp))%>%
  unnest_wider(calcKp_RR_Result) %>%
  mutate(method.partition="RodgersRowland") %>%
  select(compound, CAS, method.phys.chem.Fup, Fup, Fut, method.Fup,            
         reference.Fup, original.reference.Fup, ionization,             
         pKa1, pKa2, logP, logD, method.phys.chem, quat.nitrogen,         
         frac.unionized, frac.ionized.acid, frac.ionized.base, type,
         Kpad, Kpbo, Kpbr, Kpgu, Kphe, Kpki, 
         Kpli, Kplu, Kpmu, Kpsk, Kpsp, method.partition, MW, BP, method.BP) 

# Schmitt.partition.coefficient <-
#   full_join(insilico.Fup, phys.chem)  %>%
#   rowwise() %>%
#   mutate(calcKp_Schmitt_Result= calcKp_Schmitt(logP, pKa1, pKa2, pKa3=NA, Fup, type, dat = tissue_comp))%>%
#   unnest_wider(calcKp_Schmitt_Result) %>%
#   mutate(method.partition="Schmitt") %>%
#   select(compound, CAS, method.phys.chem.Fup, Fup, Fut, method.Fup,            
#          reference.Fup, original.reference.Fup, ionization,             
#          pKa1, pKa2, logP, logD, method.phys.chem, quat.nitrogen,         
#          frac.unionized, frac.ionized.acid, frac.ionized.base, type,
#          Kpad, Kpbo, Kpbr, Kpgu, Kphe, Kpki, 
#          Kpli, Kplu, Kpmu, Kpsk, Kpsp, method.partition)
# 
partition.coefficients<-rbind(Berezhkovskiy.partition.coefficient,
                              RodgersRowland.partition.coefficient) %>%
                              ungroup()

invivo.dose.Cmax <- read_excel("input/Invivo.xlsx", sheet = "Clinical data")

PBK.input.data<-inner_join(KaFa, partition.coefficients, by = c("compound","CAS")) %>%  
  inner_join(., clearance, by = c("compound","CAS")) %>% 
  inner_join(., BSEPIC50, by=c("compound", "CAS")) %>%
  rowwise() %>% 
  mutate(Fuinc = round(f.Fuinc(method.Clint, frac.unionized, frac.ionized.acid, frac.ionized.base, logD, logP, conc.enzymes.Clint, bound.unbound),3)) %>%
  mutate(method.Fuinc="in silico")%>% 
  mutate(Fuinh=round(f.Fuinc(method.Clint="hep", frac.unionized, frac.ionized.acid, frac.ionized.base, logD, logP, conc.enzymes.Clint=conc.enzymes.inh, bound.unbound="bound"),3)) %>%
  mutate(Ki=ifelse(method.BSEP=="SHH", yes=Ki/Fuinh, no=Ki)) %>% 
  inner_join(., invivo.dose.Cmax, by = c("compound","CAS")) %>%
  mutate(Cmax.plasma=ifelse(blood_plasma=="blood", (Cmax/BP)/MW*1000, Cmax/MW*1000)) %>%
  add_column(Cmax.plasma.units="umol/L")

write.csv(PBK.input.data, "./input/inputDataFile_final.csv")

PBK.input.data %>% group_by(compound) %>%count() %>% View()
