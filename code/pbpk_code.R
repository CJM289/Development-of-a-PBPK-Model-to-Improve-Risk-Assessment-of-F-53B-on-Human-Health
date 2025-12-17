rm(list=ls())
# Clean environment and load packages
library(deSolve)
library(dplyr)
library(openxlsx)


# Experimental parameters
dose <- 5 # mg/kg
blood_v <- 0.08 # Blood volume (L/kg)
bw <- 0.03 # Body weight (kg)
route<-'ORAL'


# ======================
# PBPK Model Function
# ======================
F53B_combined <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    
    # Blood compartment
    dCb <- (Qli* blood_f * (Ali / Pliu / Wli - Cb ) - 
              Qki * blood_f *(Cb  - Aki / Pkiu / Wki) -
              Qsp * blood_f* (Cb  - Asp / Pspu / Wsp) -
              Qte * blood_f* (Cb  - Ate / Pteu / Wte) -
              Qbr * blood_f * (Cb  - Abr / Pbru / Wbr) -
              Qfa  * blood_f* (Cb - Afa / Pfau / Wfa) -
              Qhe * blood_f * (Cb  - Ahe / Pheu / Whe) -
              Qlu * blood_f* (Cb  - Alu / Pluu / Wlu) -
              Qre * blood_f* (Cb - Are / Preu / Wre)-GFR*blood_f*Cb*FF) / Wb
    
    # Gut compartment (for oral administration)
    absorbed <- Kabs * Agut 
    excreted <- Kelm* Agut 
    dAgut <- - (absorbed + excreted)
    
    # Liver compartment (with metabolism)
    dAli <- Qli * blood_f*(Cb - Ali / Pliu / Wli) - Km * Ali + absorbed
    dAex <- excreted + Km * Ali 
    
    # Kidney compartment (with filtration and reabsorption)
    
    
    
    dAki <- Qki * blood_f*(Cb  - Aki / Pkiu / Wki) +Afi*(K_reab * Km_reab) / (Km_reab + Afi)
    dAfi<-GFR*blood_f*Cb*FF-kur*Afi- Afi*(K_reab * Km_reab) / (Km_reab + Afi)
    
    dAur <- kur*Afi
    
    # Other tissue compartments
    dAsp <- Qsp * blood_f*(Cb - Asp / Pspu / Wsp)
    dAte <- Qte * blood_f* (Cb  - Ate / Pteu / Wte)
    dAbr <- Qbr * blood_f* (Cb  - Abr / Pbru / Wbr)
    dAfa <- Qfa * blood_f* (Cb  - Afa / Pfau / Wfa)
    dAhe <- Qhe * blood_f* (Cb  - Ahe / Pheu / Whe)
    dAlu <- Qlu * blood_f* (Cb   - Alu / Pluu / Wlu)
    dAre <- Qre *blood_f* (Cb - Are / Preu / Wre)
    
    # Mass balance check
    total_mass <- Cb * Wb + Ali + Aki + Asp + Ate + Abr + Alu + Afa + Ahe + Are + Agut+Afi
    
    list(c(dCb, dAli, dAki, dAsp, dAte, dAbr, dAlu, dAfa, dAhe, dAre, dAgut, dAfi,dAex, dAur),
         Inbody_Mass = total_mass)
  })
}


params<-c(BW = 0.03, GFR = 0.01668, Wli = 0.001647, Wki = 0.000501, Wb = 0.0024, 
                       Wsp = 0.000105, Wte = 0.0001794, Wbr = 0.000495, Wfa = 0.008745, 
                       Whe = 0.00015, Wlu = 0.000219, Wre = 0.0155586, Qli = 0.19159, 
                       Qki = 0.10829, Qsp = 0.004165, Qte = 0.0071162, Qbr = 0.03927, 
                       Qfa = 0.0833, Qhe = 0.07854, Qlu = 0.00595, Qre = 0.6717788, 
                       Pliu = 13.14, Pkiu = 0.46, Pspu = 0.3, Pteu = 0.24, Pbru = 0.12, 
                       Pluu = 0.82, Pfau = 0.12, Pheu = 0.29, Kelm= 0.525, blood_f = 0.067, 
          Km = 5.35e-05, K_reab = 0.77, Kabs = 0.975, 
          FF = 0.5293, Preu = 0.53, kur=0.017,
          Km_reab = 4873
)

sim_time=seq(0,24*365,24)
force_positive <- function(t, state, parms) {
  state[state < 0] <- 0  # 
  return(state)
}

if (route == 'IV') {
  state <- c(
    Cb = dose/blood_v*1000,  # Blood concentration (ng/mL)
    Ali = 0,  # Liver amount (ng)
    Aki = 0,  # Kidney amount (ng)
    Asp = 0, Ate = 0, Abr = 0, Alu = 0, Afa = 0, Ahe = 0, Are = 0,
    Agut = 0, Afi=0,Aex = 0, Aur = 0
  )
  sim <- ode(y = state, times = sim_time, func = F53B_combined, parms = params,events = list(
    func = force_positive,
    time = sim_time) )
} else if (route=='ORAL') 
  {
  state <- c(
    Cb = 0, Ali = 0, Aki = 0, Asp = 0, Ate = 0, Abr = 0, 
    Alu = 0, Afa = 0, Ahe = 0, Are = 0,
    Agut = dose*bw*1000, Afi=0, Aex = 0, Aur = 0
  )
  sim <- ode(y = state, times = sim_time, func = F53B_combined, parms = params,
             events = list(
               func = force_positive,
               time = sim_time)
