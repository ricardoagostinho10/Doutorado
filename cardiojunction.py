import numpy as np
import scipy.integrate as inte
import math
from math import exp
import warnings
warnings.filterwarnings('ignore')
import os
import pandas as pd
import json
from scipy.integrate import solve_ivp

class cardiojunction:
    contador = 0;
    f = 0; 
    f2 = 0;
    w = 0; 
    w2 = 0; 
    t_ap = 0; 
    L = 0; 
    Delay = 0; 
    Delay2 = 0; 
    Ap = 0; 
    tap = 0; 
    tsim = 0; 
    protocol = 0; 
    A_inj = 0; 
    A_inj2 = 0; 
    BLOCKSRPUMP = 0;
    Model_Na = 0; 
    Model_Ito = 0; 
    STIMULUSRPUMP = 0; 
    BLOCKNCX = 0; 
    STIMULUSNCX = 0; 
    CAFEINA = 0; 
    BLOCKCICR = 0; 
    Model_ICaL = 0; 
    GNaL = 0; 
    PCam = 0; 
    v_resting = 0; 
    D = 0; 

    alfaneg = 0;
    betaneg = 0; 
    Lsarc = 0; 
    Model_Force = 0; 
    alfaaj = 0; 
    betaaj = 0; 
    Fmax = 0; 
    Lajust = 0; 
    Model_IKr = 0; 
    Model_IKs = 0; 
    GKs = 0; 
    cellLength = 0; 
    BLOCKIKs = 0; 
    BLOCKIKr = 0; 
    BLOCKItos = 0;  
    BLOCKItof = 0; 
    BLOCKINa = 0; 
    BLOCKICaL = 0; 
    STIMULUSICaL = 0;
    FoRT = 0; 
    Qpow = 0;  
    Nao = 0; 
    Cao = 0; 
    Ko = 0; 
    TSt = 0; 
    fit = 0; 

    GNa = 0; 
    IbarNaK = 0; 
    KmNaip = 0; 
    KmKo = 0; 
    pNaK  = 0;
    GtoSlow = 0;
    GtoFast = 0; 
    gkp = 0; 
    Fjunc = 0; 
    Fsl = 0; 
    gkr = 0; 
    kiss = 0;
    kisss = []
    Fjunc_CaL = 0; 
    Fsl_CaL = 0; 
    pCa  = 0;
    GSrLeak = 0; 
    gCai = 0; 
    gCao = 0;    

    Q10CaL = 0;
    J_na_juncsl = 0;
    Vjunc = 0; 
    GNam = 0; 
    I_scale_Na = 0; 
    GItof = 0; 
    GItos = 0; 
    fatorgks = 0; 
    GKss = 0; 
    GKur = 0;

    IbarNCX = 0;
    Q10NCX = 0;
    ksat = 0; 
    nu = 0; 
    KmNai = 0; 
    KmNao = 0; 
    KmCai = 0; 
    KmCao = 0; 
    Kdact = 0; 
    Q10SLCaP = 0; 
    IbarSLCaP = 0;
    KmPCa = 0; 
    ks = 0; 
    GCaB = 0; 
    ICa_scale = 0;
    

    Kmf = 0;
    Kmr = 0; 
    hillSRCaP = 0; 
    Q10SRCaP = 0; 
    Vmax_SRCaP = 0; 
    Frdy = 0; 
    Kp = 0; 
    Kw = 0; 
    Ke = 0; 
    Lo = 0; 
    Le = 0; 
    a = 0
    pulso = []
    sol = []
    tempo = []
    novoTempo = []
    vetoriKr = []
    repouso = 0
    step = 0
    acumuladorCorrente = []    
	
    def vetorPulso(self, largura, delay, L, t_ap):
        item = t_ap
        cont = 0
        while item <= L:
            self.pulso.append(item)
            if cont % 2 == 0:
                item = item + largura
            else:
                item = item + delay
            cont = cont + 1

                    
    def pulse_train(self,valor):
        cont = 0
        
        while cont < len(self.pulso):
            if(valor <= self.pulso[cont]):
                break
            cont+=1
   
        if cont % 2 == 0:
            
            return 0
        else:
            return 1 
            
    def calculos(self,y,t):
        ydot = np.zeros(shape=(113)); # Numero de Equacoes Diferenciais Acopladas
#################################################################################################################################################
 
## Protocolo para simulacao
        
        if self.protocol==1:
            I_app = 0;     

            
            if t <= self.t_ap:                                #[ms]
    
                A = -self.v_resting + self.Ap;    
                sig = self.pulse_train(t)
                v = A * sig + self.v_resting                     # Amplitude do pulso
                
            elif t > self.t_ap and t <= self.t_ap + self.L:    
                A = -self.v_resting + self.Ap;    
                sig = self.pulse_train(t)
                v = A * sig + self.v_resting                     # Amplitude do pulso

            else:
    
                A = -self.v_resting + self.Ap;    
                sig = self.pulse_train(t)
                v = A * sig + self.v_resting                     # Amplitude do pulso
   
            
        elif self.protocol == 2:  # voltagem clamp com voltagem controlada e injecao de corrente      
        
            v = y[38]; 
        
            if t <= self.t_ap:                                            #[ms]
    
                vv = self.v_resting;                                           #[mV]        

            elif t > self.t_ap and t <= self.t_ap + self.L:    
    
                A = -self.v_resting + self.Ap;                                            # Amplitude do pulso

                sig = self.pulse_train(t)
                vv = A * sig + self.v_resting
    
            else:
    
                vv = self.v_resting;    #[mV]                
   

            V_clamp = vv;
            R_clamp = 0.02;
            I_app = (V_clamp-y[38])/R_clamp;

      
    
        elif self.protocol==3:            		
            v = y[38];    #Injecao de corrente    
			
            if ((t > self.t_ap) and t < self.t_ap + self.D) and (np.mod((t + (1000 - self.t_ap)),(self.w + self.Delay)))<= self.w:                                
                I_app = self.A_inj;        
                			
            elif (t > self.t_ap + self.D) and (np.mod(t + (1000 - self.t_ap),(self.w2 + self.Delay2))<= self.w2):       ##ok<AND2> #corrente injetada de 10 uA/uF com frequencia de 4 Hz durante 1s                    
                I_app = self.A_inj2;                    
            else:            
                I_app = 0;
        
#################################################################################################################################################

## Parametros do modelo

        # Constantes
        self.R = 8314.0;                       # [J/kmol*K]    Constante Universal dos Gases  
        self.Frdy = 96485.0;                   # [C/mol]       Constante de Farad
        self.Temp = 310.0;                     # [K]           Temperatura
        self.FoRT = self.Frdy/self.R/self.Temp;
        self.Cmem = 1.3810e-10;              # [F] Capacitancia da Membrana
        self.Qpow = (self.Temp-310)/10;           # Fator Q10
        
        # Geometria da Celula
        cellRadius = 10.25;                         # Raio da celula [um]
        junctionLength = 160e-3;                    # Comprimento da Juncao(Cleft) [um]
        junctionRadius = 15e-3;                     # Raio do Cleft [um]
        Vcell = math.pi*cellRadius**2*self.cellLength*1e-15;   # [L]
        Vmyo = 0.65*Vcell; 
        Vsr = 0.035*Vcell; 
        Vsl = 0.02*Vcell; 
        self.Vjunc = 0.0539*.01*Vcell; 
        SAjunc = 20150*math.pi*2*junctionLength*junctionRadius;    ##ok<NASGU> # [um**2]
        SAsl = math.pi*2*cellRadius*self.cellLength;                    ##ok<NASGU> # [um**2]

        # Parametros de difusao
        J_ca_juncsl = 1/1.2134e12;                  # [L/msec] = 8.2413e-13
        J_ca_slmyo = 1/2.68510e11;                  # [L/msec] = 3.2743e-12
        self.J_na_juncsl = 0;#1/(1.6382e12/3*100);          # [L/msec] = 6.1043e-13
        J_na_slmyo = 1/(1.8308e10/3*100);           # [L/msec] = 5.4621e-11
        
        # Fractional currents in compartments
        self.Fjunc = 0;   
        self.Fsl = 1-self.Fjunc;
        self.Fjunc_CaL = 1; 
        self.Fsl_CaL = 1-self.Fjunc_CaL;

        # Fixed ion concentrations
        
        self.Cli = 15;                       # Intracellular Cl  [mM]
        self.Clo = 150;                      # Extracellular Cl  [mM]
        self.Ko = 5.4;                       # Extracellular K   [mM]
        self.Nao = 140;                      # Extracellular Na  [mM]
        self.Cao = 1.8;                      # Extracellular Ca  [mM]
        self.Mgi = 0.5;                      # Intracellular Mg  [mM]
        
        # Nernst Potentials
        ena_junc = (1/self.FoRT)*np.log(self.Nao/y[31]);     # [mV]
        ena_sl = (1/self.FoRT)*np.log(self.Nao/y[32]);       # [mV]
        ek = (1/self.FoRT)*np.log(self.Ko/y[34]);            # [mV]
        eca_junc = (1/self.FoRT/2)*np.log(self.Cao/y[35]);   # [mV]
        eca_sl = (1/self.FoRT/2)*np.log(self.Cao/y[36]);     # [mV]
        ecl = (1/self.FoRT)*np.log(self.Cli/self.Clo);            # [mV]
        
        ## Na current parameters
        
        # Na current Fast - Model Bondarenko
        self.GNam = 13*(1 - self.BLOCKINa);# 16
        self.I_scale_Na = 1; 
        
        # Na current Fast - Model HH
        self.GNa = 13*(1 - self.BLOCKINa); 
        
        # Na current Later - Model HH
        self.GNaL = 0*(1 - self.BLOCKINa);               # substitui 0.0045 por 0 de Bondarenko
        
        # Na current background - Model HH
        self.GNaB = 0.0026;                            # [mS/uF] substitui 0.297e-3 por 0.0026 de Bondarenko 
        
        
        # Na transport parameters -  Na/K Pump Current
        self.IbarNaK = 0.88;                          # [uA/uF] 1.90719 POR 0.88  
        self.KmNaip = 21;                                # [mM]  11 POR 21 DE BONDARENKO 2013
        self.KmKo = 1.5;                                 # [mM]
        
        ## K current parameters
        
        # K current itof e itos - Model HH
        self.GtoSlow = 0*(1 - self.BLOCKItos);                     # [mS/uF] 0.06 por 0.0629 Bondarenko 
        self.GtoFast = 0.4067*(1 - self.BLOCKItof);                     # [mS/uF] 0.02 por 0.4067 Bondarenko 
        
        # K current itof e itos - Model Winslow
        self.GItof = 0.4067*(1 - self.BLOCKItof);    #0.01;# 0.0775; #mS/uF substitui 0.015 por 0.4067   0.015 por 0.0798   TESTE 0.25
        self.GItos = 0*(1 - self.BLOCKItos); #4.161e-8;  substitui 1.161e-8 por 0.0629
        
        # K current IKs - Model HH
        self.pNaK = 0.01833; 
        
        # K current IKs - Model Winslow
        self.GKs = 0.00575*(1 - self.BLOCKIKs);  #mS/uF 0.0095
        
        # K current IKs - Model Severi
        self.fatorgks = 0.00575;# ou 0.0779 0.000575 BONDARENKO
        
        # K current IKr - Model Winslow e Model HH
        self.gkr = 0.078*(1 - self.BLOCKIKr)*np.sqrt(self.Ko/5.4); #0.3 BONDARENKO
        
        # K current plato - Model HH
        self.gkp = 0;#0.001;                                     # Nao modelado no Bondarenko
        
        # K current IKur - Model HH
        self.GKur = 0.160; #0 por 0.0975 Bondarenko 0.160
        
        # K current IKss - Model HH
        self.GKss = 0.0324; #0.050; # mS/uF 0 por 0.0324 Bondarenko
        
        ## Ca Current parameters
        
        # Ca current ICaL - Model HH
        self.pCa = 5.4e-4*(1 - self.BLOCKICaL)*(1 + self.STIMULUSICaL);  
        
        # Ca current ICaL - Model Mahajan
        self.PCam = 100*24.3e-6*(1 - self.BLOCKICaL)*(1 + self.STIMULUSICaL);                 # [cm/sec]
        self.Q10CaL = 1.8;
        self.ICa_scale = 1*self.Q10CaL**self.Qpow;          ###ok<NASGU>
        aff = 1;
        self.gCai = 0.0341;    
        self.gCao =0.341; 
        
        # Correntes de sodio e potassio ativadas por Ca
        self.pNa = 1.5e-8;                       # [cm/sec]
        pK = 2.7e-7;                        # [cm/sec]
        self.gKi = 0.75;
        self.gKo = 0.75;
        self.gNai = 0.75;
        self.gNao = 0.75;
        
        # Corrente de Cl ativada por Ca
        GClCa = 10;                   # [mS/uF] 0.109625 por 10
        GClB = 0;                        # [mS/uF] 9e-3 por 0 Bondarenko
        KdClCa = 10e-3;                    # [mM] 100e-3 por 10e-3 bondarendo
        
        # Transportadores de Ca
        
        # INCX Parameters                                                                  #Conferir dados do NCX Shannon vs Bondarenko
        self.IbarNCX = 9.0*(1 - self.BLOCKNCX)*(1 + self.STIMULUSNCX);     # [uA/uF] # conferir com Bondarenko
        self.KmCai = 3.59e-3;                                    # [mM]
        self.KmCao = 1.38;                                        # [mM]1.3 por 1.38 de Bondarenko
        self.KmNai = 12.29;                                      # [mM]
        self.KmNao = 87.5;                                       # [mM]
        self.ksat = 0.27;                                        # [none]  
        self.nu = 0.35;                                          # [none]
        self.Kdact = 0.256e-3;                                   # [mM] 
        self.Q10NCX = 1.57;                                      # [none]
        
        # Bomba de Ca do sarcolema
        self.IbarSLCaP = 1;                                 # [uA/uF](2.2 umol/L cytosol/sec) 0.0673 por 1 
        self.KmPCa = 0.5e-3;                                     # [mM] 
        self.Q10SLCaP = 2.35;                                    # [none]
        
        
        ## Fluxo de Ca do RS
        
        # Leak do RS
        self.GSrLeak = 5.348e-6;
        
        # Bomba de Ca do RS 
        self.Q10SRCaP =  2.6;                                                # [none]
        self.Kmf = 0.246e-3;                                                 # [mM] default
        self.Kmr = 1.7;                                                      # [mM]L cytosol
        self.hillSRCaP = 1.787;                                              # [mM]
        
        # Liberacao de Ca do SR
        self.ks = 25;                                                        # [1/ms]      
        koCa = 10;                                                      # [mM**-2 1/ms]   #default 10   modified 20
        kom = 0.06;                                                     # [1/ms]     
        kiCa = 0.5;                                                     # [1/mM/ms]
        kim = 0.005;                                                    # [1/ms]
        ec50SR = 0.45;                                                  # [mM]

        if self.CAFEINA==1:
            koCa=koCa*7.5;
            self.GCaB=0;
            self.Vmax_SRCaP=0;    
        else:
            koCa = 10;  
            self.Vmax_SRCaP = 5.3114e-3*(1 - self.BLOCKSRPUMP)*(1 + self.STIMULUSRPUMP); # 5.3114e-3 por 0.45e-3 bondarenko
            self.GCaB = 2.513e-4;


        # Buffering parameters  
        Bmax_Naj = 7.561;                                       # [mM]  # Na buffering
        Bmax_Nasl = 1.65;                                       # [mM]
        koff_na = 1e-3;                                         # [1/ms]
        kon_na = 0.1e-3;                                        # [1/mM/ms]
        Bmax_TnClow = 70e-3;                                    # [mM]                      # TnC low affinity
        koff_tncl = 19.6e-3;                                    # [1/ms] 
        kon_tncl = 32.7;                                        # [1/mM/ms]
        Bmax_TnChigh = 140e-3;                                  # [mM]                      # TnC high affinity 
        koff_tnchca = 0.032e-3;                                 # [1/ms] 
        kon_tnchca = 2.37;                                      # [1/mM/ms]
        koff_tnchmg = 3.33e-3;                                  # [1/ms] 
        kon_tnchmg = 3e-3;                                      # [1/mM/ms]
        Bmax_CaM = 24e-3;                                       # [mM]  # CaM buffering
        koff_cam = 238e-3;                                      # [1/ms] 
        kon_cam = 34;                                           # [1/mM/ms]
        Bmax_myosin = 140e-3;                                   # [mM]                      # Myosin buffering
        koff_myoca = 0.46e-3;                                   # [1/ms]
        kon_myoca = 13.8;                                       # [1/mM/ms]
        koff_myomg = 0.057e-3;                                  # [1/ms]
        kon_myomg = 0.0157;                                     # [1/mM/ms]
        Bmax_SR = 19*.9e-3;                                     # [mM] (Bers text says 47e-3) 19e-3
        koff_sr = 60e-3;                                        # [1/ms]
        kon_sr = 100;                                           # [1/mM/ms]
        Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;                        # [mM]    # SL buffering
        Bmax_SLlowj = 4.6e-3*Vmyo/self.Vjunc*0.1;                    # [mM]    
        koff_sll = 1300e-3;                                     # [1/ms]
        kon_sll = 100;                                          # [1/mM/ms]
        Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;                       # [mM] 
        Bmax_SLhighj = 1.65e-3*Vmyo/self.Vjunc*0.1;                  # [mM] 
        koff_slh = 30e-3;                                       # [1/ms]
        kon_slh = 100;                                          # [1/mM/ms]
        Bmax_Csqn = 140e-3*Vmyo/Vsr;                            # [mM] # Bmax_Csqn = 2.6;      # Csqn buffering
        koff_csqn = 65;                                         # [1/ms] 
        kon_csqn = 100;                                         # [1/mM/ms] 
        
        ## Membrane Currents
        
        # Current Late I_Na - Model HH
        aml = 0.32*(v+47.13)/(1-exp(-0.1*(v+47.13))); 
        bml = 0.08*exp(-v/11);                            
        hlinf = 1/(1+exp((v+91)/6.1));
        tauhl = 600;
        
        ydot[46] = aml*(1-y[46])-bml*y[46];
        ydot[47] = (hlinf-y[47])/tauhl;
        
        I_NaL_junc = self.Fjunc*self.GNaL*y[46]**3*y[47]*(v-ena_junc);
        I_NaL_sl = self.Fsl*self.GNaL*y[46]**3*y[47]*(v-ena_sl);
        
        # I_Na: Fast Na Current - Model Bondarenko

        if self.Model_Na == 1:
        
            alfa11 = 3.802/(0.1027*exp(-(v+2.5)/17)+ 0.2*exp(-(v+2.5)/150)); #3.802 original e alterei para 5
            alfa12 = 3.802/(0.1027*exp(-(v+2.5)/15)+ 0.23*exp(-(v+2.5)/150)); #3.802
            alfa13 = 3.802/(0.1027*exp(-(v+2.5)/12)+ 0.25*exp(-(v+2.5)/150)); #3.802
            beta11 = 0.1917*exp(-(v+2.5)/20.3);
            beta12 = 0.20*exp(-(v-2.5)/20.3); #beta12 = 0.20*exp(-(v-2.5)/20.3); no original
            beta13 = 0.22*exp(-(v-7.5)/20.3); #beta13 = 0.22*exp(-(v-7.5)/20.3); 
            alfa3 = 7e-7*exp(-(v+7)/7.7);
            beta3 = 0.0084 + 0.00002*(v+7);
            alfa2 = 1/(0.188495*exp(-(v+7)/16.6)+ 0.393956);
            beta2 = (alfa13*alfa2*alfa3)/(beta13*beta3);
            alfa4 = alfa2/1000;
            alfa5 = alfa2/1000;
            alfa5 = alfa2/95000;
            beta4 = alfa3;
            beta5 = alfa3/50; #/50
 
            ONa = 1 - (y[48] + y[49] + y[50]+ y[51] + y[52] + y[53] + y[54] + y[55]);
            
            ydot[48] = alfa3*y[51] + alfa12*y[49] + beta13*ONa - (beta3 + beta12 + alfa13)*y[48]; # CNa1 = z[0]
            ydot[49] = alfa11*y[50] + beta12*y[48] + alfa3*y[54] - (beta11 + alfa12 + beta3)*y[49];  # CNa2 = z[1]
            ydot[50] = alfa3*y[55] + beta11*y[49] - (beta3 + alfa11)*y[50]; # CNa3 = z[2]
            ydot[51] = alfa12*y[54] + beta4*y[52] + beta3*y[48] + alfa2*ONa - (beta12 + alfa4 + alfa3 + beta2)*y[51]; #IFNa = z[3]
            ydot[52] = alfa4*y[51] + beta5*y[53] - (alfa5 + beta4)*y[52];  # I1Na = z[4]
            ydot[53] = alfa5*y[52] - beta5*y[53];  # I2Na = z[5]
            ydot[54] = alfa11*y[55] + beta12*y[51] + beta3*y[49] - (beta11 + alfa12 + alfa3)*y[54];     # ICNa2 =z[6]
            ydot[55] = beta3*y[50] + beta11*y[54] - (alfa11 + alfa3)*y[55];  #ICNa3 = z[7]
            
            ENa = (1/self.FoRT)*np.log((0.9*self.Nao + 0.1*self.Ko)/(0.9*y[32] + 0.1*y[34])); #y[34]
            I_Na_sl = self.Fsl*self.GNam*ONa*(v - ENa)*self.I_scale_Na + I_NaL_sl;
            I_Na_junc =0;
            
            I_Na = I_Na_junc + I_Na_sl;                                                   ##ok<NASGU>
        
        else:
    
            #I_Na: Fast Na Current - Model HH
            am = 0.32*(v+47.13)/(1-exp(-0.1*(v+47.13)));
            bm = 0.08*exp(-v/11);
        
            if v >= -40:
                ah = 0; aj = 0;
                bh = 1/(0.13*(1+exp(-(v+10.66)/11.1)));
                bj = 0.3*exp(-2.535e-7*v)/(1+exp(-0.1*(v+32)));
            else:
                ah = 0.135*exp((80+v)/-6.8);
                bh = 3.56*exp(0.079*v)+3.1e5*exp(0.35*v);
                aj = (-1.2714e5*exp(0.2444*v)-3.474e-5*exp(-0.04391*v))*(v+37.78)/(1+exp(0.311*(v+79.23)));
                bj = 0.1212*exp(-0.01052*v)/(1+exp(-0.1378*(v+40.14)));
            
        
            ydot[0] = am*(1-y[0])-bm*y[0];
            ydot[1] = ah*(1-y[1])-bh*y[1];
            ydot[2] = aj*(1-y[2])-bj*y[2];
                
            I_Na_junc = self.Fjunc*self.GNa*y[0]**3*y[1]*y[2]*(v-ena_junc) + I_NaL_junc;
            I_Na_sl = self.Fsl*self.GNa*y[0]**3*y[1]*y[2]*(v-ena_sl)+ I_NaL_sl;


        I_Na = I_Na_junc + I_Na_sl;                                                   ##ok<NASGU>
        
        
        # I_nabk: Na Background Current - Model HH
        I_nabk_junc = self.Fjunc*self.GNaB*(v-ena_junc);
        I_nabk_sl = self.Fsl*self.GNaB*(v-ena_sl);
        I_nabk = I_nabk_junc+I_nabk_sl;                                             ##ok<NASGU>
        
        # I_nak: Na/K Pump Current 
        sigma = (exp(self.Nao/67.3)-1)/7;
        fnak = 1/(1+0.1245*exp(-0.1*v*self.FoRT)+0.0365*sigma*exp(-v*self.FoRT));
        I_nak_junc = self.Fjunc*self.IbarNaK*fnak*self.Ko /(1+(self.KmNaip/y[31])**4) /(self.Ko+self.KmKo);
        I_nak_sl = self.Fsl*self.IbarNaK*fnak*self.Ko /(1+(self.KmNaip/y[32])**4) /(self.Ko+self.KmKo);
        I_nak = I_nak_junc+I_nak_sl;
        
        # I_kr: Rapidly Activating K Current
        
        if self.Model_IKr == 1:  # Model Winslow
        
            alfa0r = 0.0171*exp(0.0330*v);# ms**-1
            beta0r = 0.0397*exp(-0.0431*v);# ms**-1
            alfa1r = 0.0206*exp(0.00262*v);# ms**-1 #0.0262
            beta1r = 0.0013*exp(-0.0269*v);# ms**-1
            alfair = 0.1067*exp(0.0057*v);# ms**-1
            betair = 0.0065*exp(-0.0454*v);# ms**-1
            alfai3r = 8.04e-5*exp(6.98e-7*v); # ms**-1 e-7
            kfr = 0.0261; # ms**-1
            kbr = 0.1483; # ms**-1
            psir = (beta1r*betair*alfai3r)/(alfa1r*alfair);
            
            OKr = 1 - (y[86] + y[87] + y[88] + y[89]);
            
            ydot[86] = beta0r*y[87] - alfa0r*y[86]; # z[22] = CK1r
            ydot[87] =  alfa0r*y[86] + kbr*y[88] - (beta0r + kfr)*y[87]; #z[23] = CK2r
            ydot[88] = kfr*y[87]+ psir*y[89] + beta1r*OKr - (kbr + alfai3r + alfa1r)*y[88]; # z[24] = CK3r
            ydot[89] = alfair*OKr + alfai3r*y[88] - (betair + psir)*y[89];# z[25] = IKr
            
            I_kr = self.gkr*OKr*(v-ek); # [uA/uF]    
            
        else:  #Model HH
            
            self.gkr = 0.078*(1 - self.BLOCKIKr)*np.sqrt(self.Ko/5.4); #0.03 por 0.078 Bondarenko
            xrss = 1/(1+exp(-(v+50)/7.5));
            tauxr = 1/(1.38e-3*(v+7)/(1-exp(-0.123*(v+7)))+6.1e-4*(v+10)/(exp(0.145*(v+10))-1));
            rkr = 1/(1+exp((v+33)/22.4));
            
            ydot[11] = (xrss-y[11])/tauxr;
            
            I_kr = self.gkr*y[11]*rkr*(v-ek);
        


        # I_ks: Slowly Activating K Current
        pcaks_junc = -np.log10(y[35])+3.0; 
        pcaks_sl = -np.log10(y[36])+3.0;  
        gks_junc = 0.07*(1 - self.BLOCKIKs)*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
        gks_sl = 0.07*(1 - self.BLOCKIKs)*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6))); 
        eks = (1/self.FoRT)*np.log((self.Ko+self.pNaK*self.Nao)/(y[34]+self.pNaK*y[33]));    
        xsss = 1/(1+exp(-(v-1.5)/16.7));
        tauxs = 1/(7.19e-5*(v+30)/(1-exp(-0.148*(v+30)))+1.31e-4*(v+30)/(exp(0.0687*(v+30))-1)); 
        
        ydot[12] = (xsss-y[12])/tauxs;
        
        if self.Model_IKs==1: # Model Winslow
        
            alfas = 7.956e-3;
            betas = 2.16e-1*exp(-2e-5*v);
            gamas = 3.97e-2;
            deltas = 7e-3*exp(-0.15*v);
            epslons = 7.67e-3*exp(0.087*v);
            omegas = 3.80e-3*exp(-0.014*v);
            
            C0Ks = y[90];
            C1Ks = y[91];
            O1Ks = y[92];
            O2Ks = 1 - (C0Ks + C1Ks + O1Ks);
            
            ydot[90] = betas*C1Ks - alfas*C0Ks;
            ydot[91] = alfas*C0Ks + deltas*O1Ks - (betas + gamas)*C1Ks; 
            ydot[92] = gamas*C1Ks + omegas*O2Ks - (deltas + epslons)*O1Ks;
            
            I_ks_junc = self.Fjunc*gks_junc*y[12]**2*(v-eks);
            I_ks_sl = self.Fsl*self.GKs*(O1Ks + O2Ks)*(v-eks); #[uA/uF]
            I_ks = I_ks_junc + I_ks_sl;
            
        elif self.Model_IKs==2: # Model Severi
            
            P1alfa = 8.7743e-3; 
            P2alfa = 0.0204; 
            P3alfa = 1.0349; 
            P1beta = 14.5254e-3;
            P2beta = 0.0399; 
            P3beta = 1.4191; 
            P1gama = 747.771e-3; 
            P2gama = -0.1112; 
            P3gama = 0.3418;
            P1delta = 93.2953e-3; 
            P2delta = -0.0695; 
            P1zeta = 55.2459e-3; 
            P2zeta = -0.0481; 
            P3zeta = 0.714;
            P4zeta = 3.8883e-3; 
            P1psi = 0.3215e-3;
            P2psi = 0.2529; 
            P1omega = 1.5922e-3; 
            P2omega =  -0.8405;
            teta = 3.0763e-3;
             
             
            alfa = P1alfa/(1+exp(-(v-P2alfa)/P3alfa*self.FoRT));       
            beta = P1beta/(1+exp((v-P2beta)/P3beta*self.FoRT));       
            gama = P1gama/(1+exp(-(v-P2gama)/P3gama*self.FoRT));       
            delta = P1delta*exp(P2delta*v*self.FoRT);        
            zeta = (P1zeta-P4zeta)/(1+exp((v-P2zeta)/P3zeta*self.FoRT))+ P4zeta;       
            psi = P1psi*exp(P2psi*v*self.FoRT);       
            omega = P1omega*exp(P2omega*v*self.FoRT);
                  
             
            gks = (1 - self.BLOCKIKs)*self.fatorgks*(1 + (0.6/(1 + (3.8e-5/y[36])**1.4)));
             
            C1 = y[93]; C2 = y[94]; C3 = y[95]; C4 = y[96]; C5 = y[97]; 
             
            C6 = y[98]; C7 = y[99]; C8 = y[100]; C9 = y[101]; C10 = y[102];
             
            C11 = y[103]; C12 = y[104]; C13 = y[105]; C14 = y[106]; C15 = y[107]; O1 = y[108]; 
            
            O2 = 1 - (O1 + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 + C11 + C12 + C13 + C14 + C15);
             
            ydot[93] = beta*C2 - 4*alfa*C1;
            ydot[94] = 2*beta*C3 + delta*C6 + 4*alfa*C1 - (beta + 3*alfa + gama)*C2;
            ydot[95] = 3*beta*C4 + delta*C7 + 3*alfa*C2 - (2*beta + 2*alfa + 2*gama)*C3;
            ydot[96] = 4*beta*C5 + delta*C8 + 2*alfa*C3 - (3*beta + alfa + 3*gama)*C4;
            ydot[97] = alfa*C4 + delta*C9 - (4*beta + 4*gama)*C5;
            ydot[98] = gama*C2 + beta*C7 - (delta + 3*alfa)*C6;
            ydot[99] = 3*alfa*C6 + 2*delta*C10 + 2*beta*C8 + 2*gama*C3 - (delta + beta + 2*alfa + gama)*C7;
            ydot[100] = 2*alfa*C7 + 2*delta*C11 + 3*beta*C9 + 3*gama*C4 - (delta + 2*beta + alfa + 2*gama)*C8;
            ydot[101] = 4*gama*C5 + 2*delta*C12 + alfa*C8 - (delta + 3*beta + 3*gama)*C9;
            ydot[102] = gama*C7 + beta*C11 - (2*delta + 2*alfa)*C10;
            ydot[103] = 2*alfa*C10 + 3*delta*C13 + 2*beta*C12 + 2*gama*C8 - (2*delta + beta + alfa + gama)*C11;
            ydot[104] = 3*gama*C9 + 3*delta*C14 + alfa*C11 - (2*delta + 2*beta + 2*gama)*C12;
            ydot[105] = gama*C11 + beta*C14 - (3*delta + alfa)*C13;
            ydot[106] = 2*gama*C12 + 4*delta*C15 + alfa*C13 - (3*delta + beta + gama)*C14;
            ydot[107] = gama*C14 + zeta*O1 - (4*delta + teta)*C15;
            ydot[108] = teta*C15 + omega*O2 - (zeta + psi)*O1;
             
            eks = (1/self.FoRT)*np.log((self.Ko+self.pNaK*y[32])/(y[34]+self.pNaK*y[32]));    
             
            I_ks_junc = self.Fjunc*gks_junc*y[12]**2*(v-eks);
            I_ks_sl = self.Fsl*gks*(O1 + O2)*(v-eks); #[uA/uF]
            I_ks = I_ks_junc + I_ks_sl;
  
        else: # Model HH    
            I_ks_junc = self.Fjunc*gks_junc*y[12]**2*(v-eks);
            I_ks_sl = self.Fsl*gks_sl*y[12]**2*(v-eks);
            I_ks = I_ks_junc+I_ks_sl;



        #I_kp: Plateau K current
        kp_kp = 1/(1+exp(7.488-v/5.98));
        I_kp_junc = self.Fjunc*self.gkp*kp_kp*(v-ek);
        I_kp_sl = self.Fsl*self.gkp*kp_kp*(v-ek);
        I_kp = I_kp_junc+I_kp_sl;
        
        
        ## IKur
        tauaur = 0.493*exp(-0.0629*v) + 2.058;
        tauiur = 1200 - (170/(1 + exp((v+45.2)/5.7)));
        ass = 1/(1 + exp(-(v+22.5)/7.7));
        iss = 1/(1 + exp((v+45.2)/5.7));
        aur = y[109];
        iur = y[110];
        
        ydot[109] = (ass - aur)/tauaur;
        ydot[110] = (iss - iur)/tauiur;
        
        IKur = self.GKur*aur*iur*(v - ek);
        
        ###IKss
        taukss = 39.3*exp(-0.0862*v) + 13.17;
        akss = y[111];
        ikss = y[112];
        
        ydot[111] = (ass - akss)/taukss;
        ydot[112] = 0;
        
        IKss = self.GKss*akss*ikss*(v-ek);
        
        ## I_to: Transient Outward K Current (slow and fast components) - Model HH
        xtoss = 1/(1+exp(-(v+3.0)/15));
        ytoss = 1/(1+exp((v+33.5)/10));
        rtoss = 1/(1+exp((v+33.5)/10));
        tauxtos = 9/(1+exp((v+3.0)/15))+0.5;
        tauytos = 3e3/(1+exp((v+60.0)/10))+30;
        taurtos = 2.8e3/(1+exp((v+60.0)/10))+220;                                   #Fei changed here!! time-dependent gating variable
        
        ydot[7] = (xtoss-y[7])/tauxtos;
        ydot[8] = (ytoss-y[8])/tauytos;
        ydot[39]= (rtoss-y[39])/taurtos;                                            #Fei changed here!! time-dependent gating variable
        
        tauxtof = 3.5*exp(-v*v/30/30)+1.5;
        tauytof = 20.0/(1+exp((v+33.5)/10))+20.0;
        
        ydot[9] = (xtoss-y[9])/tauxtof;
        ydot[10] = (ytoss-y[10])/tauytof;
        
        ## I_to: Transient Outward K Current (slow and fast components) - Model Winslow
        
        #Itof
        alfa_a = 0.675*exp(0.0255*v);
        beta_a = 0.088269*exp(-0.0883*v);
        alfa_i = 0.109566;
        beta_i = 3.03334e-4;
        f1 = 1.66120;
        f2 = 22.2463;
        f3 = 195.978;
        f4 = 181.609;
        b1 = 0.72246;
        b2 = 0.47656;
        b3 = 7.77537;
        b4 = 318.232;
        C0f = y[56];
        C1f = y[57];
        C2f = y[58];
        C3f = y[59];
        CI0f = y[60];
        CI1f = y[61];
        CI2f = y[62];
        CI3f = y[63];
        If = y[64];
        Of = 1 - (y[56]+ y[57] + y[58] + y[59] + y[60] + y[61]+ y[62]+ y[63]+ y[64]);
        
        #Itos
        alfa_as = 1.840024*exp(0.0077*v);
        beta_as = 0.010817*exp(-0.0779*v);
        alfa_is = 0.003058;
        beta_is = 2.4936e-6;
        f1s = 0.52465;
        f2s = 17.5188;
        f3s = 938.587;
        f4s = 54749.1;
        b1s = 1.00947;
        b2s = 1.17100;
        b3s = 0.63902;
        b4s = 2.12035;
        C0s = y[65];
        C1s = y[66];
        C2s = y[67];
        C3s = y[68];
        CI0s = y[69];
        CI1s = y[70];
        CI2s = y[71];
        CI3s = y[72];
        Is = y[73];
        Os = 1 - (y[65]+ y[66] + y[67] + y[68] + y[69] + y[70]+ y[71]+ y[72]+ y[73]);
        Ki = y[34];
        Nass = y[32];
        
        ydot[56] = alfa_i*CI0f + beta_a*C1f - (4*alfa_a + beta_i)*C0f;
        ydot[57] = 4*alfa_a*C0f + 2*beta_a*C2f + (alfa_i/b1)*CI1f - (beta_a + 3*alfa_a + f1*beta_i)*C1f; 
        ydot[58] = 3*alfa_a*C1f + 3*beta_a*C3f + (alfa_i/b2)*CI2f - (2*beta_a + 2*alfa_a + f2*beta_i)*C2f;
        ydot[59] = 2*alfa_a*C2f + 4*beta_a*Of + (alfa_i/b3)*CI3f - (3*beta_a + alfa_a + f3*beta_i)*C3f;
        ydot[60] = beta_i*C0f + (beta_a/f1)*CI1f - (alfa_i + b1*4*alfa_a)*CI0f;
        ydot[61] = b1*4*alfa_a*CI0f + (f1*2*beta_a/f2)*CI2f + f1*beta_i*C1f - ((beta_a/f1)+ (alfa_i/b1) + (b2*3*alfa_a/b1))*CI1f;
        ydot[62] = (b2*3*alfa_a/b1)*CI1f + (f2*3*beta_a/f3)*CI3f + f2*beta_i*C2f -((f1*2*beta_a/f2)+ (alfa_i/b2) + (b3*2*alfa_a/b2))*CI2f;
        ydot[63] = (b3*2*alfa_a/b2)*CI2f + (f3*4*beta_a/f4)*If + f3*beta_i*C3f - ((f2*3*beta_a/f3)+ (alfa_i/b3) + (b4*alfa_a/b3))*CI3f;
        ydot[64] = f4*beta_i*Of + (b4*alfa_a/b3)*CI3f - ((alfa_i/b4) + (f3*4*beta_a/f4))*If;
        ydot[65] = alfa_is*CI0s + beta_as*C1s - (4*alfa_as + beta_is)*C0s;
        ydot[66] = 4*alfa_as*C0s + 2*beta_as*C2s + (alfa_is/b1s)*CI1s - (beta_as + 3*alfa_as + f1s*beta_is)*C1s; 
        ydot[67] = 3*alfa_as*C1s + 3*beta_as*C3s + (alfa_is/b2s)*CI2s - (2*beta_as + 2*alfa_as + f2s*beta_is)*C2s;
        ydot[68] = 2*alfa_as*C2s + 4*beta_as*Os + (alfa_is/b3s)*CI3s - (3*beta_as + alfa_as + f3s*beta_is)*C3s;
        ydot[69] = beta_is*C0s + (beta_as/f1s)*CI1s - (alfa_is + b1s*4*alfa_as)*CI0s;
        ydot[70] = b1s*4*alfa_as*CI0s + (f1s*2*beta_as/f2s)*CI2s + f1s*beta_is*C1s - ((beta_as/f1s)+ (alfa_is/b1s) + (b2s*3*alfa_as/b1s))*CI1s;
        ydot[71] = (b2s*3*alfa_as/b1s)*CI1s + (f2s*3*beta_as/f3s)*CI3s + f2s*beta_is*C2s -((f1s*2*beta_as/f2s)+ (alfa_is/b2s) + (b3s*2*alfa_as/b2s))*CI2s;
        ydot[72] = (b3s*2*alfa_as/b2s)*CI2s + (f3s*4*beta_as/f4s)*Is + f3s*beta_is*C3s - ((f2s*3*beta_as/f3s)+ (alfa_is/b3s) + (b4s*alfa_as/b3s))*CI3s;
        ydot[73] = f4s*beta_is*Os + (b4s*alfa_as/b3s)*CI3s - ((alfa_is/b4s) + (f3s*4*beta_as/f4s))*Is;


        if self.Model_Ito == 1: # Model Winslow

            I_tof = self.GItof*Of*(v-ek); #[uA/uF]   
            I_tos = self.GItos*Os*(4*v*self.Frdy*self.FoRT)*((Ki*exp(v*self.FoRT) - self.Ko) /(exp(v*self.FoRT)-1) + 0.02*(Nass*exp(v*self.FoRT) - self.Nao) /(exp(v*self.FoRT)-1)); 
            #Itos = ito*(10**7)*SAsl/Cmem;

        else:  # Model HH

            I_tos = self.GtoSlow*y[7]*(y[8]+0.5*y[39])*(v-ek);                               # [uA/uF]
            I_tof = self.GtoFast*y[9]*y[10]*(v-ek);



        I_to = I_tos + I_tof;
        
        # I_ki: Time-Independent K Current # ALTERAR COMO NO PAPER DO BONDARENKO
        aki = 1.02/(1+exp(0.2385*(v-ek-59.215)));
        bki =(0.49124*exp(0.08032*(v+5.476-ek)) + exp(0.06175*(v-ek-594.31))) /(1 + exp(-0.5143*(v-ek+4.753)));
        self.kiss = aki/(aki+bki);
        
        I_ki = 0.9*np.sqrt(self.Ko/5.4)*self.kiss*(v-ek);
        
        # I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
        I_ClCa_junc = self.Fjunc*GClCa/(1+KdClCa/y[35])*(v-ecl);
        I_ClCa_sl = self.Fsl*GClCa/(1+KdClCa/y[36])*(v-ecl);
        I_ClCa = I_ClCa_junc+I_ClCa_sl;
        I_Clbk = GClB*(v-ecl);
        
        ## I_Ca L-type Current - Model Mahajan
        Pc2_LCCj_m1=y[40]; 
        Pc1_LCCj_m1=y[41]; 
        Pi1Ca_LCCj_m1=y[42];
        Pi2Ca_LCCj_m1=y[43]; 
        Pi1Ba_LCCj_m1=y[44]; 
        Pi2Ba_LCCj_m1=y[45];
        
        cajLCC = y[35]; 
        ICa_speed = 1;
        
        # LTCC Current - Fixed Parameters
        cpt = 3.75e-3;                                                      # [mM]
        cat = 7.617e-3;                                                     # [mM]
        s1o = 0.0182688;                                                    # [1/ms]
        k1o = 0.024168;                                                     # [1/ms]
        k2o = 0.000103615;                                                  # [1/ms]
        sp0 = 1.5;
        sp1 = 3;                                                            # [ms]
        sp2 = 40;                                                           # [mV]
        sp3 = 3;                                                            # [mV]
        sp4 = 4;                                                            # [mV]
        sp5 = 11.32;                                                        # [mV]
        sp6 = 15.6;                                                         # [mV]
        sp7 = 10;                                                           # [ms]
        sp8 = 4954;                                                         # [ms]
        sp9 = 78.0329;                                                      # [ms]
        sp10 = 0.1;                                                         # [ms]
        aR2 = 1;
        sR2 = -2;                                                           # [mV]
        pR2 = 0.145;                                                        # [1/mV]
        aT2 = 1;                                                            # [1/ms]
        sT2 = -1000;                                                        # [mV]
        pT2 = 0.100;                                                        # [1/mV]
        aR1 = 0.09091;     
        sR1 = -1000;                                                        # [mV]
        pR1 = 0.100;                                                        # [1/mV]
        aT1 = 0.30303;                                                      # [1/ms]
        sT1 = -1000;                                                        # [mV]
        pT1 = 0.100;                                                        # [1/mV]
        aRv2 = 0.9;
        sRv2 = -29;                                                         # [mV]
        pRv2 = 0.135;                                                       # [1/mV]
        aTv2 = 500;                                                         # [1/ms]
        sTv2 = -25;                                                         # [mV]
        pTv2 = 0.050;                                                       # [1/mV]
        aRv1 = 0.85;
        sRv1 = -180;                                                        # [mV]
        pRv1 = 0.090;                                                       # [1/mV]
        aTv1 = 270;                                                         # [1/ms]
        sTv1 = -180;                                                        # [mV]
        pTv1 = 0.100;                                                       # [1/mV]
        aTrev1 = 205.12;                                                    # [1/ms]
        sTrev1 = -65;                                                       # [mV]
        pTrev1 = 0.100;                                                     # [1/mV]
        aTrev2 = 7e8;                                                       # [1/ms]
        sTrev2 = 60;                                                        # [mV]
        pTrev2 = 0.130;                                                     # [1/mV]
        
        # Voltage and Ca-dependent Parameters
        fcp=1/(1+(cpt/cajLCC/aff)**3);                                       # Ca-dep
        tca=sp9/(1+(cajLCC*aff/cat)**4)+sp10;                                # Ca-dep
        R2=aR2/(1+exp(-(v-sR2)*pR2));
        T2=aT2/(1+exp(-(v-sT2)*pT2));
        PT=1-(1/(1+exp(-(v+sp2)/sp3)));
        R1=aR1/(1+exp(-(v-sR1)*pR1));
        T1=aT1/(1+exp(-(v-sT1)*pT1));
        RV=sp7+sp8*exp(v/sp6);
        Pr=1-(1/(1+exp(-(v+sp2)/sp4)));
        Pq=1+sp0/(1+exp(-(v+sp2)/sp4));
        TCa=Pq*((RV-tca)*Pr+tca);
        Ps=1/(1+exp(-(v+sp2)/sp5));
        Rv1=aRv1/(1+exp(-(v-sRv1)*pRv1));
        Tv1=aTv1/(1+exp(-(v-sTv1)*pTv1));
        Rv2=aRv2/(1+exp(-(v-sRv2)*pRv2));
        Tv2=aTv2/(1+exp(-(v-sTv2)*pTv2));
        Trev1=aTrev1/(1+exp(-(v-sTrev1)*pTrev1));
        Frev1=(1-Rv1)/Rv1*R1/(1-R1);
        Trev2=aTrev2/(1+exp(-(v-sTrev2)*pTrev2));
        Frev2=(1-Rv2)/Rv2*R2/(1-R2)*Rv1/(1-Rv1);
        
        # Transition Rates (20 rates)
        alphaLCC=ICa_speed*R2/T2;
        betaLCC=ICa_speed*(1-R2)/T2;
        r1=ICa_speed*R1/T1;
        r2=ICa_speed*(1-R1)/T1;
        k1=ICa_speed*k1o*fcp;
        k2=ICa_speed*k2o;
        k3=ICa_speed*PT/sp1;
        k5=ICa_speed*(1-Ps)/TCa;
        k6=ICa_speed*fcp*Ps/TCa;
        s1=ICa_speed*s1o*fcp;
        k1p=ICa_speed*Rv1/Tv1;
        k2p=ICa_speed*(1-Rv1)/Tv1;
        k3p=ICa_speed*1/(Trev2*(1+Frev2));
        k4p=Frev2*k3p;                            
        k5p=ICa_speed*(1-Rv2)/Tv2;
        k6p=ICa_speed*Rv2/Tv2;
        s1p=ICa_speed*1/(Trev1*(1+Frev1));
        s2p=Frev1*s1p;                            
        k4=k3*(alphaLCC/betaLCC)*(k1/k2)*(k5/k6); 
        s2=s1*(k2/k1)*(r1/r2);                    
            
        Po = 1 - y[40] - y[41] - y[42] - y[43] - y[44] - y[45];
        Po_LCCj_m1 = Po;
        
        # State transitions for mode-1 junctional LCCs
        ydot[40] = betaLCC*Pc1_LCCj_m1 + k5*Pi2Ca_LCCj_m1 + k5p*Pi2Ba_LCCj_m1 - (k6+k6p+alphaLCC)*Pc2_LCCj_m1;                      # C2_m1j 
        ydot[41] = alphaLCC*Pc2_LCCj_m1 + k2*Pi1Ca_LCCj_m1 + k2p*Pi1Ba_LCCj_m1 + r2*Po_LCCj_m1 - (r1+betaLCC+k1+k1p)*Pc1_LCCj_m1;   # C1_m1j 
        ydot[42] = k1*Pc1_LCCj_m1 + k4*Pi2Ca_LCCj_m1 + s1*Po_LCCj_m1 - (k2+k3+s2)*Pi1Ca_LCCj_m1;                                    # I1Ca_m1j  
        ydot[43] = k3*Pi1Ca_LCCj_m1 + k6*Pc2_LCCj_m1 - (k4+k5)*Pi2Ca_LCCj_m1;                                                       # I2Ca_m1j
        ydot[44] = k1p*Pc1_LCCj_m1 + k4p*Pi2Ba_LCCj_m1 + s1p*Po_LCCj_m1 - (k2p+k3p+s2p)*Pi1Ba_LCCj_m1;                              # I1Ba_m1j
        ydot[45] = k3p*Pi1Ba_LCCj_m1 + k6p*Pc2_LCCj_m1 - (k5p+k4p)*Pi2Ba_LCCj_m1;                                                   # I2Ba_m1j
        
        ## I_Ca: L-type Calcium Current - Model HH
        dss = 1/(1+exp(-(v+14.5)/6.0));
        taud = dss*(1-exp(-(v+14.5)/6.0))/(0.035*(v+14.5));
        fss = 1/(1+exp((v+35.06)/3.6))+0.6/(1+exp((50-v)/20));
        tauf = 1/(0.0197*exp( -(0.0337*(v+14.5))**2 )+0.02);
        
        ydot[3] = (dss-y[3])/taud;
        ydot[4] = (fss-y[4])/tauf;
        ydot[5] = 1.7*y[35]*(1-y[5])-11.9e-3*y[5]; 
        ydot[6] = 1.7*y[36]*(1-y[6])-11.9e-3*y[6]; 
                              
        fcaCaMSL=0;
        fcaCaj= 0;


        if self.Model_ICaL == 1: # Model Mahajan

            ibarca_j = self.PCam*4*(v*self.Frdy*self.FoRT)*(self.gCai*y[35]*exp(2*v*self.FoRT)-self.gCao*self.Cao) /(exp(2*v*self.FoRT)-1);
            ibarca_sl = self.PCam*4*(v*self.Frdy*self.FoRT)*(self.gCai*y[36]*exp(2*v*self.FoRT)-self.gCao*self.Cao) /(exp(2*v*self.FoRT)-1);
            I_Ca_junc = self.Fjunc_CaL*ibarca_j*Po*self.ICa_scale;
            I_Ca_sl = self.Fsl_CaL*ibarca_sl*Po*self.ICa_scale;

        else: # Model HH
            ibarca_j = self.pCa*4*(v*self.Frdy*self.FoRT)*(self.gCai*y[35]*exp(2*v*self.FoRT)-self.gCao*self.Cao) /(exp(2*v*self.FoRT)-1);
            ibarca_sl = self.pCa*4*(v*self.Frdy*self.FoRT)*(self.gCai*y[36]*exp(2*v*self.FoRT)-self.gCao*self.Cao) /(exp(2*v*self.FoRT)-1);
            I_Ca_junc = (self.Fjunc_CaL*ibarca_j*y[3]*y[4]*((1-y[5])+fcaCaj)*self.Q10CaL**self.Qpow)*0.45*1;
            I_Ca_sl = (self.Fsl_CaL*ibarca_sl*y[3]*y[4]*((1-y[6])+fcaCaMSL)*self.Q10CaL**self.Qpow)*0.45*1;



        I_Ca = I_Ca_junc + I_Ca_sl;

        #Correntes de Sodio e Potassio - ativadas por Calcio
        ibark = pK*(v*self.Frdy*self.FoRT)*(self.gKi*y[34]*exp(v*self.FoRT)-self.gKo*self.Ko) /(exp(v*self.FoRT)-1);
        ibarna_j = self.pNa*(v*self.Frdy*self.FoRT) *(self.gNai*y[31]*exp(v*self.FoRT)-self.gNao*self.Nao)/(exp(v*self.FoRT)-1);
        ibarna_sl = self.pNa*(v*self.Frdy*self.FoRT) *(self.gNai*y[32]*exp(v*self.FoRT)-self.gNao*self.Nao)/(exp(v*self.FoRT)-1);
        
        I_CaK = (ibark*y[3]*y[4]*(self.Fjunc_CaL*(fcaCaj+(1-y[5]))+self.Fsl_CaL*(fcaCaMSL+(1-y[6])))*self.Q10CaL**self.Qpow)*0.45*1;
        I_CaNa_junc =(self.Fjunc*ibarna_j*y[3]*y[4]*((1-y[5])+fcaCaj)*self.Q10CaL**self.Qpow)*0.45*1;
        I_CaNa_sl = (self.Fsl*ibarna_sl*y[3]*y[4]*((1-y[6])+fcaCaMSL)*self.Q10CaL**self.Qpow)*.45*1;
        I_CaNa = I_CaNa_junc+I_CaNa_sl;
        I_Catot = I_Ca+I_CaK+I_CaNa; ##ok<NASGU>
        
        # I_ncx: Na/Ca Exchanger flux
        Ka_junc = 1/(1+(self.Kdact/y[35])**3);
        Ka_sl = 1/(1+(self.Kdact/y[36])**3);
        s1_junc = exp(self.nu*v*self.FoRT)*y[31]**3*self.Cao;
        s1_sl = exp(self.nu*v*self.FoRT)*y[32]**3*self.Cao;
        s2_junc = exp((self.nu-1)*v*self.FoRT)*self.Nao**3*y[35];
        s3_junc = self.KmCai*self.Nao**3*(1+(y[31]/self.KmNai)**3) + self.KmNao**3*y[35]*(1+y[35]/self.KmCai)+self.KmCao*y[31]**3+y[31]**3*self.Cao+self.Nao**3*y[35];
        s2_sl = exp((self.nu-1)*v*self.FoRT)*self.Nao**3*y[36];
        s3_sl = self.KmCai*self.Nao**3*(1+(y[32]/self.KmNai)**3) + self.KmNao**3*y[36]*(1+y[36]/self.KmCai)+self.KmCao*y[32]**3+y[32]**3*self.Cao+self.Nao**3*y[36];
        
        I_ncx_junc = self.Fjunc*self.IbarNCX*self.Q10NCX**self.Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+self.ksat*exp((self.nu-1)*v*self.FoRT));
        I_ncx_sl = self.Fsl*self.IbarNCX*self.Q10NCX**self.Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+self.ksat*exp((self.nu-1)*v*self.FoRT));
        I_ncx = I_ncx_junc+I_ncx_sl; ##ok<NASGU>
        
        # I_pca: Sarcolemmal Ca Pump Current
        I_pca_junc = self.Fjunc*self.Q10SLCaP**self.Qpow*self.IbarSLCaP*y[35]**2/(self.KmPCa**2+y[35]**2);
        I_pca_sl = self.Fsl*self.Q10SLCaP**self.Qpow*self.IbarSLCaP*y[36]**2/(self.KmPCa**2+y[36]**2);
        I_pca = I_pca_junc+I_pca_sl; ##ok<NASGU>
        
        # I_cabk: Ca Background Current
        I_cabk_junc = self.Fjunc*self.GCaB*(v-eca_junc);
        I_cabk_sl = self.Fsl*self.GCaB*(v-eca_sl);
        I_cabk = I_cabk_junc+I_cabk_sl; ##ok<NASGU>
        
        ## SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
        MaxSR = 15; 
        MinSR = 1;
        kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y[30])**2.5);
        koSRCa = koCa/kCaSR;
        kiSRCa = kiCa*kCaSR;
        
        RI = 1-y[13]-y[14]-y[15];
        
        ydot[13] = (kim*RI-kiSRCa*y[35]*y[13])-(koSRCa*y[35]**2*y[13]-kom*y[14]);   # R
        ydot[14] = (koSRCa*y[35]**2*y[13]-kom*y[14])-(kiSRCa*y[35]*y[14]-kim*y[15]);# O
        ydot[15] = (kiSRCa*y[35]*y[14]-kim*y[15])-(kom*y[15]-koSRCa*y[35]**2*RI);   # I
        
        J_SRCarel = self.ks*(1 - self.BLOCKCICR)*(y[14])*(y[30]-y[35]);                                       # [mM/ms]
        
        J_serca = self.Q10SRCaP**self.Qpow*self.Vmax_SRCaP*((y[37]/self.Kmf)**self.hillSRCaP-(y[30]/self.Kmr)**self.hillSRCaP)/(1+(y[37]/self.Kmf)**self.hillSRCaP+(y[30]/self.Kmr)**self.hillSRCaP);
        
        J_SRleak = self.GSrLeak*(y[30]-y[37]);                                           #  [mM/ms]
        
        ## Forca de Contracao - Modelo de Negroni & Lascano - 2008

        if self.Model_Force == 0:
            hw = y[74];
            hp = y[75];
            TSCa = y[76];   #complexo TS + 3Ca
            TSCaw = y[77];  #complexo TS + 3Ca + miosina no estado fraco
            TSCap = y[78];  #complexo TS + 3Ca + miosina no estado forte
            TSp = y[79];    #complexo TS + miosina no estado forte
            TSw = y[80];    #complexo TS + miosina no estado fraco
            
            self.TSt = 140e-3;    #concentracao total de troponina em mM
            
            TS = self.TSt - TSCa - TSCaw - TSCap - TSw - TSp;
            
            Rneg = 15; #um**-2
            La = 1.15; #um
            
            TSCaeff = TSCa*exp(-Rneg*(self.Lsarc - La)*(self.Lsarc-La));
            
            hpr = 0.006; #um
            hwr = 0.0001; #um
            B = 0.5; #ms**-1
            Za = 0.0023; #ms**-1
            Yv = 1.5; #ms**-1
            gama = 28000; #um**2
            Yd = 0.0333; # ms**-1
            Yc = 1; #ms**-1.um**2
            Lc = 1.2; #um
            Yvd = 1.5; #ms**-1
            Yb = 0.1816e9; #(mM)**-3.ms**-1   CONFERIR
            Yp = 0.1397; # ms**-1
            Yr = 0.1397;# ms**-1
            Zb = 0.1397;# ms**-1
            Zp = 0.2095;# ms**-1
            Yq = 0.2328;# ms**-1
            Zq = 0.3724;# ms**-1
            fneg = 0.0023; # ms**-1
            Zr = 7.2626e9; #(mM)**-3.ms**-1   CONFERIR
            
            g = Za + Yv*(1 - exp(-gama*(hw - hwr)*(hw - hwr)));
            gd = Yd + Yc*(self.Lsarc - Lc)**2 + Yvd*(1 - exp(-gama*(hw - hwr)*(hw - hwr)));
            
            ydot[74] = -B*(hw - hwr);
            ydot[75] = -B*(hp - hpr);
            ydot[76]= TS*(y[37]**3)*Yb + g*TSCaw - Zb*TSCa - fneg*TSCaeff; #[TSCa3]
            ydot[77] = fneg*TSCaeff + Zp*TSCap - g*TSCaw - Yp*TSCaw;      #[TSCa3~]
            ydot[78] = Yp*TSCaw + Zr*TSp*(y[37]**3) - Zp*TSCap - TSCap*Yr; #[TSCa3*]
            ydot[79] = Yr*TSCap + Zq*TSw - Zr*TSp*(y[37]**3) - Yq*TSp;     #[TS*]
            ydot[80] = Yq*TSp - Zq*TSw - gd*TSw;                          #[TS~]
            
            # Forca de contracao
            
            self.Kp = 2700e3; #mN.mm**2.um**-1.mM**-1
            self.Kw = 540e3; #mN.mm**2.um**-1.mM**-1
            self.Lo = 0.97; #um
            self.Ke = 105000; #mN.mm**-2.um**-5
            self.Le = 10; #mN.mm**-2.um**-1
            
            Fp = self.Ke*(self.Lsarc - self.Lo)**5 + self.Le*(self.Lsarc - self.Lo);
            Fb = self.Kw*(TSCaw + TSw)*hw + self.Kp*(TSCap + TSp)*hp;
            FORCA = Fp + Fb;
            
            #Modelo para encurtamento do sarcomero
            self.alfaneg = 0.5;
            self.betaneg = 80;
            Ls = (1/self.betaneg)*np.log((FORCA/self.alfaneg)+1);
            self.fit = Ls;
            #Lsim = Lsarc - (Ls - fit); #comprimento do sarcomero, sendo Lsarc o comprimento inicial
            #CelL = cellLength *(Lsim/Lsarc);
            
            # Buffer Troponina
            Btrop = 3*(ydot[76] + ydot[77] + ydot[78]);

        else:
            ## Forca de contracao - Modelo de Rice et al (1999)
            
            #Parametros
            
            kPN = 0.045;    #ms**-1 taxa de transicao do estado permissivo para nao permissivo
            fXB = 0.10;     #ms**-1 taxa de transicao da ponte cruzada no estado fraco para o estado forte.
            gminxb = 0.14;  #ms**-1 taxa minima de desligamento da ponte cruzada forte para fraca
            f01 = 3*fXB;
            f12 = 10*fXB;
            f23 = 7*fXB;
            SLnorm = self.Lsarc - self.Lajust; #verificar
            if SLnorm < 0:
                gxbSL = gminxb*(2 - math.pow(SLnorm*-1,1.6));
            else:
                gxbSL = gminxb*(2 - math.pow(SLnorm,1.6));
            g10SL = gxbSL;
            g21SL = 2*gxbSL;
            g32SL = 3*gxbSL;
            KCa = koff_tncl/kon_tncl;
            Khalf = 1/(1 + ((KCa)/(1.5e-3 - SLnorm*1)));
            Ntm = 5 + 3*SLnorm;
            kNP = kPN*(y[18]/(Bmax_TnClow*Khalf))**Ntm;
            
            N0 = y[81];
            N1 = y[82];
            P0 = y[83];
            P1 = y[84];
            P2 = y[85];
            
            P3 = 1 - (N0 + N1 + P0 + P1 + P2);
            
            ydot[81] = kPN*P0 + g10SL*N1 - kNP*N0;
            ydot[82] = kPN*P1 -(kNP + g10SL)*N1;
            ydot[83] = g10SL*P1 + kNP*N0 - (kPN + f01)*P0;
            ydot[84] = g21SL*P2 + f01*P0 + kNP*N1 - (f12 + g10SL + kPN)*P1;
            ydot[85] = g32SL*P3 + f12*P1 - (f23 + g21SL)*P2;
            
            #Calculo da forca de contracao
            
            SOMA = gminxb*2*gminxb*3*gminxb + f01*2*gminxb*3*gminxb + f01*f12*3*gminxb + f01*f12*f23;
            
            P3max = (f01*f12*f23)/SOMA;
            P2max = (f01*f12*3*gminxb)/SOMA;
            P1max = (f01*2*gminxb*3*gminxb)/SOMA;
            self.Fmax = P1max + 2*P2max + 3*P3max;



        ## Sodium Buffering
        ydot[16] = kon_na*y[31]*(Bmax_Naj-y[16])-koff_na*y[16];                     # NaBj      [mM/ms]
        ydot[17] = kon_na*y[32]*(Bmax_Nasl-y[17])-koff_na*y[17];                    # NaBsl     [mM/ms]
        
        # Cytosolic Ca Buffers
        ydot[18] = kon_tncl*y[37]*(Bmax_TnClow-y[18])-koff_tncl*y[18];              # TnCL      [mM/ms]
        ydot[19] = kon_tnchca*y[37]*(Bmax_TnChigh-y[19]-y[20])-koff_tnchca*y[19];   # TnCHc     [mM/ms]
        ydot[20] = kon_tnchmg*self.Mgi*(Bmax_TnChigh-y[19]-y[20])-koff_tnchmg*y[20];     # TnCHm     [mM/ms]
        ydot[21] = kon_cam*y[37]*(Bmax_CaM-y[21])-koff_cam*y[21];                   # CaM       [mM/ms]
        ydot[22] = kon_myoca*y[37]*(Bmax_myosin-y[22]-y[23])-koff_myoca*y[22];      # Myosin_ca [mM/ms]
        ydot[23] = kon_myomg*self.Mgi*(Bmax_myosin-y[22]-y[23])-koff_myomg*y[23];        # Myosin_mg [mM/ms]
        ydot[24] = kon_sr*y[37]*(Bmax_SR-y[24])-koff_sr*y[24];                      # SRB       [mM/ms]

        if self.Model_Force == 0:
            J_CaB_cytosol = np.sum(ydot[20:24])- ydot[22]; # CONFERIR
        else:
            J_CaB_cytosol = np.sum(ydot[18:24]);


        # Junctional and SL Ca Buffers
        ydot[25] = kon_sll*y[35]*(Bmax_SLlowj-y[25])-koff_sll*y[25];       # SLLj      [mM/ms]
        ydot[26] = kon_sll*y[36]*(Bmax_SLlowsl-y[26])-koff_sll*y[26];      # SLLsl     [mM/ms]
        ydot[27] = kon_slh*y[35]*(Bmax_SLhighj-y[27])-koff_slh*y[27];      # SLHj      [mM/ms]
        ydot[28] = kon_slh*y[36]*(Bmax_SLhighsl-y[28])-koff_slh*y[28];     # SLHsl     [mM/ms]
        
        J_CaB_junction = ydot[25]+ydot[27];
        J_CaB_sl = ydot[26]+ydot[28];
        
        # SR - CSQ Buffers
        ydot[29] = kon_csqn*y[30]*(Bmax_Csqn-y[29])-koff_csqn*y[29];                # Csqn      [mM/ms]
        
        ## Ion concentrations
        
        # SR Ca Concentrations
        ydot[30] = J_serca-(J_SRleak*Vmyo/Vsr + J_SRCarel)- ydot[29];                                                   # Ca_sr em mM   
        
        # Sodium Concentrations
        I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3*I_ncx_junc + 3*I_nak_junc + I_CaNa_junc;                            # [uA/uF]
        I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3*I_ncx_sl + 3*I_nak_sl + I_CaNa_sl;                                        # [uA/uF]
        
        ydot[31] = -I_Na_tot_junc*self.Cmem/(self.Vjunc*self.Frdy) + self.J_na_juncsl/self.Vjunc*(y[32]-y[31])- ydot[16];                             #Na_j
        ydot[32] = -I_Na_tot_sl*self.Cmem/(Vsl*self.Frdy)+ self.J_na_juncsl/Vsl*(y[31]-y[32])+ J_na_slmyo/Vsl*(y[33]-y[32])- ydot[17];      #Na_sl
        ydot[33] = J_na_slmyo/Vmyo*(y[32]-y[33]);                                                                            #Na_cy
        
        # Potassium Concentration
        I_K_tot = I_to +I_kr +I_ks + I_ki - 2*I_nak + I_CaK + I_kp + IKur + IKss;                                                          # [uA/uF]
        
        ydot[34] = -I_K_tot*self.Cmem/(Vmyo*self.Frdy);                                                                               #K

        # Calcium Concentrations
        I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2*I_ncx_junc;                                                # [uA/uF]
        I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2*I_ncx_sl;            # [uA/uF]
        
        ydot[35] = -I_Ca_tot_junc*self.Cmem/(self.Vjunc*2*self.Frdy) + J_ca_juncsl/self.Vjunc*(y[36]-y[35])- J_CaB_junction + (J_SRCarel)*Vsr/self.Vjunc; # Ca_j
        ydot[36] = -I_Ca_tot_sl*self.Cmem/(Vsl*2*self.Frdy)+ J_ca_juncsl/Vsl*(y[35]-y[36])+ J_ca_slmyo/Vsl*(y[37]-y[36])- J_CaB_sl;        # Ca_sl
        
        if self.Model_Force == 0:
            ydot[37] = -J_serca*Vsr/Vmyo - J_CaB_cytosol + J_ca_slmyo/Vmyo*(y[36]-y[37]) + J_SRleak - Btrop;                                 # Ca_cy
        else:
            ydot[37] = -J_serca*Vsr/Vmyo - J_CaB_cytosol + J_ca_slmyo/Vmyo*(y[36]-y[37]) + J_SRleak;                                 # Ca_cy  
        
            
        # Membrane Potential
        I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;            # [uA/uF]
        I_Cl_tot = I_ClCa + I_Clbk;                        # [uA/uF]
        I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
        I_tot = I_Na_tot + I_Cl_tot +I_Ca_tot + I_K_tot;
        #if self.w < 30 and ((t > self.t_ap) and t < self.t_ap + self.D) and (np.mod((t + (1000 - self.t_ap)),(self.w + self.Delay)))<= self.w:                                
        #    print(I_tot);                   
        ydot[38] = -(I_tot-I_app);
        #print(ydot[38]);									
        return ydot;
#######################################################################################################################################################################
            
    def principal(self, protocol, Model_ICaL, Model_Na, Model_Ito, Model_IKr, Model_IKs, Model_Force, cellLength, Lsarc,
		BLOCKSRPUMP,STIMULUSRPUMP,BLOCKNCX,STIMULUSNCX,CAFEINA,BLOCKCICR,BLOCKIKs,BLOCKIKr,BLOCKItof,BLOCKItos,BLOCKINa,BLOCKICaL,STIMULUSICaL,
		t_ap,L,Ap,v_resting,tap,w,f,Delay,A_inj,tig):            
###################################################################################################################################################
# PROTOCOLO 1: - Potencial da Membrana gerado de forma artificial
#                com amplitude e frequencia controladas: w = 200, Ap = 15, v_resting = -80 

# PROTOCOLO 2: - Experimento de voltage clamp com voltagem controlada
#                e injecao de corrente dependente da voltagem aplicada: w = 200, Ap = 15, v_resting = -80 

# PROTOCOLO 3: - Corrente injetada (9.5 ou 10 uA/uF)com possibilidade de aplicacao de dois trens de pulso
#                com frequencia f e f2, largura w e w2 com duracao do primeiro trem de pulso igual a D ms

#        pasta = str(protocol)+str(Model_ICaL)+str(Model_Na)+str(Model_Ito)+str(Model_IKr)+str(Model_IKs)+str(Model_Force)+str(cellLength)+str(Lsarc)+str(BLOCKSRPUMP)+str(STIMULUSRPUMP)+str(BLOCKNCX)+str(STIMULUSNCX)+str(CAFEINA)+str(BLOCKCICR)+str(BLOCKIKs)+str(BLOCKIKr)+str(BLOCKItof)+str(BLOCKItos)+str(BLOCKINa)+str(BLOCKICaL)+str(STIMULUSICaL)+str(t_ap)+str(L)+str(Ap)+str(v_resting)+str(tap)+str(w)+str(f)+str(Delay)+str(A_inj)+str(tig)
		
        
#        del self # so apagas esta referencia para o objeto nao o objeto em si
        #self.limparClasse();		
        self.protocol = protocol;
        
        self.Model_ICaL = Model_ICaL;

        self.Model_Na = Model_Na;

        self.Model_Ito =  Model_Ito;

        self.Model_IKr = Model_IKr;

        self.Model_IKs =  Model_IKs

        self.Model_Force = Model_Force

        self.cellLength  = cellLength

        self.Lsarc = Lsarc
        

        self.alfaaj = -730.26;               # Parametro de ajuste da forca de contracao para o modelo de Rice: Fcontr[i] = alfaaj*Fcontrn[i];  
        self.betaaj = 0.8;                   # Parametro de ajuste para o encurtamento do sarcomero no modelo de Rice: SL[i] = betaaj*Fcontrn[i] + Lsarc;
        self.Lajust = 1.3;                   # Parametro de ajuste para calculo do SLnorm utilizado no kNP - modelo de Rice; SLnorm = (SL - Lajust)/(2.3 - Lajust)

        self.alfaneg = 0.5;                  # Parametro de ajuste para o encurtamento do sarcomero no modelo de Negroni Ls[i] = (1/betaneg)*np.log((FORCA[i]/alfaneg)+1);
        
        self.betaneg = 80;                   # Parametro de ajuste para o encurtamento do sarcomero no modelo de Negroni

        self.BLOCKSRPUMP = BLOCKSRPUMP
        self.STIMULUSRPUMP = STIMULUSRPUMP
        self.BLOCKNCX = BLOCKNCX
        self.STIMULUSNCX = STIMULUSNCX
        self.CAFEINA = CAFEINA
        self.BLOCKCICR = BLOCKCICR
        self.BLOCKIKs = BLOCKIKs
        self.BLOCKIKr = BLOCKIKr
        self.BLOCKItof = BLOCKItof
        self.BLOCKItos = BLOCKItos
        self.BLOCKINa = BLOCKINa
        self.BLOCKICaL = BLOCKICaL
        self.STIMULUSICaL = STIMULUSICaL

        self.t_ap = t_ap;
        
        self.L = L;

        self.L = self.L - self.t_ap;

        self.Ap = Ap;
        
        self.v_resting = v_resting;

        self.tap = tap;

#        if(w<10):
#            w=10;		
        self.w = w;

        self.f = f;
        
        self.Delay = 1000*1/self.f - self.w;

        self.A_inj = A_inj;
        
        tig = tig;

        self.D = self.L;                         # Duracao do primeiro trem de pulso
        
        self.w2 = 20;                         # Largura do pulso retangular no primeiro trem de pulso
        
        self.f2 = 1;                         # Frequencia em Hz do segundo trem de pulso
        
        self.Delay2 = 1000*1/self.f2 - self.w2;        # Delay em ms
        
        self.A_inj2 = 9.5;                   # Amplitude da corrente injetada para o segundo trem de pulso
        
        self.tsim = self.t_ap + self.L + self.tap          ##ok<NOPTS> # tempo total de simulacao em ms
        
        self.vetorPulso(self.w,self.Delay,L, self.t_ap)
    
#tig = 2000;                         # instante a partir do qual deseja visualizar os graficos

###############################################################################################################################################################
#Initial Conditions
        p=0;
        mo=1.405627e-3;
        ho= 9.867005e-1;
        jo=9.915620e-1; 
        do=7.175662e-6; 
        fo=1.000681; 
        fcaBjo=2.421991e-2;
        fcaBslo=1.452605e-2;
        xtoso=4.051574e-3;
        ytoso=9.945511e-1; 
        xtofo=4.051574e-3; 
        ytofo= 9.945511e-1; 
        xkro=8.641386e-3; 
        xkso= 5.412034e-3;
        RyRro=8.884332e-1;
        RyRoo=8.156628e-7; 
        RyRio=1.024274e-7; 
        NaBjo=3.539892;
        NaBslo=7.720854e-1; 
        TnCLo=8.773191e-3; 
        TnCHco=1.078283e-1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
        TnCHmo=1.524002e-2; 
        CaMo=2.911916e-4; 
        Myoco=1.298754e-3; 
        Myomo=1.381982e-1;
        SRBo=2.143165e-3; 
        SLLjo=9.566355e-3; 
        SLLslo=1.110363e-1; 
        SLHjo=7.347888e-3; 
        SLHslo=7.297378e-2; 
        Csqnbo= 1.242988;
        Ca_sro=5.545201e-1; 
        Najo= 8.80329; 
        Naslo=8.80733; 
        Naio=8.80853; 
        Kio=135; 
        Cajo=1.737475e-4; 
        Caslo= 1.031812e-4; 
        Caio=8.597401e-5; 
        Vmo=-8.556885e+1; 
        rtoso=0.9946;
        C2ca = 1; 
        C1ca = 0; 
        I1Ca=0; 
        I2Ca = 0; 
        I1Ba = 0; 
        I2Ba =0; 
        mol= 0; 
        hol= 0;
        CNa1 = 0;
        CNa2 = 0;
        CNa3 = 1;
        IFNa = 0;
        I1Na = 0;
        I2Na = 0;
        ICNa2 = 0;
        ICNa3 = 0;
        C0fo = 0;
        C1fo = 1;
        C2fo = 0;
        C3fo = 0;
        CI0fo = 0;
        CI1fo = 0;
        CI2fo = 0; 
        CI3fo = 0;
        Ifo = 0; 
        C0so = 0; 
        C1so = 1; 
        C2so = 0; 
        C3so = 0;
        CI0so = 0;
        CI1so = 0;
        CI2so = 0;
        CI3so = 0;
        Iso = 0;
        hw = 0.0001;
        hp = 0.006;
        TSCa = 2.6695e-05;
        TSCaw = 1.0427e-06; 
        TSCap = 4.5130e-07;
        TSw = 9.6370e-07;
        TSp =1.7690e-06;
        N0o = 0.998770;
        N1o = 0.367612e-4;
        P0o = 0.112735e-3;
        P1o = 0.148856e-3;
        P2o = 0.408484e-3;
        CK1r0 = 1;#9.967e-1;
        CK2r0 = 0; #4.341e-4;
        CK3r0 = 0; #7.634e-5;
        IKr0 = 0; #1.533e-6;
        C0Ks = 0;
        C1Ks = 1;
        O1Ks = 0;
        C1 = 0;
        C2 = 0;
        C3 = 0;
        C4 = 1; 
        C5 = 0;
        C6 = 0;
        C7 = 0;
        C8 = 0; 
        C9 = 0; 
        C10 = 0;
        C11 = 0;
        C12 = 0;
        C13 = 0; 
        C14 = 0; 
        C15 = 0;
        O1= 0;
        auro = 0.417069e-3;
        iuro = 0.998543;
        aksso = 0.417069e-3;
        iksso = 1;

## Matrix - Initial Conditions

# Gating variables      
        Y0  = [mo, ho, jo, do, fo, fcaBjo, fcaBslo, xtoso, ytoso, xtofo, ytofo, xkro, xkso,
               RyRro, RyRoo, RyRio, NaBjo, NaBslo, TnCLo, TnCHco, TnCHmo, CaMo, Myoco, Myomo,
               SRBo, SLLjo, SLLslo, SLHjo, SLHslo, Csqnbo,
               Ca_sro, Najo, Naslo, Naio, Kio, Cajo, Caslo, Caio, Vmo, rtoso,
               C2ca, C1ca, I1Ca, I2Ca, I1Ba, I2Ba, mol, hol, CNa1, CNa2, CNa3, IFNa, I1Na, I2Na, ICNa2, ICNa3,
               C0fo, C1fo, C2fo, C3fo, CI0fo, CI1fo, CI2fo, CI3fo, Ifo, C0so, C1so, C2so, C3so, CI0so, CI1so, CI2so, CI3so, Iso,
               hw, hp, TSCa, TSCaw, TSCap, TSw, TSp, N0o, N1o, P0o, P1o, P2o, CK1r0, CK2r0, CK3r0, IKr0, C0Ks, C1Ks, O1Ks,
               C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, O1, auro, iuro, aksso, iksso];  
       
   
#################################################################################################################################################
        t_start = 0.0;
        t_final = self.tsim;
        delta_t = 1e-5;
        ynit = Y0;
        iteracoes = t_final * 5;		
        t= np.zeros(iteracoes)        
        tempo= np.zeros(iteracoes)		
        tempo = np.linspace(t_start,t_final,iteracoes)
		
        if(self.protocol==3 and self.w<50):
            vetCritico = [];
            a = 0;		
            while a < self.tsim:
                if ((a > self.t_ap) and a < self.t_ap + self.D) and (np.mod((a + (1000 - self.t_ap)),(self.w + self.Delay)))<= self.w:
                    vetCritico.append(a);
                a = a +1;
            vetCritico = np.array(vetCritico);
            solver = inte.odeint (self.calculos, ynit, tempo , tcrit=vetCritico, full_output=True)
            z = np.array(solver[0]);
            t = np.array(tempo);
        else:
            solver = inte.odeint (self.calculos, ynit, tempo , full_output=True)                   
            z = np.array(solver[0]);
            t = np.array(tempo);
		
        Cac = z[:,35];
        Cass = z[:,36];
        Cacy = z[:,37];
        Nac = z[:,31];
        Nass = z[:,32];
        Nacy = z[:,33];
        Ki = z[:,34];
        tauuu = z[:,12];
        tar = z[:,11];
        y8 = z[:,7];
        y9 = z[:,8];
        y10 = z[:,9];
        y11 = z[:,10];
        y40 = z[:,39];
        y4 = z[:,3];
        y5 = z[:,4];
        y6 = z[:,5];
        y7 = z[:,6];
        y1 = z[:,0];
        y2 = z[:,1];
        y3 = z[:,2];
        y15 = z[:,14];
        CaSR = z[:,30];
        Po = 1 - z[:,40] - z[:,41] - z[:,42] - z[:,43] - z[:,44] - z[:,45];
        y47 = z[:,46];
        y48 = z[:,47];
        y39 = z[:,38];
        ONa = 1 - (z[:,48] + z[:,49] + z[:,50] + z[:,51] + z[:,52] + z[:,53] + z[:,54] + z[:,55]);
        C0f = z[:,56];
        C1f = z[:,57];
        C2f = z[:,58];
        C3f = z[:,59];
        CI0f = z[:,60];
        CI1f = z[:,61];
        CI2f = z[:,62];
        CI3f = z[:,63];
        If = z[:,64];
        C0s = z[:,65];
        C1s = z[:,66];
        C2s = z[:,67];
        C3s = z[:,68];
        CI0s = z[:,69];
        CI1s = z[:,70];
        CI2s = z[:,71];
        CI3s = z[:,72];
        Is = z[:,73];
        Of = 1 - (z[:,56]+ z[:,57] + z[:,58] + z[:,59] + z[:,60] + z[:,61]+ z[:,62]+ z[:,63]+ z[:,64]);
        Os = 1 - (z[:,65]+ z[:,66] + z[:,67] + z[:,68] + z[:,69] + z[:,70]+ z[:,71]+ z[:,72]+ z[:,73]);
              
        hw = z[:,74];
        hp = z[:,75];
        TSCa = z[:,76];   #complexo TS + 3Ca
        TSCaw = z[:,77];  #complexo TS + 3Ca + miosina no estado fraco
        TSCap = z[:,78];  #complexo TS + 3Ca + miosina no estado forte
        TSp = z[:,79];    #complexo TS + miosina no estado forte
        TSw = z[:,80];    #complexo TS + miosina no estado fraco
        #N0 = y(:,82);
        #N1 = y(:,83);
        #P0 = y(:,84);
        #P1 = y(:,85);
        #P2 = y(:,86);
        #P3 = 1 - (N0 + N1 + P0 + P1 + P2);
        C1Kr = z[:,86];
        C2Kr = z[:,87];
        C3Kr = z[:,88];
        IKr = z[:,89];
        OKr = 1 - (z[:,86] + z[:,87] + z[:,88] + z[:,89]);
        C0Ks = z[:,90];
        C1Ks = z[:,91];
        O1Ks = z[:,92];
        O2Ks = 1 - (C0Ks + C1Ks + O1Ks);
        
        C1 = z[:,93]; 
        C2 = z[:,94];
        C3 = z[:,95];
        C4 = z[:,96];
        C5 = z[:,97]; 
         
        C6 = z[:,98];
        C7 = z[:,99];
        C8 = z[:,100];
        C9 = z[:,101];
        C10 = z[:,102];
         
        C11 = z[:,103];
        C12 = z[:,104];
        C13 = z[:,105];
        C14 = z[:,106];
        C15 = z[:,107];
        O1 = z[:,108]; 
         
        O2 = 1 - (O1 + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10 + C11 + C12 + C13 + C14 + C15);
        
        aur = z[:,109]; 
        iur = z[:,110];
        akss = z[:,111];
        ikss = z[:,112];

        if self.Model_Force == 0:
            TS = self.TSt - TSCa - TSCaw - TSCap - TSw - TSp;

        N0 = z[:,81];
        N1 = z[:,82];
        P0 = z[:,83]; 
        P1 = z[:,84];
        P2 = z[:,85];
        P3 = 1 - (N0 + N1 + P0 + P1 + P2);
        
        
        I_app = np.zeros(t.size)
        v = np.zeros(t.size)
        vv= np.zeros(t.size)
        #t= np.zeros(t.size)
        V_clamp= np.zeros(t.size)
        Ka_junc= np.zeros(t.size)
        Ka_sl= np.zeros(t.size)
        s1_junc= np.zeros(t.size)
        s1_sl= np.zeros(t.size)
        s2_junc= np.zeros(t.size)
        s3_junc= np.zeros(t.size)
        s2_sl= np.zeros(t.size)
        s3_sl= np.zeros(t.size)
        I_ncx_junc= np.zeros(t.size)
        I_ncx_sl= np.zeros(t.size)
        I_ncx= np.zeros(t.size)
        ek= np.zeros(t.size)
        aki= np.zeros(t.size)
        bki= np.zeros(t.size)
        I_ki= np.zeros(t.size)
        kp_kp= np.zeros(t.size)
        I_kp_junc= np.zeros(t.size)
        I_kp_sl= np.zeros(t.size)
        I_kp= np.zeros(t.size)
                
        pcaks_junc= np.zeros(t.size)
        pcaks_sl= np.zeros(t.size)
        gks_junc= np.zeros(t.size)
        gks_sl= np.zeros(t.size)
        eks= np.zeros(t.size)
        I_ks_junc= np.zeros(t.size)
        I_ks_sl= np.zeros(t.size)
        eks= np.zeros(t.size)
        gks= np.zeros(t.size)
        I_ks_sl= np.zeros(t.size)
                
        I_ks= np.zeros(t.size)
        
        ek= np.zeros(t.size)
        rkr= np.zeros(t.size)
        I_kr= np.zeros(t.size)
            
        #IKur    
        IKur= np.zeros(t.size)
        
        IKss= np.zeros(t.size)
        
        fnak= np.zeros(t.size)
        I_nak_junc= np.zeros(t.size)
        I_nak_sl= np.zeros(t.size)
        I_nak= np.zeros(t.size)
        
        
        I_tof= np.zeros(t.size)
        I_tos= np.zeros(t.size)
                    
        fcaCaMSL= np.zeros(t.size)
        fcaCaj= np.zeros(t.size)
        ibarca_j= np.zeros(t.size)
        ibarca_sl= np.zeros(t.size)
        I_Ca_junc= np.zeros(t.size)
        I_Ca_sl= np.zeros(t.size)
                    
            
        I_Ca= np.zeros(t.size)
        
        
        ena_sl= np.zeros(t.size)
        ena_junc= np.zeros(t.size)
        I_NaL_sl= np.zeros(t.size)
        I_NaL_junc= np.zeros(t.size)
                    

        ENa= np.zeros(t.size)
        I_Na_sl= np.zeros(t.size)
        I_Na_junc= np.zeros(t.size)
            
        
               
        I_Na= np.zeros(t.size)
        
        DifNa_sl_cl= np.zeros(t.size)
        
        I_pca_junc= np.zeros(t.size)
        I_pca_sl= np.zeros(t.size)
        I_pca= np.zeros(t.size)
                
        eca_junc= np.zeros(t.size)
        eca_sl= np.zeros(t.size)
        I_cabk_junc= np.zeros(t.size)
        I_cabk_sl= np.zeros(t.size)
        I_cabk= np.zeros(t.size)
        endena_sl= np.zeros(t.size)
        J_SRCarel= np.zeros(t.size)
        J_serca= np.zeros(t.size)
        J_SRleak= np.zeros(t.size)
        
        Fb= np.zeros(t.size)
        Fp= np.zeros(t.size)
        FORCA= np.zeros(t.size)
        Ls= np.zeros(t.size)
        Lsim= np.zeros(t.size)
        CelL= np.zeros(t.size)
            

        Fcontrn= np.zeros(t.size)
        Fcontr= np.zeros(t.size)
        SL= np.zeros(t.size)
        CelL= np.zeros(t.size)
        self.kisss = np.zeros(t.size)
        ibarca_j= np.zeros(t.size)
        i=0;            
        while i < t.size:
            if (self.protocol==1):
        
                I_app[i] = 0;
            
            elif(self.protocol==2):
        
                if t[i] <= self.t_ap:                                            #[ms]
    
                    vv[i] = self.v_resting;                                          #[mV]
  
                elif (t[i] > self.t_ap) and (t[i] <= self.t_ap + self.L):    
                    A = -self.v_resting + self.Ap
                    sig = self.pulse_train(t[i])
                    vv[i] = A * sig + self.v_resting 
                else:
                    vv[i] = self.v_resting;    #[mV]
                    
                V_clamp[i] = vv[i];
                R_clamp = 0.02;
                I_app[i] = (V_clamp[i]- y39[i]/R_clamp);          
     
            elif(self.protocol==3):
                if ((t[i] > self.t_ap) and t[i] < self.t_ap + self.D) and (np.mod((t[i] + (1000 - self.t_ap)),(self.w + self.Delay)))<= self.w:                                				
                    I_app[i] = self.A_inj;    
                elif (t[i] > self.t_ap + self.D) and ((np.mod(t[i] + (1000 - self.t_ap)),(self.w2 + self.Delay2)))<= self.w2:        ##ok<AND2> #corrente injetada de 10 uA/uF com frequencia de 4 Hz durante 1s
                    I_app[i] = self.A_inj2;
                else:
                    I_app[i] = 0;
            i=i+1;
        
        i=0;              
        while i < t.size:

            
            if self.protocol==1:    
                A = -self.v_resting + self.Ap

                sig = self.pulse_train(t[i])
                v[i] = A * sig + (self.v_resting)                     

            elif(self.protocol==2):
                v[i] = z[i,38];    
            elif(self.protocol==3):
                v[i] = z[i,38];

            Ka_junc[i] = 1/(1+(self.Kdact/Cac[i])**3);
            Ka_sl[i] = 1/(1+(self.Kdact/Cass[i])**3);
            s1_junc[i] = exp(self.nu*v[i]*self.FoRT)*Nac[i]**3*self.Cao;
            s1_sl[i] = exp(self.nu*v[i]*self.FoRT)*Nass[i]**3*self.Cao;
            s2_junc[i] = exp((self.nu-1)*v[i]*self.FoRT)*self.Nao**3*Cac[i];
            s3_junc[i] = self.KmCai*self.Nao**3*(1+(Nac[i]/self.KmNai)**3) + self.KmNao**3*Cac[i]*(1+Cac[i]/self.KmCai)+self.KmCao*Nac[i]**3+Nac[i]**3*self.Cao+self.Nao**3*Cac[i];
            s2_sl[i] = exp((self.nu-1)*v[i]*self.FoRT)*self.Nao**3*Cass[i];
            s3_sl[i] = self.KmCai*self.Nao**3*(1+(Nass[i]/self.KmNai)**3) + self.KmNao**3*Cass[i]*(1+Cass[i]/self.KmCai)+self.KmCao*Nass[i]**3+Nass[i]**3*self.Cao+self.Nao**3*Cass[i];
            I_ncx_junc[i] = self.Fjunc*self.IbarNCX*self.Q10NCX**self.Qpow*Ka_junc[i]*(s1_junc[i]-s2_junc[i])/s3_junc[i]/(1+self.ksat*exp((self.nu-1)*v[i]*self.FoRT));
            I_ncx_sl[i] = self.Fsl*self.IbarNCX*self.Q10NCX**self.Qpow*Ka_sl[i]*(s1_sl[i]-s2_sl[i])/s3_sl[i]/(1+self.ksat*exp((self.nu-1)*v[i]*self.FoRT));
            I_ncx[i] = I_ncx_junc[i]+I_ncx_sl[i];
        
            #I_ki: Time-Independent K Current
            ek[i] = (1/self.FoRT)*np.log(self.Ko/Ki[i]);
            aki[i] = 1.02/(1+exp(0.2385*(v[i]-ek[i]-59.215)));
            bki[i] =(0.49124*exp(0.08032*(v[i]+5.476-ek[i])) + exp(0.06175*(v[i]-ek[i]-594.31))) /(1 + exp(-0.5143*(v[i]-ek[i]+4.753)));
            self.kisss[i] = aki[i]/(aki[i]+bki[i]);
            I_ki[i] = 0.9*np.sqrt(self.Ko/5.4)*self.kisss[i]*(v[i]-ek[i]);
        
            #I_kp: Plateau K current
            kp_kp[i] = 1/(1+exp(7.488-v[i]/5.98));
            I_kp_junc[i] = self.Fjunc*self.gkp*kp_kp[i]*(v[i]-ek[i]);
            I_kp_sl[i] = self.Fsl*self.gkp*kp_kp[i]*(v[i]-ek[i]);
            I_kp[i] = I_kp_junc[i]+I_kp_sl[i];
            
            # I_ks: Slowly Activating K Current
            pcaks_junc[i] = -np.log10(Cac[i])+3.0; 
            pcaks_sl[i] = -np.log10(Cass[i])+3.0;  
            gks_junc[i] = 0.07*(1 - self.BLOCKIKs)*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc[i])/0.6)));
            gks_sl[i] = 0.07*(1 - self.BLOCKIKs)*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl[i])/0.6))); 
            eks[i] = (1/self.FoRT)*np.log((self.Ko+self.pNaK*self.Nao)/(Ki[i]+self.pNaK*Nacy[i]));    
            I_ks_junc[i] = self.Fjunc*gks_junc[i]*tauuu[i]**2*(v[i]-eks[i]);
            if self.Model_IKs ==1:
                I_ks_sl[i] = self.Fsl*self.GKs*(O1Ks[i] + O2Ks[i])*(v[i]-eks[i]); #[uA/uF]    
            elif self.Model_IKs ==2:
                eks[i] = (1/self.FoRT)*np.log((self.Ko+self.pNaK*Nass[i])/(Ki[i]+self.pNaK*Nass[i]));    
                gks[i] = (1 - self.BLOCKIKs)*self.fatorgks*(1 + (0.6/(1 + (3.8e-5/Cass[i])**1.4)));    
                I_ks_sl[i] = self.Fsl*gks[i]*(O1[i] + O2[i])*(v[i]-eks[i]); #[uA/uF]     
            else:
                I_ks_sl[i] = self.Fsl*gks_sl[i]*tauuu[i]**2*(v[i]-eks[i]);
            
            I_ks[i] = I_ks_junc[i]+ I_ks_sl[i];
    
    # I_kr: Rapidly Activating K Current
            ek[i] = (1/self.FoRT)*np.log(self.Ko/Ki[i]);
            if self.Model_IKr==1:
                I_kr[i] = self.gkr*OKr[i]*(v[i]-ek[i]);
            else:         
                rkr[i] = 1/(1+exp((v[i]+33)/22.4));
                I_kr[i] = self.gkr*tar[i]*rkr[i]*(v[i]-ek[i]);
        
    #IKur    
            IKur[i] = self.GKur*aur[i]*iur[i]*(v[i] - ek[i]);   
    
    #IKss
            IKss[i] = self.GKss*akss[i]*ikss[i]*(v[i]-ek[i]);
    
    # I_nak: Na/K Pump Current
            sigma = (exp(self.Nao/67.3)-1)/7;
            fnak[i] = 1/(1+0.1245*exp(-0.1*v[i]*self.FoRT)+0.0365*sigma*exp(-v[i]*self.FoRT));
            I_nak_junc[i] = self.Fjunc*self.IbarNaK*fnak[i]*self.Ko/(1+(self.KmNaip/Nac[i])**4)/(self.Ko+self.KmKo);
            I_nak_sl[i] = self.Fsl*self.IbarNaK*fnak[i]*self.Ko /(1+(self.KmNaip/Nass[i])**4)/(self.Ko+self.KmKo);
            I_nak[i] = I_nak_junc[i] + I_nak_sl[i];
    
    # I_to: Transient Outward K Current (slow and fast components)
    
            if self.Model_Ito == 1:
                I_tof[i] = self.GItof*Of[i]*(v[i]-ek[i]); #[uA/uF]   
                I_tos[i] = self.GItos*Os[i]*(4*v[i]*self.Frdy*self.FoRT)*((Ki[i]*exp(v[i]*self.FoRT) - self.Ko) /(exp(v[i]*self.FoRT)-1) + 0.02*(Nass[i]*exp(v[i]*self.FoRT) - self.Nao) /(exp(v[i]*self.FoRT)-1)); 
    #Itos = ito*(10**7)*SAsl/Cmem;
            else: 
                I_tos[i] = self.GtoSlow*y8[i]*(y9[i]+0.5*y40[i])*(v[i]-ek[i]);    # [uA/uF]
                I_tof[i] = self.GtoFast*y10[i]*y11[i]*(v[i]-ek[i]);
                
            #ICal
            ibarca_j[i] = self.pCa*4*(v[i]*self.Frdy*self.FoRT)*(self.gCai*Cac[i]*exp(2*v[i]*self.FoRT)-self.gCao*self.Cao) /(exp(2*v[i]*self.FoRT)-1);
            ibarca_sl[i] = self.pCa*4*(v[i]*self.Frdy*self.FoRT)*(self.gCai*Cass[i]*exp(2*v[i]*self.FoRT)-self.gCao*self.Cao) /(exp(2*v[i]*self.FoRT)-1);   
    
            if self.Model_ICaL== 1:
                ibarca_j[i] = self.PCam*4*(v[i]*self.Frdy*self.FoRT)*(self.gCai*Cac[i]*exp(2*v[i]*self.FoRT)-self.gCao*self.Cao) /(exp(2*v[i]*self.FoRT)-1);
                ibarca_sl[i] = self.PCam*4*(v[i]*self.Frdy*self.FoRT)*(self.gCai*Cass[i]*exp(2*v[i]*self.FoRT)-self.gCao*self.Cao) /(exp(2*v[i]*self.FoRT)-1);
                I_Ca_junc[i] = self.Fjunc_CaL*ibarca_j[i]*Po[i]*self.ICa_scale;
                I_Ca_sl[i] = self.Fsl_CaL*ibarca_sl[i]*Po[i]*self.ICa_scale;
                
            else:
                fcaCaMSL[i]= 0.1/(1+(0.01/Cass[i])); 
                fcaCaj[i]= 0.1/(1+(0.01/Cac[i])); 
                ibarca_j[i] = self.pCa*4*(v[i]*self.Frdy*self.FoRT)*(self.gCai*Cac[i]*exp(2*v[i]*self.FoRT)-self.gCao*self.Cao) /(exp(2*v[i]*self.FoRT)-1);
                ibarca_sl[i] = self.pCa*4*(v[i]*self.Frdy*self.FoRT)*(self.gCai*Cass[i]*exp(2*v[i]*self.FoRT)-self.gCao*self.Cao) /(exp(2*v[i]*self.FoRT)-1);
                I_Ca_junc[i] = (self.Fjunc_CaL*ibarca_j[i]*y4[i]*y5[i]*((1-y6[i])+fcaCaj[i])*self.Q10CaL**self.Qpow)*0.45*1;
                I_Ca_sl[i] = (self.Fsl_CaL*ibarca_sl[i]*y4[i]*y5[i]*((1-y7[i])+fcaCaMSL[i])*self.Q10CaL**self.Qpow)*0.45*1;
                
        
            I_Ca[i] = I_Ca_junc[i] + I_Ca_sl[i];
    
    # I_Na: Fast Na Current
            I_NaL_sl[i] = self.Fsl*self.GNaL*y47[i]**3*y48[i]*(v[i]-ena_sl[i]);
            I_NaL_junc[i] = self.Fjunc*self.GNaL*y47[i]**3*y48[i]*(v[i]-ena_junc[i]);

            if self.Model_Na == 1:
                ENa[i] = (1/self.FoRT)*np.log((0.9*self.Nao + 0.1*self.Ko)/(0.9*Nass[i] + 0.1*Ki[i])); #Ki[i]
                I_Na_sl[i] = self.Fsl*self.GNam*ONa[i]*(v[i] - ENa[i])*self.I_scale_Na + I_NaL_sl[i];
                I_Na_junc[i] =0;
        
            else:
                ena_sl[i] = (1/self.FoRT)*np.log(self.Nao/Nass[i]);       # [mV]
                ena_junc[i] = (1/self.FoRT)*np.log(self.Nao/Nac[i]);     # [mV]

                I_Na_junc[i] = self.Fjunc*self.GNa*y1[i]**3*y2[i]*y3[i]*(v[i]-ena_junc[i])+ I_NaL_junc[i];
                I_Na_sl[i] = self.Fsl*self.GNa*y1[i]**3*y2[i]*y3[i]*(v[i]-ena_sl[i]) + I_NaL_sl[i];

                endena_sl[i] = (1/self.FoRT)*np.log(self.Nao/Nass[i]);       # [mV]
                ena_junc[i] = (1/self.FoRT)*np.log(self.Nao/Nac[i]);     # [mV]
                I_Na_junc[i] = self.Fjunc*self.GNa*y1[i]**3*y2[i]*y3[i]*(v[i]-ena_junc[i])+ I_NaL_junc[i];                    
                I_Na_sl[i] = self.Fsl*self.GNa*y1[i]**3*y2[i]*y3[i]*(v[i]-ena_sl[i]) + I_NaL_sl[i];
                

           
            I_Na[i] = I_Na_junc[i] + I_Na_sl[i];
    
    # Difusao
            DifNa_sl_cl[i] = self.J_na_juncsl/self.Vjunc*(Nass[i]-Nac[i]);
    
            # I_pca: Sarcolemmal Ca Pump Current
            I_pca_junc[i] = self.Fjunc*self.Q10SLCaP**self.Qpow*self.IbarSLCaP*Cac[i]**1.6/(self.KmPCa**1.6+Cac[i]**1.6);
            I_pca_sl[i] = self.Fsl*self.Q10SLCaP**self.Qpow*self.IbarSLCaP*Cass[i]**1.6/(self.KmPCa**1.6+Cass[i]**1.6);
            I_pca[i] = I_pca_junc[i] + I_pca_sl[i];
            
            # I_cabk: Ca Background Current
            eca_junc[i] = (1/self.FoRT/2)*np.log(self.Cao/Cac[i]);   # [mV]
            eca_sl[i] = (1/self.FoRT/2)*np.log(self.Cao/Cass[i]);     # [mV]
            I_cabk_junc[i] = self.Fjunc*self.GCaB*(v[i]- eca_junc[i]);
            I_cabk_sl[i] = self.Fsl*self.GCaB*(v[i]- eca_sl[i]);
            I_cabk[i] = I_cabk_junc[i] + I_cabk_sl[i];
            
            # SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
            J_SRCarel[i] = self.ks*(1 -self.BLOCKCICR)*y15[i]*(CaSR[i]-Cac[i]);          # [mM/ms]
            J_serca[i] = self.Q10SRCaP**self.Qpow*self.Vmax_SRCaP*((Cacy[i]/self.Kmf)**self.hillSRCaP-(CaSR[i]/self.Kmr)**self.hillSRCaP)/(1+(Cacy[i]/self.Kmf)**self.hillSRCaP+(CaSR[i]/self.Kmr)**self.hillSRCaP);
            J_SRleak[i] = self.GSrLeak*(CaSR[i]- Cacy[i]); 
    
            #Forca
            if self.Model_Force == 0:
                Fb[i] = self.Kw*(TSCaw[i] + TSw[i])*hw[i] + self.Kp*(TSCap[i] + TSp[i])*hp[i];
                Fp[i] = self.Ke*(self.Lsarc - self.Lo)**5 + self.Le*(self.Lsarc - self.Lo);
                FORCA[i] = Fp[i] + Fb[i];
                Ls[i] = (1/self.betaneg)*np.log((FORCA[i]/self.alfaneg)+1);
                Lsim[i] = self.Lsarc - (Ls[i] - self.fit);
                CelL[i] = self.cellLength *(Lsim[i]/self.Lsarc);
        
            else:
                Fcontrn[i] = -(P1[i] + N1[i] + 2*P2[i] + 3*P3[i])/self.Fmax;
                Fcontr[i] = self.alfaaj*Fcontrn[i];
                SL[i] = self.betaaj*Fcontrn[i] + self.Lsarc;
                CelL[i] = self.cellLength *(SL[i]/self.Lsarc);
    
            i=i+1;
## Difusao

        DifNa_sl_cl = self.J_na_juncsl/self.Vjunc*(z[:,32]-z[:,31]);
        C0KsC1Ks = C0Ks + C1Ks
        O1KsO2Ks = O1Ks + O2Ks
        C1KrC2KrC3Kr=C1Kr + C2Kr + C3Kr
        canalSerca = 1 - (z[:,13]+ z[:,14] + z[:,15])
        z4041=z[:,40] + z[:,41]
        z4243=z[:,42] + z[:,43]
        z4445=z[:,44] + z[:,45]
        z484950=z[:,48]+ z[:,49]+ z[:,50]
        z515253=z[:,51]+ z[:,52]+ z[:,53]
        z5455=z[:,54]+ z[:,55]
        C0fC1fC2fC3f = C0f + C1f + C2f + C3f
        CI0fCI1fCI2fCI3f = CI0f + CI1f + CI2f + CI3f	
        C0sC1sC2sC3s = C0s + C1s + C2s + C3s		
        CI0sCI1sCI2sCI3s = CI0s + CI1s + CI2s + CI3s
		
     
		
        #py.sign_in('cardiolab', 'K3sUD1esgJduklRiy0Lp') # Replace the username, and API key with your credentials.
        return [t,I_app,v,I_Ca,I_Na,z[:,37],z[:,35],z[:,36],z[:,30],I_kr,I_ks,C0KsC1Ks,O1KsO2Ks,OKr,C1KrC2KrC3Kr,I_kp,I_ki,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,O1,O2,J_SRCarel,J_serca,J_SRleak,z[:,13],z[:,15],z[:,14],canalSerca,I_ncx,I_nak,z[:,34],z[:,32],z[:,31],z[:,33],DifNa_sl_cl,I_cabk_junc,I_cabk_sl,I_cabk,I_pca_junc,I_pca_sl,I_pca,Po,z[:,40],z[:,41],z[:,42],z[:,43],z[:,44],z[:,45],z4041,z4243,z4445,z[:,3],z[:,4],ONa,z484950,z515253,z5455,        z[:,0],z[:,1],z[:,2],I_tos,I_tof,Of,If,Os,Is,C0fC1fC2fC3f,CI0fCI1fCI2fCI3f,C0sC1sC2sC3s,CI0sCI1sCI2sCI3s,z[:,7],        z[:,8],z[:,9],z[:,10],FORCA,Lsim,Cacy,CelL,Fcontr,N0,N1,P0,P1,P2,P3,SL,IKur,IKss];
        
#################################################################################################################################################
