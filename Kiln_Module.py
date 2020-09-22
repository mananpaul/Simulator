# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 15:57:35 2020

@author: paulm
"""

from chempy.chemistry import balance_stoichiometry
from chempy.chemistry import Species
from chempy.util.arithmeticdict import ArithmeticDict
from Stream import Stream
import pandas as pd

class Kiln_Module(object):

    def __init__(self,name, S_Str1 = None, S_Str2 = None, G_Str1 = None, G_Str2 = None, Rxn_S = None,Rxn_L = None,Rxn_G = None,Ext_S = None,Ext_G = None, DR = None):
        self.name = name   
        self.S_Str1 = S_Str1
        self.S_Str2 = S_Str2
        #self.L_Str1 = L_Str1
        #self.L_Str2 = L_Str2
        self.G_Str1 = G_Str1
        self.G_Str2 = G_Str2
        self.Rxn_S = Rxn_S
        #self.Rxn_L = Rxn_L
        self.Rxn_G = Rxn_G
        self.Ext_S = Ext_S
        #self.Ext_L = Ext_L
        self.Ext_G = Ext_G
        self.DR = DR
        
    def Mass_Balance(self):  
        
        Str_S_wf = ArithmeticDict(float)
        Str_G_wf = ArithmeticDict(float)
        #Combined Solid Stream
        Solid_In = self.S_Str1 + self.S_Str2
        
        #print("Solid In = {}".format(Solid_In))
        #Combined Liquid Stream        
        #Liquid_In = self.L_Str1 + self.L_Str2
        
        #Combined Gas Stream
        Gas_In = self.G_Str1 + self.G_Str2      
                     
        Solid_Out = self.Rxn_Mol_Delta(self.Rxn_S, self.Ext_S, Solid_In)  
        
        #print("Solid Out = {}".format(Solid_Out))   
        
        for k,v in Solid_Out.items():
            ph_S = Species.from_formula(k).phase_idx        
            if ph_S == 3:
                Gas_In[k] = Gas_In[k]  + v
                Solid_Out[k] = Solid_Out[k] - v
                
        Gas_Out = self.Rxn_Mol_Delta(self.Rxn_G, self.Ext_G, Gas_In) 
        
        Str_S_Wflow = self.Comp_wflow(Solid_Out)
        Str_G_Wflow = self.Comp_wflow(Gas_Out)

        print("S_Wflow = {}".format(Str_S_Wflow))  
        print("G_Wflow = {}".format(Str_G_Wflow))
        
        """        
        for k,v in Str_S_Wflow.items():
            ph_S = Species.from_formula(k).phase_idx 
            if ph_S == 1:
                Str_G_Wflow[k] = Str_G_Wflow[k] + v*(self.DR) 
                Str_S_Wflow[k] = Str_S_Wflow[k] - v*(1-self.DR)
        """
        
        Str_S_wf = self.Comp_wf(Str_S_Wflow)
        Str_G_wf = self.Comp_wf(Str_G_Wflow)        
        
        return Str_S_wf, Str_G_wf, sum(Str_S_Wflow.values()), sum(Str_G_Wflow.values())

    def Rxn_Mol_Delta(self, Rxn,Ext,Str_In):     
        Str_Out = Str_In
        global reac, prod
        delta_reac = ArithmeticDict(float)
        delta_prod = ArithmeticDict(float)
        #print("Before Rxn Str = {}".format(Str_Out))
        for i,j in zip(Rxn,Ext):          
            rxn_split = i.split("=")
            reac_l = rxn_split[0].split("+")
            prod_l = rxn_split[1].split("+")  
            param = float(j)
            reac, prod = balance_stoichiometry(reac_l, prod_l)  
            in_moles = {k: v for k,v in Str_Out.items() if k in reac}
            lim_reac = min(in_moles.values())
            lim_reac_name = list({k: v for k,v in in_moles.items() if v == lim_reac}.keys())[0]  
            lim_reac_coeff = list({k: v for k,v in reac.items() if k == lim_reac_name}.values())[0]            
            for k,v in reac.items():
                delta_reac[k] = lim_reac*param*v/lim_reac_coeff                   
            for k,v in prod.items():
                delta_prod[k] = lim_reac*param*v/lim_reac_coeff       
            for k,v in Str_Out.items():
              if k in reac:
                  Str_Out[k] = v - delta_reac[k]
              if k in prod:
                  Str_Out[k] = v + delta_prod[k]
        #print("After Rxn Str = {}".format(Str_Out))
        return Str_Out

    def Bal_Rxn(self, Rxn):                
        reac_keys = list(reac.keys())
        reac_values = list(reac.values())
        prod_keys = list(prod.keys())    
        prod_values = list(prod.values())        
        total_reac = []
        total_prod = []        
        for i,j in zip(reac_keys,reac_values):
            total_reac.append(str(j)+str(i))             
        t_reac = '+'.join(map(str, total_reac))            
        for i,j in zip(prod_keys,prod_values):
            total_prod.append(str(j)+str(i))               
        t_prod = '+'.join(map(str, total_prod))        
        t_rxn = str(t_reac) + "=" + str(t_prod)     
        return t_rxn

    # Splitting Streams into Phases 
    
    def Dusting(self, Str_S, Str_G, DR):   
        for k,v in Str_S.items():            
            phase = Species.from_formula(k).phase_idx        
            if phase == 1:
                Str_G[k] = Str_G[k] + v*DR
                Str_S[k] = Str_S[k] - v*(1-DR)          
            #if phase == 2:
            #    Str_L[k] = v
            if phase == 3:
                Str_G[k] = Str_G[k] + v 
                Str_S[k] = Str_S[k] - v 
        return Str_S, Str_G

    def Comp_wflow(self,Str_In):
         comp_wf_inter = {k: v*Species.from_formula(k).mass for (k, v) in Str_In.items()}
         comp_wfl = ArithmeticDict(float, comp_wf_inter)
         return comp_wfl
    
    def Comp_wf(self,Str_In):
         sumcomp = sum(Str_In.values())
         comp_wf_inter = {k: v/sumcomp for (k, v) in Str_In.items()}
         comp_wf = ArithmeticDict(float, comp_wf_inter)
         return comp_wf

Rxn_S_1 = "NiO(s)+C(s)=Ni(s)+CO(g)"
Rxn_S_2 = "Fe2O3(s)+C(s)=Fe3O4(s)+CO(g)"
Rxn_S_3 = "Fe3O4(s)+C(s)=FeO(s)+CO(g)"
Rxn_S_4 = "FeO(s)+C(s)=Fe(s)+CO(g)"
Rxn_S_5 = "CoO(s)+C(s)=Co(s)+CO(g)"
Rxn_S_6 = "H2O(s)=H2O(l)"
Rxn_S_7 = "C2H4(s)=C2H4(g)"
Rxn_S_8 = "H2O(l)=H2O(g)"

Rxn_G_1 = "C(s)+O2(g)=CO(g)"
Rxn_G_2 = "C2H4(s)=C2H4(g)"
Rxn_G_3 = "CO(g)+O2(g)=CO2(g)"
Rxn_G_4 = "C2H4(g)+O2(g)=CO2(g)+H2O(g)"

RxnExt_S = [0.9,1,1,0.1,0.9,0.9,1,1]  
#RxnExt_L = []
RxnExt_G = [1,1,0.99,1] 

RKF_M = {"Al2O3(s)": 0.03, "Fe(s)": 0.0, "FeO(s)": 0.0, "Fe3O4(s)": 0.0, "Fe2O3(s)": 0.24, "Ni(s)": 0.0, "NiO(s)": 0.03, "Co(s)": 0.0, "CoO(s)": 0.01, "H2O(l)":0.0, "SiO2(s)": 0.39, "MgO(s)": 0.2, "C(s)":0.0,"S(s)":0.0, "CO(g)":0.0, "CO2(g)": 0.0, "H2O(g)":0.0, "N2(g)":0.0, "O2(g)":0.0, "C2H4(s)": 0.0, "H2O(s)":0.1,"C2H4(g)":0.0}
Coal_R_M = {"Al2O3(s)": 0.1, "Fe(s)": 0.0, "FeO(s)": 0.0, "Fe3O4(s)": 0.0, "Fe2O3(s)": 0.0, "Ni(s)": 0.0, "NiO(s)": 0.0, "Co(s)": 0.0,"CoO(s)": 0.0, "H2O(l)":0.1, "SiO2(s)": 0.0, "MgO(s)": 0.0, "C(s)":0.45,"S(s)":0.0, "CO(g)":0.0, "CO2(g)": 0.0, "H2O(g)":0.0, "N2(g)":0.0, "O2(g)":0.0, "C2H4(s)": 0.35, "H2O(s)":0.0,"C2H4(g)":0.0}

Coal_B_M = {"Al2O3(s)": 0.1, "Fe(s)": 0.0, "FeO(s)": 0.0, "Fe3O4(s)": 0.0, "Fe2O3(s)": 0.0, "Ni(s)": 0.0, "NiO(s)": 0.0, "Co(s)": 0.0,"CoO(s)": 0.0, "H2O(l)":0.0, "SiO2(s)": 0.0, "MgO(s)": 0.0, "C(s)":0.45,"S(s)":0.0, "CO(g)":0.0, "CO2(g)": 0.0, "H2O(g)":0.0, "N2(g)":0.0, "O2(g)":0.0, "C2H4(s)": 0.45, "H2O(s)":0.0,"C2H4(g)":0.0}
Air_B_M = {"Al2O3(s)": 0.0, "Fe(s)": 0.0, "FeO(s)": 0.0, "Fe3O4(s)": 0.0, "Fe2O3(s)": 0.0, "Ni(s)": 0.0, "NiO(s)": 0.0, "Co(s)": 0.0,"CoO(s)": 0.0, "H2O(l)":0.0, "SiO2(s)": 0.0, "MgO(s)": 0.0, "C(s)":0.0,"S(s)":0.0, "CO(g)":0.1, "CO2(g)": 0.0, "H2O(g)":0.2, "N2(g)":0.8, "O2(g)":0.1, "C2H4(s)": 0.0, "H2O(s)":0.0,"C2H4(g)":0.0}

RKF = Stream("RKF",1,30,101.325,RKF_M)
Coal_R = Stream("Reductant Coal",0.06,30,101.325,Coal_R_M)
Coal_B = Stream("Burner Coal",0.07,30,101.325,Coal_B_M)
Air_B = Stream("Burner Air",1,30,101.325,Air_B_M)

K_S_Rxn = [Rxn_S_1,Rxn_S_2,Rxn_S_3,Rxn_S_4,Rxn_S_5,Rxn_S_6,Rxn_S_7,Rxn_S_8]
K_G_Rxn = [Rxn_G_1,Rxn_G_2,Rxn_G_3,Rxn_G_4]

K = Kiln_Module("Kiln", S_Str1 = RKF, S_Str2 = Coal_R, G_Str1 = Coal_B, G_Str2 = Air_B, Rxn_S = K_S_Rxn, Rxn_G = K_G_Rxn, Ext_S = RxnExt_S, Ext_G = RxnExt_G, DR = 0.1)

S,G, Ssum, Gsum = K.Mass_Balance()

RKP = Stream("RKP",Ssum,30,101.325, S)
RKG = Stream("RKG", Gsum,30, 101.325, G)

print("Solid In = {0}, comp = {1}, ele = {2}".format(RKF.frate,RKF.Comp_wf(),RKF.Ele_wf()))
print("Solid Out = {0}, comp = {1}, ele = {2}".format(RKP.frate,S,RKP.Ele_wf()))
print("Gas Out= {0}, comp = {1}".format(RKG.frate,RKG.Comp_wf()))

df = pd.DataFrame.from_dict(RKF.Comp_wf(), orient = "index")
df1 = pd.DataFrame.from_dict(RKF.Ele_wf(), orient = "index")
df2 = pd.DataFrame.from_dict(RKP.Comp_wf(), orient = "index")
df3 = pd.DataFrame.from_dict(RKP.Ele_wf(), orient = "index")
df4 = pd.DataFrame.from_dict(G, orient = "index")
df5 = pd.DataFrame.from_dict(RKG.Ele_wf(), orient = "index")

df.to_csv("RKF_comp.csv")
df1.to_csv("RKF_ele.csv")
df2.to_csv("RKP_comp.csv")
df3.to_csv("RKP_ele.csv")
df4.to_csv("RKG_comp.csv")
df5.to_csv("RKG_ele.csv")
