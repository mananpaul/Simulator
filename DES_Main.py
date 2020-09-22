# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:23:20 2020

@author: paulm
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 20:36:03 2019

@author: Manan
"""

import simpy as sim
import numpy as np # a maths and plotting module
import pandas as pd # more data analysis
import matplotlib.pyplot as plt # 
import seaborn as sns
import math
import time
from tkinter import *
from tkinter.ttk import *
import random

import itertools
from Kiln_Module import Kiln_Module, RKF
from Container import Container

Random_Seed = 100
Surge_Bin_Size = 40
Furnace_Bin_Size = 18 # 16+2
Container_Size = 18 # 16+2
Container_Factor = 16
T_INTER = [10,50]
RKP_RANGE = [1,2]
power = 82
GHE = 0.465
Container_Heel = 1

#os.system("Reaction_Module.py &")

output_dict = {'Event Type':[], 'Time':[], 'Surge Bin Tonnes': [], 'Container Tonnes': [], 'FBin1':[],'FBin2':[],'FBin3':[],'FBin4':[],'FBin5':[],'FBin6':[],'FBin7':[],'FBin8':[],'FBin9':[]}
output_dict1 = {'Time':[], 'RKP Feed':[], 'Surge Bin Tonnes': [], 'Container Tonnes': [], 'FBin1':[],'FBin2':[],'FBin3':[],'FBin4':[],'FBin5':[],'FBin6':[],'FBin7':[],'FBin8':[],'FBin9':[]}

event = []
time1 = []
rkpt = []
sbt = []
ct = []
fb1t = []
fb2t = []
fb3t = []
fb4t = []
fb5t = []
fb6t = []
fb7t = []
fb8t = []
fb9t = []

#Rate_Feeding_Ring = (Power / GhE) * 2 / 5 * 1 / 60
#Rate_Feeding_Center = (Power / GhE) * 3 / 5 * 1 / 120

def write_data(time1, rkpt, sbt, ct, fb1t, fb2t, fb3t,fb4t,fb5t,fb6t,fb7t, fb8t, fb9t):
    output_dict1['Time'].append(time1)
    output_dict1['RKP Feed'].append(rkpt)   
    output_dict1['Surge Bin Tonnes'].append(sbt)
    output_dict1['Container Tonnes'].append(ct)
    output_dict1['FBin1'].append(fb1t)
    output_dict1['FBin2'].append(fb2t)
    output_dict1['FBin3'].append(fb3t)
    output_dict1['FBin4'].append(fb4t)
    output_dict1['FBin5'].append(fb5t)
    output_dict1['FBin6'].append(fb6t)
    output_dict1['FBin7'].append(fb7t)
    output_dict1['FBin8'].append(fb8t)
    output_dict1['FBin9'].append(fb9t)
    
def write_data1(event, time1, sbt, ct,fb1t, fb2t, fb3t,fb4t,fb5t,fb6t,fb7t, fb8t, fb9t):
    if event == 'SBL':
        output_dict['Event Type'].append('Surge Bin Load')        
    elif event == 'SBU':
        output_dict['Event Type'].append('Surge Bin Unload')
    elif event == 'FBL':
        output_dict['Event Type'].append('Furnace Bin Load')
    elif event == 'FBU':
        output_dict['Event Type'].append('Furnace Bin Unload')
    else:
        raise Exception('Event type code not properly defined in data write')
    output_dict['Time'].append(time1)
    output_dict['Surge Bin Tonnes'].append(sbt)
    output_dict['Container Tonnes'].append(ct)
    output_dict['FBin1'].append(fb1t)
    output_dict['FBin2'].append(fb2t)
    output_dict['FBin3'].append(fb3t)
    output_dict['FBin4'].append(fb4t)
    output_dict['FBin5'].append(fb5t)
    output_dict['FBin6'].append(fb6t)
    output_dict['FBin7'].append(fb7t)
    output_dict['FBin8'].append(fb8t)
    output_dict['FBin9'].append(fb9t)
 
def Surge_Bin_Load(env, SurgeBin, CalcineContainer, RKP_Feed): 
    if (SurgeBin.level <= Surge_Bin_Size): 
        event = 'SBL'
        yield SurgeBin.put(RKP_Feed)        
        print("Time = {}".format(env.now))
        print("SBin Level = {}".format(SurgeBin.level))
        write_data1(event, env.now, SurgeBin.level, CalcineContainer.level,FBin1.level, FBin2.level,FBin3.level,FBin4.level,FBin5.level,FBin6.level, FBin7.level, FBin8.level, FBin9.level)
        
def Surge_Bin_Unload(env,SurgeBin,CalcineContainer):     
    if SurgeBin.level >= Container_Factor and (CalcineContainer.level <= 1):                      
       event = 'SBU'
       yield CalcineContainer.put(Container_Factor)
       yield SurgeBin.get(Container_Factor)
       yield env.timeout(600)  
       #Container_Up()
       print("Container Level = {}".format(CalcineContainer.level))    
       write_data1(event, env.now, SurgeBin.level, CalcineContainer.level,FBin1.level, FBin2.level,FBin3.level,FBin4.level,FBin5.level,FBin6.level, FBin7.level, FBin8.level, FBin9.level)
   
def Furnace_Bin_Load(env,CalcineContainer,BinPriority):    
    Bin_List=[FBin1, FBin2,FBin3,FBin4,FBin5,FBin6,FBin7,FBin8,FBin9] 
    for i in Bin_List:
        if (i.level) <= 1 and (CalcineContainer.level >= Container_Factor):               
           event = 'FBL'
           print("Furnace Bin Level = {}".format(i.level))           
           yield i.put(Container_Factor)
           yield CalcineContainer.get(Container_Factor)
           yield env.timeout(600)  
           #Container_Down()
           print("Calcine Container Level After FBin Load = {}".format(CalcineContainer.level))
           write_data1(event, env.now, SurgeBin.level, CalcineContainer.level,FBin1.level, FBin2.level,FBin3.level,FBin4.level,FBin5.level,FBin6.level, FBin7.level, FBin8.level, FBin9.level)
       
def Furnace_Bin_Unload(env,BinPriority):      
    Bin_Center_Feed = [FBin2,FBin5,FBin8]
    Bin_Ring_Feed = [FBin1,FBin4,FBin7,FBin3,FBin6,FBin9]         
    for i in BinPriority:
        for j in Bin_Center_Feed:
            if i == j:
                if i.level >=1:
                    event = 'FBU'
                    yield i.get(0.08)
                    yield env.timeout(6)
                    print("Center FBin Bin level ===================================={}".format(i.level))
                    write_data1(event, env.now, SurgeBin.level, CalcineContainer.level,FBin1.level, FBin2.level,FBin3.level,FBin4.level,FBin5.level,FBin6.level, FBin7.level, FBin8.level, FBin9.level)                    
        for k in Bin_Ring_Feed:
            if i == k:
                if i.level >=1:
                    event = 'FBU'
                    yield i.get(0.06)
                    yield env.timeout(6)
                    print("Ring FBin Bin level ===================================={}".format(i.level))
                    #print("Calcine Container Level After Bin = {}".format(CalcineContainer.level))
                    write_data1(event, env.now, SurgeBin.level, CalcineContainer.level,FBin1.level, FBin2.level,FBin3.level,FBin4.level,FBin5.level,FBin6.level, FBin7.level, FBin8.level, FBin9.level)
              
def Control(env, SurgeBin, CalcineContainer, BinPriority):        
    for i in itertools.count():            
        yield env.timeout(1)
        #print("time step = {}".format(env.))
        RKF_F = np.random.normal(1,0.1)
        #RKF_F = 1
        RKF_Feed = RKF.frate
        print("RKP_Feed = {}".format(RKF_Feed))
        yield env.process(Surge_Bin_Load(env,SurgeBin,CalcineContainer, RKF_Feed))
        yield env.process(Surge_Bin_Unload(env,SurgeBin,CalcineContainer))        
        yield env.process(Furnace_Bin_Load(env,CalcineContainer,BinPriority))
        yield env.process(Furnace_Bin_Unload(env,BinPriority))
        write_data(i, RKF_Feed, SurgeBin.level, CalcineContainer.level, FBin1.level, FBin2.level,FBin3.level,FBin4.level,FBin5.level,FBin6.level, FBin7.level, FBin8.level, FBin9.level)
        
#random.seed(Random_Seed)

print('Container Simulation')

env = sim.Environment()

SurgeBin = Container(env, capacity= Surge_Bin_Size, init=10)
CalcineContainer = Container(env, capacity= Container_Size, init=0)
FBin1 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin2 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin3 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin4 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin5 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin6 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin7 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin8 = Container(env, capacity= Furnace_Bin_Size, init=10)
FBin9 = Container(env, capacity= Furnace_Bin_Size, init=10)

BinPriority = [FBin2,FBin5,FBin8,FBin2,FBin5,FBin8,FBin1,FBin4,FBin7,FBin2,FBin5,FBin8,FBin2,FBin5,FBin8,FBin3,FBin6,FBin9]

env.process(Control(env,SurgeBin, CalcineContainer, BinPriority))

env.run(until=5000)





