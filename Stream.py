# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 15:24:44 2020

@author: paulm
"""
from collections import OrderedDict 
from chempy.util.arithmeticdict import ArithmeticDict
from chempy.chemistry import Species
import periodictable as pt 
#from Kiln_Module import S, G

class Stream(object):
    
    def __init__(self, name, frate, temp, pres, comp_d):
        self.name = name
        self.frate = frate
        self.temp = temp
        self.pres = pres
        self.comp = comp_d
        
    def Comp_wf(self):
        c = ArithmeticDict(float, self.comp)
        sum_comp = sum(c.values())
        first_key = list(c.keys())[0]
        if sum_comp < 1:
            c[first_key] += (1 - sum_comp)            
        if sum_comp > 1:
            c[first_key] -= (1 - sum_comp)              
        return c
    
    def Ele_wf(self):
        ele_list = []
        ele_wtfr = []
        for k,v in self.Comp_wf().items():
            compound = Species.from_formula(k).composition   
            for i,j in compound.items():
                m = "{}".format(pt.elements[i])
                h = pt.elements[i].mass*j/Species.from_formula(k).mass*v 
                #print("m = {}".format(m))
                #print("pt.elements[i].mass ={}".format(pt.elements[i].mass))
                #print("compound = {}".format(Species.from_formula(k).mass))
                ele_list.append(m)                    
                ele_wtfr.append(h)
                #print("ele_list = {}".format(ele_list))
                #print("ele_wtfr = {}".format(ele_wtfr))
        #print("length ele_List {0}, ele_wtfr {1}".format(len(ele_list), len(ele_wtfr)))
        #duplicates = len(ele_list) > len(set(ele_list))
        #if duplicates:
        ele = ArithmeticDict(float)
        for i,j in zip(ele_list, ele_wtfr):
            ele[i] += j
        #print("Ele = {}".format(ele))
        #print("Comp = {}".format(self.Comp_wf()))
        return ele  
     
    def Comp_mflow(self):
        comp_mf_inter = {k: v*self.frate/Species.from_formula(k).mass for (k, v) in self.Comp_wf().items()}
        comp_mf = ArithmeticDict(float, comp_mf_inter)
        return comp_mf
    
    def Ele_mw(self):
        ele_mf_inter = {k: v/pt.elements(k).mass for (k, v) in self.Ele_wf().items()}
        ele_mf = ArithmeticDict(float, ele_mf_inter)
        return ele_mf 
    
    def __add__(self,other):        
        return self.Comp_mflow() + other.Comp_mflow() 
    
    def __sub__(self,other):        
        return self.Comp_mflow() - other.Comp_mflow()  
    
    def __mul__(self,other):
        return self.Comp_mflow()*other
    
