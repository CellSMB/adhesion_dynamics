# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 09:44:21 2019

@author: qiliny
"""

"""
Given a force 'F', the program will output one kinetic constant for each of the
following reactions.
alpha-catenin
formin
myosin
"""
import numpy as np

#F_alpha = 7.9e-12
t = 0.1
k10z,k20z,k12z,k21z = 11, 0.14, 3, 20   # 1/s
x10,x20,x12,x21 = 0.0, 0.4e-9, 0.2e-9, 4.0e-9    # nm
kB = 1.38064852e-23                  # Boltzmann constant
T = 298.0                                 # absolute temprature 
kBT = kB*T
# alpha-catenin
# parameters setting

def alpha_dis(F_alpha):
    
    k10 = k10z*np.exp(F_alpha*x10/(kBT))
    k20 = k20z*np.exp(F_alpha*x20/(kBT))
    k12 = k12z*np.exp(F_alpha*x12/(kBT))
    k21 = k21z*np.exp(-F_alpha*x21/(kBT))
    
    b = k21+k20+k12+k10
    c = k21*k10+k10*k20+k12*k20
    lambda1 = (b+np.sqrt(b**2-4*c))/2.0
    lambda2 = (b-np.sqrt(b**2-4*c))/2.0
    B01 = (k21z*k10z)/(k21z*k10z+k12z*k20z)
    B02 = 1-B01
    
    C1 = (k21+k12+B01*k20+B02*k10-lambda1)/(lambda2-lambda1)
    C2 = 1-C1
    
#    f = C1*lambda1*np.exp(-lambda1*t)+C2*lambda2*np.exp(-lambda2*t)
    
    tao = C1/lambda1+C2/lambda2
    return tao

# Formin
# Parameters setting
#F_formin = 3.0e-12
p0 = 0.56
kon = 84.8
#C = 3.0
Cc = 0.07
sigma = 5.4e-9

def formin_elong(C,F_formin):
    
    p = 1.0/(1+((1-p0)/p0)*np.exp(-F_formin*sigma/(kB*T)))
    v_elong = kon*(C-Cc)*p
    return v_elong

# Myosin
#F_myosin_N = 6e-12  # force per myosin head

kc, ks = 176, 15
xc, xs = -2.5e-9, 0.40e-9
def myosin_dis(F_myosin_N):
    F_myosin_N = F_myosin_N*10**-12
    P = kc*np.exp(F_myosin_N*xc/(kB*T))+ks*np.exp(F_myosin_N*xs/(kB*T))
    return P

def k5r_dis_cont(con_myo_actin,init_myo_actin):
    k5r = (1.0*(1/(1+np.exp((con_myo_actin-init_myo_actin))))+0.5)*0.03
    return k5r
    
def ForceonEcad(F_cort,t):
    C = 965.9
    D = 4.105
    F_Ecad = F_cort*(D*(t/(t+C))+0.1)
    return F_Ecad



















