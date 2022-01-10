#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:56:04 2019

@author: qiliny
"""
import re
import numpy as np
import math
import random
import os
#import sys
#import time

#import matplotlib.pyplot as plt
#import matplotlib.animation as animation


from modules import *
from forcebased_k_new import *
#%%
#_________________________________________________________________________________________________
# INPUTS
#_________________________________________________________________________________________________

#-------------------------------------------------------------------------------------------------
# READING INPUT FILE
#-------------------------------------------------------------------------------------------------

#Open the input file
input_file = open("input_2.0_15.txt","r")

record_freq = 1000
#Create variable that extracts each line from text. Also shifts starting index from 0 to 1.
input_lines = input_file.readlines()
input_lines.insert(0,' ')
#Dimension of the model
number_D = 2

#Geometric Inputs
nodeinput, eleminput                                    = re.split(',|\n',input_lines[6])[0:2] #Link to nodes and elements created via triangle
nodeincenter, elemincenter			                    = re.split(',|\n',input_lines[7])[0:2] #Link to nodes and elements created via triangle
#Simulation runtime inputs
Tstep, Total_time, startT, freq_update, ODE_Timestep    = input_lines[13].split(',')           #Simulation Time Parameters
load_steps                                              = input_lines[16]                      #Output Frequence of Results
outputfreq                                              = input_lines[19]                      #Output Frequence of Results

#Reaction Diffusion Inputs
    #Actin Ecad Integrated model

Actin_Ecad                                              = input_lines[26].split('\n')[0]      #Link to Cellml model

init_E_trans, Dx_E_trans, Dy_E_trans                    = input_lines[29].split(',')          #E_trans inputs
init_actin, Dx_actin, Dy_actin                          = input_lines[32].split(',')          #actin inputs
init_actin_bran, Dx_actin_bran, Dy_actin_bran           = input_lines[35].split(',')          #actin_a inputs
init_actin_bund, Dx_actin_bund, Dy_actin_bund           = input_lines[38].split(',')          #actin_a inputs
init_Rac, Dx_Rac, Dy_Rac                                = input_lines[41].split(',')          #Rac inputs
init_Rac_a, Dx_Rac_a, Dy_Rac_a                          = input_lines[44].split(',')          #Rac_a inputs
init_Rho, Dx_Rho, Dy_Rho                                = input_lines[47].split(',')          #Rho inputs
init_Rho_a, Dx_Rho_a, Dy_Rho_a                          = input_lines[50].split(',')          #Rho_a inputs

init_myo, Dx_myo, Dy_myo                                = input_lines[53].split(',')          #myo inputs
init_myo_a, Dx_myo_a, Dy_myo_a                          = input_lines[56].split(',')          #myo_a inputs
init_E_trans_actin, Dx_E_trans_actin, Dy_E_trans_actin  = input_lines[59].split(',')          #E_trans_actin inputs
init_Avg_E_trans, Dx_Avg_E_trans, Dy_Avg_E_trans        = input_lines[62].split(',')          #E_trans inputs

init_myo_actin_bran, Dx_myo_actin_bran, Dy_myo_actin_bran   = input_lines[65].split(',')          #myo_actin inputs
init_myo_actin_bund, Dx_myo_actin_bund, Dy_myo_actin_bund   = input_lines[68].split(',')          #myo_actin inputs

init_actin_all, Dx_actin_all, Dy_actin_all              = input_lines[71].split(',')          #myo_actin inputs
init_myo_all, Dx_myo_all, Dy_myo_all                    = input_lines[74].split(',')          #myo_actin inputs
init_E_trans_all, Dx_E_trans_all, Dy_E_trans_all        = input_lines[77].split(',')          #myo_actin inputs

k1                      = input_lines[80].split('\t')[0]
k1p                     = input_lines[81].split('\t')[0]          
k2                      = input_lines[82].split('\t')[0]          
k3                      = input_lines[83].split('\t')[0]  
k3r                     = input_lines[84].split('\t')[0]     
k4                      = input_lines[85].split('\t')[0]          
k5r1                    = input_lines[86].split('\t')[0]          
k5r2                    = input_lines[87].split('\t')[0]          
k6                      = input_lines[88].split('\t')[0]
k6r                     = input_lines[89].split('\t')[0]          
k7                      = input_lines[90].split('\t')[0]           
k8                      = input_lines[91].split('\t')[0]           
k9                      = input_lines[92].split('\t')[0]          
k10                     = input_lines[93].split('\t')[0]          
k10r                    = input_lines[94].split('\t')[0]          
k11                     = input_lines[95].split('\t')[0]          
k12                     = input_lines[96].split('\t')[0]
                 
Ecad_conc               = input_lines[103]          # Concentration of the E-cad
side_len                = input_lines[105]          # Side length of a lattice(nm)
diffusivity             = input_lines[107]          # Diffusivity of the E - cad monomer and trans dimer nm ^ 2 / s
x_len                   = input_lines[109]          # The number of lattices alone the x axis
y_len                   = input_lines[111]          # The number of lattices alone the y axis
Kd_trans                = input_lines[113]         # Free energy change for trans dimer
Kd_cis                  = input_lines[115]         # Free energy change for cis

cortex_rad              = input_lines[128]

radius_um               = float(input_lines[135])
N0                      = float(input_lines[138])                   # Avogadro's number
C0                      = float(input_lines[141])
R                       = float(input_lines[144])                   # gas constant
k                       = float(input_lines[147])                   # Boltzmann constant
T                       = float(input_lines[150])                   # Temprature
h                       = float(input_lines[153])                   # height nm

F_cort                  = float(input_lines[156])                   # Cortical tension (N/m or kg/s^2) 900pN/um or 9e-4 N/m
cort_thick              = float(input_lines[159])                   # Thickness of the cortical actin (um)
con_cons                = float(input_lines[162])                   # Concentration constants
Fila_len                = float(input_lines[165])                   # Estimated Filament length 
real_radius             = float(input_lines[168])                   # radius of the real intercellular contact
F_myo                   = float(input_lines[171])                   # force generated by each myosin head

mode                    = int(input_lines[179])
fixed_angle             = float(input_lines[182])
init_angle              = np.deg2rad(float(input_lines[185]))       # initial theta calculated from Hertz equation.
T_fix_angle             = int(input_lines[188])                     # Time to fix the contact angle 

factor                  = int(input_lines[195])

#debug_mode
debug_mode              = int(input_lines[202])

input_file.close()

#-------------------------------------------------------------------------------------------------
# CONVERTING INPUTS FROM STRING TO APPROPRIATE DATA TYPES
#-------------------------------------------------------------------------------------------------

#Simulation Runtime Inputs
Tstep                       = float(Tstep)
Total_time                  = float(Total_time)           #Simulation Time Parameters
startT                      = float(startT)
freq_update                 = int(freq_update)
ODE_Timestep                = float(ODE_Timestep)

load_steps                  = int(load_steps)
outputfreq                  = int(outputfreq)
debug_mode                  = int(debug_mode)

#Reaction Diffusion Inputs
init_E_trans                = float(init_E_trans)        #E_trans inputs
Dx_E_trans                  = float(Dx_E_trans)
Dy_E_trans                  = float(Dy_E_trans)

init_actin                  = float(init_actin)        #actin inputs
Dx_actin                    = float(Dx_actin)
Dy_actin                    = float(Dy_actin)
          
init_actin_bran             = float(init_actin_bran)        #actin_a inputs
Dx_actin_bran               = float(Dx_actin_bran)
Dy_actin_bran               = float(Dy_actin_bran)    

init_actin_bund             = float(init_actin_bund)        #actin_a inputs
Dx_actin_bund               = float(Dx_actin_bund)
Dy_actin_bund               = float(Dy_actin_bund) 
    
init_Rac                    = float(init_Rac)        #Rac inputs
Dx_Rac                      = float(Dx_Rac)
Dy_Rac                      = float(Dy_Rac)    
    
init_Rac_a                  = float(init_Rac_a)        #Rac_a inputs
Dx_Rac_a                    = float(Dx_Rac_a)
Dy_Rac_a                    = float(Dy_Rac_a)    
    
init_Rho                    = float(init_Rho)        #Rho inputs
Dx_Rho                      = float(Dx_Rho)
Dy_Rho                      = float(Dy_Rho)    
    
init_Rho_a                  = float(init_Rho_a)        #Rho_a inputs
Dx_Rho_a                    = float(Dx_Rho_a)
Dy_Rho_a                    = float(Dy_Rho_a)    

init_myo                    = float(init_myo)        #myo inputs
Dx_myo                      = float(Dx_myo)
Dy_myo                      = float(Dy_myo)    
    
init_myo_a                  = float(init_myo_a)        #myo_a inputs
Dx_myo_a                    = float(Dx_myo_a)
Dy_myo_a                    = float(Dy_myo_a)    
    
init_E_trans_actin          = float(init_E_trans_actin)        #E_trans_actin inputs
Dx_E_trans_actin            = float(Dx_E_trans_actin)
Dy_E_trans_actin            = float(Dy_E_trans_actin)    

init_Avg_E_trans            = float(init_Avg_E_trans)        #Averger E_trans inputs
Dx_Avg_E_trans              = float(Dx_Avg_E_trans)
Dy_Avg_E_trans              = float(Dy_Avg_E_trans)

    
init_myo_actin_bran         = float(init_myo_actin_bran)        #myo_actin inputs
Dx_myo_actin_bran           = float(Dx_myo_actin_bran)
Dy_myo_actin_bran           = float(Dy_myo_actin_bran)

init_myo_actin_bund         = float(init_myo_actin_bund)        #myo_actin inputs
Dx_myo_actin_bund           = float(Dx_myo_actin_bund)
Dy_myo_actin_bund           = float(Dy_myo_actin_bund)

init_actin_all         = float(init_actin_all)        #myo_actin inputs
Dx_actin_all           = float(Dx_actin_all)
Dy_actin_all           = float(Dy_actin_all)

init_myo_all         = float(init_myo_all)        #myo_actin inputs
Dx_myo_all           = float(Dx_myo_all)
Dy_myo_all           = float(Dy_myo_all)

init_E_trans_all         = float(init_E_trans_all)        #myo_actin inputs
Dx_E_trans_all           = float(Dx_E_trans_all)
Dy_E_trans_all           = float(Dy_E_trans_all)

Ecad_conc                   = float(Ecad_conc)              # Concentration of the E-cad
side_len                    = float(side_len)               # Side length of a lattice(nm)
diffusivity                 = float(diffusivity)            # Diffusivity of the E - cad monomer and trans dimer nm ^ 2 / s
x_len                       = int(x_len)                  # The number of lattices alone the x axis
y_len                       = int(y_len)                  # The number of lattices alone the y axis
Kd_trans                    = float(Kd_trans)                # Free energy change for trans dimer
Kd_cis                      = float(Kd_cis)                  # Free energy change for cis

k1                      = float(k1)   
k1p                     = float(k1p)       
k2                      = float(k2)         
k3                      = float(k3)  
k3r                     = float(k3r)
k4                      = float(k4)
k5r1                    = float(k5r1)
k5r2                    = float(k5r2)
k6                      = float(k6)
k6r                     = float(k6r)
k7                      = float(k7)
k8                      = float(k8)
k9                      = float(k9)
k10                     = float(k10)
k10r                    = float(k10r)
k11                     = float(k11)
k12                     = float(k12)

cortex_rad              = float(cortex_rad)

# Time related parameters calculation
deltat = side_len ** 2 / (4 * diffusivity)          # lenght of each time step (s)
t_total_step = int(math.ceil(Total_time/deltat))    # total number of steps
Total_time_actual = t_total_step*deltat             # actual total time (s)

g_trans = (-k*N0*T*np.log(Kd_trans/1)/N0-k*T*np.log(side_len**2*h*10**-27*N0*C0*10**3))
g_cis = (-k*N0*T*np.log(Kd_cis/1)/N0-k*T*np.log(side_len**2*h*10**-27*N0*C0*10**3))

alpha_theta = np.deg2rad(69)                        # Average degree between alpha-catenin and F-actin
theta_cons = np.cos(alpha_theta)*np.sin(alpha_theta)# A constant which is derived from the angle theta
k6r_central = alpha_dis(0.0)                        # dissociation constant on central area

upper_angle_limit = np.deg2rad(270)
lower_angle_limit = np.deg2rad(90)      


# -------------------------------------------------------------------------------------------------
# MESH
# -------------------------------------------------------------------------------------------------

    # ___________________________________________________________________________________
    #                              Initializing Node File

# Inputing node file
node_file = open(nodeinput, 'r')

# Reading the values of the first line
num_nodes, num_coords, num_attributes, boundary_markers = node_file.readline().split()

# Converting the numbers to integers
num_nodes = int(num_nodes)
num_coords = int(num_coords)
num_attributes = int(num_attributes)
boundary_marker = int(boundary_markers)

# Creating variables to store node number & boundary marker
NodeNums = [[0 for m in range(2)] for n in range(num_nodes)]

# Creating variable to store x, y and z coordinates
NodeCoords = [[0 for m in range(num_coords)] for n in range(num_nodes)]

# Reading data from nodefiler
for i in range(num_nodes):
    NodeNums[i][0], NodeCoords[i][0], NodeCoords[i][1], NodeNums[i][1] = node_file.readline().split()
    # Converting from string to int/float
    NodeNums[i][0] = int(NodeNums[i][0])
    NodeNums[i][1] = int(NodeNums[i][1])
    NodeCoords[i][0] = float(NodeCoords[i][0])
    NodeCoords[i][1] = float(NodeCoords[i][1])
    
# The number of Boundary Nodes
Num_BoundaryNodes = [row[1] for row in NodeNums].count(1)

# Closing the file
node_file.close()

# ___________________________________________________________________________________
#                              Initializing Element File

# Inputing element file
eleminput = eleminput.replace("\r","")
elem_file = open(eleminput, 'r')

# Reading the values of the first line
num_elements, nodes_per_ele, ele_attributes , extra_zero = elem_file.readline().split()

# Converting the numbers to integers
num_elements = int(num_elements)
nodes_per_ele = int(nodes_per_ele)
ele_attributes = int(ele_attributes)

# Creating variable to store element map
ElemMap = [[0 for x in range(3)] for y in range(num_elements)]
Elemindex = [[0 for m in range(1)] for n in range(num_elements)]
# Reading data from elemfile
for i in range(num_elements):

    Elemindex[i][0], ElemMap[i][0], ElemMap[i][1], ElemMap[i][2] = elem_file.readline().split()

    # Converting from string to int
    Elemindex[i][0] = int(Elemindex[i][0])
    ElemMap[i][0] = int(ElemMap[i][0])
    ElemMap[i][1] = int(ElemMap[i][1])
    ElemMap[i][2] = int(ElemMap[i][2])

# Closing the file
elem_file.close()


total_side_len = side_len*(x_len-10)         # Side length of the square
radius_nm = total_side_len/2.0             # radius of the circle in nm
num_sites = x_len * y_len                   # The number of lattices on the surface

gamma_t = math.exp(g_trans/(k*T))
ass_p_t = gamma_t / (1 + gamma_t)   # Association probability trans
dis_p_t = 1 / (1 + gamma_t)         # Dissociation probability trans

gamma_c = math.exp(g_cis/(k*T))
ass_p_c = gamma_c / (1 + gamma_c)   # Association probability cis
dis_p_c = 1 / (1 + gamma_c)         # Dissociation probability cis

num_mono = np.zeros((t_total_step // 100+1, 1), dtype=int)
num_trans = np.zeros((t_total_step // 100+1, 1), dtype=int)
num_cis = np.zeros((t_total_step // 100+1, 1), dtype=int)


#Initialize position and orientations of M1 and M2
#col_0: index, col_1: monomer, col_2: trans, col_3: cis
#col_4: orientation(N: 5, E: 6, S: 7, W: 8),
#col_5: x, col_6: y, col_7: occupied
#col_8: E_trans_actin bond (bound == 1, free == 0)

# M1[:, 0] = range(num_sites) + np.ones(num_sites)



M1 = [[0 for x in range(num_sites)] for y in range(9)]  # Initialize M1 matrix
M2 = [[0 for x in range(num_sites)] for y in range(9)]  # Initialize M2 matrix

# assign index and x,y positions
for i in range(x_len*y_len):
    M1[0][i] = i
    M1[5][i] = i % x_len
    M1[6][i] = int(math.floor(i/x_len))

M2[0][:], M2[5][:], M2[6][:] = M1[0][:], M1[5][:], M1[6][:]

# Find all lattice points within the circle
par_x, par_y = [0 for x in range(num_sites)],[0 for x in range(num_sites)]
if x_len%2 != 0:
    par_x = (np.array(M1[5][:])-x_len//2)*side_len
    par_y = (np.array(M1[6][:])-y_len//2)*side_len
else:
    par_x = (np.array(M1[5][:])-(x_len//2-0.5))*side_len
    par_y = (np.array(M1[6][:])-(y_len//2-0.5))*side_len
dist_center = np.sqrt(par_x**2+par_y**2)
incircle = np.array(np.where(dist_center<radius_nm))        # Indexes of lattice points within the circle
incircle_temp = np.zeros(x_len*y_len)
incircle_temp[incircle] = 1

# edit properties of occupied sites
num_par = int(Ecad_conc*incircle.size)
incircle_list = list(map(int,incircle.ravel().tolist()))                  # turn it into a list
sites_index = range(x_len*y_len)
occupied_site_M1 = [select(incircle_list) for x in range(num_par)]  # randomly select num_par occupied sites,M1
sites_index = range(x_len*y_len)
occupied_site_M2 = [select(incircle_list) for x in range(num_par)]  # randomly select num_par occupied sites,M2
for i in range(num_par):
    M1[1][occupied_site_M1[i]] = 1
    M1[4][occupied_site_M1[i]] = random.randint(5, 8)
    M1[7][occupied_site_M1[i]] = 1
    
    M2[1][occupied_site_M2[i]] = 1
    M2[4][occupied_site_M2[i]] = random.randint(5, 8)
    M2[7][occupied_site_M2[i]] = 1

M1 = np.transpose(np.array(M1))
M2 = np.transpose(np.array(M2))
onedge = M1[np.logical_and((dist_center<radius_nm),(dist_center>cortex_rad*1000)),0]
incenter = M1[(dist_center<cortex_rad*1000),0]
area_onedge = np.pi*((radius_um)**2-(cortex_rad)**2)
area_incenter = np.pi*cortex_rad**2

col_index = np.array([1, 2, 3, 4, 7])
#----------------------------------------------------------------------------------
# Match lattices with elements of the mesh. Find index of elements surrounding a node  
#----------------------------------------------------------------------------------

boundTol = 1e-12
matchedlist = []    # Elements with corresponding lattice points      
Area_all = []       # Store area of each element
matchedfile = open('MatchLatticeElement_7.0.txt','w')
Areafile = open('AreaAll_7.0.txt','w')
# loop through all elements
for m in range(num_elements):
#    matchedfile.write(str(m+1)+'    ')
    templist = []
    # assign x,y of the three nodes of the element to new variables 
    node1_x, node1_y = NodeCoords[ElemMap[m][0]-1][0], NodeCoords[ElemMap[m][0]-1][1]
    node2_x, node2_y = NodeCoords[ElemMap[m][1]-1][0], NodeCoords[ElemMap[m][1]-1][1]
    node3_x, node3_y = NodeCoords[ElemMap[m][2]-1][0], NodeCoords[ElemMap[m][2]-1][1]
    Area = area(node1_x, node1_y, node2_x, node2_y, node3_x, node3_y)   # Area of the element
    # loop through all lattice points
    for n in range(num_sites):
        lat_x = par_x[n]/1000.0     # x of the lattice point
        lat_y = par_y[n]/1000.0     # y of the lattice point

        # Area Calculation, area is a module
        Area1 = area(lat_x, lat_y, node2_x, node2_y, node3_x, node3_y)      # Area of the triangle with vertexs P,B,C
        Area2 = area(node1_x, node1_y, lat_x, lat_y, node3_x, node3_y)      # Area of the triangle with vertexs A,P,C
        Area3 = area(node1_x, node1_y, node2_x, node2_y, lat_x, lat_y)      # Area of the triangle with vertexs A,B,P
        
        
        # if P is in the Element
        if abs(Area-Area1-Area2-Area3)<boundTol:
            templist.append(n) 
    Area_all.append(Area)
    Areafile.write(str(Area))
    Areafile.write('\n')
    
    if templist != []:
        matchedfile.write('   '.join(map(str,templist)))
    matchedfile.write('\n')             # file written
    matchedlist.append(templist)        # list stored
matchedfile.close()                     # close file
Areafile.close()
Area_all = np.array(Area_all)
# Check duplication
# result shows that the program successfully match lattice with nodes except when 
# the position of a lattice point is exactly the same as the position of a node 
seen = set()
duplicated = []
for x in matchedlist:
    for y in x:
        if y not in seen:
            seen.add(y)
        else:
            duplicated.append(y)
#------------------------------------------------------------------------------
# generate a list, find nearby elements of a node

ElemMap = np.array(ElemMap)
nodenearbyelem = []
nodenearbyelemfile = open('MatchNodeElement_7.0.txt','w')
for i in range(num_nodes):
    nodenearbyelemfile.write(str(i+1)+'    ')
    node_elem = np.array(np.where(ElemMap==(i+1)))      # find all elements nearby the node
    nodenearbyelem.append(node_elem[0,:]+1)             # append the elements to the list
    
    nodenearbyelemfile.write('   '.join(map(str,node_elem[0,:]+1)))
    nodenearbyelemfile.write('\n')

nodenearbyelemfile.close()
