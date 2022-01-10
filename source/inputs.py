# Author: Qilin Yu
# Date: April 26 2021



import math
import numpy as np
from forcebased_k_new import alpha_dis

####################################
# GEOMETRIC INPUTS
####################################

#Inputing the mesh node and ele files

nodeinput, eleminput = './ellipse.3.node','./ellipse.3.ele'
nodeincenter, elemincenter = './incenter_fine.2.node','./incenter_fine.2.ele'

####################################
# SIMULATION RUNTIME INPUTS
####################################
#Dimension of the model
number_D = 2

#simulation rd time parameters - PDE Timestep, overall simulation time (seconds), start time, update frequency of ODE, ODE Timestep
Tstep, Total_time, startT, freq_update, ODE_Timestep    = 1.0, 5, 0.0, 100, 0.00001

#simulation mech load parameters - Load Step
load_steps = 1

#results output frequency
outputfreq = 1

# record_freq
record_freq = 1000

####################################
# REACTION DIFFUSION INPUTS
####################################

#Inputing the cellml model of Holmes wave-pinning actin polymerization signalling model
Actin_Ecad = './Actin_Ecad_4.0_1.cellml'
 
# E_trans Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_E_trans, Dx_E_trans, Dy_E_trans = 0.0,0.0,0.0

# actin Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_actin_1, Dx_actin, Dy_actin  = 5119,0.1,0.1
init_actin_2 = 5119

# Branched Actin Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_actin_a_1, Dx_actin_a, Dy_actin_a = 2408,0.0,0.0
init_actin_a_2 = 2408

# Rac Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_Rac_1, Dx_Rac, Dy_Rac = 202.5,0.1,0.1
init_Rac_2 = 202.5

# Rac_a Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_Rac_a_1, Dx_Rac_a, Dy_Rac_a = 22.5,0.0,0.0
init_Rac_a_2 = 22.5

# Rho Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_RhoA_1, Dx_RhoA, Dy_RhoA = 83.7,0.1,0.1
init_RhoA_2 = 83.7

# Rho_a Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_RhoA_a_1, Dx_RhoA_a, Dy_RhoA_a = 9.3,0.0,0.0
init_RhoA_a_2 = 9.3

# myo Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_myo_1, Dx_myo, Dy_myo = 2976.67,0.1,0.1
init_myo_2 = 2976.67

# myo_a Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_myo_a_1, Dx_myo_a, Dy_myo_a = 2.0,0.0,0.0
init_myo_a_2 = 2.0

# E_trans_actin Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_E_trans_actin_1, Dx_E_trans_actin, Dy_E_trans_actin = 0.0,0.0,0.0
init_E_trans_actin_2 = 0.0

# myo_actin Initial Conc. (molecules/um^2) / Diffusivities in x (um^2/s)/ Diffusivities in y (um^2/s)
init_myo_actin_1, Dx_myo_actin, Dy_myo_actin = 190,0.0,0.0
init_myo_actin_2 = 190

###############################################################################
# Kinetic Constant
###############################################################################

k0_1  = 0.0           #k0
k1_1  = 0.00001        #k1
k2_1  = 0.0006272	    #k2
k3_1  = 0.01		    #k3
k3r_1 = 0.09		    #k3r
k4_1  = 0.03		    #k4
k5_1  = 0.001		    #k5
k5r_1 = 10.0   		#k5r
k6_1  = 0.05		    #k6	
k7_1  = 0.136         #k7
k7r_1 = 0.1           #k7r
k8_1  = 0.26          #k8		
k9_1  = 0.2           #k9		
k9r_1 = 37.6          #k9r	

k0_2  = 0.0           #k0
k1_2  = 0.00001        #k1
k2_2  = 0.0006272	    #k2
k3_2  = 0.01		    #k3
k3r_2 = 0.09		    #k3r
k4_2  = 0.03		    #k4
k5_2  = 0.001		    #k5
k5r_2 = 10.0   		#k5r
k6_2  = 0.05		    #k6	
k7_2  = 0.136         #k7
k7r_2 = 0.1           #k7r
k8_2  = 0.26          #k8		
k9_2  = 0.2           #k9		
k9r_2 = 37.6          #k9r	

####################################
# E-cad lattice model inputs
####################################

# E-cad concentration (Ecad_conc)
Ecad_conc_M1 = 0.04
Ecad_conc_M2 = 0.04

# lattice side length (side_len)
side_len = 7.0

# E-cad diffusivity	(diffusivity)
diffusivity = 28000.0

# The number of lattices alone the x axis	(x_len)
x_len = 286

# The number of lattices alone the y axis	(y_len)
y_len = 286

# Free energy change for trans dimer	(g_trans/kT) (KD)
g_trans = 6

# Free energy change for cis interaction	(g_cis/kT) (KD)
g_cis = 2

# Free energy change for trans dimer for higher g_trans region (only useful if g_trans is not uniform on the whole region)	(g_trans/kT) (KD)
g_trans_ringright = 6

####################################
# CONSTANTS INPUTS
####################################

# Radius of the contact in micron radius_um
radius_um = 0.5

# Cortex radius
cortex_rad = 0.45	
									
# Avogadro's number
N0 = 6.0221409e23
  
# Concentration constant C0                                
C0 = 1.0

# gas constant R
R = 8.3145
 
# Boltzmann constant k                                     
k = 1.38064852e-23
  
# Temprature T                               
T = 298.15

# height nm h                                         
h = 10.0
                                            
# F_cort (Cortical tension (N/m or kg/s^2) 900pN/um or 9e-4 N/m) 
F_cort = 4000
       
# cort_thick (Thickness of the cortical actin (um))  
cort_thick = 0.05

# con_cons (Concentration constants) 
con_cons = 180.667

# Fila_len (Estimated Filament length)  
Fila_len = 3.0

# real_radius (radius of the real intercellular contact) 
real_radius = 5.0

# F_myo (force generated by each myosin head) 
F_myo = 4.0
                                         
# L_actin_mono (length of each G-actin/actin monomer, 0.0027um or 2.7nm) 
L_actin_mono = 0.0027

####################################
# CELL ATTACHMENT MODE
####################################

# mode
mode = 1

# fixed_angle
fixed_angle = 30

# init_angle, initial theta calculated from Hertz equation.
init_angle = 1.825  
                    
# time to fix the contact angle (unit: second)
T_fix_angle = 60

# time to change local g_trans (unit: second)
t_increase_gtrans = 60

# If E-cad Trans dimer is affected by Force (Mode1: Yes, Mode2: No)
Force_sensitive = 1

#-------------------------------------------------------------------------------------------------
# Constants used in the simulation
#-------------------------------------------------------------------------------------------------

# Time related parameters calculation
deltat = side_len ** 2 / (4 * diffusivity)          # lenght of each time step (s)
t_total_step = int(math.ceil(Total_time/deltat))    # total number of steps
Total_time_actual = t_total_step*deltat             # actual total time (s)

alpha_theta = np.deg2rad(69)                        # Average degree between alpha-catenin and F-actin
theta_cons = np.cos(alpha_theta)*np.sin(alpha_theta)# A constant which is derived from the angle theta
k5r_central = alpha_dis(0.0)                        # dissociation constant on central area

upper_angle_limit = np.deg2rad(270)
lower_angle_limit = np.deg2rad(90)    

#-------------------------------------------------------------------------------------------------
# Define Fields' ID
#-------------------------------------------------------------------------------------------------
store_coeff = 1.0

numberOfXi = 2

CoordinateSystemUserNumber  = 10
RegionUserNumber            = 20
BasisUserNumber             = 30
PressureBasisUserNumber     = 40
GeneratedMeshUserNumber     = 50
MeshUserNumber              = 60
DecompositionUserNumber     = 70
GeometricFieldUserNumber    = 80

rd_problemUserNumber                = 100

E_trans_EquationsSetFieldUserNumber = 111
E_trans_EquationsSetUserNumber      = 121
E_trans_FieldUserNumber             = 131
E_trans_MaterialsFieldUserNumber    = 141
E_trans_SourceFieldUserNumber       = 151

actin_1_EquationsSetFieldUserNumber   = 161
actin_1_EquationsSetUserNumber        = 171
actin_1_FieldUserNumber               = 181
actin_1_MaterialsFieldUserNumber      = 191
actin_1_SourceFieldUserNumber         = 201

actin_a_1_EquationsSetFieldUserNumber = 211
actin_a_1_EquationsSetUserNumber      = 221
actin_a_1_FieldUserNumber             = 231
actin_a_1_MaterialsFieldUserNumber    = 241
actin_a_1_SourceFieldUserNumber       = 251

Rac_1_EquationsSetFieldUserNumber     = 481
Rac_1_EquationsSetUserNumber          = 491
Rac_1_FieldUserNumber                 = 501
Rac_1_MaterialsFieldUserNumber        = 511
Rac_1_SourceFieldUserNumber           = 521

Rac_a_1_EquationsSetFieldUserNumber   = 531
Rac_a_1_EquationsSetUserNumber        = 541
Rac_a_1_FieldUserNumber               = 551
Rac_a_1_MaterialsFieldUserNumber      = 561
Rac_a_1_SourceFieldUserNumber         = 571

RhoA_1_EquationsSetFieldUserNumber     = 581
RhoA_1_EquationsSetUserNumber          = 591
RhoA_1_FieldUserNumber                 = 601
RhoA_1_MaterialsFieldUserNumber        = 611
RhoA_1_SourceFieldUserNumber           = 621

RhoA_a_1_EquationsSetFieldUserNumber   = 631
RhoA_a_1_EquationsSetUserNumber        = 641
RhoA_a_1_FieldUserNumber               = 651
RhoA_a_1_MaterialsFieldUserNumber      = 661
RhoA_a_1_SourceFieldUserNumber         = 671

myo_1_EquationsSetFieldUserNumber         = 681
myo_1_EquationsSetUserNumber              = 691
myo_1_FieldUserNumber                     = 701
myo_1_MaterialsFieldUserNumber            = 711
myo_1_SourceFieldUserNumber               = 721

myo_a_1_EquationsSetFieldUserNumber       = 731
myo_a_1_EquationsSetUserNumber            = 741
myo_a_1_FieldUserNumber                   = 751
myo_a_1_MaterialsFieldUserNumber          = 761
myo_a_1_SourceFieldUserNumber             = 771

E_trans_actin_1_EquationsSetFieldUserNumber   = 781
E_trans_actin_1_EquationsSetUserNumber        = 791
E_trans_actin_1_FieldUserNumber               = 801
E_trans_actin_1_MaterialsFieldUserNumber      = 811
E_trans_actin_1_SourceFieldUserNumber         = 821

myo_actin_1_EquationsSetFieldUserNumber       = 831
myo_actin_1_EquationsSetUserNumber            = 841
myo_actin_1_FieldUserNumber                   = 851
myo_actin_1_MaterialsFieldUserNumber          = 861
myo_actin_1_SourceFieldUserNumber             = 871

k0_1_FieldUserNumber                      = 2001
k1_1_FieldUserNumber                      = 2011
k2_1_FieldUserNumber                      = 2021
k3_1_FieldUserNumber                      = 2031
k3r_1_FieldUserNumber                     = 2041
k4_1_FieldUserNumber                      = 2071
k5_1_FieldUserNumber                      = 2101
k5r_1_FieldUserNumber                     = 2111
k6_1_FieldUserNumber                      = 2121
k7_1_FieldUserNumber                      = 2131
k7r_1_FieldUserNumber                     = 2141
k8_1_FieldUserNumber                      = 2151
k9_1_FieldUserNumber                      = 2161
k9r_1_FieldUserNumber                     = 2171

actin_2_EquationsSetFieldUserNumber   = 162
actin_2_EquationsSetUserNumber        = 172
actin_2_FieldUserNumber               = 182
actin_2_MaterialsFieldUserNumber      = 192
actin_2_SourceFieldUserNumber         = 202

actin_a_2_EquationsSetFieldUserNumber = 212
actin_a_2_EquationsSetUserNumber      = 222
actin_a_2_FieldUserNumber             = 232
actin_a_2_MaterialsFieldUserNumber    = 242
actin_a_2_SourceFieldUserNumber       = 252

Rac_2_EquationsSetFieldUserNumber     = 482
Rac_2_EquationsSetUserNumber          = 492
Rac_2_FieldUserNumber                 = 502
Rac_2_MaterialsFieldUserNumber        = 512
Rac_2_SourceFieldUserNumber           = 522

Rac_a_2_EquationsSetFieldUserNumber   = 532
Rac_a_2_EquationsSetUserNumber        = 542
Rac_a_2_FieldUserNumber               = 552
Rac_a_2_MaterialsFieldUserNumber      = 562
Rac_a_2_SourceFieldUserNumber         = 572

RhoA_2_EquationsSetFieldUserNumber     = 582
RhoA_2_EquationsSetUserNumber          = 592
RhoA_2_FieldUserNumber                 = 602
RhoA_2_MaterialsFieldUserNumber        = 612
RhoA_2_SourceFieldUserNumber           = 622

RhoA_a_2_EquationsSetFieldUserNumber   = 632
RhoA_a_2_EquationsSetUserNumber        = 642
RhoA_a_2_FieldUserNumber               = 652
RhoA_a_2_MaterialsFieldUserNumber      = 662
RhoA_a_2_SourceFieldUserNumber         = 672

myo_2_EquationsSetFieldUserNumber         = 682
myo_2_EquationsSetUserNumber              = 692
myo_2_FieldUserNumber                     = 702
myo_2_MaterialsFieldUserNumber            = 712
myo_2_SourceFieldUserNumber               = 722

myo_a_2_EquationsSetFieldUserNumber       = 732
myo_a_2_EquationsSetUserNumber            = 742
myo_a_2_FieldUserNumber                   = 752
myo_a_2_MaterialsFieldUserNumber          = 762
myo_a_2_SourceFieldUserNumber             = 772

E_trans_actin_2_EquationsSetFieldUserNumber   = 782
E_trans_actin_2_EquationsSetUserNumber        = 792
E_trans_actin_2_FieldUserNumber               = 802
E_trans_actin_2_MaterialsFieldUserNumber      = 812
E_trans_actin_2_SourceFieldUserNumber         = 822

myo_actin_2_EquationsSetFieldUserNumber       = 832
myo_actin_2_EquationsSetUserNumber            = 842
myo_actin_2_FieldUserNumber                   = 852
myo_actin_2_MaterialsFieldUserNumber          = 862
myo_actin_2_SourceFieldUserNumber             = 872

k0_2_FieldUserNumber                      = 2002
k1_2_FieldUserNumber                      = 2012
k2_2_FieldUserNumber                      = 2022
k3_2_FieldUserNumber                      = 2032
k3r_2_FieldUserNumber                     = 2042
k4_2_FieldUserNumber                      = 2072
k5_2_FieldUserNumber                      = 2102
k5r_2_FieldUserNumber                     = 2112
k6_2_FieldUserNumber                      = 2122
k7_2_FieldUserNumber                      = 2132
k7r_2_FieldUserNumber                     = 2142
k8_2_FieldUserNumber                      = 2152
k9_2_FieldUserNumber                      = 2162
k9r_2_FieldUserNumber                     = 2172

cellmlUserNumber                    = 3200
cellmlModelsFieldUserNumber         = 3210
cellmlStateFieldUserNumber          = 3220
cellmlIntermediateFieldUserNumber   = 3230
cellmlParametersFieldUserNumber     = 3240
  
