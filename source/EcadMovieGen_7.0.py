# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 10:25:34 2019
Using this program to generate E-cad movies to visulize E-cads' movements
@author: qiliny
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.animation as animation
# load data
M1_data = np.loadtxt('M1_data.txt')
M2_data = np.loadtxt('M2_data.txt')
num_par_index = 0
num_par = M1_data[0,0]
x_len, y_len = 286, 286
M1_data = M1_data[0:int(1+len(M1_data)//num_par*(num_par+1)),:]
M2_data = M2_data[0:int(1+len(M1_data)//num_par*(num_par+1)),:]
num_frames = int(len(M1_data)//num_par)

M1_data = np.c_[M1_data,M1_data[:,0]%x_len]
M1_data = np.c_[M1_data,np.int_(np.floor(M1_data[:,0]/x_len))]
M2_data = np.c_[M2_data,M2_data[:,0]%x_len]
M2_data = np.c_[M2_data,np.int_(np.floor(M2_data[:,0]/x_len))]

x_y_pos = np.zeros((x_len*y_len,3),dtype = int)
x_y_pos[:,0] = np.arange(x_len**2)
x_y_pos[:,1] = x_y_pos[:,0]%x_len
x_y_pos[:,2] = x_y_pos[:,0]//y_len

######################################################
#  movie initialization 
######################################################
dpi=100     # dots per inch, which is the resolution
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter([],[], s=3, c='k')      # all Lattices
ax.axis('off')                              
M1data = ax.scatter([],[], s=1, c='r')   # M1 monomers
M2data = ax.scatter([],[], s=1, c='g')   # M2 monomers
transdata = ax.scatter([],[], s=2, c='b')   # Trans dimers
#fixlatticedata = ax.scatter([],[], s=1, c='tab:gray')   # fixed sites
#E_trans_actindata = ax.scatter([],[], s=30, c='b')   # Trans dimers

ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax.tick_params(axis='y', which='both', bottom=False, top=False, labelleft=False)
fig.set_size_inches(15, 15)   
ax.axis([-0.3, x_len+1, -0.3, y_len+1])

def update_img(n,M1,M2):
    global num_par_index
    num_par = int(M1_data[num_par_index,0])
    
    M1_temp = M1[num_par_index+1:num_par_index+num_par+1,:]
    M2_temp = M2[num_par_index+1:num_par_index+num_par+1,:]
    
    M1_mono_or_cis = (M1_temp[:,1]==1) | (M1_temp[:,3]==1) | (M1_temp[:,4]==1)
    M2_mono_or_cis = (M2_temp[:,1]==1) | (M2_temp[:,3]==1) | (M2_temp[:,4]==1)
    M1data.set_offsets(np.c_[M1_temp[np.where(M1_mono_or_cis),6].ravel()+1,M1_temp[np.where(M1_mono_or_cis),7].ravel()+1])
    
    M2data.set_offsets(np.c_[M2_temp[np.where(M2_mono_or_cis),6].ravel()+1,M2_temp[np.where(M2_mono_or_cis),7].ravel()+1])
    
    transdata.set_offsets(np.c_[M1_temp[np.where(M1_temp[:,2]==1),6].ravel()+1,M1_temp[np.where(M1_temp[:,2]==1),7].ravel()+1])
    
    num_par_index += num_par+1
#    E_trans_actindata.set_offsets(np.c_[M1_temp[np.where(M1_temp[:,4]==1),5].ravel()+1,M1_temp[np.where(M1_temp[:,4]==1),6].ravel()+1])

    return transdata,M1data,M2data#,,E_trans_actindata


ani = animation.FuncAnimation(fig,update_img,frames=num_frames,fargs=(M1_data,M2_data))
writer = animation.writers['ffmpeg'](fps=23)

ani.save('Latt_small.mp4',writer=writer,dpi=dpi)














