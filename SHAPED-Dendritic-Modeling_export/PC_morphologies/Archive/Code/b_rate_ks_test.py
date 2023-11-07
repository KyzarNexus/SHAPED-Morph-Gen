#%%
# Purpose: Use pt1 from ppfem72 to create branching rate distributions. After which, perform a ks test between the generated and original rates. Record max D score and p value for each. 

import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import itertools as iter
import pickle
import pandas as  pd
import plotly 
import plotly.graph_objs as go

def read_file(file_name): # Reads a point process file and returns a point_process_list

    point_process_list = []

    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()
    
    for i in range(0, len(lines)):
        parsed_line = lines[i].split(' ') 
        
        point_process_list.append([0] + list(map(int,parsed_line)))
    
    return point_process_list

def external_rescale_pp(n,rescale_flag = 1):
    # Rescales a point process independently of the scope of a class. 
    max_N = len(n[0])
    max_R = len(n)
    if rescale_flag == 1:
    # Rescale point processes
        for r in range(0, max_R):
            # Calculate scaling factor    
            list_of_indices = [i for i, x in enumerate(n[r]) if x == 1]
            termination = list_of_indices[-1]
            scaling_factor = (max_N-1)/termination
            
            # Rescaled point process indices
            list_of_indices_scaled = [int(x*scaling_factor) for x in list_of_indices]
            new_point_process = [0]*(max_N)
            for idx in range(0, len(list_of_indices_scaled)):
                new_point_process[list_of_indices_scaled[idx]] = 1
            new_point_process[-1] = 0
            new_point_process[-2] = 0
            n[r] = new_point_process
    return n

def ks_w_trt(n, part2_flag): # K-S Test With Time Rescaling Theorem Implemented
    # Validation of the model via KOLMOGOROV SMIRNOV GOODNESS OF FIT TEST.
    # Input: n (point process file)
    if part2_flag == 1:
        rescaled_list = []
        first_match = 0
        C = len(n)
        T = len(n[0])
        for c in range(0, C):
            # if (c == 0) or (n[c][first_match] != n[c+1][first_match]):
            if 1 == 1:
                k = 0
                old_lambda = 0
                while k < T:
                    old_lambda = old_lambda + lamb[k]
                    if n[c][k] == 1:
                        if first_match == 0:
                            first_match = k
                        tau_k = old_lambda
                        old_lambda = 0
                        # val = 1 - math.exp(-tau_k)
                        # if val not in rescaled_list:
                        rescaled_list.append(1 - math.exp(-tau_k))
                    k = k + 1
            # #only for rescaled
            # if rescale_flag == 1:        
            #     rescaled_list.append(1-math.exp(-old_lambda))

        b_n = len(rescaled_list)
        b_k = []
        for jj in range(0, b_n):
            b_k.append((jj + .05)/b_n)
        
        rescaled_list.sort()
        
        b_k_upper = [b_l + 1.36/math.sqrt(b_n) for b_l in b_k]
        b_k_lower = [b_l - 1.36/math.sqrt(b_n) for b_l in b_k]
        plt.figure()
        plt.plot(rescaled_list, b_k)
        plt.plot(b_k, b_k, color='black')
        plt.plot(b_k, b_k_upper, linestyle='dashed', color='black')
        plt.plot(b_k, b_k_lower, linestyle='dashed', color='black')
        #        while n[c][k] != 1:
        #            tau_k = n[c][k]
    else:
        return 0

# Declaring File Names
file_pp = '..\\..\\Data\\Metrics\\amaral_point_process.txt'
file_pp_apical = '..\\..\\Data\\Metrics\\amaral_point_process_apical.txt'
file_pp_basal = '..\\..\\Data\\Metrics\\amaral_point_process_basal.txt'

file_con = '..\\..\\Data\\Metrics\\amaral_conditionals.txt'
file_con_apical = '..\\..\\Data\\Metrics\\amaral_conditionals_apical.txt'
file_con_basal = '..\\..\\Data\\Metrics\\amaral_conditionals_basal.txt'

# # Import Sample Distributions
# Basal
n = read_file(file_pp_basal); n = external_rescale_pp(n); T = len(n[0])
max_N = len(n[0]); max_R = len(n)
# # Import Optimization Dictionary
pickle_in = open("optPicFull.pickle","rb"); opDict = pickle.load(pickle_in)
keys = opDict.keys()
final_keys = []
for key in keys: # Grabbing Final Values
    final_keys.append(tuple(opDict[key][-1]))
point_process_dict = dict({key : None for key in final_keys})

# # Generating branching list for each key and adding to point_process_dict
# p_max: [Description] ; 
p_max = 700; 
for key in final_keys:
    rho = key[0]; var = key[1]; beta = key[2]; mu = key[3]

    coef1 = 0.8; tau = 200
    
    max_sv = 10; sv_tau = 0.008; sv_shift = 0
    sv_filterx = np.linspace(0, T-1, T)
    sv_filter = [max_sv*np.exp(-((x-sv_shift)*sv_tau)) + 1 for x in sv_filterx]

    lamb = [0]*len(x_smooth)
    for kdx in range(0, len(x_smooth)):
        lamb[kdx] = np.exp(mu + beta*x_smooth[kdx])/sv_filter[kdx]
    lamb = [np.exp(mu + beta*x_s)*coef1 for x_s in x_smooth]
    lamb_tmp = list(lamb)
    lambda_max = max(lamb)
    lamb2 = list(lamb)
    k = 0
    max_poisson = [[]]*p_max
    T_max = max_N
    
    delta = 1
    rate_time = np.linspace(0, int(len(lamb))-1, int(len(lamb)/delta))
    rate_res = delta
    
    L = intn_flag = 0
    while poisson_flag == 0:
        max_poisson[i] = [0]*(max_N+1)    
        k = 0
        while k < T_max-1:
            u_k = np.random.uniform(0, 1)
            w_k = -math.log(u_k)/(lambda_max)
            k = k + 1 + int(w_k)
            if k >= T_max:
                max_poisson[i][-1] = 0
            else:
                max_poisson[i][k] = 1
        
        for j in range(0, len(max_poisson[i])):
            if max_poisson[i][j] == 1:
                p = lamb2[j]/(questionable*max(lamb2))
                b = np.random.binomial(1, p)
                if b == 0:
                    max_poisson[i][j] = 0
                else:
                    aa = 0
                    if (j+L) > len(lamb2):
                        lamb2[j:] *= refract[:len(lamb2)-j]
                    else:
                        lamb2[j:j+L] *= refract
        poisson_sum = sum(max_poisson[i])
        if poisson_sum < 1 or poisson_sum > 6:
            poisson_flag = 1
        else:
            poisson_flag = 0
    poisson_sums = [0]*len(max_poisson)
    for i in range(len(max_poisson)):
        poisson_sums[i] = sum(max_poisson[i])(10*tau/rate_res)
    refract = 1-np.exp(-rate_time[:L]/tau)
    questionable = 1
    weighting = 1
    for i in range(0, p_max):
        lamb2 = list(lamb)
        max_poisson[i] = [0]*(max_N+1)   
        k = 0  
    
    # Return lamb2
    
