'''
DEND_funn_classes.py
Purpose: To act as a centralized repository for class and function definitions to be called by DEND_main.py
In the intermediate stage of organization, file sources will be noted before the initial dump. 

Currently Integrated:           | I - P - F 
angle_noise_updated.py          | X - 
B_measure_tester.py             | X - 
B_measure_conditionals.py       | X - 
IBI observations.py             | X - O - X
point-process_filter_EM_72.py   | X - 
B_morph_170811_2.py             | X - 

Steps:
Initial Integration 
Removing Commented Sections
Resolving Class Disputes
Resolving Function Disputes
Handling I/O from DEND_fun_class.py

Notes:
-Check if some function definitions are redundant. (Example: There should already be a function for dot and cross products in np, hopefully.)
- File progress partially discontinued. See 'DENT_TEST.py'
'''
# Combined Imports/Classes/Functions
# # Combined Imports
from __future__ import division
import os
import random
import math
import pandas
import glob
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
from scipy.ndimage.filters import gaussian_filter
import pickle
from mpl_toolkits import mplot3d
from scipy.ndimage.measurements import label 
import matplotlib as mpl
from ast import literal_eval



# # Combined Constants

# # Combined Classes # #
'''
Status:
Compartment: Partial
Branch: N/C
Bifurcation: N/C
Morphology: Complete
Data: Partial
Stem: Partial
'''
class Compartment:
        
    def __init__(self, line):
        self.compartment_id, self.section_type, self.x, self.y, self.z, self.diameter, self.parent_id = map(float, line)
        self.diameter = 2*self.diameter
        self.parent_compartment = []
        self.daughter_compartment = []
        self.start_coordinates = (0, 0, 0)
        self.end_coordinates = (self.x, self.y, self.z)

        self.bifurcation_flag = 0

        self.branch = []

        self.length = 0
        self.euclidean_from_soma = euclidean_distance(self.start_coordinates, self.end_coordinates)

        # -----------------------------------------------------------------------

        self.start_coordinates = origin
        self.end_coordinates = end_coordinates
        self.diameter = diameter
        self.vector = [end_coordinates[i] - origin[i] for i in range(0, 3)]
        self.length = euclidean_distance(self.start_coordinates, self.end_coordinates)
        self.parent_compartment = []
        self.jitter_phi = 0
        self.jitter_theta = 0
        global c_id
        self.compartment_id = c_id
        c_id = c_id + 1
        self.terminal_flag = 0
        self.bifurcation_flag = 0
        self.stem_flag = 0
        
        #------------------------------------------------------------------------------------------

        self.compartment_id = 1
        self.parent_compartment = []
        self.start_coordinates = start_coordinates
        self.end_coordinates = end_coordinates
        self.diameter = diameter
        self.jitter_theta = 0
        self.jitter_phi = -90
        self.counter = 1
    
    def set_length(self):
        self.length = euclidean_distance(self.start_coordinates, self.end_coordinates)

class Branch:
    def __init__(self, compartment_list):
        global delta
        self.start_compartment = compartment_list[0]
        self.end_compartment = compartment_list[-1]

        self.start_coordinates = (self.start_compartment.start_coordinates)
        self.end_coordinates = (self.end_compartment.end_coordinates)

        self.vector = [self.end_coordinates[i] - self.start_coordinates[i] for i in range(0, 3)]

        self.compartment_list = compartment_list

        self.parent_bifurcation = []
        self.daughter_bifurcation = []
     

        self.fragmentation = len(compartment_list) - 1
        self.euclidean = euclidean_distance(self.start_coordinates, self.end_coordinates)
        self.pathlength = sum(compartment.length for compartment in compartment_list)
        self.euclidean_from_soma = euclidean_distance(self.start_coordinates, (0, 0, 0))
        self.pathlength_to_soma = self.pathlength
        
        if self.pathlength/delta < 1:
            print('delta too big', self.pathlength)
        self.point_process_branch = [0]*(int(self.pathlength/delta)-1) + [1]
        self.point_process = [0]*(int(self.pathlength/delta)-1) + [1]
        
        self.start_diameter = compartment_list[1].diameter
        self.end_diameter = self.end_compartment.diameter

        self.taper2 = (self.start_diameter - self.end_diameter)/self.start_diameter

        self.branch_order = 0
        self.branch_order_1 = 0
        self.branch_order_2 = 0
        self.branch_order_3 = 0
        self.branch_order_4 = 0

            # last_compartment = self.start_coordinates
        for compartment in compartment_list:
            if compartment != compartment_list[0]:
                compartment.branch  = self
                # self.pathlength = self.pathlength + math.sqrt((compartment.x-last_compartment[0])**2 + (compartment.y-last_compartment[1])**2 + (compartment.z-last_compartment[2])**2)
                # last_compartment = (compartment.x, compartment.y, compartment.z)
        if self.pathlength == 0:
            self.contraction = 0
        else:
            self.contraction = self.euclidean/self.pathlength

    def set_pathlength_to_soma(self):
        if self.daughter_bifurcation != []:
            for daughter in self.daughter_bifurcation.daughter_branch:
                daughter.pathlength_to_soma = daughter.pathlength + self.pathlength_to_soma
                daughter.point_process = self.point_process + daughter.point_process_branch
                daughter.set_pathlength_to_soma() # Recursion Error # Previously Removed for Testing

    def set_parent_bifurcation(self, parent_bifurcation):
        self.parent_bifurcation = parent_bifurcation
        parent_bifurcation.set_daughter_branch(self)

    def set_daughter_bifurcation(self, daughter_bifurcation):
        self.daughter_bifurcation = daughter_bifurcation
        daughter_bifurcation.set_parent_branch(self)
#-----------------------------------------------------------------------------
class Branch:
    def __init__(self, bifurcation, vector_remote, vector_local, normal, start_diameter):
        self.parent_bifurcation = bifurcation
        self.daughter_bifurcation = []
        self.normal = normal
        self.branch_order = bifurcation.branch_order
        self.branch_order_1 = bifurcation.branch_order_1
        self.branch_order_2 = bifurcation.branch_order_2
        self.branch_order_3 = bifurcation.branch_order_3
        self.branch_order_4 = bifurcation.branch_order_4
        
        flag = 0
        self.termination = 0
        self.pathlength = 0
        self.pathlength_proposal = 0
        self.stem_flag = 0

        while flag == 0 and self.termination == 0:
            self.termination = 0
            # self.pathlength = draw_value('pathlength')

            if bifurcation.pathlength <= 60:
                self.pathlength = draw_value('ibf_pathlength_1')
                if random.random() < 0.2 or self.branch_order_1 >= 4:
                    self.pathlength = draw_value('terminal_pathlength_1')
            elif bifurcation.pathlength <= 150:
                self.pathlength = draw_value('ibf_pathlength_2')
                if random.random() < 0.5 or self.branch_order_2 >= 3:
                    self.termination = 1                
            elif bifurcation.pathlength <= 250:
                self.pathlength = draw_value('ibf_pathlength_3')
                if random.random() < 0.8 or self.branch_order_3 >= 2:
                    self.termination = 1
            else:
                self.pathlength = draw_value('ibf_pathlength_4')
                if random.random() < 0.9 or self.branch_order_4 >= 2:
                    self.termination = 1

            if self.pathlength > 0:
                flag = 1
            # if self.pathlength + bifurcation.pathlength > neuron_pathlength:
            #     flag = 0
            #     self.termination = 1

            # if random.random() < 0.4:
            #     self.pathlength = random.gauss(40,10)
            #     self.termination = 1


        if self.branch_order >= max_branch_order:
            self.termination = 1
            # self.pathlength = random.gauss(neuron_pathlength - bifurcation.pathlength, 10)
        if neuron_pathlength - self.pathlength - bifurcation.pathlength < 0:
            self.termination = 1
        # if neuron_pathlength - draw_value('terminal_pathlength') - bifurcation.pathlength < 0:
        #     self.termination = 1
            
        self.start_diameter= start_diameter
        self.taper = -1
        self.threshold_diameter= -1
        while self.taper < 0:
            self.taper = draw_value('ibf_taper')
        self.end_diameter = self.taper*self.start_diameter
        while self.threshold_diameter < 0:
            self.threshold_diameter = draw_value('diam_threshold')
        
        # if self.end_diameter < self.threshold_diameter:
        #     self.termination = 1            
        
        while self.termination == 1:
            while self.taper < 0:
                self.taper = self.taper*2
                
                
            if bifurcation.euclidean_from_soma < 50:
                self.pathlength_proposal = draw_value('terminal_pathlength_1')
                self.termination = 2
            elif bifurcation.euclidean_from_soma < 150:
                self.pathlength_proposal = draw_value('terminal_pathlength_2')
                self.termination = 2
            elif bifurcation.euclidean_from_soma < 250:
                self.pathlength_proposal = draw_value('terminal_pathlength_3')
                self.termination = 2
            elif bifurcation.euclidean_from_soma < 300:
                # self.pathlength_proposal = neuron_pathlength - draw_value('terminal_pathlength_4')
                self.pathlength_proposal = draw_value('terminal_pathlength_4')
                self.termination = 2
            elif bifurcation.euclidean_from_soma > 300:
                self.pathlength_proposal = draw_value('terminal_pathlength_4')
                self.termination = 2
            # self.pathlength_proposal = draw_value('terminal_pathlength')
                
            # if self.pathlength_proposal + bifurcation.pathlength < neuron_pathlength-10:
            #     self.termination = 0
                
            # if self.pathlength_proposal + bifurcation.pathlength > neuron_pathlength + 10:
            #     self.termination = 1

        
            self.pathlength = self.pathlength_proposal
            # if self.pathlength_proposal + bifurcation.pathlength > neuron_pathlength:
            #     self.pathlength = random.gauss(neuron_pathlength - bifurcation.pathlength, 10)
            #     self.termination = 2
            # self.pathlength = self.pathlength_proposal
            # if bifurcation.euclidean_from_soma > 300:
            #     self.pathlength = draw_value('terminal_pathlength_4')
            #     self.termination = 2            

            if self.pathlength + bifurcation.pathlength > 270:
                self.pathlength = random.gauss(neuron_pathlength - bifurcation.pathlength, 10)
            # else:
            #     self.termination = 0
                
            # if random.random() < 0.3 and bifurcation.pathlength <=200:
                
            #     self.pathlength_proposal = draw_value('terminal_pathlength')
            #     if neuron_pathlength - bifurcation.pathlength - self.pathlength_proposal > 0:
            #         self.pathlength = -1
            #         while self.pathlength < 0: 
            #             self.pathlength = random.gauss(neuron_pathlength - bifurcation.pathlength - self.pathlength_proposal, 1)
            #     self.pathlengh = self.pathlength_proposal
            #     if random.random() < 0.5:
            #         self.termination = 0
                        
            # if self.pathlength + bifurcation.pathlength > neuron_pathlength + 10:
            #     self.pathlength = neuron_pathlength - bifurcation.pathlength
            #     self.termination = 2                

        self.euclidean = draw_value('contraction') * self.pathlength
        self.origin = bifurcation.bifurcating_compartment.end_coordinates
        self.orientation = vector_remote

        scalar_val = self.euclidean/math.sqrt(vector_remote[0]**2 +
                                              vector_remote[1]**2 +
                                              vector_remote[2]**2)
        self.remote_vector = [vector_remote[i]*scalar_val for i in range(0, 3)]
        self.fragmentation = int(draw_value('fragmentation'))
        self.contraction = draw_value('contraction')
        scalar_val = (self.euclidean/math.sqrt(vector_local[0]**2 +
                                               vector_local[1]**2 +
                                               vector_local[2]**2))/self.fragmentation
        self.local_vector = [vector_local[i] * scalar_val for i in range(0, 3)]

        self.end_coordinates = [self.origin[i] + self.remote_vector[i] for i in range(0, 3)]
        self.initial_coordinates = [self.origin[i] + self.local_vector[i] for i in range(0, 3)]


        self.compartment_list = [bifurcation.bifurcating_compartment] + make_compartments(self.origin, self.initial_coordinates, self.end_coordinates, self.fragmentation, self.contraction, self.start_diameter, self.end_diameter)
        self.compartment_list[-1].terminal_flag = 1
                
        if self.termination == 0 and flag == 1:
            self.daughter_bifurcation = Bifurcation(self)
#----------------------------------------------------------------------------
class Branch:
    def __init__(self, start_coordinates, start_vector, diameter):
        self.compartment_list = []
        self.start_coordinates = start_coordinates
        self.end_coordinates = start_coordinates
        self.diameter = diameter
        self.inertial_growth_vector = start_vector
        self.local_growth_vector = start_vector
        self.remote_growth_vector = start_vector
        self.x_axis = unit_vector([start_vector[1], -start_vector[0], 0])
        random_angle = draw_value('random_angle')
        self.random_vector = unit_vector([math.cos(random_angle), 0, math.sin(random_angle)])

        candidate = 0
        while (candidate < 0.8 or candidate > 1.0):
            candidate = random.gauss(0.9, 0.05)
        self.contraction = candidate

        self.parent_branch = []
        self.terminal_compartment = []
        self.local_GF = 0 
        self.terminal_flag = 1
        self.primary_flag = 2
        
        self.full_vector = [0,0,0]
        self.count = -1
        
        self.stem_length = random.gammavariate(stem_pathlength_coef[0], stem_pathlength_coef[2]) + stem_pathlength_coef[1]
        self.pathlength = 0
        self.pathlength_from_soma = euclidean_distance(self.start_coordinates, [0,0,0])
        self.order = -1
        self.canonical_end = [0,0,0]
        self.canonical_pathlength = 5
        self.grow_branch([])
        self.parent_bif_amp = 0

        global lamb
        global tau
        
        lamb_tmp2 = list(lamb)
        
        lambda_max = max(lamb_tmp2)
        max_N = len(lamb_tmp2)
        p_flag = 0
        scalar1 = 1
        
        delta = 1
        rate_time = np.linspace(0, int(len(lamb))-1, int(len(lamb)/delta))
        rate_res = delta
        
        L = int(10*tau/rate_res)
        refract = 1-np.exp(-rate_time[:L]/tau)
        # questionable = 1
        
        j = int(self.pathlength_from_soma)
        
        if (j+L) > len(lamb):
            lamb_tmp2[j:] *= refract[:len(lamb_tmp2)-j]
        else:
            lamb_tmp2[j:j+L] *= refract        
        
        while p_flag == 0:
            max_poisson2 = [0]*(max_N+1)   
            k = 0
            while k < max_N-1:
                u_k = np.random.uniform(0, 1)
                w_k = -math.log(u_k)/(lambda_max)
                # z_k = np.random.exponential(lambda_max)
                # k = k + 1 + int(z_k)
                k = k + 1 + int(w_k)
                if k >= max_N:
                    max_poisson2[-1] = 0
                else:
                    max_poisson2[k] = 1
            
            for j in range(0, len(max_poisson2)):
                if max_poisson2[j] == 1:
                    p = lamb_tmp2[j]/(scalar1*max(lamb_tmp2))

                    b = np.random.binomial(1, p)
                    if b == 0:
                        max_poisson2[j] = 0
                    else:
                        if (j+L) > len(lamb_tmp):
                            lamb_tmp2[j:] *= refract[:len(lamb_tmp2)-j]
                        else:
                            lamb_tmp2[j:j+L] *= refract
            max_poisson2[-1] = 1
            poisson_sum = sum(max_poisson2)

            
            if poisson_sum < 1 or poisson_sum > 7:
                p_flag = 1
            else:
                p_flag = 1
            # if sum(max_poisson2[0:75]) == 0:
            #     p_flag = 0
                

        self.point_process = max_poisson2
        # print poisson_sum



        new_coord = [0,0,0]
        end_coord = [0,0,0]
        self.count = 0
        for i in range(int(self.pathlength_from_soma),len(max_poisson2)):
            if max_poisson2[i] == 0:
                self.count = self.count + 1
            else:
                break
        global simulation_time
        
        if self.count + self.pathlength_from_soma > simulation_time:
            self.count = int(simulation_time - self.pathlength_from_soma)
 
        for i in range(0, 3):
            new_coord[i] = self.start_coordinates[i] + start_vector[i]
            end_coord[i] = self.start_coordinates[i] + self.count * start_vector[i]
        self.new_compartment_list = make_compartments(self.start_coordinates, new_coord, end_coord, self.count, 0.9, 0.5, 0.5)
        
        # self.compartment_list = self.new_compartment_list
        # self.end_coordinates =  self.new_compartment_list[-1].end_coordinates
        # self.point_process_id = 378
        
        self.canonical_end = [self.count*self.remote_growth_vector[i] + self.new_compartment_list[-1].end_coordinates[i] for i in range(0, 3)]
        self.canonical_pathlength = self.count
        
    def set_order(self, branch_order):
        self.order = branch_order
        
    def rescale_point_process(self, simulation_time):
        # Calculate scaling factor    
        list_of_indices = [i for i, x in enumerate(self.point_process) if x == 1]
        termination = list_of_indices[-1]
        scaling_factor = simulation_time/(len(self.point_process)-1)
        
        # Rescaled point process indices
        list_of_indices_scaled = [int(x*scaling_factor) for x in list_of_indices]
        new_point_process = [0]*len(self.point_process)
        for idx in range(0, len(list_of_indices_scaled)):
            new_point_process[list_of_indices_scaled[idx]] = 1

        self.point_process = new_point_process
        # print list_of_indices_scaled
        
    def grow_branch(self, branch_list):
        compartment_length = draw_value('compartment_length')
        self.pathlength = self.pathlength + compartment_length
        self.pathlength_from_soma = self.pathlength_from_soma + compartment_length
        if self.terminal_compartment != []:
            self.inertial_growth_vector = self.terminal_compartment.vector

        theta_list = []
        total_difference_vector = [0, 0, 0]
        biggest_angle = 0
        biggest_average = 0
        neighbor_gradient = 0 
        consumption_rate = 0
        distance_falloff = 0
      
        for branch2 in branch_list:
            if branch2 != self:
                r = euclidean_distance(branch2.end_coordinates, self.end_coordinates)
                difference_vector = [branch2.terminal_compartment.vector[i] - self.terminal_compartment.vector[i] for i in range(0, 3)]
                
                proj_xz = [difference_vector[0], 0, difference_vector[2]]
                # proj_xz[1] = 0
                if proj_xz[2] >= 0:
                    theta_list.append(math.acos(unit_vector(proj_xz)[0]))
                else:
                    theta_list.append(2*math.pi-math.acos(unit_vector(proj_xz)[0]))
                    
                # theta_list.append(math.acos(difference_vector[0]/r))
                                
                if r != 0:
                    total_difference_vector = [total_difference_vector[i] + difference_vector[i]/(r) for i in range(0, 3)]
                    # branch.local_GF = branch.local_GF - 300/((r**2) * (branch2.diameter**2))
                    neighbor_gradient = neighbor_gradient + consumption_rate*(1/(branch2.compartment_list[-1].diameter**2))*(distance_falloff/(r**2))
                    
                theta_list.sort()
                                
                biggest_angle = theta_list[0] + 2*math.pi - theta_list[-1]
                biggest_average = (theta_list[0] + 2*math.pi + theta_list[-1])/2
                for k in range(0, len(theta_list) - 1):
                    current_angle = theta_list[k+1] - theta_list[k]
                    if current_angle > biggest_angle:
                        biggest_angle = current_angle
                        biggest_average = (theta_list[k+1] + theta_list[k])/2

                
    
                
        self.neighbor_vector = (math.cos(biggest_average), 0, math.sin(biggest_average))
        if biggest_average == 0:
            self.neighbor_vector = [0, 1, 0]
        self.gradient_vector = [0,1,0]
        random_angle = draw_value('random_angle')

        self.random_vector = unit_vector([math.cos(random_angle), random.random(), math.sin(random_angle)])

        boundary = [0, 400, 0]
        # self.gradient_vector = unit_vector([boundary[i] - self.end_coordinates[i] for i in range(0, 3)])

        # print self.inertial_growth_vector, self.neighbor_vector, self.gradient_vector, self.random_vector

        # jitter
    
        theta_sd = 5
        phi_sd = 5
        # theta_sd = 15
        # phi_sd = 5
        growing_inertial_weight = 1
        growing_neighbor_weight = 0
        growing_gradient_weight = 0
        growing_random_weight = 0
        if self.pathlength > 50:
            growing_neighbor_weight = 0.5
        else:
            growing_neighbor_weight = 1
        growing_neighbor_weight = 0
        
        
        weight = 0.03
        if self.canonical_pathlength == 0:
            growing_inertial_weight = 1
        else:
            growing_inertial_weight = 1-weight*(float(len(self.compartment_list))/self.canonical_pathlength)
        growing_gradient_weight = 1-growing_inertial_weight
        
        if self.canonical_end == [0,0,0]:
            new_growth_vector = self.inertial_growth_vector
        else:
            new_growth_vector = unit_vector([self.canonical_end[i] - self.end_coordinates[i] for i in range(0,3)])
            
        new_growth_vector = self.inertial_growth_vector
        
        self.remote_growth_vector = new_growth_vector
        
        growing_inertial_weight = 1
        growing_neighbor_weight = 0
        growing_gradient_weight = 0
        growing_canonical_weight = 0
        
        # if self.count <= 0:
        #     growing_inertial_weight = 1
        # else:
        #     growing_inertial_weight = len(self.compartment_list)/self.count
        # growing_canonical_weight = 1 - (growing_inertial_weight)
            
        # print self.inertial_growth_vector, new_growth_vector, growing_inertial_weight, growing_gradient_weight      
        # self.total_vector = unit_vector([growing_inertial_weight*self.inertial_growth_vector[i] + growing_neighbor_weight*self.neighbor_vector[i] + growing_gradient_weight*new_growth_vector[i] + growing_random_weight*self.random_vector[i] for i in range(0, 3)])        
        self.total_vector = unit_vector([growing_inertial_weight*self.local_growth_vector[i] + growing_canonical_weight*self.remote_growth_vector[i] for i in range(0, 3)])        
        
        rand_theta = random.gauss(0, theta_sd)
        if self.compartment_list == []:
            last_phi = -90
            last_theta = 25
            self.counter = 1
        else:
            last_phi = self.compartment_list[-1].jitter_phi
            last_theta = self.compartment_list[-1].jitter_theta
            self.counter = self.compartment_list[-1].counter + 1
            
        if self.counter >= 10:
            self.counter = 0
            if last_phi >= 180:
                last_phi2 = last_phi - 180
            elif last_phi < 0:
                last_phi2 = last_phi + 180
            else:
                last_phi2 = last_phi - 180
                
            if last_theta >= 0:
                rand_theta = random.gauss(-25, theta_sd) - last_theta
            else:
                rand_theta = random.gauss(25, theta_sd) + last_theta
                        
            reverse_growth_vector = self.inertial_growth_vector
            reverse_growth_vector[0] = -reverse_growth_vector[0]
            reverse_growth_vector[2] = -reverse_growth_vector[2]
            # self.total_vector = unit_vector([growing_inertial_weight*reverse_growth_vector[i] + growing_canonical_weight*new_growth_vector[i] for i in range(0, 3)])
            rand_phi = -90
            # self.total_vector, self.x_axis = angle_noise([0,0,0], self.total_vector, self.total_vector, self.x_axis, rand_theta, rand_phi, 0)

        # self.total_vector = unit_vector([growing_inertial_weight*self.local_growth_vector[i] + growing_canonical_weight*self.remote_growth_vector[i] for i in range(0, 3)])        
        
        last_phi2 = -90
        # rand_phi = random.gauss(last_phi2, phi_sd)
        rand_phi = -90
        # rand_phi = random.uniform(-180, 180)
        rand_theta = random.gauss(0, theta_sd)
        rand_phi = random.gauss(-90, phi_sd)        
        
        # self.total_vector, self.x_axis = angle_noise([0,0,0], self.total_vector, self.total_vector, self.x_axis, rand_theta, rand_phi, 0)
        # print self.total_vector
        new_end_coordinates = [self.end_coordinates[i] + compartment_length*self.total_vector[i] for i in range(0, 3)]
        

        new_compartment = Compartment(self.end_coordinates, new_end_coordinates, self.diameter)
        # if self.start_coordinates == [0,0,0]:
                
        #     print new_end_coordinates
        new_compartment.jitter_theta = rand_theta
        new_compartment.jitter_phi = rand_phi
        if self.compartment_list != []:
            new_compartment.parent_compartment = self.compartment_list[-1]
        elif self.parent_branch != []:
            new_compartment.parent_compartment = self.parent_branch.compartment_list[-1]
        
        # if len(branch_list) > 50:
        #     self.terminal_flag == 0
        if self.terminal_flag == 1:    
            self.compartment_list.append(new_compartment)
        self.end_coordinates = new_end_coordinates
        self.terminal_compartment = self.compartment_list[-1]
        self.full_vector = [self.end_coordinates[i] - self.start_coordinates[i] for i in range(0,3)]

class Bifurcation:
    def __init__(self, bifurcating_compartment):
        self.bifurcating_compartment = bifurcating_compartment
        self.bifurcation_id = bifurcating_compartment.compartment_id
        self.diameter = bifurcating_compartment.diameter

        self.bifurcating_compartment.bifurcation_flag = 1

        self.euclidean_from_soma = euclidean_distance(bifurcating_compartment.end_coordinates, (0, 0, 0))

        self.parent_branch = []
        self.daughter_branch = []

        self.local_vector = (bifurcating_compartment.end_coordinates[i] - bifurcating_compartment.start_coordinates[i] for i in range(0, 3))
        self.remote_vector = (0, 0, 0)

        self.local_daughters = []
        self.remote_daughters = []

        self.vector1 = (0, 0, 0)
        self.vector2 = (0, 0, 0)
        self.vector1_local = (0, 0, 0)
        self.vector2_local = (0, 0, 0)

        self.projection= (0, 0, 1)

        self.daughter_ratio = 0
        self.pk = 0
        self.azimuth = -1
        self.elevation = -1
        self.bif_amp_local = -1
        self.bif_amp_remote = -1
        self.bif_tilt_local = -1
        self.bif_tilt_remote = -1
        self.bif_torque_local = -1
        self.bif_torque_remote = -1
        self.bif_twist_local = -1
        self.bif_twist_remote = -1
        self.bif_amp_vector = -1

        self.bif_midline = (0, 0, 0)
        self.bif_normal = (0, 0, 0)

        self.bif_midline_local = (0, 0, 0)
        self.bif_normal_local = (0, 0, 0)

    def set_parent_branch(self, parent_branch):
        self.parent_branch = parent_branch
        self.remote_vector = (self.bifurcating_compartment.end_coordinates[i] - parent_branch.start_compartment.end_coordinates[i] for i in range(0, 3))

    def set_daughter_branch(self, daughter_branch):
        self.daughter_branch.append(daughter_branch)
        self.remote_daughters.append(daughter_branch.end_compartment)
        self.local_daughters.append(daughter_branch.compartment_list[1])

        if len(self.daughter_branch) == 2:
            self.daughter_ratio = self.daughter_branch[1].start_diameter/self.daughter_branch[0].start_diameter
            if self.daughter_ratio < 1:
                self.daughter_ratio = 1/self.daughter_ratio
            self.pk = (self.daughter_branch[0].start_diameter**1.5 + self.daughter_branch[1].start_diameter**1.5)/(self.bifurcating_compartment.diameter**1.5)
class Bifurcation:
    def __init__(self, parent_branch):
        self.origin = parent_branch.end_coordinates
        self.orientation = parent_branch.orientation
        self.reference = parent_branch.orientation

        if parent_branch.parent_bifurcation != []:
            self.pathlength = parent_branch.pathlength + parent_branch.parent_bifurcation.pathlength
            
            self.azimuth = parent_branch.parent_bifurcation.azimuth
            self.elevation = parent_branch.parent_bifurcation.elevation
            self.orientation = parent_branch.parent_bifurcation.orientation
        else:
            self.pathlength = parent_branch.pathlength
            
            self.azimuth = 0
            self.elevation = 0
            self.orientation = 0
            
        self.bifurcating_compartment = parent_branch.compartment_list[-1]
        self.bifurcating_compartment.bifurcation_flag = 1

        self.parent_branch = parent_branch
        self.daughter_branches = []
        self.branch_order = parent_branch.branch_order + 1
        
        self.branch_order_1 = 0
        self.branch_order_2 = 0
        self.branch_order_3 = 0
        self.branch_order_4 = 0
        
        self.euclidean_from_soma = euclidean_distance(self.origin, (0, 0, 0))
        
        if self.euclidean_from_soma > 0 and self.euclidean_from_soma < 50:
            self.branch_order_1 = parent_branch.branch_order_1 + 1
        elif self.euclidean_from_soma > 50 and self.euclidean_from_soma < 150:
            self.branch_order_2 = parent_branch.branch_order_2 + 1
        elif self.euclidean_from_soma > 150 and self.euclidean_from_soma < 250:
            self.branch_order_3 = parent_branch.branch_order_3 + 1
        elif self.euclidean_from_soma > 250:
            self.branch_order_4 = parent_branch.branch_order_4 + 1


        self.amplitude = 0
        while self.amplitude <= 0:
            self.amplitude = draw_value('bif_amplitude_remote') * math.pi/180
        # if self.parent_branch.stem_flag == 1:
        #     self.amplitude = self.amplitude*2
            
        proposed_azimuth = draw_value('bif_azimuth_remote')
        # proposed_azimuth = 0
        self.azimuth = proposed_azimuth*math.pi/180 + self.azimuth
                
        proposed_elevation = 100
        proposed_orientation = 100

        while proposed_elevation * proposed_orientation > 200:        
            proposed_elevation = draw_value('bif_elevation_remote')
            if proposed_elevation < 0 and self.bifurcating_compartment.end_coordinates[2] < 0:
                proposed_elevation = draw_value('bif_elevation_remote')
            elif proposed_elevation > 0 and self.bifurcating_compartment.end_coordinates[2] > 0:
                proposed_elevation = draw_value('bif_elevation_remote')
            self.elevation = proposed_elevation*math.pi/180 + self.elevation
            
            proposed_orientation = draw_value('bif_orientation_remote')
            if proposed_orientation > 0 and self.bifurcating_compartment.end_coordinates[0] < 0:
                proposed_orientation = draw_value('bif_orientation_remote')
            elif proposed_orientation < 0 and self.bifurcating_compartment.end_coordinates[0] > 0:
                proposed_orientation = draw_value('bif_orientation_remote')
            self.orientation = proposed_orientation*math.pi/180 + self.orientation



        # print self.bifurcating_compartment.end_coordinates
        # print proposed_elevation, proposed_azimuth, proposed_orientation

        # v1_remote = np.array((math.sin(self.amplitude/2), math.cos(self.amplitude/2), 0))
        # v2_remote = np.array((math.sin(-self.amplitude/2), math.cos(-self.amplitude/2), 0))

        v1_remote = np.array((math.sin(0), math.cos(0), 0))
        v2_remote = np.array((math.sin(self.amplitude), math.cos(self.amplitude), 0))
        
        if self.branch_order <= 1:
            v1_remote = np.array((math.sin(self.amplitude/2), math.cos(self.amplitude/2), 0))
            v2_remote = np.array((math.sin(-self.amplitude/2), math.cos(-self.amplitude/2), 0))
            self.elevation = 0
            self.orientation = 0


        Rr_x = np.matrix( ((1, 0, 0), (0, math.cos(self.elevation), -math.sin(self.elevation)), (0, math.sin(self.elevation), math.cos(self.elevation))) )
        Rr_z = np.matrix( ((math.cos(self.orientation), -math.sin(self.orientation), 0), (math.sin(self.orientation), math.cos(self.orientation), 0), (0, 0, 1)) )
        Rr_y = np.matrix( ((math.cos(self.azimuth), 0, -math.sin(self.azimuth)), (0, 1, 0), (math.sin(self.azimuth), 0, math.cos(self.azimuth))) )
        
        R = Rr_y * Rr_z * Rr_x

        self.vector1_remote = np.ndarray.tolist(np.dot(v1_remote, R))[0]
        self.vector2_remote = np.ndarray.tolist(np.dot(v2_remote, R))[0]

        self.amplitude = 0
        while self.amplitude <= 0:
            self.amplitude = draw_value('bif_amplitude_local') * math.pi/180
        # if self.parent_branch.stem_flag == 1:
        #     self.amplitude = self.amplitude*2
            
        # self.azimuth = draw_value('bif_azimuth_local')*math.pi/180 + self.azimuth
        
        # self.elevation = draw_value('bif_elevation_local')*math.pi/180 + self.elevation
        # self.orientation = draw_value('bif_orientation_local')*math.pi/180 + self.orientation
        
        
        v1_local = np.array((math.sin(0), math.cos(0), 0))
        v2_local = np.array((math.sin(self.amplitude), math.cos(self.amplitude), 0))
        if self.branch_order <= 1:
            v1_local = np.array((math.sin(self.amplitude/2), math.cos(self.amplitude/2), 0))
            v2_local = np.array((math.sin(-self.amplitude/2), math.cos(-self.amplitude/2), 0))
            self.elevation = 0
            self.orientation = 0
            
        Rr_x = np.matrix( ((1, 0, 0), (0, math.cos(self.elevation), -math.sin(self.elevation)), (0, math.sin(self.elevation), math.cos(self.elevation))) )
        Rr_z = np.matrix( ((math.cos(self.orientation), -math.sin(self.orientation), 0), (math.sin(self.orientation), math.cos(self.orientation), 0), (0, 0, 1)) )
        Rr_y = np.matrix( ((math.cos(self.azimuth), 0, -math.sin(self.azimuth)), (0, 1, 0), (math.sin(self.azimuth), 0, math.cos(self.azimuth))) )
        
        R = Rr_y * Rr_z * Rr_x

        self.vector1_local = np.ndarray.tolist(np.dot(v1_local, R))[0]
        self.vector2_local = np.ndarray.tolist(np.dot(v2_local, R))[0]
        
        # self.vector1_local = self.vector1_remote
            # self.vector2_local = self.vector2_remote

            # self.midline = np.ndarray.tolist(np.dot(np.asarray((0,0,1)), R))[0]
            

            #     self.bif_amplitude_local = draw_value('bif_amplitude_local')
            #     self.bif_torque_local = draw_value('bif_torque_local')
            #     self.bif_tilt_local = draw_value('bif_tilt_local')

            # self.bif_amplitude_local = self.bif_amplitude_remote
            # self.bif_torque_local = self.bif_torque_remote
            # self.bif_tilt_local = self.bif_tilt_remote

            # v1_remote = np.array((math.sin(self.theta1_remote), math.sin(self.phi_remote), math.cos(self.theta1_remote)*math.cos(self.phi_remote)))
            # v2_remote = np.array((math.sin(self.theta2_remote), math.sin(self.phi_remote), math.cos(self.theta2_remote)*math.cos(self.phi_remote)))

            # v1_local = np.array((math.sin(self.theta1_local), math.sin(self.phi_local), math.cos(self.theta1_local)*math.cos(self.phi_local)))
            # v2_local = np.array((math.sin(self.theta2_local), math.sin(self.phi_local), math.cos(self.theta2_local)*math.cos(self.phi_local)))

            # R_v = np.matrix(((math.cos(self.psi), -math.sin(self.psi), 0), (math.sin(self.psi), math.cos(self.psi), 0), (0, 0, 1)))

            # v1_remote = np.dot(v1_remote, R_v)
            # v2_remote = np.dot(v2_remote, R_v)

            # v1_local = np.dot(v1_local, R_v)
            # v2_local = np.dot(v2_local, R_v)

            # u = self.orientation

            # norm = np.asarray(parent_branch.normal)
            # norm = [norm[i]/math.sqrt(magnitude_product(norm, norm)) for i in range(0, 3)]
            # n1 = (1, 0, 0)
            # n2 = (0, 1, 0)

            # w1 = [dot_product(u, n1)/magnitude_product(n1, n1), 0, 0]
            # w2 = [0, dot_product(u, n2)/magnitude_product(n2, n2), 0]

            # proj_yz = [self.orientation[i] - w1[i] for i in range(0, 3)]
            # proj_xz = [self.orientation[i] - w2[i] for i in range(0, 3)]

            # n = (0, 0, 1)

            # rot_x = math.acos(dot_product(proj_yz, n)/magnitude_product(proj_yz, n))
            # rot_y = math.acos(dot_product(proj_xz, n)/magnitude_product(proj_xz, n))
            # if proj_yz[1] < 0:
            #     rot_x = -rot_x
            # if proj_xz[0] < 0:
            #     rot_y = -rot_y

            # Rr_x = np.matrix(((1, 0, 0), (0, math.cos(-rot_x), -math.sin(-rot_x)), (0, math.sin(-rot_x), math.cos(-rot_x))))
            # Rr_y = np.matrix(((math.cos(-rot_y), 0, -math.sin(-rot_y)), (0, 1, 0), (math.sin(-rot_y), 0, math.cos(-rot_y))))

            # Rr_xy = Rr_x*Rr_y

            # new_norm = np.ndarray.tolist(np.dot(norm, Rr_xy))[0]
            # n_norm = (0, 0, 1)
            # rot_z = math.acos(dot_product(new_norm, n_norm)/magnitude_product(new_norm, n_norm))
            # if new_norm[2] > 0:
            #     rot_z = -rot_z

            # tropism = draw_value('tropism')

            # if self.branch_order >= 1:
            #     rot_x = tropism*rot_x
            #     rot_y = tropism*rot_y

            # R_x = np.matrix(((1, 0, 0), (0, math.cos(rot_x), -math.sin(rot_x)), (0, math.sin(rot_x), math.cos(rot_x))))
            # R_y = np.matrix(((math.cos(rot_y), 0, -math.sin(rot_y)), (0, 1, 0), (math.sin(rot_y), 0, math.cos(rot_y))))
            # R_z = np.matrix(((math.cos(rot_z), -math.sin(rot_z), 0), (math.sin(rot_z), math.cos(rot_z), 0), (0, 0, 1)))

            # R = R_x * R_y * R_z

            # self.vector1_remote = np.ndarray.tolist(np.dot(v1_remote, R))[0]
            # self.vector2_remote = np.ndarray.tolist(np.dot(v2_remote, R))[0]

            # self.vector1_local = np.ndarray.tolist(np.dot(v1_local, R))[0]
            # self.vector2_local = np.ndarray.tolist(np.dot(v2_local, R))[0]

            # if abs(self.vector1_remote[1]) < 0.000000000001:
            #     self.vector1_remote[1] = 0
            # if abs(self.vector2_remote[1]) < 0.000000000001:
            #     self.vector2_remote[1] = 0

        self.normal = cross_product(self.vector1_remote, self.vector2_remote)
        self.normal = [-self.normal[i] for i in range (0, 3)]

        self.parent_diameter = self.bifurcating_compartment.diameter
        self.parent_daughter_ratio = draw_value('parent_daughter_ratio')
        self.daughter_ratio = draw_value('daughter_ratio')

        if parent_branch.stem_flag == 1:
            self.parent_diameter = -1
            while self.parent_diameter <= 0:
                self.parent_diameter = draw_value('ibf_diameter')

        self.diameter1 = self.parent_daughter_ratio * self.parent_diameter
        self.diameter2 = self.diameter1/self.daughter_ratio


        branch1 = Branch(self, self.vector1_remote, self.vector1_local, self.normal, self.diameter1)
        branch2 = Branch(self, self.vector2_remote, self.vector2_local, self.normal, self.diameter2)

        self.daughter_branches.append(branch1)
        self.daughter_branches.append(branch2)

        parent_branch.daughter_bifurcation = self

class Morphology:
    def __init__(self, compartment_list, branch_list, bifurcation_list, terminal_list):
        self.compartment_list = compartment_list
        self.branch_list = branch_list
        self.bifurcation_list = bifurcation_list
        self.terminal_list = terminal_list

        self.num_stems = len([branch for branch in branch_list if branch.branch_order == 0])
        self.total_dendritic_length = 0
        self.total_surface_area = 0
        self.normalized_total_surface_area = 0

        # self.major_coordinates = [compartment.x for compartment in terminal_list]
        # self.minor_coordinates = [compartment.z for compartment in terminal_list]
       
        # self.tmp_major_axis = abs(max(self.major_coordinates) - min(self.major_coordinates))
        # self.tmp_minor_axis = abs(max(self.minor_coordinates) - min(self.minor_coordinates))
       
        original_coordinates = [np.array((compartment.x, compartment.z)) for compartment in terminal_list]
        
        self.major_axis = 0
        self.minor_axis = 0
        self.max_euclidean = 0
        
        self.GF_list = []
        self.GF_x_raw_list = []
        self.GF_y_raw_list = []
        self.GF_z_raw_list = []
        self.GF_neighbor_raw_list = []
        self.branch_flag_list = []

        
        
        max_y = 0
        
        for compartment in self.compartment_list:
            local_GF, GF_x_raw, GF_y_raw, GF_z_raw, GF_neighbor_raw = get_local_GF(compartment.end_coordinates, self.branch_list)
            self.GF_list.append(local_GF)
            self.GF_x_raw_list.append(GF_x_raw)
            self.GF_y_raw_list.append(GF_y_raw)
            self.GF_z_raw_list.append(GF_z_raw)
            self.GF_neighbor_raw_list.append(GF_neighbor_raw)
            if compartment.euclidean_from_soma > self.max_euclidean:
                self.max_euclidean = compartment.euclidean_from_soma
            if compartment.y > max_y:
                max_y = compartment.y



            self.branch_flag_list.append(compartment.bifurcation_flag)

        for angle in range(0,180):
            self.tmp_angle = angle *math.pi/180
            R = np.matrix(((math.cos(self.tmp_angle), -math.sin(self.tmp_angle)), (math.sin(self.tmp_angle), math.cos(self.tmp_angle))))
            
            rotated_coordinates = [np.ndarray.tolist(np.dot(entry, R))[0] for entry in original_coordinates]
            
            self.major_coordinates = [coordinates[0] for coordinates in rotated_coordinates]
            self.minor_coordinates = [coordinates[1] for coordinates in rotated_coordinates]
            
            self.tmp_major_axis = abs(max(self.major_coordinates) - min(self.major_coordinates))
            self.tmp_minor_axis = abs(max(self.minor_coordinates) - min(self.minor_coordinates))
            
            if self.tmp_major_axis > self.major_axis:
                self.major_axis = self.tmp_major_axis
                self.minor_axis = self.tmp_minor_axis
        
        # print(self.max_euclidean) # Testing purposes
        self.tdl_distribution = [0]*int(round(self.max_euclidean))        
        for branch in branch_list:
            self.total_dendritic_length = self.total_dendritic_length + branch.pathlength
        for compartment in compartment_list:
            tdl_extent = int(compartment.euclidean_from_soma)
            for tdl_idx in range(tdl_extent, len(self.tdl_distribution)):
                self.tdl_distribution[tdl_idx] = self.tdl_distribution[tdl_idx] + compartment.length

        self.length_scale_factor = self.total_dendritic_length/100    
        for compartment in compartment_list:
            self.total_surface_area = self.total_surface_area + math.pi*compartment.diameter*compartment.length
            self.normalized_total_surface_area = self.normalized_total_surface_area + math.pi*compartment.diameter*compartment.length/self.length_scale_factor
        for termination in terminal_list:
            self.total_surface_area = self.total_surface_area + math.pi*(0.5*termination.diameter)**2
        self.num_bifurcations = len(bifurcation_list)
        
        normed_bins = np.linspace(0, self.max_euclidean, len(sholl_bins))
        ynormed_bins = np.linspace(0, max_y, len(sholl_bins))
        self.sholl_counts = [0]*len(sholl_bins)
        self.ysholl_counts = [0]*len(sholl_bins)
        self.sholl_normed_counts = [0]*len(sholl_bins)
        self.ysholl_normed_counts = [0]*len(sholl_bins)
        # print self.total_surface_area
        for compartment in compartment_list:
            if compartment.parent_compartment != []:
                for i in range(0, len(sholl_bins)):
                    if compartment.euclidean_from_soma >= sholl_bins[i] and compartment.parent_compartment.euclidean_from_soma < sholl_bins[i]:
                        self.sholl_counts[i] = self.sholl_counts[i] + 1
                    if compartment.y >= sholl_bins[i] and compartment.parent_compartment.y < sholl_bins[i]:
                        self.ysholl_counts[i] = self.ysholl_counts[i] + 1
                    if compartment.euclidean_from_soma >= normed_bins[i] and compartment.parent_compartment.euclidean_from_soma < normed_bins[i]:
                        self.sholl_normed_counts[i] = self.sholl_normed_counts[i] + 1
                    if compartment.y >= ynormed_bins[i] and compartment.parent_compartment.y < ynormed_bins[i]:
                        self.ysholl_normed_counts[i] = self.ysholl_normed_counts[i] + 1
        self.sholl_counts = [x for x in self.sholl_normed_counts]         
        
        bif_angle_list = [x.bif_amp_remote for x in self.bifurcation_list]
        
        # bin_edges=np.histogram(bif_angle_list, bins=50)[1] #get the bin edges
        # plt.hist(bif_angle_list, alpha=0.7, normed=True, bins=bin_edges)
        # plt.plot(sholl_bins, self.sholl_counts)
        # plt.scatter(self.major_coordinates, self.minor_coordinates)
        # plt.show()     

class Data:
    def __init__(self):
        self.section_types = []
        self.diameters = []
        self.taper = []
        self.branch_order = []
        self.euclidean = []
        self.pathlength = []
        self.contraction = []
        self.fragmentation = []
        self.daughter_ratio = []
        self.pk = [] 
        self.local_bif_angle = []
        self.bif_branch_order = []

class Stem:
    def __init__(self, orientation):
        self.daughter_bifurcation = []
        self.parent_bifurcation = []

        self.orientation = orientation
        self.branch_order = 0
        
        self.branch_order_1 = 0
        self.branch_order_2 = 0
        self.branch_order_3 = 0
        self.branch_order_4 = 0

        flag = 0
        while flag == 0:
            self.pathlength = draw_value('stem_pathlength')

            if self.pathlength > 0:
                flag = 1
        self.euclidean = draw_value('contraction')*self.pathlength
        self.stem_flag = 1

        self.origin = [0, 0, 0]
        self.end_coordinates = [self.euclidean*orientation[0], self.euclidean*orientation[1], self.euclidean*orientation[2]]
        self.fragmentation = int(draw_value('fragmentation'))
        self.initial_coordinates = [0, self.euclidean/self.fragmentation, 0]
        self.contraction = draw_value('contraction')
        self.normal = (0, 0, 1)
        
        self.start_diameter = draw_value('stem_diameter')
        self.taper = draw_value('ibf_taper')
        self.end_diameter = self.start_diameter - self.taper*self.start_diameter

        self.compartment_list = make_compartments(self.origin, self.initial_coordinates, self.end_coordinates, self.fragmentation, self.contraction, self.start_diameter, self.end_diameter)
        self.compartment_list[0].stem_flag = 1

# # Combined Functions # # 
# Basic Operations
def euclidean_distance(u,v):
    ans = math.sqrt((u[0] - v[0])**2 + (u[1] - v[1])**2 + (u[2] - v[2])**2)
    return ans
def cross_product(u,v):
    s1 = u[1]*v[2] - u[2]*v[1]
    s2 = u[2]*v[0] - u[0]*v[2]
    s3 = u[0]*v[1] - u[1]*v[0]
    return (s1, s2, s3)
def dot_product(u,v):
    s = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    return s
def magnitude_product(u,v):
    s1 = math.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    s2 = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return s1*s2
def unit_vector(u):
    magnitude = euclidean_distance((u[0],u[1],u[2]), (0, 0, 0))
    
    if magnitude != 0:  
        return [u[0]/magnitude, u[1]/magnitude, u[2]/magnitude]
    else:
        return [0, 0, 0]
def rotation_angles(u, norm):
    # n1 = (1, 0, 0)
    # n2 = (0, 1, 0)

    # w1 = [dot_product(u,n1)/magnitude_product(n1,n1), 0, 0]
    # w2 = [0, dot_product(u,n2)/magnitude_product(n2,n2), 0]

    # proj_yz = [u[i] - w1[i] for i in range(0,3)]
    # proj_xz = [u[i] - w2[i] for i in range(0,3)]

    proj_yz = [0, u[1], u[2]]
    
    n = (0,1,0)
    
    if proj_yz == [0, 0, 0]:
        rot_x = 0
    else:
        rot_x = math.acos(dot_product(proj_yz, n)/magnitude_product(proj_yz, n))
        
    if proj_yz[2] > 0:
        rot_x = -rot_x
    
    Rr_x = np.matrix( ((1, 0, 0), (0, math.cos(-rot_x), -math.sin(-rot_x)), (0, math.sin(-rot_x), math.cos(-rot_x))) )
    
    new_u = np.ndarray.tolist(np.dot(u, Rr_x))[0]
    
    if new_u == [0,0,0]:
        rot_z = 0
    else:
        rot_z = math.acos(dot_product(new_u, n)/magnitude_product(new_u, n))
    if new_u[0] < 0:
        rot_z = -rot_z
    
    Rr_z = np.matrix( ((math.cos(-rot_z), -math.sin(-rot_z), 0), (math.sin(-rot_z), math.cos(-rot_z), 0), (0, 0, 1)) )
    
    Rr_xz = Rr_x*Rr_z
    
    
    if norm == (0, 0, 0):
        norm = (0, 0, 1)
    
    norm = [norm[i]/math.sqrt(magnitude_product(norm, norm)) for i in range(0,3)]
    new_norm = np.ndarray.tolist(np.dot(norm, Rr_xz))[0]
    n_norm = [0, 0, 1]
    
    rot_y = math.acos(dot_product(new_norm,n_norm)/magnitude_product(new_norm, n_norm))
    if new_norm[0] < 0:
        rot_y = -rot_y
    # Rr_y = np.matrix( ((math.cos(-rot_y), 0, -math.sin(-rot_y)), (0, 1, 0), (math.sin(-rot_y), 0, math.cos(-rot_y))) )

    # R_total = Rr_x*Rr_z*Rr_y
    
    # u_inv = np.ndarray.tolist(np.dot(u, R_total))[0]
    # norm_inv = np.ndarray.tolist(np.dot(norm, R_total))[0]
    # u_plot = zip((0,0,0), u)
    # norm_plot = zip((0,0,0), norm_inv)
    
    # print u_inv
    
    # ax = plt.axes(projection='3d')
    # ax.set_xlim([-50,50])
    # ax.set_ylim([-50,50])
    # ax.set_zlim([0,100])
    
    # print u_plot
    # plt.plot(u_plot[0], u_plot[1], u_plot[2])
    # plt.plot(norm_plot[0], norm_plot[1], norm_plot[2])
        
    # plt.show()

    # print rot_x
    return rot_x, rot_y, rot_z
def calculate_remote_angles(bifurcation_list):
    for bifurcation in bifurcation_list:       
        start_coordinates = bifurcation.bifurcating_compartment.end_coordinates
        end_coordinates1 = bifurcation.remote_daughters[0].end_coordinates
        end_coordinates2 = bifurcation.remote_daughters[1].end_coordinates

        bifurcation.vector1 = [end_coordinates1[i] - start_coordinates[i] for i in range(0, 3)]
        bifurcation.vector2 = [end_coordinates2[i] - start_coordinates[i] for i in range(0, 3)]  
        
        vector1 = bifurcation.vector1
        vector2 = bifurcation.vector2
        

        
        bif_norm = cross_product(vector1, vector2)
        bif_norm_mag = math.sqrt(sum(i**2 for i in bif_norm))

        if bif_norm_mag == 0:
            bifurcation.bif_normal = (0,1,0)
        else:
            bifurcation.bif_normal = (bif_norm[0]/bif_norm_mag, bif_norm[1]/bif_norm_mag, bif_norm[2]/bif_norm_mag)

        bifurcation.bif_amp_remote = math.acos(dot_product(vector1, vector2)/magnitude_product(vector1, vector2))*180/math.pi

        bifurcation.bif_midline = ((vector1[0]+vector2[0])/2, (vector1[1]+vector2[1])/2, (vector1[2]+vector2[2])/2)

    for bifurcation in bifurcation_list:
        if bifurcation.parent_branch.parent_bifurcation != []:
            rot_x1, rot_y1, rot_z1 = rotation_angles(bifurcation.bif_midline, bifurcation.bif_normal)
            rot_x2, rot_y2, rot_z2 = rotation_angles(bifurcation.parent_branch.vector, bifurcation.parent_branch.parent_bifurcation.bif_normal)
        
            rot_x_diff = rot_x2 - rot_x1
            rot_y_diff = rot_y2 - rot_y1
            rot_z_diff = rot_z2 - rot_z1

            
            if rot_y_diff > math.pi:
                rot_y_diff = rot_y_diff - 2*math.pi
            elif rot_y_diff < -math.pi:
                rot_y_diff = rot_y_diff + 2*math.pi


            rot_x = rot_x2
            rot_y = rot_y2
            rot_z = rot_z2

            R_x = np.matrix( ((1, 0, 0), (0, math.cos(rot_x), -math.sin(rot_x)), (0, math.sin(rot_x), math.cos(rot_x))) )
            R_y = np.matrix( ((math.cos(rot_y), 0, -math.sin(rot_y)), (0, 1, 0), (math.sin(rot_y), 0, math.cos(rot_y))) )
            R_z = np.matrix( ((math.cos(rot_z), -math.sin(rot_z), 0), (math.sin(rot_z), math.cos(rot_z), 0), (0, 0, 1)) )

            R = R_x * R_y * R_z

            v1_remote = np.asarray(bifurcation.vector1)
            v2_remote = np.asarray(bifurcation.vector2)

            vector1_remote = np.ndarray.tolist(np.dot(v1_remote, R))[0]
            vector2_remote = np.ndarray.tolist(np.dot(v2_remote, R))[0]

            bif_midline2 = ((vector1_remote[0]+vector2_remote[0])/2, (vector1_remote[1]+vector2_remote[1])/2, (vector1_remote[2]+vector2_remote[2])/2)
            bif_norm2 = cross_product(vector1_remote, vector2_remote)

            bifurcation.bif_torque_remote, bifurcation.bif_tilt_remote, bifurcation.bif_twist_remote = rotation_angles(bif_midline2, bif_norm2)

            # bifurcation.bif_tilt_remote = bifurcation.bif_tilt_remote * 180/math.pi
            # bifurcation.bif_torque_remote = bifurcation.bif_torque_remote * 180/math.pi
            # bifurcation.bif_twist_remote = bifurcation.bif_twist_remote * 180/math.pi
                
            bifurcation.bif_tilt_remote = rot_x_diff * 180/math.pi
            bifurcation.bif_torque_remote = rot_y_diff * 180/math.pi
            bifurcation.bif_twist_remote = rot_z_diff * 180/math.pi
            
    return
def calculate_local_angles(bifurcation_list):
    for bifurcation in bifurcation_list:
        start_coordinates = bifurcation.bifurcating_compartment.end_coordinates
        end_coordinates1 = bifurcation.local_daughters[0].end_coordinates
        end_coordinates2 = bifurcation.local_daughters[1].end_coordinates

        bifurcation.vector1_local = [end_coordinates1[i] - start_coordinates[i] for i in range(0, 3)]
        bifurcation.vector2_local = [end_coordinates2[i] - start_coordinates[i] for i in range(0, 3)]     

        vector1 = bifurcation.vector1_local
        vector2 = bifurcation.vector2_local

        bif_norm = cross_product(vector1, vector2)
        bif_norm_mag = math.sqrt(sum(i**2 for i in bif_norm))

        if bif_norm_mag == 0:
            bifurcation.bif_normal_local = (0, 1, 0)
        else:       
            bifurcation.bif_normal_local = (bif_norm[0]/bif_norm_mag, bif_norm[1]/bif_norm_mag, bif_norm[2]/bif_norm_mag)

        bifurcation.bif_amp_local = math.acos(dot_product(vector1, vector2)/magnitude_product(vector1, vector2))*180/math.pi

        bifurcation.bif_midline_local = ((vector1[0]+vector2[0])/2, (vector1[1]+vector2[1])/2, (vector1[2]+vector2[2])/2)
    
    for bifurcation in bifurcation_list:
        if bifurcation.parent_branch.parent_bifurcation != []:
            rot_x1, rot_y1, rot_z1 = rotation_angles(bifurcation.bif_midline, bifurcation.bif_normal)
            rot_x2, rot_y2, rot_z2 = rotation_angles(bifurcation.parent_branch.vector, bifurcation.parent_branch.parent_bifurcation.bif_normal)
        
            rot_x_diff = rot_x2 - rot_x1
            rot_y_diff = rot_y2 - rot_y1
            rot_z_diff = rot_z2 - rot_z1

            
            if rot_y_diff > math.pi:
                rot_y_diff = rot_y_diff - 2*math.pi
            elif rot_y_diff < -math.pi:
                rot_y_diff = rot_y_diff + 2*math.pi


            rot_x = rot_x2
            rot_y = rot_y2
            rot_z = rot_z2

            R_x = np.matrix( ((1, 0, 0), (0, math.cos(rot_x), -math.sin(rot_x)), (0, math.sin(rot_x), math.cos(rot_x))) )
            R_y = np.matrix( ((math.cos(rot_y), 0, -math.sin(rot_y)), (0, 1, 0), (math.sin(rot_y), 0, math.cos(rot_y))) )
            R_z = np.matrix( ((math.cos(rot_z), -math.sin(rot_z), 0), (math.sin(rot_z), math.cos(rot_z), 0), (0, 0, 1)) )

            R = R_x * R_z * R_y

            v1_local = np.asarray(bifurcation.vector1)
            v2_local = np.asarray(bifurcation.vector2)

            vector1_local = np.ndarray.tolist(np.dot(v1_local, R))[0]
            vector2_local = np.ndarray.tolist(np.dot(v2_local, R))[0]

            bif_midline_local_2 = ((vector1_local[0]+vector2_local[0])/2, (vector1_local[1]+vector2_local[1])/2, (vector1_local[2]+vector2_local[2])/2)
            bif_norm_local_2 = cross_product(vector1_local, vector2_local)

            bifurcation.bif_torque_local, bifurcation.bif_tilt_local, bifurcation.bif_twist_local = rotation_angles(bif_midline_local_2, bif_norm_local_2)

            # bifurcation.bif_tilt_local = bifurcation.bif_tilt_local * 180/math.pi
            # bifurcation.bif_torque_local = bifurcation.bif_torque_local * 180/math.pi
            # bifurcation.bif_twist_local = bifurcation.bif_twist_local * 180/math.pi

            bifurcation.bif_tilt_local = rot_x_diff * 180/math.pi
            bifurcation.bif_torque_local = rot_y_diff * 180/math.pi
            bifurcation.bif_twist_local = rot_z_diff * 180/math.pi 

    return
def angle_noise(origin, end_point, y_axis, x_axis, x_rot, y_rot, z_rot):
    
    origin_array = np.asarray(origin)
    end_array = np.asarray(end_point)
    x_axis_array = np.asarray(x_axis)
    
    vector_real = end_array - origin_array
    vector = np.asarray(y_axis)
    mag = np.sqrt(vector.dot(vector))
    
    vector_unit = vector/mag
    vector_list = np.ndarray.tolist(vector_unit)
    
    ux = vector_list[0]
    uy = vector_list[1]
    uz = vector_list[2]
    
    perp_vector_list = x_axis
    
    vx = perp_vector_list[0]
    vy = perp_vector_list[1]
    vz = perp_vector_list[2]
    
    perp_vector_list2 = cross_product(vector_list, perp_vector_list)
    
    wx = perp_vector_list2[0]
    wy = perp_vector_list2[1]
    wz = perp_vector_list2[2]


    rot_x = x_rot*math.pi/180
    rot_y = y_rot*math.pi/180
    rot_z = z_rot*math.pi/180


    # print rot_x, rot_y, rot_z

    # s = math.sin(rot_x)
    # c = math.cos(rot_x)
    # omc = 1-c

    R_x = np.matrix(( (math.cos(rot_x) + vx*vx*(1-math.cos(rot_x)), vx*vy*(1-math.cos(rot_x))-vz*math.sin(rot_x), vx*vz*(1-math.cos(rot_x))+vy*math.sin(rot_x) ),
                      (vy*vx*(1-math.cos(rot_x))+vz*math.sin(rot_x), math.cos(rot_x)+vy*vy*(1-math.cos(rot_x)), vy*vz*(1-math.cos(rot_x))-vx*math.sin(rot_x)),
                      (vz*vx*(1-math.cos(rot_x))-vy*math.sin(rot_x), vz*vy*(1-math.cos(rot_x))+vx*math.sin(rot_x), math.cos(rot_x)+vz*vz*(1-math.cos(rot_x))) ))
    
    
    R_y = np.matrix( ((math.cos(rot_y) + ux*ux*(1-math.cos(rot_y)), ux*uy*(1-math.cos(rot_y))-uz*math.sin(rot_y), ux*uz*(1-math.cos(rot_y))+uy*math.sin(rot_y)),
                      (uy*ux*(1-math.cos(rot_y))+uz*math.sin(rot_y), math.cos(rot_y)+uy*uy*(1-math.cos(rot_y)), uy*uz*(1-math.cos(rot_y))-ux*math.sin(rot_y)),
                      (uz*ux*(1-math.cos(rot_y))-uy*math.sin(rot_y), uz*uy*(1-math.cos(rot_y))+ux*math.sin(rot_y), math.cos(rot_y)+uz*uz*(1-math.cos(rot_y)))) )
   
    
    R_z = np.matrix( ((math.cos(rot_z) + wx*wx*(1-math.cos(rot_z)), wx*wy*(1-math.cos(rot_z))-wz*math.sin(rot_z), wx*wz*(1-math.cos(rot_z))+wy*math.sin(rot_z)),
                      (wy*wx*(1-math.cos(rot_z))+wz*math.sin(rot_z), math.cos(rot_z)+wy*wy*(1-math.cos(rot_z)), wy*wz*(1-math.cos(rot_z))-wx*math.sin(rot_z)),
                      (wz*wx*(1-math.cos(rot_z))-wy*math.sin(rot_z), wz*wy*(1-math.cos(rot_z))+wx*math.sin(rot_z), math.cos(rot_z)+wz*wz*(1-math.cos(rot_z)))) )

    # print R_x
    # print R_y
    # print R_z

    R = R_x * R_y * R_z
    
    v_array = np.asarray(vector_real)

    new_vector = np.ndarray.tolist(np.dot(v_array, R))[0]
    new_x_axis = np.ndarray.tolist(np.dot(x_axis_array, R))[0]
    
    new_end = [origin[i] + new_vector[i] for i in range(0, 3)]
    return new_end, new_x_axis
def calc_angle(vector1, vector2):
    if sum([(vector1[i] - vector2[i])**2 for i in range(0, 3)]) < 0.001:
        angle = 0
    elif sum([(-vector1[i] - vector2[i])**2 for i in range(0, 3)]) < 0.001:
        angle = 180
    else:
        angle = math.acos(dot_product(vector1, vector2)/magnitude_product(vector1, vector2))*180/math.pi
    return angle

# Visualization
def plot_branches(branch_list, graph_title = 'Blank Title'):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for branch in branch_list:
        if branch.end_compartment.section_type != 1:
            for compartment in branch.compartment_list:
                x_list.append(compartment.end_coordinates[0])
                y_list.append(compartment.end_coordinates[1])
                z_list.append(compartment.end_coordinates[2])
        plt.plot(x_list, z_list, y_list, color = 'black')
        x_list = []
        y_list = []
        z_list = []
    # Title Labeling 
    plt.title(graph_title)

    ax.set_xlim([-200,200])
    ax.set_ylim([-200,200])
    ax.set_zlim([0,400])
    plt.axis('off')
    plt.grid(b=None)
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.view_init(elev=0, azim=90)
    plt.show()
def plot_bifurcations(bifurcation_list):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for bifurcation in bifurcation_list:
        for branch in bifurcation.daughter_branch:
            for compartment in branch.compartment_list:
                x_list.append(compartment.x)
                y_list.append(compartment.y)
                z_list.append(compartment.z)
            plt.plot(x_list, y_list, z_list, color = color)
            x_list = []
            y_list = []
            z_list = []    

    # ax.view_init(elev=20, azim=-90)
    plt.show()    
def plot_bifurcation_feature(bifurcation_list):
    fig  = plt.figure()

    distribution = [bif.bif_tilt_remote for bif in bifurcation_list if bif.bif_tilt_remote != -1]
    # bins = np.linspace(-180, 180, 36)
    plt.hist(distribution, normed=True, bins=50)
    # print distribution
    # xt = plt.xticks()[0]
    # xmin, xmax = min(xt), max(xt)
    # lnspc = np.linspace(xmin, xmax, len(distribution))
    # be, ce = stats.expon.fit(distribution)
    # pdf_expon = stats.expon.pdf(lnspc, be, ce)
    # plt.plot(lnspc, pdf_expon, label = "Exp")
    # m, s = stats.norm.fit(distribution)
    # ag, bg, cg = stats.gamma.fit(distribution)
    # pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
    # pdf_norm = stats.norm.pdf(lnspc, m, s)
    # plt.plot(lnspc, pdf_gamma, label="Gamma")
    # plt.plot(lnspc, pdf_norm, label="Norm")
   
    # print("Exponential Coef", be, ce)
    # print("Normal Coef", m, s)
    # print("Gamma Coef", ag, bg, cg)
    
    plt.show()
    return sum(distribution)/len(distribution)
def plot_compartment_feature(compartment_list):
    fig  = plt.figure()

    distribution = [comp.diameter for comp in compartment_list]
    plt.hist(distribution, normed=True, bins=50)
    
    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(distribution))
    be, ce = stats.expon.fit(distribution)
    # pdf_expon = stats.expon.pdf(lnspc, be, ce)
    # plt.plot(lnspc, pdf_expon, label = "Exp")
    m, s = stats.norm.fit(distribution)
    ag, bg, cg = stats.gamma.fit(distribution)
    # pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
    pdf_norm = stats.norm.pdf(lnspc, m, s)
    # plt.plot(lnspc, pdf_gamma, label="Gamma")
    plt.plot(lnspc, pdf_norm, label="Norm")
    
    print("Exponential Coef", be, ce)
    print("Normal Coef", m, s)
    print("Gamma Coef", ag, bg, cg)

    # ks_check(distribution, 'stem_diameter')

    plt.show()
    return sum(distribution)/len(distribution)
def plot_pathlength_feature(branch_list):
    fig = plt.figure()

    distribution = [branch.start_compartment.diameter for branch in branch_list]

    plt.hist(distribution, normed=True, bins = 50)

    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(distribution))
    be, ce = stats.expon.fit(distribution)
    pdf_expon = stats.expon.pdf(lnspc, be, ce)
    plt.plot(lnspc, pdf_expon, label = "Exp")
    m, s = stats.norm.fit(distribution)
    ag, bg, cg = stats.gamma.fit(distribution)
    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
    pdf_norm = stats.norm.pdf(lnspc, m, s)
    plt.plot(lnspc, pdf_gamma, label="Gamma")
    plt.plot(lnspc, pdf_norm, label="Norm")
    
    print("Exponential Coef", be, ce)
    print("Normal Coef", m, s)
    print("Gamma Coef", ag, bg, cg)

    # ks_check(distribution, 'pathlength')

    plt.show()
def plot_pathlength_to_soma(branch_list):
    fig = plt.figure()

    distribution = [branch.start_compartment.diameter for branch in branch_list]

    plt.hist(distribution, normed=True, bins = 50)
    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(distribution))
    be, ce = stats.expon.fit(distribution)
    pdf_expon = stats.expon.pdf(lnspc, be, ce)
    plt.plot(lnspc, pdf_expon, label = "Exp")
    m, s = stats.norm.fit(distribution)
    ag, bg, cg = stats.gamma.fit(distribution)
    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
    pdf_norm = stats.norm.pdf(lnspc, m, s)
    plt.plot(lnspc, pdf_gamma, label="Gamma")
    plt.plot(lnspc, pdf_norm, label="Norm")
    
    print("Exponential Coef", be, ce)
    print("Normal Coef", m, s)
    print("Gamma Coef", ag, bg, cg)

    plt.show()    
def plot_euclidean_from_soma(branch_list):
    fig = plt.figure()

    distribution = [branch.euclidean_from_soma for branch in branch_list]

    plt.hist(distribution, normed=True, bins = 50)
    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(distribution))
    be, ce = stats.expon.fit(distribution)
    pdf_expon = stats.expon.pdf(lnspc, be, ce)
    plt.plot(lnspc, pdf_expon, label = "Exp")
    m, s = stats.norm.fit(distribution)
    ag, bg, cg = stats.gamma.fit(distribution)
    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
    pdf_norm = stats.norm.pdf(lnspc, m, s)
    plt.plot(lnspc, pdf_gamma, label="Gamma")
    plt.plot(lnspc, pdf_norm, label="Norm")
    

    # ks_check(distribution, 'euclidean_from_soma')
    #     print("Exponential Coef", be, ce)
    print("Normal Coef", m, s)
    print("Gamma Coef", ag, bg, cg)
    plt.show()   
def plot_branch_feature(branch_list):
    fig = plt.figure()

    distribution = [branch.taper2 for branch in branch_list]

    plt.hist(distribution, normed=True, bins = 50)

    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(distribution))
    be, ce = stats.expon.fit(distribution)
    pdf_expon = stats.expon.pdf(lnspc, be, ce)
    plt.plot(lnspc, pdf_expon, label = "Exp")
    m, s = stats.norm.fit(distribution)
    ag, bg, cg = stats.gamma.fit(distribution)
    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
    pdf_norm = stats.norm.pdf(lnspc, m, s)
    plt.plot(lnspc, pdf_gamma, label="Gamma")
    plt.plot(lnspc, pdf_norm, label="Norm")
    
    print("Exponential Coef", be, ce)
    print("Normal Coef", m, s)
    print("Gamma Coef", ag, bg, cg)

    plt.show()   

    return sum(distribution)/len(distribution)
def plot_morphology_feature(distribution):
    fig = plt.figure()

    plt.hist(distribution, normed=True, bins = 50)

    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    # lnspc = np.linspace(xmin, xmax, len(distribution))
    # be, ce = stats.expon.fit(distribution)
    # pdf_expon = stats.expon.pdf(lnspc, be, ce)
    # plt.plot(lnspc, pdf_expon, label = "Exp")
    # m, s = stats.norm.fit(distribution)
    # ag, bg, cg = stats.gamma.fit(distribution)
    # pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
    # pdf_norm = stats.norm.pdf(lnspc, m, s)
    # plt.plot(lnspc, pdf_gamma, label="Gamma")
    # plt.plot(lnspc, pdf_norm, label="Norm")
    
    # print("Exponential Coef", be, ce)
    # print("Normal Coef", m, s)
    # print("Gamma Coef", ag, bg, cg)

    plt.show()

    # ks_check(distribution, 'total_dendritic_length')
def plot_pairs(pair_list):
    fig2  = plt.figure()

    distribution1 = []
    distribution2 = []

    for pair in pair_list:
        item1 = pair[0].parent_branch.rot_z
        item2 = pair[0].parent_branch.end_coordinates[1]
        if item1 != -1 and item2 != -1:
            distribution1.append(item1)
            distribution2.append(item2)

    # H, xedges, yedges = np.histogram2d(distribution1, distribution2, bins=(6,50))
    # H = H.T    
    # plt.imshow(H, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='cool', vmin=-100, vmax=200)

    plt.scatter(distribution1, distribution2, marker='.')
    # plt.xlabel('Euclidean Distance from the Soma (um)')
    # plt.ylabel('Branch Pathlength (um)')
    plt.show()
    # plt.colorbar()
    # plt.xlim(0, 300)
    # plt.ylim(0,300)
    # plt.show()
def plot_euclidean_orders(orders_list):
    fig = plt.figure()

    num_orders = len(orders_list)
    for i in range(0,num_orders):
        plt.subplot(num_orders, 1, i + 1)
        distribution = [branch.pathlength for branch in orders_list[i]]

        plt.hist(distribution, normed=True, bins = 50)

        xt = plt.xticks()[0]
        xmin, xmax = min(xt), max(xt)
        lnspc = np.linspace(xmin, xmax, len(distribution))
        be, ce = stats.expon.fit(distribution)
        pdf_expon = stats.expon.pdf(lnspc, be, ce)
        plt.plot(lnspc, pdf_expon, label = "Exp")
        m, s = stats.norm.fit(distribution)
        ag, bg, cg = stats.gamma.fit(distribution)
        pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
        pdf_norm = stats.norm.pdf(lnspc, m, s)
        plt.plot(lnspc, pdf_gamma, label="Gamma")
        plt.plot(lnspc, pdf_norm, label="Norm")
    
        print("Exponential Coef", be, ce)
        print("Normal Coef", m, s)
        print("Gamma Coef", ag, bg, cg)

    plt.show()
    #    parameters = {
    #        'stem_diameter':(stem_diameter_constants[0], stem_diameter_constants[1], stem_diameter_constants[2], 'gamma'),
    #        'taper2': (taper_2_constants[0], taper_2_constants[1], taper_2_constants[2], 'gamma'),
    #        'soma_diameter': (soma_diameter_mean, soma_diameter_sd, 0,'gaussian'),
    #        'taper': (taper_2_mean, taper_2_sd, 0,'gaussian'),
    #        'daughter_ratio': (daughter_ratio_min, daughter_ratio_max, 0,'uniform'),
    #        'parent_daughter_ratio': (0.9,0,0,'gaussian'),
    #        'diam_threshold': (diam_threshold_mean, diam_threshold_sd, 0,'gaussian'),
    #        'fragmentation': (fragmentation_mean, fragmentation_sd, 0,'gaussian'),
    #        'pathlength': (branch_length_constants[0], branch_length_constants[1], branch_length_constants[2],'gamma'),
    #        'contraction': (0.9,0, 0,'gaussian'),
    #        'azimuth': (0.1,0, 0,'gaussian'),
    #        'elevation': (0.1,0, 0,'gaussian'),
    #        'bif_amplitude_remote': (bif_amplitude_remote_mean, bif_amplitude_remote_sd, 0,'gaussian'),
    #        'bif_tilt_remote': (bif_tilt_remote_mean, bif_tilt_remote_sd, 0,'gaussian'),
    #        'bif_torque_remote': (bif_torque_remote_mean, bif_torque_remote_sd, 0,'gaussian'),
    #        'bif_amplitude_local': (bif_amplitude_local_mean, bif_amplitude_local_sd, 0,'gaussian'),
    #        'bif_tilt_local': (bif_tilt_local_mean, bif_tilt_local_sd, 0,'gaussian'),
    #        'bif_torque_local': (bif_torque_local_mean, bif_torque_local_sd, 0,'gaussian'),
    #        'bif_twist': (bif_twist_min, bif_twist_max, 0,'uniform'),
    #        'tropism': (tropism_mean, tropism_sd, 0,'gaussian'),
    #        'num_bifurcations': (num_bifurcations_constants[0], num_bifurcations_constants[1], num_bifurcations_constants[2], 'gamma'),
    #        'total_dendritic_length': (total_dendritic_length_constants[0], total_dendritic_length_constants[1], total_dendritic_length_constants[2], 'gamma'),
    #        'branch_order': (branch_order_constants[0], branch_order_constants[1], 0, 'gaussian'),
    #        'pathlength_to_soma': (pathlength_to_soma_constants[0], pathlength_to_soma_constants[1], 0, 'gaussian'),
    #        'euclidean_from_soma': (euclidean_from_soma_constants[0], euclidean_from_soma_constants[1], 0, 'norm'),
    #        'diameter': (diameter_constants[0], diameter_constants[1], 0, 'gaussian'),

    #        'branch_order_0': (branch_order_0_constants[0], branch_order_0_constants[1], 0,'exponential'),
    #        'branch_order_1': (branch_order_1_constants[0], branch_order_1_constants[1], branch_order_1_constants[2],'gamma'),
    #        'branch_order_2': (branch_order_2_constants[0], branch_order_2_constants[1], 0,'exponential'),
    #        'branch_order_3': (branch_order_3_constants[0], branch_order_3_constants[1], branch_order_3_constants[2],'gamma'),
    #        'branch_order_4': (branch_order_4_constants[0], branch_order_4_constants[1], branch_order_4_constants[2],'gamma'),
    #        'branch_order_5': (branch_order_5_constants[0], branch_order_5_constants[1], branch_order_5_constants[2],'gamma'),
    #        'branch_order_6': (branch_order_6_constants[0], branch_order_6_constants[1], branch_order_6_constants[2],'gamma')    
    #    }

    #    distribution_type = parameters[request][-1]
    #    if distribution_type == 'gamma':
    #        fit_params = (parameters[request][0], parameters[request][1], parameters[request][2])
    #    else:
    #        fit_params = (parameters[request][0], parameters[request][1])

    #    (d_val, p_val) = stats.kstest(distribution, distribution_type, args=fit_params)
    #    print(d_val, p_val)
def plot_heatmap(bifurcation_list):
    x_distribution = [bif.bifurcating_compartment.x for bif in bifurcation_list]
    y_distribution = [bif.bifurcating_compartment.y for bif in bifurcation_list]
    angle_distribution = [bif.bif_twist_remote for bif in bifurcation_list]
    
    plt.scatter(x_distribution, y_distribution, c=angle_distribution, s = 50)
    plt.gray()
    
    plt.show()
def plot_avg_data(avg_data):

    plt.hist(avg_data, normed=True, bins = 50)
    plt.show()
def plot_terminal_feature(terminal_list):
    
    distribution = [termination.branch.pathlength for termination in terminal_list]

    plt.hist(distribution, normed=True, bins = 50)
    plt.show()

    return sum(distribution)/len(distribution)
def raster(event_times_list, **kwargs):
    """
    Creates a raster plot
    Parameters
    ----------
    event_times_list : iterable
                       a list of event time iterables
    color : string
            color of vlines
    Returns
    -------
    ax : an axis containing the raster plot
    """
    ax = plt.gca()
    for ith, trial in enumerate(event_times_list):
        plt.vlines(trial, ith + .5, ith + 1.5, **kwargs)
    plt.ylim(.5, len(event_times_list) + .5)
    return ax
def draw_value(request):
    parameters = {

        'compartment_length': 1,
        'GF_threshold': 1,
        'bif_amplitude_local': random.gammavariate(bif_amplitude_local_coef[0], bif_amplitude_local_coef[2]) + bif_amplitude_local_coef[1],
        'bif_amplitude_remote': random.gammavariate(bif_amplitude_remote_coef[0], bif_amplitude_remote_coef[2]) + bif_amplitude_remote_coef[1],
        # 'random_angle': random.random()*2*math.pi
        'random_angle': 0
    }

    return parameters[request]

# I/O
def read_file_swc(file_name):
    # Returns lists of all compartments, branches, bifurcations, and termination objects in the chosen swc file
    all_compartments = []
    compartment_dict = {}
    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()

    # Initialize last line check
    last_line = (-1, -1, -1)
    id_offset = 0
    # Compartment Initialization
    for i in range(0, len(lines)):
        parsed_line = lines[i].split(' ')
        line = [entry for entry in parsed_line if entry != '']  #ignores empty strings in the case of double spaces
        # Ignore last line if it is a duplicate of the second to last line
        if len(line) > 0:
            if (line[0] != '#') and line[0][0] != '#':
                if last_line != (line[2], line[3], line[4]):
                    # if line[1] == '3' or line[-1] == '-1':                    
                    new_compartment = Compartment(line)
                    if new_compartment.section_type == 3 or new_compartment.section_type == 4 or new_compartment.section_type == 1 :
                        if new_compartment.compartment_id != 1:
                            new_compartment.compartment_id = new_compartment.compartment_id - id_offset
                        if new_compartment.parent_id != 1:
                            new_compartment.parent_id = new_compartment.parent_id - id_offset

                        all_compartments.append(new_compartment)
                        compartment_dict[new_compartment.compartment_id] = new_compartment
                        if new_compartment.parent_id >= 1:
                            start_compartment = compartment_dict[new_compartment.parent_id]
                            new_compartment.parent_compartment = start_compartment
                            new_compartment.start_coordinates = start_compartment.end_coordinates
        
                            start_compartment.daughter_compartment.append(new_compartment)
                    else:
                        id_offset = id_offset + 1
                    last_line = (line[2], line[3], line[4])              
    for i in range(0, len(all_compartments)):
        all_compartments[i].set_length()

    # Bifurcation Initialization

    all_bifurcations = [Bifurcation(compartment) for compartment in all_compartments if len(compartment.daughter_compartment) > 1 and compartment.compartment_id != 1.0]
    
    # Breakpoint Initialization
    breakpoints = []
    for ii in range(len(all_compartments) - 2, -1, -1):
        if all_compartments[ii +1].parent_id != all_compartments[ii].compartment_id:
            if len(all_compartments[ii + 1].parent_compartment.daughter_compartment) == 1:
                breakpoints.append(all_compartments[ii + 1])

    # Termination Initialization
    all_terminations = [compartment for compartment in all_compartments if len(compartment.daughter_compartment) == 0]

    # Branch Initialization
    all_branches = divide_into_branches(all_compartments, all_bifurcations, breakpoints)

    calculate_remote_angles(all_bifurcations)
    calculate_local_angles(all_bifurcations)
    
    return all_compartments, all_branches, all_bifurcations, all_terminations
def read_file_pp(file_name):

    point_process_list = []

    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()
    
    for i in range(0, len(lines)):
        parsed_line = lines[i].split(' ') 
        
        point_process_list.append([0] + list(map(int,parsed_line)))
    
    return point_process_list
def read_file_conditionals(file_name):

    conditional_names1 = []
    conditional_names2 = []
    xedges = []
    yedges = []
    H = []
    
    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()
    
    for i in range(0, len(lines)):
        # for i in range(0, 1):

        parsed_line = lines[i].split(';')
        # print parsed_line
        
        tmp_name = parsed_line[0].split('_')
        conditional_names1.append(' '.join(tmp_name))
        
        del parsed_line[0]
        
        tmp_name2 = parsed_line[0].split('_')
        conditional_names2.append(' '.join(tmp_name2))
        triple_parsed = []
        sholl_append = []
        
        del parsed_line[0]
     
        xedges.append(list(map(float, parsed_line[0].split(' '))))
        del parsed_line[0]
        yedges.append(list(map(float, parsed_line[0].split(' '))))
        del parsed_line[0]
        
        # print parsed_line
        # distribution_types.append(parsed_line[0])
        # del parsed_line[0]
        # if parsed_line[0] != '':
        parsed_line = parsed_line[0].split(' ')
        # del parsed_line[0]
        # del parsed_line[-1]
        for entry in parsed_line:
            if entry[0] == '[':
                double_parsed = []
                double_parsed.append(entry[1:-1])
            elif entry[-1] == ']':
                double_parsed.append(entry[:-1])
                triple_parsed.append(double_parsed)
            else:
                double_parsed.append(entry[:-1])
        for j in range(0, len(triple_parsed)):
            sholl_append.append([float(x) for x in triple_parsed[j]])
        H.append(sholl_append)
            

        # else:
        #     distribution_entries.append([0])
   
    return conditional_names1, conditional_names2, xedges, yedges, H
def divide_into_branches(compartment_list, bifurcation_list, breakpoints):

    current_branch = []
    branches = []
    bifurcations = bifurcation_list
    bifurcation_compartments = [bif.bifurcating_compartment for bif in bifurcations]
    broken_branches = []
    
    last_parent_id = compartment_list[-1].compartment_id
    soma_added = 0
    broken_flag = 0
    # print breakpoints
    for i in range(0,len(compartment_list)):
        compartment = compartment_list[len(compartment_list)-i-1]
        broken_flag = 0
        if compartment != []:
            if len(broken_branches) != 0:
                # print broken_branches
                for jj in range(0, len(broken_branches)):
                    # print compartment.compartment_id, broken_branches[jj][-1].parent_id
                    if compartment.compartment_id == broken_branches[jj][-1].parent_id:
                        print('yes')
            for kk in range(0, len(breakpoints)):
                if compartment == breakpoints[kk]:
                    current_branch.append(compartment_list[int(last_parent_id)-1])
                    broken_branches.append(current_branch)
                    
                    current_branch = []
                    current_branch.append(compartment)
                    broken_flag = 1
                    break
            if broken_flag != 1:
                
                if compartment.compartment_id == last_parent_id and compartment not in bifurcation_compartments and last_parent_id != 1:
                    current_branch.append(compartment_list[int(last_parent_id)-1])                              
                else:                        
                    # print compartment.compartment_id
                    current_branch.append(compartment_list[int(last_parent_id)-1])
                    current_branch.reverse()
                    new_branch = Branch(current_branch)
                    branches.append(new_branch)
                    if compartment not in bifurcations:
                        if last_parent_id == 1:
                            if soma_added == 0:
                                soma_added = 1
                            branches[-1].branch_order = 0
    
                    # for i in range(0, len(bifurcation_list)):
    
                    if last_parent_id != 1:
                        bif_match = [bif for bif in bifurcation_list if int(bif.bifurcation_id) == int(last_parent_id)]
                        # print last_parent_id
                        # for bif in bifurcation_list:
                        #     print bif.bifurcation_id
                        winner = bif_match[0]   
                        # if bif_match = []:
                            
                        # new_branch.set_parent_bifurcation(bif_match[0])
                        
    
                    current_branch = []
                    current_branch.append(compartment)

            last_parent_id = compartment.parent_id
    for branch in branches:
                
        start_match = [bif for bif in bifurcation_list if int(bif.bifurcation_id) == int(branch.start_compartment.compartment_id)]
        if start_match != []:
            branch.set_parent_bifurcation(start_match[0])
        end_match = [bif for bif in bifurcation_list if int(bif.bifurcation_id) == int(branch.end_compartment.compartment_id)]
        if end_match != []:
            branch.set_daughter_bifurcation(end_match[0])
    for branch in branches:
        branch_order = 0
        branch_order_4 = 0
        branch_order_3 = 0
        branch_order_2 = 0
        branch_order_1 = 0
        
        parent_bif = branch.parent_bifurcation

        while parent_bif != []:
            if parent_bif.euclidean_from_soma > 250:
                branch_order_4 = branch_order_4 + 1
            elif parent_bif.euclidean_from_soma > 150:
                branch_order_3 = branch_order_3 + 1
            elif parent_bif.euclidean_from_soma > 50:
                branch_order_2 = branch_order_2 + 1
            else:
                branch_order_1 = branch_order_1 + 1
            parent_bif = parent_bif.parent_branch.parent_bifurcation
            branch_order = branch_order + 1

        branch.branch_order = branch_order
        branch.branch_order_4 = branch_order_4
        branch.branch_order_3 = branch_order_3
        branch.branch_order_2 = branch_order_2
        branch.branch_order_1 = branch_order_1
        
        current_compartment = branch.end_compartment
        
        while current_compartment.compartment_id != branch.start_compartment.compartment_id:

            current_compartment.branch_order = branch_order
            current_compartment = current_compartment.parent_compartment
    for branch in branches:
        if branch.branch_order == 0:
            branch.set_pathlength_to_soma()  
            
    return branches
def merge_swc(ap_input,ba_input, merge_flag):
    # Goal: take two swc file names for apical and basal datasets and merge into one swc file. This assumes both files have been pre-formatted.
    # Input: Arrays containing strings for the apical and basal file names. 
    # Action: Loops through files and outputs a file that is a merger of the two with the basal flipped on the vertical axis
    if merge_flag == 1 and len(ap_input) == len(ba_input):
        for i in range(len(ap_input)):
            output_title = '../Data/Morphologies/test_output/full_test_morphology_{}.swc'.format(i+1)
            # Apical
            with open(ap_input[i]) as ap_file:
                ap_data = ap_file.readlines()
            # Basal
            with open(ba_input[i]) as ba_file:
                ba_data = ba_file.readlines()
                if ba_data[0].split()[0] == '1':
                    del ba_data[0]
                for j in range(len(ba_data)):
                    split = ba_data[j].split()
                    split[3] = str(-float(split[3])) # Flipping y axis
                    if split[0] == '2':
                        split[0] = str(j + len(ap_data) + 1)
                        split[6] = '1'
                    else:
                        split[0] = str(int(split[0]) + len(ap_data) - 1)
                        split[6] = str(int(split[6]) + len(ap_data) - 1)

                    ba_data[j] = ' '.join(split) + '\n'
            # Merge
            full_data = ap_data + ba_data
            # Saving
            with open(output_title,'w+') as full_out:
                full_out.writelines(full_data)
        print('All Files Merged')
    elif len(ap_input) != len(ba_input):
        print('File list dimensions do not match.')
        return 0 
    else:
        print('No Merge Performed; Flag not 1')
        return 0
def list_aggregate(file_list):
    # Purpose
    # Output: the following lists:  compartment_list, branch_list, bifurcation_list, terminal_list, morphologies_list. The general, apical, and basal aggregates are outputted to separate lists
    # Initialize Empty arrays (not entirely necessary)
    general_aggregate = []
    apical_aggregate = []
    basal_aggregate = []
    # General
    all_morphologies_list = []
    all_compartments_list = []
    all_branches_list = []
    all_bifurcations_list = []
    all_terminations_list = []
    all_total_dendritic_length = []
    all_bifurcation_counts = []
    # Apical
    all_apical_arbor_list = []
    all_compartments_list_apical = []
    all_branches_list_apical = []
    all_bifurcations_list_apical = []
    all_terminations_list_apical = []
    all_total_dendritic_length_apical = []
    all_bifurcation_counts_apical = []
    # Basal
    all_basal_arbor_list = []
    all_compartments_list_basal = []
    all_branches_list_basal = []
    all_bifurcations_list_basal = []
    all_terminations_list_basal = []
    all_total_dendritic_length_basal = []
    all_bifurcation_counts_basal = []


    for i in range (0, len(file_list)-1):
        pair_table = []
        id_pairs = []
        ordered_pairs = []

        print(file_list[i])

        general_read = read_file(file_list[i])
        
        compartment_list = general_read[0]
        branch_list = general_read[1]
        bifurcation_list = general_read[2]
        termination_list = general_read[3]
    
        compartment_list, branch_list, bifurcation_list, termination_list = read_file(file_list[i])
        
        # # Apical and Basal Partitioning (Done earlier in read_file() )
        compartment_list_apical = [comp for comp in compartment_list if comp.section_type == 4 or comp.section_type == 1]
        branch_list_apical = [branch for branch in branch_list if branch.compartment_list[-1].section_type == 4 or branch.compartment_list[-1].section_type == 1]
        bifurcation_list_apical = [bif for bif in bifurcation_list if bif.bifurcating_compartment.section_type == 4 or bif.bifurcating_compartment.section_type == 1]
        termination_list_apical = [comp for comp in compartment_list_apical if len(comp.daughter_compartment) == 0]
        compartment_list_basal = [comp for comp in compartment_list if comp.section_type == 3 or comp.section_type == 1]
        branch_list_basal = [branch for branch in branch_list if branch.compartment_list[-1].section_type == 3 or branch.compartment_list[-1].section_type == 1]
        bifurcation_list_basal = [bif for bif in bifurcation_list if bif.bifurcating_compartment.section_type == 3 or bif.bifurcating_compartment.section_type == 1]
        termination_list_basal = [comp for comp in compartment_list_basal if len(comp.daughter_compartment) == 0]
        
        if i > 50:
            continue
        if compartment_list == [] and branch_list == []:
            continue
        # General
        current_morphology = Morphology(compartment_list, branch_list, bifurcation_list, termination_list)
        all_morphologies_list.append(current_morphology)
        all_compartments_list = all_compartments_list  + compartment_list
        all_branches_list = all_branches_list + branch_list
        all_bifurcations_list = all_bifurcations_list + bifurcation_list
        all_terminations_list = all_terminations_list + termination_list
        # Apical
        current_apical_arbor = Morphology(compartment_list_apical,branch_list_apical,bifurcation_list_apical,termination_list_apical)
        all_apical_arbor_list.append(current_apical_arbor)
        all_compartments_list_apical = all_compartments_list_apical + compartment_list_apical
        all_branches_list_apical = all_branches_list_apical + branch_list_apical
        all_bifurcations_list_apical = all_bifurcations_list_apical + bifurcation_list_apical
        all_terminations_list_apical = all_terminations_list_apical + termination_list_apical
        # Basal
        current_basal_arbor = Morphology(compartment_list_basal,branch_list_basal,bifurcation_list_basal,termination_list_basal)
        all_basal_arbor_list.append(current_basal_arbor)
        all_compartments_list_basal = all_compartments_list_basal + compartment_list_basal
        all_branches_list_basal = all_branches_list_basal + branch_list_basal
        all_bifurcations_list_basal = all_bifurcations_list_basal + bifurcation_list_basal
        all_terminations_list_basal = all_terminations_list_basal + termination_list_basal
        
        total_dendritic_length = 0

        for branch in branch_list:
            total_dendritic_length = total_dendritic_length + branch.pathlength
        all_total_dendritic_length.append(total_dendritic_length)

        all_bifurcation_counts.append(len(bifurcation_list))
    

    plot_branches(all_branches_list,'General')
    plot_branches(all_branches_list_apical,'Apical')
    plot_branches(all_branches_list_basal,'Basal')
    
    # Packaging for function return
    general_aggregate = [all_compartments_list,all_branches_list,all_bifurcations_list,all_terminations_list,all_morphologies_list]
    apical_aggregate = [all_compartments_list_apical,all_branches_list_apical,all_bifurcations_list_apical,all_terminations_list_apical,all_apical_arbor_list]
    basal_aggregate = [all_compartments_list_basal,all_branches_list_basal,all_bifurcations_list_basal,all_terminations_list_basal,all_basal_arbor_list]
    return general_aggregate, apical_aggregate, basal_aggregate

# Analysis
def save_distributions(cumulative_list):
    # Input: compartment_list, birfuction_list, terminal_list, morphologies_list, all inside another list.
    # Operation: Processes the data contained in each list to return a metric. This metric is then appended to the output list
    # Output: distribution_list (name), distribution_type_list (type), entries_list (data)
    
    # Extracting data from function input
    compartment_list = cumulative_list[0]
    branch_list = cumulative_list[1]
    bifurcation_list = cumulative_list[2] 
    terminal_list = cumulative_list[3]
    morphologies_list = cumulative_list[4]

    # Declaring empty lists to be filled
    distribution_list = []
    distribution_type_list = []
    entries_list = []

    stem_list = [stem for stem in branch_list if stem.branch_order == 0]
    ibf_list = [ibf.parent_branch for ibf in bifurcation_list if ibf.parent_branch.branch_order != 0]
    terminal_branches = [branch for branch in branch_list if branch.daughter_bifurcation == []]
    
    distribution_list.append('Stem_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([stem.pathlength for stem in stem_list]) 
    
    # distribution_list.append('All_pathlength_GCL_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 50])
   
    # distribution_list.append('All_pathlength_IML_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 150 and ibf.euclidean_from_soma > 50])
   
    # distribution_list.append('All_pathlength_MML_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([ibf.pathlength for ibf in branch_list if (ibf.euclidean_from_soma > 150 and ibf.euclidean_from_soma <= 250)])

    # distribution_list.append('All_pathlength_OML_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma > 250])    
   
    # distribution_list.append('Terminal_Pathlength_GCL_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma <= 50])

    # distribution_list.append('Terminal_Pathlength_IML_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 50 and termination.branch.euclidean_from_soma <= 150])
   
    # distribution_list.append('Terminal_pathlength_MML_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([termination.branch.pathlength for termination in terminal_list if (termination.branch.euclidean_from_soma > 150 and termination.branch.euclidean_from_soma <= 250)])

    # distribution_list.append('Terminal_pathlength_OML_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 250])  

    distribution_list.append('Interbifurcation_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([branch.pathlength for branch in branch_list])
    
    # # Apical Data
    # distribution_list_apical.append('Interbifurcation_Pathlength_Apical(um)')
    # distribution_type_list_apical.append('basic')
    # entries_list_apical.append([branch.pathlength for branch in branch_list_apical])
    # # Basal Data
    # distribution_list_basal.append('Interbifurcation_Pathlength_Basal(um)')
    # distribution_type_list_basal.append('basic')
    # entries_list_basal.append([branch.pathlength for branch in branch_list_basal])
    
    pathlength_distribution = [branch.pathlength for branch in branch_list]
    branch_distance_distribution = [branch.pathlength_to_soma - branch.pathlength for branch in branch_list]
    max_distance = max(branch_distance_distribution)
    
    branch_distance_distribution = []
    
    pathlength_distribution = []
    rotx_distribution = []
    rotz_distribution = []
    previous_branch_length = []
    previous_bif_amp = []
    previous_rot_x = []
    previous_rot_z = []
    previous_x = []
    previous_z = []
    bif_amp_distribution2 = []
    dummy = []
    for branch in branch_list:
        if branch.parent_bifurcation != []:
            branch_distance_distribution.append(branch.pathlength_to_soma - branch.pathlength)
            pathlength_distribution.append(branch.pathlength)
            rotx_distribution.append(branch.rot_x)
            rotz_distribution.append(branch.rot_z)
            previous_branch_length.append(branch.parent_bifurcation.parent_branch.pathlength)
            # previous_bif_amp.append(branch.parent_bifurcation.bif_amp_remote)
            previous_rot_x.append(branch.parent_bifurcation.parent_branch.rot_x)
            previous_rot_z.append(branch.parent_bifurcation.parent_branch.rot_z)
            previous_x.append(branch.start_coordinates[0])
            previous_z.append(branch.start_coordinates[2])
            bif_amp_distribution2.append(branch.parent_bifurcation.bif_amp_remote)
            dummy.append(random.random())
            if branch.parent_bifurcation.parent_branch.parent_bifurcation != []:
                previous_bif_amp.append(branch.parent_bifurcation.parent_branch.parent_bifurcation.bif_amp_remote)
            else:
                previous_bif_amp.append(0)
        else:
            # previous_branch_length.append(0)
            # previous_bif_amp.append(0)
            # previous_rot_x.append(0)
            # previous_rot_z.append(0)
            # previous_x.append(0)
            # previous_z.append(0)
            # bif_amp_distribution2.append(0)
            dummy.append(random.random())
            
    normalized_distance = [x*1000/max_distance for x in branch_distance_distribution]

    # distribution_list.append('Stem_Diameter_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([stem.start_diameter for stem in stem_list])

    # distribution_list.append('IBF_Diameter_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([ibf.start_diameter for ibf in ibf_list])

    # distribution_list.append('Terminal_Diameter_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([termination.branch.start_diameter for termination in terminal_list])

    # distribution_list.append('Compartment_Diameter_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([compartment.diameter for compartment in compartment_list if compartment.section_type == 3])
    
    # distribution_list.append('Post-bifurcation_Diameter_1_(um)')
     # distribution_type_list.append('basic')
     # entries_list.append([max(bifurcation.daughter_branch[0].start_diameter, bifurcation.daughter_branch[1].start_diameter) for bifurcation in bifurcation_list])
    
    # distribution_list.append('IBF_Taper')
     # distribution_type_list.append('basic')
     # entries_list.append([ibf.taper2 for  ibf in ibf_list])

    # distribution_list.append('Terminal_Taper')
     # distribution_type_list.append('basic')
     # entries_list.append([termination.branch.taper2 for  termination in terminal_list])
   
    # distribution_list.append('Daughter_Ratio')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.daughter_ratio for bifurcation in bifurcation_list])   

    # distribution_list.append('Parent_Daughter_Ratio')
     # distribution_type_list.append('basic')
     # entries_list.append([max(bifurcation.daughter_branch[0].start_diameter, bifurcation.daughter_branch[1].start_diameter)/bifurcation.diameter for bifurcation in bifurcation_list])   


    distribution_list.append('Bifurcation_Amplitude_Remote(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_amp_remote for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    distribution_list.append('Bifurcation_Amplitude_Vector(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_amp_vector for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    bif_amp_distribution = [bifurcation.bif_amp_vector for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1]
    bif_distance_distribution = [bifurcation.daughter_branch[0].pathlength_to_soma for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1]

    # distribution_list.append('Elevation_Remote(deg)')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.bif_tilt_remote for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    # distribution_list.append('Azimuth_Remote(deg)')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.bif_torque_remote for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
       
    # distribution_list.append('Orientation_Remote(deg)')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.bif_twist_remote for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])

    # distribution_list.append('Bifurcation_Amplitude_Local(deg)')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.bif_amp_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    # distribution_list.append('Elevation_Local(deg)')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.bif_tilt_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    # distribution_list.append('Azimuth_Local(deg)')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.bif_torque_local for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
       
    # distribution_list.append('Orientation_Local(deg)')
     # distribution_type_list.append('basic')
     # entries_list.append([bifurcation.bif_twist_local for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])

    distribution_list.append('rot_x')
    distribution_type_list.append('basic')
    entries_list.append([branch.rot_x for branch in branch_list])

    # rotx_distribution = [branch.rot_x for branch in branch_list]

    distribution_list.append('rot_z')
    distribution_type_list.append('basic')
    entries_list.append([branch.rot_z for branch in branch_list])
    
    # rotz_distribution = [branch.rot_z for branch in branch_list]
    
    distribution_list.append('x_coordinate')
    distribution_type_list.append('conditional')
    entries_list.append([bifurcation.bifurcating_compartment.x for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])
    
    distribution_list.append('z_coordinate')
    distribution_type_list.append('conditional')
    entries_list.append([bifurcation.bifurcating_compartment.z for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])    

    # distribution_list.append('Terminal_Euclidean_to_soma_(um)')
     # distribution_type_list.append('conditional')
     # entries_list.append([termination.branch.euclidean_from_soma for termination in terminal_list])

    # distribution_list.append('IBF_Euclidean_to_soma_(um)')
     # distribution_type_list.append('conditional')
     # entries_list.append([ibf.euclidean_from_soma for ibf in ibf_list])

    # distribution_list.append('Pre-bifurcation_Diameter_(um)')
     # distribution_type_list.append('conditional')
     # entries_list.append([bifurcation.bifurcating_compartment.diameter for bifurcation in bifurcation_list])
   
    # distribution_list.append('Pre-bifurcation_Pathlength_(um)')
     # distribution_type_list.append('conditional')
     # entries_list.append([bifurcation.parent_branch.pathlength for bifurcation in bifurcation_list])         

    distribution_list.append('Total_Dendritic_Length_(um)')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.total_dendritic_length for morphology in morphologies_list])
    
    distribution_list.append('Branch_Order')
    distribution_type_list.append('emergent')
    entries_list.append([branch.branch_order for branch in branch_list])

    
    distribution_list.append('Number_of_Bifurcations')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.num_bifurcations for morphology in morphologies_list])

    # distribution_list.append('Total_Surface_Area_(um^2)')
     # distribution_type_list.append('emergent')
     # entries_list.append([morphology.total_surface_area for morphology in morphologies_list])
    
     # distribution_list.append('Normalized_Total_Surface_Area_(um^2)')
     # distribution_type_list.append('emergent')
     # entries_list.append([morphology.normalized_total_surface_area for morphology in morphologies_list])    
    
    distribution_list.append('Major_axis_(um)')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.major_axis for morphology in morphologies_list])
    
    distribution_list.append('Minor_axis_(um)')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.minor_axis for morphology in morphologies_list])
    
    distribution_list.append('sholl_bins')
    distribution_type_list.append('emergent')
    entries_list.append(sholl_bins)

    distribution_list.append('sholl_counts')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.sholl_normed_counts for morphology in morphologies_list])

    # distribution_list.append('max_euclidean')
     # distribution_type_list.append('emergent')
     # entries_list.append([morphology.max_euclidean for morphology in morphologies_list])

    distribution_list.append('tdl_distribution')
    distribution_type_list.append('emergent')
    
    max_max_euclidean = 0
    for morphology in morphologies_list:
        if morphology.max_euclidean > max_max_euclidean:
            max_max_euclidean = morphology.max_euclidean
    
    max_max_euclidean = int(math.ceil(max_max_euclidean))
    for morphology in morphologies_list:
        while len(morphology.tdl_distribution) < max_max_euclidean:
            morphology.tdl_distribution.append(morphology.tdl_distribution[-1])
    entries_list.append([morphology.tdl_distribution for morphology in morphologies_list])

    tdl_bins = np.linspace(0, max_max_euclidean-1, max_max_euclidean)

    tmp_list = [morphology.tdl_distribution for morphology in morphologies_list]
    plt.plot(tdl_bins, [sum(x)/len(morphologies_list) for x in zip(*tmp_list)], label='Generated TDL')    
    plt.title('TDL for Generated')
    plt.legend()

    # print(len(morphologies_list)) # Testing
    
    # CROSS CORRELATIONS
    rotx_mean = stats.circmean(rotx_distribution)
    rotz_mean = stats.circmean(rotz_distribution)
    
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for idx in range(1, len(rotx_distribution)):
        a = rotx_distribution[idx]
        b = rotz_distribution[idx-1]
        a_mean = rotx_mean
        b_mean = rotz_mean
        sum1 = sum1 + math.sin(a - a_mean)*math.sin(b - b_mean)
        sum2 = sum2 + math.sin(a - a_mean)**2
        sum3 = sum3 + math.sin(b - b_mean)**2
    
    ra = sum1/math.sqrt(sum2*sum3)    
    
    a = np.corrcoef(branch_distance_distribution, pathlength_distribution)
    b = np.corrcoef(bif_distance_distribution, bif_amp_distribution)
    c = np.corrcoef(branch_distance_distribution, rotx_distribution)
    d = np.corrcoef(rotx_distribution, rotz_distribution)
    
    
    previous_bif_amp = [x for x in previous_bif_amp]
    previous_rot_x = [x*180/(math.pi) for x in previous_rot_x]
    previous_rot_z = [x*180/(math.pi) for x in previous_rot_z]
    bif_amp_distribution2 = [x for x in bif_amp_distribution2]

    rotx_distribution = [x*180/(math.pi) for x in rotx_distribution]
    rotz_distribution = [x*180/(math.pi) for x in rotz_distribution]
    
    data = [branch_distance_distribution, previous_bif_amp, previous_rot_x, previous_rot_z, previous_x, previous_z, pathlength_distribution, bif_amp_distribution2, rotx_distribution, rotz_distribution]
    # data = [normalized_distance, previous_bif_amp, previous_rot_x, previous_rot_z, previous_x, previous_z, pathlength_distribution, bif_amp_distribution2, rotx_distribution, rotz_distribution]
    
    data_length = len(data)
    data_length2 = len(data[0])
    
    
    names = ['branch_distance_distribution', 'previous_bif_amp', 'previous_rot_x', 'previous_rot_z', 'previous_x', 'previous_z', 'pathlength_distribution', 'bif_amp_distribution2', 'rotx_distribution', 'rotz_distribution']
    # datfra = pandas.DataFrame(np.array(data).reshape(data_length2, data_length), columns = names)
    # correlations = datfra.corr()
    
    correlation_list = [x[:] for x in [[0]*data_length] * data_length]
    for i in range(0, data_length):
        for j in range(0, data_length):
            correlation_list[i][j] = np.corrcoef(data[i], data[j])[0][1]
    
    correlation_matrix = np.array(correlation_list)
    
    correlations = correlation_matrix
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(correlations, vmin=0, vmax=1, cmap='viridis')
    fig.colorbar(cax)
    ticks = np.arange(0,len(names),1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(names, rotation=90)
    ax.set_yticklabels(names)
    ax.set_xticklabels([x+1 for x in ticks])
    ax.set_yticklabels([x+1 for x in ticks])  
    
    
    fig = plt.figure()
    plt.scatter(branch_distance_distribution, pathlength_distribution)
    plt.title('pathlength')

    
    
    for j in range(6,7):
        for i in range(0,1):
            
            dist1 = data[i]
            dist2 = data[j]
            
            distribution1 = []
            distribution2 = []
            
            for idx in range(0, len(dist1)):
                if dist1[idx] < 550:
                    distribution1.append(dist1[idx])
                    distribution2.append(dist2[idx])
            
            
            distribution1 = [x for x in distribution1 if x < 550]
            
            fig = plt.figure()
            H, xedges, yedges = np.histogram2d(distribution1, distribution2, bins=(8,25), normed=True)
            # print(xedges, yedges) # Possibly for testing purposes
            H = H.T
            H_max = np.amax(H)
            column_max = np.amax(H, axis=0)

            weights = []
            for x in column_max:
                if x == 0:
                    y = 0
                else:
                    y = H_max/x
                weights.append(y)

            # if column_max == 0: # Bodged the value here. Another value might be needed
            #     weights = [0 for x in column_max]
            # else:
            #     weights = [H_max/x for x in column_max]
            for k in range(len(weights)):
                H[:,k] *= weights[k]
                
            x_edge_min = 0
            x_edge_max = 600
            y_edge_min =0
            y_edge_max = 500
            vmin_val = 0
            vmax_val = 0.00005
            plt.imshow(H, interpolation='nearest', aspect='auto', origin='low', extent=[x_edge_min, x_edge_max, y_edge_min, y_edge_max], cmap='viridis', vmin=vmin_val, vmax=H_max)


            plt.colorbar()
            # plt.title(names[i] + ' vs ' + names[j])

    return distribution_list, distribution_type_list, entries_list
def save_distributions(cumulative_list, conditionals_file):
    # Input Format: [compartments_list, branches_list, bifurcations_list, terminations_list, morphologies_list, distribution_list, distribution_type_list, entries_list]
    # cumulative_list
    compartment_list = cumulative_list[0]
    branch_list = cumulative_list[1]
    bifurcations_list = cumulative_list[2]
    terminations_list = cumulative_list[3]
    morphologies_list = cumulative_list[4]
    distribution_list = cumulative_list[5]
    distribution_type_list = cumulative_list[6]
    entries_list = cumulative_list[7]

    stem_list = [stem for stem in branch_list if stem.branch_order == 0]
    ibf_list = [ibf.parent_branch for ibf in bifurcation_list if ibf.parent_branch.branch_order != 0]
    terminal_branches = [branch for branch in branch_list if branch.daughter_bifurcation == []]
    
    distribution_list.append('Stem_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([stem.pathlength for stem in stem_list]) 
    
    # distribution_list.append('All_pathlength_GCL_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 50])
    
    # distribution_list.append('All_pathlength_IML_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 150 and ibf.euclidean_from_soma > 50])
    
    # distribution_list.append('All_pathlength_MML_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([ibf.pathlength for ibf in branch_list if (ibf.euclidean_from_soma > 150 and ibf.euclidean_from_soma <= 250)])

    # distribution_list.append('All_pathlength_OML_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma > 250])    
    
    # distribution_list.append('Terminal_Pathlength_GCL_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma <= 50])

    # distribution_list.append('Terminal_Pathlength_IML_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 50 and termination.branch.euclidean_from_soma <= 150])
    
    # distribution_list.append('Terminal_pathlength_MML_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([termination.branch.pathlength for termination in terminal_list if (termination.branch.euclidean_from_soma > 150 and termination.branch.euclidean_from_soma <= 250)])

    # distribution_list.append('Terminal_pathlength_OML_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 250])  

    distribution_list.append('Interbifurcation_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([branch.pathlength for branch in branch_list])
    
    pathlength_distribution = [branch.pathlength for branch in branch_list]
    branch_distance_distribution = [branch.pathlength_to_soma - branch.pathlength for branch in branch_list]
    max_distance = max(branch_distance_distribution)
    
    branch_distance_distribution = []
    
    pathlength_distribution = []
    rotx_distribution = []
    rotz_distribution = []
    previous_branch_length = []
    previous_bif_amp = []
    previous_rot_x = []
    previous_rot_z = []
    previous_x = []
    previous_z = []
    bif_amp_distribution2 = []
    rotz_local_distribution = []
    rotx_local_distribution = []
    dummy = []
    for branch in branch_list:
        if branch.parent_bifurcation != []:
            branch_distance_distribution.append(branch.pathlength_to_soma - branch.pathlength)
            pathlength_distribution.append(branch.pathlength)
            rotx_distribution.append(branch.rot_x)
            rotz_distribution.append(branch.rot_z)
            rotx_local_distribution.append(branch.rot_x_local)
            rotz_local_distribution.append(branch.rot_z_local)
            previous_branch_length.append(branch.parent_bifurcation.parent_branch.pathlength)
            # previous_bif_amp.append(branch.parent_bifurcation.bif_amp_remote)
            previous_rot_x.append(branch.parent_bifurcation.parent_branch.rot_x)
            previous_rot_z.append(branch.parent_bifurcation.parent_branch.rot_z)
            previous_x.append(branch.start_coordinates[0])
            previous_z.append(branch.start_coordinates[2])
            bif_amp_distribution2.append(branch.parent_bifurcation.bif_amp_remote)
            dummy.append(random.random())
            if branch.parent_bifurcation.parent_branch.parent_bifurcation != []:
                previous_bif_amp.append(branch.parent_bifurcation.parent_branch.parent_bifurcation.bif_amp_remote)
            else:
                previous_bif_amp.append(0)
        else:
            # previous_branch_length.append(0)
            # previous_bif_amp.append(0)
            # previous_rot_x.append(0)
            # previous_rot_z.append(0)
            # previous_x.append(0)
            # previous_z.append(0)
            # bif_amp_distribution2.append(0)
            dummy.append(random.random())
            
    normalized_distance = [x*1000/max_distance for x in branch_distance_distribution]
    
    

    distribution_list.append('Stem_Diameter_(um)')
    distribution_type_list.append('basic')
    entries_list.append([stem.start_diameter for stem in stem_list])

    distribution_list.append('IBF_Diameter_(um)')
    distribution_type_list.append('basic')
    entries_list.append([ibf.start_diameter for ibf in ibf_list])

    # distribution_list.append('TerminaXl_Diameter_(um)')
    # distribution_type_list.append('basic')
    # entries_list.append([termination.branch.start_diameter for termination in terminal_list])

    distribution_list.append('Compartment_Diameter_(um)')
    distribution_type_list.append('basic')
    entries_list.append([compartment.diameter for compartment in compartment_list if compartment.section_type == 3])
        
    distribution_list.append('Post-bifurcation_Diameter_1_(um)')
    distribution_type_list.append('basic')
    entries_list.append([max(bifurcation.daughter_branch[0].start_diameter, bifurcation.daughter_branch[1].start_diameter) for bifurcation in bifurcation_list])
        
    distribution_list.append('IBF_Taper')
    distribution_type_list.append('basic')
    entries_list.append([ibf.taper2 for  ibf in ibf_list])

    # distribution_list.append('Terminal_Taper')
    # distribution_type_list.append('basic')
    # entries_list.append([termination.branch.taper2 for  termination in terminal_list])
    
    distribution_list.append('Daughter_Ratio')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.daughter_ratio for bifurcation in bifurcation_list])   

    distribution_list.append('Parent_Daughter_Ratio')
    distribution_type_list.append('basic')
    entries_list.append([max(bifurcation.daughter_branch[0].start_diameter, bifurcation.daughter_branch[1].start_diameter)/bifurcation.diameter for bifurcation in bifurcation_list])   


    distribution_list.append('Bifurcation_Amplitude_Remote(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_amp_remote for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    distribution_list.append('Bifurcation_Amplitude_Vector(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_amp_vector for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    bif_amp_distribution = [bifurcation.bif_amp_vector for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1]
    bif_distance_distribution = [bifurcation.daughter_branch[0].pathlength_to_soma for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1]

    # distribution_list.append('Elevation_Remote(deg)')
    # distribution_type_list.append('basic')
    # entries_list.append([bifurcation.bif_tilt_remote for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    # distribution_list.append('Azimuth_Remote(deg)')
    # distribution_type_list.append('basic')
    # entries_list.append([bifurcation.bif_torque_remote for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
        
    # distribution_list.append('Orientation_Remote(deg)')
    # distribution_type_list.append('basic')
    # entries_list.append([bifurcation.bif_twist_remote for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])

    # distribution_list.append('Bifurcation_Amplitude_Local(deg)')
    # distribution_type_list.append('basic')
    # entries_list.append([bifurcation.bif_amp_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    # distribution_list.append('Elevation_Local(deg)')
    # distribution_type_list.append('basic')
    # entries_list.append([bifurcation.bif_tilt_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    # distribution_list.append('Azimuth_Local(deg)')
    # distribution_type_list.append('basic')
    # entries_list.append([bifurcation.bif_torque_local for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
        
    # distribution_list.append('Orientation_Local(deg)')
    # distribution_type_list.append('basic')
    # entries_list.append([bifurcation.bif_twist_local for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])

    distribution_list.append('rot_x')
    distribution_type_list.append('basic')
    entries_list.append([branch.rot_x for branch in branch_list])

    # rotx_distribution = [branch.rot_x for branch in branch_list]

    distribution_list.append('rot_z')
    distribution_type_list.append('basic')
    entries_list.append([branch.rot_z for branch in branch_list])
    
    # rotz_distribution = [branch.rot_z for branch in branch_list]
    
    distribution_list.append('x_coordinate')
    distribution_type_list.append('conditional')
    entries_list.append([bifurcation.bifurcating_compartment.x for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])
    
    distribution_list.append('z_coordinate')
    distribution_type_list.append('conditional')
    entries_list.append([bifurcation.bifurcating_compartment.z for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])    

    # distribution_list.append('Terminal_Euclidean_to_soma_(um)')
    # distribution_type_list.append('conditional')
    # entries_list.append([termination.branch.euclidean_from_soma for termination in terminal_list])

    # distribution_list.append('IBF_Euclidean_to_soma_(um)')
    # distribution_type_list.append('conditional')
    # entries_list.append([ibf.euclidean_from_soma for ibf in ibf_list])

    # distribution_list.append('Pre-bifurcation_Diameter_(um)')
    # distribution_type_list.append('conditional')
    # entries_list.append([bifurcation.bifurcating_compartment.diameter for bifurcation in bifurcation_list])
    
    # distribution_list.append('Pre-bifurcation_Pathlength_(um)')
    # distribution_type_list.append('conditional')
    # entries_list.append([bifurcation.parent_branch.pathlength for bifurcation in bifurcation_list])         

    distribution_list.append('Total_Dendritic_Length_(um)')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.total_dendritic_length for morphology in morphologies_list])
    
    distribution_list.append('Branch_Order')
    distribution_type_list.append('emergent')
    entries_list.append([branch.branch_order for branch in branch_list])

    # distribution_list.append('Branch_Order_1_(um)')
    # distribution_type_list.append('emergent')
    # entries_list.append([branch.branch_order_1 for branch in branch_list if branch.branch_order_1 != 0])
    
    # distribution_list.append('Branch_Order_2_(um)')
    # distribution_type_list.append('emergent')
    # entries_list.append([branch.branch_order_2 for branch in branch_list if branch.branch_order_2 != 0])
    
    # distribution_list.append('Branch_Order_3_(um)')
    # distribution_type_list.append('emergent')
    # entries_list.append([branch.branch_order_3 for branch in branch_list if branch.branch_order_3 != 0])
    
    # distribution_list.append('Branch_Order_4_(um)')
    # distribution_type_list.append('emergent')
    # entries_list.append([branch.branch_order_4 for branch in branch_list if branch.branch_order_4 != 0])
    
    
    distribution_list.append('Number_of_Bifurcations')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.num_bifurcations for morphology in morphologies_list])

    # distribution_list.append('Total_Surface_Area_(um^2)')
    # distribution_type_list.append('emergent')
    # entries_list.append([morphology.total_surface_area for morphology in morphologies_list])
        
    # distribution_list.append('Normalized_Total_Surface_Area_(um^2)')
    # distribution_type_list.append('emergent')
    # entries_list.append([morphology.normalized_total_surface_area for morphology in morphologies_list])    
    
    
    distribution_list.append('Major_axis_(um)')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.major_axis for morphology in morphologies_list])
    
    distribution_list.append('Minor_axis_(um)')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.minor_axis for morphology in morphologies_list])
    
    distribution_list.append('sholl_bins')
    distribution_type_list.append('emergent')
    entries_list.append(sholl_bins)

    distribution_list.append('sholl_counts')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.sholl_counts for morphology in morphologies_list])
    
    # CROSS CORRELATIONS
    rotx_mean = stats.circmean(rotx_distribution)
    rotz_mean = stats.circmean(rotz_distribution)
    
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for idx in range(1, len(rotx_distribution)):
        a = rotx_distribution[idx]
        b = rotz_distribution[idx-1]
        a_mean = rotx_mean
        b_mean = rotz_mean
        sum1 = sum1 + math.sin(a - a_mean)*math.sin(b - b_mean)
        sum2 = sum2 + math.sin(a - a_mean)**2
        sum3 = sum3 + math.sin(b - b_mean)**2
    
    # ra = sum1/math.sqrt(sum2*sum3)    
    
    a = np.corrcoef(branch_distance_distribution, pathlength_distribution)
    b = np.corrcoef(bif_distance_distribution, bif_amp_distribution)
    c = np.corrcoef(branch_distance_distribution, rotx_distribution)
    d = np.corrcoef(rotx_distribution, rotz_distribution)
    
    previous_bif_amp = [x for x in previous_bif_amp]
    previous_rot_x = [x*180/(math.pi) for x in previous_rot_x]
    previous_rot_z = [x*180/(math.pi) for x in previous_rot_z]
    bif_amp_distribution2 = [x for x in bif_amp_distribution2]

    rotx_distribution = [x*180/(math.pi) for x in rotx_distribution]
    rotz_distribution = [x*180/(math.pi) for x in rotz_distribution]
    
    rotx_local_distribution = [x*180/(math.pi) for x in rotx_local_distribution]
    rotz_local_distribution = [x*180/(math.pi) for x in rotz_local_distribution]    
    
    data = [branch_distance_distribution, previous_bif_amp, previous_rot_x, previous_rot_z, previous_x, previous_z, pathlength_distribution, bif_amp_distribution2, rotx_distribution, rotz_distribution, rotx_local_distribution, rotz_local_distribution]
    # data = [normalized_distance, previous_bif_amp, previous_rot_x, previous_rot_z, previous_x, previous_z, pathlength_distribution, bif_amp_distribution2, rotx_distribution, rotz_distribution]
    
    data_length = len(data)
    data_length2 = len(data[0])
    
    # for i in range(5):
    #     for j in range(data_length2):
    #         data[i][j] = random.random()
    
    names = ['branch_distance_distribution', 'previous_bif_amp', 'previous_rot_x', 'previous_rot_z', 'previous_x', 'previous_z', 'pathlength_distribution', 'bif_amp_distribution2', 'rotx_distribution', 'rotz_distribution', 'rotx_local', 'rotz_local']
    # datfra = pandas.DataFrame(np.array(data).reshape(data_length2, data_length), columns = names)
    # correlations = datfra.corr()
    
    correlation_list = [x[:] for x in [[0]*data_length] * data_length]
    for i in range(0, data_length):
        for j in range(0, data_length):
            # print i, j
            correlation_list[i][j] = np.corrcoef(data[i], data[j])[0][1]
    
    correlation_matrix = np.array(correlation_list)
    
    correlations = correlation_matrix
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(correlations, vmin=0, vmax=1, cmap='viridis')
    fig.colorbar(cax)
    ticks = np.arange(0,len(names),1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(names, rotation=90)
    ax.set_yticklabels(names)
    # ax.set_xticklabels([x+1 for x in ticks])
    # ax.set_yticklabels([x+1 for x in ticks])  
    
    
    fig = plt.figure()
    plt.scatter(branch_distance_distribution, pathlength_distribution)
    plt.title('pathlength')
    
    H_list = []
    edge_list = []
    conditional_names1 = []
    conditional_names2 = []
    
    for j in range(6, 12):
        for i in range(0, 6):
            distribution1 = data[i]
            distribution2 = data[j]
            fig = plt.figure()
            H, xedges, yedges = np.histogram2d(distribution1, distribution2, bins=(8,25), normed=True)
            H = H.T
            H_max = np.amax(H)
            column_max = np.amax(H, axis=0)

            weights = []
            for x in column_max:
                if x == 0:
                    y = 0
                else:
                    y = H_max/x
                weights.append(y)

            # weights = [H_max/x for x in column_max]
            for k in range(len(weights)):
                H[:,k] *= weights[k]
            plt.imshow(H, interpolation='nearest', aspect='auto', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='viridis', vmin=0, vmax=H_max)
            plt.colorbar()
            plt.title(names[i] + ' vs ' + names[j])
            plt.show()
            H_list.append(np.ndarray.tolist(H))
            edge_list.append([np.ndarray.tolist(xedges), np.ndarray.tolist(yedges)])
            conditional_names1.append(names[i])
            conditional_names2.append(names[j])
               
    # # Saving to File Temporarily Suspended Due to Testing Purposes           
    # with open(conditionals_file, 'w') as f:
    #     for i in range(0, len(H_list)):
    #         xedges = ' '.join(str(x) for x in edge_list[i][0])
    #         yedges = ' '.join(str(x) for x in edge_list[i][1])
    #         Hentry = ' '.join(str(x) for x in H_list[i])
    #         entry = xedges + ';' + yedges + ';' + Hentry
    #         strs = conditional_names1[i] + ';' + conditional_names2[i] + ';' + entry
    #         f.write(strs + '\n')
    #     f.closed
    # #

    # fig = plt.figure()
    # plt.scatter(previous_rot_x, rotx_distribution)
    # plt.title('rotx')
    # fig = plt.figure()
    # plt.scatter(previous_rot_z, rotz_distribution)
    # plt.title('rotz')

    return distribution_list, distribution_type_list, entries_list

# Generation
def conditional_probability(index, conditions):4
    conditionals_file = '../Data/Metrics/claiborne_conditionals.txt'   

    conditional_names1, conditional_names2, xedges, yedges, H = read_conditional_file(conditionals_file) 
        
    heatmap = H[index]
    for idx in range(len(xedges[index])):
        x_bin = idx - 1
        if conditions < xedges[index][idx]:
            break
    if x_bin == -1:
        x_bin = 0
    conditional_distribution = []
    for ii in range(len(H[index])):
        conditional_distribution.append(H[index][ii][x_bin])
    cond_sum = sum(conditional_distribution)
    conditional_distribution2 = [0]*len(conditional_distribution)
    cum_sum = 0
    for ii in range(len(conditional_distribution)):
        cum_sum = cum_sum + conditional_distribution[ii]
        conditional_distribution2[ii] = cum_sum/cond_sum
    random_draw = random.random()
    
    for idx in range(len(conditional_distribution2)):
        cond_bin = idx
        if random_draw < conditional_distribution2[idx]:
            break
        
    return random.uniform(yedges[index][cond_bin], yedges[index][cond_bin + 1])
def make_compartments(origin, initial_coordinates, end_coordinates, fragmentation, contraction, start_diameter, end_diameter):
    compartment_list = []
    new_origin = origin
    # local_coordinates = [new_origin[i] + initial_coordinates]
    compartment_list.append(Compartment(origin, initial_coordinates, start_diameter))
    new_origin = initial_coordinates

    new_fragmentation = fragmentation
    new_end_coordinates = end_coordinates
    vector = (end_coordinates[0]-new_origin[0], end_coordinates[1]-new_origin[1], end_coordinates[2]-new_origin[2])
    # euclidean = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

    vector = (initial_coordinates[0] - origin[0], initial_coordinates[1] - origin[1], initial_coordinates[2] - origin[2])
    rho = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

    if abs(vector[2]-rho) < 0.0000000001:
        parent_phi = math.acos(1)
    else:
        parent_phi = math.acos(vector[2]/rho)

    if abs(vector[1] - rho*math.sin(parent_phi)) < 0.0000000001:
        parent_theta = math.asin(1)
    else:
        if abs(abs(vector[1]) - abs(rho*math.sin(parent_phi))) < 0.0000000001:
            parent_theta = math.asin(-1)
        else:
            if math.sin(parent_phi) == 0:
                parent_theta = 0
            else:
                parent_theta = math.asin(vector[1]/(rho*math.sin(parent_phi)))
            if vector[0] < 0:
                parent_theta = math.pi - parent_theta

    for i in range(1, fragmentation):
        new_vector = [(end_coordinates[j] - new_origin[j])/new_fragmentation for j in range(0, 3)]
        new_end_coordinates = [new_origin[j]+new_vector[j] for j in range(0, 3)]
        new_diameter = start_diameter - (start_diameter-end_diameter)*i/fragmentation

        vector = (end_coordinates[0] - new_origin[0], end_coordinates[1] - new_origin[1], end_coordinates[2] - new_origin[2])
        rho = math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

        if abs(vector[2]-rho) < 0.0000000001:
            phi = math.acos(1)
        else:
            phi = math.acos(vector[2]/rho)

        if abs(vector[1] - rho*math.sin(phi)) < 0.0000000001:
            theta = math.asin(1)
        else:
            if abs(abs(vector[1]) - abs(rho*math.sin(phi))) < 0.0000000001:
                theta = math.asin(-1)
            else:
                if math.sin(phi) == 0:
                    theta = 0
                else:
                    theta = math.asin(vector[1]/(rho*math.sin(phi)))
                if vector[0] < 0:
                    theta = math.pi - theta

        jitter_deviation = 0.5
        phi_adjust = random.gauss(0, jitter_deviation)
        theta_adjust = random.gauss(0, jitter_deviation)

        target_weight = 1-weight*(float(new_fragmentation)/fragmentation)

        # if new_fragmentation < fragmentation/2:
        #     target_weight = 1
        # else:
        #     target_weight = 0

        new_theta = theta*target_weight + parent_theta*(1-target_weight) + (1-target_weight)*theta_adjust
        new_phi = phi*target_weight + parent_phi*(1-target_weight) + (1-target_weight)*phi_adjust

        parent_theta = new_theta
        parent_phi = new_phi

        dx = (rho/new_fragmentation) * math.cos(new_theta) * math.sin(new_phi)
        dy = (rho/new_fragmentation) * math.sin(new_theta) * math.sin(new_phi)
        dz = (rho/new_fragmentation) * math.cos(new_phi)

        new_end_coordinates = [new_origin[0]+dx, new_origin[1]+dy, new_origin[2]+dz]

        compartment_list.append(Compartment(new_origin, new_end_coordinates, new_diameter))

        new_origin = new_end_coordinates
        new_fragmentation = new_fragmentation-1

    compartment_list.append(Compartment(new_origin, end_coordinates, end_diameter))
    compartment_list[-1].parent_compartment = compartment_list[-2]

    return compartment_list
def simulate_growth(branch_list, arg1, arg2, arg3, arg4, arg5, t, simulation_time,conditionals_file):
    global max_prob
    theta_sd = 0
    phi_sd = 0
    for branch in branch_list:

        # if len(branch_list) > 100:
        #     branch.terminal_flag = 0
        if branch.terminal_flag == 1:
            total_difference_vector = [0, 0, 0]
            boundary = [0, max_N, 0]
            gradient_vector = [0, 1, 0]
            
            rand_theta = random.gauss(0, theta_sd)
            rand_phi = random.gauss(0, phi_sd)
            random_vector = angle_noise([0, 0, 0], branch.inertial_growth_vector, branch.inertial_growth_vector, branch.x_axis, rand_theta, rand_phi, 0)
    
            
            r_tip = euclidean_distance([branch.end_coordinates[0], 0, branch.end_coordinates[2]], [0, 0, 0])
            
            x = branch.end_coordinates[0]
            y = branch.end_coordinates[1]
            z = branch.end_coordinates[2]
            x_max = 500
            y_max  = 500
            z_max = 500
            consumption_rate = 1
            distance_falloff = 10
    
            GF_weight = 0
            neighbor_weight = arg4
            
            GF_gradient = (y_max - y)/y_max + abs(x)/x_max + abs(z)/z_max
            neighbor_gradient = 0
            neighbor_max = 0
    
            GF_scaling_y = arg1
            GF_scaling_x = arg2
            GF_scaling_z = arg3
            GF_scaling_n = arg4
    
            branch.local_GF = (GF_scaling_y*((y_max-branch.end_coordinates[1])/y_max) + GF_scaling_x*((abs(branch.end_coordinates[0]))/x_max) + GF_scaling_z*((abs(branch.end_coordinates[2]))/z_max))
            
            theta_list = []
            
            biggest_angle = 0
            biggest_average = 0
            
            proj_xz = [branch.terminal_compartment.vector[0], 0, branch.terminal_compartment.vector[2]]
            if proj_xz[2] >= 0:
                theta_list.append(math.acos(unit_vector(proj_xz)[0]))
            else:
                theta_list.append(2*math.pi-math.acos(unit_vector(proj_xz)[0]))
        
            for branch2 in branch_list:
                if branch2 != branch:
                    r = euclidean_distance(branch2.end_coordinates, branch.end_coordinates)
                    difference_vector = [branch2.terminal_compartment.vector[i] - branch.terminal_compartment.vector[i] for i in range(0, 3)]
                    
                    proj_xz = [difference_vector[0], 0, difference_vector[2]]
                    # proj_xz[1] = 0
                    if proj_xz[2] >= 0:
                        theta_list.append(math.acos(unit_vector(proj_xz)[0]))
                    else:
                        theta_list.append(2*math.pi-math.acos(unit_vector(proj_xz)[0]))
                        
                    # theta_list.append(math.acos(difference_vector[0]/r))
                                    
                    if r != 0:
                        total_difference_vector = [total_difference_vector[i] + difference_vector[i]/(r) for i in range(0, 3)]
                        # branch.local_GF = branch.local_GF + GF_scaling_n*neighbor_weights/neighbor_max
                        # branch.local_GF = branch.local_GF - neighbor_weight/((r**2) * (branch2.diameter**-2))
        
                        # neighbor_gradient = neighbor_gradient + consumption_rate*(1/(branch2.compartment_list[-1].diameter**2))*(distance_falloff/(r**2))
                        neighbor_gradient = neighbor_gradient + 1/((1+r)**2)
                        neighbor_max = neighbor_max + 1
    
                    else:
                        branch.local_GF = 0
                    
                    
                    theta_list.sort()
                                    
                    biggest_angle = theta_list[0] + 2*math.pi - theta_list[-1]
                    biggest_average = (theta_list[0] + 2*math.pi + theta_list[-1])/2
                    
                    # print biggest_angle, biggest_average
                    for k in range(0, len(theta_list) - 1):
                        current_angle = theta_list[k+1] - theta_list[k]
                        if current_angle > biggest_angle:
                            biggest_angle = current_angle
                            biggest_average = (theta_list[k+1] + theta_list[k])/2
        
            neighbor_subtotal = 0
            if neighbor_max != 0:
                neighbor_subtotal = neighbor_gradient/neighbor_max
            branch.local_GF = branch.local_GF + GF_scaling_n*neighbor_subtotal
                    
            total_difference_vector = [math.cos(biggest_average), 0, math.sin(biggest_average)]
            if biggest_average == 0:
                total_difference_vector = [0, 1, 0]
            
            jndex = int(branch.pathlength_from_soma)
            if branch.point_process[jndex] == 1:
                branching_probability = 1
            else:
                branching_probability = 0
            
            
            
            conditional_names1, conditional_names2, xedges, yedges, H = read_conditional_file(conditionals_file)   
            if branching_probability == 1 and (simulation_time-jndex) > 5 and len(branch_list) < 60:
                    
                conditional_names1, conditional_names2, xedges, yedges, H = read_conditional_file(conditionals_file)   

                branch.terminal_flag = 0
                new_vector = unit_vector([inertial_weight*branch.inertial_growth_vector[idx] for idx in range(0, 3)])
                
                parent_vector = unit_vector([branch.inertial_growth_vector[i] for i in range(0, 3)])
                rand_theta = random.gauss(0, 1)
                rand_phi = random.uniform(-180, 180)
             
                rand_theta = 0
                rand_phi = 0
                
                first_vector = new_vector
                
                if branch.diameter > 0.75:
                    taper = random.gauss(0.85, 0.03)
                    new_diameter = branch.diameter*taper
                    branch.diameter  = branch.diameter*taper
                else:
                    new_diameter = branch.diameter
                    branch.diameter = branch.diameter
                    # if (unit_vector(first_vector)[1] > 0.9999):
                if (branch.order == -2):
                    # bif_amp = draw_value('bif_amplitude_remote')
                    bif_amp = random.gauss(40, 10)
                    if rand_theta > 0:
                        rand_theta2 = rand_theta - bif_amp/2
                        rand_theta = rand_theta + bif_amp/2
                    else:
                        rand_theta2 = rand_theta + bif_amp/2
                        rand_theta = rand_theta - bif_amp/2
                    rand_phi = random.gauss(-90,10)
                    rand_phi2 = -90
                    # new_rand_vector = angle_noise([0,0,0], [-new_vector[0], new_vector[1], -new_vector[2]], rand_theta, rand_phi)
                    new_rand_vector, new_x_axis = angle_noise([0,0,0], first_vector, first_vector, branch.x_axis, rand_theta2, rand_phi2, 0)
                    new_vector, new_x_axis2 = angle_noise([0,0,0], first_vector, first_vector, branch.x_axis, rand_theta, rand_phi, 0)
                    new_branch = Branch(branch.end_coordinates, new_rand_vector, new_diameter)
                    new_branch.x_axis = new_x_axis
                else:
                    # bif_amp = draw_value('bif_amplitude_remote')
                    # heatmap = H[6]
                    # for idx in range(len(xedges[6])):
                    #     x_bin = idx - 1
                    #     if branch.pathlength_from_soma < xedges[6][idx]:
                    #         break
                    # if x_bin == -1:
                    #     x_bin = 0
                    # conditional_distribution = []
                    # for ii in range(len(H[6])):
                    #     conditional_distribution.append(H[6][ii][x_bin])
                    # cond_sum = sum(conditional_distribution)
                    # conditional_distribution2 = [0]*len(conditional_distribution)
                    # cum_sum = 0
                    # for ii in range(len(conditional_distribution)):
                    #     cum_sum = cum_sum + conditional_distribution[ii]
                    #     conditional_distribution2[ii] = cum_sum/cond_sum
                    # random_draw = random.random()
                    
                    # for idx in range(len(conditional_distribution2)):
                    #     cond_bin = idx
                    #     if random_draw < conditional_distribution2[idx]:
                    #         break
                        
                    bif_amp = conditional_probability(6, branch.pathlength_from_soma)
                    if branch.parent_bif_amp != 0:    
                        bif_amp_alt = conditional_probability(6, branch.parent_bif_amp)
                    else:
                        bif_amp_alt = bif_amp
                    
                    # bif_amp = 30
                    # bif_amp = (bif_amp + bif_amp_alt)/2
                    new_vector1, new_x_axis1 = angle_noise([0,0,0], new_vector, new_vector, branch.x_axis, 0, 0, bif_amp)
                    new_vector2, new_x_axis2 = angle_noise([0,0,0], new_vector, new_vector, branch.x_axis, 0, 0, -bif_amp)
                    
                    vert_check1 = new_vector1[1]
                    vert_check2 = new_vector2[1]
                       
                    if vert_check1 > vert_check2:
                        new_vector = new_vector1
                        new_x_axis = new_x_axis1
                        
                    else:
                        new_vector = new_vector2
                        new_x_axis = new_x_axis2  
                    
                    
                    norm_vector = unit_vector(cross_product(parent_vector, [0,1,0]))
                    if norm_vector == [0,0,0]:
                        norm_vector = [0,0,1]
                    s2 = math.sqrt(magnitude_product(new_vector, new_vector))
                    s1 = s2*math.cos(bif_amp*math.pi/180)
                    i_vector = [s1*parent_vector[i] for  i in range(0, 3)]
                    
                    p_vector = unit_vector([new_vector[i] - i_vector[i] for i in range(0,3)])

                    c_vector1 = unit_vector([p_vector[0], 0, p_vector[2]])
                    c_vector2 = unit_vector([-branch.end_coordinates[0], 0, -branch.end_coordinates[2]])
                    if c_vector2 == [0,0,0]:
                        c_vector2 = [0,0,-1]  
                    
                    # v[0] = [new_vector[0], 0, new_vector[2]]
                    # v[1] = [0, 0, 1]
    
                    # if first_vector[0] > 0 and first_vector[1] > 0.5:
                    #     rand_theta2 = rand_theta - bif_amp
                    # elif first_vector[0] < 0 and first_vector[1] < 0.5:
                    #     rand_theta2 = rand_theta - bif_amp
                    # else:
                    #     rand_theta2 = rand_theta + bif_amp

                    # rand_theta2_local = rand_theta2
    
                    # rand_phi2 = calc_angle(v[0], v[1])
                    # if v[0][0] > 0:
                    #     rand_phi2 = -rand_phi2
                    # rand_phi2 = random.gauss(-90, 20)
                    
                    previous_x = branch.end_coordinates[0]
                    previous_z = branch.end_coordinates[2]
                    
                    rand_theta2 = conditional_probability(14, previous_z)
                    rand_phi2 = conditional_probability(13, previous_x)
                    # print conditional_names1, conditional_names2
                    rand_theta2_local = conditional_probability(14, previous_z)
                    rand_phi2_local = conditional_probability(13, previous_x)
                    
                    # rand_theta2 = random.gauss(0, 10)
                    # rand_theta2_local = random.gauss(0, 10)
    

                    c_angle_remote = conditional_probability(41, previous_z)
                
                    if rand_theta2 > rand_theta2_local:
                        rand_theta2 = rand_theta2_local
                        rand_phi2 = rand_phi2_local
                    
                    azimuth = rand_phi2_local
                    azimuth = rand_phi2
                    # azimuth = 0
                    azimuth = -c_angle_remote
                    # azimuth = random.uniform(0, 180)
                    # print norm_vector, p_vector
                    # c_angle = 90 - calc_angle(norm_vector, p_vector)
                    c_angle = calc_angle(norm_vector, p_vector) + azimuth
                    # print 'bif amp =', bif_amp, ', c-angle =', c_angle
                    # c_angle = -calc_angle(norm_vector, p_vector)
                    remote_angle = calc_angle(norm_vector, p_vector) + rand_theta2
                    
                    new_vector1, new_x_axis1 = angle_noise([0,0,0], new_vector, parent_vector, new_x_axis, 0, c_angle, 0)
                    new_remote_vector1, new_remote_x_axis1 = angle_noise([0,0,0], new_vector, parent_vector, new_x_axis, 0, remote_angle, 0)
                                
                    new_vector2, new_x_axis2 = angle_noise([0,0,0], new_vector, parent_vector, new_x_axis, 0, -c_angle, 0)
                    new_remote_vector2, new_remote_x_axis2 = angle_noise([0,0,0], new_vector, parent_vector, new_x_axis, 0, remote_angle, 0)
                    
                    c_vectorA = unit_vector([new_vector1[0], 0, new_vector1[2]])
                    c_vectorB = unit_vector([new_vector2[0], 0, new_vector2[2]])
                    
                    # ang_check1 = calc_angle(c_vectorA, c_vector2)
                    # ang_check2 = calc_angle(c_vectorB, c_vector2)
                    
                    # for vertical check
                    vert_check1 = new_vector1[1]
                    vert_check2 = new_vector2[1]
                    
                    # ang_check1 = 0
                    # ang_check2 = 1
            
                    # print c_angle, c_vector2, c_vector1, c_vectorA, c_vectorB, ang_check1, ang_check2
        
                    if vert_check1 > vert_check2:
                        new_vector = new_vector1
                        new_x_axis = new_x_axis1
                        new_remote = new_remote_vector1
                        
                    else:
                        new_vector = new_vector2
                        new_x_axis = new_x_axis2      
                        new_remote = new_remote_vector2
                    
                    new_branch = Branch(branch.end_coordinates, new_vector, new_diameter)
                
                new_branch.rescale_point_process(simulation_time)
                new_branch.local_growth_vector = new_vector
                new_branch.remote_growth_vector = new_remote
                new_branch.canonical_end = [new_branch.count*new_branch.remote_growth_vector[i] + new_branch.new_compartment_list[-1].end_coordinates[i] for i in range(0, 3)]

                
                    # new_branch = Branch(branch.end_coordinates, new_vector, new_diameter)
                    
                new_branch2 = Branch(branch.end_coordinates, branch.inertial_growth_vector, new_diameter)
                new_branch2.point_process = branch.point_process
                new_branch2.local_growth_vector = branch.remote_growth_vector
                new_branch2.remote_growth_vector = branch.remote_growth_vector
                new_branch2.canonical_end = [new_branch2.count*new_branch2.remote_growth_vector[i] + new_branch2.new_compartment_list[-1].end_coordinates[i] for i in range(0, 3)]
            
                new_branch2.point_process[jndex] = 0
                if branch_list != []:
                    # new_branch2 = Branch(branch.end_coordinates, branch.inertial_growth_vector, branch.diameter)
    
                    new_branch.compartment_list[-1].parent_compartment = branch.terminal_compartment
                    new_branch.parent_branch = branch_list[-1]
    
                    new_branch2.compartment_list[-1].parent_compartment = branch.terminal_compartment
                    new_branch2.parent_branch = branch_list[-1]                
                    
                    if new_branch.parent_branch.compartment_list != []:
                        new_branch.parent_compartment = new_branch.parent_branch.compartment_list[-1]
                    else:
                        new_branch.parent_compartment = 1
                        
                    if new_branch2.parent_branch.compartment_list != []:
                        new_branch2.parent_compartment = new_branch2.parent_branch.compartment_list[-1]
                    else:
                        new_branch2.parent_compartment = 1                    
                        
                new_branch2.set_order(branch.order + 1)
                new_branch2.primary_flag = 1
                new_branch2.parent_bif_amp = bif_amp
                branch_list.append(new_branch2)
    
                new_branch.set_order(branch.order+1)
                new_branch.parent_bif_amp = bif_amp
                branch_list.append(new_branch)
                branch.terminal_compartment.vector = new_vector
            if branch.terminal_flag == 1:
                branch.grow_branch(branch_list)
    return branch_list
def external_rescale_pp(n,rescale_flag):
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
def part_1(n):
    # Generates a branching rate estimate & performs spike thinning from the inputted point process file. 
    # Input: n (point process file)

    # Establishing Constants        
    print('Part 1 Start')
    alpha = 0.25
    alpha1 = 0.5
    alpha2 = 0.25
    alpha3 = 0.5
    alpha4 = 0.5
    x_smooth = []
    T = len(n[0])
    max_EM_iterations = 5  # maximum number of iterations for EM algorithm
    
    rho_new = 0.993
    state_variance_new = 0.11949
    beta_new = 0.9999
    mu_new = -11.09

    delta = 1

    tolerance = 0.001
    max_iter2 = 100

    for ii in range(max_EM_iterations):
        # E-Step
        rho = rho_new
        state_variance = state_variance_new
        # sv_filterx = np.linspace(0, T-1, T)
        # max_sv = 10
        # sv_tau = 0.008
        # sv_shift =100
        # sv_filter = [max_sv*np.exp(-((x-sv_shift)*sv_tau)) + 1 for x in sv_filterx]
        # plt.plot(sv_filterx, sv_filter)
            
        # state_variance2 = [state_variance_new * x for x in sv_filter]
        beta = beta_new
        mu = mu_new
        print(rho, state_variance, beta, mu)
        
        if x_smooth == []:
            x_0_new = 0.1
        else:
            x_0_new = rho_new*x_smooth[1]
        var_0_new = state_variance_new/(1-rho_new**2)
        
        x1 = [0] * (T)
        sig2 = [0] * (T)
        x_smooth = [0] * (T)
        x_kkm1 = [0] * (T)
        sig2_kkm1 = [0] * (T)
        sig2_kN = [0] * (T) 
        sig_kp1kN = [0] * (T) 
        Wkkp1 = [0] * (T)
        W = [0] * (T)
        x_kIK = [0]*(T)
        var_kIK = [0]*(T)
        
        x1[0] = x_0_new
        sig2[0] = var_0_new

        for k in range(0, T):
            
            if k == 0:
                x_kkm1[0] = x_0_new
                sig2_kkm1[0] = var_0_new
            else:
                x_kkm1[k] = rho*x1[k-1]
                sig2_kkm1[k] = (rho**2)*sig2[k-1] + state_variance
                # sig2_kkm1[k] = (rho**2)*sig2[k-1] + state_variance2[k]
            
            n_sum = 0
            C = len(n)
            for c in range(0, C):
                n_sum = n_sum + beta*(n[c][k] - np.exp(mu + beta*x_kkm1[k]))
            
            x_new = x_kkm1[k]
            d_sum = 0
            f_sum = 0
            for j in range(max_iter2):
                x_old = x_new
                tolerance3 = 0.001
                delta3 = 0.001
                f_sum = 0            
                for c in range(0, C):
                    f_sum = f_sum + beta*(n[c][k] - np.exp(mu + beta*x_old))
                
                # print f_sum
                f = x_kkm1[k] + sig2_kkm1[k]*f_sum - x_old
                d_sum = -sig2_kkm1[k]*C*(beta**2)*delta*np.exp(mu + beta*x_old) - 1
                x_new = x_old - f/d_sum
                if d_sum != 0:
                    x_new = x_old - alpha*f/d_sum
                else:
                    x_new = x_old
                if abs(x_new - x_old) < tolerance3:
                    break
            
            # if k > 666:
            #     print n[c][k]
            # print k, f_sum
            # x[k] = x_kkm1[k] + sig2_kkm1[k]*n_sum
            x1[k] = x_new
            # print x_new, x_kkm1[k] + sig2_kkm1[k]*n_sum
            sig2[k] = 1/((1/sig2_kkm1[k]) + C*(beta**2)*np.exp(mu + beta*x1[k])*delta)
        x_smooth[T-1] = x1[T-1]
        W[T-1] = sig2[T-1] + x_smooth[T-1]**2
        sig2_kN[T-1] = sig2[T-1]
        for k in range(T-2, -1, -1):
            s_k = rho*sig2[k]/sig2_kkm1[k+1]
            x_smooth[k] = x1[k] + s_k*(x_smooth[k+1] - x_kkm1[k+1])
            x_kIK[k] = x_smooth[k]
            sig2_kN[k] = sig2[k] + (s_k**2)*(sig2_kN[k+1] - sig2[k+1])
            var_kIK[k] = sig2_kN[k]
            if k == T-2:
                sig_kp1kN[k] = s_k*sig2_kN[T-1]
            else:
                # sig_kp1kN[k] = s_k*(sig_kp1kN[k+1])
                sig_kp1kN[k] = s_k*sig2_kN[k+1]
            
            Wkkp1[k] = sig_kp1kN[k] + x_smooth[k]*x_smooth[k+1]
            W[k] = sig2_kN[k] + x_smooth[k]**2
        # print W[-2]
        W0 = W[0]
        Wkkp10 = Wkkp1[0]        
        # M-Step
        # u, B, N_c, k2, x, var = symbols('u B N_c k2 x var')
        N_c = np.sum(n[c])
        tolerance2 = 0.001
        # M-Step Newton's method
        B = [0]*max_iter2
        B[0] = beta
        f_array = [0]*max_iter2
        for i in range(max_iter2-1):
            f_sum = 0
            f2_sum = 0
            u_sum = 0
            n_sum = 0
            u2_sum = 0
            du_sum = 0
            delta2 = 0.001
            du_sum2 = 0
            # B[i] = i/10
            # print x_kIK
            for k in range(0, T):
                u_sum = u_sum + np.exp(B[i]*x_kIK[k] + (1/2)*(B[i]**2)*var_kIK[k])*delta
                # u2_sum = u2_sum + np.exp((B[i]+delta2)*x_kIK[k] + (1/2)*((B[i]+delta2)**2)*var_kIK[k])*delta
                f_sum = f_sum + (np.exp(B[i]*x_kIK[k] + (1/2)*(B[i]**2)*var_kIK[k]))*(x_kIK[k] + B[i]*var_kIK[k])*delta
                # f2_sum = f2_sum + (np.exp((B[i]+delta2)*x_kIK[k] + (1/2)*(B[i]+delta2)*var_kIK[k]))*(x_kIK[k] + (B[i]+delta2)*var_kIK[k])*delta
                n_sum = n_sum + n[c][k]*x_kIK[k]
                du_sum = du_sum + (np.exp(x_kIK[k] + B[i]*var_kIK[k]))*delta
                du_sum2 = du_sum2 + x_kIK[k] + B[i]*var_kIK[k]
                # print B[i], x_kIK[k], var_kIK[k]
            
            u = np.log(N_c) - np.log(u_sum)
            # u2 = np.log(N_c) - np.log(u2_sum)
            # f = np.exp(u*(f_sum*(x_kIK[k] + B[i]*var_kIK[k])*delta)) - n_sum
            # f2 = np.exp(u*(f2_sum*(x_kIK[k] + (B[i]+delta2)*var_kIK[k])*delta)) - n_sum
            f = np.exp(u)*f_sum - n_sum
            # f2 = np.exp(u[1])*f2_sum - n_sum
            f_array[i] = f
            # df1 = (f2-f)/delta2
            
            # du = -1/(u_sum)
            du = -du_sum/u_sum
            # du[1] = -1/u[1]_sum
            d_f_sum = 0
            for k in range(0, T):
                left_side = (np.exp(B[i]*x_kIK[k] + (1/2)*B[i]*var_kIK[k]))
                right_side = (x_kIK[k] + B[i]*var_kIK[k])*delta
                d_left_side = (x_kIK[k] + var_kIK[k]/2)*(np.exp(B[i]*x_kIK[k] + (1/2)*B[i]*var_kIK[k]))
                d_right_side = var_kIK[k]*delta
                d_f_sum = d_f_sum + d_left_side*right_side + left_side*d_right_side
            
            df = (-(N_c*du_sum2)/u_sum)*f_sum + (N_c/u_sum)*d_f_sum
            # df = df1
            # print df1, df
                
            # print ii, i, B[i]
            # print u_sum
            if abs(df) > 1e-8:
                B[i+1] = B[i] - alpha*f/df
            else:
                B[i+1] = B[i]
            if abs(B[i+1] - B[i]) < tolerance2:
                break
        
        # plt.plot(B, f_array)
        # print B[i]
        
        
        #     likelihood calculation
        # for k in range(0, T):
        #     n[c]
        
        sum_a = 0
        sum_b = 0
        sum_c = 0
        for k in range(0, T):
            if k == 0:
                sum_b = sum_b + W0
                sum_a = sum_a + Wkkp10
            else:
                sum_a = sum_a + Wkkp1[k-1]
                sum_b = sum_b + W[k-1]
            sum_c = sum_c + W[k]
        rho_new = sum_a/sum_b
        state_variance_new = (1/T)*(sum_c + (rho_new**2)*sum_b - 2*rho_new*sum_a + W0*(1-rho_new**2))
        beta_new = B[i+1]
        mu_new = u
        
        
        rho_diff = rho_new - rho
        state_diff = state_variance_new - state_variance
        beta_diff = beta_new - beta
        mu_diff = mu_new - mu
        
        # if mu_diff*alpha4 < tolerance2*10:
        #     alpha4 = alpha4*0.5
        
        rho_new = rho + alpha1*rho_diff
        state_variance_new = state_variance + alpha2*state_diff
        beta_new = beta + alpha3*beta_diff
        mu_new = mu + alpha4*mu_diff
        # rho_new = 0.946
        # state_variance_new = 0.26
        # beta_new = 2.485
        # mu_new = -9.53 
        
        rho_new = 0.993
        mu_new = -11.09
        
        diff1 = abs(rho_new - rho)
        diff2 = abs(state_variance_new - state_variance)
        diff3 = abs(beta_new - beta)
        diff4 = abs(mu_new - mu)
        
        if diff1 < tolerance and diff2 < tolerance and diff3 < tolerance and diff4 < tolerance:
            break
        
        x_0_new = rho_new*x_kIK[1]
        var_0_new = state_variance_new/(1-rho_new**2)

    p_max = 700
    
    global max_poisson
    
    # lamb = [np.exp(mu + gain*x_s) for x_s in x_smooth]
    global lamb
    global lamb_tmp
    global tau


    coef1 = 0.8
    tau = 200
    
    sv_filterx = np.linspace(0, T-1, T)
    max_sv = 10
    sv_tau = 0.008
    sv_shift = 0
    sv_filter = [max_sv*np.exp(-((x-sv_shift)*sv_tau)) + 1 for x in sv_filterx]
    # plt.plot(sv_filterx, sv_filter) 
    
    # coef1 = [1/x for x in sv_filter]
    # coef1 = [1]*len(sv_filter)
    lamb = [0]*len(x_smooth)
    for kdx in range(0, len(x_smooth)):
        # sv_filter[kdx] = 1
        lamb[kdx] = np.exp(mu_new + beta_new*x_smooth[kdx])/sv_filter[kdx]
    lamb = [np.exp(mu_new + beta_new*x_s)*coef1 for x_s in x_smooth]
    lamb_tmp = list(lamb)
    lambda_max = max(lamb)
    lambda_min = min(lamb)
    # lamb = lamb/lambda_max
    lamb2 = list(lamb)
    k = 0
    max_poisson = [[]]*p_max
    T_max = max_N
    
    delta = 1
    # tau = 250
    rate_time = np.linspace(0, int(len(lamb))-1, int(len(lamb)/delta))
    rate_res = delta
    
    L = int(10*tau/rate_res)
    refract = 1-np.exp(-rate_time[:L]/tau)
    questionable = 1
    weighting = 1
    for i in range(0, p_max):
        lamb2 = list(lamb)
        max_poisson[i] = [0]*(max_N+1)   
        k = 0
        poisson_flag = 0
        while poisson_flag == 0:
            max_poisson[i] = [0]*(max_N+1)    
            k = 0
            while k < T_max-1:
                u_k = np.random.uniform(0, 1)
                w_k = -math.log(u_k)/(lambda_max)
                # z_k = np.random.exponential(lambda_max)
                # k = k + 1 + int(z_k)
                k = k + 1 + int(w_k)
                if k >= T_max:
                    max_poisson[i][-1] = 0
                else:
                    max_poisson[i][k] = 1
            
            for j in range(0, len(max_poisson[i])):
                if max_poisson[i][j] == 1:
                    # p = lamb2[j]/lambda_max
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
            # print poisson_sum
            if poisson_sum < 1 or poisson_sum > 6:
                poisson_flag = 1
            else:
                poisson_flag = 0
    poisson_sums = [0]*len(max_poisson)
    for i in range(len(max_poisson)):
        poisson_sums[i] = sum(max_poisson[i])
    plt.figure()
    plt.plot(lamb2, color='black')
    plt.title("Estimated Branching Rate", fontsize=18)
    plt.xlabel("Pathlength from Soma (um)", fontsize=18)
    plt.ylabel("Branching Rate (branch/um)", fontsize=18)
    #raster plot of generated point processes  
    plt.figure()
    n_subplot = p_max-1
    to_be_rastered = []
    for i in range(0,n_subplot):
        event_index_list = [x for x, y in enumerate(max_poisson[i]) if y == 1]
        to_be_rastered.append(event_index_list)
    ax = raster(to_be_rastered)
    plt.xlabel('Distance of Bifurcation Occurrence (um)', fontsize=18)
    plt.ylabel('Trial #', fontsize=18)
    plt.show()

    # Returning max_poisson pp list of pps
    # return lamb2
def part_2(n, part2_flag):
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
def part_3(conditionals_file, part3_flag,type_marker=''):
    # # Morphologies constructed using conidtional branching angles.
    if part3_flag == 1:
        # Establishing Constants
        arg1 = 0
        arg2 = 0
        arg3 = 0
        arg4 = 0
        arg5 = 0
        global c_id
        global max_prob
        c_id = 1

        num_morphologies = 1
        plot_flag = 1
        flag_val = 0
        weight_bonus = 0.05
        exponent_val = 3
        max_stems = 3
        fail_count = 0

        ## 
        view1 = 15
        view2 = 45

        simulation_time = int(300 + random.random()*130)
        stem_pathlength_coef = (0.49, 1.49, 1.8)

        inertial_weight = 1
        neighbor_weight = 0
        gradient_weight = 0.02
        random_weight = 1
        max_prob = 0

        growing_inertial_weight = 1
        growing_neighbor_weight = 0
        growing_gradient_weight = 0.01
        growing_random_weight = 1

        bif_amplitude_local_coef = (5.51, -0.83, 9.275)
        bif_amplitude_remote_coef = (3.51, 0.17, 9.44)

        origin = [0, 0, 0]
        stem_diameter = 1
        stem_vector = unit_vector([0, 1, 0])
        branch_list = []

        theta_sd = 10
        phi_sd = 180

        output_titles = []
        
        for j in range(0, num_morphologies):
            flag = 1
            while flag == 1:
                c_id = 1
                simulation_time = -1
                while simulation_time > max_N or simulation_time < 0:
                    simulation_time = int(250 + random.gauss(100,50))   
                    
                n_stem = int(random.random()*max_stems) + 1
                print(('morphology #', j+1))
            
                if n_stem == 1:
                    origin = [0, 0, 0]
                    
                    theta = random.uniform(-theta_sd, theta_sd)
                    phi = random.uniform(-phi_sd, phi_sd)
            
                    stem_diameter = random.gauss(1.5, 0.2)
                    stem_vector = unit_vector([0, 1, 0])
                    stem_x_axis = [1,0,0]
                    new_stem_vector, new_x_axis = angle_noise(origin, stem_vector, stem_vector, stem_x_axis, theta, phi, 0)
                    branch_list = []
                    stem = Branch(origin, new_stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)
                    branch_list.append(stem)
                elif n_stem == 2:
                    origin = [0, 0, 0]
                    stem_diameter = 1.5
                    bif_amp = random.gauss(40, 20)
                    rand_theta = random.gauss(0, theta_sd)
                    rand_theta = 0
                    if rand_theta > 0:
                        rand_theta2 = rand_theta - bif_amp/2
                        rand_theta = rand_theta + bif_amp/2
                    else:
                        rand_theta2 = rand_theta + bif_amp/2
                        rand_theta = rand_theta - bif_amp/2
                        
                    rand_phi = -90
                    
                    rand_phi2 = rand_phi
                    stem_x_axis = [1,0,0]
            
                    stem_vector = unit_vector([0, 1, 0])
                    new_stem_vector, new_x_axis = angle_noise(origin, stem_vector, stem_vector, stem_x_axis, 0, -90, rand_theta)
                    branch_list = []
                    stem = Branch(origin, new_stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
            
                    stem_vector = unit_vector([0, 1, 0])
                    new_stem_vector, new_x_axis, = angle_noise(origin, stem_vector, stem_vector, stem_x_axis, 0, -90, -rand_theta2)
                    stem = Branch(origin, new_stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem) 
                elif n_stem == 3:
                    origin = [0, 0, 0]
                    stem_diameter = 1.5
                    
                    rand_theta = random.gauss(0, theta_sd)
                    rand_phi = random.uniform(-phi_sd, phi_sd) - 90

                    bif_amp = draw_value('bif_amplitude_local')
                    bif_amp = random.gauss(40, 20)

                    if rand_theta == 0:
                        theta_diff = bif_amp
                    else:
                        theta_diff = 180*bif_amp/(2*rand_theta)
                    if rand_theta > 0:
                        rand_theta2 = rand_theta - bif_amp
                    else:
                        rand_theta2 = rand_theta + bif_amp
                    rand_phi2 = -90
                    
                    bif_amp2 = draw_value('bif_amplitude_local')
                    bif_amp2 = random.gauss(40, 20)

                    if rand_theta > 0:
                        rand_theta3 = rand_theta - bif_amp2
                    else:
                        rand_theta3 = rand_theta + bif_amp2
                    rand_phi3 = rand_phi
                    stem_x_axis = [1,0,0]
                    
                    stem_vector = unit_vector([0, 1, 0])
                    new_stem_vector, new_x_axis = angle_noise(origin, stem_vector, stem_vector, stem_x_axis, rand_theta, rand_phi, 0)
                    branch_list = []
                    stem = Branch(origin, new_stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
                    
                    stem_vector = unit_vector([0, 1, 0])
                    new_stem_vector, new_x_axis = angle_noise(origin, stem_vector, stem_vector, stem_x_axis, rand_theta2, rand_phi2, 0)
                    stem = Branch(origin, new_stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
            
                    stem_vector = unit_vector([0, 1, 0])
                    new_stem_vector, new_x_axis = angle_noise(origin, stem_vector, stem_vector, stem_x_axis, rand_theta3, rand_phi3, 0)
                    stem = Branch(origin, new_stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
                elif n_stem == 4:
                    origin = [0, 0, 0]
                    stem_diameter = 1.5
                    theta = random.uniform(-theta_sd, theta_sd)
                    phi = random.uniform(-phi_sd, phi_sd)    
                    stem_vector = unit_vector([0.5, 1, 0.5])
                    branch_list = []
                    stem = Branch(origin, stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
            
                    stem_vector = unit_vector([0.5, 1, -0.5])
                    stem = Branch(origin, stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
            
                    stem_vector = unit_vector([-0.5, 1, -0.5])
                    stem = Branch(origin, stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
                    
                    stem_vector = unit_vector([-0.5, 1, 0.5])
                    stem = Branch(origin, stem_vector, stem_diameter)
                    stem.x_axis = new_x_axis
                    stem.rescale_point_process(simulation_time)

                    branch_list.append(stem)
                for t in range(0, simulation_time):
                    branch_list = simulate_growth(branch_list, arg1, arg2, arg3, arg4, arg5, t, simulation_time-1,conditionals_file)
                b_list = branch_list
                
                sholl_bins = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0]
                sholl_avgs = [0.35931748573555294, 0.5462821744859561, 0.7400286882429741, 0.8722411536487171, 0.8449991865433044, 0.6591798626987704, 0.25765960230245943, 0.06795795241173393, 0.02768380079304449, 0.017316017316017316]
                sholl_counts = [0]*len(sholl_bins)
                for branch in b_list:
                    for compartment in branch.compartment_list:
                        if compartment.parent_compartment != []:
                            for i in range(0, len(sholl_bins)):
                                efs1 = math.sqrt(sum([compartment.end_coordinates[k]**2 for k in range(0, 3)]))
                                efs2 = math.sqrt(sum([compartment.start_coordinates[k]**2 for k in range(0, 3)]))
                                if efs1 >= sholl_bins[i] and efs2 < sholl_bins[i]:
                                    sholl_counts[i] = sholl_counts[i] + 1
                sholl_adjusted = [x/max(sholl_counts) for x in sholl_counts]
                sholl_error = sum((sholl_avgs[i] - sholl_adjusted[i])**2 for i in range(0, len(sholl_counts)))
                
                bif_num_min = 12
                bif_num_max = 45
                bif_num_mid = (bif_num_min + bif_num_max)/2
                bif_num_range = bif_num_max - bif_num_min
                
                bif_num = len(b_list)
                z = (bif_num - bif_num_mid)
                sigma = 10
                gauss_max = np.exp(-(0)**2/2)/(sigma*math.sqrt(2*math.pi))
                gaussian_filter = np.exp(-(-z/sigma)**2/2)/(sigma*math.sqrt(2*math.pi))/gauss_max
                if random.uniform(0,1) > gaussian_filter:
                    # print gaussian_filter
                    flag = flag_val
                    fail_count = fail_count + 1
                elif sholl_error > 2:
                    # print sholl_error
                    flag = flag_val
                    fail_count = fail_count + 1
                else:
                    flag = 0
            
            sum_vector = [0, 0, 0]
            for branch in branch_list:
                if branch.terminal_flag != -1:
                    for i in range(0, 3):
                        sum_vector[i] = sum_vector[i] + branch.full_vector[i]
            sum_array = np.asarray(sum_vector)
            sum_mag = np.sqrt(sum_array.dot(sum_array))
            sum_unit = np.ndarray.tolist(sum_array/sum_mag)
            rot_x, rot_y, rot_z = rotation_angles(sum_unit, [0, 0, 1])
            rot_x = rot_x*180/math.pi
            rot_y = rot_y*180/math.pi
            rot_z = rot_z*180/math.pi

            rotated_sum_vector, throw_away_axis = angle_noise([0,0,0], sum_vector, sum_vector, [1,0,0], -rot_x, rot_y, rot_z)
            
            for branch in b_list:
                for compartment in branch.compartment_list:
                    if compartment.start_coordinates != [0,0,0]:
                        compartment.start_coordinates, throw_away_axis = angle_noise([0, 0, 0], compartment.start_coordinates, [0,1,0], [1,0,0], -rot_x, rot_y, rot_z)
                    compartment.end_coordinates, throw_away_axis = angle_noise([0, 0, 0], compartment.end_coordinates, [0,1,0], [1,0,0], -rot_x, rot_y, rot_z)
                branch.start_coordinates = branch.compartment_list[0].start_coordinates
                branch.end_coordinates = branch.compartment_list[-1].end_coordinates
            
            max_index = 0
            max_distance = 0        
            for ii in range(0, len(b_list)):

                if b_list[ii].terminal_flag != -1:
                    if b_list[ii].end_coordinates[0]**2 + b_list[ii].end_coordinates[2]**2 > max_distance:
                        max_distance = b_list[ii].end_coordinates[0]**2 + b_list[ii].end_coordinates[2]**2
                        max_index = ii
                        # print max_index, max_distance
            ref_coord = b_list[max_index].end_coordinates
            max_index = 0
            max_distance = 0    
            for jj in range(0, len(b_list)):

                if b_list[jj].terminal_flag != -1:
                    if (ref_coord[0] - b_list[jj].end_coordinates[0])**2 + (ref_coord[2] - b_list[jj].end_coordinates[2])**2 > max_distance:
                        max_distance = (ref_coord[0] - b_list[jj].end_coordinates[0])**2 + (ref_coord[2] - b_list[jj].end_coordinates[2])**2
                        max_index = jj
                    
            ref_coord2 = b_list[max_index].end_coordinates
            
            ref_difference = [ref_coord[kk] - ref_coord2[kk] for kk in range(0,3)]
            ref_difference[1] = 0

            final_rot_x, final_rot_y, final_rot_z = rotation_angles(ref_difference, [0,0,1])
            
            final_rot_z = final_rot_z * 180/math.pi

            for branch in b_list:
                for compartment in branch.compartment_list:
                    if compartment.start_coordinates != [0,0,0]:
                        compartment.start_coordinates, throwaway_axis = angle_noise([0, 0, 0], compartment.start_coordinates, [0,1,0], [1,0,0], 0, -final_rot_z, 0 )
                    compartment.end_coordinates, throwaway_axis = angle_noise([0, 0, 0], compartment.end_coordinates, [0,1,0], [1,0,0], 0, -final_rot_z, 0)
                branch.start_coordinates = branch.compartment_list[0].start_coordinates
                branch.end_coordinates = branch.compartment_list[-1].end_coordinates    

            output_title = '../Data/Morphologies/test_output/{}{}{}{}'.format('morphology_tester_output_',type_marker,'_', j+1, '.swc')
            output_titles.append(output_title) # Potentially for saving separate files. 
            soma_diameter = random.gauss(2.8, 0.3)
            
            with open(output_title, 'w') as f:
                true_origin = b_list[0].compartment_list[0]
                line = (true_origin.compartment_id, 1, true_origin.start_coordinates[0], true_origin.start_coordinates[1], true_origin.start_coordinates[2], soma_diameter, -1)
                parent_id = 1
                strs=' '.join(str(x) for x in line)
                f.write(strs + '\n') # Writes origin point to file
                
                # Post Note: Make sure you have the else-if syntax corret when you get back. 
                section_type = 0 # Accounting for either case
                if type_marker == 'apical':   
                    section_type = 4
                elif type_marker == 'basal':
                    section_type = 3
                else:
                    section_type = 0

                delete_flag = 0
                for branch in b_list:
                    for compartment in branch.compartment_list:
                        if compartment != []:
                            compartment.compartment_id = c_id+1
                            c_id = c_id + 1
                            if compartment.parent_compartment != []:
                                parent_id = compartment.parent_compartment.compartment_id
                            else:
                                parent_id = 1
                            new_line = (compartment.compartment_id, section_type, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter, parent_id)
                            if new_line == line:
                                delete_flag = 1
                            line = (compartment.compartment_id, section_type, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter, parent_id)
                            strs=' '.join(str(x) for x in line)
                            f.write(strs + '\n') # Writes each compartment to the file. (Needs correction to account for apical and basal types (3 or 4))
                        parent_id = compartment.compartment_id
            f.closed
        
        if  delete_flag == 1:
            os.remove(output_title)
        print(('Success Rate:', num_morphologies/(num_morphologies + fail_count)))
        
        if plot_flag == 1:
            plot_branches(b_list)
            plot_branches(branch_list)
       
        return output_titles
    else:
        return 0
