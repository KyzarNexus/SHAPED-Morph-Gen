#%%
from __future__ import division
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import os
import pickle

from angle_noise_updated import angle_noise
from B_morph_170811_2 import make_compartments

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
    Rr_y = np.matrix( ((math.cos(-rot_y), 0, -math.sin(-rot_y)), (0, 1, 0), (math.sin(-rot_y), 0, math.cos(-rot_y))) )

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

def read_file(file_name):

    point_process_list = []

    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()
    
    for i in range(0, len(lines)):
        parsed_line = lines[i].split(' ')
        point_process_list.append([0] + list(map(int,parsed_line)))
    
    return point_process_list

def read_conditional_file(file_name):

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
     
        xedges.append(map(float, parsed_line[0].split(' ')))
        del parsed_line[0]
        yedges.append(map(float, parsed_line[0].split(' ')))
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

def euclidean_distance(u, v):
    ans = math.sqrt((u[0] - v[0])**2 + (u[1] - v[1])**2 + (u[2] - v[2])**2)
    return ans

def unit_vector(u):
    magnitude = euclidean_distance((u[0], u[1], u[2]), (0, 0, 0))
    
    if magnitude != 0:  
        return [u[0]/magnitude, u[1]/magnitude, u[2]/magnitude]
    else:
        return [0, 0, 0]

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

def plot_branches(branch_list):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for branch in branch_list:
        x_list.append(branch.start_coordinates[0])
        y_list.append(branch.start_coordinates[1])
        z_list.append(branch.start_coordinates[2])        
        for compartment in branch.compartment_list:
            x_list.append(compartment.end_coordinates[0])
            y_list.append(compartment.end_coordinates[1])
            z_list.append(compartment.end_coordinates[2])
        plt.plot(x_list, z_list, y_list, color='black')
        
        x_list = []
        y_list = []
        z_list = []
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_xlim([-200,200])
    ax.set_ylim([-200,200])
    ax.set_zlim([0,400])
  
    # ax.set_xlim([-20,20])
    # ax.set_ylim([-20,20])
    # ax.set_zlim([0,75])
    
    ax.view_init(elev=view1, azim=view2)
    plt.show()
      
class Compartment:
    def __init__(self, start_coordinates, end_coordinates, diameter):
        
        # global c_id
        self.compartment_id = 1
        self.parent_compartment = []
        # c_id = c_id + 1
        
        self.start_coordinates = start_coordinates
        self.end_coordinates = end_coordinates
        self.diameter = diameter
        self.jitter_theta = 0
        self.jitter_phi = -90
        self.counter = 1
        
        self.vector = [end_coordinates[i] - start_coordinates[i] for i in range(0, 3)]
        self.length = euclidean_distance(self.start_coordinates, self.end_coordinates)

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
        

        
        # global p_max
        # global max_poisson
        # p_flag = 0
        # while p_flag == 0:
        #     self.point_process_id = random.randrange(0, p_max)
        #     if sum(max_poisson[self.point_process_id]) != 0:
        #         if sum(max_poisson[self.point_process_id][0:75]) != -2:           
        #             p_flag = 1


            
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
        questionable = 1
        
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
                z_k = np.random.exponential(lambda_max)
                # k = k + 1 + int(z_k)
                k = k + 1 + int(w_k)
                if k >= max_N:
                    max_poisson2[-1] = 0
                else:
                    max_poisson2[k] = 1
            
            for j in range(0, len(max_poisson2)):
                if max_poisson2[j] == 1:
                    # p = lamb2[j]/lambda_max
                    p = lamb_tmp2[j]/(scalar1*max(lamb_tmp2))

                    b = np.random.binomial(1, p)
                    if b == 0:
                        max_poisson2[j] = 0
                    else:
                        aa = 0
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
            if sum(max_poisson2[0:75]) == 0:
                p_flag = 1
                

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

def calc_angle(vector1, vector2):
    if sum([(vector1[i] - vector2[i])**2 for i in range(0, 3)]) < 0.001:
        angle = 0
    elif sum([(-vector1[i] - vector2[i])**2 for i in range(0, 3)]) < 0.001:
        angle = 180
    else:
        angle = math.acos(dot_product(vector1, vector2)/magnitude_product(vector1, vector2))*180/math.pi
    return angle

def conditional_probability(index, conditions):
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

def simulate_growth(branch_list, arg1, arg2, arg3, arg4, arg5, t, simulation_time):
    global max_prob
    theta_sd = 0
    phi_sd = 0
    for branch in branch_list:

        
        
        # if len(branch_list) > 100:
        #     branch.terminal_flag = 0
        if branch.terminal_flag == 1:
            total_difference_vector = [0, 0, 0]
            boundary = [0, max_N, 0]
            # gradient_vector = unit_vector([boundary[i] - branch.end_coordinates[i] for i in range(0, 3)])
            gradient_vector = [0, 1, 0]
            
            # random_angle = draw_value('random_angle')
            # random_vector = unit_vector([math.cos(random_angle), 0, math.sin(random_angle)])
     
            rand_theta = random.gauss(0, theta_sd)
            rand_phi = random.gauss(0, phi_sd)
            # rand_phi = random.uniform(-phi_sd, phi_sd)
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
            # neighbor_weight = 0.06
    
            GF_gradient = (y_max - y)/y_max + abs(x)/x_max + abs(z)/z_max
            neighbor_gradient = 0
            neighbor_max = 0
    
            GF_scaling_y = arg1
            GF_scaling_x = arg2
            GF_scaling_z = arg3
            GF_scaling_n = arg4
    
            branch.local_GF = (GF_scaling_y*((y_max-branch.end_coordinates[1])/y_max) + GF_scaling_x*((abs(branch.end_coordinates[0]))/x_max) + GF_scaling_z*((abs(branch.end_coordinates[2]))/z_max))
            # branch.local_GF = (GF_scaling_y*((branch.end_coordinates[1])/y_max) + GF_scaling_x*((abs(branch.end_coordinates[0]))/x_max) + GF_scaling_z*((abs(branch.end_coordinates[2]))/z_max))
    
            # print branch.end_coordinates[1], branch.local_GF
            theta_list = []
            
            biggest_angle = 0
            biggest_average = 0
            
            proj_xz = [branch.terminal_compartment.vector[0], 0, branch.terminal_compartment.vector[2]]
            # proj_xz[1] = 0
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
                    
            # branch.local_GF = GF_gradient*GF_weight + neighbor_gradient*neighbor_weight
            # print branch.local_GF, neighbor_gradient
                
            total_difference_vector = [math.cos(biggest_average), 0, math.sin(biggest_average)]
            if biggest_average == 0:
                total_difference_vector = [0, 1, 0]
            # total_difference_vector = unit_vector(total_difference_vector)            
            #     if branch.local_GF > draw_value('GF_threshold') and random.random() < (branch.local_GF - 5)/1.99:
        
            #     print 1+branch.local_GF
            #     if random.random() < float(branch.diameter**2)/(3000*branch.local_GF):
            #     print 1+branch.local_GF*(branch.diameter**-2)
                
    
    
            # index = branch.point_process_id
            jndex = int(branch.pathlength_from_soma)
            if branch.point_process[jndex] == 1:
                branching_probability = 1
            else:
                branching_probability = 0
            # max_prob = weight_bonus*branch.local_GF/(GF_scaling_x + GF_scaling_y + GF_scaling_z + GF_scaling_n)
            # lamb = np.exp(x_smooth[int(branch.pathlength)])
            # branching_probability = weight_bonus * (lamb)/(2**branch.order)
            # if max_prob > 0.05:
            #     max_prob = 0.05
                
            # if branching_probability > max_prob:
            #     max_prob = branching_probability
            conditionals_file = '../Data/Metrics/claiborne_conditionals.txt'   
 
            conditional_names1, conditional_names2, xedges, yedges, H = read_conditional_file(conditionals_file)   
            # if random.random() < branching_probability:
            if branching_probability == 1 and (simulation_time-jndex) > 5 and len(branch_list) < 60:
                # global lamb_tmp
                # global lamb
                
                # delta = 1
                # tau = 200
                # rate_time = np.linspace(0, int(len(lamb))-1, int(len(lamb)/delta))
                # rate_res = delta
                
                # L = int(10*tau/rate_res)
                # refract = 1-np.exp(-rate_time[:L]/tau)
                # questionable = 1                
                # j = simulation_time

                # lamb_tmp = list(lamb)

                # if (j+L) > len(lamb_tmp):
                #     lamb_tmp[j:] *= refract[:len(lamb_tmp)-j]
                # else:
                #     lamb_tmp[j:j+L] *= refract
                    
                conditional_names1, conditional_names2, xedges, yedges, H = read_conditional_file(conditionals_file)   

                branch.terminal_flag = 0
                # print index, jndex, max_poisson[index][jndex]
                # if (branch.local_GF < 2) and (random.random() < 0.01):
                
                new_vector = unit_vector([inertial_weight*branch.inertial_growth_vector[idx] for idx in range(0, 3)])
                # rand_theta = random.gauss(20, theta_sd)
                parent_vector = unit_vector([branch.inertial_growth_vector[i] for i in range(0, 3)])
                rand_theta = random.gauss(0, 1)
                rand_phi = random.uniform(-180, 180)
             
                rand_theta = 0
                rand_phi = 0
                
                first_vector = new_vector
                # new_vector = angle_noise([0,0,0], [new_vector[0], new_vector[1], new_vector[2]], rand_theta, rand_phi, 0)
    
                # print new_vector
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
                    
                    v[0] = [new_vector[0], 0, new_vector[2]]
                    v[1] = [0, 0, 1]
    
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
                    
                    
                    
                    # rand_theta2_local = 45
                    # rand_phi2_local = -90

                    # rand_theta2 = rand_theta2_local
                    # rand_phi2 = rand_phi2_local
                    # rand_theta2_local = rand_theta2
                    # rand_phi2_local = rand_phi2

                    # rand_theta2 = random.gauss(previous_z/4, 5)
                    # rand_phi2 = random.gauss(previous_x/5, 10)
                        
                    # rand_phi2 = random.gauss(-90, 5/(branch.end_coordinates[1]/100))
                    # new_rand_vector = angle_noise([0,0,0], [-new_vector[0], new_vector[1], -new_vector[2]], rand_theta, rand_phi)
                    # new_rand_vector = angle_noise([0,0,0], [new_vector[0], new_vector[1], new_vector[2]], rand_theta2, rand_phi2, 0)
                    # new_local_vector = angle_noise([0,0,0], [new_vector[0], new_vector[1], new_vector[2]], rand_theta2_local, rand_phi2_local, 0)
                    
                    # ang_check1 = calc_angle(new_rand_vector, [0, 1, 0])
                    # ang_check2 = calc_angle(first_vector, [0, 1, 0])
                    
                    # local_ang_check1 = calc_angle(new_local_vector, [0, 1, 0])
                    
                    # if local_ang_check1 > ang_check2:
                    #     new_local_vector = angle_noise([0,0,0], [new_vector[0], new_vector[1], new_vector[2]], rand_theta2_local, -rand_phi2_local, 0)
                    #     new_rand_vector = angle_noise([0,0,0], [new_vector[0], new_vector[1], new_vector[2]], rand_theta2, -rand_phi2, 0)
    
                    # if branch.parent_branch != []:
                    #     second_flag = 1
                    # elif branch.parent_branch.primary_flag == 2:
                    #     second_flag = 1
                    # else:
                    #     second_flag = 0
                    # if ang_check1 > ang_check2:
                    #     new_rand_vector = angle_noise([0,0,0], [new_vector[0], new_vector[1], new_vector[2]], rand_theta2, -rand_phi2, 0)
                    
                    new_branch = Branch(branch.end_coordinates, new_vector, new_diameter)
                # print bif_amp, c_angle, vert_check1, vert_check2, branch.inertial_growth_vector[1], branch.end_coordinates[1]
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
    
                
    
    
                # bif_angle = calc_angle(new_rand_vector, new_vector)
    
                new_branch.set_order(branch.order+1)
                new_branch.parent_bif_amp = bif_amp
                branch_list.append(new_branch)
                branch.terminal_compartment.vector = new_vector
            #         print new_vector, branch.terminal_compartment.vector
                        
            #         if unit_vector(branch.terminal_compartment.vector)[1] > 0.95:
            #             rand_theta = random.gauss(0, theta_sd)
            #             rand_phi = random.gauss(0, phi_sd)
            # #                new_rand_vector = angle_noise([0,0,0], [-new_vector[0], new_vector[1], -new_vector[2]], rand_theta, rand_phi)
            #             new_rand_vector = angle_noise([0,0,0], [0, 1, 0], rand_theta, rand_phi)
            
            #             branch.terminal_compartment.vector = new_rand_vector
            #     else:
            #         if t == simulation_time:
            #             branch_list.remove(branch)
            
            #     if branch.terminal_flag == 1:
            if branch.terminal_flag == 1:
                branch.grow_branch(branch_list)
            # else:
            #     new_branch2.grow_branch(branch_list)
            # branch.grow_branch(branch_list)
    return branch_list


file_name = '..\\..\\Data\\Metrics\\amaral_point_process_basal.txt'

#file_name = '/home/zzchou/Desktop/claiborne_point_process.txt'

#file_name = '\Users\Zane\Desktop\Morphologies\Codes\claiborne_point_process.txt'
#file_name = '\Users\Zane\Desktop\Morphologies\Codes\claiborne_point_process.txt'

n = read_file(file_name)
max_N = len(n[0])
max_R = len(n)

# Rescale point processes
for r in range(0, max_R):
    # Calculate scaling factor

    list_of_indices = [i for i, x in enumerate(n[r]) if x == 1]
    termination = list_of_indices[-1]
    scaling_factor = (max_R-1)/termination
    
    # Rescaled point process indices
    list_of_indices_scaled = [int(x*scaling_factor) for x in list_of_indices]
    new_point_process = [0]*(max_R)
    for idx in range(0, len(list_of_indices_scaled)):
        new_point_process[list_of_indices_scaled[idx]] = 1
    n[r] = new_point_process

# n = [[0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
#     [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
#     [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
#     [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
#     [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
#     [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
#     [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
#     [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0]]

one_count = 0
for row in n:
    for column in row:
        if column == 1:
            one_count = one_count+1
total_count = len(n) * (len(n[0]) - 1)
mean_branching_rate = one_count/total_count

fail_count = 0
state_noise_variance = 0.015
gain = 1
history_gain = 1
self_gain = -1
mu = np.log(mean_branching_rate)
max_S = 30

num_morphologies = 1
plot_flag = 1
flag_val = 0
weight_bonus = 0.05
exponent_val = 3
max_stems = 3

x_00 = 0
sig2_00 = 0

max_N = len(n[0])
max_R = len(n)

part1 = 1
part2 = 0

C = len(n)
T = len(n[0])

delta = 1
alpha = 1

ones_const = np.ones(T)
r_old = np.zeros(T)
r_new = 0.5*r_old

closeEnough = 1e-6

rho_new = 0.9
state_variance_new = 0.1
beta_new = 1
mu_new = -13

# rho_new = 0.3
# state_variance_new = 0.5
# beta_new = -2
# mu_new = -5

tolerance = 0.001
max_iter = 1
max_iter2 = 100

x_0_new = 0.1
var_0_new = state_variance_new/(1-rho_new**2)


part1_flag = 1
part2_flag = 0
part3_flag = 0

# file_name_pp = '../Data/Metrics/claiborne_point_process.txt'
file_name_pp = '..\\..\\Data\\Metrics\\amaral_point_process_basal.txt'

#file_name_pp = '/home/zzchou/Desktop/claiborne_point_process.txt'

#file_name_pp = '\Users\Zane\Desktop\Morphologies\Codes\claiborne_point_process.txt'


# Flags
rescale_flag = 1    # rescales all point processes to terminate at max R
unique_flag = 0   # removes recounted point processes
one_per_stem_flag = 0 # one point process per stem
median_per_stem_flag = 0 # one point process per stem with the median number of bifurcations (rounded up)
raster_flag = 1     # visualizes loaded point process data
EM_flag = 1         # executes expectation maximization algorithm
KS_flag = 1         # executes KS goodness of fit test
simulate_flag = 1   # simulates point processes
# Variables
max_EM_iterations = 3  # maximum number of iterations for EM algorithm

# Read point processes from file
n = read_file(file_name_pp)

max_N = len(n[0])
max_R = len(n)

if one_per_stem_flag == 1:
    current_indices = [0]
    all_indices = []
    new_n = []
    for r in range(0, max_R):
        list_of_indices = [i for i, x in enumerate(n[r]) if x == 1]
        if list_of_indices[0] != current_indices[0]:
            current_indices = list_of_indices
            new_n.append(n[r])
    n = new_n
max_R = len(n)

if median_per_stem_flag == 1:
    current_indices = [0]
    all_indices = []
    new_n = []
    for r in range(0, max_R):
        list_of_indices = [i for i, x in enumerate(n[r]) if x == 1]
        if list_of_indices[0] != current_indices[0]:
            all_indices = []
            current_indices = list_of_indices
            r2 = r
            list_of_indices2 = [current_indices[0]]
            while list_of_indices2[0] == current_indices[0] and r2 < max_R:
                list_of_indices2 = [i for i, x in enumerate(n[r2]) if x == 1]
                all_indices.append(list_of_indices2)
                r2 = r2+1
            len_of_indices = [len(x) for x in all_indices]
            median_len = int(np.median(len_of_indices))
            median_indices = [i for i, x in enumerate(len_of_indices) if x == median_len]
            ii = 0
            while median_indices == []:
                median_indices = [i for i, x in enumerate(len_of_indices) if x == (median_len + ii)]
                ii = ii + 1
            new_n.append(n[r + median_indices[0]])
    n = new_n
max_R = len(n)
    
if unique_flag == 1:
    current_indices = [0]
    for r in range(0, max_R):
        list_of_indices = [i for i, x in enumerate(n[r]) if x == 1]
        if list_of_indices[0] != current_indices[0]:
            current_indices = []
        
        new_point_process = [0]*max_R
        
        for j in range(len(list_of_indices)):
            index = list_of_indices[j]
            if index not in current_indices:
                new_point_process[index] = 1
                current_indices.append(index)
        
        n[r] = new_point_process

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
        


    
# VISUAL: Raster plot of rescaled point processes
if raster_flag == 1:
    plt.figure()
    to_be_rastered = []
    for r in range(0,max_R):
        event_index_list = [x for x, y in enumerate(n[r]) if y == 1]
    
        to_be_rastered.append(event_index_list)
    ax = raster(to_be_rastered)
    plt.xlabel('Rescaled Distance of Bifurcation Occurrence (um)', fontsize=18)
    plt.ylabel('Point Process #', fontsize=18)
    plt.show()

T = len(n[0])

x_smooth = []
alpha = 0.25
alpha1 = 0.5
alpha2 = 0.25
alpha3 = 0.5
alpha4 = 0.5

if part1_flag == 1:

    
    # pickle_in = open("convergence_filtered_keys.pickle","rb")
    # filtered_keys = pickle.load(pickle_in); keys = filtered_keys[0]; #print('Filtered Key:',keys)
    # # rho_new = keys[0] ; state_variance_new = keys[1] ; beta_new = keys[2] ; mu_new = keys[3]
    rho_new = 0.9
    state_variance_new = 0.1
    beta_new = 1
    mu_new = -13
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
                    # print('Slice:',mu,beta,x_old)
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
                u2_sum = u2_sum + np.exp((B[i]+delta2)*x_kIK[k] + (1/2)*((B[i]+delta2)**2)*var_kIK[k])*delta
                f_sum = f_sum + (np.exp(B[i]*x_kIK[k] + (1/2)*(B[i]**2)*var_kIK[k]))*(x_kIK[k] + B[i]*var_kIK[k])*delta
                f2_sum = f2_sum + (np.exp((B[i]+delta2)*x_kIK[k] + (1/2)*(B[i]+delta2)*var_kIK[k]))*(x_kIK[k] + (B[i]+delta2)*var_kIK[k])*delta
                n_sum = n_sum + n[c][k]*x_kIK[k]
                du_sum = du_sum + (np.exp(x_kIK[k] + B[i]*var_kIK[k]))*delta
                du_sum2 = du_sum2 + x_kIK[k] + B[i]*var_kIK[k]
                # print B[i], x_kIK[k], var_kIK[k]
            u = np.log(N_c) - np.log(u_sum)
            u2 = np.log(N_c) - np.log(u2_sum)
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
    
    # print rho_new, state_variance_new, beta_new, mu_new
    global p_max
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
    plt.plot(sv_filterx, sv_filter)
    plt.xlabel('sv_filterx')
    plt.ylabel('sv_filter')
    
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
                # print(lambda_max)
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
                poisson_flag = 1
    
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

# KOLMOGOROV SMIRNOV GOODNESS OF FIT TEST
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

global c_id
global max_prob
c_id = 1

view1 = 15
view2 = 45

#simulation_time = 430
# num_morphologies = 50
simulation_time = int(300 + random.random()*130)
stem_pathlength_coef = (0.49, 1.49, 1.8)

#inertial_weight = 2
#neighbor_weight = 0.25
#gradient_weight = 0.5
#random_weight = 1
#max_prob = 0
#
#growing_inertial_weight = 2
#growing_neighbor_weight = 1
#growing_gradient_weight = 1.5 + random.gauss(0, 0.1)
#growing_random_weight = 0.15

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



# MAIN

# 1 stem
    
#arg1 = 0.08
#arg2 = 0.01
#arg3 = 0.997
#arg4 = 0.9
#arg5 = 0.1
    
#arg1 = 0.08
#arg2 = 0.12
#arg3 = 0.06
#arg4 = 0.35
#arg5 = 0.05

arg1 = 0
arg2 = 0
arg3 = 0
arg4 = 0
arg5 = 0

#arg1 = 2.536
#arg2 = 3.747
#arg3 = 1.904
#arg4 = 62.336
#arg5 = 0.277


#arg1 = 1.566
#arg2 = 8.473
#arg3 = 0.846
#arg4 = 28.208
#arg5 = 0.185
# 1-2 stems

#arg1 = 0.068
#arg2 = 0.306
#arg3 = 0.2066
#arg4 = 0.884
#arg5 = 0.04

#simulation_time = 300
origin = [0, 0, 0]
stem_diameter = 1
stem_vector = unit_vector([0, 1, 0])
branch_list = []

theta_sd = 10
phi_sd = 180

if part3_flag == 1:
   
    for j in range(0, num_morphologies):
        flag = 1
        while flag == 1:
            c_id = 1
            simulation_time = -1
            while simulation_time > max_N or simulation_time < 0:
                simulation_time = int(250 + random.gauss(100,50))
            # simulation_time = 50
            # print 'simulation_time =',simulation_time    
                
            n_stem = int(random.random()*max_stems) + 1
            # n_stem = 1
            print('morphology #', j+1)
        
            if n_stem == 1:
                origin = [0, 0, 0]
                
                theta = random.uniform(-theta_sd, theta_sd)
                phi = random.uniform(-phi_sd, phi_sd)
                
                # theta = 0
                # phi = 0
                
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
                # rand_theta = random.gauss(0, theta_sd)
                # rand_phi = random.uniform(-phi_sd, phi_sd)     
                rand_phi = -90
                # bif_amp = draw_value('bif_amplitude_local')
                # if rand_theta == 0:
                #     theta_diff = bif_amp
                # else:
                #     theta_diff = 180*bif_amp/(2*rand_theta)
                # if rand_theta > 0:
                #     rand_theta2 = rand_theta - bif_amp/2
                #     rand_theta = rand_theta + bif_amp
                # else:
                #     rand_theta2 = rand_theta + bif_amp/2
                #     rand_theta = rand_theta - bif_amp
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
                branch_list = simulate_growth(branch_list, arg1, arg2, arg3, arg4, arg5, t, simulation_time-1)
            
            # for branch in branch_list:
            #     print branch.pathlength_from_soma - branch.pathlength
            # b_list = [b for b in branch_list if (b.terminal_flag == 0 or b.pathlength_from_soma - b.pathlength < simulation_time - 100)]
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
            # print sholl_error
            # print '# of branches', len(b_list)

            bif_num_min = 12
            bif_num_max = 45
            bif_num_mid = (bif_num_min + bif_num_max)/2
            bif_num_range = bif_num_max - bif_num_min
            
            bif_num = len(b_list)
            z = (bif_num - bif_num_mid)
            sigma = 10
            gauss_max = np.exp(-(0)**2/2)/(sigma*math.sqrt(2*math.pi))
            gaussian_filter = np.exp(-(-z/sigma)**2/2)/(sigma*math.sqrt(2*math.pi))/gauss_max
            
            # if len(b_list) < 0 or len(b_list) > 48:
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
        # plt.plot(sholl_bins, sholl_avgs)
        # plt.plot(sholl_bins, sholl_adjusted)    
        print('done')
        
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
        
    
        # rot_z = 0
        # print rot_x, rot_y, rot_z
    
        rotated_sum_vector, throw_away_axis = angle_noise([0,0,0], sum_vector, sum_vector, [1,0,0], -rot_x, rot_y, rot_z)
    
        # if plot_flag == 1:
        # #        plot_branches(b_list)
        #     plot_branches(branch_list)
        #     blah = zip([0,0,0],sum_vector)
        #     blah2 = zip([0,0,0], rotated_sum_vector)
        #     plt.plot(blah[0], blah[2], blah[1])
        #     plt.plot(blah2[0], blah2[2], blah2[1])
        
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
                    
                    # print jj, ref_coord[0], branch_list[jj].end_coordinates[0], ref_coord[2], branch_list[jj].end_coordinates[2]
    
        ref_coord2 = b_list[max_index].end_coordinates
        
        ref_difference = [ref_coord[kk] - ref_coord2[kk] for kk in range(0,3)]
        ref_difference[1] = 0
        # print ref_coord, ref_coord2
        final_rot_x, final_rot_y, final_rot_z = rotation_angles(ref_difference, [0,0,1])
        
        final_rot_z = final_rot_z * 180/math.pi
    
        for branch in b_list:
            for compartment in branch.compartment_list:
                if compartment.start_coordinates != [0,0,0]:
                    compartment.start_coordinates, throwaway_axis = angle_noise([0, 0, 0], compartment.start_coordinates, [0,1,0], [1,0,0], 0, -final_rot_z, 0 )
                compartment.end_coordinates, throwaway_axis = angle_noise([0, 0, 0], compartment.end_coordinates, [0,1,0], [1,0,0], 0, -final_rot_z, 0)
            branch.start_coordinates = branch.compartment_list[0].start_coordinates
            branch.end_coordinates = branch.compartment_list[-1].end_coordinates    
        #file_name = '/home/zzchou/Documents/Codes/claiborne_point_process.txt'
        # output_title = 'C:/Users/Zane/Desktop/Morphologies/Data/test_output/{}{}{}'.format('morphology_tester_output_', j+1, '.swc')
        output_title = '../Data/Morphologies/test_output/{}{}{}'.format('morphology_tester_output_', j+1, '.swc')
        soma_diameter = random.gauss(2.8, 0.3)
        with open(output_title, 'w') as f:
            true_origin = b_list[0].compartment_list[0]
            line = (true_origin.compartment_id, 1, true_origin.start_coordinates[0], true_origin.start_coordinates[1], true_origin.start_coordinates[2], soma_diameter, -1)
            parent_id = 1
            strs=' '.join(str(x) for x in line)
            f.write(strs + '\n')
    
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
                        new_line = (compartment.compartment_id, 3, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter, parent_id)
                        if new_line == line:
                           delete_flag = 1
                        line = (compartment.compartment_id, 3, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter, parent_id)
                        strs=' '.join(str(x) for x in line)
                        f.write(strs + '\n')
                    parent_id = compartment.compartment_id
        f.closed
    if  delete_flag == 1:
        os.remove(output_title)
    print('Success Rate:', num_morphologies/(num_morphologies + fail_count))
    if plot_flag == 1:
        plot_branches(b_list)
        # plot_branches(branch_list)



