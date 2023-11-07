#%% B_measure_tester.py
from __future__ import division
import os
import glob
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import random; np.random.seed(8675309); random.seed(8675309)
from sys import platform
import pickle

# Importing Custom Functions
from function_repo import vertical_align

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

    def set_length(self):
        self.length = euclidean_distance(self.start_coordinates, self.end_coordinates)

class Branch:
    def __init__(self, compartment_list):
        self.start_compartment = compartment_list[0]
        self.end_compartment = compartment_list[-1]

        self.start_coordinates = (self.start_compartment.start_coordinates)
        self.end_coordinates = (self.end_compartment.end_coordinates)

        self.vector = [self.end_coordinates[i] - self.start_coordinates[i] for i in range(0, 3)]

        self.compartment_list = compartment_list

        self.parent_bifurcation = []
        self.daughter_bifurcation = []
        self.rot_x = 0
        self.rot_y = 0
        self.rot_z = 0     

        self.fragmentation = len(compartment_list) - 1
        self.euclidean = euclidean_distance(self.start_coordinates, self.end_coordinates)
        self.pathlength = sum(compartment.length for compartment in compartment_list)
        self.euclidean_from_soma = euclidean_distance(self.start_coordinates, (0, 0, 0))
        self.pathlength_to_soma = self.pathlength  
        
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
                daughter.vector = [daughter.end_coordinates[i] - daughter.start_coordinates[i] for i in range(0, 3)]
                daughter.rot_x, daughter.rot_y, daughter.rot_z = rotation_angles(daughter.vector, [0,0,1])                
                # Recursion Error
                daughter.set_pathlength_to_soma() # Removed for testing

    def set_parent_bifurcation(self, parent_bifurcation):
        self.parent_bifurcation = parent_bifurcation
        parent_bifurcation.set_daughter_branch(self)
        # print(parent_bifurcation.daughter_branch) # For testing

    def set_daughter_bifurcation(self, daughter_bifurcation):
        self.daughter_bifurcation = daughter_bifurcation
        daughter_bifurcation.set_parent_branch(self)

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

class Morphology:
    def __init__(self, compartment_list, branch_list, bifurcation_list, terminal_list):
        global sholl_bins
        self.sholl_bins = sholl_bins
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

        
        
        self.max_x = 0; self.max_y = 0; self.max_z = 0
        for compartment in self.compartment_list:
            local_GF, GF_x_raw, GF_y_raw, GF_z_raw, GF_neighbor_raw = get_local_GF(compartment.end_coordinates, self.branch_list)
            self.GF_list.append(local_GF)
            self.GF_x_raw_list.append(GF_x_raw)
            self.GF_y_raw_list.append(GF_y_raw)
            self.GF_z_raw_list.append(GF_z_raw)
            self.GF_neighbor_raw_list.append(GF_neighbor_raw)
            if compartment.euclidean_from_soma > self.max_euclidean:
                self.max_euclidean = compartment.euclidean_from_soma
            if abs(compartment.x) > abs(self.max_x):
                self.max_x = compartment.x
            if abs(compartment.y) > abs(self.max_y):
                self.max_y = compartment.y
            if abs(compartment.z) > abs(self.max_z):
                self.max_z = compartment.z

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
        ynormed_bins = np.linspace(0, self.max_y, len(sholl_bins))
        self.sholl_counts = [0]*len(sholl_bins)
        self.ysholl_counts = [0]*len(sholl_bins)
        self.sholl_normed_counts = [0]*len(sholl_bins)
        self.ysholl_normed_counts = [0]*len(sholl_bins)
        # print self.total_surface_area
        for compartment in compartment_list:
            if compartment.parent_compartment != []:
                for i in range(0, len(sholl_bins)):
                    if compartment.euclidean_from_soma >= sholl_bins[i] and compartment.parent_compartment.euclidean_from_soma < sholl_bins[i]:
                        self.sholl_counts[i] = self.sholl_counts[i] + 1 # Primarily for Sholl analysis
                    if compartment.y >= sholl_bins[i] and compartment.parent_compartment.y < sholl_bins[i]:
                        self.ysholl_counts[i] = self.ysholl_counts[i] + 1
                    if compartment.euclidean_from_soma >= normed_bins[i] and compartment.parent_compartment.euclidean_from_soma < normed_bins[i]:
                        self.sholl_normed_counts[i] = self.sholl_normed_counts[i] + 1
                    if compartment.y >= ynormed_bins[i] and compartment.parent_compartment.y < ynormed_bins[i]:
                        self.ysholl_normed_counts[i] = self.ysholl_normed_counts[i] + 1
        self.sholl_counts = [x for x in self.sholl_normed_counts]         
        
        # bif_angle_list = [x.bif_amp_remote for x in self.bifurcation_list]
        # bin_edges = np.histogram(bif_angle_list, bins=50)[1] #get the bin edges

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

def euclidean_distance(u,v):
    ans = math.sqrt((u[0] - v[0])**2 + (u[1] - v[1])**2 + (u[2] - v[2])**2)
    return ans

def get_local_GF(coordinates, branch_list):
    x = coordinates[0]
    y = coordinates[1]
    z = coordinates[2]

    x_max = 500
    y_max = 500
    z_max = 500

    # consumption_rate = 1
    # distance_falloff = 10

    # GF_weight = 0.5


    # neighbor_weight = 300

    GF_scaling_y = 0.01
    GF_scaling_x = 0.3
    GF_scaling_z = 0.3
    GF_scaling_n = 0.3

    # GF_y_raw = (y_max - y)/y_max
    GF_y_raw = abs(y_max - y)/y_max
    GF_x_raw = abs(x)/x_max
    GF_z_raw = abs(z)/z_max
    GF_neighbor_raw = 0

    local_GF = GF_scaling_y*((y_max-y)/y_max) + GF_scaling_x*(abs(x)/x_max) + GF_scaling_z*(abs(z)/z_max)
    neighbor_max = 10


    # GF_gradient = (y_max - y)/y_max + abs(x)/x_max + abs(z)/z_max
    neighbor_gradient = 0

    for branch in branch_list:
        if branch.start_coordinates[1] < y and branch.end_coordinates[1] > y:
            i = 0
            neighbor= branch.compartment_list[0]
            while (branch.compartment_list[i].end_coordinates[1] < y) and (i < len(branch.compartment_list)):
                neighbor = branch.compartment_list[i]
                i = i + 1

            r = euclidean_distance(neighbor.end_coordinates, coordinates)
            if r != 0:
                # neighbor_gradient = neighbor_gradient + consumption_rate*(1/neighbor.diameter**2)*(distance_falloff/(r**2))
                neighbor_gradient = neighbor_gradient + 1/((1+r)**2)
                # neighbor_max = neighbor_max + 1

                # local_GF = local_GF - neighbor_weight/((r**2) * (neighbor.diameter**2))
                # GF_neighbor_raw = GF_neighbor_raw + 1/((r**2) * (neighbor.diameter**2))
                # neighbor_vector = [neighbor.end_coordinates[i] - coordinates[i] for i in range(0, 3)]

    neighbor_subtotal = 0
    if neighbor_max != 0:
        neighbor_subtotal = neighbor_gradient/neighbor_max
    GF_neighbor_raw = neighbor_subtotal
    local_GF = local_GF + GF_scaling_n*neighbor_subtotal
    local_GF = local_GF / (GF_scaling_x + GF_scaling_y + GF_scaling_z + GF_scaling_n)

    # print GF_gradient, neighbor_gradient
    # GF_score = GF_gradient*GF_weight + neighbor_gradient*neighbor_weight

    # print coordinates, GF_score
    # return GF_score
    return local_GF, GF_x_raw, GF_y_raw, GF_z_raw, GF_neighbor_raw

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
        proj1 = [vector1[0] - vector2[0], 0, vector1[2] - vector2[2]]
        proj2 = [0,0,1]
        bifurcation.bif_amp_vector = math.acos(dot_product(proj1, proj2)/magnitude_product(proj1, proj2))*180/math.pi

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
        # print compartment
        broken_flag = 0
        if compartment != []:
            if len(broken_branches) != 0:
                # print broken_branches
                # print(f'Broken Branch Detected: ({len(broken_branches)})')
                for jj in range(0, len(broken_branches)):
                    # print(f'cID: {compartment.compartment_id}\npID: {broken_branches[jj][-1].parent_id}')
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
                        # bif_match = [bif for bif in bifurcation_list if int(bif.bifurcation_id) == int(last_parent_id)]
                        
                        # print last_parent_id
                        # for bif in bifurcation_list:
                        #     print bif.bifurcation_id
                        # winner = bif_match[0]   
                        # if bif_match = []:
                            
                        # new_branch.set_parent_bifurcation(bif_match[0])
                        pass

                        
    
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

def plot_branches(branch_list, graph_title = 'Blank Title'):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

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
        plt.plot(x_list,  y_list, z_list,color = 'black')
        x_list = []
        y_list = []
        z_list = []
    # # Title Labeling 
    plt.title(graph_title)

    ax.set_xlim([-200,200])
    ax.set_ylim([-200,200])
    ax.set_zlim([0,400])
    # plt.axis('off')
    plt.grid(b=None)
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # ax.view_init(elev=0, azim=90)
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
    pdf_expon = stats.expon.pdf(lnspc, be, ce)
    # plt.plot(lnspc, pdf_expon, label = "Exp")
    m, s = stats.norm.fit(distribution)
    ag, bg, cg = stats.gamma.fit(distribution)
    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
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
    # print("Exponential Coef", be, ce)
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

def read_file(file_path, plt_morph_flag = False, title_addendum = ''):
    # Returns lists of all compartments, branches, bifurcations, and termination objects in the chosen swc file
    
    all_compartments = []
    compartment_dict = {}
    file_handle = open(file_path, 'r')
    lines = file_handle.read().splitlines()


    if platform == "linux" or platform == "linux2" or platform == 'darwin':
        splitChar = '/'
    elif platform == "win32":
        splitChar = '\\'
    name_split = file_path.split(f'{splitChar}'); file_name = name_split[-1]; print(file_name)

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

    # --------------- #
    # Vertically Transforming Compartments with respect to apical
    # WARNING: Only designed for pyramidal neurons for now
    apical_coords = np.empty(3) # Empty array for loop
    basal_coords = np.empty(3)
    for comp in all_compartments: # Getting Coordinates List
        if comp.section_type == 4 or comp.section_type == 1: 
            apical_coords = np.vstack((apical_coords,list(comp.start_coordinates)))
            apical_coords = np.vstack((apical_coords,list(comp.end_coordinates)))
        elif comp.section_type == 3:
            basal_coords = np.vstack((basal_coords,list(comp.start_coordinates)))
            basal_coords = np.vstack((basal_coords,list(comp.end_coordinates)))
        else:
            continue 
    apical_coords = np.delete(apical_coords,0,axis=0) # Removing the empty array
    basal_coords = np.delete(basal_coords,0,axis=0)
    # Transform
    apical_coords, basal_coords = vertical_align(apical_coords,basal_coords, plot_flag= False)
    # Reassigning
    for comp in all_compartments:
        if comp.section_type == 4 or comp.section_type == 1:
            comp.x,comp.y,comp.z = apical_coords[0][0:3]
            comp.start_coordinates = tuple(apical_coords[0])
            comp.end_coordinates = tuple(apical_coords[1])
            apical_coords = np.delete(apical_coords,[0,1], axis = 0) # Deletes most recent entry so same indexing can be used
        elif comp.section_type == 3:
            comp.x,comp.y,comp.z = basal_coords[0][0:3]#; comp.z = -comp.z
            comp.start_coordinates = tuple(basal_coords[0])
            comp.end_coordinates = tuple(basal_coords[1])
            basal_coords = np.delete(basal_coords,[0,1], axis = 0)
        else:
            continue
    # --------------- #


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

    if plt_morph_flag == True:
        plot_branches(all_branches, graph_title = ' '.join([file_name, title_addendum]))
    
    return all_compartments, all_branches, all_bifurcations, all_terminations

def morph_aggregate(file_list, pyramidal_flag = False):
    # Extracts the morphological features of cells in a dataset. Features extracted include the following:
    # compartment_list, branch_list, bifurcation_list, terminal_list, morphologies_list

    plot_read_flag = True # Plots each morphology that is read from the file list 

    if pyramidal_flag == False:
        # Pre-allocating for appending
        feature_aggregate = []
        all_morphologies_list = []
        all_compartments_list = []
        all_branches_list = []
        all_bifurcations_list = []
        all_terminations_list = []
        all_total_dendritic_length = []
        all_bifurcation_counts = []
        
        for i in range(len(file_list)):
            compartment_list, branch_list, bifurcation_list, termination_list = read_file(file_list[i],plt_morph_flag = plot_read_flag)
            
            current_morphology = [Morphology(compartment_list,branch_list,bifurcation_list,termination_list)]
            all_compartments_list = all_compartments_list + compartment_list
            all_branches_list = all_branches_list + branch_list
            all_bifurcations_list = all_bifurcations_list + bifurcation_list
            all_terminations_list = all_terminations_list + termination_list
            all_morphologies_list = all_morphologies_list + current_morphology


            total_dendritic_length = 0
            for branch in branch_list:
                total_dendritic_length = total_dendritic_length + branch.pathlength
            all_total_dendritic_length.append(total_dendritic_length)

            all_bifurcation_counts.append(len(bifurcation_list))

        # Packaging and returning
        feature_aggregate  = [all_compartments_list,all_branches_list,all_bifurcations_list,all_terminations_list,all_morphologies_list]
        return feature_aggregate
    elif pyramidal_flag == True:
        # Pre-allocating for appending

        # Apical
        feature_aggregate_apical = []
        all_apical_arbor_list = []
        all_compartments_list_apical = []
        all_branches_list_apical = []
        all_bifurcations_list_apical = []
        all_terminations_list_apical = []
        all_total_dendritic_length_apical = []
        all_bifurcation_counts_apical = []

        # Basal
        feature_aggregate_basal = []
        all_basal_arbor_list = []
        all_compartments_list_basal = []
        all_branches_list_basal = []
        all_bifurcations_list_basal = []
        all_terminations_list_basal = []
        all_total_dendritic_length_basal = []
        all_bifurcation_counts_basal = []
        
        for i in range(len(file_list)):
            try:
                compartment_list, branch_list, bifurcation_list, termination_list = read_file(file_list[i],plt_morph_flag = plot_read_flag)

                # Section: Apical
                compartment_list_apical = [comp for comp in compartment_list if comp.section_type == 4 or comp.section_type == 1]
                branch_list_apical = [branch for branch in branch_list if branch.compartment_list[-1].section_type == 4 or branch.compartment_list[-1].section_type == 1]
                bifurcation_list_apical = [bif for bif in bifurcation_list if bif.bifurcating_compartment.section_type == 4 or bif.bifurcating_compartment.section_type == 1]
                termination_list_apical = [comp for comp in compartment_list_apical if len(comp.daughter_compartment) == 0]

                total_dendritic_length_apical = 0

                for branch in branch_list_apical:
                    total_dendritic_length_apical = total_dendritic_length_apical + branch.pathlength
                all_total_dendritic_length_apical.append(total_dendritic_length_apical)

                all_bifurcation_counts_apical.append(len(bifurcation_list_apical))

                global sholl_bins
                # # Appending to aggregate
                all_compartments_list_apical = all_compartments_list_apical + compartment_list_apical
                # Dataset-Dependent Sholl Bin Setup
                comp_euclideans = [comp.euclidean_from_soma for comp in all_compartments_list_apical]; comps_max_euclidean = max(comp_euclideans)
                sholl_bins = np.linspace(0, comps_max_euclidean, 10)
                all_branches_list_apical = all_branches_list_apical + branch_list_apical
                all_bifurcations_list_apical = all_bifurcations_list_apical + bifurcation_list_apical
                all_terminations_list_apical = all_terminations_list_apical + termination_list_apical
                current_apical = [Morphology(compartment_list_apical,branch_list_apical,bifurcation_list_apical,termination_list_apical)]
                all_apical_arbor_list = all_apical_arbor_list + current_apical

                # ----------------------------------------------------------------- #

                # Section: Basal
                compartment_list_basal = [comp for comp in compartment_list if comp.section_type == 3 or comp.section_type == 1]
                branch_list_basal = [branch for branch in branch_list if branch.compartment_list[-1].section_type == 3 or branch.compartment_list[-1].section_type == 1]
                bifurcation_list_basal = [bif for bif in bifurcation_list if bif.bifurcating_compartment.section_type == 3 or bif.bifurcating_compartment.section_type == 1]
                termination_list_basal = [comp for comp in compartment_list_basal if len(comp.daughter_compartment) == 0]

                total_dendritic_length_basal = 0

                for branch in branch_list_basal:
                    total_dendritic_length_basal = total_dendritic_length_basal + branch.pathlength
                all_total_dendritic_length_basal.append(total_dendritic_length_basal)

                all_bifurcation_counts_basal.append(len(bifurcation_list_basal))

                
                
                all_compartments_list_basal = all_compartments_list_basal + compartment_list_basal
                # Dataset-Dependent Sholl Bin Setup
                comp_euclideans = [comp.euclidean_from_soma for comp in all_compartments_list_basal]; comps_max_euclidean = max(comp_euclideans)
                sholl_bins = np.linspace(0, comps_max_euclidean, 10)
                all_branches_list_basal = all_branches_list_basal + branch_list_basal
                all_bifurcations_list_basal = all_bifurcations_list_basal + bifurcation_list_basal
                all_terminations_list_basal = all_terminations_list_basal + termination_list_basal
                current_basal = [Morphology(compartment_list_basal,branch_list_basal,bifurcation_list_basal,termination_list_basal)]
                all_basal_arbor_list = all_basal_arbor_list + current_basal


                # ------------------------------------------------------------------ #
                plot_flag = False
                if plot_flag == True:
                    title_addendum = ''
                    plot_branches(branch_list_apical, graph_title = 'apical ')
                    plot_branches(branch_list_basal, graph_title = 'basal')
            except AttributeError:
                print(f'Attribute ERROR: {file_list[i]}')
            except ValueError:
                print(f'Value ERROR: {file_list[i]}')
            # Packing and Returning
            feature_aggregate_apical = [all_compartments_list_apical,all_branches_list_apical,all_bifurcations_list_apical,all_terminations_list_apical,all_apical_arbor_list]
            feature_aggregate_basal = [all_compartments_list_basal,all_branches_list_basal,all_bifurcations_list_basal,all_terminations_list_basal,all_basal_arbor_list]
            

        return feature_aggregate_apical, feature_aggregate_basal
    else:
        print('Error: Invalid Pyramidal Flag')

def distribution_gen(cumulative_list):
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

    distribution_list.append('Interbifurcation_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([branch.pathlength for branch in branch_list])

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
    
    distribution_list.append('max_y')
    distribution_type_list.append('basic')
    entries_list.append([morphology.max_y for morphology in morphologies_list])

    distribution_list.append('max_z')
    distribution_type_list.append('basic')
    entries_list.append([morphology.max_z for morphology in morphologies_list])

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
    entries_list.append(morphologies_list[0].sholl_bins) # Taken from the first morph since they should all be the same

    distribution_list.append('sholl_counts')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.sholl_normed_counts for morphology in morphologies_list])

    distribution_list.append('max_euclidean')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.max_euclidean for morphology in morphologies_list])

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

    # tdl_bins = np.linspace(0, max_max_euclidean-1, max_max_euclidean)

    # tmp_list = [morphology.tdl_distribution for morphology in morphologies_list]
    # tdl = [sum(x)/len(morphologies_list) for x in zip(*tmp_list)]
    # plt.plot(tdl_bins, tdl, label='Generated TDL')    
    # plt.title('TDL for Generated')
    # plt.legend()
    # plt.show()
    # print(f'{len(tdl_bins)} | {len(tdl)}')
    
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
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # cax = ax.matshow(correlations, vmin=0, vmax=1, cmap='viridis')
    # fig.colorbar(cax)
    # ticks = np.arange(0,len(names),1)
    # ax.set_xticks(ticks)
    # ax.set_yticks(ticks)
    # ax.set_xticklabels(names, rotation=90)
    # ax.set_yticklabels(names)
    # ax.set_xticklabels([x+1 for x in ticks])
    # ax.set_yticklabels([x+1 for x in ticks])  
    
    
    # fig = plt.figure()
    # plt.scatter(branch_distance_distribution, pathlength_distribution)
    # plt.title('pathlength')
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
            # plt.imshow(H, interpolation='nearest', aspect='auto', origin='low', extent=[x_edge_min, x_edge_max, y_edge_min, y_edge_max], cmap='viridis', vmin=vmin_val, vmax=H_max)


            # plt.colorbar()
            # plt.title(names[i] + ' vs ' + names[j])


    # Packing & Returning
    dist_pack = [distribution_list, distribution_type_list, entries_list]
    return dist_pack

def save_distributions(dist_file, dist_pack):
    # Saves distributions into a separate text file for later use
    dist_list = dist_pack[0]
    dist_type_list = dist_pack[1]
    entries_list = dist_pack[2]
    
    # Saving to pickle to conserve data types
    pickle_file = dist_file.split('.'); pickle_file = pickle_file[0] + '.pickle' # Changing file extenstion
    dist_dict = {}
    for i in range(len(dist_list)):
        dict_entry = {dist_list[i] : [dist_type_list[i], entries_list[i]]}
        dist_dict.update(dict_entry)
    with open(pickle_file,'wb') as f:
        pickle.dump(dist_dict,f)

    # Saving to text
    with open(dist_file, 'w') as f:
        for i in range(len(dist_list)):
            entries = ' '.join(str(x) for x in entries_list[i])
            line = f'{dist_list[i]};{dist_type_list[i]};{entries}\n'
            f.write(line)
    print('Distribution Saved')

def analyze_dataset(dataset_identifier, synthetic_flag, save_flag, pyramidal_flag):
    
    if dataset_identifier == 'claiborne':
        source_files = glob.glob(os.path.abspath(f'../../Data/Morphologies/{dataset_identifier}/*.swc'))
        dist_file = os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_.txt')

        morph_cumulative = morph_aggregate(source_files, pyramidal_flag = pyramidal_flag)
        dist_pack = distribution_gen(morph_cumulative)
        if save_flag == True:
            save_distributions(dist_file,dist_pack)
   
    elif dataset_identifier == 'amaral': 
        pyramidal_flag = True
        if synthetic_flag == False:
            source_files = glob.glob(os.path.abspath(f'../../Data/Morphologies/{dataset_identifier}/CA3/*.swc'))
            dist_file_apical = os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_apical.txt')
            dist_file_basal = os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_basal.txt')
        elif synthetic_flag == True:
            source_files = glob.glob(os.path.abspath(f'../../Data/Morphologies/Saved-Output/Full/*.swc'))
            dist_file_apical = os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_apical_SYNTHETIC.txt')
            dist_file_basal = os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_basal_SYNTHETIC.txt')
        else:
            print('Error: Invalid Synthetic Flag')

        morph_cumulative_apical, morph_cumulative_basal = morph_aggregate(source_files, pyramidal_flag = pyramidal_flag)

        dist_pack_apical = distribution_gen(morph_cumulative_apical)
        if save_flag == True:
            save_distributions(dist_file_apical,dist_pack_apical)

        dist_pack_basal = distribution_gen(morph_cumulative_basal)
        if save_flag == True:
            save_distributions(dist_file_basal,dist_pack_basal)
    
    else:
        print('Error: Invalid Dataset Identifier')




#%%
# # B_measure_tester.py
if __name__ == '__main__':
    pyramidal_flag = False # Flag set to False by default for any dataset unless specified otherwise (Do not change)
    
    dataset_identifier = 'amaral' # str to designate the file in which the database is located. Also designates which procedure will be followed 
    synthetic_flag = False # Designates switching between organic and synthetic datasets 
    save_flag = True # Whether to save the distributions
    
    analyze_dataset(dataset_identifier, synthetic_flag, save_flag, pyramidal_flag)

    # compare_datasets(dataset_identifier,'apical')










    


