from __future__ import division

import glob
import os
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from angle_noise_updated import angle_noise
import random

class Compartment:
    def __init__(self, line):
        self.compartment_id, self.section_type, self.x, self.y, self.z, self.diameter, self.parent_id = map(float, line)
        self.diameter = 2*self.diameter
        self.parent_compartment = []
        self.daughter_compartment = []
        self.start_coordinates = (0, 0, 0)
        
        x_sd = 0
        z_sd = 0
        
        self.jitter_x = 0
        self.jitter_z = 0
        
        self.end_coordinates = (self.x + self.jitter_x, self.y, self.z + self.jitter_z)

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

#         last_compartment = self.start_coordinates
        for compartment in compartment_list:
            if compartment != compartment_list[0]:
                compartment.branch  = self
#             self.pathlength = self.pathlength + math.sqrt((compartment.x-last_compartment[0])**2 + (compartment.y-last_compartment[1])**2 + (compartment.z-last_compartment[2])**2)
#             last_compartment = (compartment.x, compartment.y, compartment.z)
        if self.pathlength == 0:
            self.contraction = 0
        else:
            self.contraction = self.euclidean/self.pathlength

    def set_pathlength_to_soma(self):
        if self.daughter_bifurcation != []:
            for daughter in self.daughter_bifurcation.daughter_branch:
                daughter.pathlength_to_soma = daughter.pathlength + self.pathlength_to_soma
                daughter.vector = [daughter.end_coordinates[i] - daughter.start_coordinates[i] for i in range(0, 3)]
#                print daughter.vector
                daughter.rot_x, daughter.rot_y, daughter.rot_z = rotation_angles(daughter.vector, [0,1,0])
#                print daughter.rot_x, daughter.rot_y

                daughter.set_pathlength_to_soma()

    def set_parent_bifurcation(self, parent_bifurcation):
        self.parent_bifurcation = parent_bifurcation
        parent_bifurcation.set_daughter_branch(self)

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

sholl_bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

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
#        self.major_coordinates = [compartment.x for compartment in terminal_list]
#        self.minor_coordinates = [compartment.z for compartment in terminal_list]
#        
#        self.tmp_major_axis = abs(max(self.major_coordinates) - min(self.major_coordinates))
#        self.tmp_minor_axis = abs(max(self.minor_coordinates) - min(self.minor_coordinates))
#        
        original_coordinates = [np.array((compartment.x, compartment.z)) for compartment in terminal_list]
        
        self.major_axis = 0
        self.minor_axis = 0
        
        self.GF_list = []
        self.GF_x_raw_list = []
        self.GF_y_raw_list = []
        self.GF_z_raw_list = []
        self.GF_neighbor_raw_list = []
        self.branch_flag_list = []

        for compartment in self.compartment_list:
            local_GF, GF_x_raw, GF_y_raw, GF_z_raw, GF_neighbor_raw = get_local_GF(compartment.end_coordinates, self.branch_list)
            self.GF_list.append(local_GF)
            self.GF_x_raw_list.append(GF_x_raw)
            self.GF_y_raw_list.append(GF_y_raw)
            self.GF_z_raw_list.append(GF_z_raw)
            self.GF_neighbor_raw_list.append(GF_neighbor_raw)



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
        
        for branch in branch_list:
            self.total_dendritic_length = self.total_dendritic_length + branch.pathlength

        self.length_scale_factor = self.total_dendritic_length/100    
        for compartment in compartment_list:
            self.total_surface_area = self.total_surface_area + math.pi*compartment.diameter*compartment.length
            self.normalized_total_surface_area = self.normalized_total_surface_area + math.pi*compartment.diameter*compartment.length/self.length_scale_factor
        for termination in terminal_list:
            self.total_surface_area = self.total_surface_area + math.pi*(0.5*termination.diameter)**2
        self.num_bifurcations = len(bifurcation_list)
        
        self.sholl_counts = [0]*len(sholl_bins)
        
        for compartment in compartment_list:
            if compartment.parent_compartment != []:
                for i in range(0, len(sholl_bins)):
                    if compartment.euclidean_from_soma >= sholl_bins[i] and compartment.parent_compartment.euclidean_from_soma < sholl_bins[i]:
                        self.sholl_counts[i] = self.sholl_counts[i] + 1
                        
#        plt.scatter(sholl_bins, self.sholl_counts)
#        plt.scatter(self.major_coordinates, self.minor_coordinates)
#        plt.show()
        

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


def read_file(file_name):

    all_compartments = []
    compartment_dict = {}

    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()

#     Initialize last line check
    last_line = (-1, -1, -1)
    id_offset = 0
    for i in range(0, len(lines)):
        parsed_line = lines[i].split(' ')
        line = [entry for entry in parsed_line if entry != '']  #ignores empty strings in the case of double spaces
        

#         Ignore last line if it is a duplicate of the second to last line
        if len(line) > 0:
            if (line[0] != '#') and line[0][0] != '#':
                if last_line != (line[2], line[3], line[4]):
    #                 if line[1] == '3' or line[-1] == '-1':                    
                    new_compartment = Compartment(line)
                    if new_compartment.section_type == 3 or new_compartment.section_type == 1:
                        if new_compartment.compartment_id != 1:
                            new_compartment.compartment_id = new_compartment.compartment_id - id_offset
                        if new_compartment.parent_id != 1:
                            new_compartment.parent_id = new_compartment.parent_id - id_offset

                        all_compartments.append(new_compartment)
                        compartment_dict[new_compartment.compartment_id] = new_compartment
                        if new_compartment.parent_id >= 1:
                            try:
                                start_compartment = compartment_dict[new_compartment.parent_id]
                            except:
                                continue
                            new_compartment.parent_compartment = start_compartment
                            new_compartment.start_coordinates = start_compartment.end_coordinates
        
                            start_compartment.daughter_compartment.append(new_compartment)
                    else:
                        id_offset = id_offset + 1
                    last_line = (line[2], line[3], line[4])
                
    for i in range(0, len(all_compartments)):
        all_compartments[i].set_length()
    breakpoints = []    
    all_bifurcations = [Bifurcation(compartment) for compartment in all_compartments if len(compartment.daughter_compartment) > 1 and compartment.compartment_id != 1.0]
    for ii in range(len(all_compartments) - 2, -1, -1):
        if all_compartments[ii +1].parent_id != all_compartments[ii].compartment_id:
            try:
                aa = all_compartments[ii + 1].parent_compartment.daughter_compartment
            except:
                return [], [], [], []
            if len(all_compartments[ii + 1].parent_compartment.daughter_compartment) == 1:
                breakpoints.append(all_compartments[ii + 1])
#                print breakpoints
    all_terminations = [compartment for compartment in all_compartments if len(compartment.daughter_compartment) == 0]

    all_branches = divide_into_branches(all_compartments, all_bifurcations, breakpoints)
    
#    for branch in all_branches:
#        print branch.branch_order, branch.pathlength, branch.pathlength_to_soma

    calculate_remote_angles(all_bifurcations)
    calculate_local_angles(all_bifurcations)

    return all_compartments, all_branches, all_bifurcations, all_terminations

def cross_product((u1, u2, u3),(v1, v2, v3)):
    s1 = u2*v3 - u3*v2
    s2 = u3*v1 - u1*v3
    s3 = u1*v2 - u2*v1
    return (s1, s2, s3)


def dot_product((u1, u2, u3),(v1, v2, v3)):
    s = u1*v1 + u2*v2 + u3*v3
    return s


def magnitude_product((u1, u2, u3),(v1, v2, v3)):
    s1 = math.sqrt(u1**2 + u2**2 + u3**2)
    s2 = math.sqrt(v1**2 + v2**2 + v3**2)
    return s1*s2


def euclidean_distance((u1, u2, u3), (v1, v2, v3)):
    ans = math.sqrt((u1 - v1)**2 + (u2 - v2)**2 + (u3 - v3)**2)
    return ans

def get_local_GF(coordinates, branch_list):
    x = coordinates[0]
    y = coordinates[1]
    z = coordinates[2]

    x_max = 500
    y_max = 500
    z_max = 500

    consumption_rate = 1
    distance_falloff = 10

    GF_weight = 0.5


    neighbor_weight = 300

    GF_scaling_y = 0.01
    GF_scaling_x = 0.3
    GF_scaling_z = 0.3
    GF_scaling_n = 0.3

#    GF_y_raw = (y_max - y)/y_max
    GF_y_raw = abs(y_max - y)/y_max
    GF_x_raw = abs(x)/x_max
    GF_z_raw = abs(z)/z_max
    GF_neighbor_raw = 0

    local_GF = GF_scaling_y*((y_max-y)/y_max) + GF_scaling_x*(abs(x)/x_max) + GF_scaling_z*(abs(z)/z_max)
    neighbor_max = 10


    GF_gradient = (y_max - y)/y_max + abs(x)/x_max + abs(z)/z_max
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
                neighbor_vector = [neighbor.end_coordinates[i] - coordinates[i] for i in range(0, 3)]

    neighbor_subtotal = 0
    if neighbor_max != 0:
        neighbor_subtotal = neighbor_gradient/neighbor_max
    GF_neighbor_raw = neighbor_subtotal
    local_GF = local_GF + GF_scaling_n*neighbor_subtotal
    local_GF = local_GF / (GF_scaling_x + GF_scaling_y + GF_scaling_z + GF_scaling_n)

    # print GF_gradient, neighbor_gradient
    GF_score = GF_gradient*GF_weight + neighbor_gradient*neighbor_weight

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
        else:
            rot_x1, rot_y1, rot_z1 = rotation_angles(bifurcation.bif_midline, bifurcation.bif_normal)
            rot_x2, rot_y2, rot_z2 = rotation_angles(bifurcation.parent_branch.vector, [0,0,1])            
        
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

#            bifurcation.bif_tilt_remote = bifurcation.bif_tilt_remote * 180/math.pi
#            bifurcation.bif_torque_remote = bifurcation.bif_torque_remote * 180/math.pi
#            bifurcation.bif_twist_remote = bifurcation.bif_twist_remote * 180/math.pi
        
        bifurcation.bif_tilt_remote = rot_x_diff * 180/math.pi
        bifurcation.bif_torque_remote = rot_y_diff * 180/math.pi
        bifurcation.bif_twist_remote = rot_z_diff * 180/math.pi
            
    return

def rotation_angles(u, norm):
#    n1 = (1, 0, 0)
#    n2 = (0, 1, 0)
#
#    w1 = [dot_product(u,n1)/magnitude_product(n1,n1), 0, 0]
#    w2 = [0, dot_product(u,n2)/magnitude_product(n2,n2), 0]
#
#    proj_yz = [u[i] - w1[i] for i in range(0,3)]
#    proj_xz = [u[i] - w2[i] for i in range(0,3)]

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

#    R_total = Rr_x*Rr_z*Rr_y
#    
#    u_inv = np.ndarray.tolist(np.dot(u, R_total))[0]
#    norm_inv = np.ndarray.tolist(np.dot(norm, R_total))[0]
#    u_plot = zip((0,0,0), u)
#    norm_plot = zip((0,0,0), norm_inv)
#    
#    print u_inv
#    
#    ax = plt.axes(projection='3d')
#    ax.set_xlim([-50,50])
#    ax.set_ylim([-50,50])
#    ax.set_zlim([0,100])
#    
#    print u_plot
#    plt.plot(u_plot[0], u_plot[1], u_plot[2])
#    plt.plot(norm_plot[0], norm_plot[1], norm_plot[2])
    
#    plt.show()
#
#    print rot_x

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

#            bifurcation.bif_tilt_local = bifurcation.bif_tilt_local * 180/math.pi
#            bifurcation.bif_torque_local = bifurcation.bif_torque_local * 180/math.pi
#            bifurcation.bif_twist_local = bifurcation.bif_twist_local * 180/math.pi

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
#    print breakpoints
    for i in range(0,len(compartment_list)):
        compartment = compartment_list[len(compartment_list)-i-1]
#        print compartment
        broken_flag = 0
        if compartment != []:
            if len(broken_branches) != 0:
#                print broken_branches
                for jj in range(0, len(broken_branches)):
#                    print compartment.compartment_id, broken_branches[jj][-1].parent_id
                    if compartment.compartment_id == broken_branches[jj][-1].parent_id:
                        print 'yes'
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
#                    print compartment.compartment_id
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
    #                    print last_parent_id
    #                    for bif in bifurcation_list:
    #                        print bif.bifurcation_id
                        winner = bif_match[0]                    
#                        if bif_match = []:
                            
                        new_branch.set_parent_bifurcation(bif_match[0])
                        
    
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

def plot_branches(branch_list):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for branch in branch_list:
#        if branch.end_compartment.section_type != 1:
        if branch.parent_bifurcation == []:
            for compartment in branch.compartment_list:
                x_list.append(compartment.end_coordinates[0])
                y_list.append(compartment.end_coordinates[1])
                z_list.append(compartment.end_coordinates[2])
        elif branch.branch_order <= 2:
            for compartment in branch.compartment_list:
                x_list.append(compartment.end_coordinates[0])
                y_list.append(compartment.end_coordinates[1])
                z_list.append(compartment.end_coordinates[2])                
                
        plt.plot(x_list, z_list, y_list, color = 'black')
        x_list = []
        y_list = []
        z_list = []
#    ax.set_xlim([-200,200])
#    ax.set_ylim([-200,200])
#    ax.set_zlim([0,400])
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.view_init(elev=20, azim=130)
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
#    bins = np.linspace(-180, 180, 36)
    plt.hist(distribution, normed=True, bins=50)
#    print distribution
#    xt = plt.xticks()[0]
#    xmin, xmax = min(xt), max(xt)
#    lnspc = np.linspace(xmin, xmax, len(distribution))
#    be, ce = stats.expon.fit(distribution)
#    pdf_expon = stats.expon.pdf(lnspc, be, ce)
#    plt.plot(lnspc, pdf_expon, label = "Exp")
#    m, s = stats.norm.fit(distribution)
#    ag, bg, cg = stats.gamma.fit(distribution)
#    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
#    pdf_norm = stats.norm.pdf(lnspc, m, s)
#    plt.plot(lnspc, pdf_gamma, label="Gamma")
#    plt.plot(lnspc, pdf_norm, label="Norm")
#    
#    print("Exponential Coef", be, ce)
#    print("Normal Coef", m, s)
#    print("Gamma Coef", ag, bg, cg)
    
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
#    plt.plot(lnspc, pdf_gamma, label="Gamma")
    plt.plot(lnspc, pdf_norm, label="Norm")
    
    print("Exponential Coef", be, ce)
    print("Normal Coef", m, s)
    print("Gamma Coef", ag, bg, cg)

#    ks_check(distribution, 'stem_diameter')

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

#    ks_check(distribution, 'pathlength')

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
    

#    ks_check(distribution, 'euclidean_from_soma')
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
#    lnspc = np.linspace(xmin, xmax, len(distribution))
#    be, ce = stats.expon.fit(distribution)
#    pdf_expon = stats.expon.pdf(lnspc, be, ce)
#    plt.plot(lnspc, pdf_expon, label = "Exp")
#    m, s = stats.norm.fit(distribution)
#    ag, bg, cg = stats.gamma.fit(distribution)
#    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
#    pdf_norm = stats.norm.pdf(lnspc, m, s)
#    plt.plot(lnspc, pdf_gamma, label="Gamma")
#    plt.plot(lnspc, pdf_norm, label="Norm")
#    
#    print("Exponential Coef", be, ce)
#    print("Normal Coef", m, s)
#    print("Gamma Coef", ag, bg, cg)

    plt.show()

#    ks_check(distribution, 'total_dendritic_length')

def plot_pairs(pair_list):
    fig2  = plt.figure()

    distribution1 = []
    distribution2 = []

    for pair in pair_list:
        item1 = pair[0].bif_amp_vector
        item2 = pair[0].bif_amp_remote
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

#def ks_check(distribution, request):
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
#
#        'branch_order_0': (branch_order_0_constants[0], branch_order_0_constants[1], 0,'exponential'),
#        'branch_order_1': (branch_order_1_constants[0], branch_order_1_constants[1], branch_order_1_constants[2],'gamma'),
#        'branch_order_2': (branch_order_2_constants[0], branch_order_2_constants[1], 0,'exponential'),
#        'branch_order_3': (branch_order_3_constants[0], branch_order_3_constants[1], branch_order_3_constants[2],'gamma'),
#        'branch_order_4': (branch_order_4_constants[0], branch_order_4_constants[1], branch_order_4_constants[2],'gamma'),
#        'branch_order_5': (branch_order_5_constants[0], branch_order_5_constants[1], branch_order_5_constants[2],'gamma'),
#        'branch_order_6': (branch_order_6_constants[0], branch_order_6_constants[1], branch_order_6_constants[2],'gamma')    
#    }
#
#    distribution_type = parameters[request][-1]
#    if distribution_type == 'gamma':
#        fit_params = (parameters[request][0], parameters[request][1], parameters[request][2])
#    else:
#        fit_params = (parameters[request][0], parameters[request][1])
#
#
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

def reconstruct_morphology(compartment_list, branch_list, bifurcation_list, terminal_list):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for j in range(len(branch_list)):
        branch = branch_list[len(branch_list) - j - 1]
        
        if branch.parent_bifurcation == []:        
            coordinate_list = [branch.start_compartment, branch.end_compartment]
            for coord in coordinate_list:
                x_list.append(coord.end_coordinates[0])
                y_list.append(coord.end_coordinates[1])
                z_list.append(coord.end_coordinates[2])
            
    
            plt.plot(x_list, z_list, y_list, color = 'black')
            x_list = []
            y_list = []
            z_list = []
            
            x_list2 = []
            y_list2 = []
            z_list2 = []
        elif branch.branch_order <= 10:
            bif = branch.parent_bifurcation
#            x_rot = bif.bif_tilt_remote
#            y_rot = bif.bif_torque_remote
#            z_rot = bif.bif_twist_remote
             
            x_rot = branch.rot_x*180/math.pi
            y_rot = branch.rot_y*180/math.pi
            z_rot = branch.rot_z*180/math.pi
#            x_rot = 0
#            y_rot = 0
#            z_rot = 0
#            bif_amp = bif.bif_amp_remote
            
#            origin_array = np.asarray(branch.start_compartment.end_coordinates)
#            end_array = np.asarray(branch.end_compartment.end_coordinates)
#            
#            vector = end_array - origin_array
#            mag = np.sqrt(vector.dot(vector))
            parent_vector = bif.parent_branch.vector
            parent_array = np.asarray(parent_vector)
            parent_mag = np.sqrt(parent_array.dot(parent_array))
            extend_coord = [0]*3
#            print branch.vector, branch.rot_x, branch.rot_y, branch.rot_z
            for i in range(0,3):
                extend_coord[i] = branch.start_compartment.end_coordinates[i] + branch.pathlength*parent_vector[i]/parent_mag
            new_coord = angle_noise(branch.start_compartment.end_coordinates, extend_coord, x_rot, y_rot, z_rot)
            coordinate_list = [branch.start_compartment.end_coordinates, new_coord]
            print x_rot, y_rot, z_rot
            for coord in coordinate_list:
                x_list.append(coord[0])
                y_list.append(coord[1])
                z_list.append(coord[2])
            
            plt.plot(x_list, z_list, y_list, color = 'blue')

            for compartment in branch.compartment_list:
                x_list2.append(compartment.end_coordinates[0])
                y_list2.append(compartment.end_coordinates[1])
                z_list2.append(compartment.end_coordinates[2])
            plt.plot(x_list2, z_list2, y_list2, color = 'black')
            x_list = []
            y_list = []
            z_list = []

            x_list2 = []
            y_list2 = []
            z_list2 = []

    ax.set_xlim([-200,200])
    ax.set_ylim([-200,200])
    ax.set_zlim([0,400])
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.view_init(elev=20, azim=130)
    plt.show()

        
        
    
def save_distributions(compartment_list, branch_list, bifurcation_list, terminal_list, morphologies_list, distribution_list, distribution_type_list, entries_list):
    
    stem_list = [stem for stem in branch_list if stem.branch_order == 0]
    ibf_list = [ibf.parent_branch for ibf in bifurcation_list if ibf.parent_branch.branch_order != 0]
    terminal_branches = [branch for branch in branch_list if branch.daughter_bifurcation == []]
    
    distribution_list.append('Stem_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([stem.pathlength for stem in stem_list])   
    
    distribution_list.append('All_pathlength_GCL_(um)')
    distribution_type_list.append('basic')
    entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 50])
    
    distribution_list.append('All_pathlength_IML_(um)')
    distribution_type_list.append('basic')
    entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 150 and ibf.euclidean_from_soma > 50])
    
    distribution_list.append('All_pathlength_MML_(um)')
    distribution_type_list.append('basic')
    entries_list.append([ibf.pathlength for ibf in branch_list if (ibf.euclidean_from_soma > 150 and ibf.euclidean_from_soma <= 250)])

    distribution_list.append('All_pathlength_OML_(um)')
    distribution_type_list.append('basic')
    entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma > 250])    
    
    distribution_list.append('Terminal_Pathlength_GCL_(um)')
    distribution_type_list.append('basic')
    entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma <= 50])

    distribution_list.append('Terminal_Pathlength_IML_(um)')
    distribution_type_list.append('basic')
    entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 50 and termination.branch.euclidean_from_soma <= 150])
    
    distribution_list.append('Terminal_pathlength_MML_(um)')
    distribution_type_list.append('basic')
    entries_list.append([termination.branch.pathlength for termination in terminal_list if (termination.branch.euclidean_from_soma > 150 and termination.branch.euclidean_from_soma <= 250)])

    distribution_list.append('Terminal_pathlength_OML_(um)')
    distribution_type_list.append('basic')
    entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 250])  

    distribution_list.append('All_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([branch.pathlength for branch in branch_list])
    
    

#    distribution_list.append('Stem_Diameter_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([stem.start_diameter for stem in stem_list])
#
#    distribution_list.append('IBF_Diameter_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([ibf.start_diameter for ibf in ibf_list])
#
#    distribution_list.append('TerminaXl_Diameter_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([termination.branch.start_diameter for termination in terminal_list])
#
#    distribution_list.append('Compartment_Diameter_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([compartment.diameter for compartment in compartment_list if compartment.section_type == 3])
    
#    distribution_list.append('Post-bifurcation_Diameter_1_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([max(bifurcation.daughter_branch[0].start_diameter, bifurcation.daughter_branch[1].start_diameter) for bifurcation in bifurcation_list])
    
#    distribution_list.append('IBF_Taper')
#    distribution_type_list.append('basic')
#    entries_list.append([ibf.taper2 for  ibf in ibf_list])
#
#    distribution_list.append('Terminal_Taper')
#    distribution_type_list.append('basic')
#    entries_list.append([termination.branch.taper2 for  termination in terminal_list])
#    
#    distribution_list.append('Daughter_Ratio')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.daughter_ratio for bifurcation in bifurcation_list])   
#
#    distribution_list.append('Parent_Daughter_Ratio')
#    distribution_type_list.append('basic')
#    entries_list.append([max(bifurcation.daughter_branch[0].start_diameter, bifurcation.daughter_branch[1].start_diameter)/bifurcation.diameter for bifurcation in bifurcation_list])   


    distribution_list.append('Bifurcation_Amplitude_Remote(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_amp_remote for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    distribution_list.append('Bifurcation_Amplitude_Vector(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_amp_vector for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    distribution_list.append('Elevation_Remote(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_tilt_remote for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    distribution_list.append('Azimuth_Remote(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_torque_remote for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
        
    distribution_list.append('Orientation_Remote(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_twist_remote for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])


    distribution_list.append('Bifurcation_Amplitude_Local(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_amp_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    distribution_list.append('Elevation_Local(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_tilt_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

    distribution_list.append('Azimuth_Local(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_torque_local for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
        
    distribution_list.append('Orientation_Local(deg)')
    distribution_type_list.append('basic')
    entries_list.append([bifurcation.bif_twist_local for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])

    distribution_list.append('x_coordinate')
    distribution_type_list.append('conditional')
    entries_list.append([bifurcation.bifurcating_compartment.x for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])
    
    distribution_list.append('z_coordinate')
    distribution_type_list.append('conditional')
    entries_list.append([bifurcation.bifurcating_compartment.z for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])    

#    distribution_list.append('Terminal_Euclidean_to_soma_(um)')
#    distribution_type_list.append('conditional')
#    entries_list.append([termination.branch.euclidean_from_soma for termination in terminal_list])
#
#    distribution_list.append('IBF_Euclidean_to_soma_(um)')
#    distribution_type_list.append('conditional')
#    entries_list.append([ibf.euclidean_from_soma for ibf in ibf_list])

#    distribution_list.append('Pre-bifurcation_Diameter_(um)')
#    distribution_type_list.append('conditional')
#    entries_list.append([bifurcation.bifurcating_compartment.diameter for bifurcation in bifurcation_list])
#    
#    distribution_list.append('Pre-bifurcation_Pathlength_(um)')
#    distribution_type_list.append('conditional')
#    entries_list.append([bifurcation.parent_branch.pathlength for bifurcation in bifurcation_list])         

    distribution_list.append('Total_Dendritic_Length_(um)')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.total_dendritic_length for morphology in morphologies_list])
    
    distribution_list.append('Branch_Order')
    distribution_type_list.append('emergent')
    entries_list.append([branch.branch_order for branch in branch_list])

#    distribution_list.append('Branch_Order_1_(um)')
#    distribution_type_list.append('emergent')
#    entries_list.append([branch.branch_order_1 for branch in branch_list if branch.branch_order_1 != 0])
#    
#    distribution_list.append('Branch_Order_2_(um)')
#    distribution_type_list.append('emergent')
#    entries_list.append([branch.branch_order_2 for branch in branch_list if branch.branch_order_2 != 0])
#    
#    distribution_list.append('Branch_Order_3_(um)')
#    distribution_type_list.append('emergent')
#    entries_list.append([branch.branch_order_3 for branch in branch_list if branch.branch_order_3 != 0])
#    
#    distribution_list.append('Branch_Order_4_(um)')
#    distribution_type_list.append('emergent')
#    entries_list.append([branch.branch_order_4 for branch in branch_list if branch.branch_order_4 != 0])
    
    
    distribution_list.append('Number_of_Bifurcations')
    distribution_type_list.append('emergent')
    entries_list.append([morphology.num_bifurcations for morphology in morphologies_list])

#    distribution_list.append('Total_Surface_Area_(um^2)')
#    distribution_type_list.append('emergent')
#    entries_list.append([morphology.total_surface_area for morphology in morphologies_list])
    
#    distribution_list.append('Normalized_Total_Surface_Area_(um^2)')
#    distribution_type_list.append('emergent')
#    entries_list.append([morphology.normalized_total_surface_area for morphology in morphologies_list])    
#    
    
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

# Main
#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/claiborne/CNG version/*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/claiborne_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/bausch/CNG version/*.swc')
#distributions_file = '/home/zzchou/Documents/bausch_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/epsztein/CNG version/*.swc')
#distributions_file = '/home/zzchou/Documents/epsztein_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/claiborne/CNG version/*5199202*.swc')
#distributions_file = '/home/zzchou/Documents/claiborne_3319202_distributions.txt'

# file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/beining/CNG version/*Mature*supra*.swc')
# distributions_file = '/home/zzchou/Documents/Codes/beining_distributions.txt'

# file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/turner/CNG version/*.swc')
# distributions_file = '/home/zzchou/Documents/Codes/turner_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/pierce/CNG version/*.swc')
#distributions_file = '/home/zzchou/Documents/pierce_distributions.txt'
## #
#file_list = glob.glob('/home/zzchou/Documents/Data/test_output/morphology_tester_output*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/zmorph_tester_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/300 morphologies/*.swc')
#potential_outlier = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/claiborne/CNG version/*3319202*.swc')

#file_list = glob.glob('C:\Users\Zane\Desktop\LNeuron_workspace\output\output*.swc')
#distributions_file = 'C:\Users\Zane\Desktop\Morphologies\Codes\LNeuron_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/test_output/300zmorph_3/*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/zmorph_tester_distributions.txt'
#
#file_list = glob.glob('../Data/Morphologies/generator/*.swc')
#distributions_file = '../Data/Metrics/zmorph_tester_distributions.txt'
#
#file_list = glob.glob('../Data/Morphologies/l_neuron/*.swc')
#distributions_file = '../Data/Metrics/LNeuron_recreated_distributions.txt'
#
file_list = glob.glob('../Data/Morphologies/flat2/*.swc')
distributions_file = '../Data/Metrics/flat_distributions.txt'

#file_list = glob.glob('../Data/Morphologies/test_output/*.swc')
#distributions_file = '../Data/Metrics/zmorph_tester_distributions.txt'

#recreated_files = glob.glob('C:\Users\Zane\Desktop\LNeuron_workspace\output\*recreated*.swc')
#recreated_files = glob.glob('/home/zzchou/Documents/Data/test_output/300zmorph_3/*recreated*.swc')
#recreated_files = glob.glob('../Data/Morphologies/generator/*recreated*.swc')
#recreated_files = glob.glob('../Data/Morphologies/l_neuron/*recreated*.swc')
recreated_files = glob.glob('../Data/Morphologies/flat2/*recreated*.swc')
#recreated_files = glob.glob('../Data/Morphologies/test_output/*recreated*.swc')


for r_file in recreated_files:
    os.remove(r_file)

all_compartments_list = []
all_branches_list = []
all_bifurcations_list = []
all_bifurcationcounts = []
all_total_dendritic_length = []
all_terminations_list = []
all_morphologies_list = []
all_surface_area = []

pair_table = {}
id_pairs = []
ordered_pairs = []

avg_data = []

major_axis_mean = 380
major_axis_sd = 45
minor_axis_mean = 180
minor_axis_sd = 40
remodel_count = 1

max_file_size = 15000
min_file_size = 5
for i in range (0, len(file_list)):
    print file_list[i]

    file_size = int(os.path.getsize(file_list[i]))
    
#    if file_size > max_file_size or file_size < min_file_size:
#        print file_size
#        continue

    compartment_list, branch_list, bifurcation_list, termination_list = read_file(file_list[i])

#    if len(bifurcation_list) < 8:
#        continue

    bif_num_min = 8
    bif_num_max = 22
    bif_num_mid = (bif_num_min + bif_num_max)/2
    bif_num_range = bif_num_max - bif_num_min
    
    tdl_num_min = 2100
    tdl_num_max = 5200
    tdl_num_mid = (tdl_num_min + tdl_num_max)/2
    tdl_num_range = tdl_num_max - tdl_num_min
    
    
    bif_num = len(bifurcation_list)
    z = (bif_num - bif_num_mid)
    sigma = 3
    gauss_max = np.exp(-(0)**2/2)/(sigma*math.sqrt(2*math.pi))
    g_filter = np.exp(-(-z/sigma)**2/2)/(sigma*math.sqrt(2*math.pi))/gauss_max
    
#            if len(b_list) < 0 or len(b_list) > 48:
    print bif_num, g_filter
    if random.uniform(0,1) > g_filter:
#                print gaussian_filter
        continue

    if compartment_list == []:
        continue
    all_morphologies_list.append(Morphology(compartment_list, branch_list, bifurcation_list, termination_list))

#    tdl_num = all_morphologies_list[-1].total_dendritic_length
#    z = tdl_num - tdl_num_mid
#    sigma = 550
#    gauss_max = np.exp(-(0)**2/2)/(sigma*math.sqrt(2*math.pi))
#    g_filter = np.exp(-(-z/sigma)**2/2)/(sigma*math.sqrt(2*math.pi))/gauss_max
#    
#    if random.uniform(0, 1) > g_filter:
#        continue

    all_compartments_list = all_compartments_list  + compartment_list
    all_branches_list = all_branches_list + branch_list
    all_bifurcations_list = all_bifurcations_list + bifurcation_list
    all_terminations_list = all_terminations_list + termination_list
    total_dendritic_length = 0;

    for bif in bifurcation_list:
        if bif.parent_branch != []:
            if bif.parent_branch.parent_bifurcation != []:
                pair_table[bif.parent_branch.parent_bifurcation] = bif

    for bif in pair_table:
        id_pairs.append((bif.bifurcating_compartment.compartment_id, pair_table[bif].bifurcating_compartment.compartment_id))
        ordered_pairs.append((bif, pair_table[bif]))

    for branch in branch_list:
        total_dendritic_length = total_dendritic_length + branch.pathlength
    all_total_dendritic_length.append(total_dendritic_length)

    all_bifurcationcounts.append(len(bifurcation_list))
    if i >= 0 and i < 3:
  ##        avg_data.append(plot_compartment_feature(compartment_list))
        plot_branches(branch_list)
    

    rand_major_axis = 0
    rand_minor_axis = 0
    major_axis_min = 220
    major_axis_max = 540
    minor_axis_min = 105
    minor_axis_max = 270
  
    avg_branches = 37
    branch_diff = 26
    num_branches = len(branch_list)
    terms_x = [blah.x for blah in termination_list]
    terms_z = [blah.z for blah in termination_list]
    
    
    mj_axis = max(terms_x) - min(terms_x)
    mn_axis = max(terms_z) - min(terms_z)
    
    
    
    mj_axis = all_morphologies_list[-1].major_axis
    mn_axis = all_morphologies_list[-1].minor_axis
    
    rand_major_axis = mj_axis
    rand_minor_axis = mn_axis
    
    shift_axis1 = (num_branches - avg_branches)*((major_axis_max - major_axis_min)/(2*branch_diff))
    shift_axis2 = (num_branches - avg_branches)*((minor_axis_max - minor_axis_min)/(2*branch_diff))
   
#    print num_branches - avg_branches
    
    while rand_major_axis < major_axis_min or rand_major_axis > major_axis_max:
        rand_major_axis = random.gauss(major_axis_mean + shift_axis1, major_axis_sd)
    while rand_minor_axis < minor_axis_min or rand_minor_axis > minor_axis_max:    
        rand_minor_axis = random.gauss(minor_axis_mean + shift_axis2, minor_axis_sd)

    if all_morphologies_list[-1].major_axis == 0:
        major_scale = 1
    else:
        major_scale = rand_major_axis/all_morphologies_list[-1].major_axis
    if all_morphologies_list[-1].minor_axis == 0:
        minor_scale = 1
    else:
        minor_scale = rand_minor_axis/all_morphologies_list[-1].minor_axis
        
#    major_scale = rand_major_axis/all_morphologies_list[-1].major_axis
#    minor_scale = rand_minor_axis/all_morphologies_list[-1].minor_axis
    
    if major_scale == 1 and minor_scale ==1:
        print 'no remodeling necessary:', remodel_count
        remodel_count += 1
    
#    major_scale = 1
#    minor_scale = 1
    x_sd = 0.18
    z_sd = 0.18
    jitter_count = 5
    
#    output_title = '/home/zzchou/Documents/Data/test_output/300zmorph_3/{}{}{}'.format('morphology_recreated_', i+1, '.swc')
#    output_title = '../Data/Morphologies/generator/{}{}{}'.format('morphology_recreated_', i+1, '.swc')
    output_title = '../Data/Morphologies/l_neuron/{}{}{}'.format('morphology_recreated_', i+1, '.swc')
    output_title = '../Data/Morphologies/flat2/{}{}{}'.format('morphology_recreated_', i+1, '.swc')
#    output_title = '../Data/Morphologies/test_output/{}{}{}'.format('morphology_recreated_', i+1, '.swc')
#    output_title = 'C:\Users\Zane\Desktop\LNeuron_workspace\output\{}{}{}'.format('morphology_recreated_', i+1, '.swc')
    
    with open(output_title, 'w') as f:
#        true_origin = branch_list[0].compartment_list[0]
#        line = (true_origin.compartment_id, 1, true_origin.start_coordinates[0], true_origin.start_coordinates[1], true_origin.start_coordinates[2], stem_diameter, -1)
#        parent_id = 1
#        strs=' '.join(str(x) for x in line)
#        f.write(strs + '\n')

        
        c_id = 0
        compartment = branch_list[-1].compartment_list[0]
        compartment.compartment_id = c_id + 1
        line = (compartment.compartment_id, 1, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter/2, -1)
        strs=' '.join(str(x) for x in line)
        f.write(strs + '\n')  
        c_id = c_id + 1
        old_x_jitter = 0
        old_z_jitter = 0
        x_jitter = 0
        z_jitter = 0
        jitter_counter = 1
        for j in range(len(branch_list)):
#            print 'new branch'
            old_x_jitter = old_x_jitter + (jitter_counter-1)*x_jitter
            old_z_jitter = old_z_jitter + (jitter_counter-1)*z_jitter            
            x_jitter = random.gauss(0, x_sd)
            z_jitter = random.gauss(0, z_sd)         

            jitter_counter = 0

            x_jitter = random.gauss(0, x_sd)
            z_jitter = random.gauss(0, z_sd)
#            x_jitter = 1
#            z_jitter = 1
            branch = branch_list[len(branch_list) - j - 1]
            old_x_jitter = branch.compartment_list[0].jitter_x
            old_z_jitter = branch.compartment_list[0].jitter_z
            for compartment in branch.compartment_list:
#                jitter_counter = jitter_counter + 1
                if compartment != branch.compartment_list[0]:
                    if compartment != []:
                        compartment.compartment_id = c_id+1
                        c_id = c_id + 1
                        if compartment.parent_compartment != []:
                            parent_id = compartment.parent_compartment.compartment_id
                        else:
                            parent_id = -1

                        if jitter_counter >= jitter_count:
                            old_x_jitter = old_x_jitter + (jitter_counter)*x_jitter
                            old_z_jitter = old_z_jitter + (jitter_counter)*z_jitter
                            x_jitter = random.gauss(0, x_sd)
                            z_jitter = random.gauss(0, z_sd)
#                            x_jitter = 1
#                            z_jitter = 1
                            jitter_counter = 0
                            
#                        jitter_counter = jitter_counter + 1

#                        line = (compartment.compartment_id, 3, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter/2, parent_id)
#                        print old_x_jitter + jitter_counter*x_jitter, x_jitter
                        line = (compartment.compartment_id, 3, compartment.end_coordinates[0]*major_scale + old_x_jitter + jitter_counter*x_jitter, compartment.end_coordinates[1], compartment.end_coordinates[2]*minor_scale + old_z_jitter + jitter_counter*z_jitter, compartment.diameter/2, parent_id)
#                        line = (compartment.compartment_id, 3, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter/2, parent_id)
                        strs=' '.join(str(x) for x in line)
                        f.write(strs + '\n')
                        jitter_counter = jitter_counter + 1
                    parent_id = compartment.compartment_id
                if compartment == branch.compartment_list[-1]:
                    if branch.daughter_bifurcation != []:
                        for d_branch in branch.daughter_bifurcation.daughter_branch:
                            d_branch.compartment_list[0].jitter_x = old_x_jitter + (jitter_counter-1)*x_jitter
                            d_branch.compartment_list[0].jitter_z = old_z_jitter + (jitter_counter-1)*z_jitter
                       
    f.closed

#    reconstruct_morphology(compartment_list, branch_list, bifurcation_list, termination_list)


#print 'avgs'
#plot_avg_data(avg_data)


ibf_branches = [branch for branch in all_branches_list if branch.daughter_bifurcation != []]
terminal_branches = [branch for branch in all_branches_list if branch.daughter_bifurcation == []]
stem_compartments = []
for branch in all_branches_list:
    if branch.branch_order == 0:
        stem_compartments = stem_compartments + branch.compartment_list

branch_order_0 = [branch for branch in all_branches_list if branch.branch_order == 0]

euclidean_order_1 = []
euclidean_order_2 = []
euclidean_order_3 = []
euclidean_order_4 = []
euclidean_order_5 = []
euclidean_order_6 = []

for branch in all_branches_list:
    if branch.euclidean_from_soma - branch.euclidean < 50:
        euclidean_order_1.append(branch)
    elif branch.euclidean_from_soma - branch.euclidean < 100:
        euclidean_order_2.append(branch)
    elif branch.euclidean_from_soma - branch.euclidean < 150:
        euclidean_order_3.append(branch)
    elif branch.euclidean_from_soma - branch.euclidean < 200:
        euclidean_order_4.append(branch)    
    elif branch.euclidean_from_soma - branch.euclidean < 250:
        euclidean_order_5.append(branch)
    elif branch.euclidean_from_soma - branch.euclidean >= 250:
        euclidean_order_6.append(branch)
        
distribution_list = []
distribution_type_list = []
entries_list = []
       
all_surface_area= [morphology.total_surface_area for morphology in all_morphologies_list]

all_bifurcation_flag_list = []
all_GF_list = []
all_GF_x_raw_list = []
all_GF_y_raw_list = []
all_GF_z_raw_list = []
all_GF_neighbor_raw_list = []
GF_tuple_list = []

for morphology in all_morphologies_list:
    for entry in morphology.GF_list:
        all_GF_list.append(entry)
    for entry2 in morphology.branch_flag_list:
        all_bifurcation_flag_list.append(entry2)
    for entry3 in morphology.GF_x_raw_list:
        all_GF_x_raw_list.append(entry3)    
    for entry4 in morphology.GF_y_raw_list:
        all_GF_y_raw_list.append(entry4)
    for entry5 in morphology.GF_z_raw_list:
        all_GF_z_raw_list.append(entry5)
    for entry6 in morphology.GF_neighbor_raw_list:
        all_GF_neighbor_raw_list.append(entry6)

for i in range(0, len(all_GF_list)):
    GF_tuple_list.append([all_GF_list[i], all_bifurcation_flag_list[i], all_GF_x_raw_list[i], all_GF_y_raw_list[i], all_GF_z_raw_list[i], all_GF_neighbor_raw_list[i]])

#sorted_GF_tuple_list = sorted(GF_tuple_list)
#sorted_GF_list = [x[0] for x in sorted_GF_tuple_list]
#sorted_bif_flag_list = [x[1] for x in sorted_GF_tuple_list]
#
#sorted_GF_tuple_list2 = sorted(GF_tuple_list, key=lambda x: x[2])
#sorted_GF_x_list = [x[2] for x in sorted_GF_tuple_list2]
#sorted_bif_flag_list2 = [x[1] for x in sorted_GF_tuple_list2]
#
#sorted_GF_tuple_list3 = sorted(GF_tuple_list, key=lambda x: x[3])
#sorted_GF_y_list = [x[3] for x in sorted_GF_tuple_list3]
#sorted_bif_flag_list3 = [x[1] for x in sorted_GF_tuple_list3]
#
#sorted_GF_tuple_list4 = sorted(GF_tuple_list, key=lambda x: x[4])
#sorted_GF_z_list = [x[4] for x in sorted_GF_tuple_list4]
#sorted_bif_flag_list4 = [x[1] for x in sorted_GF_tuple_list4]
#
#sorted_GF_tuple_list5 = sorted(GF_tuple_list, key=lambda x: x[5])
#sorted_GF_neighbor_list = [x[5] for x in sorted_GF_tuple_list5]
#sorted_bif_flag_list5 = [x[1] for x in sorted_GF_tuple_list5]
#
#sorted_GF_tuple_list2 = sorted(GF_tuple_list, key=lambda x: x[2])
#sorted_GF_x_list = [x[2] for x in sorted_GF_tuple_list2]
#sorted_bif_flag_list2 = [x[1] for x in sorted_GF_tuple_list2]
#sorted_GF_y_list = [x[3] for x in sorted_GF_tuple_list2]
#sorted_GF_z_list = [x[4] for x in sorted_GF_tuple_list2]
#sorted_GF_neighbor_list = [x[5] for x in sorted_GF_tuple_list2]
#
#GF_x_max = max(sorted_GF_x_list)
#GF_y_max = max(sorted_GF_y_list)
#GF_z_max = max(sorted_GF_z_list)
#GF_neighbor_max = max(sorted_GF_neighbor_list)
#
#GF_bin1 = 0
#GF_bin2 = 0
#GF_bin3 = 0
#GF_bin4 = 0
#GF_bin5 = 0
#BF_bin1 = 0
#BF_bin2 = 0
#BF_bin3 = 0
#BF_bin4 = 0
#BF_bin5 = 0

#for i in range(0, len(all_GF_list)):
#    if all_GF_neighbor_raw_list[i] > 0 and all_GF_neighbor_raw_list[i]< 0.01:
#        GF_bin1 = GF_bin1 + 1
#        if all_bifurcation_flag_list[i] == 1:
#            BF_bin1 = BF_bin1 + 1
#    elif all_GF_neighbor_raw_list[i] < 0.02:
#        GF_bin2 = GF_bin2 + 1
#        if all_bifurcation_flag_list[i] == 1:
#            BF_bin2 = BF_bin2 + 1
#    elif all_GF_neighbor_raw_list[i] < 0.03:
#        GF_bin3 = GF_bin3 + 1
#        if all_bifurcation_flag_list[i] == 1:
#            BF_bin3 = BF_bin3 + 1
#    elif all_GF_neighbor_raw_list[i] < 0.04:
#        GF_bin4 = GF_bin4 + 1
#        if all_bifurcation_flag_list[i] == 1:
#            BF_bin4 = BF_bin4 + 1
#    elif all_GF_neighbor_raw_list[i] < 0.05:
#        GF_bin5 = GF_bin5 + 1
#        if all_bifurcation_flag_list[i] == 1:
#            BF_bin5 = BF_bin5 + 1   
#
#print 'bin1', (BF_bin1/GF_bin1)
#print 'bin2', (BF_bin2/GF_bin2)
#print 'bin3', (BF_bin3/GF_bin3)
#print 'bin4', (BF_bin4/GF_bin4)
#print 'bin5', (BF_bin5/GF_bin5)
#print all_GF_list

#N = 1000
#Nmax = 1000
#
##plt.figure()
#nd_array = np.ndarray(shape=(Nmax,Nmax,Nmax,Nmax), dtype=float)
#for i in range(0, len(sorted_GF_x_list)):
#    index1 = int(sorted_GF_x_list[i]*Nmax)
#    index2 = int(sorted_GF_y_list[i]*Nmax)
#    index3 = int(sorted_GF_z_list[i]*Nmax)
#    index4 = int(sorted_GF_neighbor_list[i]*Nmax)
#    
#    nd_array[index1, index2, index3, index4] = sorted_bif_flag_list2
#    
#test_array = np.asarray(sorted_bif_flag_list)
#test_out = np.convolve(nd_array, np.ones((N,N,N,N,))/N)[(N-1):]
#plt.plot(sorted_GF_list, test_out)
#plt.xlim([0,1])
#plt.show()

#plt.figure()
#test_array = np.asarray(sorted_bif_flag_list)
#test_out = np.convolve(test_array, np.ones((N,))/N)[(N-1):]
#plt.plot(sorted_GF_list, test_out)
#plt.xlim([0,1])
#plt.show()
#
#plt.figure()
#test_array = np.asarray(sorted_bif_flag_list2)
#test_out = np.convolve(test_array, np.ones((N,))/N)[(N-1):]
## print len(test_array), len(test_out)
#
#plt.plot(sorted_GF_x_list, test_out)
#plt.xlim([0, 1])
#plt.show()
#
#plt.figure()
#test_array = np.asarray(sorted_bif_flag_list3)
#test_out = np.convolve(test_array, np.ones((N,))/N)[(N-1):]
#plt.plot(sorted_GF_y_list, test_out)
#plt.xlim([0,1])
#plt.show()
#
#plt.figure()
#test_array = np.asarray(sorted_bif_flag_list4)
#test_out = np.convolve(test_array, np.ones((N,))/N)[(N-1):]
#plt.plot(sorted_GF_z_list, test_out)
#plt.xlim([0,1])
#plt.show()
#
#plt.figure()
#test_array = np.asarray(sorted_bif_flag_list5)
#test_out = np.convolve(test_array, np.ones((N,))/N)[(N-1):]
#plt.plot(sorted_GF_neighbor_list, test_out)
#plt.xlim([0,1])
#plt.show()



    
#save_distributions(all_compartments_list, all_branches_list, all_bifurcations_list, all_terminations_list, all_morphologies_list, distribution_list, distribution_type_list, entries_list)
 
#plot_euclidean_orders([euclidean_order_1, euclidean_order_2, euclidean_order_3, euclidean_order_4, euclidean_order_5, euclidean_order_6])

#plot_pathlength_feature(all_branches_list)

#plot_pathlength_feature(terminal_branches)

#plot_compartment_feature(all_compartments_list)

# plot_morphology_feature(all_surface_area)

#plot_morphology_feature(all_total_dendritic_length)

#plot_branch_feature(all_branches_list)

#plot_pathlength_to_soma(terminal_branches)

#plot_euclidean_from_soma(terminal_branches)

#plot_bifurcations(bifurcation_list)

#plot_bifurcation_feature(all_bifurcations_list)

#plot_branches(all_branches_list)

#plot_pairs(ordered_pairs)



#with open(distributions_file, 'w') as f:
#    for i in range(0, len(distribution_list)):
#        entries = ' '.join(str(x) for x in entries_list[i])
#        strs = distribution_list[i] + ' ' + distribution_type_list[i] + ' ' + entries
#        f.write(strs + '\n')
#    f.closed