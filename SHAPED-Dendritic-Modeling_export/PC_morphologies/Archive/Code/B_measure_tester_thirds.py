from __future__ import division

import glob
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import pandas
import random

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
                daughter.rot_x, daughter.rot_y, daughter.rot_z = rotation_angles(daughter.vector, [0,0,1])                
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
        
        print self.max_euclidean
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
#        print self.total_surface_area
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
        bin_edges=np.histogram(bif_angle_list, bins=50)[1] #get the bin edges

#        plt.hist(bif_angle_list, alpha=0.7, normed=True, bins=bin_edges)
#        plt.plot(sholl_bins, self.sholl_counts)
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
        if branch.end_compartment.section_type != 1:
            for compartment in branch.compartment_list:
                x_list.append(compartment.end_coordinates[0])
                y_list.append(compartment.end_coordinates[1])
                z_list.append(compartment.end_coordinates[2])
        plt.plot(x_list, z_list, y_list, color = 'black')
        x_list = []
        y_list = []
        z_list = []

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
    
def save_distributions(compartment_list, branch_list, bifurcation_list, terminal_list, morphologies_list, distribution_list, distribution_type_list, entries_list):
    
    stem_list = [stem for stem in branch_list if stem.branch_order == 0]
    ibf_list = [ibf.parent_branch for ibf in bifurcation_list if ibf.parent_branch.branch_order != 0]
    terminal_branches = [branch for branch in branch_list if branch.daughter_bifurcation == []]
    
    distribution_list.append('Stem_Pathlength_(um)')
    distribution_type_list.append('basic')
    entries_list.append([stem.pathlength for stem in stem_list]) 
    
#    distribution_list.append('All_pathlength_GCL_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 50])
#    
#    distribution_list.append('All_pathlength_IML_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma <= 150 and ibf.euclidean_from_soma > 50])
#    
#    distribution_list.append('All_pathlength_MML_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([ibf.pathlength for ibf in branch_list if (ibf.euclidean_from_soma > 150 and ibf.euclidean_from_soma <= 250)])
#
#    distribution_list.append('All_pathlength_OML_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([ibf.pathlength for ibf in branch_list if ibf.euclidean_from_soma > 250])    
#    
#    distribution_list.append('Terminal_Pathlength_GCL_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma <= 50])
#
#    distribution_list.append('Terminal_Pathlength_IML_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 50 and termination.branch.euclidean_from_soma <= 150])
#    
#    distribution_list.append('Terminal_pathlength_MML_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([termination.branch.pathlength for termination in terminal_list if (termination.branch.euclidean_from_soma > 150 and termination.branch.euclidean_from_soma <= 250)])
#
#    distribution_list.append('Terminal_pathlength_OML_(um)')
#    distribution_type_list.append('basic')
#    entries_list.append([termination.branch.pathlength for termination in terminal_list if termination.branch.euclidean_from_soma > 250])  

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
    dummy = []
    for branch in branch_list:
        if branch.parent_bifurcation != []:
            branch_distance_distribution.append(branch.pathlength_to_soma - branch.pathlength)
            pathlength_distribution.append(branch.pathlength)
            rotx_distribution.append(branch.rot_x)
            rotz_distribution.append(branch.rot_z)
            previous_branch_length.append(branch.parent_bifurcation.parent_branch.pathlength)
#            previous_bif_amp.append(branch.parent_bifurcation.bif_amp_remote)
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
#            previous_branch_length.append(0)
#            previous_bif_amp.append(0)
#            previous_rot_x.append(0)
#            previous_rot_z.append(0)
#            previous_x.append(0)
#            previous_z.append(0)
#            bif_amp_distribution2.append(0)
            dummy.append(random.random())
            
    normalized_distance = [x*1000/max_distance for x in branch_distance_distribution]

    

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

    bif_amp_distribution = [bifurcation.bif_amp_vector for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1]
    bif_distance_distribution = [bifurcation.daughter_branch[0].pathlength_to_soma for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1]

#    distribution_list.append('Elevation_Remote(deg)')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.bif_tilt_remote for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])
#
#    distribution_list.append('Azimuth_Remote(deg)')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.bif_torque_remote for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
#        
#    distribution_list.append('Orientation_Remote(deg)')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.bif_twist_remote for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])


#    distribution_list.append('Bifurcation_Amplitude_Local(deg)')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.bif_amp_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])

#    distribution_list.append('Elevation_Local(deg)')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.bif_tilt_local for bifurcation in bifurcation_list if bifurcation.bif_tilt_remote != -1])
#
#    distribution_list.append('Azimuth_Local(deg)')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.bif_torque_local for bifurcation in bifurcation_list if bifurcation.bif_torque_remote != -1])
#        
#    distribution_list.append('Orientation_Local(deg)')
#    distribution_type_list.append('basic')
#    entries_list.append([bifurcation.bif_twist_local for bifurcation in bifurcation_list if bifurcation.bif_twist_remote != -1])

    distribution_list.append('rot_x')
    distribution_type_list.append('basic')
    entries_list.append([branch.rot_x for branch in branch_list])

#    rotx_distribution = [branch.rot_x for branch in branch_list]

    distribution_list.append('rot_z')
    distribution_type_list.append('basic')
    entries_list.append([branch.rot_z for branch in branch_list])
    
#    rotz_distribution = [branch.rot_z for branch in branch_list]
    
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
    entries_list.append([morphology.sholl_normed_counts for morphology in morphologies_list])

#    distribution_list.append('max_euclidean')
#    distribution_type_list.append('emergent')
#    entries_list.append([morphology.max_euclidean for morphology in morphologies_list])

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

    print tdl_bins
    tmp_list = [morphology.tdl_distribution for morphology in morphologies_list]
    plt.plot(tdl_bins, [sum(x)/len(morphologies_list) for x in zip(*tmp_list)], label='Generated TDL')    
    plt.title('TDL for Generated')
    plt.legend()

    print len(morphologies_list)








#
#    tmp_list = [morphology.sholl_counts for morphology in morphologies_list]
#    tmp_list2 = [morphology.ysholl_counts for morphology in morphologies_list]
#    plt.plot(sholl_bins, [sum(x)/len(morphologies_list) for x in zip(*tmp_list)], label='Generated Sholl')    
#    plt.plot(sholl_bins, [sum(x)/len(morphologies_list) for x in zip(*tmp_list2)], label='Generated YSholl')
#    plt.title('Sholl for Generated')
#    plt.legend()
    
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
#    data = [normalized_distance, previous_bif_amp, previous_rot_x, previous_rot_z, previous_x, previous_z, pathlength_distribution, bif_amp_distribution2, rotx_distribution, rotz_distribution]
    
    data_length = len(data)
    data_length2 = len(data[0])
    
#    for i in range(5):
#        for j in range(data_length2):
#            data[i][j] = random.random()
    
    names = ['branch_distance_distribution', 'previous_bif_amp', 'previous_rot_x', 'previous_rot_z', 'previous_x', 'previous_z', 'pathlength_distribution', 'bif_amp_distribution2', 'rotx_distribution', 'rotz_distribution']
#    datfra = pandas.DataFrame(np.array(data).reshape(data_length2, data_length), columns = names)
#    correlations = datfra.corr()
    
    correlation_list = [x[:] for x in [[0]*data_length] * data_length]
    for i in range(0, data_length):
        for j in range(0, data_length):
#            print i, j
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
    
#    i = 0
#    j = 6
#    
#    print len(data[i])
#    for idx in range(len(data[i])):
#        if data[i][len(data[i]) - idx -1] > 550:
#            print data[i][len(data[i]) - idx -1]
#            del data[i][len(data[i])-idx -1]                    
#            del data[j][len(data[i])-idx -1]
#        
#    print len(data[i])
    
    
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
            print xedges, yedges
            H = H.T
            H_max = np.amax(H)
            column_max = np.amax(H, axis=0)
            if x == 0:
                weights = [0 for x in column_max]
            else:
                weights = [H_max/x for x in column_max]
            for k in range(len(weights)):
                H[:,k] *= weights[k]
            plt.imshow(H, interpolation='nearest', aspect='auto', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap='viridis', vmin=0, vmax=H_max)
            plt.colorbar()
#            plt.title(names[i] + ' vs ' + names[j])
    
#    fig = plt.figure()
#    plt.scatter(previous_rot_x, rotx_distribution)
#    plt.title('rotx')
#    fig = plt.figure()
#    plt.scatter(previous_rot_z, rotz_distribution)
#    plt.title('rotz')
    
# Main
#file_list = glob.glob('../Data/Morphologies/claiborne/CNG version/*.swc')
#distributions_file = '../Data/Metrics/claiborne_distributions.txt'
##
#file_list = glob.glob('../Data/Morphologies/generator/*recreated*.swc')
#distributions_file = '../Data/Metrics/zmorph_tester_distributions.txt'

file_list = glob.glob('../Data/Morphologies/l_neuron/*recreated*.swc')
distributions_file = '../Data/Metrics/LNeuron_recreated_distributions.txt'

#file_list = glob.glob('../Data/Morphologies/flat2/*recreated*.swc')
#distributions_file = '../Data/Metrics/flat_distributions.txt'



#file_list = glob.glob('/home/zzchou/Desktop/codes2018/Data/neuron_nmo/claiborne/stem1/*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/claiborne_1stem_distributions.txt'

#file_list = glob.glob('/home/zzchou/Desktop/codes2018/Data/neuron_nmo/claiborne/stem2/*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/claiborne_2stem_distributions.txt'
#
#file_list = glob.glob('/home/zzchou/Desktop/codes2018/Data/neuron_nmo/claiborne/stem34/*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/claiborne_34stem_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/bausch/CNG version/*.swc')
#distributions_file = '/home/zzchou/Documents/bausch_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/epsztein/CNG version/*.swc')
#distributions_file = '/home/zzchou/Documents/epsztein_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/claiborne/CNG version/*124894*.swc')
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

#file_list = glob.glob('/home/zzchou/Desktop/181227/c300t200/morphology_tester_output*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/zmorph_tester_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/test_output/morphology_recreated*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/zmorph_tester_distributions.txt'

#file_list = glob.glob('/home/zzchou/Desktop/codes2018/Data/tau5/*output*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/zmorph_tester_distributions_t5.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/testing/*output*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/LNeuron_recreated_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/testing/*recreated*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/LNeuron_recreated_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/test_output/claiborne_recreated*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/claiborne_distributions.txt'

#file_list = glob.glob('/home/zzchou/Documents/Data/300 morphologies/*.swc')
#potential_outlier = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/claiborne/CNG version/*3319202*.swc')

#file_list = glob.glob('/home/zzchou/Downloads/output.swc')
#distributions_file = '/home/zzchou/Downloads/output_distributions.txt'
            
#file_list = glob.glob('/home/zzchou/Documents/Data/test_output/300zmorph_3/morphology_recreated*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/zmorph_tester_distributions.txt'
            
#file_list = glob.glob('/home/zzchou/Documents/Data/test_output/300zmorph_2/morphology_tester_output*.swc')
#distributions_file = '/home/zzchou/Documents/Codes/zmorph_tester_distributions.txt'


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

for i in range (0, len(file_list)):
    print file_list[i]
    compartment_list, branch_list, bifurcation_list, termination_list = read_file(file_list[i])

    if i > 50:
        continue

    if compartment_list == [] and branch_list == []:
        continue

    all_morphologies_list.append(Morphology(compartment_list, branch_list, bifurcation_list, termination_list))

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
    if i >= 0 and i < 40:
  ##        avg_data.append(plot_compartment_feature(compartment_list))
        print i, all_morphologies_list[-1].num_stems
        plot_branches(branch_list)
#        plt.grid(None)

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

#for morphology in all_morphologies_list:
#    for entry in morphology.GF_list:
#        all_GF_list.append(entry)
#    for entry2 in morphology.branch_flag_list:
#        all_bifurcation_flag_list.append(entry2)
#    for entry3 in morphology.GF_x_raw_list:
#        all_GF_x_raw_list.append(entry3)    
#    for entry4 in morphology.GF_y_raw_list:
#        all_GF_y_raw_list.append(entry4)
#    for entry5 in morphology.GF_z_raw_list:
#        all_GF_z_raw_list.append(entry5)
#    for entry6 in morphology.GF_neighbor_raw_list:
#        all_GF_neighbor_raw_list.append(entry6)
#
#for i in range(0, len(all_GF_list)):
#    GF_tuple_list.append([all_GF_list[i], all_bifurcation_flag_list[i], all_GF_x_raw_list[i], all_GF_y_raw_list[i], all_GF_z_raw_list[i], all_GF_neighbor_raw_list[i]])
#
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

save_distributions(all_compartments_list, all_branches_list, all_bifurcations_list, all_terminations_list, all_morphologies_list, distribution_list, distribution_type_list, entries_list)

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



with open(distributions_file, 'w') as f:
    for i in range(0, len(distribution_list)):
        entries = ' '.join(str(x) for x in entries_list[i])
        strs = distribution_list[i] + ' ' + distribution_type_list[i] + ' ' + entries
        f.write(strs + '\n')
    f.closed