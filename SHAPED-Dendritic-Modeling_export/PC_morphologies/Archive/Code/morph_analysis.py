from __future__ import division
import time
import glob
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import pandas
import random

from angle_noise_updated import angle_noise
from angle_noise_updated import angle_noise_old
from B_morph_170811_2 import make_compartments

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

def unit_vector((u1, u2, u3)):
    magnitude = euclidean_distance((u1, u2, u3), (0, 0, 0))
    
    if magnitude != 0:  
        return [u1/magnitude, u2/magnitude, u3/magnitude]
    else:
        return [0, 0, 0]
    
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
#        self.branch_order_1 = 0
#        self.branch_order_2 = 0
#        self.branch_order_3 = 0
#        self.branch_order_4 = 0

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

    def set_branch_order(self):
        if self.daughter_bifurcation != []:
            for daughter in self.daughter_bifurcation.daughter_branch:
                daughter.branch_order = self.branch_order + 1
                daughter.set_branch_order()
                
    def set_pathlength_to_soma(self):
        if self.daughter_bifurcation != []:
            for daughter in self.daughter_bifurcation.daughter_branch:
                daughter.pathlength_to_soma = daughter.pathlength + self.pathlength_to_soma
                daughter.vector = [daughter.end_coordinates[i] - daughter.start_coordinates[i] for i in range(0, 3)]
                daughter.rot_x, daughter.rot_y, daughter.rot_z = rotation_angles_old(daughter.vector, [0,0,1])                
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
        
        max_euclidean = 0
        max_y = 0
        
        for compartment in self.compartment_list:
            local_GF, GF_x_raw, GF_y_raw, GF_z_raw, GF_neighbor_raw = get_local_GF(compartment.end_coordinates, self.branch_list)
            self.GF_list.append(local_GF)
            self.GF_x_raw_list.append(GF_x_raw)
            self.GF_y_raw_list.append(GF_y_raw)
            self.GF_z_raw_list.append(GF_z_raw)
            self.GF_neighbor_raw_list.append(GF_neighbor_raw)
            if compartment.euclidean_from_soma > max_euclidean:
                max_euclidean = compartment.euclidean_from_soma
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
        
        for branch in branch_list:
            self.total_dendritic_length = self.total_dendritic_length + branch.pathlength

        self.length_scale_factor = self.total_dendritic_length/100    
        for compartment in compartment_list:
            self.total_surface_area = self.total_surface_area + math.pi*compartment.diameter*compartment.length
            self.normalized_total_surface_area = self.normalized_total_surface_area + math.pi*compartment.diameter*compartment.length/self.length_scale_factor
        for termination in terminal_list:
            self.total_surface_area = self.total_surface_area + math.pi*(0.5*termination.diameter)**2
        self.num_bifurcations = len(bifurcation_list)
        
        normed_bins = np.linspace(0, max_euclidean, len(sholl_bins))
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
        self.sholl_counts = [x for x in self.sholl_counts]         
        
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

def calculate_absolute_angles(bifurcation_list):
#    print bifurcation_list
    for bifurcation in bifurcation_list:
#        print bifurcation.parent_branch, bifurcation.parent_branch.parent_bifurcation
        start_coordinates = bifurcation.bifurcating_compartment.end_coordinates
        end_coordinates1 = bifurcation.remote_daughters[0].end_coordinates
        end_coordinates2 = bifurcation.remote_daughters[1].end_coordinates

        bifurcation.vector1 = [end_coordinates1[i] - start_coordinates[i] for i in range(0, 3)]
        bifurcation.vector2 = [end_coordinates2[i] - start_coordinates[i] for i in range(0, 3)]  
        
        vector1 = unit_vector(bifurcation.vector1)
        vector2 = unit_vector(bifurcation.vector2)

        bif_norm = cross_product(vector1, vector2)
        bif_norm_mag = math.sqrt(sum(i**2 for i in bif_norm))

        if bif_norm_mag == 0:
            bifurcation.bif_normal = (0,1,0)
        else:
            bifurcation.bif_normal = (bif_norm[0]/bif_norm_mag, bif_norm[1]/bif_norm_mag, bif_norm[2]/bif_norm_mag)
            
        parent_vector = bifurcation.parent_branch.vector
        soma_vector = [bifurcation.bifurcating_compartment.x, bifurcation.bifurcating_compartment.y, bifurcation.bifurcating_compartment.z]
        soma_vector[1] = 0
        
#        print soma_vector

        bifurcation.bif_amp_remote = math.acos(dot_product(vector1, vector2)/magnitude_product(vector1, vector2))*180/math.pi
        bifurcation.bif_midline = ((vector1[0]+vector2[0])/2, (vector1[1]+vector2[1])/2, (vector1[2]+vector2[2])/2)

        if bifurcation.parent_branch.parent_bifurcation != []:
#            print bifurcation.bif_midline, bifurcation.bif_normal
            bifurcation.bif_normal = [0,0,1]
            rot_x1, rot_y1, rot_z1 = rotation_angles_old(bifurcation.bif_midline, bifurcation.bif_normal)
#            print rot_x1, rot_y1, rot_z1
            rot_x2, rot_y2, rot_z2 = rotation_angles_old(bifurcation.parent_branch.vector, bifurcation.parent_branch.parent_bifurcation.bif_normal)
        
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
            
#            bifurcation.bif_torque_remote, bifurcation.bif_tilt_remote, bifurcation.bif_twist_remote = rotation_angles(bif_midline2, bif_norm2)
            bifurcation.bif_tilt_remote = rot_x1*180/math.pi
            bifurcation.bif_torque_remote = rot_y1*180/math.pi
            bifurcation.bif_twist_remote = rot_z1*180/math.pi
            
#            print bifurcation.bif_torque_remote, bifurcation.bif_tilt_remote, bifurcation.bif_twist_remote
        else:
#            print 'yes'
            rot_x1, rot_y1, rot_z1 = rotation_angles_old(bifurcation.bif_midline, [0,0,1])
            
#            print rot_x1, rot_y1, rot_z1
            
            bifurcation.bif_tilt_remote = rot_x1*180/math.pi
            bifurcation.bif_torque_remote = rot_y1*180/math.pi
            bifurcation.bif_twist_remote = rot_z1*180/math.pi            
            
            
        
    return

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
            
            print bifurcation.bifurcating_compartment.end_coordinates, bifurcation.bif_amp_remote, bifurcation.bif_tilt_remote, bifurcation.bif_torque_remote, bifurcation.bif_twist_remote
    return

def rotation_angles_old(u, norm):
#    n1 = (1, 0, 0)
#    n2 = (0, 1, 0)
#
#    w1 = [dot_product(u,n1)/magnitude_product(n1,n1), 0, 0]
#    w2 = [0, dot_product(u,n2)/magnitude_product(n2,n2), 0]
#
#    proj_yz = [u[i] - w1[i] for i in range(0,3)]
#    proj_xz = [u[i] - w2[i] for i in range(0,3)]

    proj_yz = u
    n = (0,1,0)
    
    n2 = (0,0,1)
    
    if proj_yz == [0, 0, 0]:
        rot_x = 0
    else:
        rot_x = math.acos(dot_product(proj_yz, n)/magnitude_product(proj_yz, n))
        
    if proj_yz[2] > 0:
        rot_x = -rot_x
        
    Rr_x = np.matrix( ((1, 0, 0), (0, math.cos(rot_x), -math.sin(rot_x)), (0, math.sin(rot_x), math.cos(rot_x))) )
    
    new_u = np.ndarray.tolist(np.dot(n, Rr_x))[0]
    proj_xz = [u[0], 0, u[2]]
    proj_new_u = [new_u[0], 0, new_u[2]]
#    print new_u
    if new_u == [0,0,0]:
        rot_y = 0
    else:
        rot_y = math.acos(dot_product(proj_new_u, proj_xz)/magnitude_product(proj_new_u, proj_xz))
    if new_u[0] < 0:
        rot_y = -rot_y

    
    Rr_y = np.matrix( ((math.cos(-rot_y), -math.sin(-rot_y), 0), (0, 1, 0), (math.sin(-rot_y), math.cos(-rot_y), 0)) )
    
    Rr_xy = Rr_x*Rr_y
    
    
    if norm == (0, 0, 0):
        norm = (0, 0, 1)
    
    norm = [norm[i]/math.sqrt(magnitude_product(norm, norm)) for i in range(0,3)]
    new_norm = np.ndarray.tolist(np.dot(n2, Rr_xy))[0]
    n_norm = [0, 0, 1]
    
    rot_z = math.acos(dot_product(n2, norm)/magnitude_product(n2, norm))
    if new_norm[0] < 0:
        rot_z = -rot_z
    Rr_z = np.matrix( ((math.cos(-rot_z), 0, -math.sin(-rot_z)), (math.sin(-rot_z), 0, math.cos(-rot_z)), (0, 0, 1)) )

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
#        start_match = [bif for bif in bifurcation_list if int(bif.bifurcation_id) == int(branch.start_compartment.compartment_id)]
#        if start_match != []:
#            branch.set_parent_bifurcation(start_match[0])
            
        end_match = [bif for bif in bifurcation_list if int(bif.bifurcation_id) == int(branch.end_compartment.compartment_id)]
        if end_match != []:
            branch.set_daughter_bifurcation(end_match[0])


            
#    for branch in branches:
#        branch_order = 0
#        branch_order_4 = 0
#        branch_order_3 = 0
#        branch_order_2 = 0
#        branch_order_1 = 0
#        
#        parent_bif = branch.parent_bifurcation
#
#        while parent_bif != []:
#            if parent_bif.euclidean_from_soma > 250:
#                branch_order_4 = branch_order_4 + 1
#            elif parent_bif.euclidean_from_soma > 150:
#                branch_order_3 = branch_order_3 + 1
#            elif parent_bif.euclidean_from_soma > 50:
#                branch_order_2 = branch_order_2 + 1
#            else:
#                branch_order_1 = branch_order_1 + 1
#            parent_bif = parent_bif.parent_branch.parent_bifurcation
#            branch_order = branch_order + 1
#
#        branch.branch_order = branch_order
#        branch.branch_order_4 = branch_order_4
#        branch.branch_order_3 = branch_order_3
#        branch.branch_order_2 = branch_order_2
#        branch.branch_order_1 = branch_order_1
#        
#        current_compartment = branch.end_compartment
#        
#        while current_compartment.compartment_id != branch.start_compartment.compartment_id:
#
#            current_compartment.branch_order = branch_order
#            current_compartment = current_compartment.parent_compartment
#
    for branch in branches:
        if branch.branch_order == 0:
            branch.set_pathlength_to_soma()
                
    return branches

def plot_dendrogram(branch_list):
    fig = plt.figure()
#    terminations = 0
#    branch_order_counts = [0]*10
#    
#    for branch in branch_list:
#        if branch.daughter_bifurcation == []:
#            terminations = terminations + 1
#        branch_order_counts[branch.branch_order] += 1
    
#    i_max = max(branch_order_counts)
#    for i in range(0, len(branch_order_counts)):
#        j_max = branch_order_counts[i]
#        for j in range(0, branch_order_counts[i]):
#            plt.scatter(j/(j_max)+ 0.5/j_max, i)
#    plt.plot(branch_order_counts)

    

    index = 0
    dendro_x = []
    dendro_y = []
    for branch in branch_list:
        if branch.daughter_bifurcation == []:
            dendro_x.append(index)
            dendro_y.append(branch.pathlength_to_soma - branch.pathlength)
            x_list = [index, index]
            index += 1
            y_list = [branch.pathlength_to_soma - branch.pathlength, branch.pathlength_to_soma]
            plt.plot(x_list, y_list, color = 'black')
            
    tmp_flag = 0
    tmp_list = branch_list
    new_tmp_list = branch_list
    while tmp_flag == 0:
        new_tmp_list = []
        for branch in tmp_list:
            if branch.daughter_bifurcation != []:
                ind_matches = [jj for jj, x in enumerate(dendro_y) if round(x, 5) == round(branch.pathlength_to_soma, 5)]
#                print len(branch_list), len(dendro_y)

                if len(ind_matches) == 2:
                    cross_bar_x = [dendro_x[ii] for ii in ind_matches]
                    cross_bar_y = [branch.pathlength_to_soma, branch.pathlength_to_soma]
                                    
    
                
                    x_ind = np.mean([dendro_x[ii] for ii in ind_matches])
                    x_list = [x_ind, x_ind]
                    
                    dendro_x.append(x_ind)
                    dendro_y.append(branch.pathlength_to_soma - branch.pathlength)
                    
                    index += 1
                    y_list = [branch.pathlength_to_soma - branch.pathlength, branch.pathlength_to_soma]
                    plt.plot(x_list, y_list, color = 'black')   
                    plt.plot(cross_bar_x, cross_bar_y, color = 'black')
                elif len(ind_matches) == 1:
                    x_ind = np.mean([dendro_x[ii] for ii in ind_matches])
                    dendro_x.append(x_ind)
                    dendro_y.append(branch.pathlength_to_soma - branch.pathlength)                
                    x_list = [x_ind, x_ind]
                    y_list = [branch.pathlength_to_soma - branch.pathlength, branch.pathlength_to_soma]
                    plt.plot(x_list, y_list, color = 'black')
                elif len(ind_matches) == 0:
                    new_tmp_list.append(branch)
#                else:
#                    print dendro_y, branch.pathlength_to_soma, ind_matches
        tmp_list = new_tmp_list
#        print new_tmp_list
        if len(tmp_list) == 0:
            tmp_flag = 1        
    plt.show()

def plot_sholl(branch_list):
    plt.figure()

    sholl_bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

    sholl_counts = [0]*len(sholl_bins)
    normed_counts = [0]*len(sholl_bins)

    max_euclidean = 0
    for i in range(0, len(branch_list)):
        branch = branch_list[i]
        branch_end_euclidean = branch.euclidean_from_soma + branch.euclidean
        for j in range(0, len(sholl_bins)):
            if branch.euclidean_from_soma < sholl_bins[j] and branch_end_euclidean > sholl_bins[j]:
                sholl_counts[j] = sholl_counts[j] + 1
        if branch_end_euclidean > max_euclidean:
            max_euclidean = branch_end_euclidean

    normed_bins = np.linspace(0, max_euclidean, len(sholl_bins))

    for i in range(0, len(branch_list)):
        branch = branch_list[i]
        branch_end_euclidean = branch.euclidean_from_soma + branch.euclidean
        for j in range(0, len(sholl_bins)):
            if branch.euclidean_from_soma < normed_bins[j] and branch_end_euclidean > normed_bins[j]:
                normed_counts[j] = normed_counts[j] + 1
 
#    plt.plot(sholl_bins, sholl_counts, label='Sholl')
    plt.plot(normed_bins, normed_counts, label='Normed')

    plt.show()
    
def plot_branches(branch_list):
    
    xy_flag = 1
    xz_flag = 1
    branch_order_zoom = 20
    
    if xy_flag == 1:
    
        fig = plt.figure()
        ax = plt.axes(projection='3d')
    
        x_list = []
        y_list = []
        z_list = []
    
        color = next(ax._get_lines.prop_cycler)['color']
    
    #    for branch in branch_list:
    #        if branch.end_compartment.section_type != 1:
    #            for compartment in branch.compartment_list:
    #                x_list.append(compartment.end_coordinates[0])
    #                y_list.append(compartment.end_coordinates[1])
    #                z_list.append(compartment.end_coordinates[2])
    #        plt.plot(x_list, z_list, y_list, color = 'black')
    #        x_list = []
    #        y_list = []
    #        z_list = []
            
        for branch in branch_list:
            if branch.end_compartment.section_type != 1 and branch.branch_order <= branch_order_zoom:
                x_list.append(branch.start_compartment.end_coordinates[0])
                y_list.append(branch.start_compartment.end_coordinates[1])
                z_list.append(branch.start_compartment.end_coordinates[2])
                x_list.append(branch.end_compartment.end_coordinates[0])
                y_list.append(branch.end_compartment.end_coordinates[1])
                z_list.append(branch.end_compartment.end_coordinates[2])
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
        
        ax.set_aspect('equal')
        
        ax.view_init(elev=0, azim=-90)
        plt.show()
    
    if xz_flag == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x_list = []
        y_list = []
        z_list = []
        
        
        plt.scatter(0, 0, linewidths=5)
        for branch in branch_list:
            if branch.end_compartment.section_type != 1 and branch.branch_order <= branch_order_zoom:
                x_list.append(branch.start_compartment.end_coordinates[0])
                y_list.append(branch.start_compartment.end_coordinates[1])
                z_list.append(branch.start_compartment.end_coordinates[2])
                x_list.append(branch.end_compartment.end_coordinates[0])
                y_list.append(branch.end_compartment.end_coordinates[1])
                z_list.append(branch.end_compartment.end_coordinates[2])
    
            plt.plot(x_list, z_list, color = 'black')
            x_list = []
            y_list = []
            z_list = []
        ax.set_aspect('equal')    
        plt.show()
        
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
    all_terminations = [compartment for compartment in all_compartments if len(compartment.daughter_compartment) == 0]

    stem_count = 0
    all_branches = divide_into_branches(all_compartments, all_bifurcations, breakpoints)
    for branch in all_branches:
        if branch.start_coordinates == (0,0,0):
            stem_count = stem_count + 1
    for branch in all_branches:
        if branch.start_coordinates == (0,0,0):
            if stem_count > 1:
                branch.branch_order = 1
            branch.set_branch_order()

#    calculate_remote_angles(all_bifurcations)
#    calculate_local_angles(all_bifurcations)
    calculate_absolute_angles(all_bifurcations)

    return all_compartments, all_branches, all_bifurcations, all_terminations


#file_list = glob.glob('/home/zzchou/Documents/Data/neuron_nmo/claiborne/CNG version/*.swc')
file_list = glob.glob('../Data/Morphologies/claiborne/CNG version/*124894*.swc')
#file_list = glob.glob('../Data/Morphologies/generator/*recreated*.swc')
#file_list = glob.glob('../Data/Morphologies/l_neuron/*recreated*.swc')

start_time = time.time()

for i in range (0, len(file_list)):
    elapsed_time = time.time() - start_time
    print '{} {}% complete'.format(file_list[i], round(100*(i+1)/len(file_list),2))
    compartment_list, branch_list, bifurcation_list, termination_list = read_file(file_list[i])
    
    if i >= 0 and i < 50:
#        avg_data.append(plot_compartment_feature(compartment_list))
        plot_dendrogram(branch_list)
        aa = 0
#        plot_branches(branch_list)
#        plot_star(branch_list)
#        plot_sholl(branch_list)
  
    reconstruct_flag = 0
    if reconstruct_flag == 1:
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        branch_order_zoom = 2
        x_list = []
        y_list = []
        z_list = []
    
        color = next(ax._get_lines.prop_cycler)['color']
    
    #    for branch in branch_list:
    #        if branch.end_compartment.section_type != 1:
    #            for compartment in branch.compartment_list:
    #                x_list.append(compartment.end_coordinates[0])
    #                y_list.append(compartment.end_coordinates[1])
    #                z_list.append(compartment.end_coordinates[2])
    #        plt.plot(x_list, z_list, y_list, color = 'black')
    #        x_list = []
    #        y_list = []
    #        z_list = []
            
        for branch in branch_list:
            if branch.end_compartment.section_type != 1 and branch.branch_order <= branch_order_zoom:
                x_list.append(branch.start_compartment.end_coordinates[0])
                y_list.append(branch.start_compartment.end_coordinates[1])
                z_list.append(branch.start_compartment.end_coordinates[2])
                x_list.append(branch.end_compartment.end_coordinates[0])
                y_list.append(branch.end_compartment.end_coordinates[1])
                z_list.append(branch.end_compartment.end_coordinates[2])
            plt.plot(x_list, z_list, y_list, color = 'black')
            x_list = []
            y_list = []
            z_list = []
    
        ax.set_xlim([-50,50])
        ax.set_ylim([-50,50])
        ax.set_zlim([0,100])
        ax.grid(False)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        
#        ax.set_aspect('equal')
        
        ax.view_init(elev=0, azim=-90)
        
#        plt.figure()
#        ax = plt.axes(projection='3d')
        
        for bifurcation in bifurcation_list:
            
            if bifurcation.parent_branch.branch_order <= branch_order_zoom -1:
        #        print bifurcation.bif_amp_remote, bifurcation.bif_tilt_remote, bifurcation.bif_torque_remote, bifurcation.bif_twist_remote
        #        print bifurcation.daughter_branch
        #        for daughter in bifurcation.daughter_branch:
        #            print daughter.euclidean
                
        
                
                for daughter in bifurcation.daughter_branch:
        
                    x_rot = bifurcation.bif_tilt_remote
                    y_rot = bifurcation.bif_torque_remote
                    z_rot = bifurcation.bif_twist_remote
                    bif_amp = bifurcation.bif_amp_remote*math.pi/180
#                    print x_rot, y_rot, z_rot
#                    x_rot = 0
#                    y_rot = 0
#                    z_rot = 0           
                    if daughter == bifurcation.daughter_branch[0]:
                        x_proj = math.sin(bif_amp/2)*daughter.euclidean
                        y_proj = math.cos(bif_amp/2)*daughter.euclidean
                        z_proj = 0
                        daughter_vector = [x_proj, y_proj, z_proj]
                    else:
                        x_proj = -math.sin(bif_amp/2)*daughter.euclidean
                        y_proj = math.cos(bif_amp/2)*daughter.euclidean
                        z_proj = 0
                        daughter_vector = [x_proj, y_proj, z_proj]
                        
                    new_vector, rotated_axis = angle_noise([0,0,0], daughter_vector, [0,1,0], [0,0,1], x_rot, y_rot, z_rot)
                    start_coordinates = bifurcation.bifurcating_compartment.end_coordinates
                    end_coordinates = [start_coordinates[0] + new_vector[0], start_coordinates[1] + new_vector[1], start_coordinates[2] + new_vector[2]]
                    
                    x_list = []
                    y_list = []
                    z_list = []
                    
#                    x_list.append(start_coordinates[0])
#                    y_list.append(start_coordinates[1])
#                    z_list.append(start_coordinates[2])
#                    x_list.append(end_coordinates[0])
#                    y_list.append(end_coordinates[1])
#                    z_list.append(end_coordinates[2])
                    
                    parent_vector = bifurcation.parent_branch.vector
                    parent_array = np.asarray(parent_vector)
                    parent_mag = np.sqrt(parent_array.dot(parent_array))
                    extend_coord = [0]*3
                    branch = daughter
        #            print branch.vector, branch.rot_x, branch.rot_y, branch.rot_z
                    for i in range(0,3):
                        extend_coord[i] = branch.start_compartment.end_coordinates[i] + branch.pathlength*parent_vector[i]/parent_mag
                    
#                    extend_coord = [0, branch.start_compartment.end_coordinates[1] + branch.pathlength, 0]
                    
                    new_coord = angle_noise_old(branch.start_compartment.end_coordinates, extend_coord, x_rot, y_rot, z_rot)

                    origin_array = np.asarray(branch.start_compartment.end_coordinates)
                    end_array = np.asarray(extend_coord)
                    
                    vector = end_array - origin_array
                    mag = np.sqrt(vector.dot(vector))
                    print vector, mag
#                    print branch.start_compartment.end_coordinates, extend_coord, x_rot, y_rot, z_rot
                    
#                    print x_rot, y_rot, z_rot
                    coordinate_list = [branch.start_compartment.end_coordinates, new_coord]
#                    coordinate_list = [branch.start_compartment.end_coordinates, end_coordinates]

#                    print x_rot, y_rot, z_rot
                    for coord in coordinate_list:
                        x_list.append(coord[0])
                        y_list.append(coord[1])
                        z_list.append(coord[2])
            
#            plt.plot(x_list, z_list, y_list, color = 'blue')
                    plt.plot(x_list, z_list, y_list)
            
#        ax.grid(False)
#        ax.xaxis.pane.fill = False
#        ax.yaxis.pane.fill = False
#        ax.zaxis.pane.fill = False
#        
#        ax.set_aspect('equal')
#        
#        ax.view_init(elev=90, azim=-90)
        plt.show()
            
    