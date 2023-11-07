from __future__ import division

import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
        
        self.branch_flag_list = []

        
        
        max_y = 0
        
        for compartment in self.compartment_list:
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
        
#        print self.max_euclidean
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
            bifurcation.bif_normal = (0,0,1)
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

            bifurcation.bif_tilt_remote, bifurcation.bif_torque_remote, bifurcation.bif_twist_remote = rotation_angles(bif_midline2, bif_norm2)

            bifurcation.bif_tilt_remote = bifurcation.bif_tilt_remote * 180/math.pi
            bifurcation.bif_torque_remote = bifurcation.bif_torque_remote * 180/math.pi
            bifurcation.bif_twist_remote = bifurcation.bif_twist_remote * 180/math.pi
            
#            bifurcation.bif_tilt_remote = rot_x_diff * 180/math.pi
#            bifurcation.bif_torque_remote = rot_y_diff * 180/math.pi
#            bifurcation.bif_twist_remote = rot_z_diff * 180/math.pi
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
#        rot_x = math.acos(dot_product(proj_yz, n)/magnitude_product(proj_yz, n))
        rot_x = math.acos(dot_product(u, n)/magnitude_product(u, n))
        
#    if proj_yz[2] > 0:
#        rot_x = -rot_x
    
    Rr_x = np.matrix( ((1, 0, 0), (0, math.cos(-rot_x), -math.sin(-rot_x)), (0, math.sin(-rot_x), math.cos(-rot_x))) )
    Rn_x = np.matrix( ((1, 0, 0), (0, math.cos(rot_x), -math.sin(rot_x)), (0, math.sin(rot_x), math.cos(rot_x))) )    
#    new_u = np.ndarray.tolist(np.dot(u, Rr_x))[0]
    
    proj_xz = [u[0], 0, u[2]]
    n2 = [0,0,1]
    if proj_xz == [0,0,0]:
        rot_y = 0
    else:
        rot_y = math.acos(dot_product(proj_xz, n2)/magnitude_product(proj_xz, n2))
    if u[0] < 0:
        rot_y = -rot_y
    Rr_y = np.matrix( ((math.cos(-rot_y), 0, -math.sin(-rot_y)), (0, 1, 0), (math.sin(-rot_y), 0, math.cos(-rot_y))) )
    Rn_y = np.matrix( ((math.cos(rot_y), 0, -math.sin(rot_y)), (0, 1, 0), (math.sin(rot_y), 0, math.cos(rot_y))) )
    
#    if new_u == [0,0,0]:
#        rot_z = 0
#    else:
#        rot_z = math.acos(dot_product(new_u, n)/magnitude_product(new_u, n))
#    if new_u[0] < 0:
#        rot_z = -rot_z
#    
#    Rr_z = np.matrix( ((math.cos(-rot_z), -math.sin(-rot_z), 0), (math.sin(-rot_z), math.cos(-rot_z), 0), (0, 0, 1)) )
    
    Rr_yx = Rr_y*Rn_x
    Rn_xy = Rr_x*Rn_y
    Rn_2 = Rr_y*Rr_x
    new_u = np.ndarray.tolist(np.dot(u, Rr_yx))[0]
    tmp_u = np.dot(u, Rr_y)
    tmp_u2 = np.ndarray.tolist(np.dot(tmp_u, Rn_x))[0]
    print tmp_u2
#    proj_xz2 = [new_u[0], 0, new_u[2]]
#    n3 = [1,0,0]
#    if proj_xz2 == [0,0,0]:
#        rot_y2 = 0
#    else:
#        rot_y2 = math.acos(dot_product(proj_xz, n3)/magnitude_product(proj_xz2, n3))
#    if new_u[2] < 0:
#        rot_y2 = -rot_y2
    if norm == (0, 0, 0):
        norm = (0, 0, 1)
    
#    new_norm = np.ndarray.tolist(np.dot(norm, Rn_xy))[0]
    new_norm = np.ndarray.tolist(np.dot(norm, Rr_yx))[0]
#    new_norm[1] = 0

    n_norm = [0, 0, -1]

    print 'norm=', norm
    print 'new_norm=', new_norm

    if magnitude_product(new_norm, n_norm) == 0:
        rot_z = 0
    else:
        rot_z = math.acos(dot_product(new_norm,n_norm)/magnitude_product(new_norm, n_norm))
#    if new_norm[0] < 0:
#        rot_z = -rot_z
    Rr_z = np.matrix( ((math.cos(-rot_z), -math.sin(-rot_z), 0), (math.sin(-rot_z), math.cos(-rot_z), 0), (0, 0, 1)) )

#    if norm[2] <= 0 and u[2] <= 0:
#        rot_z = -rot_z
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

def unit_vector((u1, u2, u3)):
    magnitude = euclidean_distance((u1, u2, u3), (0, 0, 0))
    
    if magnitude != 0:  
        return [u1/magnitude, u2/magnitude, u3/magnitude]
    else:
        return [0, 0, 0]
    
def angle_noise_old(origin, daughter_mag, half_angle, theta, phi, psi):
    half_angle = half_angle*math.pi/180
    if half_angle >= 0:
        new_end2 = [math.sin(half_angle)*daughter_mag, abs(math.cos(half_angle))*daughter_mag, 0]
    else:
        new_end2 = [math.sin(half_angle)*daughter_mag, abs(math.cos(half_angle))*daughter_mag, 0]
    end_point2 = [new_end2[ii] + origin[ii] for ii in range(3)]
    
    end_point = [0, daughter_mag, 0]
    origin_array = np.asarray(origin)
    end_array = np.asarray(end_point)
    
    vector = end_array
    mag = daughter_mag
    
    vector2 = new_end2
    
    vector_unit = vector/mag
    vector_list = np.ndarray.tolist(vector_unit)
    
    ux = vector_list[0]
    uy = vector_list[1]
    uz = vector_list[2]
    
    R_pre = np.matrix(( (0, -1, 0), (1, 0, 0), (0, 0, 1) ))
    
    perp_vector_list = np.ndarray.tolist(np.dot(vector_list, R_pre))[0]
    
    vx = perp_vector_list[0]
    vy = perp_vector_list[1]
    vz = perp_vector_list[2]

    R_pre2 = np.matrix(( (1, 0, 0), (0, 0, -1), (0, 1, 0) ))
    
    perp_vector_list2 = np.ndarray.tolist(np.dot(vector_list, R_pre2))[0]
    
    wx = perp_vector_list2[0]
    wy = perp_vector_list2[1]
    wz = perp_vector_list2[2]




    rot_x = theta*math.pi/180
    rot_y = phi*math.pi/180
    rot_z = psi*math.pi/180

#    s = math.sin(rot_x)
#    c = math.cos(rot_x)
#    omc = 1-c

#    R_x = np.matrix(( (math.cos(rot_x) + vx*vx*(1-math.cos(rot_x)), vx*vy*(1-math.cos(rot_x))-vz*math.sin(rot_x), vx*vz*(1-math.cos(rot_x))+vy*math.sin(rot_x) ),
#                      (vy*vx*(1-math.cos(rot_x))+vz*math.sin(rot_x), math.cos(rot_x)+vy*vy*(1-math.cos(rot_x)), vy*vz*(1-math.cos(rot_x))-vx*math.sin(rot_x)),
#                      (vz*vx*(1-math.cos(rot_x))-vy*math.sin(rot_x), vz*vy*(1-math.cos(rot_x))+vx*math.sin(rot_x), math.cos(rot_x)+vz*vz*(1-math.cos(rot_x))) ))
#    
#    
#    R_y = np.matrix( ((math.cos(rot_y) + ux*ux*(1-math.cos(rot_y)), ux*uy*(1-math.cos(rot_y))-uz*math.sin(rot_y), ux*uz*(1-math.cos(rot_y))+uy*math.sin(rot_y)),
#                      (uy*ux*(1-math.cos(rot_y))+uz*math.sin(rot_y), math.cos(rot_y)+uy*uy*(1-math.cos(rot_y)), uy*uz*(1-math.cos(rot_y))-ux*math.sin(rot_y)),
#                      (uz*ux*(1-math.cos(rot_y))-uy*math.sin(rot_y), uz*uy*(1-math.cos(rot_y))+ux*math.sin(rot_y), math.cos(rot_y)+uz*uz*(1-math.cos(rot_y)))) )
#   
#    
#    R_z = np.matrix( ((math.cos(rot_z) + wx*wx*(1-math.cos(rot_z)), wx*wy*(1-math.cos(rot_z))-wz*math.sin(rot_z), wx*wz*(1-math.cos(rot_z))+wy*math.sin(rot_z)),
#                      (wy*wx*(1-math.cos(rot_z))+wz*math.sin(rot_z), math.cos(rot_z)+wy*wy*(1-math.cos(rot_z)), wy*wz*(1-math.cos(rot_z))-wx*math.sin(rot_z)),
#                      (wz*wx*(1-math.cos(rot_z))-wy*math.sin(rot_z), wz*wy*(1-math.cos(rot_z))+wx*math.sin(rot_z), math.cos(rot_z)+wz*wz*(1-math.cos(rot_z)))) )
    
    
    Rp_x = np.matrix(( (1, 0, 0 ),
                      (0, math.cos(rot_x), -math.sin(rot_x)),
                      (0, math.sin(rot_x), math.cos(rot_x)) ))
    Rn_x = np.matrix(( (1, 0, 0 ),
                      (0, math.cos(-rot_x), -math.sin(-rot_x)),
                      (0, math.sin(-rot_x), math.cos(-rot_x)) ))
    
    Rp_y = np.matrix(( (math.cos(rot_y), 0, math.sin(rot_y) ),
                      (0, 1, 0),
                      (-math.sin(rot_y), 0, math.cos(rot_y)) ))
    Rn_y = np.matrix(( (math.cos(-rot_y), 0, math.sin(-rot_y) ),
                      (0, 1, 0),
                      (-math.sin(-rot_y), 0, math.cos(-rot_y)) ))
   
    rot_z = -rot_z
    R_z = np.matrix(( (math.cos(rot_z), -math.sin(rot_z), 0 ),
                      (math.sin(rot_z), math.cos(rot_z), 0),
                      (0, 0, 1) ))

#    print R_x
#    print R_y
#    print R_z

    R = R_z * Rn_y * Rn_x
    R4 = Rn_x * Rn_y * R_z
    R2 = Rn_y*Rn_x
    R3 = Rn_x*Rn_y

    
    v_array = np.asarray(vector)
    w_array = np.dot(v_array, R3)
    
    w = np.ndarray.tolist(w_array)[0]

    w = unit_vector(w)
    wx = w[0]
    wy = w[1]
    wz = w[2]

    Rp_z = np.matrix( ((math.cos(rot_z) + wx*wx*(1-math.cos(rot_z)), wx*wy*(1-math.cos(rot_z))-wz*math.sin(rot_z), wx*wz*(1-math.cos(rot_z))+wy*math.sin(rot_z)),
                      (wy*wx*(1-math.cos(rot_z))+wz*math.sin(rot_z), math.cos(rot_z)+wy*wy*(1-math.cos(rot_z)), wy*wz*(1-math.cos(rot_z))-wx*math.sin(rot_z)),
                      (wz*wx*(1-math.cos(rot_z))-wy*math.sin(rot_z), wz*wy*(1-math.cos(rot_z))+wx*math.sin(rot_z), math.cos(rot_z)+wz*wz*(1-math.cos(rot_z)))) )

    Rn_z = np.matrix( ((math.cos(-rot_z) + wx*wx*(1-math.cos(-rot_z)), wx*wy*(1-math.cos(-rot_z))-wz*math.sin(-rot_z), wx*wz*(1-math.cos(-rot_z))+wy*math.sin(-rot_z)),
                      (wy*wx*(1-math.cos(-rot_z))+wz*math.sin(-rot_z), math.cos(-rot_z)+wy*wy*(1-math.cos(-rot_z)), wy*wz*(1-math.cos(-rot_z))-wx*math.sin(-rot_z)),
                      (wz*wx*(1-math.cos(-rot_z))-wy*math.sin(-rot_z), wz*wy*(1-math.cos(-rot_z))+wx*math.sin(-rot_z), math.cos(-rot_z)+wz*wz*(1-math.cos(-rot_z)))) )

    rot_z = -rot_z
    R_y2 = np.matrix(( (math.cos(rot_z), 0, math.sin(rot_z) ),
                      (0, 1, 0),
                      (-math.sin(rot_z), 0, math.cos(rot_z)) ))

    R5 = R_y2 * Rn_x * Rn_y

#    rot_x = -rot_x
#    R_x = np.matrix(( (1, 0, 0 ),
#                      (0, math.cos(rot_x), -math.sin(rot_x)),
#                      (0, math.sin(rot_x), math.cos(rot_x)) ))
    
    vector2 = unit_vector(vector2)
    v2_array = np.asarray(vector2)
#    u_array = np.dot(v2_array, R3)
    u_array = np.dot(v2_array, R_y2)
#    u2_array = np.dot(u_array, Rn_z)
    u2_array = np.dot(u_array, R3)
#    u_array = np.asarray([0,daughter_mag, 0])
#    new_vector = np.ndarray.tolist(np.dot(u_array, R3))[0]
    new_vector = np.ndarray.tolist(u2_array)[0]
    new_vector = [x*daughter_mag for x in new_vector]
    
    print 'v=', vector
    print 'v2=', vector2
    print 'w=', w
    print 'u=', u_array
    print 'u2=', u2_array
    
    new_end = [origin[i] + new_vector[i] for i in range(0, 3)]
    return new_end

def plot_bifurcations(bifurcation_list):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for bifurcation in bifurcation_list:
        for daughter in bifurcation.daughter_branch:
            compartment = bifurcation.bifurcating_compartment
            x_list.append(compartment.end_coordinates[0])
            y_list.append(compartment.end_coordinates[1])
            z_list.append(compartment.end_coordinates[2])
            
            daughter_mag = daughter.euclidean
            if daughter == bifurcation.daughter_branch[0]:
                half_angle = bifurcation.bif_amp_remote/2
            else:
                half_angle = -bifurcation.bif_amp_remote/2
            theta = bifurcation.bif_tilt_remote
            phi = bifurcation.bif_torque_remote
            psi = bifurcation.bif_twist_remote
            
#            print theta, phi, psi
            
            rotated_end = angle_noise_old(compartment.end_coordinates, daughter_mag, half_angle, theta, phi, 0)
#            print new_end, rotated_end
            x_list.append(rotated_end[0])
            y_list.append(rotated_end[1])
            z_list.append(rotated_end[2])        
            plt.plot(x_list, z_list, y_list, color = 'black')
            x_list = []
            y_list = []
            z_list = []

    ax.set_xlim([-200,200])
    ax.set_ylim([-200,200])
    ax.set_zlim([0,400])
    plt.axis('off')
    plt.grid(None)
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.view_init(elev=0, azim=90)
    plt.show()
    
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
    plt.grid(None)
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.view_init(elev=0, azim=90)
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
#                print breakpoints
    all_terminations = [compartment for compartment in all_compartments if len(compartment.daughter_compartment) == 0]

    all_branches = divide_into_branches(all_compartments, all_bifurcations, breakpoints)
    
#    for branch in all_branches:
#        print branch.branch_order, branch.pathlength, branch.pathlength_to_soma

    calculate_remote_angles(all_bifurcations)
    calculate_local_angles(all_bifurcations)

    return all_compartments, all_branches, all_bifurcations, all_terminations