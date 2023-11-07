import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np

global c_id
c_id = 1

n_stem_min = 1
n_stem_max = 1
daughter_ratio_min = 1
daughter_ratio_max = 1
bif_twist_min = -180
bif_twist_max = 180

max_branch_order = 8
jitter_deviation = 0.3
weight = 0.5

tropism_mean = 0
tropism_sd = 0
neuron_pathlength = random.gauss(400, 20)

num_morphologies = 50

# Emergent screening
total_dendritic_length_range = (2000, 6000)
number_of_bifurcations_range = (10, 25)
total_surface_area_range = (0, 12000)
major_axis_range = (200, 550)
minor_axis_range = (100, 300)

# Gaussian Distribution Constants
soma_diameter_mean = 2
soma_diameter_sd = 0

stem_diameter_mean = 9
stem_diameter_sd = 2

taper_2_mean = 0.1
taper_2_sd = 0

n_stem_coef = (1,2)
soma_diameter_coef = (2, 0)

stem_pathlength_coef = (30, 10)
ibf_pathlength_coef = (0.98, 2.99, 71.67)
terminal_pathlength_coef = (3.93, 138.1)

ibf_pathlength_1_coef = (0.98, 2.99, 71.67)
ibf_pathlength_2_coef = (1.95, 0.56, 45.9)
ibf_pathlength_3_coef = (1.2, 6.81, 59.33)
ibf_pathlength_4_coef = (0.16, 13.86, 21.1)

terminal_pathlength_1_coef = (333.6, -1006.2, 3.776)
terminal_pathlength_2_coef = (43.48, -238.77, 10.50)
terminal_pathlength_3_coef = (3.82, -6.76, 32.35)
terminal_pathlength_4_coef = (2.426, -0.42, 29.4)


stem_diameter_coef = (5.101, 0.7)
ibf_diameter_coef = (1.446, 0.672)
terminal_diameter_coef = (0.538, 0.196)

ibf_taper_coef = (7.978, -0.258, 0.057)
terminal_taper_coef = (260.844, -2.521, 0.011)

daughter_ratio_coef = (1.150, 0.194)
parent_daughter_ratio_coef = (1.041, 0.063)

bif_amplitude_remote_coef = (3.51, 0.17, 9.44)
bif_elevation_remote_coef = (0, 18.4)
bif_azimuth_remote_coef = (-180, 180)
bif_orientation_remote_coef = (0, 14.54)

diam_threshold_coef = (0.441, 0.2)

# bif_amplitude_remote_coef = (22.727, 20)


bif_amplitude_local_coef = (5.51, -0.83, 9.275)
bif_elevation_local_coef = (0, 18.4)
bif_azimuth_local_coef = (-180, 180)
bif_orientation_local_coef = (0, 13.232)

# bif_elevation_coef = (0, 0)
# bif_azimuth_coef = (0, 0)
# bif_orientation_coef = (0, 0)

ibf_branch_pathlength_mean = 80
ibf_branch_pathlength_sd = 20
diam_threshold_mean = 0.15
diam_threshold_sd = 0.1
term_branch_pathlength_mean = 150
term_branch_pathlength_sd = 20
bif_amplitude_remote_mean = 50
bif_amplitude_remote_sd = 10
bif_amplitude_local_mean = 50
bif_amplitude_local_sd = 10
bif_elevation_mean = 0
bif_elevation_sd = 0
bif_azimuth_min = 0
bif_azimuth_max = 0
bif_orientation_mean = 0
bif_orientation_sd = 0



# bif_tilt_local_mean = 0
# bif_tilt_local_sd = 10
# bif_tilt_remote_mean = 0
# bif_tilt_remote_sd = 5
# bif_torque_local_mean = 0
# bif_torque_local_sd = 10
# bif_torque_remote_mean = 0
# bif_torque_remote_sd = 5
# bif_twist_mean = 90
# bif_twist_sd = 30
# elevation_mean = math.pi/2
# elevation_sd = 0.1

fragmentation_mean = 17
fragmentation_sd = 0

view1 = 20
view2 = -135

def cross_product(u, v):
    s1 = u[1] * v[2] - u[2] * v[1]
    s2 = u[2] * v[0] - u[0] * v[2]
    s3 = u[0] * v[1] - u[1] * v[0]
    return [s1, s2, s3]

def dot_product(u, v):
    s = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    return s

def magnitude_product(u, v):
    s1 = math.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    s2 = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return s1*s2

def euclidean_distance(u,v):
    ans = math.sqrt((u[0] - v[0])**2 + (u[1] - v[1])**2 + (u[2] - v[2])**2)
    return ans

class Compartment:
    def __init__(self, origin, end_coordinates, diameter):
        self.start_coordinates = origin
        self.vector = [end_coordinates[i] - origin[i] for i in range(0, 3)]
        self.parent_compartment = []
        self.end_coordinates = end_coordinates
        self.jitter_phi = 0
        self.jitter_theta = 0
        global c_id
        self.compartment_id = c_id
        c_id = c_id + 1
        
        self.diameter = diameter
        self.terminal_flag = 0
        self.bifurcation_flag = 0
        self.stem_flag = 0
        self.length = euclidean_distance(self.start_coordinates, self.end_coordinates)

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

def draw_value(request):
    parameters = {

        'n_stem': (n_stem_coef, 'uniform'),
        
        'soma_diameter': (soma_diameter_coef, 'gaussian'),
        'stem_diameter': (stem_diameter_coef, 'gaussian'),
        'ibf_diameter': (ibf_diameter_coef, 'gaussian'),
        'terminal_diameter': (terminal_diameter_coef, 'gaussian'),
        
        'ibf_taper': (ibf_taper_coef, 'gamma'),
        'terminal_taper': (terminal_taper_coef, 'gamma'),
        
        'daughter_ratio': (daughter_ratio_coef, 'gaussian'),
        'parent_daughter_ratio': (parent_daughter_ratio_coef, 'gaussian'),
        
        'diam_threshold': (diam_threshold_coef, 'gaussian'),
        
        'fragmentation': ((fragmentation_mean, fragmentation_sd), 'gaussian'),
        
        'ibf_pathlength_1': (ibf_pathlength_1_coef, 'gamma'),
        'ibf_pathlength_2': (ibf_pathlength_2_coef, 'gamma'),        
        'ibf_pathlength_3': (ibf_pathlength_3_coef, 'gamma'),
        'ibf_pathlength_4': (ibf_pathlength_4_coef, 'gamma'),
        
        'terminal_pathlength_1': (terminal_pathlength_1_coef, 'gamma'),
        'terminal_pathlength_2': (terminal_pathlength_2_coef, 'gamma'),        
        'terminal_pathlength_3': (terminal_pathlength_3_coef, 'gamma'),        
        'terminal_pathlength_4': (terminal_pathlength_4_coef, 'gamma'),
        
        'stem_pathlength': (stem_pathlength_coef, 'gaussian'),
        'terminal_pathlength': (terminal_pathlength_coef, 'exponential'),
        
        'contraction': ((0.9, 0), 'gaussian'),

        'bif_amplitude_remote': (bif_amplitude_remote_coef, 'gamma'),
        'bif_elevation_remote': (bif_elevation_remote_coef, 'gaussian'),
        'bif_azimuth_remote': (bif_azimuth_remote_coef, 'uniform'),
        'bif_orientation_remote': (bif_orientation_remote_coef, 'gaussian'),

        'bif_amplitude_local': (bif_amplitude_local_coef, 'gamma'),
        'bif_elevation_local': (bif_elevation_local_coef, 'gaussian'),
        'bif_azimuth_local': (bif_azimuth_local_coef, 'uniform'),
        'bif_orientation_local': (bif_orientation_local_coef, 'gaussian'),
        
        'tropism': ((tropism_mean, tropism_sd), 'gaussian')
    }

    (coef, distribution_type) = parameters[request]

    if distribution_type == 'uniform':
        return random.randint(coef[0], coef[1])
    elif distribution_type == 'gaussian':
        return random.gauss(coef[0], coef[1])
    elif distribution_type == 'gamma':
        return random.gammavariate(coef[0], coef[2]) + coef[1]
    elif distribution_type == 'exponential':
        return random.expovariate(1/coef[1]) + coef[0]
        
    else:
        return random.gauss(coef[0], coef[1])
                            
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

def plot_branches(branch_list):
    # fig = plt.figure()
    ax = plt.axes(projection='3d')

    x_list = []
    y_list = []
    z_list = []

    color = next(ax._get_lines.prop_cycler)['color']

    for branch in branch_list:
        for compartment in branch.compartment_list:
            x_list.append(compartment.end_coordinates[0])
            y_list.append(compartment.end_coordinates[1])
            z_list.append(compartment.end_coordinates[2])
        plt.plot(x_list, z_list, y_list, color=color)
        x_list = []
        y_list = []
        z_list = []

    ax.set_xlim([-150,150])
    ax.set_ylim([-150,150])
    ax.set_zlim([0,300])
    ax.view_init(elev=view1, azim=view2)
    plt.show()

def add_branch_to_list(branch, all_branch_list):
    all_branch_list.append(branch)
    if branch.daughter_bifurcation != []:
        add_branch_to_list(branch.daughter_bifurcation.daughter_branches[0], all_branch_list)
        add_branch_to_list(branch.daughter_bifurcation.daughter_branches[1], all_branch_list)

# # Main
    # print n_stem

    # all_branch_list = []

    # for index in range(0, n_stem):
    #    stem = Stem()
    #    Bifurcation(stem)

    #    add_branch_to_list(stem, all_branch_list)




    # for j in range(0, num_morphologies):
    #    global c_id
    #    global neuron_pathlength_mean
    
    #    redo_flag = 1
    
    

    #    while redo_flag == 1:
    #        n_stem = int(draw_value('n_stem'))
    # #        print n_stem
        
    #        if n_stem == 1:
    #            orientation_vectors = [(0,1,0)]
    #        elif n_stem == 2:
    #            stem_angle = 30*math.pi/180
    #            orientation_vectors = [(math.cos(stem_angle), math.sin(stem_angle), 0), (-math.cos(stem_angle), math.sin(stem_angle), 0)]
    #        elif n_stem == 3:
    #            orientation_vectors = [(0,1,0.5), (0.5,1,-0.3), (-0.5,1,-0.3)]

    #        neuron_pathlength = random.gauss(400, 20)
    #        c_id = 1

    #        all_branch_list = []
    #        all_compartment_list = []
    #        all_terminal_list = []
    #        all_bifurcation_list= []
    
    #        for index in range(0, n_stem):    
    #            if index != 0:
    #                c_id = c_id - 1            
    #            stem = Stem(orientation_vectors[index])

    #            Bifurcation(stem)
                
        
    #            add_branch_to_list(stem, all_branch_list)
            
    #        for branch in all_branch_list:
    #            for compartment in branch.compartment_list:
    #                if compartment != branch.compartment_list[0] or compartment.compartment_id == 1:
    #                    all_compartment_list.append(compartment)
    #                if compartment.terminal_flag == 1:
    #                    all_terminal_list.append(compartment)
    #        for compartment in all_compartment_list:
    #                if compartment.bifurcation_flag == 1:
    #                    all_bifurcation_list.append(compartment)
    #        total_dendritic_length = 0
    #        total_surface_area = 0
    #        normalized_total_surface_area = 0        
    #        major_coordinates = [compartment.end_coordinates[0] for compartment in all_terminal_list]
    #        minor_coordinates = [compartment.end_coordinates[2] for compartment in all_terminal_list]
    #        num_bifurcations = len(all_bifurcation_list)
        
    #        original_coordinates = [np.array((compartment.end_coordinates[0], compartment.end_coordinates[2])) for compartment in all_terminal_list]
        
    #        major_axis = 0
    #        minor_axis = 0
        
    #        for angle in range(0,180):
    #            tmp_angle = angle *math.pi/180
    #            R = np.matrix(((math.cos(tmp_angle), -math.sin(tmp_angle)), (math.sin(tmp_angle), math.cos(tmp_angle))))
            
    #            rotated_coordinates = [np.ndarray.tolist(np.dot(entry, R))[0] for entry in original_coordinates]
            
    #            major_coordinates = [coordinates[0] for coordinates in rotated_coordinates]
    #            minor_coordinates = [coordinates[1] for coordinates in rotated_coordinates]
            
    #            tmp_major_axis = abs(max(major_coordinates) - min(major_coordinates))
    #            tmp_minor_axis = abs(max(minor_coordinates) - min(minor_coordinates))
            
    #            if tmp_major_axis > major_axis:
    #                major_axis = tmp_major_axis
    #                minor_axis = tmp_minor_axis
                
                
    #        for branch in all_branch_list:
    #            total_dendritic_length = total_dendritic_length + branch.pathlength
    
    #        length_scale_factor = total_dendritic_length/100    
    #        for compartment in all_compartment_list:
    #            total_surface_area = total_surface_area + math.pi*compartment.diameter*compartment.length
    #            normalized_total_surface_area = normalized_total_surface_area + math.pi*compartment.diameter*compartment.length/length_scale_factor
    #        for termination in all_terminal_list:
    #            total_surface_area = total_surface_area + math.pi*(0.5*termination.diameter)**2
            
    #        if total_dendritic_length < total_dendritic_length_range[1] and total_dendritic_length > total_dendritic_length_range[0]:
    #            if major_axis < major_axis_range[1] and major_axis > major_axis_range[0]:
    #                if minor_axis < minor_axis_range[1] and minor_axis > minor_axis_range[0]:
    #                    if num_bifurcations < number_of_bifurcations_range[1] and num_bifurcations > number_of_bifurcations_range[0]: 
    #                        redo_flag = 0
    

    # #    print num_bifurcations        
    #    output_title = '/home/zzchou/Documents/Data/test_output/{}{}{}'.format('morphology_output_', j+1, '.swc')
    #    with open(output_title, 'w') as f:
    #        true_origin = all_branch_list[0].compartment_list[0]
    #        line = (1, 1, true_origin.end_coordinates[0], true_origin.end_coordinates[1], true_origin.end_coordinates[2], draw_value('soma_diameter'), -1)
    #        parent_id = 1
    #        strs=' '.join(str(x) for x in line)
    #        f.write(strs + '\n')

    #        for branch in all_branch_list:
    #            for compartment in branch.compartment_list:
    #                if compartment.stem_flag == 1:
    #                    parent_id = 1
    #                if compartment != branch.compartment_list[0]:
                
    #                    line = (compartment.compartment_id, 3, compartment.end_coordinates[0], compartment.end_coordinates[1], compartment.end_coordinates[2], compartment.diameter, parent_id)
    #                    strs=' '.join(str(x) for x in line)
    #                    f.write(strs + '\n')
    #                parent_id = compartment.compartment_id
    #    f.closed

    #    c_id = 1

    #    neuron_pathlength_mean = random.gauss(400, 20)


    # plot_branches(all_branch_list)


