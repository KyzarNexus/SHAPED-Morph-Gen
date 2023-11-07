from __future__ import division

import math
import numpy as np
import morphology_classes as mc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

daughter1 = [1,1,-2]
daughter_mag1 = mc.euclidean_distance([0,0,0], daughter1)
daughter_unit1 = mc.unit_vector(daughter1)
daughter2 = [-1,1,2]
daughter_mag2 = mc.euclidean_distance([0,0,0], daughter2)
daughter_unit2 = mc.unit_vector(daughter2)



bif_midline = [(daughter_unit1[i] + daughter_unit2[i])/2 for i in range(0,3)]
bif_midline = mc.unit_vector(bif_midline)

bif_amp = math.acos(mc.dot_product(daughter_unit1, daughter_unit2)/mc.magnitude_product(daughter_unit1, daughter_unit2))
norm = mc.cross_product(daughter_unit1, daughter_unit2)
norm = mc.unit_vector(norm)
if norm[1] < 0:
    norm = mc.cross_product(daughter_unit2, daughter_unit1)
    norm = mc.unit_vector(norm)




u = [1,1,0]
theta, phi, psi = mc.rotation_angles(bif_midline, norm)

print 'daughter1 = ', daughter_unit1
print 'bif_midline = ', bif_midline

if daughter_unit1[2]*bif_midline[2] > 0:
    if abs(daughter_unit1[2])<abs(bif_midline[2]):
        psi = psi + math.pi
elif bif_midline[2] == 0:
    if daughter_unit1[2] < 0:
        psi = psi + math.pi
else:
    psi = psi + math.pi
    
if daughter_unit1[0]*bif_midline[0] > 0:
    if abs(daughter_unit1[0])<abs(bif_midline[0]):
        psi = psi + math.pi
elif bif_midline[0] == 0:
    if daughter_unit1[0] < 0:
        psi = psi + math.pi
else:
    psi = psi + math.pi
    
theta = theta*180/math.pi
phi = phi*180/math.pi
psi = psi*180/math.pi

origin = [0,0,0]
daughter_mag = mc.euclidean_distance([0,0,0], daughter1)
half_angle = bif_amp*180/(math.pi*2)
print 'inputs:', origin, daughter_mag, half_angle, theta, phi, psi
new_end = mc.angle_noise_old(origin, daughter_mag, half_angle, theta, phi, psi)


print 'coordinates = ', new_end
