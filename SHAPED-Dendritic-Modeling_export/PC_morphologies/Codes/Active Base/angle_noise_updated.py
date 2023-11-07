from __future__ import division

import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def euclidean_distance(u,v):
    ans = math.sqrt((u[0] - v[0])**2 + (u[1] - v[1])**2 + (u[2] - v[2])**2)
    return ans

def cross_product(u,v):
    s1 = u[1]*v[2] - u[2]*v[1]
    s2 = u[2]*v[0] - u[0]*v[2]
    s3 = u[0]*v[1] - u[1]*v[0]
    return (s1, s2, s3)

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

def angle_noise_old(origin, end_point, theta, phi, psi):
    
    origin_array = np.asarray(origin)
    end_array = np.asarray(end_point)
    
    vector = end_array - origin_array
    mag = np.sqrt(vector.dot(vector))
    
    vector = [0, 1, 0]
    
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
    
    v_array = np.asarray(vector)

    new_vector = np.ndarray.tolist(np.dot(v_array, R))[0]
    # print origin, v_array, mag
    
    new_end = [origin[i] + new_vector[i]*mag for i in range(0, 3)]
    return new_end    

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
    
def rotation_angles_old(u, norm):
    # n1 = (1, 0, 0)
    # n2 = (0, 1, 0)

    # w1 = [dot_product(u,n1)/magnitude_product(n1,n1), 0, 0]
    # w2 = [0, dot_product(u,n2)/magnitude_product(n2,n2), 0]

    # proj_yz = [u[i] - w1[i] for i in range(0,3)]
    # proj_xz = [u[i] - w2[i] for i in range(0,3)]

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
    # print new_u
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
    # n_norm = [0, 0, 1]
    
    rot_z = math.acos(dot_product(n2, norm)/magnitude_product(n2, norm))
    if new_norm[0] < 0:
        rot_z = -rot_z
    # Rr_z = np.matrix( ((math.cos(-rot_z), 0, -math.sin(-rot_z)), (math.sin(-rot_z), 0, math.cos(-rot_z)), (0, 0, 1)) )

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
 