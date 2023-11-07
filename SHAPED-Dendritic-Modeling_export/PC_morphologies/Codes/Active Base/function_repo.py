import os
import pickle
import math
from mpl_toolkits import mplot3d
import numpy as np
from scipy.ndimage.measurements import label 
import matplotlib as mpl
import matplotlib.pyplot as plt
from ast import literal_eval

def vertical_align(source_coords,secondary_coords = 0, plot_flag = False):
    # Purpose: Align a 3d point cloud toward the vertical (Z) axis
    # Input: A 3xN np array of 3D coordinate data
    # Method: Using the average point cloud vector, compute the transformation matrix for the dataset. Then apply the transform to the data.
    # Requirements: numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, mpl_toolkits import mplot3d
    
    xdata = [x[0] for x in source_coords]; xavg = np.mean(xdata)
    ydata = [x[1] for x in source_coords]; yavg = np.mean(ydata)
    zdata = [x[2] for x in source_coords]; zavg = np.mean(zdata)
    avgVect = np.array([0,0,0, xavg,yavg,zavg]) # Origin included for plotting purposes

    # Z Axis Transform To X Axis Alignment
    phi = np.arctan(avgVect[4]/avgVect[3])
    # Rotation Around Z Axis
    rotZ = np.array([
        [np.cos(phi),-np.sin(phi),0],
        [np.sin(phi),np.cos(phi),0],
        [0,0,1]]); rotZ = np.transpose(rotZ) # Array initialized w/o transpose to preserve readability
    rotZvect = np.matmul(rotZ,avgVect[3:]) # Original vector after rotation
    rotZvect = np.append(np.array([0,0,0]),rotZvect)

    # Y Axis Transform to Z Axis Alignment
    theta = np.arctan(rotZvect[3]/rotZvect[5])
    # Rotation Around Y Axis
    rotY = np.array([
        [np.cos(theta),0,np.sin(theta)],
        [0,1,0],
        [-np.sin(theta),0,np.cos(theta)]]); rotY = np.transpose(rotY)
    rotYvect = np.matmul(rotY,avgVect[3:]) # Original vector after rotation
    rotYvect = np.append(np.array([0,0,0]),rotYvect)

    # Computing and executing transform
    transform = np.matmul(rotY,rotZ)
    vert_coords = np.matmul(transform,np.transpose(source_coords)); vert_coords = np.transpose(vert_coords)
    xnew = [x[0] for x in vert_coords]; xavgnew = np.mean(xnew)
    ynew = [x[1] for x in vert_coords]; yavgnew = np.mean(ynew)
    znew = [x[2] for x in vert_coords]; zavgnew = np.mean(znew)
    avgVectNew = np.array([0,0,0, xavgnew,yavgnew,zavgnew]) # Used to verify transform
    invert_flag = False
    if zavgnew < 0: # Inverting if the dataset is upside-down
        invert_flag = True # Used for basal flip
        for coord_index in range(len(vert_coords)):
            vert_coords[coord_index][2] = -vert_coords[coord_index][2] 




    if type(secondary_coords) == np.ndarray:
        second_vert = np.matmul(transform,np.transpose(secondary_coords)); second_vert = np.transpose(second_vert)
        if invert_flag == True:
            for coord_index in range(len(second_vert)):
                second_vert[coord_index][2] = -second_vert[coord_index][2] 


    if plot_flag == True:
        # Axis Vectors
        refvects = np.array([0,0,0, 0,0,1])
        refvects = np.vstack((refvects,[0,0,0, 0,1,0])) 
        refvects = np.vstack((refvects,[0,0,0, 1,0,0])) 

        # Plotting
        fig = plt.figure
        ax = plt.axes(projection='3d')
        ax.set_title('Vertical Transform')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        # Point Clouds

        # if type(secondary_coords) == np.ndarray:
        #     xAll = [x[0] for x in source_coords] + [x[0] for x in secondary_coords]
        #     yAll = [x[1] for x in source_coords] + [x[1] for x in secondary_coords]
        #     zAll = [x[2] for x in source_coords] + [x[2] for x in secondary_coords]
        #     xAllNew = [x[0] for x in vert_coords] + [x[0] for x in second_vert]
        #     yAllNew = [x[1] for x in vert_coords] + [x[1] for x in second_vert]
        #     zAllNew = [x[2] for x in vert_coords] + [x[2] for x in second_vert]
        #     ax.scatter3D(xAll,yAll,zAll, color = 'b', alpha = 0.5)
        #     ax.scatter3D(xAllNew,yAllNew,zAllNew, color = 'g')
        # else:
        #     pass
        ax.scatter3D(xdata,ydata,zdata, color = 'b', alpha = 0.5)
        ax.scatter3D(xnew,ynew,znew, color = 'g')
        # Vectors
        X,Y,Z,U,V,W = zip(*refvects)
        ax.quiver(X,Y,Z,U,V,W, color = 'k', alpha = 0.75) # Axis Vectors
        ax.quiver(avgVect[0],avgVect[1],avgVect[2],avgVect[3],avgVect[4],avgVect[5], color = 'b', alpha = 0.5) # Source Vector
        ax.quiver(avgVectNew[0],avgVectNew[1],avgVectNew[2],avgVectNew[3],avgVectNew[4],avgVectNew[5], color = 'g') # Transformed Vector
        plt.show

    if type(secondary_coords) == np.ndarray:

        # Flipping along y axis (to conform to pre-existing code)
        for i in range(len(vert_coords)):
            y = vert_coords[i][1]
            vert_coords[i][1] = vert_coords[i][2]
            vert_coords[i][2] = y
        for i in range(len(second_vert)):
            y = second_vert[i][1]
            second_vert[i][1] = second_vert[i][2]
            second_vert[i][2] = y
        return vert_coords, second_vert
    else:
        for i in range(len(vert_coords)):
            y = vert_coords[i][1]
            vert_coords[i][1] = vert_coords[i][2]
            vert_coords[i][2] = y
        return vert_coords
