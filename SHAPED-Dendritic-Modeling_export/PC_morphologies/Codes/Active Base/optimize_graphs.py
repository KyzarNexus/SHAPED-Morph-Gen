#%%
# Import
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
import itertools as iter
import pickle
import pandas as  pd
import plotly 
import plotly.graph_objs as go

#%%
# Downloading Picked Files
# pickle_in = open("optPicFull.pickle","rb"); opDict = pickle.load(pickle_in)
pickle_in = open('ksDict.pickle','rb'); lambDict = pickle.load(pickle_in) #[rescaled_list,b_k,b_k_upper,b_k_lower]

lambDictPlotFlag = 1
if lambDictPlotFlag == 1:
    keys = lambDict.keys()
    count = 0
    for key in keys:
        count+=1
        rescaled_list, b_k, b_k_upper, b_k_lower = lambDict[key]
        plt.figure(count)
        label = 'ID: '+str(count)+' | '+str([round(x,3) for x in key])
        plt.title(label)
        plt.plot(b_k, rescaled_list)
        plt.plot(b_k, b_k, color='black')
        plt.plot(b_k, b_k_upper, linestyle='dashed', color='black')
        plt.plot(b_k, b_k_lower, linestyle='dashed', color='black')


#%% Plotting Branching Rate per initial condition
pickle_in = open('lambPic.pickle','rb'); lambDict = pickle.load(pickle_in)
ratePlotFlag = 1
if ratePlotFlag == 1:
    keys = lambDict.keys()
    count = 0
    for key in keys:
        count+=1
        lamb = lambDict[key]
        xspan = list(range(len(lamb)))
        plt.figure(count)
        label = 'ID: '+ str(count) +' | '+ str([round(x,3) for x in key])
        plt.title(label)
        plt.plot(xspan,lamb)





#%%
# Setting Up for Plotting

#%% Individual Trend (Outdated)
opDictPlotFlag = 0
if opDictPlotFlag == 1: 
    keys = opDict.keys()
    for key in keys:
        dictLength = len(opDict[key])
        iters = range(1,dictLength+1) # X Values
        rhoVals = []; varVals = []; betaVals = []; muVals = [] # Y Values
        for i in range(dictLength):
            rhoVals.append(opDict[key][i][0]); varVals.append(opDict[key][i][1]); betaVals.append(opDict[key][i][2]); muVals.append(opDict[key][i][3])
        
        # Plotting subplots for dictEntry
        fig, axs = plt.subplots(2,2)
        fig.suptitle(str(key))
        axs[0,0].plot(iters,rhoVals) ; axs[0,0].set_title('Rho Trend')
        axs[0,1].plot(iters,varVals) ; axs[0,1].set_title('Variance Trend')
        axs[1,0].plot(iters,betaVals) ; axs[1,0].set_title('Beta Trend')
        axs[1,1].plot(iters,muVals) ; axs[1,1].set_title('Mu Trend')
        
        for ax in axs.flat:
            ax.set(xlabel='Iterations', ylabel='Values')
        for ax in axs.flat:
            ax.label_outer()
# %% Compressing to Relative Error and 'Convergence' Data
convFlag = 0
if convFlag == 1:
    relE_Dict = dict({key : None for key in keys}); con_Dict = dict({key : None for key in keys})
    for key in keys:
        if opDict[key][-2] == 0:
            relE_Dict[key] = opDict[key][-1]
        else:
            relE_Dict[key] = abs(np.divide(np.subtract(opDict[key][-1],opDict[key][-2]),opDict[key][-2]))
        con_Dict[key] = np.nanmean(relE_Dict[key]) # # Legend: Rho, Variance, Beta, Mu, 'Convergence'
    #%%


    convThresh = 0.01; threshList = []; threshKeys = [] # Sorting for Highest Convergence Values
    rhoList = []; varList = []; betaList = []; muList = []
    for key in keys:
        if con_Dict[key] <= convThresh:
            threshKeys.append(key)
            threshList.append(opDict[key][-1] + [con_Dict[key]])
            rhoList.append(threshList[-1][0]); varList.append(threshList[-1][1]); betaList.append(threshList[-1][2]); muList.append(threshList[-1][3])

    #%% Saving filtered list for further use.
    #%%
    pickle_out = open("convergence_filtered_keys.pickle","wb")
    pickle.dump(threshKeys,pickle_out); pickle_out.close()
#%%
plot_flag2 = 0
if plot_flag2 == 1:
    data = pd.DataFrame({
        'Rho': rhoList,
        'Variance': varList,
        'Beta': betaList,
        'Mu': muList})
    # # Removing Outliers (Assuming Single Optimization End Value. Comment out if otherwise)
    # data = data[data['Rho'].between(data['Rho'].quantile(.10), data['Rho'].quantile(.90))]
    # data = data[data['Variance'].between(data['Variance'].quantile(.10), data['Variance'].quantile(.90))]
    # data = data[data['Beta'].between(data['Beta'].quantile(.10), data['Beta'].quantile(.90))]
    data = data[data['Mu'].between(data['Mu'].quantile(.10), data['Mu'].quantile(.90))]
    # # Popping Non-relevant Values
    data = data[data['Rho'] != 0]; data = data[data['Beta'] != 0]
    # Marker Properties
        #Set marker properties
    markercolor = data['Variance'].values 

    #Make Plotly figure
    fig1 = go.Scatter3d(x=data['Beta'].values,
                        y=data['Mu'].values,
                        z=data['Rho'].values,
                        marker=dict(size=5,
                                    color=markercolor,
                                    opacity= .9,
                                    reversescale=True,
                                    colorscale='plasma'),
                        line=dict (width=0.02),
                        mode='markers')

    #Make Plot.ly Layout
    mylayout = go.Layout(scene=dict(xaxis=dict( title="Beta"),
                                    yaxis=dict( title="Mu"),
                                    zaxis=dict(title="Rho")),)

    #Plot and save html
    plotly.offline.plot({"data": [fig1],
                        "layout": mylayout},
                        auto_open=True,
                        filename=("5D Plot.html"))
#%%
plot_flag1 = 0 # Plotting convergence data for each parameter
if plot_flag1 == 1:
    
    # # # Plotting Section # # #
    #%% 5D Array Generation
    key_list = np.array([list(key) for key in keys])
    full_list = [key_list[:,0],key_list[:,1],key_list[:,2],key_list[:,3],np.array([con_Dict[key] for key in keys])]
    # # full_list Legend: Rho, Variance, Beta, Mu, 'Convergence'
    # Graph Planning:
    # Beta - Axis(x) # Mu - Axis(y) # Rho - Axis(z) # Variance - Color Gradient (Label in Graph) # Conv - Opacity
    data = pd.DataFrame({
        'Beta': full_list[2], # x
        'Mu': full_list[3], # y
        'Rho': full_list[0], # z
        'Variance': full_list[1], 
        'Convergence': full_list[4]})

    #%% Plotting

    #Set marker properties
    markersize = data['Variance'].values * 10
    markercolor = data['Convergence'].values 

    #Make Plotly figure
    fig1 = go.Scatter3d(x=data['Beta'].values,
                        y=data['Mu'].values,
                        z=data['Rho'].values,
                        marker=dict(size=markersize,
                                    color=markercolor,
                                    opacity=0.9,
                                    reversescale=True,
                                    colorscale='plasma'),
                        line=dict (width=0.02),
                        mode='markers')

    #Make Plot.ly Layout
    mylayout = go.Layout(scene=dict(xaxis=dict( title="Beta"),
                                    yaxis=dict( title="Mu"),
                                    zaxis=dict(title="Rho")),)

    #Plot and save html
    plotly.offline.plot({"data": [fig1],
                        "layout": mylayout},
                        auto_open=True,
                        filename=("5D Plot.html"))





    # %%
