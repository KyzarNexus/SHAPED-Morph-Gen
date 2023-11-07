#%% 
import os
import pickle
import math
from mpl_toolkits import mplot3d
import numpy as np
from scipy.ndimage.measurements import label #; np.random.seed(65465555)
import matplotlib as mpl
import matplotlib.pyplot as plt
from ast import literal_eval
from function_repo import *


#%% Current File Purpose: Reading Distribution Files, Error Propogation


# if __name__ == "__main__":
    #     # dataset_identifier = 'amaral'
    #     # type_identifier = 'apical'

    #     # file_path = os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_{type_identifier}.txt')
    #     # file_name = f'{dataset_identifier}_distributions_{type_identifier}'
    #     # pickle_path = os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_{type_identifier}.pickle')

    #     # with open(pickle_path,'rb') as f:
    #     #     distribution_dict = pickle.load(f)

    #     # y = distribution_dict['max_y'][1:][0]
    #     # y_mean = 


    #     # sholl_bins = distribution_dict['sholl_bins'][1:][0] 
    #     # sholl_counts = distribution_dict['sholl_counts'][1:][0] 
    #     # sholl_avg = np.mean(sholl_counts,0)
        

    #     # # Error
    #     # sholl_error_abs = [np.divide(np.abs(np.subtract(sholl_avg,sholl_counts[i])),sholl_avg) for i in range(len(sholl_counts))]
    #     # sholl_error_abs = np.nan_to_num(sholl_error_abs)
    #     # sholl_error = [np.subtract(sholl_avg,sholl_counts[i]) for i in range(len(sholl_counts))]
    #     # # sholl_error = np.nan_to_num(sholl_error)
    #     # # sholl_error = (np.subtract(sholl_counts, sholl_avg))**2


    #     # # Error Summary
    #     # error_sums_abs = []
    #     # for i in range(len(sholl_error_abs)):
    #     #     error_sums_abs.append(sum(sholl_error_abs[i]))

    #     # error_sums = []
    #     # for i in range(len(sholl_error)):
    #     #     error_sums.append(sum(sholl_error[i]))

    #     # # Plotting Error Summary
    #     # # fig, ax = plt.subplots()
    #     # # ax.hist(error_sums_abs, bins = 15)
    #     # # ax.set_title('Error Distribution (Abs): ' + file_name )
    #     # # ax.set_xlabel('Magnitude')
    #     # # ax.set_ylabel('Frequency')

        
    
    #     # fig, ax = plt.subplots()
    #     # ax.hist(error_sums, bins = 15)
    #     # ax.set_title('Error Distribution: ' + file_name )
    #     # ax.set_xlabel('Magnitude')
    #     # ax.set_ylabel('Frequency')
    #     # mean = np.mean(error_sums)
    #     # median = np.median(error_sums)
    #     # sigma = np.std(error_sums)
    #     # textstr = f'\u03bc: {mean:.2f}\nMedian: {median:.2f}\n\u03C3: {sigma:.2f}'
    #     # ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    #     #     verticalalignment='top')
    #     # #plt.savefig(file_name + '_Error' + '.png')

    #     # # Plotting sholl counts
    #     # fig, ax = plt.subplots()
    #     # ax.set_xlabel('Sholl Bins')
    #     # ax.set_ylabel('Sholl Counts')
    #     # ax.set_title(file_name + ' All')

    #     # for i in range(len(sholl_counts)):
    #     #     ax.plot(sholl_bins,sholl_counts[i], color = 'b', alpha = 0.3, linewidth = 1)
    #     # ax.plot(sholl_bins, sholl_avg, color = 'k', linewidth = 2, label = 'Average')
    #     # ax.legend()
    #     # plt.savefig(file_name + '_All' + '.png')

    #     # fig, ax = plt.subplots()
    #     # ax.set_xlabel('Sholl Bins')
    #     # ax.set_ylabel('Sholl Counts')
    #     # ax.set_title(file_name + ' Average')
    #     # ax.plot(sholl_bins, sholl_avg)

    #     # plt.savefig(file_name + '_Avg' + '.png')

    #     # Plotting Sholl Error Over each sholl count

    #     # for i in range(len(sholl_counts)):
    #     #     fig, ax = plt.subplots()
    #     #     ax.set_xlabel('Sholl Bins')
    #     #     ax.set_ylabel('Sholl Counts')
    #     #     ax.set_title(file_name + ' ' + str(i))
    #     #     ax.plot(sholl_bins,sholl_counts[i], label = str(i))
    #     #     ax.plot(sholl_bins,sholl_avg, label = 'Average')
    #     #     ax.legend()

    #     #     ax2 = ax.twinx()
    #     #     color = 'tab:red'
    #     #     ax2.bar(sholl_bins,sholl_error_abs[i], width = 25, alpha = 0.3),
    #     #     ax2.set_ylim([0,2])
    #     #     ax2.set_ylabel('Error')

    #     #     plt.savefig(file_name + '_' + str(i) + '.png')
def compare_datasets(dataset_identifier, type_marker = '', synth_compare_flag = False):
    # Designed to co-plot two different distribution types. If the compare flag is set to false, only a single dataset will be shown.

    with open(os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_{type_marker}.pickle'),'rb') as f:
        dist_dict = pickle.load(f)
    if synth_compare_flag == True:
        with open(os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_{type_marker}_SYNTHETIC.pickle'),'rb') as f:
            synth_dict = pickle.load(f)

    # Plot: Normalized Frequency of Branch Pathlength(um)
    pathlength1 = dist_dict['Interbifurcation_Pathlength_(um)'][1:][0]
    pl_hist1 = np.histogram(pathlength1, 20)
    fig, ax = plt.subplots()
    ax.plot(pl_hist1[1],pl_hist1[0])

    # Plot: Normalized Frequency of Bifurcation Amplitude (deg)
        # Bifurcation_Amplitude_Vector(deg) or Bifurcation_Amplitude_Remote(deg) [Ask Gene for Clarification]
    bif_amp1 = dist_dict['Bifurcation_Amplitude_Remote(deg)'][1:][0]
    # Plot: Normalized Frequency of Bifurcation Azimuth (deg)
        # Doesn't explicitly exist in the distributions file. Possibly one of the prior mentioned dists
    bif_az1 = dist_dict['Bifurcation_Amplitude_Vector(deg)'][1:][0] # Temporary Dict ref
    # Plot: Normalized Frequency of Number of Bifurcations
    bif_num1 = dist_dict['Number_of_Bifurcations'][1:][0]
    # Plot: Normalized Frequency of Branch Order
    bif_order1 = dist_dict['Branch_Order'][1:][0]
    # Plot: Normalized Frequency of Total Dendritic Length
    tdl1 = dist_dict['tdl_distribution'][1:][0]

    # Plot: Avg Total Dendritic Length (um) vs. Normalized Distance from Soma
    # Plot: Number of Intersections vs. Normalized Distance from Soma
    
    
    
    
    pass

if __name__ == "__main__":
    dataset_identifier = 'amaral'
    type_marker = 'apical'
    synth_compare_flag = False
    compare_datasets(dataset_identifier,type_marker,synth_compare_flag)
    

# # # TDL:
    # tdl_distribution = distribution_dict['tdl_distribution'][1:][0]
    # max_euclidean = distribution_dict['max_euclidean'][1:][0] 
    # supermax_euclidean = max(max_euclidean); supermax_euclidean = int(math.ceil(supermax_euclidean))
    # tdl_bins = np.linspace(0, supermax_euclidean-1, supermax_euclidean)

    # tdl_distribution = [sum(x)/ len(tdl_distribution) for x in zip(*tdl_distribution)]
    # # [sum(x)/len(morphologies_list) for x in zip(*tmp_list)]
    # # Dimensions check
    # print(f'{len(tdl_distribution)} | {len(tdl_bins)}')

    # # tmp_list = [morphology.tdl_distribution for morphology in morphologies_list]
    # # tdl = [sum(x)/len(morphologies_list) for x in zip(*tmp_list)]
    # # plt.plot(tdl_bins, tdl, label='Generated TDL')  
    # # Plotting
    # fig, ax = plt.subplots()
    # ax.plot(tdl_bins, tdl_distribution, label = type_marker)
    # ax.set_xlabel('Euclidean Distance From Soma (\u03bcm)')
    # ax.set_ylabel('Total Dendritic Length (\u03bcm)')
    # # ax.set_title(f'Total Dendritic Length | {type_marker}')
    # # Temporary to test overlay
    # type_marker = 'basal'
    # with open(os.path.abspath(f'../../Data/Metrics/{dataset_identifier}_distributions_{type_marker}.pickle'),'rb') as f:
    #     distribution_dict = pickle.load(f)
    # tdl_distribution = distribution_dict['tdl_distribution'][1:][0]
    # max_euclidean = distribution_dict['max_euclidean'][1:][0] 
    # supermax_euclidean = max(max_euclidean); supermax_euclidean = int(math.ceil(supermax_euclidean))
    # tdl_bins = np.linspace(0, supermax_euclidean-1, supermax_euclidean)
    # tdl_distribution = [sum(x)/ len(tdl_distribution) for x in zip(*tdl_distribution)]
    # ax.plot(tdl_bins, tdl_distribution, label = type_marker)
        # ax.legend()









# %%
