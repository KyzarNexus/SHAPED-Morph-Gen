from __future__ import division

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
import itertools

def read_file(file_name):

    distribution_names = []
    distribution_types = []
    distribution_entries = []

    file_handle = open(file_name, 'r')
    lines = file_handle.read().splitlines()
    
    for i in range(0, len(lines)):
        parsed_line = lines[i].split(' ')

        tmp_name = parsed_line[0].split('_')

        distribution_names.append(' '.join(tmp_name))
        triple_parsed = []
        sholl_append = []
        
        del parsed_line[0]
        
        distribution_types.append(parsed_line[0])
        del parsed_line[0]
        if parsed_line[0] != '':
            if parsed_line[0][0] == '[':
    #            del parsed_line[0]
    #            del parsed_line[-1]
                for entry in parsed_line:
                    if entry[0] == '[':
                        double_parsed = []
                        double_parsed.append(entry[1:-1])
                    elif entry[-1] == ']':
                        double_parsed.append(entry[:-1])
                        triple_parsed.append(double_parsed)
                    else:
                        double_parsed.append(entry[:-1])
                for j in range(0, len(triple_parsed)):
                    sholl_append.append([float(x) for x in triple_parsed[j]])
                distribution_entries.append(sholl_append)
            else:
                distribution_entries.append(map(float,parsed_line))
        else:
            distribution_entries.append([0])
    
    return distribution_names, distribution_types, distribution_entries

def plot_2_distributions(distribution_1, distribution_2, distribution_name):
    fig = plt.figure()
    bin_num = 25
    both_distributions = distribution_1 + distribution_2
    [hist_peaks, bin_edges] =np.histogram(both_distributions, bins=bin_num) #get the bin edges
    [hist_peaks1, bin_edges] =np.histogram(distribution_1, bins=bin_num)
    [hist_peaks2, bin_edges] =np.histogram(distribution_2, bins=bin_num)
    
#    print bin_edges
    bin_edges = np.ndarray.tolist(bin_edges)
    
    del bin_edges[-1]
    plt.hist(distribution_1, alpha=0.7, normed=True, bins=bin_edges, label=legend_name_1)
    plt.hist(distribution_2, alpha=0.7, normed=True, bins=bin_edges, label=legend_name_2)
#    plt.plot(bin_edges, hist_peaks1, color='black', linestyle='solid', label=legend_name_1)
#    plt.plot(bin_edges, hist_peaks2, color='black', linestyle='dashed',label=legend_name_2)
    plt.xlabel(distribution_name, fontsize=18)
    plt.ylabel('Frequency', fontsize=18)

    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(distribution_1))
#    print lnspc
#    be, ce = stats.expon.fit(distribution_1)
#    pdf_expon = stats.expon.pdf(lnspc, be, ce)
#    plt.plot(lnspc, pdf_expon, label = "Exp")
#    m, s = stats.norm.fit(distribution_1)
#    ag, bg, cg = stats.gamma.fit(distribution_1)
#    pdf_gamma = stats.gamma.pdf(lnspc, ag, bg, cg)
#    pdf_norm = stats.norm.pdf(lnspc, m, s)
#    plt.plot(lnspc, pdf_gamma, label="Fitted Distribution", color='red')
#    plt.plot(lnspc, pdf_norm, label="Norm")
    
#    print("Exponential Coef", be, ce)
#    print("Normal Coef", m, s)
#    print("Gamma Coef", ag, bg, cg)


    distribution_1.sort()
    distribution_2.sort()
    percentile1 = [0]*len(distribution_1)
    percentile2 = [0]*len(distribution_2)
    
    delta = 1
    max1 = math.ceil(distribution_1[-1])
    max2 = math.ceil(distribution_2[-1])
    index1 = np.linspace(0, max1, int(max1/delta)).tolist()
    index2 = np.linspace(0, max2, int(max2/delta)).tolist()
    cumulative1 = [0]*int(max1/delta)
    cumulative2 = [0]*int(max2/delta)
    
    for i in range(0, int(max1/delta)):
        matches = [ii for ii, x in enumerate(distribution_1) if x > (i+1)*delta]
#        print matches
        if matches == []:
            cumulative1[i] = 1
        else:
            cumulative1[i] = (matches[0]+1)/(len(distribution_1)+1)

    for i in range(0, int(max2/delta)):
        matches = [ii for ii, x in enumerate(distribution_2) if x > (i+1)*delta]
        if matches == []:
            cumulative2[i] = 1
        else:
            cumulative2[i] = (matches[0]+1)/(len(distribution_2)+1)

    for i in range(0, len(distribution_1)):
        percentile1[i] = (i+1)/(len(distribution_1)+1)
        
    for j in range(0, len(distribution_2)):
        percentile2[j] = (j+1)/(len(distribution_2)+1)


#    fig = plt.figure()
#    plt.plot(index1, cumulative1, color='red')
#    plt.plot(index2, cumulative2, color='blue')
#    plt.show()
#    
#    fig = plt.figure()
#    plt.plot(distribution_1, percentile1, color='red')
#    plt.plot(distribution_2, percentile2, color='blue')
##    plt.yscale('log')
##    plt.xscale('log')
#    plt.show()
    
    KS_stat, p_val = stats.ks_2samp(distribution_1, distribution_2)
        
    print '\nKS stat:', KS_stat, '\n', 'p value:', p_val


    #plot KS-plot
    
#
    plt.legend()
    plt.show()

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
    
def plot_3_distributions(distribution_1, distribution_2, distribution_3, distribution_name):
    fig = plt.figure()
    ax = plt.axes()
    bin_num = 8
    both_distributions = distribution_1 + distribution_2 + distribution_3
    
#    print both_distributions, bin_num
    [hist_peaks, bin_edges] =np.histogram(both_distributions, normed=True, bins=bin_num) #get the bin edges
    [hist_peaks1, bin_edges1] =np.histogram(distribution_1, normed=True, bins=bin_num)
    [hist_peaks2, bin_edges2] =np.histogram(distribution_2, normed=True, bins=bin_num)
    [hist_peaks3, bin_edges3] =np.histogram(distribution_3, normed=True, bins=bin_num)
    
#    print bin_edges
    bin_edges = np.ndarray.tolist(bin_edges)
    bin_edges1 = np.ndarray.tolist(bin_edges1)
    bin_edges2 = np.ndarray.tolist(bin_edges2)
    bin_edges3 = np.ndarray.tolist(bin_edges3)
    
    
    
    
    del bin_edges[-1]
    del bin_edges1[-1]
    del bin_edges2[-1]
    del bin_edges3[-1]
   
#    plt.hist(distribution_1, alpha=0.6, normed=True, bins=bin_edges, label=legend_name_1)
#    plt.hist(distribution_2, alpha=0.6, normed=True, bins=bin_edges, label=legend_name_2)
#    plt.hist(distribution_3, alpha=0.6, normed=True, bins=bin_edges, label=legend_name_3)
    
    smoothing_factor = 3
    hist_peaks1 = smooth(hist_peaks1, smoothing_factor)
    hist_peaks2 = smooth(hist_peaks2, smoothing_factor)
    hist_peaks3 = smooth(hist_peaks3, smoothing_factor)
    
    for i in range(0, len(hist_peaks1)):
        if hist_peaks1[i] == 0:
            hist_peaks1[i] = 0.00000001
    for i in range(0, len(hist_peaks2)):
        if hist_peaks2[i] == 0:
            hist_peaks2[i] = 0.00000001
    for i in range(0, len(hist_peaks3)):
        if hist_peaks3[i] == 0:
            hist_peaks3[i] = 0.00000001

    p1 = [hp * (bin_edges1[1] - bin_edges1[0]) for hp in hist_peaks1]
    p2 = [hp * (bin_edges2[1] - bin_edges2[0]) for hp in hist_peaks2]
    p3 = [hp * (bin_edges3[1] - bin_edges3[0]) for hp in hist_peaks3]


    p1 = [max(hist_peaks1)]
    p2 = [max(hist_peaks2)]
    p3 = [max(hist_peaks3)]

    sf1 = 1
    sf2 = 1
    sf3 = 1

    sf1 = 1/sum(p1)
    sf2 = 1/sum(p2)
    sf3 = 1/sum(p3)
    
#    hist_peaks1 = [hp * sf1 for hp in hist_peaks1]
#    hist_peaks2 = [hp * sf2 for hp in hist_peaks2]
#    hist_peaks3 = [hp * sf3 for hp in hist_peaks3]

    palette = itertools.cycle(sns.color_palette())

    plt.plot(bin_edges1, hist_peaks1, color=next(palette), linestyle='solid', label=legend_name_1)
    color=next(palette)
    plt.plot(bin_edges2, hist_peaks2, color=next(palette), linestyle='solid',label=legend_name_2)
    color=next(palette)
    color=next(palette)
    color=next(palette)

    plt.plot(bin_edges3, hist_peaks3, color=next(palette), linestyle='solid',label=legend_name_3)    



    smallest_distribution = min(len(distribution_1), len(distribution_2), len(distribution_3))

    distribution_1_trimmed = distribution_1[0:smallest_distribution]
    distribution_2_trimmed = distribution_2[0:smallest_distribution]
    distribution_3_trimmed = distribution_3[0:smallest_distribution]

#    print len(hist_peaks1), len(hist_peaks2)

    kl1 = stats.entropy(hist_peaks1, hist_peaks2)
    kl2 = stats.entropy(hist_peaks1, hist_peaks3)
    
#    print kl1, kl2

#    plt.title('zmorph KL:{} lneuron KL:{}'.format(kl1, kl2), fontsize=18)
    plt.xlabel(distribution_name)
    plt.ylabel('Normalized Frequency')

#    plt.legend()
    sns.despine()
    plt.show()    
    
def plot_4_distributions(distribution_1, distribution_2, distribution_3, distribution_4, distribution_name):
    fig = plt.figure()
    ax = plt.axes()
    bin_num = 8
    both_distributions = distribution_1 + distribution_2 + distribution_3 + distribution_4
    
#    print both_distributions, bin_num
    [hist_peaks, bin_edges] =np.histogram(both_distributions, normed=True, bins=bin_num) #get the bin edges
    [hist_peaks1, bin_edges1] =np.histogram(distribution_1, normed=True, bins=bin_edges)
    [hist_peaks2, bin_edges2] =np.histogram(distribution_2, normed=True, bins=bin_edges)
    [hist_peaks3, bin_edges3] =np.histogram(distribution_3, normed=True, bins=bin_edges)
    [hist_peaks4, bin_edges4] = np.histogram(distribution_4, normed=True, bins=bin_edges)
    
#    print bin_edges
    bin_edges = np.ndarray.tolist(bin_edges)
    bin_edges1 = np.ndarray.tolist(bin_edges1)
    bin_edges2 = np.ndarray.tolist(bin_edges2)
    bin_edges3 = np.ndarray.tolist(bin_edges3)
    bin_edges4 = np.ndarray.tolist(bin_edges4)
    
    bin_edges1 = list(bin_edges)
    bin_edges2 = list(bin_edges)
    bin_edges3 = list(bin_edges)
    bin_edges4 = list(bin_edges)
    
    
    
    del bin_edges[-1]
    del bin_edges1[-1]
    del bin_edges2[-1]
    del bin_edges3[-1]
    del bin_edges4[-1]
   
#    plt.hist(distribution_1, alpha=0.6, normed=True, bins=bin_edges, label=legend_name_1)
#    plt.hist(distribution_2, alpha=0.6, normed=True, bins=bin_edges, label=legend_name_2)
#    plt.hist(distribution_3, alpha=0.6, normed=True, bins=bin_edges, label=legend_name_3)
    
    smoothing_factor = 3
    hist_peaks1 = smooth(hist_peaks1, smoothing_factor)
    hist_peaks2 = smooth(hist_peaks2, smoothing_factor)
    hist_peaks3 = smooth(hist_peaks3, smoothing_factor)
    hist_peaks4 = smooth(hist_peaks4, smoothing_factor)
    
#    both_distributions = distribution_1 + distribution_2 + distribution_3 + distribution_4
#    bin_edges=np.histogram(both_distributions, bins=50)[1] #get the bin edges

    
    for i in range(0, len(hist_peaks1)):
        if hist_peaks1[i] == 0:
            hist_peaks1[i] = 0.00000001
    for i in range(0, len(hist_peaks2)):
        if hist_peaks2[i] == 0:
            hist_peaks2[i] = 0.00000001
    for i in range(0, len(hist_peaks3)):
        if hist_peaks3[i] == 0:
            hist_peaks3[i] = 0.00000001
    for i in range(0, len(hist_peaks4)):
        if hist_peaks4[i] == 0:
            hist_peaks4[i] = 0.000000001

    p1 = [hp * (bin_edges1[1] - bin_edges1[0]) for hp in hist_peaks1]
    p2 = [hp * (bin_edges2[1] - bin_edges2[0]) for hp in hist_peaks2]
    p3 = [hp * (bin_edges3[1] - bin_edges3[0]) for hp in hist_peaks3]
    p4 = [hp * (bin_edges4[1] - bin_edges4[0]) for hp in hist_peaks4]


    p1 = [max(hist_peaks1)]
    p2 = [max(hist_peaks2)]
    p3 = [max(hist_peaks3)]
    p4 = [max(hist_peaks4)]

    sf1 = 1
    sf2 = 1
    sf3 = 1
    sf4 = 1

    sf1 = 1/sum(p1)
    sf2 = 1/sum(p2)
    sf3 = 1/sum(p3)
    sf4 = 1/sum(p4)
    
#    hist_peaks1 = [hp * sf1 for hp in hist_peaks1]
#    hist_peaks2 = [hp * sf2 for hp in hist_peaks2]
#    hist_peaks3 = [hp * sf3 for hp in hist_peaks3]

    err1 = 0
    err2 = 0
    err3 = 0
    cum_total = 0
    for zz in range(len(hist_peaks1)):
        err1 += abs(hist_peaks1[zz] - hist_peaks2[zz])*(bin_edges[1] - bin_edges[0])
        err2 += abs(hist_peaks1[zz] - hist_peaks3[zz])*(bin_edges[1] - bin_edges[0])
        err3 += abs(hist_peaks1[zz] - hist_peaks4[zz])*(bin_edges[1] - bin_edges[0])
        cum_total += hist_peaks1[zz]*(bin_edges[1] - bin_edges[0])
    
    print cum_total
    print err1/cum_total, err2/cum_total, err3/cum_total

    palette = itertools.cycle(sns.color_palette())

    plt.plot(bin_edges1, hist_peaks1, color=next(palette), linestyle='solid', label=legend_name_1)
    color=next(palette)
    plt.plot(bin_edges2, hist_peaks2, color=next(palette), linestyle='solid',label=legend_name_2)
    plt.plot(bin_edges3, hist_peaks3, color=next(palette), linestyle='solid',label=legend_name_3)

#    color=next(palette)
    color=next(palette)
    color=next(palette)

    plt.plot(bin_edges4, hist_peaks4, color=next(palette), linestyle='solid',label=legend_name_4)

#    plt.legend()
    plt.xlabel(distribution_name)
    plt.ylabel('Normalized Frequency')
    sns.despine()
    plt.savefig('../Figures/{}_compare4.svg'.format(distribution_name), format='svg', dpi=1200)

    plt.show()        
    
def plot_conditional_relationships(distribution_names, distributions, basic_parameters, conditional_parameters):
    for i in range(0, len(basic_parameters)):
        for j in range(0, len(conditional_parameters)):
            basic_index = basic_parameters[i]
            conditional_index = conditional_parameters[j]
            if len(distributions[conditional_index]) == len(distributions[basic_index]):
                
                print (distribution_names[basic_index] + ' depending on ' + distribution_names[conditional_index])
                plt.scatter(distributions[conditional_index], distributions[basic_index])
                plt.show()

def compare_conditional_relationships(distribution_names, distributions_1, distributions_2, basic_parameters, conditional_parameters):
    for i in range(0, len(basic_parameters)):
        for j in range(0, len(conditional_parameters)):
            basic_index = basic_parameters[i]
            conditional_index = conditional_parameters[j]
            if len(distributions_1[conditional_index]) == len(distributions_1[basic_index]):
                
                print (distribution_names[basic_index] + ' depending on ' + distribution_names[conditional_index])
                plt.scatter(distributions_1[conditional_index], distributions_1[basic_index], alpha=0.6)
                plt.scatter(distributions_2[conditional_index], distributions_2[basic_index], alpha=0.6)
                plt.show()

file_name_1 = '../Data/Metrics/claiborne_distributions.txt'
file_name_2 = '../Data/Metrics/zmorph_tester_distributions.txt'
file_name_3 = '../Data/Metrics/LNeuron_recreated_distributions.txt'
file_name_4 = '../Data/Metrics/flat_distributions.txt'

plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
#plt.rc('xlabel', labelsize=18)
#plt.rc('ylabel', labelsize=18)
legend_name_1 = file_name_1.split('/')[-1].split('_')[0]
legend_name_2 = file_name_2.split('/')[-1].split('_')[0]
legend_name_3 = file_name_3.split('/')[-1].split('_')[0]
legend_name_4 = file_name_4.split('/')[-1].split('_')[0]

legend_name_1 = 'Real'
legend_name_2 = 'L-ARBOR'
legend_name_3 = 'L-Neuron'
legend_name_4 = 'Flat Branch Rate'
sns.set_style("ticks")
#sns.set_context("poster", font_scale = 2.5, rc={"lines.linewidth":10.5})
sns.set_context("talk")
sns.set_palette("deep")

distribution_names_1, distribution_type, distributions_1 = read_file(file_name_1)
distribution_names_2, distribution_type, distributions_2 = read_file(file_name_2)
distribution_names_3, distribution_type, distributions_3 = read_file(file_name_3)
distribution_names_4, distribution_type, distributions_4 = read_file(file_name_4)
basic_parameters = []
conditional_parameters = []

for x in range(0, len(distribution_type)):
    if distribution_type[x] == 'basic':
        basic_parameters.append(x)
    elif distribution_type[x] == 'conditional':
        conditional_parameters.append(x)
# plot_conditional_relationships(distribution_names_1, distributions_1, basic_parameters, conditional_parameters)
#compare_conditional_relationships(distribution_names_1, distributions_1, distributions_2, basic_parameters, conditional_parameters)

for i in range(0, min(len(distribution_names_1), len(distribution_names_2))):
    if distribution_names_1[i] == 'sholl bins':
        sholl_bins = distributions_1[i]
    elif distribution_names_1[i] == 'sholl counts':
        sholl_avgs_1 = [0]*len(sholl_bins)
        sholl_avgs_2 = [0]*len(sholl_bins)
        sholl_avgs_3 = [0]*len(sholl_bins)
        sholl_avgs_4 = [0]*len(sholl_bins)
        normed_counts = [0]*len(sholl_bins)
        for j in range(0, len(sholl_bins)):
            for entry in distributions_1[i]:
                sholl_avgs_1[j]  = sholl_avgs_1[j] + entry[j]/len(distributions_1[i])
            for entry2 in distributions_2[i]:
                sholl_avgs_2[j] = sholl_avgs_2[j] + entry2[j]/len(distributions_2[i])
            for entry3 in distributions_3[i]:
                sholl_avgs_3[j] = sholl_avgs_3[j] + entry3[j]/len(distributions_3[i])
            for entry4 in distributions_4[i]:
                sholl_avgs_4[j] = sholl_avgs_4[j] + entry4[j]/len(distributions_4[i])

#        sholl_avgs_1 = [x/max(sholl_avgs_1) for x in sholl_avgs_1]
#        sholl_avgs_2 = [x/max(sholl_avgs_2) for x in sholl_avgs_2]
#        sholl_avgs_3 = [x/max(sholl_avgs_3) for x in sholl_avgs_3]
#        sholl_avgs_4 = [x/max(sholl_avgs_4) for x in sholl_avgs_4]
        sholl_avgs_1 = [x for x in sholl_avgs_1]
        sholl_avgs_2 = [x for x in sholl_avgs_2]
        sholl_avgs_3 = [x for x in sholl_avgs_3]
        sholl_avgs_4 = [x for x in sholl_avgs_4]
        sholl_bins = [b/sholl_bins[-1] for b in sholl_bins]
        
        
        print '\n\n\nSholl Analysis'
        palette = itertools.cycle(sns.color_palette())

        plt.plot(sholl_bins, sholl_avgs_1, marker='', color=next(palette), label=legend_name_1)
#        print sholl_bins, sholl_avgs_1
        color = next(palette)
        plt.plot(sholl_bins, sholl_avgs_2, marker='', color = next(palette), label=legend_name_2)
        plt.plot(sholl_bins, sholl_avgs_3, marker='', color = next(palette), label=legend_name_3)
#        color = next(palette)
        color = next(palette)
        color = next(palette)
        plt.plot(sholl_bins, sholl_avgs_4, marker='', color = next(palette), label=legend_name_4)
#        plt.plot(sholl_bins, sholl_avgs_4, marker='o', label=legend_name_4)
        plt.xlabel(u"Normalized Distance from Soma")
        plt.ylabel('Number of intersections')
        sns.despine()
#        plt.legend()
        plt.savefig('../Figures/sholl_compare4.svg', format='svg', dpi=1200)

        plt.show()
    elif distribution_names_1[i] == 'tdl distribution':
        
        min_len = 10000000000
        for distribution in distributions_1[i]:
            if len(distribution) < min_len:
                min_len = len(distribution)
        for distribution in distributions_1[i]:
            if len(distribution) > min_len:
                del distribution[-1]
        min_len = 10000000000
        for distribution in distributions_2[i]:
            if len(distribution) < min_len:
                min_len = len(distribution)
        for distribution in distributions_2[i]:
            if len(distribution) > min_len:
                del distribution[-1]
        min_len = 10000000000
                
        for distribution in distributions_3[i]:
            if len(distribution) < min_len:
                min_len = len(distribution)
        for distribution in distributions_3[i]:
            if len(distribution) > min_len:
                del distribution[-1]
                
        min_len = 10000000000
                
        for distribution in distributions_4[i]:
            if len(distribution) < min_len:
                min_len = len(distribution)
        for distribution in distributions_4[i]:
            if len(distribution) > min_len:
                del distribution[-1]                
        
        tdl_avgs_1 = distributions_1[i][0]
        tdl_avgs_2 = distributions_2[i][0]
        tdl_avgs_3 = distributions_3[i][0]
        tdl_avgs_4 = distributions_4[i][0]

        tdl_bins1 = np.ndarray.tolist(np.linspace(0, len(tdl_avgs_1)-1, len(tdl_avgs_1)))
        tdl_bins2 = np.ndarray.tolist(np.linspace(0, len(tdl_avgs_2)-1, len(tdl_avgs_2)))
        tdl_bins3 = np.ndarray.tolist(np.linspace(0, len(tdl_avgs_3)-1, len(tdl_avgs_3)))       
        tdl_bins4 = np.ndarray.tolist(np.linspace(0, len(tdl_avgs_4)-1, len(tdl_avgs_4)))
        
        
#        print len(distributions_3[i])
#
        gcl_c1= []
        gcl_c2 = []
        first_third_c1 = []
        second_third_c1 = []
        first_third_c2 = []
        second_third_c2 = []
        first_third_c3 = []
        second_third_c3 = []
        for tdl_dist in distributions_1[i]:
            for k in range(1, len(tdl_dist)):
                if tdl_dist[len(tdl_dist)-k] != tdl_dist[-1]:
                    break
            max_idx = len(tdl_dist) - k
            gcl = max_idx - 350
            if gcl < 0:
                gcl = 0
            first_third = int(117 + gcl)
            second_third = int(233 + gcl)
            total_third = max_idx
        
            first_third_c1.append(tdl_dist[first_third])
            second_third_c1.append(tdl_dist[second_third])
        for tdl_dist2 in distributions_2[i]:
            for k in range(1, len(tdl_dist2)):
                if tdl_dist2[len(tdl_dist2)-k] != tdl_dist2[-1]:
                    break
            max_idx2 = len(tdl_dist2) - k
            gcl2 = max_idx2 - 350
            if gcl2 < 0:
                gcl2 = 0
            first_third2 = int(117+ gcl)
            second_third2 = int(233 + gcl)
            total_third2 = max_idx2
        
            first_third_c2.append(tdl_dist2[first_third2])
            second_third_c2.append(tdl_dist2[second_third2])
        for tdl_dist3 in distributions_3[i]:
            for k in range(1, len(tdl_dist3)):
                if tdl_dist3[len(tdl_dist3)-k] != tdl_dist3[-1]:
                    break
#            print k
            max_idx3 = len(tdl_dist3) - k
            gcl3 = max_idx3 - 350
            if gcl3 < 0:
                 gcl3 = 0
            first_third3 = int(117+ gcl)
            second_third3 = int(233 + gcl)
            total_third3 = max_idx3
        
            first_third_c3.append(tdl_dist3[first_third3])
            second_third_c3.append(tdl_dist3[second_third3])
            
        print 'cdf1', np.mean(first_third_c1)/tdl_dist[max_idx], np.mean(second_third_c1)/tdl_dist[max_idx] - np.mean(first_third_c1)/tdl_dist[max_idx]
        print 'cdf2', np.mean(first_third_c2)/tdl_dist2[max_idx2], np.mean(second_third_c2)/tdl_dist2[max_idx2] - np.mean(first_third_c2)/tdl_dist2[max_idx2]
        print 'cdf3', np.mean(first_third_c3)/tdl_dist3[max_idx3], np.mean(second_third_c3)/tdl_dist3[max_idx3] - np.mean(first_third_c3)/tdl_dist3[max_idx3]
        
        cdf1 = [sum(x)/len(distributions_1[i]) for x in zip(*distributions_1[i])]
        pdf1 = np.diff(cdf1)
#        del tdl_bins1[-1]
        cdf2 = [sum(x)/len(distributions_2[i]) for x in zip(*distributions_2[i])]
        pdf2 = np.diff(cdf2)
#        del tdl_bins2[-1]
        cdf3 = [sum(x)/len(distributions_3[i]) for x in zip(*distributions_3[i])]
        pdf3 = np.diff(cdf3)
#        del tdl_bins3[-1]
        cdf4 = [sum(x)/len(distributions_4[i]) for x in zip(*distributions_4[i])]

        
#        tdl_bins1 = [b/tdl_bins1[-1] for b in tdl_bins1]
#        tdl_bins2 = [b/tdl_bins2[-1] for b in tdl_bins2]
#        tdl_bins3 = [b/tdl_bins3[-1] for b in tdl_bins3]
        
        smoothing_factor = 200
        pdf1 = smooth(pdf1, smoothing_factor)
        pdf2 = smooth(pdf2, smoothing_factor)
        pdf3 = smooth(pdf3, smoothing_factor)        
        
        tdl_bins2.append(600)
        cdf2.append(cdf2[-1]*1.0)
#        tdl_bins3.append(600)
#        cdf3.append(cdf3[-1]*1.0)
        tdl_bins4.append(600)
        cdf4.append(cdf4[-1]*1.0)        

        plt.figure()
        palette = itertools.cycle(sns.color_palette())

        plt.plot(tdl_bins1, cdf1, color=next(palette))
        color=next(palette)
        plt.plot(tdl_bins2, cdf2, color= next(palette))
        plt.plot(tdl_bins3, cdf3, color= next(palette))  
#        color=next(palette)
        color=next(palette)
        color=next(palette)
        plt.plot(tdl_bins4, cdf4, color=next(palette))   
#        plt.title('TDL Analysis')
        plt.xlabel(u"Y-Coordinate (\u03bcm)")
        plt.ylabel(u"Avg Total Dendritic Length (\u03bcm)")
        plt.xlim([0,600])
#        plt.ylim()
        sns.despine()
#        plt.legend()
        plt.savefig('../Figures/tdl_compare4.svg', format='svg', dpi=1200)
        plt.show()

#        tdl_avgs_1 = [0]*len(sholl_bins)
#        tdl_avgs_2 = [0]*len(sholl_bins)
#        tdl_avgs_3 = [0]*len(sholl_bins)
#        for j in range(0, len(distributions_1[i][0])):
#            for entry in distributions_1[i]:
#                tdl_avgs_1[j]  = sholl_avgs_1[j] + entry[j]/len(distributions_1[i])
#            for entry2 in distributions_2[i]:
#                tdl_avgs_2[j] = sholl_avgs_2[j] + entry2[j]/len(distributions_2[i])
#            for entry3 in distributions_3[i]:
#                tdl_avgs_3[j] = sholl_avgs_3[j] + entry3[j]/len(distributions_3[i])
#
#        print '\n\n\nTDL Analysis'
#        plt.plot(tdl_bins, tdl_avgs_1, marker='', color='C0', label=legend_name_1)
#        print tdl_bins, tdl_avgs_1
#        plt.plot(tdl_bins, tdl_avgs_2, marker='', color = 'C1', label=legend_name_2)
#        plt.plot(tdl_bins, tdl_avgs_3, marker='', color = 'C2', label=legend_name_3)               
                
   
    else:
        if i >=0 and i <= 200:
#            print '\n\n\n' + distribution_names_1[i]
#            plot_2_distributions(distributions_1[i], distributions_2[i], distribution_names_1[i])
#            plot_3_distributions(distributions_1[i], distributions_2[i], distributions_3[i], distribution_names_1[i])
            plot_4_distributions(distributions_1[i], distributions_2[i], distributions_3[i], distributions_4[i], distribution_names_1[i])

        
    