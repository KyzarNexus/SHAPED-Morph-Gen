from __future__ import division
import glob
import morphology_classes as mc

# Load real morphologies

file_list = glob.glob('../Data/Morphologies/claiborne/CNG version/*.swc')
distributions_file = '../Data/Metrics/claiborne_distributions.txt'

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
    if i > 1:
        continue
    
    print file_list[i]
    compartment_list, branch_list, bifurcation_list, termination_list = mc.read_file(file_list[i])

    if compartment_list == [] and branch_list == []:
        continue    

    all_morphologies_list.append(mc.Morphology(compartment_list, branch_list, bifurcation_list, termination_list))

    all_compartments_list = all_compartments_list  + compartment_list
    all_branches_list = all_branches_list + branch_list
    all_bifurcations_list = all_bifurcations_list + bifurcation_list
    all_terminations_list = all_terminations_list + termination_list
    total_dendritic_length = 0

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
    
    if i >= 0 and i < 4:
#        mc.plot_branches(branch_list)
        mc.plot_bifurcations(bifurcation_list)