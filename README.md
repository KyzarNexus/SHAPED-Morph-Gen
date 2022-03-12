# SHAPED_Morph_Gen
## Description
Dendritic arbor modeling using the SHAPED algorithm. This code was adapted from the prior work of Zane Chou in modeling realistic granule cell dendritic arbors in the dentate gyrus. The application of this repository is for modeling pyramidal dendritic morphologies in the CA3 region of the rat hippocampus. Due to pyramidal cells having an apical and basal arbor, their morphological data was separated and processed separately. All source data comes from neuromorpho.org under the archive name 'Amaral'. 
## Pipeline Description
The raw .swc files downloaded from neuromorpho.org are passed through b_measure_tester.py and b_measure_conditionals.py to create distributions of morphological attributes that will be used in the morphology generation step. After running b_measure_conditionals.py, run IBI_observations.py to create a point proces file that will be used in the morphology generation step. At which point, run ppf_em_73.py to create a branching rate estimate and generate novel morphologies. After the synthetic morphology dataset has been created, run the synthetic dataset through b_measure_tester.py to assess generation accuracy to the organic dataset. 
## Raw Data & 3D Reconstruction
Each dataset can be split into three separate classes describing its morphology: compartment, branch, & morphology. Compartments are created from the raw coordinates seen in the .swc file. Branches are a collection of compartments that starts at the beginning of a bifurcation and ends at the next bifurcation (or termination). Morphologies are collections of the resulting list of branches. 
## File Description
### B_measure_tester.py
The purpose of this file is to parse through the raw .swc data, reconstruct a 3D representation of the morphology, and create distributions of morphological characteristics such as max branch length, branching angle, etc. These distributions are then stored in '[dataset]_distributions_[arbor type].txt'. To make retrieval easier in later steps, the object storing the distributions was pickled into a .pickle file. 
### B_measure_conditionals.py
The purpose of this file is similar to b_measure_conditionals, but creates distributions that are dependent on other morphological values. For example, to create a conditional distribution of branching angles with respect to distance from the soma, the script creates a distribution of branching angles for each section of distance away from the soma. The resulting conditional file, '[dataset]_conditionals_[arbor_type].txt' is integral in the morphology generation step. 
### IBI_observations.py
The purpose of this file is to create a point processes ([0 1 0 0 1 0 etc.]) that describe the branching pattern for each terminal branch. For each terminal branch, a list of compartments is tested for bifurcations. If a bifurcation occurs at a compartment, then that places in the array is '1'. Likewise, if no bufurcation occurs, it is denoted as '0'. This collection of point processes is used to create a branching rate estimate that can then be used in the morphology generation process. 
### ppf_em_73.py
'Point Process Filter / Estimation Maximization (Iteration 73)'
This file preforms multiple functions and must be outlined separately.
#### Point Process Filter
Due to the point processes created in IBI_obervations having overlapping recorded bifurcations, the resulting file must be filtered to obtain a better branching estimate in the next step. 
#### Estimation Maximization
This step is to directly create a branching rate estimate by optimizing the parameters used in the function that creates the branching rate estimate. A parameter sweep of initial conditions is performed and the set of initial conditions that converges to the best parameters to use for the rate generation are passed through. NOTE: This section takes on the order of a week to perform, depending on the speed of the computer and the size of the initial conditions passed into the EM algorithm. Due to this, is is only necessary to perform this step once to obtain the branching rate estimate that is then saved as '[dataset]-[arbor type]-lambPic.pickle'. Optimization history of this step is also saved in [dataset]-[arbor type]-opDict.pickle' & [dataset]-[arbor type]-Rel_E_Dict.pickle' accordingly
#### Morphology Generation
This next step uses the conditional, distributions, and branching rate files to create a new dataset of synthetic morphologies. 


## Final Note
This codebase is still a bit glitchy and may require some debugging to work correctly. 
