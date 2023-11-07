from mpi4py import MPI
from neuron import h
pc = h.ParallelContext() # Need to import for parallel functions

id_pc = int(pc.id()) # Get id of the core
nhost = int(pc.nhost()) # Get the total number of cores being used

print("I am %i out of %i"%(id_pc,nhost))#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 19:44:39 2019

@author: zzchou
"""

