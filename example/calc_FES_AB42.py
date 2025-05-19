#!/usr/bin/python3

import numpy as np
import sys
import glob

from DRIDmetric import DRID

import pyemma
import pyemma.coordinates as pyc

sys.path.append("/home/moritz/Software/scripts/freeEnergyCalculation/")
import freenet


############################
###   Input Parameters   ###
############################

sysname = "DRID_AB42"
drid_dir = "./data/"

cluster_cutoff = 0.02

############################
###     Clustering       ###
############################
"""
Clusters the state trajectories and calculates the corresponding Trainsition Matrix
"""

M = freenet.clusterStates(
    prefix = sysname, 
    directory = drid_dir, 
    cutoff = cluster_cutoff
)


############################
###   branching Prob.    ###
############################
"""
Calculate equilibrium probabilities of each minimum/state,
a  s well branching probabilities for transitions between connected state.
"""
state_probabilities, probability_Matrix = freenet.calcProbabilities(
    tm=M, sysname=sysname, save=True, remove=True
)

############################
###         FES          ###
############################
"""
Calculate free energies of minimas and transition states (barrier hight)
and output them in in pathsample style
"""
gE = freenet.getEnergy(
    state_probabilities,
    probability_Matrix,
    temperature=300,
    timescale=20e-12,
)

gE.run(mts=True)


############################
###        TEST          ###
############################

from legacy.getEnergy import getEnergy as getEnergy_legacy

n_min = len(state_probabilities)

legacy = getEnergy_legacy(
    "./legacy/DRID_AB42_equiProb.txt",
    "./legacy/DRID_AB42_branchingProbabilities.txt",
    n_min,
    temperature=300,
    timescale=20e-12,
    mts=True,
    dic="./legacy/DRID_AB42_TM-Dict.txt",
)

legacy_state_probabilities = legacy.loadEqui("./legacy/DRID_AB42_equiProb.txt")
legacy_probability_Matrix = legacy.loadMatrix(
    "./legacy/DRID_AB42_branchingProbabilities.txt"
)

dif_M = probability_Matrix - legacy_probability_Matrix
dif_eq = state_probabilities - legacy_state_probabilities

print("---------------------------------------------------------")
print("Difference between legacy matrix: ", np.sum(np.abs(dif_M)))
print("Difference between legacy eqi. prob.: ", np.sum(np.abs(dif_M)))
print("---------------------------------------------------------")
print("---------------------------------------------------------")
