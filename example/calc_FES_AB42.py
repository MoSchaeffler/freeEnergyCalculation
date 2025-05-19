#!/usr/bin/python3

import freenet


############################
###   Input Parameters   ###
############################

prefix = "DRID_AB42"
drid_dir = "./data/"

cluster_cutoff = 0.02

############################
###     Clustering       ###
############################
"""
Clusters the state trajectories and calculates the corresponding Trainsition Matrix
"""

M = freenet.clusterStates(
    prefix = prefix, 
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
    tm=M, sysname=prefix, save=True, remove=True
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

