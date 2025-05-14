#!/usr/bin/python3

import numpy as np

verbose = False


class calcProbabilities:
    """
    Calculates the rquilibrium probabilities of each minimum, 
    as well branching probabilities for transitions between connected state

    Parameters
    ----------
    traj : Trajectory
        MD trajectory object (e.g., from MDTraj or MDAnalysis).
    cutoff : float
        Clustering cutoff distance.
    temperature : float
        Simulation temperature in Kelvin.

    Attributes
    ----------
    nodes : list of Node
        Identified minima.
    edges : list of Edge
        Transition-state connections.
    """

    def __init__(self, sysname, tm_path, remove=False):

        # Load transition matrix
        self.tm = np.load(tm_path)
        # nstates
        self.nmin = len(self.tm[0])

        ############################
        ###  equi. Probabilities ###
        ############################

        # population of each node by summing over all incoming flow
        pop_m = np.sum(self.tm, axis=0)
        # probabiliy
        peq_m = pop_m / np.sum(pop_m)

        ############################
        ###   Branching Prob.    ###
        ############################

        # norm matrix so that rows sum to 0 -> right stochastic matrix
        rm = np.zeros((self.nmin, self.nmin))

        # if mts == True delete transitions that are not bidirectional and delete disconnected states
        if remove == False:
            T = self.tm

        if remove == True:
            print("Processing Transition Matrix")
            T = np.zeros((self.nmin, self.nmin))
            for i in range(self.nmin):
                for j in range(self.nmin):

                    t = self.tm[i, j]

                    if i == j:
                        T[i, j] = t
                    else:
                        if self.tm[j, i] != 0:
                            T[i, j] = t

            print("Removing non bidirectional transition states")
            i = 0
            s = 1
            while i < self.nmin:

                l = np.count_nonzero(T[i])

                if l == 1:

                    T = np.delete(T, i, axis=0)
                    T = np.delete(T, i, axis=1)

                    self.nmin -= 1

                    if verbose == True:
                        print("removed state {}".format(s))

                if l != 1:
                    i += 1

                s += 1

            print("Done")

        """
        for i in range(self.nmin):
            row = T[i]
            n = np.sum(row)
            rm[i] = row/n
        """
        i = 0
        delEntries = []

        print("---------------------------------------------------------")
        print("Start Normalization")
        rm = np.zeros((self.nmin, self.nmin))
        print("Dimension Probability Matrix: ", rm.shape)
        while i < self.nmin:
            if i == 0:
                rm = np.zeros((self.nmin, self.nmin))
                if verbose == True:
                    print("Dimension Probability Matrix: ", rm.shape)

            row = T[i]
            n = np.sum(row)
            if n == 0:
                if verbose == True:
                    print("state {} has no outgoing transitions".format(i))
                    print("delete entry and reinitialize")
                T = np.delete(T, i, 0)
                T = np.delete(T, i, 1)
                delEntries.append(i)
                i = 0
                self.nmin = self.nmin - 1
                # update prob eq list
                peq_m = peq_m[:-1]

            else:
                rm[i] = row / n
                # print(np.sum(rm[i]))
                i += 1

        print("End Normalization")
        print("Dimension Probability Matrix: ", rm.shape)

        # write everyting
        self.writeProbMatrix(rm, sysname + "_branchingProbabilities.txt")

        self.writeProbabilities(peq_m, sysname + "_equiProb.txt")

    ############################
    ###   Helper Functions   ###
    ############################

    def writeProbMatrix(self, T, outname):

        with open(outname, "w") as f:

            for i in range(self.nmin):
                row = str(T[i, 0])
                for j in range(1, self.nmin):
                    row += "\t " + str(T[i, j])

                row += "\n"
                f.write(row)

    def writeProbabilities(self, P, outname):

        with open(outname, "w") as f:

            prob = str(P[0])
            for i in range(1, self.nmin):
                prob += "\t" + str(P[i])

            prob += "\n"
            f.write(prob)
