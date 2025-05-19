#!/usr/bin/python3

from typing import Dict, Callable, Optional
import numpy as np

verbose = False


def calcProbabilities(
    tm: np.array,
    remove: bool = True,
    save: bool = False,
    sysname: str = "system",
) -> (np.array, np.array):
    """
    Computes the equilibrium probabilities of each minimum,
    as well branching probabilities for transitions between connected state

    Parameters
    ----------
    tm : Transition Matrix
        2D numpy object containing number of transitions between each state
    sysname : str
        System name for saving
    remove : bool
        use algorithm removing non biderctional transitions between states
    save: bool
        save branching probabilities and equilibrium state probabilities to file

    Returns
    -------
    state_probabilities : ndarray of shape (n_states,)
        The euilibrium probability of occupying each free‐energy minimum.
    probability_matrix : ndarray of shape (n_states, n_states)
        The transition‐probability matrix between minima, where
        entry `(i, j)` is the probability of moving from state `i` to `j`.
    """

    # n states
    n_min = len(tm[0])

    ############################
    ###  equi. Probabilities ###
    ############################

    # population of each node by summing over all incoming flow
    pop_m = np.sum(tm, axis=0)
    # normalize
    peq_m = pop_m / np.sum(pop_m)

    ############################
    ###   Branching Prob.    ###
    ############################

    # norm matrix so that rows sum to 0 -> right stochastic matrix
    rm = np.zeros((n_min, n_min))

    # if mts == True delete transitions that are not bidirectional and delete disconnected states
    if remove == False:
        T = tm

    if remove == True:
        print("Processing Transition Matrix")
        T = np.zeros((n_min, n_min))
        for i in range(n_min):
            for j in range(n_min):

                t = tm[i, j]

                if i == j:
                    T[i, j] = t
                else:
                    if tm[j, i] != 0:
                        T[i, j] = t

        print("Removing non bidirectional transition states")
        i = 0
        s = 1
        while i < n_min:

            l = np.count_nonzero(T[i])

            if l == 1:

                T = np.delete(T, i, axis=0)
                T = np.delete(T, i, axis=1)

                n_min -= 1

                if verbose == True:
                    print("removed state {}".format(s))

            if l != 1:
                i += 1

            s += 1

        print("Done")

    """
    for i in range(n_min):
        row = T[i]
        n = np.sum(row)
        rm[i] = row/n
    """
    i = 0
    delEntries = []

    print("---------------------------------------------------------")
    print("Start Normalization")
    rm = np.zeros((n_min, n_min))
    print("Dimension Probability Matrix: ", rm.shape)
    while i < n_min:
        if i == 0:
            rm = np.zeros((n_min, n_min))
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
            n_min = n_min - 1
            # update prob eq list
            peq_m = peq_m[:-1]

        else:
            rm[i] = row / n
            # print(np.sum(rm[i]))
            i += 1

    print("End Normalization")
    print("Dimension Probability Matrix: ", rm.shape)

    # write everyting
    if save == True:
        writeProbMatrix(rm, sysname + "_branchingProbabilities.txt")

        writeProbabilities(peq_m, sysname + "_equiProb.txt")

    return peq_m, rm


############################
###   Helper Functions   ###
############################


def writeProbMatrix(M: np.array, outname: str) -> None:
    """
    Write matrix of branching probabilities to .txt file

    Parameters
    ----------
    M : np.array
        matrix of branching probabilities
    outname: str
        name of out file
    """
    nmin = len(M)

    with open(outname, "w") as f:

        for i in range(nmin):
            row = str(M[i, 0])
            for j in range(1, nmin):
                row += "\t " + str(M[i, j])

            row += "\n"
            f.write(row)


def writeProbabilities(P, outname):
    """
    Write minima equilibrium probabilities to .txt file

    Parameters
    ----------
    P : np.array
        matrix of brequilibrium probabilities
    outname: str
        name of out file
    """
    nmin = len(P)

    with open(outname, "w") as f:

        prob = str(P[0])
        for i in range(1, nmin):
            prob += "\t" + str(P[i])

        prob += "\n"
        f.write(prob)
