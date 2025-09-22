#!/usr/bin/python3

from typing import Dict, Callable, Optional
import numpy as np
import sys
import glob
from pathlib import Path

import pickle as pkl

#import pyemma
#import pyemma.coordinates as pyc

from deeptime.clustering import RegularSpace

def clusterStates(
    prefix: str,
    directory: Optional[str] = "./",
    cutoff: Optional[float] = 0.02,
    max_centers: Optional[int] = 2000,
    outname: Optional[str] = False,
    verbose: Optional[bool] = False
) -> np.array:
    """
    description

    Parameters
    ----------
    directory : str
        path to data directory
    prefix: str
        prefix of all state trajectories
    cutoff: float
        cutoff for clusteriong in statespace
    max_centers: int
        maximum number of clusters
    outname : str
        name of output matrix

    Returns
    -------
    None:
        Saves all data directly as numpy files
    """

    # get data
    pattern = str(Path(directory) / f"{prefix}_*")
    files = glob.glob(pattern)
    nf = len(files)

    # construct trajectory of all states
    traj_full = np.load(files[0])
    for i in range(1,nf):
        traj_full = np.concatenate(
            (traj_full, np.load(files[i]))
        )

    np.save(str(Path(directory) / f"fullTrajectory_{prefix}.npy"),traj_full)

    #### Cluster full dataset
    # stack dimensions into 2d array of size (n_frames, n_dims)
    n_frames = traj_full.shape[0]
    data = traj_full.reshape(n_frames,-1)

    print(traj_full.shape,data.shape)
    print(traj_full[0])
    print(data[0])

    model = RegularSpace(dmin=cutoff, max_centers=max_centers)
    c_rS = model.fit(data).fetch_model()

    # save model
    with open(f'cluster_model_{prefix}.pkl', 'wb') as f:
        pkl.dump(c_rS, f, protocol=pkl.HIGHEST_PROTOCOL)

    N_cl = c_rS.n_clusters

    if verbose:
        print("Number of Clusters: {}".format(c_rS.n_clusters))


    ########## assign original data and build Transition Matrix
    # loading the pyemma mapping for each trajectory individually avoids artifacts due to concatinating
    M = np.zeros((N_cl, N_cl))

    for i in range(nf):
        traj = np.load(files[i])
        n_frames = traj.shape[0]
        data = traj.reshape(n_frames,-1)

        t = c_rS.transform(data)

        # add to transition matrix
        transitions = zip(t[:-1], t[1:])

        for trans in transitions:
            # print(trans)
            M[trans[0], trans[1]] += 1

        # write state trj as text file 
        with open(str(Path(directory) / f"state-trj_{prefix}_{i}.txt"), "w") as f:

            f.write(f"@ {files[i]}\n")
            f.write("@ frame\t state\n")

            for j in range(len(t)):
                f.write("{}\t {}\t\n".format(j, t[j]))

        # save as np.array
        np.save(str(Path(directory) / f"state-trj_{prefix}_{i}.npy"), t)


    # save transition matrix
    if outname == False:
        np.save(f"{prefix}_TransitionMatrix.npy", M)
    else:
        np.save(outname, M)

    """
    # write dict
    with open(sysname + "_TM-Dict.txt", "w") as f:

        f.write("@ Transition Matrix Dictionary\n")
        f.write("@ state population\n")

        for i in range(N_cl):
            f.write("{}\t {}\t\n".format(i, i))
    """
    return M