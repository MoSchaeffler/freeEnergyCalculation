# Example: Free Energy Surface Calculation

This example demonstrates how to compute a free energy surface (FES) from molecular‐dynamics trajectories using the `calc_FES_AB42.py` script included in the `freenet` module. In this workflow we go through the following steps:

- **Clustering** of the data to define discrete microstates  
- **Transition‐matrix estimation** between microstates  
- **Free‐energy database** generation (minima and transition‐state free energies) in a format readable by PATHSAMPLE  

All steps are shown in the running running example of the `calc_FES_AB42.py` script

---

## Step by step guide

### Step 1: Prepare input data

Assume you have saved molecular dynamics trajectory descriptors (e.g., DRID vectors) as `.npy` files with a common prefix:

```
data/traj_0.npy
data/traj_1.npy
data/traj_2.npy
...
```

Each file corresponds to one trajectory segment. The entries are per-frame feature vectors of shape (n_frames, n_features) or (n_frames, n_features, n_dim).

---

### Step 2: Cluster trajectories into states

```python
import freenet as fn

# Cluster states with a cutoff of 0.02 and a maximum of 2000 clusters
M = fn.clusterStates(
    prefix="traj",
    directory="./data",
    cutoff=0.02,
    max_centers=2000,
    verbose=True
)
```

This step:
- Reads all files matching `traj*npy` in `./data/`.
- Concatenates them into one trajectory.
- Clusters them into discrete states using `RegularSpace` clustering.
- Writes per-trajectory state assignments (`state-trj_traj_i.txt`, `.npy`).
- Builds and saves the transition matrix `traj_TransitionMatrix.npy`.

---

### Step 3: Compute equilibrium and branching probabilities

```python
peq, probM = fn.calcProbabilities(
    M,
    remove=True,
    save=True,
    sysname="mySystem"
)
```

This step:
- Computes the equilibrium occupation probability of each state.
- Normalizes and symmetrizes the transition matrix.
- Optionally saves results to `mySystem_branchingProbabilities.txt` and `mySystem_equiProb.txt`.

---

### Step 4: Calculate free energies

```python
energy = fn.getEnergy(
    equiProb=peq,
    probMatrix=probM,
    temperature=300,        # Kelvin
    timescale=100e-12       # 100 ps in seconds
)

# Run the workflow with the mean-transition-state algorithm enabled
energy.run(mts=True)
```

This step:
- Validates data consistency (`checkData`).
- Converts transition probabilities to rates (`rateMatrix`).
- Computes minima free energies and transition state barriers (`freeEenergy`).
- Writes output files:
  - `min.data` : minima free energies
  - `ts.data`  : averaged transition state free energies
  - `mapping.data` : state index remapping if disconnected states are removed

These outputs are compatible with disconnectivity graph analysis tools (e.g., PATHSAMPLE / OPTIM).

---

### Full Example Script

The example script `calc_FES_AB42.py` is a working example which uses two MD trajectories, in reduced DRID space. 
For more details regarding dimensionaltity reduction via the DRID metric see: https://github.com/MoSchaeffler/DRIDmetric

The example performs all previous described steps at once similar to this:
```python
import freenet as fn

# Cluster trajectories
M = fn.clusterStates("traj", "./data", cutoff=0.02, verbose=True)

# Probabilities
peq, probM = fn.calcProbabilities(M, remove=True, save=True, sysname="mySystem")

# Free energy
energy = fn.getEnergy(peq, probM, temperature=300, timescale=100e-12)
energy.run(mts=True)
```

Run this script after preparing your trajectory data, and you will obtain the free energy landscape in `min.data` and `ts.data`.
