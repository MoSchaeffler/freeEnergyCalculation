# freenet

`freenet` is a Python package to cluster molecular dynamics trajectories into discrete states and calculate the free energy surface (FES) based on the resulting states. It is designed to provide tools for understanding conformational landscapes, equilibrium probabilities, and transition state free energies.

---

## Theory

The free energy surface (FES) of a biomolecule encodes its structural and dynamical properties. To derive the FES from molecular dynamics (MD) simulations, the trajectory is discretized into microstates by clustering in a suitable metric space (e.g., DRID). Each state corresponds to a free-energy minimum.

The free energy of state $i$ is computed from its equilibrium probability $p_i$ as

$$
F_i = -k_B T \log(p_i),
$$

where $k_B$ is the Boltzmann constant and $T$ the temperature.

Transition state free energies between minima $j$ and $k$ are obtained from the Eyring–Polanyi formulation:

$$
F_{jk} = F_k - k_B T \log(k_{jk}) + k_B T \log\left(\frac{k_B T}{h}\right),
$$

where $k_{jk}$ is the transition rate and $h$ is Planck’s constant. To minimize finite-sampling errors, forward and backward estimates are averaged:

$$
F^{ts}_{jk} = \frac{F_{jk} + F_{kj}}{2}.
$$

This framework ensures that the equilibrium distribution and kinetic rates are faithfully reproduced on the network level.

---

## Installation

Clone the repository and install with pip:

```bash
git clone https://github.com/MoSchaeffler/freeEnergyCalculation.git
cd freenet
pip install .
```

Dependencies include:
- numpy
- networkx
- deeptime
- MDAnalysis
- PyYAML

---

## Key Functions and Classes

### `clusterStates` (from `clustering.py`)
Clusters trajectories into discrete states and builds the transition matrix.
Trajectories are expected to be numpy arrays of size (n_frames, n_features) or (n_frames, n_features, n_dim).

**Parameters:**
- `prefix` (str): Prefix of input trajectory *npy files.
- `directory` (str, default `"./"`): Path to data directory.
- `cutoff` (float, default `0.02`): Minimum distance cutoff for clustering.
- `max_centers` (int, default `2000`): Maximum number of clusters.
- `outname` (str, optional): Output file name for transition matrix.
- `verbose` (bool, default `False`): Print details if True.

**Returns:**  
Transition matrix as a NumPy array.

---

### `calcProbabilities` (from `branchingProbabilities.py`)
Computes equilibrium and branching probabilities.

**Parameters:**
- `tm` (np.array): Transition matrix.
- `remove` (bool, default `True`): Remove non-bidirectional transitions and disconnected states.
- `save` (bool, default `False`): Save results to file.
- `sysname` (str, default `"system"`): Prefix for saved files.

**Returns:**
- `state_probabilities`: Equilibrium probabilities of minima.
- `probability_matrix`: Transition probability matrix.

---

### `getEnergy` (main class in `getEnergy.py`)
Calculates free energies and transition states.

**Initialization Parameters:**
- `equiProb` (array): Equilibrium probabilities.
- `probMatrix` (array): Transition probability matrix.
- `temperature` (float, default `300`): Temperature in Kelvin.
- `timescale` (float, default `100e-12`): Time between MD frames (s).

**Main Methods:**
- `run(mts=False, stateDictionary=None)`: Execute workflow, write minima and transition state files.
- `checkData(prob, trans)`: Validate probability normalization.
- `rateMatrix(probMatrix, timescale)`: Compute continuous-time rate matrix.
- `freeEenergy(K, prob, T)`: Compute free energies of minima and barriers.
- `effectiveBarrier(...)`: Compute symmetrized effective barriers.


---

## Basic Usage

A typical workflow:

```python
import freenet as fn

# Step 1: Cluster trajectories
M = fn.clusterStates(prefix="traj", directory="./data/", cutoff=0.02, verbose=True)

# Step 2: Calculate equilibrium and branching probabilities
peq, probM = fn.calcProbabilities(M, remove=True, save=True, sysname="mySystem")

# Step 3: Compute free energies
energy = fn.getEnergy(peq, probM, temperature=300, timescale=100e-12)
energy.run(mts=True)
```

This produces `min.data`, `ts.data`, and related files suitable for disconnectivity graph analysis.

---

## References

- Schäffler, M., Wales, D.J., Strodel, B. *The energy landscape of Aβ42: a funnel to disorder for the monomer becomes a folding funnel for self-assembly*. ChemComm, 2024.  
- Krivov, S.V. & Karplus, M. *Hidden complexity of free energy surfaces for peptide (protein) folding*. J. Chem. Phys. **117**, 10894 (2002).
