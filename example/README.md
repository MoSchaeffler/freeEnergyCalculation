# Example: Free Energy Surface Calculation

This example demonstrates how to compute a free energy surface (FES) from molecular‐dynamics trajectories using the `calc_FES.py` script included in the `freenet` module. In this workflow we:

1. **Clustering** of the data to define discrete microstates  
2. **Transition‐matrix estimation** between microstates  
3. **Free‐energy database** generation (minima and transition‐state free energies) in a format readable by PATHSAMPLE  

---

## Data Format

- **Input**: NumPy arrays of shape `(N_f, N_s)`, where `N_f` = number of frames and `N_s` = reduced dimensionality (DRID coordinates).  
- In this example, two MD trajectories of the Aβ-42 peptide have been processed with the DRID package to produce low-dimensional representations for each frame.  
- Files must share a common prefix (e.g. `DRID_AB42`).

For more details regarding dimensionaltity reduction via the DRID metric see: https://github.com/MoSchaeffler/DRIDmetric