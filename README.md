# Block Copolymer Chain Builder for LAMMPS (Pizza.py Style)

**Status**: 🚧 Active Development  
**Author**: Claire Murphy DePompa  
**License**: [MIT License](LICENSE)  
**Citation**: See [Attribution](#attribution--citation)

This toolkit builds bead–spring polymer melts (e.g., ABA triblocks, diblocks) for use in LAMMPS simulations.  
It’s a modern, lightweight update to the classic `Pizza.py` FENE chain builder—designed for quick input generation in block copolymer DPD/coarse-grained studies.


---

## 📁 Repository Contents

### 🧱 Chain Builder Scripts (`Simple FENE Chain Builder/`)

| File               | Description                                                                                                   |
|--------------------|---------------------------------------------------------------------------------------------------------------|
| `chain_builder.py` | Main chain builder class. Creates atoms, bonds, and simulation box.                                           |
| `data.py`          | Helper for writing LAMMPS `.data` files in proper format.                                                     |
| `example.py`       | **Example driver script** – define your system and run this to create your `.data` file.                      |

### 📈 Post-Processing (`Post-Processing Analysis/`)

Includes a prototype Jupyter notebook for RDF analysis using [MDAnalysis](https://www.mdanalysis.org/).  
More tools for analyzing microphase structure are in development.

---

## 🚀 Quickstart: Creating DPD Input Configurations
Adoption and modification of these scripts is welcome and encouraged for researchers and students
who are studying block copolymers via DPD simulations. To get started, open the file
`example.py` in a text editor or IDE. At the top, look for the section labeled **USER INPUT**:

```python
# ─── USER INPUT ─────────────────────────────────────────────────────────────
nchains   = 100                    # number of chains
block_def = [(1, 30),              # 30 beads of type‑1  (block A)
             (2, 70)]              # 70 beads of type‑2  (block B)
rho_star  = 3.0                    # number density ρ*
seed      = 12345                  # RNG seed for reproducibility
outfile   = "N100_n100_AB_fA3.data"    # output filename
# ────────────────────────────────────────────────────────────────────────────
```
**Key parameters:**
* `nchains`: How many chains go in the box
* `block_def`: The layout of one polymer chain (see next section)
* `rho_star`: The number density (e.g. 3.0 for DPD)
* `seed`: A random seed for initial configuration (reuse the same seed when creating replicate experiments)
* `outfile`: File name of the .data file to create. Choose something descriptive and keep track of your experiments somewhere handy.
  * Note: These scripts include logic to print your input parameters at the top of the .data file that's created for better
    recordkeeping. So if you forget to keep track of the experiments you are running, or lose your records to the 
    ever-vigilant gremlins of science chaos, you'll have a breadcrumb trail to bring you back. ;)

Use `block_def` to describe the bead types and block lengths; each entry is a tuple with:
`(bead_type, number_of_beands_in_block)`

**Supported Polymer Types**

Feel free to copy/paste and try simulating!

| Polymer Type     | `block_def` Example           |
| ---------------- | ----------------------------- |
| Homopolymer      | `[(1, 100)]`                  |
| Diblock (A–B)    | `[(1, 40), (2, 60)]`          |
| Triblock (A–B–A) | `[(1, 15), (2, 70), (1, 15)]` |
| Multiblock (AB)₅ | `[(1,10), (2,10)] * 5`        |

### How to run the script and make your input configurations:
1. Open a terminal (Command Prompt, Anaconda Prompt, etc.)
2. Navigate to the folder with the scripts
3. Run the script: `python example.py`
4. You will then see an output in the terminal: 
   `Simulation box: 14.9 x 14.9 x 14.9`
   `Finished writing example.data`

### Output `.data` file contains:
* LAMMPS-compatible header
* Simulation box size
* Masses section
* Atoms section (atom ID, molecule ID, atom type, x, y, z)
* Bonds section (bond ID, bond type, atom1, atom2)
* Your original parameters are printed as # comment lines at the top for traceability.

### Advanced usage
Some advanced tweaks for creating multiple experiments:
* You can sweep compositions by looping over `block_def` entries
* Change bond type, bond length, or chain growth rules inside `chain_builder.py`
* Add `output/` or a full path to `outfile` to control where files are saved

Here's an example composition sweep code you can use instead of the provided `example.py`.
It loops through four compositions and outputs fours separate data files across the composition
ranges
```
nchains  = 400
n_total  = 100            # beads in ONE chain
rho_star = 3.0

for fA in [0.10, 0.30, 0.50, 0.70]:
    nA = int(round(n_total * fA))   # beads of type 1 (A)
    nB = n_total - nA               # beads of type 2 (B)
    pattern = [1]*nA + [2]*nB       # e.g. [1,1,…,2,2,…]

    cb = ChainBuilder(N=nchains * n_total, rho_star=rho_star)
    cb.build(nchains=nchains, nper=n_total, pattern=pattern)
    cb.write(f"A{nA}_B{nB}.data")   # files: A10_B90.data, A30_B70.data, …
    print(f"Created A{nA}_B{nB}.data")
```
---
## Post-processing of your experiments
The included Jupyter notebook is intended to be run on an HPC. It uses the aforementioned MDAnalysis package,
which defaults to serial execution in Jupyter.
Consult your facility's documentation for instructions on initializing and running compute nodes on your HPC.
Once you are comfortable with that process, request one node and high amounts of memory  for the default execution.
You can modify the notebook to explicitly parallelize MDAnalysis, in which case you would adjust the node/memory requests
accordingly. It's hard to provide general resource recommendations without seeing your system, but the 
guidelines below should get you in the ballpark for serial execution. It's up to you to consult your HPC's
documentation and determine how to monitor and scale your resource usage.   

**Some things to keep in mind:**
* Trajectory reading and frame iteration are typically I/O-bound, especially from ASCII LAMMPS dumps.
  * Use frame streaming for large trajectory files to reduce resource usage and avoid crashes
  * Avoid `MemoryReader` for big trajectories
* RDF calculations scale ∝ N × ⟨neighbors⟩ × F, and memory load can spike if multiple frames or selection masks are held at once.

**Rough HPC Resource Requests for Serial RDF (g(r)) Post-Processing on DPD Block Copolymer Trajectories**  

| System size (beads) | Frames | Traj Size (binary, GB) | Est. RAM Needed (GB) | Notes                                        |
| ------------------- | ------ | ---------------------- | -------------------- | -------------------------------------------- |
| 50,000              | 1,000  | \~1.5                  | 6–8                  | Fast, safe on most modern laptops            |
| 100,000             | 1,000  | \~3.0                  | 10–16                | Should still run fine in Jupyter             |
| 250,000             | 1,000  | \~7.5                  | 24–32                | Recommend >32 GB RAM to avoid slowdowns      |
| 500,000             | 1,000  | \~15.0                 | 48–64                | HPC strongly advised                         |
| 1,000,000           | 1,000  | \~30.0                 | 96–128               | Likely requires batch job, Jupyter may crash |

These assume you stream the trajectory with a frame iterator, and don’t store all pairwise distances in memory. If you store all distances, double/triple RAM.

**Options for Parallelization in MDAnalysis**

| Method                      | Parallel? | Notes                       |
| --------------------------- | --------- | --------------------------- |
| `InterRDF` from MDAnalysis  | ❌         | Fast for small systems      |
| `multiprocessing`, `joblib` | ✅         | Frame-parallel manual split |
| Dask + dask-jobqueue        | ✅         | Scalable HPC solution       |

