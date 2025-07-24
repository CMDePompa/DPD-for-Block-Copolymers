"""
EXAMPLE DRIVER SCRIPT
________________________________________________________________________________________________________________________
 - Builds an AB diblock copolymer melt* for Dissipative‑Particle‑Dynamics
   (or any bead‑spring) simulations in LAMMPS.
 - You choose:
       1.  number of chains in the box (`nchains`)
       2.  block layout per chain (`block_def` --> edit this list to customize)
       3.  target number density (`rho_star`)
       4.  output file name (`outfile`)
 - Script writes a single compressed‑format LAMMPS data file (e.g. “ab.data”)
   containing the box, atoms, bonds, and Masses section.
________________________________________________________________________________________________________________________
 HOW TO EDIT FOR A NEW SYSTEM:
 - Change only the values in the USER INPUT block, leave the rest
 - Change the tuples in `block_def`:
       (atom_type, number_of_beads_in_block)
   Example below makes an **A15‑B70‑A15** triblock (A = type‑1, B = type‑2).
 - Keep `nchains * total_beads_per_chain` reasonable (less than 50000) for quick tests.
"""
import os
from chain_builder import ChainBuilder
from data          import data     # <- patched version with d.comments support

# ─── USER INPUT ─────────────────────────────────────────────────────────────
nchains   = 100                    # number of chains
block_def = [(1, 30),              # 35 beads of type‑1  (block A)
             (2, 70)]              # 70 beads of type‑2  (block B)
rho_star  = 3.0                    # number density ρ*
seed      = 12345                  # RNG seed for reproducibility
outfile   = "N100_n100_AB_fA3.data"    # output filename
# ────────────────────────────────────────────────────────────────────────────


# DO NOT EDIT THE CODE BELOW
# Expand block_def --> full per‑bead pattern list, e.g. [1,1,…,2,2,…,1,1,…]
pattern = [atype for atype, L in block_def for _ in range(L)]

# Build the system
builder = ChainBuilder(N=nchains * len(pattern), rho_star=rho_star)
builder.seed = seed
builder.build(nchains=nchains, nper=len(pattern), pattern=pattern)

# Assemble a data() object and attach build parameters as comments
d = data()
d.title = f"LAMMPS FENE chain data : {os.path.basename(outfile)}"
d.comments = [
    f"nchains   = {nchains}",
    f"block_def = {block_def}",
    f"rho*      = {rho_star}",
    f"seed      = {seed}",
]

# Required header dictionary
d.headers = {
    "atoms":       len(builder.atoms),
    "bonds":       len(builder.bonds),
    "atom types":  max(a[2] for a in builder.atoms),
    "bond types":  1,
    "xlo xhi": (builder.xlo, builder.xhi),
    "ylo yhi": (builder.ylo, builder.yhi),
    "zlo zhi": (builder.zlo, builder.zhi),
}

# Masses, Atoms, Bonds sections
d.sections["Masses"] = ["1 1.0\n", "2 1.0\n"]
d.sections["Atoms"]  = [
    f"{aid} {mol} {atype} {x:.6f} {y:.6f} {z:.6f}\n"
    for aid, mol, atype, x, y, z in builder.atoms
]
d.sections["Bonds"]  = [
    f"{bid} 1 {a1} {a2}\n"
    for bid, _, a1, a2 in builder.bonds
]

# Write the file
d.write(outfile)
print("Finished writing", outfile)