# This file is part of the Block Copolymer Chain Builder
# Copyright (C) 2025 Claire Murphy DePompa
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# chain_builder.py
import math, random
from data import data

# --- deterministic Park–Miller RNG for old runs -----------------------------
IM, AM, IA, IQ, IR = 2147483647, 1.0/2147483647, 16807, 127773, 2836
def pm_rand(state):
    k = state // IQ
    state = IA*(state - k*IQ) - IR*k
    if state < 0: state += IM
    return AM*state, state
# ---------------------------------------------------------------------------

class ChainBuilder:
    def __init__(self, N, rho_star, xasp=1, yasp=1, zasp=1, rng="pm"):
        self.N, self.rho_star = N, rho_star
        self.seed = 12345
        self.rng_kind = rng        # "pm" or "py"
        self.mtype, self.btype = 1, 1
        self.blen, self.dmin = 0.97, 1.02
        self.id = "chain"
        self.atoms, self.bonds = [], []

        # --- one canonical box calculation (was duplicated) ----------------
        volume = N / rho_star
        side = (volume / (xasp*yasp*zasp))**(1/3)
        self.xprd, self.yprd, self.zprd = xasp*side, yasp*side, zasp*side
        self.xlo, self.ylo, self.zlo = -self.xprd/2, -self.yprd/2, -self.zprd/2
        self.xhi, self.yhi, self.zhi =  self.xprd/2,  self.yprd/2,  self.zprd/2
        # -------------------------------------------------------------------
        print(f"Simulation box: {self.xprd:.3g} × {self.yprd:.3g} × {self.zprd:.3g}")

    # -----------------------------------------------------------------------
    def _rand(self):
        if self.rng_kind == "py":
            return random.random()
        r, self.seed = pm_rand(self.seed)
        return r
    # -----------------------------------------------------------------------

    def _pbc(self, x, y, z):
        if x < self.xlo: x += self.xprd
        elif x >= self.xhi: x -= self.xprd
        if y < self.ylo: y += self.yprd
        elif y >= self.yhi: y -= self.yprd
        if z < self.zlo: z += self.zprd
        elif z >= self.zhi: z -= self.zprd
        return x, y, z

    # ------------- unified build routine -----------------------------------
    def build(self, nchains, nper, pattern=None):
        if pattern is None:
            pattern = [self.mtype]*nper
        if len(pattern) != nper:
            raise ValueError("len(pattern) must equal nper")
        for ic in range(nchains):
            id_atom_prev = self.atoms[-1][0] if self.atoms else 0
            id_bond_prev = self.bonds[-1][0] if self.bonds else 0
            mol_id_start = self.atoms[-1][1] if self.atoms else 0

            for i in range(nper):
                idatom = id_atom_prev + i + 1
                atype  = pattern[i]
                if i == 0:                     # first bead: random in box
                    x = self.xlo + self._rand()*self.xprd
                    y = self.ylo + self._rand()*self.yprd
                    z = self.zlo + self._rand()*self.zprd
                else:                          # grow freely‑jointed chain
                    rsq = 2.0
                    while rsq > 1.0:
                        dx = 2*self._rand()-1; dy = 2*self._rand()-1; dz = 2*self._rand()-1
                        rsq = dx*dx + dy*dy + dz*dz
                    r = math.sqrt(rsq)
                    dx, dy, dz = dx/r, dy/r, dz/r
                    x = self.atoms[-1][3] + dx*self.blen
                    y = self.atoms[-1][4] + dy*self.blen
                    z = self.atoms[-1][5] + dz*self.blen
                x, y, z = self._pbc(x, y, z)
                mol_id = mol_id_start + ic + 1 if self.id == "chain" else i+1
                self.atoms.append([idatom, mol_id, atype, x, y, z])
                if i:
                    self.bonds.append([id_bond_prev+i, self.btype, idatom-1, idatom])
    # -----------------------------------------------------------------------

    def write(self, fname):
        if len(self.atoms) != self.N:
            raise RuntimeError(f"{len(self.atoms)} beads generated; expected {self.N}")
        atypes = max(a[2] for a in self.atoms)
        btypes = max((b[1] for b in self.bonds), default=0)

        d = data()
        d.title = "LAMMPS FENE chain data file"
        d.headers = {"atoms": len(self.atoms), "bonds": len(self.bonds),
                     "atom types": atypes, "bond types": btypes,
                     "xlo xhi": (self.xlo, self.xhi),
                     "ylo yhi": (self.ylo, self.yhi),
                     "zlo zhi": (self.zlo, self.zhi)}
        d.sections["Masses"] = [f"{i+1} 1.0\n" for i in range(atypes)]
        d.sections["Atoms"]  = [f"{i} {j} {k} {x:.6f} {y:.6f} {z:.6f}\n"
                                for i,j,k,x,y,z in self.atoms]
        d.sections["Bonds"]  = [f"{i} {t} {a} {b}\n" for i,t,a,b in self.bonds]
        d.write(fname)
