#!/usr/bin/env python3
"""
Pizza.py toolkit - data tool
Read, write, and manipulate LAMMPS data files

Usage:
    d = data("data.poly")       # read a LAMMPS data file (gzipped files accepted)
    d = data()                  # create an empty data file
    d.write("data.new")         # write out a LAMMPS data file
"""

from os import popen

class data:

    def __init__(self, *list_args):
        self.nselect = 1
        
        if len(list_args) == 0:
            self.title = "LAMMPS data file"
            self.names = {}
            self.headers = {}
            self.sections = {}
            self.comments = []
            return

        file = list_args[0]
        if file[-3:] == ".gz":
            f = popen("%s -c %s" % (PIZZA_GUNZIP, file), 'r')
        else:
            f = open(file, 'r')
        
        self.title = f.readline().rstrip("\n")
        self.names = {}
        headers = {}
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if len(line) == 0:
                continue
            found = False
            for keyword in hkeywords:
                if keyword in line:
                    found = True
                    words = line.split()
                    if keyword in ["xlo xhi", "ylo yhi", "zlo zhi"]:
                        headers[keyword] = (float(words[0]), float(words[1]))
                    elif keyword == "xy xz yz":
                        headers[keyword] = (float(words[0]), float(words[1]), float(words[2]))
                    else:
                        headers[keyword] = int(words[0])
                    break
            if not found:
                break

        sections = {}
        while True:
            found = False
            for pair in skeywords:
                keyword, length = pair[0], pair[1]
                if keyword == line.strip():
                    found = True
                    if length not in headers:
                        raise Exception("data section %s has no matching header value" % keyword)
                    # Read and discard one line (possibly a blank line)
                    f.readline()
                    lines = []
                    for i in range(headers[length]):
                        lines.append(f.readline())
                    sections[keyword] = lines
                    break
            if not found:
                raise Exception("invalid section %s in data file" % line.strip())
            # Read the blank line between sections
            f.readline()
            line = f.readline()
            if not line:
                break
            line = line.strip()
        
        f.close()
        self.headers = headers
        self.sections = sections

    def map(self, *pairs):
        if len(pairs) % 2 != 0:
            raise Exception("data map() requires pairs of mappings")
        for i in range(0, len(pairs), 2):
            j = i + 1
            self.names[pairs[j]] = pairs[i] - 1

    def get(self, *list_args):
        if len(list_args) == 1:
            field = list_args[0]
            array = []
            lines = self.sections[field]
            for line in lines:
                words = line.split()
                values = list(map(float, words))
                array.append(values)
            return array
        elif len(list_args) == 2:
            field = list_args[0]
            n = list_args[1] - 1
            vec = []
            lines = self.sections[field]
            for line in lines:
                words = line.split()
                vec.append(float(words[n]))
            return vec
        else:
            raise Exception("invalid arguments for data.get()")

    def reorder(self, name, *order):
        natoms = len(self.sections[name])
        oldlines = self.sections[name]
        newlines = ["" for _ in range(natoms)]
        for index in order:
            for i in range(natoms):
                words = oldlines[i].split()
                newlines[i] += words[index - 1] + " "
        for i in range(natoms):
            newlines[i] += "\n"
        self.sections[name] = newlines

    def replace(self, name, icol, vector):
        lines = self.sections[name]
        newlines = []
        j = icol - 1
        for i in range(len(lines)):
            words = lines[i].split()
            words[j] = str(vector[i])
            newline = ' '.join(words) + '\n'
            newlines.append(newline)
        self.sections[name] = newlines

    def newxyz(self, dm, ntime):
        nsnap = dm.findtime(ntime)
        dm.sort(ntime)
        x, y, z = dm.vecs(ntime, "x", "y", "z")
        self.replace("Atoms", self.names['x'] + 1, x)
        self.replace("Atoms", self.names['y'] + 1, y)
        self.replace("Atoms", self.names['z'] + 1, z)
        if "ix" in dm.names and "ix" in self.names:
            ix, iy, iz = dm.vecs(ntime, "ix", "iy", "iz")
            self.replace("Atoms", self.names['ix'] + 1, ix)
            self.replace("Atoms", self.names['iy'] + 1, iy)
            self.replace("Atoms", self.names['iz'] + 1, iz)

    def delete(self, keyword):
        if keyword in self.headers:
            del self.headers[keyword]
        elif keyword in self.sections:
            del self.sections[keyword]
        else:
            raise Exception("keyword not found in data object")

    def write(self, file):
        with open(file, "w") as f:
            print(self.title, file=f)
            # --- user‑supplied comments (each becomes '# ...') ---------------
            for line in getattr(self, "comments", []):
                print(f"# {line}", file=f)
            if getattr(self, "comments", []):
                print(file=f)  # blank line after the comment block
            # -----------------------------------------------------------------

            for keyword in hkeywords:
                if keyword in self.headers:
                    if keyword in ["xlo xhi", "ylo yhi", "zlo zhi"]:
                        pair = self.headers[keyword]
                        print(f"{pair[0]} {pair[1]} {keyword}", file=f)
                    elif keyword == "xy xz yz":
                        triple = self.headers[keyword]
                        print(f"{triple[0]} {triple[1]} {triple[2]} {keyword}", file=f)
                    else:
                        print(f"{self.headers[keyword]} {keyword}", file=f)
            # Ensure "Masses" section exists before writing
            if "Masses" not in self.sections:
                self.sections["Masses"] = [
                    f"{atype} 1.0\n" for atype in range(1, self.headers["atom types"] + 1)
                ]
            
            # Write "Masses" section before other sections
            print("\nMasses\n", file=f)
            for line in self.sections["Masses"]:
                f.write(line)
            
            for pair in skeywords:
                keyword = pair[0]
                if keyword in self.sections:
                    print("\n%s\n" % keyword, file=f)
                    for line in self.sections[keyword]:
                        f.write(line)
    
    def iterator(self, flag):
        if flag == 0:
            return 0, 0, 1
        return 0, 0, -1

    def findtime(self, n):
        if n == 0:
            return 0
        raise Exception("no step %d exists" % n)

    def viz(self, isnap):
        if isnap:
            raise Exception("cannot call data.viz() with isnap != 0")
        
        id_idx = self.names["id"]
        type_idx = self.names["type"]
        x_idx = self.names["x"]
        y_idx = self.names["y"]
        z_idx = self.names["z"]

        xlohi = self.headers["xlo xhi"]
        ylohi = self.headers["ylo yhi"]
        zlohi = self.headers["zlo zhi"]
        box = [xlohi[0], ylohi[0], zlohi[0], xlohi[1], ylohi[1], zlohi[1]]

        atoms = []
        atomlines = self.sections["Atoms"]
        for line in atomlines:
            words = line.split()
            atoms.append([int(words[id_idx]), int(words[type_idx]),
                          float(words[x_idx]), float(words[y_idx]), float(words[z_idx])])
        bonds = []
        if "Bonds" in self.sections:
            bondlines = self.sections["Bonds"]
            for line in bondlines:
                words = line.split()
                bid, btype = int(words[0]), int(words[1])
                atom1, atom2 = int(words[2]), int(words[3])
                atom1words = atomlines[atom1 - 1].split()
                atom2words = atomlines[atom2 - 1].split()
                bonds.append([bid, btype,
                              float(atom1words[x_idx]), float(atom1words[y_idx]), float(atom1words[z_idx]),
                              float(atom2words[x_idx]), float(atom2words[y_idx]), float(atom2words[z_idx]),
                              float(atom1words[type_idx]), float(atom2words[type_idx])])
        tris = []
        lines = []
        return 0, box, atoms, bonds, tris, lines

    def maxbox(self):
        xlohi = self.headers["xlo xhi"]
        ylohi = self.headers["ylo yhi"]
        zlohi = self.headers["zlo zhi"]
        return [xlohi[0], ylohi[0], zlohi[0], xlohi[1], ylohi[1], zlohi[1]]

    def maxtype(self):
        return self.headers["atom types"]

# Keywords used in the data file
hkeywords = ["atoms", "ellipsoids", "lines", "triangles", "bodies",
             "bonds", "angles", "dihedrals", "impropers",
             "atom types", "bond types", "angle types", "dihedral types",
             "improper types", "xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"]

skeywords = [["Masses", "atom types"],
             ["Atoms", "atoms"],
             ["Ellipsoids", "ellipsoids"],
             ["Lines", "lines"],
             ["Triangles", "triangles"],
             ["Bodies", "bodies"],
             ["Bonds", "bonds"],
             ["Angles", "angles"],
             ["Dihedrals", "dihedrals"],
             ["Impropers", "impropers"],
             ["Velocities", "atoms"],
             ["Pair Coeffs", "atom types"],
             ["Bond Coeffs", "bond types"],
             ["Angle Coeffs", "angle types"],
             ["Dihedral Coeffs", "dihedral types"],
             ["Improper Coeffs", "improper types"],
             ["BondBond Coeffs", "angle types"],
             ["BondAngle Coeffs", "angle types"],
             ["MiddleBondTorsion Coeffs", "dihedral types"],
             ["EndBondTorsion Coeffs", "dihedral types"],
             ["AngleTorsion Coeffs", "dihedral types"],
             ["AngleAngleTorsion Coeffs", "dihedral types"],
             ["BondBond13 Coeffs", "dihedral types"],
             ["AngleAngle Coeffs", "improper types"],
             ["Molecules", "atoms"],
             ["Tinker Types", "atoms"]]
