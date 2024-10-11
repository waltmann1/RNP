from __future__ import division
import numpy as np
from Quaternion import QuaternionBetween
import numpy.linalg as la


class CGMolyAbs(object):

    def __init__(self):

        self.positions = []
        self.f_weights = []
        self.x_weights = []
        self.name = None
        self.bonds = None
        self.site_indexes = []
        self.domain_starts = []
        self.bond_types = []
        self.atom_types = None
        self.ucg_atom_types = None
        self.angles = None
        self.angle_types = None

    def write_data_file(self, name=None):

        if self.atom_types is None:
            self.set_atom_types()
        if name is None:
            name = self.name + ".data"

        f = open(name, "w")
        f.write("LAMMPS data file via mscg, version 20 Sep 2021, timestep = 0\n\n")
        f.write(str(len(self.site_indexes)) + " atoms\n")
        f.write(str(len(self.bonds)) + " bonds\n")
        if self.angles is not None:
            f.write(str(len(self.angles)) + " angles\n")
        else:
            f.write("0 angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n")
        f.write("\n")


        f.write(str(self.n_atom_types()) + " atom types\n")
        f.write(str(self.n_bond_types()) + " bond types\n")
        f.write(str(self.n_angle_types()) + " angle types\n")
        f.write("\n")

        f.write("0 2000 xlo xhi\n")
        f.write("0 2000 ylo yhi\n")
        f.write("0 2000 zlo zhi\n")
        f.write("\n")

        f.write("Masses\n\n")
        for t in self.get_unique_atoms():
            f.write(str(int(t[0])) + " " + str(t[1]) + "\n")
        f.write("\n")

        f.write("Bond Coeffs #harmonic K R0\n\n")
        uni = self.get_unique_bonds()
        for t in uni:
            f.write(str(int(t[0])) + " " + str(t[2]) + " " + str(t[1]) + "\n")
        f.write("\n")

        if self.angles is not None:
            f.write("Angle Coeffs #harmonic K (energy/radian^2) theta_0 (degrees)\n\n")
            uni = self.get_unique_angles()
            for t in uni:
                f.write(str(int(t[0])) + " " + str(t[2]) + " " + str(t[1]) + "\n")
            f.write("\n")

        f.write("Atoms #full\n\n")
        for i in range(len(self.positions)):
            f.write(str(i + 1) + " 1 " + str(self.atom_types[i][0]) + " 0 " + str(self.positions[i][0]) + " " + str(self.positions[i][1]) +
                    " " + str(self.positions[i][2]) + "\n")
        f.write("\n")

        f.write("Bonds\n\n")
        for i in range(len(self.bonds)):
            f.write(str(i+1) + " " + str(int(self.bond_types[i][0])) + " " + str(int(self.bonds[i][0])) + " " + str(int(self.bonds[i][1])) + "\n")
        f.write("\n")

        if self.angles is not None:
            f.write("Angles\n\n")
            for i in range(len(self.angles)):
                f.write(str(i + 1) + " " + str(int(self.angle_types[i][0])) + " " + str(int(self.angles[i][0])) + " " + str(
                    int(self.angles[i][1])) + " " + str(int(self.angles[i][2])) + "\n")
            f.write("\n")

    def write_yaml(self, name=None):

        if name is None:
            name = self.name
        f = open(name + ".yaml", 'w')

        f.write("site-types:\n")

        for ind in range(len(self.site_indexes)):
            f.write("    " + self.name + "_" + str(ind + 1) + ":\n")
            f.write("        " + "index:    " + str(self.site_indexes[ind]) + "\n")
            f.write("        " + "x-weight: " + str(self.x_weights[ind]) + "\n")
            f.write("        " + "f-weight: " + str(self.f_weights[ind]) + "\n")
        f.write("system:\n")
        f.write(" - anchor: 0\n")
        f.write("   repeat: 1\n")
        f.write("   offset: 0\n")
        f.write("   sites:\n")
        for i in range(len(self.domain_starts)):
            f.write("     - [" + self.name + "_" + str(i + 1) + ", " + str(self.domain_starts[i]) + "]\n")

    def read_domain_file(self, fname):
        f = open(fname, 'r')
        data = f.readlines()
        return [list(range(int(line.split()[0]), int(line.split()[1]) + 1))  for line in data]

    def get_angles(self):

        if self.angles is None:
            return []
        else:
            return self.angles

    def get_x_weight(self, name):
        real_name = ""
        for thing in list(name):
            if str(thing).isalpha():
                real_name += thing
        letters = ["C", "H", "O", "N", "S"]
        weights = [12.0, 1.0, 16.0, 14.0, 32.0]
        return weights[letters.index(str(real_name[0]))]

    def get_bonds(self, bond_file):
        #i j # t k r0
        f = open(bond_file)

        data = f.readlines()[1:]
        n = len(data)
        bond_matrix = np.zeros((n,2))
        bond_types = np.zeros((n,3))
        for i in range(n):
            line = data[i].split()
            bond_matrix[i][0] = int(line[0])
            bond_matrix[i][1] = int(line[1])
            bond_types[i][0] = int(i + 1)
            bond_types[i][1] = float(line[2])
            bond_types[i][2] = float(line[3])
        self.bonds = bond_matrix
        self.bond_types = bond_types

    def shift(self, vector):

        self.positions = np.add(self.positions, vector)

    def get_positions_from_lammpstrj(self, tfile, frame=0):

        from mscg import Trajectory

        cgt = Trajectory(tfile, fmt="lammpstrj")
        cgt.read_frame()
        self.positions = cgt.x

    def get_positions_from_pdb(self, filename):

        f = open(filename)
        data = f.readlines()
        count = 0
        for ind, line in enumerate(data):
            s = line.split()
            if s[0] == "ATOM":
                if s[6].count(".") == 2:
                    s[8] = s[7]
                    s[7] = s[6][7:]
                    s[6] = s[6][:7]
                if s[6].count(".") == 3:
                    s[8] = s[6][14:]
                    s[7] = s[6][7:14]
                    s[6] = s[6][:7]
                self.positions[count][0] = float(s[6])
                if s[7].count(".") == 2:
                    s[8] = s[7][7:]
                    s[7] = s[7][:7]
                self.positions[count][1] = float(s[7])
                self.positions[count][2] = float(s[8])
                count += 1

    def dump_pdb(self, dump_file=None, site_name="CA"):

        filename = dump_file
        if dump_file is None:
            filename = self.name + '_cg.pdb'
        f = open(filename, 'w')
        atom_number = 1
        string = ""
        for ind in range(len(self.site_indexes)):
            string = "ATOM"
            for i in range(6 - len(str(atom_number))):
                string += " "
            string += str(atom_number)
            atom_number += 1
            name = site_name
            for i in range(4 - len(str(name))):
                string += " "
            string += name
            r_name = "ALA"
            for i in range(6 - len(str(r_name))):
                string += " "
            string += r_name
            string += " "
            string += chr(ord("A"))
            #ind2 = 1
            for i in range(4 - len(str(ind+1))):
                string += " "
            string += str(ind+1)
            string += "      "
            pos = str(str(self.positions[ind][0])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "
            pos = str(str(self.positions[ind][1])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "
            pos = str(str(self.positions[ind][2])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "

            string += "1.00 10.00"
            for i in range(11):
                string += " "
            string += name[:1]
            string += "  \n"
            f.write(string)
        f.close()

    def dump_moldy(self, dump_file=None):

        filename = dump_file
        if dump_file is None:
            filename = self.name + '_cg.mol'
        f = open(filename, 'w')

        f.write(str(len(self.site_indexes)) + "\n")
        f.write("1 1 1\n")
        f.write("1000 0 0\n0 1000 0\n0 0 1000\n")
        string = ""
        for ind in range(len(self.site_indexes)):
            print(ind)
            pos = str(str(self.positions[ind][0])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "
            pos = str(str(self.positions[ind][1])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "
            pos = str(str(self.positions[ind][2])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + " 1  1"
            string += "  \n"
        f.write(string)
        f.close()

    def read_positions(self, filename):

        f = open(filename, "r")
        positions = []
        data = f.readlines()
        for line in data:
            if len(line) > 0:
                s = line.split()
                positions.append([float(s[0]), float(s[1]), float(s[2])])
        return positions

    def set_atom_types(self):

        self.atom_types = []
        for i in range(len(self.site_indexes)):
            self.atom_types.append([int(i+1), np.sum(self.x_weights[i])])

    def abs_path(self, string):
        return "/home/cwaltmann/PycharmProjects/Builders/" + string

    def n_atom_types(self):

        if len(self.atom_types) == 1:
            return 1
        return len(np.unique(self.atom_types, axis=0))

    def n_bond_types(self):

        if len(self.bond_types) == 1:
            return 1
        return len(np.unique(self.bond_types, axis=0))

    def n_angle_types(self):

        if self.angles is None:
            return 0
        elif len(self.angle_types) == 1:
            return 1
        else:
            return len(np.unique(self.angle_types, axis=0))

    def get_unique_atoms(self):

        if len(self.atom_types) == 1:
            return self.atom_types
        return np.unique(self.atom_types, axis=0)

    def get_unique_bonds(self):

        if len(self.bond_types) < 2:
            return self.bond_types
        return np.unique(self.bond_types, axis=0)

    def get_unique_angles(self):

        if self.angles is None:
            return []
        if len(self.angle_types) == 1:
            return self.angle_types
        return np.unique(self.angle_types, axis=0)

    def align(self, vec):
        q = QuaternionBetween(self.chain_vector(), vec)
        for x in range(len(self.positions)):
            self.positions[x] = q.orient(self.positions[x])

    def align_to_q(self, q):
        for x in range(len(self.positions)):
            self.positions[x] = q.orient(self.positions[x])

    def chain_vector(self):

        return np.subtract(self.positions[-1], self.positions[0])
