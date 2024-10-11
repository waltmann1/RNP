from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class CAGrimer(CGMolyAbs):

    def __init__(self):
        super(CAGrimer, self).__init__()
        mass = 10
        self.site_indexes = [[] for _ in range(268)]
        self.f_weights = [[1] for _ in range(268)]
        self.x_weights = [[10] for _ in range(268)]
        self.name = "CAGrimer"
        self.positions = self.read_positions(self.abs_path("CAGrimer_pos1.txt"))
        self.atom_types = [[i+1, mass] for i in range(134)]
        self.atom_types.extend(self.atom_types)
        self.positions.extend(self.read_positions(self.abs_path("CAGrimer_pos2.txt")))
        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        #self.positions = np.add(self.positions, [50, 50, 50])
        self.get_bonds(self.abs_path("CAGrimer_bonds.txt"))
        more_bonds = []
        for bond in self.bonds:
            new = cp.deepcopy(bond)
            new[0] += 134
            new[1] += 134
            more_bonds.append(new)
        np.array(list(self.bonds).extend(more_bonds))
        self.bonds = np.concatenate((self.bonds, more_bonds), axis=0)
        self.bond_types = np.concatenate((self.bond_types, self.bond_types), axis=0)
        temp = cp.deepcopy(self.bonds)
        temp_types = cp.deepcopy(self.bond_types)

        self.get_bonds(self.abs_path("CAGrimer_bonds_inter.txt"))
        #print(temp)
        for t in self.bond_types:
            t[0] += len(np.unique(temp_types, axis=0))
        self.bonds = np.concatenate((temp, self.bonds), axis=0)
        self.bond_types = np.concatenate((temp_types, self.bond_types), axis=0)

        self.ucg_atom_types = [[i, 10] for i in range(135, 139)]

        #np.array(list(self.bonds).extend(temp))
        #np.array(list(self.bond_types).extend(temp_types))
        #print(self.bond_types)
        #print(len(self.bond_types))
        #quit()

    def n_atom_types(self):

        return len(np.unique(self.atom_types, axis=0)) + len(np.unique(self.ucg_atom_types, axis=0))

    def get_unique_atoms(self):

        a = np.unique(self.atom_types, axis=0)
        b = np.unique(self.ucg_atom_types, axis=0)

        return np.concatenate((a, b), axis=0)


