from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class Crowder(CGMolyAbs):

    def __init__(self):
        super(Crowder, self).__init__()
        self.site_indexes = [[] for _ in range(42)]
        self.f_weights = [[1] for _ in range(42)]
        self.x_weights = [[22.62] for _ in range(42)]
        self.name = "Crowder"
        self.positions = self.read_positions(self.abs_path("Crowder_pos.txt"))
        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        #self.positions = np.add(self.positions, [50, 50, 50])
        self.atom_types = [[int(1), 22.62] for _ in range(42)]
        self.get_bonds(self.abs_path("Crowder_bonds.txt"))
        for t in self.bond_types:
            if t[1] == self.bond_types[0][1]:
                t[0] = 2
            else:
                t[0] = 1
