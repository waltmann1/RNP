from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class SingleCA(CGMolyAbs):

    def __init__(self):
        super(SingleCA, self).__init__()
        mass = 10
        self.site_indexes = [[] for _ in range(134)]
        self.f_weights = [[1] for _ in range(134)]
        self.x_weights = [[10] for _ in range(134)]
        self.name = "SingleCA"
        self.positions = self.read_positions(self.abs_path("CAGrimer_pos1.txt"))
        self.atom_types = [[i+1, mass] for i in range(134)]

        ave_pos = np.mean(self.positions, axis=0)
        self.positions = np.subtract(self.positions, ave_pos)
        #self.positions = np.add(self.positions, [50, 50, 50])
        self.get_bonds(self.abs_path("CAGrimer_bonds.txt"))