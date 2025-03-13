from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class IP10(CGMolyAbs):

    def __init__(self, weight=36, unique_types=False):
        super(IP10, self).__init__()
        self.site_indexes = [[2], [1], [2]]
        if unique_types:
            self.site_indexes = [[1],[2],[3]]
        self.f_weights = [[1] for _ in range(1)]
        self.x_weights = [[weight] for _ in range(3)]
        self.name = "IP10"
        mono_length = 4.5
        k=10
        self.positions = [[-mono_length, 0, 0], [0, 0, 0], [mono_length, 0, 0]]
        self.atom_types = [[2, weight], [1, weight], [2, weight]]
        if unique_types:
            self.atom_types =  [[1, weight], [2, weight], [3, weight]]
        self.bonds = [[i, i + 1] for i in range(1, 3)]
        self.bond_types = [[1, mono_length, k] for _ in range(3 - 1)]
        self.bond_length = mono_length

        angles_A = [[i, i + 1, i + 2] for i in range(1, 2)]
        #angles_B = np.add(angles_A, mono_length)
        self.angles = angles_A
        #self.angles.extend(angles_B)
        theta = 180
        k = 200
        # k=0
        self.angle_types = [[1, theta, k] for _ in range(len(self.angles))]