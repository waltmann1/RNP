from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class PolygonMolecule(CGMolyAbs):

    def __init__(self, N, scale=1, weight=36, unique_types=False):
        super(PolygonMolecule, self).__init__()
        if not N > 0 and type(N)  == type(int(1)):
            raise ValueError("N must be a positive integer")
        self.site_indexes = [[1]]
        rest = [[2] for _ in range(N)]
        self.site_indexes.extend(rest)
        if unique_types:
            self.site_indexes = [[i] for i in range(1,N+2)]
        self.f_weights = [[1] for _ in range(N+1)]
        self.x_weights = [[weight] for _ in range(N+1)]
        self.name = "Polygon" + str(N)
        mono_length = scale
        k=10
        k_angle = 10

        self.atom_types = [[1, weight]]
        self.atom_types.extend([[2,weight] for _ in range(N)])
        if unique_types:
            self.atom_types =  [[i, weight] for i in range(N+2)]

        dAngle = 2 * np.pi / N
        self.positions = [[0, 0, 0], [np.cos(0), np.sin(0), 0]]
        self.positions.extend([[np.cos(dAngle * i), np.sin(dAngle*i), 0] for i in range(1, N)])

        self.bonds = [[1, i] for i in range(2, N+2)]
        self.bond_types = [[1, mono_length, k] for _ in range(2, N+2)]
        if N < 1:
            bond_length_2 = np.linalg.norm(np.subtract(self.positions[1], self.positions[2]))
            self.bonds.extend([[i, i+1] for i in range(1, N)])
            self.bond_types.extend([[2, bond_length_2, k] for _ in range(1, N)])
            self.angles = [[i, 0, i + 1] for i in range(1, N)]
            self.angle_types = [[1, dAngle, k_angle] for _ in range(1, N)]
            if N < 2:
                self.bonds.extend([[1, N]])
                self.bond_types.extend([[2, bond_length_2, k]])
                self.angles.extend([[1, 0, N]])
                self.angle_types.extend([[1, dAngle, k_angle]])
        #self.bond_length = mono_length