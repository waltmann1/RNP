from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class OneBeadMoly(CGMolyAbs):

    def __init__(self, weight=36):
        super(OneBeadMoly, self).__init__()
        self.site_indexes = [[] for _ in range(1)]
        self.f_weights = [[1] for _ in range(1)]
        self.x_weights = [[weight] for _ in range(1)]
        self.name = "OneBead"
        self.positions = [[0,0,0]]
        self.atom_types = [[int(1), weight] for _ in range(1)]
        self.bonds = np.array([])
        self.bond_types = np.array([])
