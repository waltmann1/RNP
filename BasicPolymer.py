from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import copy as cp
import math
class BasicPolymer(CGMolyAbs):

    def __init__(self, nbeads, bond_length, mass=10, k=10, separation=None, lattice_points=False):
        super(BasicPolymer, self).__init__()

        self.name = "poly" + str(nbeads)
        self.site_indexes = [[] for _ in range(nbeads)]
        self.f_weights = [[1] for _ in range(nbeads)]
        self.x_weights = [[mass] for _ in range(nbeads)]
        sep = bond_length
        if separation is not None:
            sep = separation
        if lattice_points:
            self.positions = self.lattice_points(nbeads, bond_length)
        else:
            self.positions = self.spiral_points(nbeads, arc=bond_length, separation=sep)
        self.atom_types = [[1, mass] for _ in range(nbeads)]
        self.bonds = [[i, i+1] for i in range(1, nbeads)]
        self.bond_types = [[1, bond_length, k] for _ in range(nbeads - 1)]
        self.bond_length = bond_length
        self.k = k
    def spiral_points(self,n, arc=.5, separation=4):
        """generate points on an Archimedes' spiral
        with `arc` giving the length of arc between two points
        and `separation` giving the distance between consecutive
        turnings
        - approximate arc length with circle arc at given distance
        - use a spiral equation r = b * phi
        """

        def p2c(r, phi):
            """polar to cartesian
            """
            return [r * math.cos(phi), r * math.sin(phi), 0]

        # yield a point at origin
        points=  [[0,0,0]]

        # initialize the next point in the required distance
        r = arc
        b = separation / (2 * math.pi)
        # find the first phi to satisfy distance of `arc` to the second point
        phi = float(r) / b
        count = 0
        while count < n - 1:
            points.append(p2c(r, phi))
            # advance the variables
            # calculate phi that will give desired arc length at current radius
            # (approximating with circle)
            phi += float(arc) / r
            r = b * phi
            count += 1
        return points

    def lattice_points(self, n, bond_length):

        points = np.array([[0,0,0]])
        options = [0,1,2]
        options_2 = [-1, 1]
        for i in range(n-1):
            coord = np.random.choice(options)
            dir = np.random.choice(options_2)
            attempt = cp.deepcopy(points[-1])
            attempt[coord] += dir
            count = 0
            while np.any(np.all(attempt == points, axis=1)):
                count += 1
                coord = np.random.choice(options)
                dir = np.random.choice(options_2)
                attempt = cp.deepcopy(points[-1])
                attempt[coord] += dir
                if count > 10:
                    attempt[coord] += dir

            #attempt = np.array([points])
            count = 0
            print(coord,dir, attempt, points[-1], len(points))
            points = np.concatenate((points, [attempt]))

        return np.multiply(points, bond_length)
