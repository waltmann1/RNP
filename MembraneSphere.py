from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp
from MDAnalysis.analysis import distances


class MembraneSphere(CGMolyAbs):

    def __init__(self, radius, particle_radius):
        super(MembraneSphere, self).__init__()
        sa = np.pi * 4 * radius * radius
        area = np.pi * particle_radius**2
        length = int(sa / area)
        self.positions = np.multiply(self.unit_sphere(length), radius)
        self.positions = np.add(self.positions, radius)
        self.site_indexes = [[] for _ in range(length)]
        self.f_weights = [[1] for _ in range(length)]
        self.x_weights = [[22.62] for _ in range(length)]
        self.site_indexes = [[] for _ in range(length)]
        self.f_weights = [[1] for _ in range(length)]
        self.x_weights = [[22.62] for _ in range(length)]
        self.name = "memb_" + str(radius) + "nm"
        self.atom_types = [[int(1), 22.62] for _ in range(length)]
        self.bonds = np.array([])
        self.bond_types = np.array([])

    def unit_sphere(self, n):

        points = []
        offset = 2 / n
        increment = np.pi * (3 - np.sqrt(5))

        for i  in range(n):
            y = ((i * offset) - 1) + offset/2
            r = np.sqrt(1 - pow(y,2))
            phi = i * increment
            x = np.cos(phi) * r
            z = np.sin(phi) * r
            points.append([x,y,z])
        return points

    def get_bond_info(self, points):

        points = np.array(points)
        neighbors = []
        real_lengths = []
        matrix = distances.distance_array(points, points)
        for i in range(len(points)):
            my_neighbors = []
            my_distances = []
            temp = cp.deepcopy(matrix[i])
            #temp.sort(reverse=True)
            temp = sorted(temp)
            temp = temp[1:]
            for j in range(6):
                my_distances.append(temp[j])
                my_neighbors.append(list(matrix[i]).index(temp[j]))
            pair = []
            for ind, n in enumerate(my_neighbors):
                if n > i:
                    pair = [i+ 1, n+1]
                else:
                    pair = [n+1, i+1]
                if pair not in neighbors:
                    neighbors.append(pair)
                    real_lengths.append(my_distances[ind])
        return neighbors, real_lengths
