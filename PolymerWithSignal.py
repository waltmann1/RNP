import numpy.linalg as la
import copy as cp
import math
from BasicPolymer import BasicPolymer


class PolymerWithSignal(BasicPolymer):

    def __init__(self, nbeads, bond_length, mass=10, k=10, separation=None, lattice_points=False):

        super(PolymerWithSignal, self).__init__(nbeads, bond_length, mass=mass, k=k, separation=separation,
                                                lattice_points=lattice_points)

        signal_fraction = .04
        signal_number = 240
        #for i in range(int(nbeads * (signal_fraction))):
        for i in range(signal_number):
            self.atom_types[i] = [2, mass]