import numpy as np
import copy as cp
import math
from Template import Template

class TemplateWithPatch(Template):

    def __init__(self, radius, bond_k=100, phi_limit=5):

        super(TemplateWithPatch, self).__init__(radius, bond_k=bond_k)

        for i in range(len(self.positions)):
            pos = self.positions[i]
            phi = np.rad2deg(np.arccos(pos[2]/np.linalg.norm(pos)))
            if phi < phi_limit:
                self.atom_types[i] = [2, self.atom_types[i][1]]