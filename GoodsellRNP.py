from __future__ import division
from CGMolyComplex import CGMolyComplex
from CGMolyAbs import CGMolyAbs
from BasicPolymer import BasicPolymer
from MDAnalysis.analysis import distances
import numpy as np
import numpy.linalg as la
import Bio.PDB
import copy as cp


class GoodsellRNP(CGMolyComplex):

    def __init__(self, PDBName, model_number=1, mass=None, l_integrase=False):

        molys = None
        parser = Bio.PDB.PDBParser()
        self.struct = parser.get_structure(PDBName[:-4], PDBName)
        chains = [chain.get_atoms() for chain in self.struct.get_chains()]
        positions = []
        for ind, chain in enumerate(chains[(model_number-1)*5: (model_number) * 5]):
            pos = [np.multiply(atom.get_coord(),10) for atom in chain]
            positions.append(pos)
        dimer_pos = positions[0]
        dimer_pos.extend(positions[1])
        rna = GoodsellRNADimer(dimer_pos)
        grace_pos = self.read_integrase(pdbname=PDBName, modelnumber=model_number)
        #print(grace_pos[:4])
        #quit()
        indiv_pos = []
        in_rna_crosslinks = []
        in_unit_counter = 0
        #bond_types = rna.n_bond_types()
        for i in range(int(len(grace_pos)/4)):
            indiv_pos.append(grace_pos[i*4:(i+1)*4])
            ave_pos = np.average(grace_pos[i*4:(i+1)*4], axis=0)
            in_rna_dist = distances.distance_array(np.array(grace_pos[i*4:(i+1)*4]), np.array(dimer_pos))
            ave_rna_dist = distances.distance_array(np.array(ave_pos), np.array(dimer_pos))
            ave_rna_exclude = 5000 * (ave_rna_dist[0] < 46)

            #in_rna_dist = np.abs(np.subtract(in_rna_dist, 23))
            for k in range(4):
                in_rna_dist[k] = np.add(in_rna_dist[k], ave_rna_exclude)

            rna_indices = [np.argmin(in_rna_dist[i])  for i in range(4)]
            #bond_type = bond_types + 1
            for i in range(0,4):
                #for j in range(i+1,4):
                    #link_dist = 6.5
                 #   if (i == 0 and j == 1) or (i == 2) and (j == 3):
                  #      link_dist = 92
                  #      bond_type = bond_types + 2
                  #  else:
                  #      bond_type = bond_types + 1
                  #      link_dist = 65

                # add bonds from in to rna
                in_rna_crosslinks.append([rna_indices[i], len(dimer_pos) + in_unit_counter +1])
                in_unit_counter += 1
                if l_integrase and in_unit_counter % 7 == 4:
                    in_unit_counter += 3

        integrase = GoodsellINTetramer()
        if l_integrase:
            integrase = LyumkisINTetramer()
        molys = [rna, integrase]
        super(GoodsellRNP, self).__init__(molys, numbers=[1, len(indiv_pos)])
        # set integrase positions
        self.reset_master_positions()
        for i in range(len(indiv_pos)):
            self.set_moly_positions(integrase.get_positions(indiv_pos[i]), 1, i)
        #add RNA-IN crosslinks
        self.intra_bonds = in_rna_crosslinks
        length = 23
        if l_integrase:
            length= 15
        self.intra_bond_types = [[1, length, 1] for _ in range(len(self.intra_bonds))]

    def read_integrase(self, pdbname, modelnumber):

        f = open(pdbname, "r")
        data = f.readlines()
        positions = []
        on = False
        for ind, line in enumerate(data):
            s = line.split()
            if len(s) > 3:
                if s[3] == "IN" and on:
                    #print(s)
                    pos = [float(s[6]) * 10, float(s[7]) * 10, float(s[8]) * 10]
                    #print(pos)
                    #quit()
                    positions.append(pos)
            if s[0] == "MODEL":
                if s[1] == str(modelnumber):
                    on =True
                else:
                    on = False
        return positions


class GoodsellRNADimer(BasicPolymer):

    def __init__(self, positions, diff_types=True, mass=None):
        #mass includes NC


        if mass is not None:
            mass =10
        else:
            mass = 330 * 3 + 1800

        super(GoodsellRNADimer, self).__init__(len(positions),10.0, k=5, mass=mass)
        self.bonds.remove([3058,3059])
        self.bond_types = self.bond_types[:-1]
        self.positions = positions
        mono_length = int(len(self.positions)/2)

        if diff_types:
            for i in range(112, mono_length):
                self.atom_types[i][0] = 2
                self.atom_types[i+ mono_length][0] = 2

        crosslinks = [[1, 19], [2,18], [3,17], [4,16], [6,15], [7,14], [9, 13], [10,12], [20,35],
                      [21,34], [22,32], [23,31], [24,30], [25,29], [36,112], [37,111], [38, 110],
                      [41,74], [42,73], [43, 72], [44,59], [45,58], [46,57], [48,53], [49, 52],
                      [76,109], [77, 108], [78,93], [80,92], [81,91], [82,89], [83,88]]
        c2 = np.add(crosslinks, 3058)

        crosslinks.extend(c2)
        crosslinks.append([85, 85 + 3058])
        crosslinks.append([86, 86 + 3058])
        crosslinks.append([87, 87 + 3058])

        for crosslink in crosslinks:
            self.bonds.append(crosslink)
            self.bond_types.append([2, self.bond_length, self.k])

        angles_A = [[i, i + 1, i + 2] for i in range(1, mono_length-1)]
        angles_B = np.add(angles_A, mono_length)
        self.angles = angles_A
        self.angles.extend(angles_B)
        theta = 180
        k= 2.5 * .6
        #k=0
        self.angle_types = [[1, theta, k] for _ in range(len(self.angles))]
class GoodsellINTetramer(CGMolyAbs):

    def __init__(self, mass=None):

        super(GoodsellINTetramer, self).__init__()
        #print("hi")
        positions = [[0,0,0] for i in range(4)]
        print(positions)
        #print(distances.distance_array(positions, positions))
        nbeads = len(positions)

        if mass is not None:
            mass=10
        else:
            mass = 30000

        #mass = 10
        self.name = "IN" + str(nbeads)
        self.site_indexes = [[] for _ in range(nbeads)]
        self.f_weights = [[1] for _ in range(nbeads)]
        self.x_weights = [[mass] for _ in range(nbeads)]
        self.atom_types = [[1, mass] for _ in range(nbeads)]
        self.bonds = [[1, 2], [3,4], [1,3], [1,4], [2,3], [2,4]]
        k=5
        self.bond_types = [[1, 46, k], [1, 46, k],[2, 32.5, k], [2,32.5,k], [2, 32.5, k], [2,32.5,k]]

    def get_positions(self, g_int_pos):

        return g_int_pos


class LyumkisINTetramer(CGMolyAbs):

    def __init__(self, mass=None):

        super(LyumkisINTetramer, self).__init__()

        positions = [[0, 0, 0] for i in range(7)]
        print(positions)
        # print(distances.distance_array(positions, positions))
        nbeads = len(positions)

        if mass is not None:
            mass = 10
        else:
            mass = 30000

        # mass = 10
        self.name = "IN" + str(nbeads)
        self.site_indexes = [[] for _ in range(nbeads)]
        self.f_weights = [[1] for _ in range(nbeads)]
        self.x_weights = [[mass] for _ in range(nbeads)]
        self.atom_types = [[1, mass] for _ in range(4)] + [[2, mass]] + [[3,mass], [3, mass]]
        self.bonds = [[1, 2], [3, 4],
                      [1, 3], [1, 4], [2, 3], [2, 4],
                      [1, 5], [2, 5], [3, 5], [4, 5],
                      [1, 6], [3, 6], [4, 6],
                      [2, 7], [3, 7], [4, 7],
                      [5, 6], [5, 7],
                      [6,7]]

        k = 500
        long = 70
        short = 70 / np.sqrt(2)
        bottom = 46.1
        top = 20
        top_far = np.sqrt(53) * 10
        last = np.sqrt(37.25) * 10
        self.bond_types = [[1, long, k], [1, long, k],
                           [2, short, k], [2, short, k], [2, short, k], [2, short, k],
                           [3, bottom, k], [3, bottom, k], [3, bottom, k], [3, bottom, k],
                           [4, top, k], [5, top_far, k], [5, top_far, k],
                           [4, top, k], [5, top_far, k], [5, top_far, k],
                           [6, last, k], [6, last, k],
                           [7, long, k]]


    def get_positions(self, g_int_pos):

        new_int_pos = []
        one_two_vector = np.subtract(g_int_pos[1], g_int_pos[0])
        three_four_vector= np.subtract(g_int_pos[3], g_int_pos[2])

        one_two_vector = np.divide(one_two_vector, np.linalg.norm(one_two_vector))
        three_four_vector = np.divide(three_four_vector, np.linalg.norm(three_four_vector))
        new_1 = np.add(g_int_pos[0], np.multiply(one_two_vector, -12))
        #new_1 = g_int_pos[0]
        new_2 = np.add(g_int_pos[1], np.multiply(one_two_vector, 12))
        #new_2 = g_int_pos[1]
        new_3 = np.add(g_int_pos[2], np.multiply(three_four_vector, -12))
        #new_3 = g_int_pos[2]
        new_4 = np.add(g_int_pos[3], np.multiply(three_four_vector, 12))
        #new_4 = g_int_pos[3]

        triangle = [new_1, new_2, new_3]
        v1 = np.subtract(triangle[0], triangle[1])
        v2 = np.subtract(triangle[2], triangle[1])
        n_long = np.cross(v1, v2)
        n_unit = np.divide(n_long, np.linalg.norm(n_long))
        center = np.average(g_int_pos, axis=0)

        new_5 = np.add(center, np.multiply(-30, n_unit))
        #new_5 = center
        new_6 = np.add(new_1, np.multiply(20, n_unit))
        #new_6 = new_1
        new_7 = np.add(new_2, np.multiply(20, n_unit))
        #new_7 = new_2
        return [new_1, new_2, new_3, new_4, new_5, new_6, new_7]




