from __future__ import division
import numpy as np
import random
from CGMolyAbs import CGMolyAbs
import copy as cp


class CGMolyComplex(object):

    def __init__(self, molys, name=None, numbers=None, box=[1000, 1000, 1000]):

        self.molys = molys
        self.master_positions = None
        self.numbers = [1 for _ in range(len(molys))]
        if numbers is not None:
            self.numbers = numbers

        self.box = box

        self.name = "complex"
        self.intra_bonds = None
        self.intra_bond_types = None
        if name is not None:
            self.name = name

    def to_cg_moly(self):

        new_moly = CGMolyAbs()
        new_moly.positions = self.master_positions
        new_moly.angles = []
        new_moly.bonds = []
        new_moly.angle_types = []
        new_moly.atom_types = []
        #print("bond tpyes", len(new_moly.bond_types))

        atom_counter = 0
        atom_type_counter = 0
        bond_type_counter = 0
        angle_type_counter = 0
        for ind, moly in enumerate(self.molys):
            new_atypes = cp.deepcopy(moly.atom_types)
            for type in new_atypes:
                type[0] += atom_type_counter
            if moly.bond_types is not None:
                new_botypes = cp.deepcopy(moly.bond_types)
                #print("new bo types 1", len(new_botypes))
                for type in new_botypes:
                    type[0] += bond_type_counter
                #print("new bo types 2", len(new_botypes))
            if moly.angle_types is not None:
                new_antypes = cp.deepcopy(moly.angle_types)
                for type in new_antypes:
                    type[0] += angle_type_counter
            #print("atom tpyes", len(new_moly.atom_types))
            for number in range(self.numbers[ind]):
                new_moly.f_weights.extend(moly.f_weights)
                new_moly.f_weights.extend(moly.x_weights)
                if moly.bonds is not None:
                    new_moly.bonds.extend(np.add(moly.bonds, atom_counter))
                    #print("bond tpyes 3", len(new_moly.bond_types))
                    new_moly.bond_types.extend(new_botypes)
                    #print("bond tpyes 4", len(new_moly.bond_types))
                if moly.angles is not None:
                    new_moly.angles.extend(np.add(moly.angles, atom_counter))
                    new_moly.angle_types.extend(new_antypes)
                if moly.atom_types is None:
                    moly.set_atom_types()

                new_moly.atom_types.extend(new_atypes)
                new_moly.site_indexes.extend(moly.site_indexes)

                atom_counter += len(moly.site_indexes)
            #print("bond tpyes 5", len(new_moly.bond_types))

            bond_type_counter += moly.n_bond_types()
            angle_type_counter += moly.n_angle_types()
            atom_type_counter += moly.n_atom_types()


        if self.intra_bond_types is not None:
            new_ibotypes = cp.deepcopy(self.intra_bond_types)
            for type in new_ibotypes:
                type[0] += bond_type_counter
            new_moly.bond_types.extend(new_ibotypes)
            new_moly.bonds.extend(self.intra_bonds)

        #print(self.intra_bonds)
        #print(self.intra_bond_types)
        #quit()
        return new_moly


    def increase_numbers(self, moly_index, num_to_add):


        insert_index = np.sum([len(self.molys[i].site_indexes) * self.numbers[i]  for i in range(moly_index+1)])
        before = list(self.master_positions[:insert_index])
        after = list(self.master_positions[insert_index:])

        insert_pos = []
        new_moly = cp.deepcopy(self.molys[moly_index])
        for _ in range(num_to_add):
            ave_pos = np.average(new_moly.positions, axis=0)
            new_moly.shift(-ave_pos)
            random_pos = [self.box[i] * np.random.random_sample() for i in range(3)]
            new_moly.shift(random_pos)
            insert_pos.extend(new_moly.positions)

        insert_pos = np.array(insert_pos)
        before.extend(insert_pos)
        before.extend(after)
        self.master_positions = np.array(before)
        self.numbers[moly_index] += num_to_add

    def write_yaml(self, name=None):

        master_site_indexes = self.molys[0].site_indexes
        master_x_weights = self.molys[0].x_weights
        master_f_weights = self.molys[0].f_weights
        master_site_names = []
        for i in range(len(self.molys[0].site_indexes)):
            master_site_names.append(self.molys[0].name + "_" + str(i + 1))
        master_domain_starts = self.molys[0].domain_starts
        current = self.molys[0].domain_starts[-1] + len(self.molys[0].site_indexes[-1])
        for moly in self.molys[1:]:
            master_site_indexes.extend(moly.site_indexes)
            master_f_weights.extend(moly.f_weights)
            master_x_weights.extend(moly.x_weights)
            master_domain_starts.extend(np.add(moly.domain_starts, current))
            current += moly.domain_starts[-1] + len(moly.site_indexes[-1])
            site_names = []
            for i in range(len(moly.site_indexes)):
                site_names.append(moly.name + "_" + str(i+1))
            master_site_names.extend(site_names)

        if name is None:
            name = self.name
        f = open(name + ".yaml", 'w')

        f.write("site-types:\n")
        for ind in range(len(master_site_indexes)):
            f.write("   " + master_site_names[ind] + ":\n")
            f.write("        " + "index:    " + str(master_site_indexes[ind]) + "\n")
            f.write("        " + "x-weight: " + str(master_x_weights[ind]) + "\n")
            f.write("        " + "f-weight: " + str(master_f_weights[ind]) + "\n")
        f.write("system:\n")
        f.write(" - anchor: 0\n")
        f.write("   repeat: 1\n")
        f.write("   offset: 0\n")
        f.write("   sites:\n")
        for i in range(len(master_site_names)):
            f.write("     - [" + master_site_names[i] + ", " + str(master_domain_starts[i]) + "]\n")

    def get_master_positions_from_lammpstrj(self, tfile, frame=0, keep_box=True):

        from mscg import Trajectory

        cgt = Trajectory(tfile, fmt="lammpstrj")
        for _ in range(frame+1):
            cgt.read_frame()
        self.master_positions = cgt.x
        if keep_box:
            self.box = cgt.box

    def get_master_positions_from_pdb(self, filename):

        f = open(filename)
        self.reset_master_positions()
        data = f.readlines()
        count = 0
        for ind, line in enumerate(data):
            s = line.split()
            if s[0] == "ATOM":
                if s[6].count(".") == 2:
                    s[8] = s[7]
                    s[7] = s[6][7:]
                    s[6] = s[6][:7]
                if s[6].count(".") == 3:
                    s[8] = s[6][14:]
                    s[7] = s[6][7:14]
                    s[6] = s[6][:7]
                self.master_positions[count][0] = float(s[6])
                if s[7].count(".") == 2:
                    s[8] = s[7][7:]
                    s[7] = s[7][:7]
                self.master_positions[count][1] = float(s[7])
                self.master_positions[count][2] = float(s[8])
                count += 1
    def reset_master_positions(self):

        natoms = np.sum([len(moly.site_indexes) * self.numbers[ind] for ind, moly in enumerate(self.molys)])
        self.master_positions= np.zeros((natoms, 3))
        counter = 0
        for ind, moly in enumerate(self.molys):
            for _ in range(self.numbers[ind]):
                for pos in moly.positions:
                    self.master_positions[counter] = pos
                    counter += 1

    def write_data_file(self, name=None):

        if name is None:
            name = self.name + ".data"

        for moly in self.molys:
            if moly.atom_types is None:
                moly.set_atom_types()

        f = open(name, "w")
        f.write("LAMMPS data file via mscg, version 20 Sep 2021, timestep = 0\n\n")
        natomtypes = np.sum([moly.n_atom_types() for moly in self.molys])
        natoms = np.sum([len(moly.site_indexes) * self.numbers[ind] for ind, moly in enumerate(self.molys)])
        f.write(str(natoms) + " atoms\n")
        nbondtypes = np.sum([moly.n_bond_types() for moly in self.molys]) + len(self.get_unique_intra_bonds())
        nangletypes = np.sum([moly.n_angle_types() for moly in self.molys])
        nbonds = np.sum([len(moly.bonds) * self.numbers[ind] for ind, moly in enumerate(self.molys)])
        if self.intra_bonds is not None:
            nbonds += len(self.intra_bonds)

        nangles = np.sum([len(moly.get_angles()) * self.numbers[ind] for ind, moly in enumerate(self.molys)])
        f.write(str(int(nbonds)) + " bonds\n")
        f.write(str(int(nangles)) + " angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n")
        f.write("\n")

        f.write(str(natomtypes) + " atom types\n")
        f.write(str(nbondtypes) + " bond types\n")
        if nangles != 0:
            f.write(str(nangletypes) + " angle types\n")
        f.write("\n")

        f.write("0 " + str(self.box[0]) +" xlo xhi\n")
        f.write("0 " + str(self.box[1]) + " ylo yhi\n")
        f.write("0 " + str(self.box[2]) + " zlo zhi\n")
        f.write("\n")

        f.write("Masses\n\n")
        f.write(self.masses_string())

        f.write("Bond Coeffs #harmonic K R0\n\n")
        counter = 0
        for moly in self.molys:
            uni_bond_types = moly.get_unique_bonds()
            #print("uni bond types", moly.name)
            #print(uni_bond_types)
            for i in range(len(uni_bond_types)):
                f.write(str(int(counter + uni_bond_types[i][0])) + " " + str(uni_bond_types[i][2]) + " " + str(uni_bond_types[i][1]) + "\n")
            counter += len(uni_bond_types)
        if self.intra_bond_types is not None:
            for uni_bond in self.get_unique_intra_bonds():
                f.write(str(int(counter + 1)) + " " + str(uni_bond[2]) + " " + str(uni_bond[1]) + "\n")
                counter += 1
        f.write("\n")

        if nangles > 0:
            f.write("Angle Coeffs #harmonic K theta_0\n\n")
            counter = 0
            for moly in self.molys:
                uni_angle_types = moly.get_unique_angles()
                for i in range(len(uni_angle_types)):
                    f.write(str(int(counter + uni_angle_types[i][0])) + " " + str(uni_angle_types[i][2]) + " " + str(
                    uni_angle_types[i][1]) + "\n")
                counter += len(uni_angle_types)
        f.write("\n")

        f.write("Atoms #full\n\n")
        f.write(self.atoms_string())

        f.write("Bonds\n\n")
        f.write(self.bonds_string())

        if nangles != 0:
            f.write("Angles\n\n")
            f.write(self.angles_string())


    #cut in angstroms
    def calc_intracomplex_bonds(self, cut=30):

        intra_bonds = []
        offsets = []
        counter = 0
        for moly in self.molys:
            offsets.append(counter)
            counter += len(moly.positions)
        for i, molyi in enumerate(self.molys):
            for j in range(i + 1, len(self.molys)):
                distances = self.look_for_connections(molyi.positions, self.molys[j].positions, cut=cut)
                for x in range(len(molyi.positions)):
                    for y in range(len(self.molys[j].positions)):
                        d = distances[x][y]
                        if d < cut:
                            ind1 = offsets[i] + x + 1
                            ind2 = offsets[j] + y + 1
                            intra_bonds.append([ind1, ind2, d])

        self.intra_bonds = intra_bonds

    def look_for_connections(self, positionsi, positionsj, cut=30):

        from MDAnalysis.analysis import distances

        dist_array = distances.distance_array(positionsi, positionsj)
        return dist_array

    def masses_string(self):

        sring = ""
        counter = 0
        for ind, moly in enumerate(self.molys):
            if moly.atom_types is None:
                moly.set_atom_types()
            uni_atoms = moly.get_unique_atoms()
            for atom in uni_atoms:
                sring += str(counter + 1) + " " + str(atom[1]) + "\n"
                counter += 1
        sring += "\n"
        return sring

    def atoms_string(self):

        sring = ""
        counter = 0
        atoms_offset = 0
        mol_tag = 1
        for ind, moly in enumerate(self.molys):
            for x in range(self.numbers[ind]):
                for i in range(len(moly.site_indexes)):
                    #print(ind, moly, x, i, counter, mol_tag)
                    #print(moly.atom_types[i][0], self.master_positions[counter])
                    sring += (str(counter + 1) + " " + str(mol_tag) + " " + str(atoms_offset + moly.atom_types[i][0]) + " 0 " + str(self.master_positions[counter][0]) + " " +
                              str(self.master_positions[counter][1]) + " " + str(self.master_positions[counter][2]) + "\n")
                    counter += 1
                mol_tag += 1
            atoms_offset += moly.n_atom_types()
        sring += "\n"
        return sring

    def bonds_string(self):

        sring = ""
        counter = 0
        offset = 0
        bond_offset = 0
        for ind, moly in enumerate(self.molys):
            for x in range(self.numbers[ind]):
                for i in range(len(moly.bonds)):
                    sring += str(counter + 1) + " " + str(int(bond_offset + moly.bond_types[i][0])) + " " + str(int(offset + moly.bonds[i][0])) + " " + str(
                        int(offset + moly.bonds[i][1])) + "\n"
                    counter += 1
                offset += len(moly.site_indexes)
            bond_offset += moly.n_bond_types()
        if self.intra_bonds is not None:
            for i, bond in enumerate(self.intra_bonds):
                sring += str(counter + 1) + " " + str(bond_offset + self.intra_bond_types[i][0]) + " " + str(int(bond[0])) + " " + str(
                    int(bond[1])) + "\n"
                counter += 1
        sring += "\n"
        return sring

    def angles_string(self):

        sring = ""
        counter = 0
        offset = 0
        angle_offset = 0
        for ind, moly in enumerate(self.molys):
            for x in range(self.numbers[ind]):
                for i in range(len(moly.get_angles())):
                    sring += (str(counter + 1) + " " + str(int(angle_offset + moly.angle_types[i][0])) + " " + str(int(offset + moly.angles[i][0])) + " " +
                              str(int(offset + moly.angles[i][1])) + " " + str(int(offset + moly.angles[i][2])) + "\n")
                    counter += 1
                offset += len(moly.site_indexes)
            angle_offset += moly.n_angle_types()
        sring += "\n"
        return sring

    def get_start_end(self, molyindex1, molyindex2):

        start = np.sum([len(self.molys[mi1].site_indexes) * self.numbers[mi1] for mi1 in range(molyindex1)])
        start += len(self.molys[molyindex1].site_indexes) * (molyindex2)
        end = start + len(self.molys[molyindex1].site_indexes)
        return  int(start), int(end)

    def shift_moly(self, vector, molyindex1, molyindex2):
        start, end = self.get_start_end(molyindex1, molyindex2)
        self.master_positions[start:end] = np.add(self.master_positions[start:end], vector)


    def extend_master_positions(self, molyindex1):

        moly = self.molys[molyindex1]
        #print(moly.name)
        moly_posititons = np.array(moly.positions)
        #print(moly.positions.shape)
        #print("master", self.master_positions.shape)
        #print(self.master_positions)
        self.master_positions = np.append(np.array(self.master_positions), moly_posititons, axis=0)

    def set_moly_positions(self, positions, molyindex1, molyindex2):
        start, end = self.get_start_end(molyindex1, molyindex2)
        self.master_positions[start:end] = positions

    def center_complex(self):

        center = np.divide(self.box,2)
        pos = np.average(self.master_positions, axis=0)
        self.master_positions = np.add(self.master_positions, np.subtract(center, pos))

    def center_moly(self, molyindex1, molyindex2):

        start, end = self.get_start_end(molyindex1, molyindex2)
        #center = np.mean(self.master_positions[start:end])
        center = np.divide(self.box,2)
        self.master_positions[start:end] = np.add(self.master_positions[start:end], center)

    def randomize_positions(self):

        for ind1 in range(len(self.molys)):
            for ind2 in range(self.numbers[ind1]):
                vector = np.multiply(self.box, .05 + .9 * np.random.rand(1,3)[0])
                self.center_moly(ind1, ind2)
                self.shift_moly(vector, ind1, ind2)

    def pad_box(self, padding):

        self.master_positions = np.add(self.master_positions, padding)

        self.box = np.add(self.box, np.multiply(padding, 2))

    def shift_master(self, vector):

        self.master_positions = np.add(self.master_positions, vector)

    def reset_box(self, box):

        self.box=box

    def make_random_lattice(self, spacing=100, insert_polymer=False, size=0, center_poly=True):

        self.reset_master_positions()
        for ind, moly in enumerate(self.molys):
            average_pos = np.average(moly.positions, axis=0)
            for j in range(self.numbers[ind]):
                self.shift_moly(-average_pos,ind,j)
                print(average_pos, ind, j)


        extras = 0
        blocked_length = 0
        center_factor = 0
        if insert_polymer:
            blocked_length = np.ceil(size/spacing)
            extras = np.power(blocked_length, 3)


        rows = int(np.ceil(np.power(np.sum(self.numbers) + extras, 1 / 3)))
        if insert_polymer and center_poly:
            center_factor = int( (rows-blocked_length)/2 )
        print(np.sum(self.numbers), rows, blocked_length)

        lattice_points = [[i,j,k] for i in range(rows) for j in range(rows) for k in range(rows)
                          if i >= blocked_length or j >= blocked_length or k >= blocked_length]
        if np.sum(self.numbers) > len(lattice_points):
            rows += 1
            lattice_points = [[i, j, k] for i in range(rows) for j in range(rows) for k in range(rows)
                              if i >= blocked_length or j >= blocked_length or k >= blocked_length]

        if center_poly:
            lattice_points = [[(point + center_factor) % rows  for point in vector] for vector in lattice_points]




        box_l = spacing * rows
        self.box= [box_l, box_l, box_l]
        #final_lattice
        print(len(lattice_points), np.sum(self.numbers))
        keep = random.sample(range(len(lattice_points)), np.sum(self.numbers))
        offset = [spacing/2, spacing/2, spacing/2]
        final_lattice = [np.add(np.multiply(lattice_points[k], spacing), offset) for k in keep]
        if insert_polymer:
            final_lattice[-1] = np.add(offset, np.multiply(.5 * blocked_length, [spacing, spacing, spacing]))
            if center_poly:
                final_lattice[-1] = np.add(final_lattice[-1], np.multiply(center_factor, [spacing, spacing, spacing]))
            print(final_lattice[-1])

        counter = 0
        for ind, moly in enumerate(self.molys):
            for n in range(self.numbers[ind]):
                #if n == self.numbers[ind] -1 and ind == len(self.molys) - 1:
                #    self.center_moly(ind, n)
                self.shift_moly(final_lattice[counter], ind, n)
                counter += 1

    def get_unique_intra_bonds(self):

        if self.intra_bond_types is None:
            return []
        if len(self.intra_bond_types) == 1:
            return self.intra_bond_types
        return np.unique(self.intra_bond_types, axis=0)

    def dump_pdb(self, dump_file=None, site_name="CA"):

        filename = dump_file
        if dump_file is None:
            filename = self.name + '_cg.pdb'
        f = open(filename, 'w')
        atom_number = 1
        string = ""
        for ind in range(len(self.master_positions)):
            string = "ATOM"
            for i in range(6 - len(str(atom_number))):
                string += " "
            string += str(atom_number)
            atom_number += 1
            name = site_name
            for i in range(4 - len(str(name))):
                string += " "
            string += name
            r_name = "ALA"
            for i in range(6 - len(str(r_name))):
                string += " "
            string += r_name
            string += " "
            string += chr(ord("A"))
            #ind2 = 1
            for i in range(4 - len(str(ind+1))):
                string += " "
            string += str(ind+1)
            string += "      "
            pos = str(str(self.master_positions[ind][0])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "
            pos = str(str(self.master_positions[ind][1])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "
            pos = str(str(self.master_positions[ind][2])[:6])
            for i in range(6 - len(pos)):
                pos += " "
            string += pos + "  "

            string += "1.00 10.00"
            for i in range(11):
                string += " "
            string += name[:1]
            string += "  \n"
            f.write(string)
        f.close()














