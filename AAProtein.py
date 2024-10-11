import math as m
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import align
from MDAnalysis.analysis import rms
import MDAnalysis.transformations as trans
import copy as cp
import yaml


class AAProtein(object):

    def __init__(self, pdbfile, yaml=None, trajectory=None):

        if trajectory is None:
            self.u = mda.Universe(pdbfile)
        else:
            self.u = mda.Universe(pdbfile, trajectory)
        self.backbone = self.u.select_atoms("name CA and resid 0-221")
        self.backbone = self.backbone[:221]
        self.all = self.u.select_atoms("all")
        self.cg_map = None
        if not yaml:
            self.all = self.all[:int(len(self.all)/2)]
        else:
            self.cg_map = self.read_cg_yaml(yaml)
        dim = [10000,10000,10000, 90, 90, 90]
        transform = mda.transformations.boxdimensions.set_dimensions(dim)
        self.u.trajectory.add_transformations(transform)

        #print(self.backbone)
        #print(len(self.backbone.positions))


    def get_cg_positions(self):

        if self.cg_map is None:
            positions = self.backbone.positions[self.read_cg_map()]
        else:
            shifted_map = [np.subtract(cg_map_i, 0) for cg_map_i in self.cg_map]
            positions = [np.mean(self.all.positions[site], axis=0) for site in shifted_map]

        return positions

    def get_non_cg_positions(self):

        if self.cg_map is None:
            cg_positions = self.read_cg_map()
            non_cg_positions = [i for i in range(221) if i not in cg_positions]
            positions = self.backbone.positions[non_cg_positions]
            return positions

    def set_frame(self, frame_number):

        frame = self.u.trajectory[frame_number]
        return frame

    def read_cg_yaml(self, yname):

        with open(yname, 'r') as stream:
            data_loaded = yaml.safe_load(stream)

        indexes = []

        for pair in data_loaded['system'][0]['sites']:
            indexes.append(np.add(data_loaded['site-types'][pair[0]]['index'], pair[1]))

        return indexes


    def read_cg_map(self, name="definition.dat"):

        if self.cg_map:
            return self.cg_map

        f = open(name)

        data = f.readlines()
        resid = []
        for line in data[:130]:
            s = line.split()
            id = int(s[2][2:])
            resid.append(id)

        for line in data[131:135]:
            s = line.split()
            id = int(s[2][3:])
            resid.append(id)

        return resid

    def create_cg_group(self, cg_positions):
        if self.cg_map is None:
            shifted_map = np.subtract(np.array(self.read_cg_map()), 1)
            cg_group = self.backbone[shifted_map]
            cg_group = mda.Merge(cg_group).select_atoms("name CA")
        else:
            #indexes = [self.cg_map[i][int(len(self.cg_map[i])//2)] for i in range(len(self.cg_map))]
            indexes = [self.cg_map[i][0] for i in range(len(self.cg_map))]
            cg_group = self.all[indexes]
            cg_group = mda.Merge(cg_group).select_atoms("all")

        for i in range(len(cg_positions)):
            cg_group[i].position = cg_positions[i]
            #cg_group[i].segment.segid = "B"
            #print(cg_group[i].segment)
        #quit()
        return cg_group


    def write_cg_group(self,cg_positions, name=None):

        if name is None:
            name = "cg_group"
        name += ".pdb"
        cg_group = self.create_cg_group(cg_positions)
        cg_group.write(name)

    def write_selection(self, sele_string, name="selected"):

        current_group = mda.Merge(self.all).select_atoms(sele_string)
        #current_group.atoms[0].name = current_group.atoms[0].resname
        current_group.write(name)

    def one_letter_code(self, three_letter_code):

        one = ["A", "G", "I","L","P","V","F","W","Y","D","E","R","H","K","S","T","C","M","N","Q"]
        three = ["ALA","GLY","ILE","LEU","PRO","VAL","PHE","TRP","TYR","ASP","GLU","ARG","HIS","LYS",
                 "SER","THR","CYS","MET","ASN","GLN"]
        return one[three.index(three_letter_code)]

    def write_mstool_names(self, name="mstool"):

        current_group = mda.Merge(self.all).select_atoms("name CA")
        for i in range(len(current_group.atoms)):
            current_group.atoms[i].name = "BB"#self.one_letter_code(current_group.atoms[i].resname)

        current_group.write(name)


    def cg_selection_string(self):

        resids = self.read_cg_map()
        sring = " and ("
        for resid in resids[:-1]:
            sring += " resid " + str(resid) + " or"
        sring += "  resid " + str(resids[-1]) + ")"

        return sring

    def non_cg_selection_string(self):

        resids = self.read_cg_map()
        resids = [i for i in range(1,222) if i not in resids]
        sring = " and ("
        for resid in resids[:-1]:
            sring += " resid " + str(resid) + " or"
        sring += "  resid " + str(resids[-1]) + ")"

        return sring

    def custom_selection_string(self, resids):

        sring = " and ("
        for resid in resids[:-1]:
            sring += " resid " + str(resid) + " or"
        sring += "  resid " + str(resids[-1]) + ")"

        return sring


    def get_best_alignment(self, cg_positions):

        ref = self.create_cg_group(cg_positions)
        if self.cg_map is None:
            ref = ref[:130]
        min_rmsd = 100
        min_positions = []
        min_frame = -1
        for f in range(len(self.u.trajectory)):
            self.u.trajectory[f]
            current_group = mda.Merge(self.all).select_atoms("all")
            if not self.cg_map:
                rmsd, aligned_group = self.alignto(ref=ref, current_group=current_group)
            else:
                rmsd, aligned_group = self.alignto_yaml(ref=ref, current_group=current_group)
            #print(rmsd)
            if rmsd < min_rmsd:
                min_rmsd = rmsd
                min_positions = aligned_group.positions
                min_frame = f
        return min_rmsd, min_positions, min_frame

    def get_worst_alignment(self, cg_positions):

        ref = self.create_cg_group(cg_positions)[:130]
        max_rmsd = 0
        max_positions = []
        max_frame = -1
        for f in range(len(self.u.trajectory)):
            self.u.trajectory[f]
            current_group = mda.Merge(self.all).select_atoms("all")
            rmsd, aligned_group = self.alignto(ref=ref, current_group=current_group)
            # print(rmsd)
            if rmsd > max_rmsd:
                max_rmsd = rmsd
                max_positions = aligned_group.positions
                max_frame = f
        return max_rmsd, max_positions, max_frame

    def best_alignment_group(self, cg_positions, worst=False):

        if not worst:
            min_rmsd, min_positions, min_frame = self.get_best_alignment(cg_positions)
        else:
            min_rmsd, min_positions, min_frame = self.get_worst_alignment(cg_positions)
        print(min_rmsd, min_frame)

        copied = mda.Merge(self.all).atoms
        for ind, atom in enumerate(copied):
            atom.position = min_positions[ind]

        return copied

    def output_best_aligned_ca_group(self, cg_pos_array, name=None):

        group_array = []
        for cg_positions in cg_pos_array:
            group_array.append(self.best_alignment_group(cg_positions))

        current_thing = mda.Merge(group_array[0])
        for i in range(1,len(group_array)):
            current_thing = mda.Merge(current_thing.atoms, group_array[i])

        #for i in range(int(len(current_thing.atoms)/3457)):
        #    new  = current_thing.add_Segment(segid=self.group_id(i))
        #    self.u.residues[ int(i*221) : int((i+1) * 221) ].segments = new
        current_thing = current_thing.atoms
        #for atom in current_thing:
        #    chain_number = int(atom.index /3457)
        #    atom.segment.segid = self.group_id(chain_number)
        if name is None:
            current_thing.write("all.pdb")
        else:
            current_thing.write(name + ".pdb")

    def combine(self, aaprotein, name=None, segids=True):

        import time
        combined = self.u

        if type(aaprotein) == type(["l", "i", "s", "t"]):
            if segids:
                for ind, thing in enumerate(aaprotein):
                    thing.u.add_TopologyAttr('chainID', [self.get_segid(ind)] * len(thing.u.atoms))
                    print(thing.u.atoms.chainIDs)
            for ind, aap in enumerate(aaprotein):
                now = time.time()
                combined = mda.Merge(combined.atoms, aap.u.atoms)
                later = time.time()
                print(ind, later - now)
        else:
            combined = mda.Merge(self.u.atoms, aaprotein.u.atoms)

        if name is None:
            name = "combined"
        name += ".pdb"
        combined.atoms.write(name)

    def get_segid(self, index, include_numbers=False):

        num = index //26

        if num == 0:
            num = ""

        out = str(chr(index % 26 + 97)).upper()

        if include_numbers:
            out = out + str(num)

        return  out

    def combine_binary(self, lizt, name=None, segids=True, include_numbers=False):
        import time
        lizt.append(self)

        if segids:
            for ind, thing in enumerate(lizt):
                thing.u.add_TopologyAttr('chainID', [self.get_segid(ind, include_numbers=include_numbers)]*len(thing.u.atoms))
                print(thing.u.atoms.chainIDs)

        odd = len(lizt) % 2 != 0
        half = int(len(lizt) / 2)
        save = ""
        if odd:
            save = lizt[-1].u
        lizt = [mda.Merge(lizt[i].u.atoms, lizt[i + half].u.atoms) for i in range(half)]
        if odd:
            lizt.append(save)
        while len(lizt) != 1:
            odd = len(lizt) % 2 != 0
            if odd:
                save = lizt[-1]
            half = int(len(lizt)/2)
            now = time.time()
            lizt = [mda.Merge(lizt[i].atoms, lizt[i + half].atoms) for i in range(half)]
            later = time.time()
            print(later - now)
            if odd:
                lizt.append(save)

        if name is None:
            name = "combined"
        name += ".pdb"
        lizt[0].atoms.write(name)


    def output_worst_aligned_ca_group(self, cg_pos_array, name=None):

        group_array = []
        for cg_positions in cg_pos_array:
            group_array.append(self.best_alignment_group(cg_positions, worst=True))

        current_thing = mda.Merge(group_array[0])
        for i in range(1, len(group_array)):
            current_thing = mda.Merge(current_thing.atoms, group_array[i])

        for i in range(int(len(current_thing.atoms) / 3457)):
            new = current_thing.add_Segment(segid=self.group_id(i))
            self.u.residues[int(i * 221): int((i + 1) * 221)].segments = new
        current_thing = current_thing.atoms
        # for atom in current_thing:
        #    chain_number = int(atom.index /3457)
        #    atom.segment.segid = self.group_id(chain_number)
        if name is None:
            current_thing.write("all_worst.pdb")
        else:
            current_thing.write(name + ".pdb")

    def group_id(self, number):

        ids = ["A", "B", "C", "D", "E", "F", 'G', 'H', 'I', 'J', 'K', 'L' 'M',
               'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        if number < len(ids):
            return ids[number]
        else:
            first = ids[int(number%len(ids))]
            second = ids[int(number / len(ids))]
            return first + second


    def alignto(self, cg_positions=None, ref=None, current_group=None):

        if current_group is None:
            copied = mda.Merge(self.all).select_atoms("all")
        else:
            copied = current_group

        if ref is None:
            ref = self.create_cg_group(cg_positions)[:130]

        shifted_map = np.subtract(np.array(self.read_cg_map()), 1)[:130]

        align.alignto(copied, ref, select="name CA" + self.cg_selection_string())
        rmsd = rms.rmsd(copied.select_atoms("name CA")[shifted_map].positions, ref.positions)

        #copied.write("after.pdb")
        return rmsd, copied

    def align_cg_compare_non_cg(self, AAprotein, current_group=None, custom_list=None):


        if current_group is None:
            copied = mda.Merge(self.all).select_atoms("all")
        else:
            copied = current_group

        ref = mda.Merge(AAprotein.all).select_atoms("all")

        old, new = align.alignto(copied, ref, select="name CA" + self.cg_selection_string())
        print(old,new)
        print(self.cg_selection_string())
        sele_string = "name CA" + self.cg_selection_string()

        if custom_list is not None:
            sele_string = "name CA" + self.custom_selection_string(custom_list)


        pos1 = mda.Merge(copied).select_atoms(sele_string).positions
        #print(copied)
        #print(pos1, sele_string)
        pos2 = mda.Merge(ref).select_atoms(sele_string).positions

        rmsd = np.sqrt(np.average([np.linalg.norm(np.subtract(pos1[i], pos2[i])) for i in range(len(pos1))]))
        return rmsd





    def alignto_yaml(self, ref, current_group):

        copied = current_group

        #indexes = [self.cg_map[i][0] for i in range(len(self.cg_map))]
        #indexes = [self.cg_map[i][int(len(self.cg_map[i]) // 2)] for i in range(len(self.cg_map))]
        shifted_map = [np.subtract(cg_map_i, 0) for cg_map_i in self.cg_map ]
        #print(self.cg_map)
        #quit()
        #cg_string = 'index '
        #for thing in shifted_map:
        #    cg_string += str(thing) + " "
        #print(copied[indexes], len(indexes))
        #print(copied.select_atoms(cg_string))
        #align.alignto(copied, ref, select=(cg_string, "all"))
        #rmsd = rms.rmsd(copied.select_atoms("all")[shifted_map].positions, ref.positions)

        new_positions = [np.mean(self.all.positions[site], axis=0) for site in shifted_map]

        new_copied = mda.Universe.empty(len(new_positions), trajectory=True)

        new_copied.atoms.positions = np.array(new_positions)
        new_copied.add_TopologyAttr('segid', ["cg"])

        merged = mda.Merge(new_copied.atoms, copied)
        #print(merged.select_atoms("segid cg"))
        #print(ref)
        align.alignto(merged, ref, select=("segid cg", "all"))
        #align.alignto(copied, ref)
        unmerged = merged.select_atoms("not segid cg")

        check_positions = [np.mean(unmerged.positions[site], axis=0) for site in shifted_map]

        rmsd = rms.rmsd(check_positions, ref.positions)

        #unmerged.write("after.pdb")
        return rmsd, unmerged

    def alignTraj_to_CG(self):

        copied = self.u.trajectory

        ref = self.create_cg_group(self.read_cg_map())[:130]

        shifted_map = np.subtract(np.array(self.read_cg_map()), 1)[:130]

        align.AlignTraj(copied, ref, select="name CA and index 0-3458" + self.cg_selection_string()).run()
        rmsd = rms.rmsd(copied.select_atoms("name CA")[shifted_map].positions, ref.positions)

        # copied.write("after.pdb")
        return rmsd, copied



    def output_cg(self):

        shifted_map = np.subtract(np.array(self.read_cg_map()), 1)[:130]

        positions = self.backbone.positions[shifted_map]

        self.backbone[self.read_cg_map()].write("cg.pdb")

        return positions
