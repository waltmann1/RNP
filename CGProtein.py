from __future__ import division
from CGMolyAbs import CGMolyAbs
import numpy as np
import numpy.linalg as la
import Bio.PDB


class CGProtein(CGMolyAbs):

    def __init__(self, PDBName, domainfile=None, res_domains=None):

        super(CGProtein, self).__init__()
        parser = Bio.PDB.PDBParser()
        self.struct = parser.get_structure(PDBName[:-4], PDBName)
        self.name = PDBName[:-4]
        #self.res_domains = res_domains

        if res_domains is None:
            res_domains = self.read_domain_file(domainfile)

        residues = [residue for residue in self.struct.get_residues()]
        self.domain_starts.append(0)

        for domain_ind, domain in enumerate(res_domains):
            site_index = 0
            domain_f_weights = []
            domain_x_weights = []
            domain_site_indexes = []
            for res_index in domain:
                #print(residues)
                #print(residues[res_index].get_resname())
                #print(domain_ind, res_index)
                atoms = residues[res_index].get_atoms()
                atom_names = [atom.name for atom in atoms]
                domain_site_indexes += list(range(site_index, site_index + len(atom_names)))
                site_index += len(atom_names)
                domain_x_weights += [self.get_x_weight(name) for name in atom_names]
                domain_f_weights += [1 for name in atom_names]
            self.domain_starts.append(site_index + self.domain_starts[-1])
            self.site_indexes.append(domain_site_indexes)
            self.f_weights.append(domain_f_weights)
            self.x_weights.append(domain_x_weights)
        self.domain_starts = self.domain_starts[:-1]


class CGProteinComplex(CGMolyAbs):

    def __init__(self, PDBs, domainfile=None, res_domains=None, name=None):

        super(CGProteinComplex, self).__init__()
        parser = Bio.PDB.PDBParser()
        structs = [parser.get_structure(PDBName[:-4], PDBName) for PDBName in PDBs]
        if name is None:
            self.name = "ProteinComplex"
        else:
            self.name = name
        #self.res_domains = res_domains

        if res_domains is None:
            res_domains = self.read_domain_file(domainfile)

        residues = [residue for struct in structs for residue in struct.get_residues()]
        self.domain_starts.append(0)

        for domain_ind, domain in enumerate(res_domains):
            site_index = 0
            domain_f_weights = []
            domain_x_weights = []
            domain_site_indexes = []
            for res_index in domain:
                #print(residues)
                #print(residues[res_index].get_resname())
                #print(domain_ind, res_index)
                atoms = residues[res_index].get_atoms()
                atom_names = [atom.name for atom in atoms]
                domain_site_indexes += list(range(site_index, site_index + len(atom_names)))
                site_index += len(atom_names)
                domain_x_weights += [self.get_x_weight(name) for name in atom_names]
                domain_f_weights += [1 for name in atom_names]
            self.domain_starts.append(site_index + self.domain_starts[-1])
            self.site_indexes.append(domain_site_indexes)
            self.f_weights.append(domain_f_weights)
            self.x_weights.append(domain_x_weights)
        self.domain_starts = self.domain_starts[:-1]








