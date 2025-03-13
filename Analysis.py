from __future__ import division
import networkx as nx
import numpy as np
import numpy.linalg as la
import os.path
import copy as cp
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 22})
from mpl_toolkits.mplot3d import Axes3D
import math as m
import MDAnalysis as mda
from MDAnalysis.analysis import distances
#from scipy.stats import skew
#from scipy.stats import kurtosis

class Analysis(object):

    def __init__(self, data_file, traj_file, n_dimers=600, ca_length=134, no_ip6 = False, no_rnp=False, format=None):

        if format is None:
            self.u = mda.Universe(data_file, traj_file)
        else:
            self.u = mda.Universe(data_file, traj_file, format=format)
        #self.u = mda.Universe(data_file)
        n_dimers=0
        for ind, atom in enumerate(self.u.atoms):
            if atom.type =="139" and n_dimers == 0:
                n_dimers = int(ind/(ca_length*2))


        self.CAs = [[i * ca_length + j for j in range(ca_length)] for i in range(n_dimers * 2)]
        self.mer_names = []
        self.mer_name_indices = []
        self.cargo = []
        self.frames = []
        self.colors = ["green"]
        self.names = ["CA"]

        self.colors.append("orange")
        self.names.append("IP6")

        self.colors.append("red")
        self.names.append("RNC")

        protein = self.u.atoms
        self.box = self.u.dimensions
        self.ip6_indices = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="140"]
        self.rnc_indices = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="141" or atom.type == "142" or atom.type=="143"]
        self.signal_indices = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="142"]
        self.integrase = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="143"]
        self.integrase_top = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="145"]
        self.integrase_bottom = []
        if len(self.integrase_top) > 0:
            self.integrase_bottom = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="144"]
            self.memb = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="146"]
        else:
            self.memb = [ind for ind, atom in enumerate(self.u.atoms) if atom.type=="144"]

        if no_ip6:
            self.ip6_indices = []
            self.rnc_indices = [ind for ind,atom in enumerate(self.u.atoms) if atom.type == "140" or atom.type=="141" or atom.type=="142"]
            self.signal_indices = [ind for ind,atom in enumerate(self.u.atoms) if atom.type=="141"]
            self.integrase = [ind for ind,atom in enumerate(self.u.atoms) if atom.type =="142"]
            self.integrase_top = [ind for ind, atom in enumerate(self.u.atoms) if atom.type == "144"]
            self.integrase_bottom = []
            if len(self.integrase_top) > 0:
                self.integrase_bottom = [ind for ind, atom in enumerate(self.u.atoms) if atom.type == "143"]
                self.memb = [ind for ind, atom in enumerate(self.u.atoms) if atom.type == "145"]
            else:
                self.memb = [ind for ind, atom in enumerate(self.u.atoms) if atom.type == "143"]

        if no_rnp:
            self.ip6_indices = [ind for ind, atom in enumerate(self.u.atoms) if atom.type == "140"]
            self.memb = [ind for ind, atom in enumerate(self.u.atoms) if atom.type == "141"]

        self.base_graph = nx.Graph()
        for i in range(len(self.CAs)):
            name = "CA"
            self.base_graph.add_node(i, name=name, color=self.colors[self.names.index(name)])

        self.graphs = []
        self.graph_frames = []
        self.traj_file = traj_file


    def get_ca_positions(self, frame_number, ca_index):

        frame = self.u.trajectory[frame_number]

        positions = np.array(self.continuous_poly_positions(frame.positions[self.CAs[ca_index]]))

        #positions = np.subtract(positions, np.average(positions))
        return positions
    """
    def write_ca_positions(self, frame_number, ca_index, name=None):
        
        positions = self.get_ca_positions(frame_number, ca_index)
        
        if name is None:
            name = "positions"
        name += ".pdb"    
            
        f = open(name, "w")
        for pos in positions:
            f.write(str(pos[0]) + " " +str(pos[1]) + " " + str(pos[2]) + " ")
        
        """

    def add_dimer_edges(self, graph):

        for i in range(len(self.CAs)):
            if i % 2 == 1:
                graph.add_edge(i-1, i)

    def remove_dimer_edges(self, graph):

        for i in range(len(self.CAs)):
            if i % 2 == 1:
                graph.remove_edge(i-1, i)

    def network_grapher(self, G, name="g"):

        fig = plt.figure(figsize=(10, 10))
        ax1 = fig.add_subplot(111)
        ax1.set_title(name)
        color_map = [G.nodes[i]['color'] for i in G.nodes]
        nx.draw_networkx(G, with_labels=True, node_color=color_map)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
        for pos in ['right', 'top', 'bottom', 'left']:
            plt.gca().spines[pos].set_visible(False)

        plt.savefig("network_" + str(name) + ".png")
        plt.show()

    def add_two_edges(self, frame_number, graph, cut=12):

        frame = self.u.trajectory[frame_number]

        two_positions = np.array([frame.positions[ca[1]] for ca in self.CAs])

        all_distances = distances.distance_array(two_positions, two_positions, box=self.box)

        for i in range(len(self.CAs)):
            for j in range(i + 1, len(self.CAs)):
                if all_distances[i][j] < cut:
                    graph.add_edge(i, j)


    def add_thirty_edges(self, frame_number, graph, cut=20):

        frame = self.u.trajectory[frame_number]

        positions_132 = np.array([frame.positions[ca[131]] for ca in self.CAs])
        positions_133 = np.array([frame.positions[ca[132]] for ca in self.CAs])
        positions_31 = np.array([frame.positions[ca[30]] for ca in self.CAs])
        positions_33 = np.array([frame.positions[ca[32]] for ca in self.CAs])

        distances_31_132 = distances.distance_array(positions_31, positions_132, box=self.box)
        distances_33_133 = distances.distance_array(positions_33, positions_133, box=self.box)

        for i in range(len(self.CAs)):
            for j in range(len(self.CAs)):
                if i !=j and (distances_31_132[i][j] < cut or distances_33_133[i][j]< cut):
                    graph.add_edge(i, j)

    def get_all_in_cycles(self, frame_number, dimers=True, separate_adsorbed=True):

        new_graph = cp.deepcopy(self.base_graph)
        self.add_thirty_edges(frame_number, new_graph)
        if dimers:
            self.add_dimer_edges(new_graph)

        cycles = list(nx.simple_cycles(new_graph, length_bound=6))
        cycles = [cycle for cycle in cycles if len(cycle) > 4]
        ads_cycle = []
        un_cycle = []

        if not separate_adsorbed:
            return cycles

        ads = self.get_all_adsorbed(frame_number)[0]
        for cycle in cycles:
            adsorbed = False
            for node in cycle:
                if node in ads:
                    adsorbed = True
            if adsorbed:
                ads_cycle.append(cycle)
            else:
                un_cycle.append(cycle)

        return ads_cycle, un_cycle


    def largest_subgraph(self, graph, min=12, ads_color=None, one=False):

        groups = list(nx.connected_components(graph))

        keeps = []
        if not one:
            for group in groups:
                if len(group) > min:
                    keeps.extend(list(group))
        else:
            keeps = groups[0]
            max = len(keeps)
            for group in groups[1:]:
                if len(group) > max:
                    keeps = group
                    max = len(group)

        new_graph = cp.deepcopy(graph)

        sub = new_graph.subgraph(keeps)

        return sub

    def color_nodes(self, graph, nodes, color="red"):

        for i in nodes:
            graph.nodes[i]['color'] = color



    def display_graph(self, frame_number, color_adsorbed =False, largest_sub=False, dimer_edges=False, IP6=False,
                      rnc=False, thirty_cut=20, seed=None, twos=False):

        new_graph = self.get_graph_by_frame_number(frame_number)

        if dimer_edges:
            self.add_dimer_edges(new_graph)
        if IP6:
            self.add_ip6_connections(frame_number, new_graph)
        if rnc:
            self.add_polymer_connections(frame_number, new_graph)
        if twos:
            self.add_two_edges(frame_number, new_graph)
        else:
            self.add_thirty_edges(frame_number, new_graph, cut=thirty_cut)

        if color_adsorbed:
            ads = self.get_all_adsorbed(frame_number)[0]
            self.color_nodes(new_graph, ads)
        if seed is not None:
            if type(seed) is not list:
                seed = [seed]
            self.color_nodes(new_graph, seed, color="yellow")

        if largest_sub:
            new_graph = self.largest_subgraph(new_graph)

        self.network_grapher(new_graph)

    def add_ip6_connections(self, frame_number, graph, cut=15):

        if graph.number_of_nodes() == len(self.CAs):
            for ind in self.ip6_indices:
                name = "IP6"
                graph.add_node(graph.number_of_nodes(), name=name, color=self.colors[self.names.index(name)])
        else:
            return

        frame = self.u.trajectory[frame_number]

        two_positions = np.array([frame.positions[ca[1]] for ca in self.CAs])

        ip6_positions = np.array([frame.positions[ind] for ind in self.ip6_indices])

        all_distances = distances.distance_array(ip6_positions, two_positions, box=self.box)

        for i, index in enumerate(ip6_positions):
            for j in range(len(two_positions)):
                if all_distances[i][j] < cut:
                    graph.add_edge(len(self.CAs) + i, j)


    def add_polymer_connections(self, frame_number, graph, cut=15, integrase_cut=36, top_cut=31):


        poly_number = graph.number_of_nodes()
        graph.add_node(poly_number, name="RNP", color="red")

        frame = self.u.trajectory[frame_number]

        #pos_positions = np.array([frame.positions[ca[133]] for ca in self.CAs])
        pos_positions = np.array([frame.positions[ca[114]] for ca in self.CAs])
        rnc_positions = []
        if len(self.integrase_top) > 0:
            rnc_positions = np.array([frame.positions[ind] for ind in self.integrase_top])
            cut = top_cut
        elif len(self.integrase) > 0:
            rnc_positions = np.array([frame.positions[ind] for ind in self.integrase])
            cut = integrase_cut
        else:
            rnc_positions = np.array([frame.positions[ind] for ind in self.rnc_indices])
        all_distances = distances.distance_array(pos_positions, rnc_positions, box=self.box)

        for i, index in enumerate(pos_positions):
                if np.any(all_distances[i] < cut):
                    graph.add_edge(i, poly_number)

        #if len(self.integrase) > 0:
         #   int_positions = np.array([frame.positions[ind] for ind in self.integrase])
         #   int_distances = distances.distance_array(pos_positions, int_positions, box=self.box)

          #  for i, index in enumerate(pos_positions):
          #      if np.any(int_distances[i] < integrase_cut):
          #          graph.add_edge(i, poly_number)

    def get_all_in_clusters(self, frame_number, min=4, dimers=True):

        new_graph = cp.deepcopy(self.base_graph)
        self.add_thirty_edges(frame_number, new_graph)
        if dimers:
            self.add_dimer_edges(new_graph)

        groups = list(nx.connected_components(new_graph))
        string = ""
        nodes = []
        for ind, group in enumerate(groups):
            if len(group) > min:
                nodes.append(group)
                for thing in group:
                    for i in range(134*thing, 134*(thing+1) + 1, 33):
                        string += str(i) + " "
        string = "index " + string
        f = open(self.abs_path("clusters.txt"), "w")
        f.write(string)
        return nodes

    def get_all_ads_clusters(self, frame_number, min=2, dimers=True):

        clusters = self.get_all_in_clusters(frame_number, min=min, dimers=dimers)

        ads = self.get_all_adsorbed(frame_number)[0]

        groups = []

        for cluster in clusters:
            cluster = list(cluster)
            for thing in cluster:
                if thing in ads and cluster not in groups:
                    groups.append(cluster)

        string = ""
        for ind, group in enumerate(groups):
            if len(group) > min:
                for thing in group:
                    for i in range(134*thing, 134*(thing+1) + 1, 33):
                        string += str(i) + " "
        string = "index " + string
        f = open(self.abs_path("ads_clusters.txt"), "w")
        f.write(string)
        return groups

    def get_all_heptamers(self, frame_number):
        new_graph = cp.deepcopy(self.base_graph)
        self.add_thirty_edges(frame_number, new_graph)

        groups = list(nx.connected_components(new_graph))
        string = ""
        nodes = []
        for ind, group in enumerate(groups):
            if len(group) == 7:
                neighbors = [len(list(new_graph.neighbors(thing))) for thing in group]
                neighbors = np.array(neighbors)
                if np.sum(neighbors) == 14 and len(set(neighbors)) == 1:
                    nodes.append(group)
                    options = [101, 30, 32, 131, 132]
                    for thing in group:
                        # for i in options:
                        for i in range(134 * thing, 134 * (thing + 1) + 1, 33):
                            string += str(i) + " "

        string = "index " + string
        f = open(self.abs_path("hepts.txt"), "w")
        f.write(string)
        return nodes


    def get_all_hexamers(self, frame_number):

        new_graph = cp.deepcopy(self.base_graph)
        self.add_thirty_edges(frame_number, new_graph)

        groups = list(nx.connected_components(new_graph))
        string = ""
        nodes = []
        for ind, group in enumerate(groups):
            #if 491 in group:
            #    print("found it", group)
            #    print([len(list(new_graph.neighbors(thing))) for thing in group])
            #    print([list(new_graph.neighbors(thing)) for thing in group])
            if len(group) == 6:
                neighbors = [len(list(new_graph.neighbors(thing))) for thing in group]
                #print(neighbors)
                if np.sum(neighbors) == 12:
                    nodes.append(group)
                    options = [101, 30, 32, 131, 132]
                    for thing in group:
                        #for i in options:
                        for i in range(134 * thing, 134 * (thing + 1) + 1, 33):
                            string += str(i) + " "
            elif len(group) == 7:
                neighbors = [len(list(new_graph.neighbors(thing))) for thing in group]
                neighbors= np.array(neighbors)
                if np.sum(neighbors) == 14:
                    print("neighbors", neighbors)
                    print("find one", neighbors==1, np.sum(neighbors==1), np.sum(neighbors==1) == 1 )
                    print("find three", neighbors==3, np.sum(neighbors==3)==1, np.sum(neighbors==1) == 1)
                if np.sum(neighbors) == 14 and np.sum(neighbors==1) == 1 and np.sum(neighbors==3)==1:
                    print("going")
                    group = [thing for ind, thing in enumerate(group) if neighbors[ind] != 1]
                    nodes.append(group)
                    options = [101, 30, 32, 131, 132]
                    for thing in group:
                        # for i in options:
                        for i in range(134 * thing, 134 * (thing + 1) + 1, 33):
                            string += str(i) + " "

        string = "index " + string
        f = open(self.abs_path("hexes.txt"), "w")
        f.write(string)
        return nodes

    def get_all_hexamers_with_ip6(self, frame_number):

        hexes = self.get_all_hexamers(frame_number)
        bound = []
        for hexamer in hexes:
            if self.check_bound_ip6(hexamer, frame_number) != -1:
                bound.append(hexamer)
        return bound

    def get_all_pentamers_with_ip6(self, frame_number, strong=True):

        pents = self.get_all_pentamers(frame_number, strong=strong)
        bound = []
        for pentamer in pents:
            if self.check_bound_ip6(pentamer, frame_number) != -1:
                bound.append(pentamer)
        return bound


    def get_ip6_position(self, frame_number, ip6_index):

        frame = self.u.trajectory[frame_number]
        return frame.positions[self.ip6_indices[ip6_index]]



    def check_bound_ip6(self, cas,  frame_number, cut=15):

        frame = self.u.trajectory[frame_number]
        cas = [self.CAs[index] for index in cas]

        two_positions = np.array([frame.positions[ca[1]] for ca in cas])
        ip6_positions = np.array([frame.positions[ind] for ind in self.ip6_indices])
        all_distances = distances.distance_array(ip6_positions, two_positions, box=self.box)

        for i, index in enumerate(ip6_positions):
            truth = [int(all_distances[i][j] < cut) for j in range(len(two_positions))]
            if np.sum(truth) == len(two_positions):
                    print(all_distances[i])
                    return i
        return -1

    def get_all_adsorbed_hexamers(self, frame_number):

        hexes = self.get_all_hexamers(frame_number)
        ads_clust = self.get_all_ads_clusters(frame_number)
        ads = [thing for clust in ads_clust for thing in clust]
        ads_hex = []
        for hex in hexes:
            if list(hex)[0] in ads:
                ads_hex.append(hex)
        return ads_hex

    def get_all_pentamers(self, frame_number, strong=True):

        new_graph = cp.deepcopy(self.base_graph)
        self.add_thirty_edges(frame_number, new_graph)

        groups = list(nx.connected_components(new_graph))
        string = ""
        nodes = []
        for ind, group in enumerate(groups):
            if len(group) == 5:
                neighbors = [len(list(new_graph.neighbors(thing))) for thing in group]
                #print(neighbors)
                if np.sum(neighbors) == 10 or not strong:
                    nodes.append(group)
                    for thing in group:
                        for i in range(134 * thing, 134 * (thing + 1) + 1, 33):
                            string += str(i) + " "
        string = "index " + string
        f = open(self.abs_path("pents.txt"), "w")
        f.write(string)
        return nodes

    def get_all_adsorbed(self, frame_number,cut=15):

        new_graph = cp.deepcopy(self.base_graph)
        self.add_polymer_connections(frame_number, new_graph, cut=cut)

        groups = list(nx.connected_components(new_graph))
        string = ""
        nodes = []
        for ind, group in enumerate(groups):
            if len(group) > 1:
                group = list(group)
                group.sort()
                group = group[:-1]
                nodes.append(group)
                for thing in group:
                    for i in range(134 * thing, 134 * (thing + 1) + 1, 33):
                        string += str(i) + " "
        string = "index " + string
        f = open(self.abs_path("all_adsorbed.txt"), "w")
        f.write(string)
        return nodes

    def get_graph_by_frame_number(self, frame_number):

        if frame_number not in self.graph_frames:
            self.graph_frames.append(frame_number)
            new_graph = cp.deepcopy(self.base_graph)
            self.graphs.append(new_graph)
        else:
            new_graph = self.graphs[self.graph_frames.index(frame_number)]
        return new_graph

    def abs_path(self, string):
        return "/home/cwaltmann/PycharmProjects/Builders/" + string


    def write_spherical_xyz(self, frame_number, name="out.xyz", add_rnc=True, add_adsorbed=True, add_clustered=True,
                            add_ip6=True, cont=True, cont_ca=True, strong_pent=True, memb=False):

        frame = self.u.trajectory[frame_number]
        if add_adsorbed:
            on = self.get_all_adsorbed(frame_number)
            if len(on) > 0:
                on = on[0]
        else:
            on = []

        split_hexes = self.get_all_hexamers(frame_number)
        hexes = [thing for group in split_hexes for thing in group]
        split_pents = self.get_all_pentamers(frame_number, strong=strong_pent)
        pents = [thing for group in split_pents for thing in group]

        clustered = []
        if add_clustered:
            clusters = self.get_all_in_clusters(frame_number)
            keep = []
            for group in clusters:
                for thing in group:
                    if thing in hexes or thing in pents and group not in keep:
                        keep.append(group)
            for group in keep:
                for thing in group:
                    if thing not in hexes and thing not in pents and thing not in clustered:
                        clustered.append(thing)

        indexes_hexes = []
        indexes_pents = []
        indexes_on = []
        indexes_clustered = []
        for thing in hexes:
            for i in range(134 * thing, 134 * (thing + 1) + 1, 33):
                indexes_hexes.append(i)
        for thing in pents:
            for i in range(134 * thing, 134 * (thing + 1) + 1, 33):
                indexes_pents.append(i)
        if add_adsorbed:
            for thing in on:
                #if thing not in hexes and thing not in pents:
                for i in list(range(134 * thing, 134 * (thing + 1) + 1, 33)):
                    indexes_on.append(i)
        if add_clustered:
            for thing in clustered:
                for i in range(134*thing, 134 * (thing + 1), 33):
                    indexes_clustered.append(i)

        hex_positions = np.array([frame.positions[thing] for thing in indexes_hexes])

        pent_positions = np.array([frame.positions[thing] for thing in indexes_pents])

        on_positions = np.array([frame.positions[thing] for thing in indexes_on])

        clustered_positions = np.array([frame.positions[thing] for thing in indexes_clustered])

        box_edges = [0, self.box[0]]
        box_positions = []
        for i in box_edges:
            for j in box_edges:
                for k in box_edges:
                    box_positions.append([i, j, k])

        poly_positions = []
        signal_positions = []
        integrase_positions = []
        top_positions = []
        bottom_positions = []
        rnc_indices_normal = [index for index in self.rnc_indices if index not in self.signal_indices and index not in self.integrase]

        if add_rnc:
            poly_positions = [frame.positions[index] for index in rnc_indices_normal if index % 5==0]
            signal_positions = [frame.positions[index] for index in self.signal_indices if index % 5==0]
            integrase_positions = [frame.positions[index] for index in self.integrase]
            top_positions = [frame.positions[index] for index in self.integrase_top]
            bottom_positions = [frame.positions[index] for index in self.integrase_bottom]
            if cont:
                combined = cp.deepcopy(poly_positions)
                combined.extend(signal_positions)
                combined = self.continuous_poly_positions(combined)
                hex_positions = self.compute_min_image_positions(combined, hex_positions)
                pent_positions = self.compute_min_image_positions(combined, pent_positions)
                on_positions = self.compute_min_image_positions(combined, on_positions)
                clustered_positions = self.compute_min_image_positions(combined, clustered_positions)
                integrase_positions = self.compute_min_image_positions(combined, integrase_positions)
                top_positions = self.compute_min_image_positions(combined, top_positions)
                bottom_positions = self.compute_min_image_positions(combined, bottom_positions)
                poly_positions= combined[:len(poly_positions)]
                signal_positions = combined[len(poly_positions):]

        if cont_ca and len(clustered) > 0:
            seed_ca = clustered[0]
            for thing in on:
                if thing in clustered or thing in hexes or thing in pents:
                    seed_ca = thing

            cont_pos = self.compute_continuous_ca_positions(frame_number, seed_ca, combined=combined)
            hex_positions = np.array([cont_pos[thing] for thing in indexes_hexes])
            pent_positions = np.array([cont_pos[thing] for thing in indexes_pents])
            clustered_positions = np.array([cont_pos[thing] for thing in indexes_clustered])


        ip6_positions = []
        if add_ip6:
            for oligomer in split_hexes + split_pents:
                possible = self.check_bound_ip6(oligomer, frame_number)
                if possible != -1:
                    ip6_positions.append(frame.positions[self.ip6_indices[possible]])
            if cont and not cont_ca:
                ip6_positions = self.compute_min_image_positions(combined, ip6_positions)
            if cont_ca:
                temp = list(cp.deepcopy(hex_positions))
                temp.extend(list(pent_positions))
                ip6_positions = self.compute_min_image_positions(temp,ip6_positions)
        memb_positions = []

        if memb:
            memb_positions = np.array([frame.positions[i] for i in self.memb if i%10==0])



        total = (((len(hex_positions) + len(pent_positions)+ len(on_positions) +
                 len(box_positions)) + len(poly_positions) + len(clustered_positions) + len(ip6_positions)
                 + len(signal_positions) + len(integrase_positions)) + len(memb_positions) + len(bottom_positions)
                 + len(top_positions))
        f = open(name, "w")
        f.write(str(total) + "\n\n")
        count = 0

        if cont and add_rnc:
            center = np.average(poly_positions, axis=0)
            box_center = np.divide(self.box[:3], 2)
            shift = np.subtract(box_center, center)
            if memb:
                shift = [0,0,0]
            poly_positions = self.safe_add(shift, poly_positions)
            hex_positions = self.safe_add(shift, hex_positions)
            pent_positions = self.safe_add(shift, pent_positions)
            on_positions = self.safe_add(shift, on_positions)
            clustered_positions = self.safe_add(shift, clustered_positions)
            ip6_positions = self.safe_add(shift, ip6_positions)
            signal_positions = self.safe_add(shift, signal_positions)
            integrase_positions = self.safe_add(shift, integrase_positions)
            top_positions = self.safe_add(shift, top_positions)
            bottom_positions = self.safe_add(shift, bottom_positions)

        for pos in hex_positions:
            f.write(self.xyz_string(count, pos,"H"))
            count += 1
        for pos in pent_positions:
            f.write(self.xyz_string(count, pos, "P"))
            count += 1
        for pos in on_positions:
            f.write(self.xyz_string(count, pos, "O"))
            count += 1
        for pos in clustered_positions:
            f.write(self.xyz_string(count, pos, "C"))
            count += 1
        for pos in box_positions:
            f.write(self.xyz_string(count, pos, "B"))
            count += 1
        for pos in poly_positions:
            f.write(self.xyz_string(count, pos, "R"))
            count += 1
        for pos in ip6_positions:
            f.write(self.xyz_string(count, pos, "I"))
            count += 1
        for pos in signal_positions:
            f.write(self.xyz_string(count, pos, "S"))
            count += 1
        for pos in integrase_positions:
            f.write(self.xyz_string(count, pos, "N"))
            count += 1
        for pos in top_positions:
            f.write(self.xyz_string(count, pos, "Nt"))
            count += 1
        for pos in bottom_positions:
            f.write(self.xyz_string(count, pos, "Nb"))
            count += 1
        for pos in memb_positions:
            f.write(self.xyz_string(count, pos, "M"))
            count += 1
        f.close()

    def safe_add(self, shift, positions):

        if len(positions) == 0:
            return positions
        else:
            return np.add(shift, positions)

    def xyz_string(self, count, pos, type):
        for ind, dig in enumerate(pos):
            if dig == 0:
                dig = 0.00
            pos[ind] = round(dig, 2)
        return str(count) + " " + str(type) + " " + str(pos[0])+ " " + str(pos[1]) + " " + str(pos[2]) + "\n"

    def get_center(self, frame_number):

        frame = self.u.trajectory[frame_number]
        rnc_positions = np.array([frame.positions[ind] for ind in self.rnc_indices])
        origin = np.array([300,500,300])
        #origin = rnc_positions[100]
        min_positions = self.compute_min_image_positions(origin, rnc_positions)
        return np.average(min_positions, axis=0)

    def compute_min_image_positions(self, origin_positions, positions):
        if len(positions) == 0:
            return positions
        min_images = [[0, 0, 0] for _ in range(len(positions))]
        #print(origin)
        #print(positions)
        positions = np.array(positions)
        origin_positions = np.array(origin_positions)
        min_dist = distances.distance_array(positions, origin_positions)
        min_dist = [np.min(dists) for dists in min_dist]
        for i in range(-1, 2, 1):
            for j in range(-1, 2, 1):
                for k in range(-1, 2, 1):
                    image = [self.box[0] * i, self.box[1] * j, self.box[2] * k]
                    new_positions = np.array([np.add(positions[i], image) for i in range(len(positions))])
                    all_distances = distances.distance_array(new_positions, origin_positions)
                    for l in range(len(positions)):
                        #print(min_dist[l], all_distances[l])
                        if np.min(all_distances[l]) < min_dist[l]:
                            #print(all_distances[l], min_dist[l], "update", i, j, k)
                            min_images[l] = [i, j, k]
                            min_dist[l] = np.min(all_distances[l])
        return [np.add(positions[z], np.multiply(self.box[:3], min_images[z]) ) for z in range(len(positions))]

    def lattice_eccentricity(self, frame, normalized=True):

        g = self.get_graph_by_frame_number(frame)
        self.add_thirty_edges(frame, g)
        self.add_dimer_edges(g)

        g = self.largest_subgraph(g, one=True)
        ecc = nx.eccentricity(g)
        if normalized:
            ecc = list(ecc.values())
            ecc = np.divide(ecc, np.max(ecc))
        return ecc

    def graph_eccentricity(self, frame, normalized=True):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ecc = self.lattice_eccentricity(frame, normalized=normalized)
        bin_width = .033
        n_bins = int((max(ecc) - min(ecc)) / bin_width)
        hist, bin_edges = np.histogram(ecc, bins=n_bins)

        hist = np.divide(hist, np.sum(hist))

        bin_middles = [(bin_edges[i] + bin_edges[i+1])/2 for i in range(len(hist))]

        ax1.bar(bin_middles, hist, width=0.8 * bin_width)

        print("kurtosis:", kurtosis(hist))
        print("skew:", skew(hist))
        #ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        #ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("ecc.png", bbox_inches="tight")
        plt.show()

    def compute_average_min_image_positions(self, origin_positions, positions):
        if len(positions) == 0:
            return positions
        min_images = [[0, 0, 0] for _ in range(len(positions))]
        #print(origin)
        #print(positions)
        positions = np.array(positions)
        origin_positions = np.array(origin_positions)
        orig_dist = distances.distance_array(positions, origin_positions)
        min_dist = [np.min(dists) for dists in orig_dist]
        min_dist = np.average(min_dist)
        min_image = [0,0,0]
        for i in range(-1, 2, 1):
            for j in range(-1, 2, 1):
                for k in range(-1, 2, 1):
                    image = [self.box[0] * i, self.box[1] * j, self.box[2] * k]
                    new_positions = np.array([np.add(positions[x], image) for x in range(len(positions))])
                    all_distances = distances.distance_array(new_positions, origin_positions)
                    averge_min_dist = np.average([np.min(dists) for dists in all_distances])
                    if averge_min_dist < min_dist:
                        min_image = [i, j, k]
                        min_dist = averge_min_dist
        return [np.add(positions[z], np.multiply(self.box[:3], min_image)) for z in range(len(positions))]

    def compute_continuous_ca_positions(self, frame_number, seed_ca, combined=False):

        frame = self.u.trajectory[frame_number]
        new_graph = cp.deepcopy(self.base_graph)
        self.add_thirty_edges(frame_number, new_graph)
        self.add_dimer_edges(new_graph)
        ca_length = len(self.CAs[0])
        new_positions = cp.deepcopy(frame.positions[:int(len(self.CAs) * ca_length)])

        if combined:
            new_positions[ca_length * seed_ca: ca_length * (seed_ca + 1)] = (
                self.compute_min_image_positions(combined, new_positions[ca_length * seed_ca: ca_length * (seed_ca + 1)]))
            new_positions[ca_length * seed_ca: ca_length * (seed_ca + 1)] = self.continuous_poly_positions(
                new_positions[ca_length * seed_ca: ca_length * (seed_ca + 1)])


        queue = []
        visited = [False for _ in range(len(self.CAs))]
        queue.append(seed_ca)
        visited[seed_ca] = True
        while queue:
            s = queue.pop(0)
            for neigh in new_graph.neighbors(s):
                if visited[neigh] == False:
                    queue.append(neigh)
                    visited[neigh] = True
                    new_positions[ca_length*neigh] = self.periodic_fix(new_positions[ca_length * s], new_positions[ca_length * neigh])
                    new_positions[ca_length * neigh: ca_length *(neigh+1)] = self.continuous_poly_positions(new_positions[ca_length * neigh: ca_length *(neigh+1)])
        return new_positions



    def continuous_poly_positions(self, poly_positions):

        pos = cp.deepcopy(poly_positions)
        for i in range(1, len(poly_positions)):
            #print(pos[i])
            pos[i] = self.periodic_fix(pos[i-1], pos[i])
            #print("now", pos[i])
        return pos

    def periodic_fix(self, origin_pos, check_pos):

        c_pos = cp.deepcopy(check_pos)
        diff = np.subtract(origin_pos, c_pos)
        #print(diff)
        for i in range(3):
            if diff[i] > self.box[i]/2:
                c_pos[i] += self.box[i]
            elif diff[i] < - self.box[i]/2:
                c_pos[i] -= self.box[i]
        return c_pos

    def write_trajectory_frames(self, frames):

        for frame_number in frames:
            frame = self.u.trajectory[frame_number]
            with mda.Writer("clean_" + self.traj_file) as W:
                W.write(self.u.atoms)

class ReadOut(object):

    def __init__(self, frames, path="out"):

        self.path=path
        self.letters = ["H", "P", "O", "C", "B", "R", "I", "S", "N", "M", "Nb", "Nt"]
        self.frame_numbers = frames
        self.frames = [self.read_frame(frame_number) for frame_number in frames]
        self.beads_per_ca = 5

    def read_frame(self, frame_number):

        #print(self.path)
        f = open(self.path+"/out_" + str(frame_number) + ".xyz")
        positions = [[] for _ in range(len(self.letters))]
        data = f.readlines()[2:]
        for line in data:
            s = line.split()
            positions[self.letters.index(s[1])].append([float(s[2]), float(s[3]), float(s[4])])

        return positions

    def get_number_hexamers(self, frame):

        return int(len(self.frames[self.frame_numbers.index(frame)][self.letters.index("H")]) / (6 * self.beads_per_ca))

    def get_number_pentamers(self, frame):
        return int(len(self.frames[self.frame_numbers.index(frame)][self.letters.index("P")]) / (5 * self.beads_per_ca))

    def get_number_edge_middle_pentamers(self, frame):
        em_pents = self.get_edge_middle_pentamers(frame)
        return [len(em_pents[0]), len(em_pents[1])]

    def get_num_ip6(self, frame):

        return len(self.frames[self.frame_numbers.index(frame)][self.letters.index("I")])

    def get_hexamers(self, frame):

        all_ca_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("H")]
        total = int(len(all_ca_pos) / (6 * self.beads_per_ca))
        hexamers = [[] for _ in range(total)]

        for ind, ca in enumerate(all_ca_pos):
            hexamers[int(ind / (6 * self.beads_per_ca))].append(ca)
        return hexamers

    def get_pentamers(self, frame):

        all_ca_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("P")]
        total = int(len(all_ca_pos) / (5 * self.beads_per_ca))
        pentamers = [[] for _ in range(total)]

        for ind, ca in enumerate(all_ca_pos):
            pentamers[int(ind / (5 * self.beads_per_ca))].append(ca)
        return pentamers

    def get_edge_middle_pentamers(self, frame, cut = 135):

        pents = self.get_pentamers(frame)
        hexes = self.get_hexamers(frame)

        edges = []
        middles = []

        cent_pent = [np.average(pents[i], axis=0) for i in range(len(pents))]
        cent_hexes = [np.average(hexes[i], axis=0) for i in range(len(hexes))]

        if len(pents) > 0:
            #print(cent_pent)
            #print(len(pents[0]),pents)
            all_distances = distances.distance_array(np.array(cent_pent), np.array(cent_hexes+cent_pent))

            for i in range(len(cent_pent)):
                if (all_distances[i] < cut).sum() == 6:
                    middles.append(pents[i])
                else:
                    edges.append(pents[i])
        return edges, middles

    def tetrahedron_volume(self, a, b , c, d):

        return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d)))/6

    def get_capsid_volume(self, frame, sa_correction=True):

        all_hex_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("H")]
        #hex_pos = [np.average(all_hex_pos[6*i: 6*(i+1)], axis=0)  for i in list(range(int(len(all_hex_pos)/(6 * self.beads_per_ca))))]
        hex_pos = [np.average(all_hex_pos[6 * self.beads_per_ca * i: 6 * self.beads_per_ca * (i + 1)], axis=0) for i in
                   list(range(int(len(all_hex_pos) / (6 * self.beads_per_ca))))]


        all_pent_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("P")]
        pent_pos = [np.average(all_pent_pos[5 * self.beads_per_ca* i: 5 * self.beads_per_ca * (i + 1)], axis=0) for i in
                    list(range(int(len(all_pent_pos)/( 5 * self.beads_per_ca))))]

        #hex_pos.extend(pent_pos)
        points = hex_pos + pent_pos

        #from scipy.spatial import Delaunay
        from scipy.spatial import ConvexHull

        ch = ConvexHull(points)
        return ch.volume - (int(sa_correction) * 10 * ch.area)

    def get_center_distance_distribution_capsid(self, frame):

        all_hex_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("H")]
        # hex_pos = [np.average(all_hex_pos[6*i: 6*(i+1)], axis=0)  for i in list(range(int(len(all_hex_pos)/(6 * self.beads_per_ca))))]
        hex_pos = [np.average(all_hex_pos[6 * self.beads_per_ca * i: 6 * self.beads_per_ca * (i + 1)], axis=0) for i in
                   list(range(int(len(all_hex_pos) / (6 * self.beads_per_ca))))]

        all_pent_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("P")]
        pent_pos = [np.average(all_pent_pos[5 * self.beads_per_ca * i: 5 * self.beads_per_ca * (i + 1)], axis=0) for i
                    in
                    list(range(int(len(all_pent_pos) / (5 * self.beads_per_ca))))]
        points = hex_pos + pent_pos

        center = np.average(points, axis=0)
        distances = np.linalg.norm(np.subtract(points, center), axis=1)

        return points, distances

    def write_capsid_shape(self, frame_number):

        points, dists = self.get_center_distance_distribution_capsid(frame_number)
        f = open("out/center_dist.xyz", "w")
        f.write(str(len(dists)) + "\n\n")
        for i in range(len(points)):
            f.write(self.xyz_string(i, points[i], int(dists[i])))
        f.close()

    def write_cylinder(self):

        points, dists = self.get_cylinder_points()
        f = open("out/cylinder.xyz", "w")
        f.write(str(len(dists)) + "\n\n")
        for i in range(len(points)):
            f.write(self.xyz_string(i, points[i], int(dists[i])))
        f.close()

    def get_cylinder_points(self):
        radius = 200
        length = 500
        spacing  = 1
        circle_points = [[radius * np.cos(np.deg2rad( 15*i)), radius * np.sin(np.deg2rad(15*i)), 0]  for i in range(24)]
        points=cp.deepcopy(circle_points)
        for i in range(1,int(1 + length/spacing)):
            points.extend(np.add(circle_points, [0,0,i*spacing]))
            points.extend(np.add(circle_points, [0,0,-i*spacing]))
        return points, np.linalg.norm(points, axis=1)

    def xyz_string(self, count, pos, type):
        for ind, dig in enumerate(pos):
            if dig == 0:
                dig = 0.00
            pos[ind] = round(dig, 2)
        return str(count) + " " + str(type) + " " + str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + "\n"


    def graph_cylinder(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        positions, dists = self.get_cylinder_points()
        dists = np.divide(dists, 10)
        # print(dists)
        bin_width = 2
        n_bins = int((max(dists) - min(dists)) / bin_width)
        hist, bin_edges = np.histogram(dists, bins=n_bins)
        hist = np.divide(hist, np.sum(hist))

        bin_middles = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(hist))]

        ax1.bar(bin_middles, hist, width=0.8 * bin_width)

        # print("kurtosis:", kurtosis(hist))
        # print("skew:", skew(hist))
        # ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        # ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("cylinder.png", bbox_inches="tight")
        plt.show()
    def graph_capsid_shape(self, frame):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        positions, dists = self.get_center_distance_distribution_capsid(frame)
        dists = np.divide(dists,10)
        #print(dists)
        bin_width = 2
        n_bins = int((max(dists) - min(dists)) / bin_width)
        hist, bin_edges = np.histogram(dists, bins=n_bins)
        hist = np.divide(hist, np.sum(hist))

        bin_middles = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(hist))]

        ax1.bar(bin_middles, hist, width=0.8 * bin_width)

        #print("kurtosis:", kurtosis(hist))
        #print("skew:", skew(hist))
        # ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        # ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("shape_"+ str(frame) +".png", bbox_inches="tight")
        plt.show()


    def circular_variance(self, vectors):

        vectors = [np.divide(v, np.linalg.norm(v)) for v in vectors]
        v = vectors[0]
        for i in range(1, len(vectors)):
            v = np.add(v, vectors[i])
        return 1 - (np.linalg.norm(v) /len(vectors))

    def integrase_cv(self, frame):

        total = 0
        in_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("N")]
        in_pos = [in_pos[i*4] for i in range(len(in_pos)) if i * 4 < len(in_pos)]

        if len(self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]) > 0:
            in_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]
            #in_pos = [in_pos[i * 2] for i in range(len(in_pos)) if i * 2 < len(in_pos)]
        #print(len(in_pos))
        rnp_pos = (self.frames[self.frame_numbers.index(frame)][self.letters.index("R")] +
                   self.frames[self.frame_numbers.index(frame)][self.letters.index("S")])

        all_distances = distances.distance_array(np.array(in_pos), np.array(rnp_pos))
        cut_off = 100
        cvs = []
        for i in range(len(in_pos)):
            #print(in_pos[i])
            #print(all_distances[i])
            #print(np.argwhere(all_distances[i] < cut_off))
            smaller = np.argwhere(all_distances[i] < cut_off)
            cv = self.circular_variance(np.subtract([rnp_pos[s[0]] for s in smaller], in_pos[i]))
            cvs.append(cv)
        return cvs

    def graph_circular_variance(self, frame):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        cvs = self.integrase_cv(frame)
        #print(dists)
        bin_width = .1
        n_bins = int(1 / bin_width)
        hist, bin_edges = np.histogram(cvs, bins=n_bins, range=(0,1))
        print(np.sum(hist[:6]), np.sum(hist))

        bin_middles = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(hist))]

        ax1.bar(bin_middles, hist, width=0.8 * bin_width)
        plt.legend()
        plt.savefig("int_cv_"+ str(frame) +".png", bbox_inches="tight")
        plt.show()

    def write_integrase_cv(self, frame_number):

        dists = self.integrase_cv(frame_number)
        points = self.frames[self.frame_numbers.index(frame_number)][self.letters.index("Nt")]
        f = open("out/rnp_cv.xyz", "w")
        f.write(str(len(dists)) + "\n\n")
        for i in range(len(points)):
            f.write(self.xyz_string(i, points[i], round(dists[i], 2) * 100))
        f.close()

    def graph_cv_integrase(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        hexamers = [np.sum(np.array(self.integrase_cv(frame)) < .6) for frame in self.frame_numbers]
        #print(hexamers)
        time = np.multiply(self.frame_numbers, 10)


        ax1.plot(time, hexamers, label="Surface Integrase")
        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("cv_sin.png", bbox_inches="tight")
        plt.show()

    def get_box(self, frame):

        box = cp.deepcopy(self.frames[self.frame_numbers.index(frame)][self.letters.index("B")])[-1]
        box.extend([0,0,0])
        return np.array(box)

    def get_rnp_rg(self, frame):

        rnc_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("R")]
        rnc_pos = self.min_image_poly(rnc_pos, frame)
        mean = np.average(rnc_pos, axis=0)
        all_distances = distances.distance_array(np.array(rnc_pos), np.array(mean))#, box=self.get_box(frame))
        dists = np.sum(np.multiply(all_distances, all_distances))
        #print(all_distances)
        #dists = 0
        #for i in range(len(rnc_pos)):
        #    for j in range(i, len(rnc_pos)):
        #        dists += all_distances[i][j] * all_distances[i][j]
        #return np.sqrt(dists / (2 * len(rnc_pos) * len(rnc_pos)))
        return np.sqrt(dists / len(rnc_pos))


    def calc_interacting_integrase(self, frame):

        total = 0
        cut = 50
        in_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("N")]
        in_pos = [in_pos[i * 4] for i in range(len(in_pos)) if i * 4 < len(in_pos)]
        if len(self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]) > 0:
            in_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]
            #in_pos = [in_pos[i * 2] for i in range(len(in_pos)) if i * 2 < len(in_pos)]
        # print(len(in_pos))
        on_pos = (self.frames[self.frame_numbers.index(frame)][self.letters.index("H")] +
                  self.frames[self.frame_numbers.index(frame)][self.letters.index("C")] +
                  self.frames[self.frame_numbers.index(frame)][self.letters.index("P")])

        all_distances = distances.distance_array(np.array(in_pos), np.array(on_pos))
        # print(len(rnp_pos))

        in_on_dist = mda

        for ind, pos in enumerate(in_pos):
            if np.any(all_distances[ind] < cut):
                total += 1

        return total, len(in_pos)




    def calc_surface_integrase(self, frame):

        total = 0
        in_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("N")]
        in_pos = [in_pos[i*4] for i in range(len(in_pos)) if i * 4 < len(in_pos)]

        if len(self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]) > 0:
            in_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]
            in_pos = [in_pos[i * 2] for i in range(len(in_pos)) if i * 2 < len(in_pos)]
        #print(len(in_pos))
        rnp_pos = (self.frames[self.frame_numbers.index(frame)][self.letters.index("R")] +
                   self.frames[self.frame_numbers.index(frame)][self.letters.index("S")])
        #print(len(rnp_pos))
        both_pos = in_pos + rnp_pos
        mean = np.average(both_pos, axis=0)
        centered_in_pos = np.subtract(in_pos, mean)
        #centered_both_pos = np.subtract(both_pos, mean)
        for ind, pos in enumerate(centered_in_pos):
            dotted = np.dot(pos, pos)
            #print(ind, np.sqrt(dotted[ind]))
            #print(dotted.shape)
            #biggest = np.max(dotted)
            #for lk, dot in enumerate(dotted):
                #print(lk, dot/ np.abs(dot) * np.sqrt(np.abs(dot)))
            #quit()
            #if list(dotted).index(biggest) == ind:
            if np.sqrt(dotted) > 220 and len(self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]) == 0:
                total += 1
            elif np.sqrt(dotted) > 230 and len(self.frames[self.frame_numbers.index(frame)][self.letters.index("Nt")]) > 0:
                total += 1
        return total, len(in_pos)

    def min_image_poly(self, poly, frame):

        for ind, pos in enumerate(poly[:-1]):
            diff = np.subtract(pos, poly[ind + 1])
            #print(ind, diff)
            for i in range(3):
                if diff[i] > self.get_box(frame)[i]/2 or diff[i] < -self.get_box(frame)[i]/2:
                    #print(ind, pos, poly[ind+1])
                    poly[ind + 1][i] += self.get_box(frame)[i] * diff[i] / np.abs(diff[i])
                    #print(pos, poly[ind+1])
        #print(poly)
        #quit()
        return poly

    def get_number_ip6_bound_hexamers(self, frame, cut=15):

        hexamers = self.get_hexamers(frame)
        ip6_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("I")]
        hex_pos = [hexamer[0] for hexamer in hexamers]
        if len(hex_pos) == 0 or len(ip6_pos)==0:
            return 0
        all_distances = distances.distance_array(np.array(hex_pos), np.array(ip6_pos), box=self.get_box(frame))

        count = 0
        for i in range(len(hexamers)):
            if np.any(all_distances[i] < cut):
                count += 1
        return count

    def get_number_ip6_bound_pentamers(self, frame, cut=15):

        pentamers = self.get_pentamers(frame)
        ip6_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("I")]
        pent_pos = [pentamer[0] for pentamer in pentamers]
        if len(pent_pos) == 0 or len(ip6_pos) == 0:
            return 0
        all_distances = distances.distance_array(np.array(pent_pos), np.array(ip6_pos), box=self.get_box(frame))
        count = 0
        for i in range(len(pentamers)):
            if np.any(all_distances[i] < cut):
                count += 1
        return count

    def get_number_ip6_bound_edge_middle_pentamers(self, frame, cut=15):

        edges, middles = self.get_edge_middle_pentamers(frame)

        ip6_pos = self.frames[self.frame_numbers.index(frame)][self.letters.index("I")]

        pent_pos = [pentamer[0] for pentamer in edges]
        count_edge = 0
        if len(pent_pos) == 0 or len(ip6_pos) == 0:
            count_edge = 0
        else:
            all_distances = distances.distance_array(np.array(pent_pos), np.array(ip6_pos), box=self.get_box(frame))

            for i in range(len(edges)):
                if np.any(all_distances[i] < cut):
                    count_edge += 1

        count_middle = 0
        pent_pos = [pentamer[0] for pentamer in middles]
        if len(pent_pos) == 0 or len(ip6_pos) == 0:
            count_middle = 0
        else:
            all_distances = distances.distance_array(np.array(pent_pos), np.array(ip6_pos), box=self.get_box(frame))
            for i in range(len(middles)):
                if np.any(all_distances[i] < cut):
                    count_middle += 1

        return count_edge, count_middle


    def graph_growth(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        hexamers = [self.get_number_hexamers(frame) for frame in self.frame_numbers]
        pentamers = [self.get_number_pentamers(frame) for frame in self.frame_numbers]
        time = np.multiply(self.frame_numbers, 10)



        ax1.plot(time, hexamers, label="Hexamers")
        ax1.plot(time, pentamers, label="Pentamers")
        ax1.set_ylim(0, 160)
        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("growth.png", bbox_inches="tight")
        plt.show()

        index = 0

        while hexamers[index] < 50:
            index += 1

        print("time to 50 hexamers:", 10 * (index + 1))

    def graph_surface_integrase(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        hexamers = [self.calc_surface_integrase(frame)[0] for frame in self.frame_numbers]
        time = np.multiply(self.frame_numbers, 10)

        ax1.plot(time, hexamers, label="Surface Integrase")
        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("sin.png", bbox_inches="tight")
        plt.show()

    def graph_interacting_integrase(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        hexamers = [self.calc_interacting_integrase(frame)[0] for frame in self.frame_numbers]
        time = np.multiply(self.frame_numbers, 10)

        ax1.plot(time, hexamers, label="Interacting Integrase")
        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("inin.png", bbox_inches="tight")
        plt.show()

    def graph_pentamers(self, count_edges=True, hex_x=False):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        time = np.multiply(self.frame_numbers, 10)
        if hex_x:
            time = [self.get_number_hexamers(frame) for frame in self.frame_numbers]
        ipentamers = [self.get_number_ip6_bound_pentamers(frame) for frame in self.frame_numbers]
        pentamers = [self.get_number_pentamers(frame) for frame in self.frame_numbers]
        ax1.plot(time, pentamers, label="Pentamers")
        ax1.plot(time, ipentamers, label="IP6 Bound Pentamers")
        if count_edges:
            #ipentamers = [self.get_number_ip6_bound_edge_middle_pentamers(frame) for frame in self.frame_numbers]
            pentamers = [self.get_number_edge_middle_pentamers(frame) for frame in self.frame_numbers]
            #edge_ipentamers = [ip[0] for ip in ipentamers]
            #middle_ipentamers = [ip[1] for ip in ipentamers]
            edge_pentamers = [p[0] for p in pentamers]
            middle_pentamers = [p[1] for p in pentamers]
            ax1.plot(time, edge_pentamers, label="Edge Pentamers")
            ax1.plot(time, middle_pentamers, label="Middle Pentamers")
            #ax1.plot(time, edge_ipentamers, label="IP6 Bound Edge Pentamers")
            #ax1.plot(time, middle_ipentamers, label="IP6 Bound Middle Pentamers")


        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Number")
        if hex_x:
            ax1.set_xlabel("Number of Hexamers")
        plt.legend()
        fname = "pent_growth.png"
        if hex_x:
            fname = "pent_growth_hex.png"
        plt.savefig(fname, bbox_inches="tight")
        plt.show()

    def graph_rg(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        rgs = [self.get_rnp_rg(frame) for frame in self.frame_numbers]
        time = np.multiply(self.frame_numbers, 10)

        ax1.plot(time, rgs)
        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("$R_{g}$ ($\AA$)")
        plt.legend()
        plt.savefig("rg_growth.png", bbox_inches="tight")
        plt.show()

    def graph_hexamers(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        ihexamers = [self.get_number_ip6_bound_hexamers(frame) for frame in self.frame_numbers]
        hexamers = [self.get_number_hexamers(frame) for frame in self.frame_numbers]

        time = np.multiply(self.frame_numbers, 10)

        ax1.plot(time, hexamers, label="Hexamers")
        ax1.plot(time, ihexamers, label="IP6 Bound Hexamers")
        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("hex_growth.png", bbox_inches="tight")
        plt.show()

    def graph_bound_growth(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        hexamers = [self.get_number_ip6_bound_hexamers(frame) for frame in self.frame_numbers]
        pentamers = [self.get_number_ip6_bound_pentamers(frame) for frame in self.frame_numbers]
        time = np.multiply(self.frame_numbers, 10)

        ax1.plot(time, hexamers, label="IP6 Bound Hexamers")
        ax1.plot(time, pentamers, label="IP6 Bound Pentamers")
        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Units")
        plt.legend()
        plt.savefig("bound_growth.png", bbox_inches="tight")
        plt.show()

    def max_capsid_size(self, frame):

        penta = self.get_pentamers(frame)
        hexa = self.get_hexamers()
        penta_pos = [pent[0] for pent in penta]
        hexa_pos = [h[0] for h in hexa]
        total = penta_pos + hexa_pos
        all_distances = distances.distance_array(np.array(total), np.array(total), box=self.get_box(frame))

        return np.max(all_distances)


class ReadSet(object):

    def __init__(self, frames_list, path_list):

        self.ace = []

        for i in range(len(frames_list)):
            #print(i, frames_list[i], path_list[i])
            self.ace.append(ReadOut(frames_list[i], path_list[i]))

    def graph_rg(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for a in self.ace:
            rgs = [a.get_rnp_rg(frame) for frame in a.frame_numbers]
            time = np.multiply(a.frame_numbers, 10)
            ax1.plot(time, rgs)

        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("$R_{g}$ ($\AA$)")
        plt.legend()
        plt.savefig("rg_growth.png", bbox_inches="tight", dpi=800)
        plt.show()

    def graph_pentamers(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        for a in self.ace:
            #ipentamers = [a.get_number_ip6_bound_pentamers(frame) for frame in a.frame_numbers]
            pentamers = [a.get_number_pentamers(frame) for frame in a.frame_numbers]
            time = np.multiply(a.frame_numbers, 10)
            ax1.plot(time, pentamers)#, label="Pentamers")
            #ax1.plot(time, ipentamers)#, label="IP6 Bound Pentamers")

        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Pentamers")
        plt.legend()
        plt.savefig("pent_growth.png", bbox_inches="tight", dpi=800)
        plt.show()

    def graph_growth(self):

        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        for a in self.ace:
            hexamers = [a.get_number_hexamers(frame) for frame in a.frame_numbers]
            #pentamers = [a.get_number_pentamers(frame) for frame in a.frame_numbers]
            time = np.multiply(a.frame_numbers, 10)
            ax1.plot(time, hexamers)#, label="Hexamers")
            #ax1.plot(time, pentamers)#, label="Pentamers")

        ax1.set_xlabel("Time ($10^6$ CG Timesteps)")
        ax1.set_ylabel("Hexamers")
        plt.legend()
        plt.savefig("growth.png", bbox_inches="tight", dpi=800)
        plt.show()
























