import numpy as np
import sys

from pedigree import Pedigree
from networkx import DiGraph,Graph


"""
NOTE: it's assumed that the variables are named using the following
convention: i_s/i_d for i in range(# of variables). This should 
probably be more general.
"""


class Linkage(Pedigree):

    def __init__(self, mrk_file=None, ped_file=None, *args):
        if ped_file and mrk_file:
            Pedigree.__init__(self, ped_file)
        else:
            raise ValueError('Linkage.__init__: Linkage object requires a pedigree')

        if len(args) <= 0:
            args = [1,2]

        self.allele_types = set()

        if mrk_file:
            self.build_allele_assignments(mrk_file, *args)
            self.build_allele_network()
            self.build_allele_graph()
            self.triangulate()

        self.allele_types = sorted(self.allele_types)


    def print_allele_network(self, allele_network=None, sorting=lambda x: x):
        allele_network = allele_network if allele_network else self.allele_network

        for p,d in sorting(allele_network.items()):
            print('{}_0: {}; {}_1: {} -> {}'.format(p, d['s'], p, d['d'], d['cs']))


    def build_allele_network(self, allele_assigns=None):
        allele_assigns = allele_assigns if allele_assigns else self.allele_assigns
        allele_network = {}

        for p,cs in self.pedigree.items():
            s = allele_assigns[p]['s'] if p in allele_assigns else None
            d = allele_assigns[p]['d'] if p in allele_assigns else None

            allele_network[p] = {
                's': s,
                'd': d,
                'cs': cs['cs'],
                'sex': cs['sex']
            }

        self.allele_network = allele_network
        return allele_network


    def build_allele_assignments(self, mrk_file, *args):
        allele_assigns = {}
        for r in self._file_iter(mrk_file):
            id = r[0]
            s_allele,d_allele = self.get_alleles(r, *args)
            self.allele_types |= set(self.get_alleles(r, *args))
            allele_assigns[id] = {'s': s_allele, 'd': d_allele}

        self.allele_assigns = allele_assigns
        return allele_assigns

        
    def build_allele_graph(self):
        G = DiGraph()

        for k,v in self.allele_network.items():
            G.add_node(self.mk_nd(k, 'd'))
            G.add_node(self.mk_nd(k, 's'))

        for k,v in self.allele_network.items():
            for c in v['cs']:
                if v['sex'] == 'm':
                    G.add_edge(self.mk_nd(k, 'd'), self.mk_nd(c, 's'))
                    G.add_edge(self.mk_nd(k, 's'), self.mk_nd(c, 's'))
                else:
                    G.add_edge(self.mk_nd(k, 'd'), self.mk_nd(c, 'd'))
                    G.add_edge(self.mk_nd(k, 's'), self.mk_nd(c, 'd'))

        self.allele_graph = G
        return G


    def triangulate(self):
        self.tri_graph = self.allele_graph.to_undirected()
        for i in self.pedigree:
            self.tri_graph.add_edge(self.mk_nd(i, 's'), self.mk_nd(i, 'd'))

        return self.tri_graph
        

    def get_alleles(self, r, *args):
        return tuple(r[i] for i in args)


    def mk_nd(self, x, w):
        return str(x) + '_' + str(w)
        
        
    def index(self, x): 
        return (int(x[:-2]) * 2 + (0 if x[-1] == 'd' else 1)) - 1


if __name__ == '__main__':

    mrk_file = sys.argv[1] if len(sys.argv) > 1 else '../data/test_qtl.txt'
    ped_file = sys.argv[2] if len(sys.argv) > 2 else '../data/test_data.txt'

    l = Linkage(mrk_file=mrk_file, ped_file=ped_file)
