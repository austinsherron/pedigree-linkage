import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys

from networkx import DiGraph,Graph
from networkx import topological_sort
from pedigree import Pedigree


"""
NOTE: it's assumed that the variables are named using the following
convention: i_s/i_d for i in range(# of variables). This should 
probably be more general.
"""


class QTL:

    def __init__(self, ped, qtl_file=None, allele_file=None, num_alleles=30):
        self.ped = ped

        if qtl_file:
            self.build_allele_assignments(qtl_file, num_alleles)
    
        if allele_file:
            self.get_allele_info(allele_file, num_alleles)
            self.build_allele_graph()


    def print_allele_network(self):
        G = self.allele_graph
        nx.draw(G)
        plt.show()


    def print_bayes_factors(self):
        Gs = topological_sort(self.allele_graph)
        G = self.allele_graph

        for n in Gs:
            p = 'P(' + n
            parents = G.predecessors(n)
            if parents:
                p += '|' + ','.join(parents)
            p += ')'

            print(p)


    def get_allele_info(self, allele_file, lim=30):
        allele_info = []
        for i,r in enumerate(self._file_iter(allele_file)):
            if i == lim:
                break

            allele_info.append([a.split(':')[1] for a in r[4:]])

        self.allele_info = allele_info
        return allele_info


    def build_allele_assignments(self, qtl_file, lim=30):
        allele_assigns = {}
        for r in self._file_iter(qtl_file):
            id,l = r[0],r[1:lim + 1]
            allele_assigns[id] = [(l[i],l[i + 1]) for i in range(0, len(l) - 1, 2)]

        self.allele_assigns = allele_assigns
        return allele_assigns

        
    def build_allele_graph(self):
        G = DiGraph()

        for k,v in self.ped.pedigree.items():
            G.add_node(self.mk_nd(k, 'd'))
            G.add_node(self.mk_nd(k, 's'))

        for k,v in self.ped.pedigree.items():
            for c in v['cs']:
                type = 's' if v['sex'] == 'm' else 'd'
                G.add_edge(self.mk_nd(k, 'd'), self.mk_nd(c, type))
                G.add_edge(self.mk_nd(k, 's'), self.mk_nd(c, type))
                G.add_edge(self.mk_nd('S', k, c), self.mk_nd(c, type))

        self.allele_graph = G
        return G


    def mk_nd(self, *args):
        return '_'.join(map(str, args))
        
        
    def index(self, x): 
        return (int(x[:-2]) * 2 + (0 if x[-1] == 'd' else 1)) - 1


    def _file_iter(self, file):
        with open(file) as f:
            next(f)
            for line in f:
                l = line.split()
                yield l


if __name__ == '__main__':

    qtl_file = sys.argv[1] if len(sys.argv) > 1 else '../data/test_qtl.txt'
    allele_file = sys.argv[2] if len(sys.argv) > 2 else '../data/test_alleles.txt'
    ped_file = sys.argv[3] if len(sys.argv) > 3 else '../data/test_data.txt'

    ped = Pedigree(ped_file)
    qtl = QTL(ped, qtl_file=qtl_file, allele_file=allele_file)
#   print(qtl.allele_assigns)
#    print(qtl.allele_info)
    qtl.print_bayes_factors()
    qtl.print_allele_network()

