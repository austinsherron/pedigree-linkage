import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys

from networkx import DiGraph,Graph
from networkx import connected_component_subgraphs,topological_sort
from pedigree import Pedigree


class QTL:

    def __init__(self, ped, qtl_file, allele_file, num_alleles=30):
        self.ped = ped
        self.len = num_alleles

        self.build_allele_assignments(qtl_file, self.len)
        self.get_allele_info(allele_file, self.len)
        self.build_allele_graph(self.len)


    def show_bayes_net(self):
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


    def get_allele_info(self, allele_file, num_alleles=None):
        num_alleles = num_alleles if num_alleles else self.len
        allele_info = []
        for i,r in enumerate(self._file_iter(allele_file)):
            if i == num_alleles:
                break

            allele_info.append([a.split(':')[1] for a in r[4:]])

        self.allele_info = allele_info
        return allele_info


    def build_allele_assignments(self, qtl_file, num_alleles=None):
        num_alleles = num_alleles if num_alleles else self.len
        allele_assigns = {}
        for r in self._file_iter(qtl_file):
            id,l = r[0],r[1:(num_alleles * 2) + 1]
            allele_assigns[id] = [(l[i],l[i + 1]) for i in range(0, len(l) - 1, 2)]

        self.allele_assigns = allele_assigns
        return allele_assigns

        
    def build_allele_graph(self, num_alleles=None):
        num_alleles = num_alleles if num_alleles else self.len
        G = DiGraph()

        for i in range(num_alleles):
            for k,v in self.ped.pedigree.items():
                G.add_node(self.mk_nd(k, 1, i))
                G.add_node(self.mk_nd(k, 0, i))

            for k,v in self.ped.pedigree.items():
                for c in v['cs']:
                    type = 0 if v['sex'] == 'm' else 1 
                    G.add_edge(self.mk_nd(k, 1, i), self.mk_nd(c, type, i))
                    G.add_edge(self.mk_nd(k, 0, i), self.mk_nd(c, type, i))
                    G.add_edge(self.mk_nd('S', k, c, i), self.mk_nd(c, type, i))
                    
                    if i > 0:
                        G.add_edge(self.mk_nd('S', k, c, i - 1), self.mk_nd('S', k, c, i))

        self.vars = topological_sort(G)
        self.var_idxs = dict([(v,i) for i,v in enumerate(self.vars)])
        self.allele_graph = G
        return G


    def mk_nd(self, *args):
        return '_'.join(map(str, args))


    def _file_iter(self, file):
        with open(file) as f:
            next(f)
            for line in f:
                l = line.split()
                yield l


    def __getitem__(self, var):
        if type(var) is tuple:
            split_var = var[0].split('_')
            var,i = var
            return split_var[i]
        else:
            return var.split('_')


    def sort_vars(self, vars):
        srtd = sorted(vars, key=lambda x: self.var_idxs[x])
        for i,var in enumerate(srtd):
            if self[var,0] == 'S':
                    srtd.pop(i)
                    srtd.insert(0, var)
        return srtd


    def index(self, var):
        return self.var_idxs[var]


    def card(self, var):
        if self[var,0] == 'S':
            return 2
        else:
            return len(self.allele_info[int(self[var,2])])


    def allele_val(self, var, index=False):
        id,prnt,depth = tuple(map(int, self[var]))
        val = self.allele_assigns[str(id)][depth][prnt]
        if index:
            return int(val) - 1
        else:
            return val


if __name__ == '__main__':

    qtl_file = sys.argv[1] if len(sys.argv) > 1 else '../data/test_qtl.txt'
    allele_file = sys.argv[2] if len(sys.argv) > 2 else '../data/test_alleles.txt'
    ped_file = sys.argv[3] if len(sys.argv) > 3 else '../data/test_data.txt'

    ped = Pedigree(ped_file)
    qtl = QTL(ped, qtl_file=qtl_file, allele_file=allele_file, num_alleles=1)
    qtl.print_bayes_factors()

