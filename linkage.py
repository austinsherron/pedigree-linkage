import numpy as np

from pedigree import Pedigree
from networkx import adjacency_matrix,topological_sort,DiGraph,Graph
from networkx import enumerate_all_cliques,find_cliques,graph_number_of_cliques


"""
NOTE: it's assumed that the variables are named using the following
convention: i_s/i_d for i in range(# of variables). This should 
probably be more general.
"""


class Linkage(Pedigree):

    def __init__(self, mrk_file=None, ped_file=None, *args):
        self.allele_types = set()

        if ped_file:
            Pedigree.__init__(self, ped_file)
        else:
            raise ValueError('Linkage.__init__: Linkage object requires a pedigree')

        if len(args) <= 0:
            args = [1,2]

        if mrk_file:
            self.build_allele_assignments(mrk_file, *args)
            self.build_allele_network()
            self.build_allele_graph()


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
        

    def get_alleles(self, r, *args):
        return tuple(r[i] for i in args)


    def mk_nd(x, w):
        return str(x) + '_' + str(w)
        
        
    def index(x): 
        return (int(x[:-2]) * 2 + (0 if x[-1] == 'd' else 1)) - 1

        
    def build_allele_graph(self, trace=False):
        G = DiGraph()
        mk_nd = Linkage.mk_nd

        for k,v in self.allele_network.items():
            G.add_node(mk_nd(k, 'd'))
            G.add_node(mk_nd(k, 's'))

        for k,v in self.allele_network.items():
            for c in v['cs']:
                if v['sex'] == 'm':
                    G.add_edge(mk_nd(k, 'd'), mk_nd(c, 's'))
                    G.add_edge(mk_nd(k, 's'), mk_nd(c, 's'))
                    if trace:
                        print(k + '_d', '->', c + '_s')
                        print(k + '_s', '->', c + '_s')
                else:
                    G.add_edge(mk_nd(k, 'd'), mk_nd(c, 'd'))
                    G.add_edge(mk_nd(k, 's'), mk_nd(c, 'd'))
                    if trace:
                        print(k + '_d', '->', c + '_d')
                        print(k + '_s', '->', c + '_d')

        self.allele_graph = G
        return G


    def triangulate(self):
        mk_nd = Linkage.mk_nd
        self.tri_graph = self.allele_graph.to_undirected()
        for i in self.pedigree:
            self.tri_graph.add_edge(mk_nd(i, 's'), mk_nd(i, 'd'))


    def write_UAI(self):
        """
        NOTE: it's assumed that valid cliques are of size 3
        """
        G = self.allele_graph

        print('MARKOV')
        print(len(self.allele_graph))
        print(' '.join([str(len(self.allele_types)) for v in self.allele_graph]))
        self._print_clique_assigns()
        self._print_factors()


    def _print_clique_assigns(self):
        uG = self.tri_graph
        index = Linkage.index
        cliques = 0
        out = []
        for clique in find_cliques(uG):
            if len(clique) == 3:
                vars = [str(index(c)) for c in clique]
                vars.sort(key=lambda x: int(x))
                out.append(((len(clique), ' '.join(vars))))
                cliques += 1

        print(cliques)
        for clique_len,clique in out:
            print(clique_len, clique)


    def _print_factors(self):
        uG = self.tri_graph
        for clique in find_cliques(uG):
            cl = len(clique) 
            if cl == 3:
                print(len(self.allele_types) ** cl)
                self._build_factor(cl)


    def _build_factor(self, cl):
        factor = np.zeros((len(self.allele_types),) * cl)
        at = len(self.allele_types)

        for i in range(at):
            for j in range(at):
                factor[i][j][i] += 0.5
                factor[i][j][j] += 0.5

        for i in range(at):
            for j in range(at):
                for k in range(at):
                    print(factor[i][j][k], end=' ')

        print()


if __name__ == '__main__':

    l = Linkage(mrk_file='../data/test_qtl.txt', ped_file='../data/test_data.txt')
    l.triangulate()
    l.write_UAI()
