import numpy as np
import sys

from linkage import Linkage
from networkx import dfs_tree,find_cliques,topological_sort
from random import random as rand


class UAI(Linkage):

    def __init__(self, mrk_file=None, ped_file=None):
        Linkage.__init__(self, mrk_file, ped_file)


    def write(self):
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
        cliques = 0
        out = []
        for clique in find_cliques(uG):
            if len(clique) == 3:
                vars = [str(self.index(c)) for c in clique]
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


    def observe(self, prob=0.5):
        observed = {}
        nodes = topological_sort(self.allele_graph)
        num_nodes = len(nodes)
        for i,node in enumerate(nodes, 1):
            if rand() < prob * (i / num_nodes):
                if self.allele_network[node[:-2]][node[-1]] == None:
                    continue

                if self.allele_graph.in_degree(node) == 0 and rand() < 0.05:
                    val = self.allele_network[node[:-2]][node[-1]]
                    observed[self.index(node)] = val
                elif self.allele_graph.in_degree(node) > 0:
                    val = self.allele_network[node[:-2]][node[-1]]
                    observed[self.index(node)] = val

        print(len(observed), end=' ')
        for var,val in observed.items():
            print(var, val, end=' ')
        print()

        
    def _dfs(n, G, depth, visited=set()):
        visited.add(n)
        print('node', n, 'at depth', depth)
        for succ in G.successors(n):
            if succ not in visited:
                UAI._dfs(succ, G, depth + 1, visited)


if __name__ == '__main__':

    action   = sys.argv[1] if len(sys.argv) > 1 else 'gen'
    mrk_file = sys.argv[2] if len(sys.argv) > 2 else '../data/test_qtl.txt'
    ped_file = sys.argv[3] if len(sys.argv) > 3 else '../data/test_data.txt'

    uai = UAI(mrk_file=mrk_file, ped_file=ped_file)

    if action == 'gen':
        uai.write()
    elif action == 'observe':
        uai.observe(0.8)
