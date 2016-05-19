import argparse as ap
import numpy as np
import sys

from pedigree import Pedigree
from qtl import QTL
from networkx import dfs_tree,find_cliques,topological_sort
from random import random as rand


class UAI:

    def __init__(self, qtl):
        """
        Constructor for the UAI class. This is the class 
        """
        self.qtl = qtl


    def write(self):
        G = self.qtl.allele_graph

        print('MARKOV')
        print(len(G))
        self._print_cardinalities()
        print(len(G))
        self._print_cliques(G)
        self._print_factors(G)


    def _print_cardinalities(self):
        vars = self.qtl.vars
        for var in vars:
            print(self.qtl.card(var), end=' ')
        print()


    def _print_cliques(self, G):
        vars = self.qtl.vars
        for i,var in enumerate(vars):
            clique = self.qtl.sort_vars([var] + G.predecessors(var))
            print(len(clique), end=' ')
            print(' '.join([str(self.qtl.index(v)) for v in clique]))


    def _print_factors(self, G):
        vars = self.qtl.vars

        for i,var in enumerate(vars):
            clique = self.qtl.sort_vars([var] + G.predecessors(var))

            if len(clique) == 1 and self.qtl[var,0] != 'S':
                self._print_founder(clique[0])
            elif len(clique) == 1 and self.qtl[var,0] == 'S':
                self._print_seg_node()
            elif len(clique) == 4:
                self._print_non_founder(clique)
            elif len(clique) == 2 and self.qtl[var,0] == 'S':
                self._print_seg_clique()
            else:
                raise AssertionError('UAI._print_factors: unrecognized clique type')


    def _print_founder(self, var):
        allele_idx = int(self.qtl[var,2])
        print(len(self.qtl.allele_info[allele_idx]))
        for prob in self.qtl.allele_info[allele_idx]:
            print(prob, end=' ')
        print()


    def _print_seg_node(self):
        print(2)
        print(0.5, 0.5)


    def _print_seg_clique(self):
        print(4)
        print(0.5, 0.5, 0.5, 0.5)


    def _print_non_founder(self, clique):
        al     = self.qtl.allele_info[int(self.qtl[clique[1],2])]
        al_len = len(al)
        father = 1 if int(self.qtl[clique[1],1]) == 0 else 2
        factor = self._build_factor(father, al_len)
        print(len(factor))
        print(' '.join(map(str, factor)))


    def _build_factor(self, father, al_len):
        factor = []
        for i in range(2):
            for j in range(al_len):
                for k in range(al_len):
                    for l in range(al_len):
                        if i == 0 and father == 1 and j == l:
                            factor.append(1.0)
                        elif i == 0 and father == 2 and k == l:
                            factor.append(1.0)
                        elif i == 1 and father == 2 and j == l:
                            factor.append(1.0)
                        elif i == 1 and father == 1 and k == l:
                            factor.append(1.0)
                        else:
                            factor.append(0.0)

        return factor


    def observe(self, prob=0.5, see=None):
        observed = self._build_observed(prob, see)

        print(len(observed), end=' ')
        for var,val in observed.items():
            print(var, val, end=' ')
        print()

    
    def _build_observed(self, prob=0.5, see=None):
        observed = {}
        if not see:
            see = lambda *args: rand() < args[0] * (args[1] / args[2])

        for i,var in enumerate(self.qtl.vars):
            if not see(prob, i, len(self.qtl.vars)):
                continue

            if self.qtl[var,0] == 'S':
                continue

            observed[self.qtl.index(var)] = self.qtl.allele_val(var, index=True)

        return observed


    def _observe_seg_var(self, var):
        prnt         = self.qtl[var,1]
        id           = self.qtl[var,2]
        sex          = self.qtl.ped.pedigree[prnt]['sex']
        idx          = 0 if sex == 'm' else 1
        child_allele = self.allele_assigns[id][idx]
        parent_alls  = self.allele_assigns[prnt]

        
if __name__ == '__main__':

    parser = ap.ArgumentParser()
    parser.add_argument('qtl_file')
    parser.add_argument('allele_file')
    parser.add_argument('ped_file')
    parser.add_argument('-no', '--network_output_file')
    parser.add_argument('-eo', '--evid_output_file')
    parser.add_argument('-na', '--num_alleles')
    parser.add_argument('-po', '--percent_to_observe')
    args = parser.parse_args()

    ped = Pedigree(args.ped_file)
    num = int(args.num_alleles) if args.num_alleles else 1
    qtl = QTL(ped, args.qtl_file, args.allele_file, num)
    uai = UAI(qtl)

    net_out_file = args.network_output_file if args.network_output_file else 'out.uai'
    with open(net_out_file, 'w') as f:
        sys.stdout = f
        uai.write()

    evid_out_file = args.evid_output_file if args.evid_output_file else 'out.uai.evid'
    with open(evid_out_file, 'w') as f:
        sys.stdout = f
        uai.observe()

    with open('var_assigns', 'w') as f:
        sys.stdout = f
        print(qtl.vars)
        print(qtl.allele_assigns)
