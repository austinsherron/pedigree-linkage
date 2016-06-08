#!/usr/bin/python

import argparse as ap
import numpy as np
import sys

from pedigree import Pedigree
from qtl import QTL
from query import Query
from networkx import dfs_tree,find_cliques,topological_sort
from random import random as rand


class UAI:

    def __init__(self, qtl):
        """
        Constructor for the UAI class. This is the class that generates UAI
        from a QTL class with a fully formed segregation network.

        Parameters
        ----------
        qtl : QTL object
            Fully constructed QTL object that contains a well formed 
            segregation graph. 
        """
        self.qtl = qtl


    def write(self):
        """
        Container method responsible for writing UAI model files.
        """
        G = self.qtl.allele_graph

        # header prefix (model type)
        print('MARKOV')
        # print number of variables
        print(len(G))
        # cardinalities
        self._print_cardinalities()
        # number of factors (== to number of variables)
        print(len(G))
        # variable -> clique assignments
        self._print_cliques(G)
        # factor tables
        self._print_factors(G)


    def _print_cardinalities(self):
        """
        Method that prints the cardinalities of variables in the segregation 
        network.
        """
        vars = self.qtl.vars
        for var in vars:
            print(self.qtl.card(var), end=' ')
        print()


    def _print_cliques(self, G):
        """
        Method that prints the variables clique assignments.

        Parameters
        ----------
        G : networkx DiGraph
            nx.DiGraph representation of a segregation network.
        """
        vars = self.qtl.vars
        for i,var in enumerate(vars):
            # vars in clique are sorted according to special ordering defined
            # in QTL.sort_vars method; special ordering involves segregation
            # nodes coming before allele nodes; this is to make printing factor
            # tables easier
            clique = self.qtl.sort_vars([var] + G.predecessors(var))
            print(len(clique), end=' ')
            print(' '.join([str(self.qtl.index(v)) for v in clique]))


    def _print_factors(self, G):
        """
        Method that prints factor tables of cliques.

        Parameters
        ----------
        G : networkx DiGraph
            nx.DiGraph representation of a segregation network.
        """
        vars = self.qtl.vars

        for i,var in enumerate(vars):
            clique = self.qtl.sort_vars([var] + G.predecessors(var))

            # founder allele 
            if len(clique) == 1 and self.qtl[var,0] != 'S':
                self._print_founder(clique[0])
            # founder segregation node
            elif len(clique) == 1 and self.qtl[var,0] == 'S':
                self._print_seg_node()
            # non-founder clique
            elif len(clique) == 4:
                self._print_non_founder(clique)
            # non-founder segregation clique
            elif len(clique) == 2 and self.qtl[var,0] == 'S':
                self._print_seg_clique(clique)
            else:
                raise AssertionError('UAI._print_factors: unrecognized clique type of size ' + str(len(clique)))


    def _print_founder(self, var):
        """
        Method that prints the probability table for a founder node (clique).

        Parameter
        ---------
        var : str
            String name of a founder node in a segregation network.         
        """
        # get index of allele
        allele_idx = int(self.qtl[var,2])
        print(len(self.qtl.allele_info[allele_idx]))
        # print empirical probabilities for allele assignments
        for prob in self.qtl.allele_info[allele_idx]:
            print(prob, end=' ')
        print()


    def _print_seg_node(self, prob=0.5):
        """
        Method that prints probabilities for founder segregation node.

        Parameters
        ----------
        prob : float (optional)
            Probability that sire allele segregates.
        """
        print(2)
        print(prob, 1 - prob)


    def _print_seg_clique(self, clique):
        """
        Method that prints probabilities for non-founder segregation clique.
        """
        print(4)
        allele_1 = self.qtl[clique[0],3]
        allele_2 = self.qtl[clique[1],3]
        recomb_frac = self.qtl.recomb_prob(allele_1, allele_2)
        for i in range(2):
            for j in range(2):
                if i == j:
                    print(1 - recomb_frac, end=' ')
                else:
                    print(recomb_frac, end=' ')
        print()


    def _print_non_founder(self, clique):
        """
        Method that prints probability table for non-founder clique.

        Parameters
        ----------
        clique : [str]
            List of string var names in segregation network.
        """
        # figure out which allele this clique describes
        al     = self.qtl.allele_info[int(self.qtl[clique[1],2])]
        al_len = len(al)
        # index of allele inherited from father in factor table depends on 
        # whether this clique describes a sire allele or a dam allele
        father = 1 if int(self.qtl[clique[1],1]) == 0 else 2
        factor = self._build_factor(father, al_len)
        print(len(factor))
        print(' '.join(map(str, factor)))


    def _build_factor(self, father, al_len):
        """
        Method that builds factors table for non-founder clique.

        Parameters
        ----------
        father : int
            Index of parent allele inherited from father.
        al_len : int
            Number of possible allele assignments.

        Returns
        -------
        factor : [scalar]
            Factor table for a non-founder clique.
        """
        factor = []
        # for each possible assignment of the segregation node
        #   i = 0 means allele inherited from father segregated
        #   i = 1 means same for allele inherited from dam
        for i in range(2):
            # for each possible assignment of one parent allele
            for j in range(al_len):
                # for each possible assignment of the other parent allele
                for k in range(al_len):
                    # for each possible assignment of the inherited allele
                    for l in range(al_len):
                        # determine probability of assignment...
                        # inherited sire allele
                        if i == 0 and father == 1 and j == l:
                            factor.append(1.0)
                        # inherited sire allele
                        elif i == 0 and father == 2 and k == l:
                            factor.append(1.0)
                        # inherited dam allele
                        elif i == 1 and father == 2 and j == l:
                            factor.append(1.0)
                        # inherited dam allele
                        elif i == 1 and father == 1 and k == l:
                            factor.append(1.0)
                        else:
                            factor.append(0.0)

        return factor


    def observe(self, prob=0.5, see=None):
        """
        Container method responsible for writing UAI evidence files.

        Parameters
        ----------
        prob : float (optional)
            Partial probability of observing a node.  Real probablity starts
            at zero and increases to value of prob as _build_observed moves
            through a topological sort of the seg graph's nodes. This founder
            nodes from being observed too often.
        see : function (optional)
            Function that determines real probability of observing nodes. Function
            takes three paramters:
                prob : float
                    Base probability of observing node.
                i : int
                    Index of variable being considered for observation in a 
                    topological sort. 
                total : int
                    Total number of nodes in seg graph.
            If no function is provided, the following is used:
                prob to be observed = prob * (i / total)
            Nodes "higher" in a toplogical sort are less likely to be observed.
        """
        observed = self._build_observed(prob, see)

        print(len(observed), end=' ')
        for var,val in observed.items():
            print(var, val, end=' ')
        print()

    
    def _build_observed(self, prob=0.5, see=None):
        """
        Method that build map of variables and their real values. See observe
        doc string for descriptions of parameters.

        Returns
        -------
        observed = {literal: literal}
            Dictionary that maps variables to their true assignments.
        """
        observed = {}
        if not see:
            see = lambda *args: rand() < args[0] * (args[1] / args[2])

        # for each variable in a topological sort of the seg graph
        for i,var in enumerate(self.qtl.vars):
            # determine whether or not to observe the variable
            if not see(prob, i, len(self.qtl.vars)):
                continue

            # TODO: can't currently observe segregation nodes
            if self.qtl[var,0] == 'S':
                continue

            # find true value of var and record it
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
    parser.add_argument('-qtl', '--qtl_file', type=str, required=True)
    parser.add_argument('-frq', '--allele_frequency_file', type=str, required=True)
    parser.add_argument('-ped', '--ped_file', type=str, required=True)
    parser.add_argument('-pos', '--loci_pos_file', type=str, required=True)
    parser.add_argument('-no', '--network_output_file', type=str)
    parser.add_argument('-eo', '--evid_output_file', type=str)
    parser.add_argument('-na', '--num_alleles', type=int)
    parser.add_argument('-po', '--prob_to_observe', type=float, default=0.5)
    parser.add_argument('-va', '--var_assigns', type=str)
    parser.add_argument('-nf', '--num_founders', type=int)
    parser.add_argument('-ng', '--num_gens', type=int)
    parser.add_argument('-mo', '--max_offspring', type=int)
    parser.add_argument('-sg', '--start_gen', type=int)
    args = parser.parse_args()

    ped = Pedigree(args.ped_file)
    if None not in [args.num_founders, args.num_gens, args.max_offspring, args.start_gen]:
        extract_args = {
            'num_founders': args.num_founders,
            'gens': args.num_gens,
            'max_offspring': args.max_offspring,
            'start_gen': args.start_gen
        }
        ped.extract_sub_ped(**extract_args)
    num = args.num_alleles if args.num_alleles else 1
    wva = args.var_assigns
    qtl = QTL(ped, args.qtl_file, args.allele_frequency_file, args.loci_pos_file, num)
    uai = UAI(qtl)

    net_out_file = args.network_output_file if args.network_output_file else 'out.uai'
    with open(net_out_file, 'w') as f:
        sys.stdout = f
        uai.write()

    evid_out_file = args.evid_output_file if args.evid_output_file else 'out.uai.evid'
    with open(evid_out_file, 'w') as f:
        sys.stdout = f
        uai.observe(prob=args.prob_to_observe)

    if wva:
        with open('var_assigns', 'w') as f:
            sys.stdout = f
            qtl.print_var_info()
