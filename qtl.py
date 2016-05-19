import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys

from networkx import DiGraph,Graph
from networkx import connected_component_subgraphs,topological_sort
from pedigree import Pedigree


class QTL:

    def __init__(self, ped, qtl_file, allele_file, num_alleles=30):
        """
        Constructor for QTL class. The QTL class relates QTL (quantitative
        trait loci) data to pedigree data: this class relates allele assigments
        to the individuals specified in a pedigree. The QTL class ultimately
        builds a segregation network from the pedigree in 'ped'.

        Parameters
        ----------
        ped : Pedigree object
            A Pedigree object that contains a dict representation  of a 
            pedigree in ped.pedigree.
        qtl_file : str
            The name of a file that relates individuals to their allele 
            assignments. QTL files have rows with the following formats:

                ID (int > 0) allele_0_a (int) allele_0_b ... allele_n_a allele_n_b

            If there are 30 alleles, the list produced by splitting one line 
            of a qtl file should be of lenght 61: 1 ID, and 2 * 30 allele
            assignments.
        allele_file : str
            The name of a file that specifies allele distributions in the
            first generation of the provided pedigree. This file is necessary
            for building the Bayes factors for founder alleles. Allele files
            have rows with the following formats:

               ID (str) Gen (0) Chr (int) Var (float) allele_n:freq (int:float) 

            Gen refers to the generation number (only 0 gen allele frequencies
            are needed). Chr refers to the chromosome the allele is found on.
            For each row, there should be a probability for each potential
            allele assignment.
        num_alleles : int (optional)
            The number of alleles to include. Should not be larger than the
            number of alleles available in the specified file.
        """
        self.ped = ped
        self.len = num_alleles

        self.build_allele_assignments(qtl_file, self.len)
        self.get_allele_info(allele_file, self.len)
        self.build_allele_graph(self.len)


    def show_bayes_net(self):
        """
        Method that draw the segregation network (require networkx).
        """
        G = self.allele_graph
        nx.draw(G)
        plt.show()


    def print_bayes_factors(self):
        """
        Method that draws bayes factors of the segregation network.
        """
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
        """
        Method that parses allele files and stores allele assignments
        distribution info. See constructor for parameter descriptions.

        Returns
        -------
        allele_info : [[str representation of floats]]
            List of lists that contain allele assignment statistics.
            allele_info[n][m] refers to the mth allele assignment of the
            nth locus.
        """
        num_alleles = num_alleles if num_alleles else self.len
        allele_info = []
        for i,r in enumerate(self._file_iter(allele_file)):
            # loop will stop if num_alleles is reached or if the number of 
            # alleles in allele_file < num_alleles
            if i == num_alleles:
                break

            # allele stats start at index 4
            allele_info.append([a.split(':')[1] for a in r[4:]])

        self.allele_info = allele_info
        return allele_info


    def build_allele_assignments(self, qtl_file, num_alleles=None):
        """
        Method that builds assignment of allele assignments to individuals
        in the provided pedigree. See constructor for parameter descriptions.

        Returns
        -------
        allele_assigns : {literal: [(literal,literal)]}
            Dictionary that maps IDs of individuals in pedigree to their allele
            assignments. Loci are represented as 2-tuples in a list.
            allele_assigns['1'][0][0] is the paternal allele at the 0th loci 
            for individual '1'. allele_assigns['1'][0][1] is the maternal
            allele at the same loci for the same individual.
            
        """
        num_alleles = num_alleles if num_alleles else self.len
        allele_assigns = {}
        for r in self._file_iter(qtl_file):
            # grab id and alleles
            id,l = r[0],r[1:(num_alleles * 2) + 1]
            # build allele list with comprehension
            allele_assigns[id] = [(l[i],l[i + 1]) for i in range(0, len(l) - 1, 2)]

        self.allele_assigns = allele_assigns
        return allele_assigns

        
    def build_allele_graph(self, num_alleles=None):
        """
        Method that builds segregation graph from the given pedigree. See the
        constructor for descriptions of parameters.

        Returns
        -------
        G : networkx.DiGraph
            NetworkX DiGraph representation of the segregation network defined
            by the given pedigree and allele assignments.  Nodes of the graph
            are named using the following conventions:
                allele node = id_parent_locus
                    ex.: 1_0_0 -> paternal allele of locus 0 for ID = 1
                    ex.: 2_1_3 -> maternal allele of locus 3 for ID = 2
                segregation node = S_parent_child_locus
                    ex.: S_1_6_1 -> segrgation indicator of parent ID = 1 
                    to child ID = 6 for locus 1
        """
        num_alleles = num_alleles if num_alleles else self.len
        G = DiGraph()

        # for each allele
        for i in range(num_alleles):
            # for each k = parent, v = [children, sex]
            for k,v in self.ped.pedigree.items():
                # make parent nodes, one for each allele, at locus i
                G.add_node(self.mk_nd(k, 1, i))
                G.add_node(self.mk_nd(k, 0, i))

            # for each k = parent, v = [children, sex]
            for k,v in self.ped.pedigree.items():
                # for each child of k
                for c in v['cs']:
                    # determine if parent is male or female
                    type = 0 if v['sex'] == 'm' else 1 
                    # add edges from parent alleles to children alleles
                    G.add_edge(self.mk_nd(k, 1, i), self.mk_nd(c, type, i))
                    G.add_edge(self.mk_nd(k, 0, i), self.mk_nd(c, type, i))
                    # add segregation node edge
                    G.add_edge(self.mk_nd('S', k, c, i), self.mk_nd(c, type, i))
                    
                    # if applicable, add edge between segregation indicators
                    if i > 0:
                        G.add_edge(self.mk_nd('S', k, c, i - 1), self.mk_nd('S', k, c, i))

        # for UAI: vars are indexed by topological sort
        self.vars = topological_sort(G)
        self.var_idxs = dict([(v,i) for i,v in enumerate(self.vars)])
        self.allele_graph = G
        return G


    def mk_nd(self, *args):
        """
        Helper method that makes graph node name from its arguments. This 
        method is really just a join method that uses the '_' (underscore)
        character as glue.

        Parameters
        ----------
        args : mixed
            Anything that can be converted to a str with the str function.

        Returns
        -------
        str
            str representation of args concatenated with the '_' character.
            ex.: args = [1,2,3] -> output = '1_2_3'
        """
        return '_'.join(map(str, args))


    def _file_iter(self, file):
        """
        See Pedigree._file_iter.
        """
        with open(file) as f:
            next(f)
            for line in f:
                l = line.split()
                yield l


    def __getitem__(self, var):
        """
        Method that overloads QTL indexing (QTL[...]). This method is generally
        used to extract information from names of nodes in the segregation
        network.
        
        Parameters
        ----------
        var : (str,) or (str,int)
            A tuple that contains either a single str or a str and an int.

        Returns
        -------
        str
            If var contains a single str, the list containing the substrings 
            of that str split on '_' will be returned; ex.: var = ('1_0_0',) ->
            ['1', '0', '0']. If var contains a str and an index, the index i of
            the list in the last example will be returned; ex.: var = 
            ('1_0_0',0) -> '1'.
        """
        if type(var) is tuple:
            split_var = var[0].split('_')
            var,i = var
            return split_var[i]
        else:
            return var.split('_')


    def sort_vars(self, vars):
        """
        This method sorts node names according to the following rules: all 'S' 
        (segregation nodes) come first, all other variable are sorted according
        to their positions in self.var_idxs (a topological sort of the 
        segregation graph). This method is used for sorting the factors of the 
        segregation network for UAI generation. The segregation nodes are made
        to come first purely for convenience sake.

        Parameters
        ----------
        vars : [str]
            List of node names to be sorted.

        Returns
        -------
        srted : [str]
            Sorted list of node names.
        """
        srtd = sorted(vars, key=lambda x: self.var_idxs[x])
        for i,var in enumerate(srtd):
            if self[var,0] == 'S':
                srtd.pop(i)
                srtd.insert(0, var)
        return srtd


    def index(self, var):
        """
        Method that returns the index of var in the topological sort performed
        during the construction of the QTL object.

        Parameters
        ----------
        var : str
            Name of variable/node.

        Returns
        -------
        int
            Index of var in topological ordering of variables.
        """
        return self.var_idxs[var]


    def card(self, var):
        """
        Method that returns the cardinality of var.

        Parameters
        ----------
        var : str
            str name of variable in segregation network.

        Returns
        -------
        int
            Cardinality of var.
        """
        # segregation nodes are binary
        if self[var,0] == 'S':
            return 2
        # the number of possible allele assignments for alleles at the locus
        # of the allele that var corresponds to
        else:
            return len(self.allele_info[int(self[var,2])])


    def allele_val(self, var, index=False):
        """
        Method that returns the actual value of the allele that corresponds to
        var.

        Parameters
        ----------
        var : str
            str name of variable in segregation network.
        index : boolean
            Returns the actual value if false, otherwise the index of the value
            (val - 1). TODO: this might have to be changed at some point, 
            depending on how allele values are represented.

        Returns
        -------
        val : literal
            str that is the value of the allele var corresponds to, or
            int(return val) - 1 (if index=True).
        """
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
#qtl.print_bayes_factors()
    print(qtl.allele_info)

