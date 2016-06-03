import matplotlib.pyplot as plt
import networkx as nx

from collections import defaultdict
from random import sample


class Pedigree:

    def __init__(self, file=None, draw=False):
        """
        Constructor for pedigree class.  Constructor parses pedigree file into
        the format that is used by the QTL class. 
        
        Parameters
        ----------
        file : str
            Name of file that contains pedigree. Pedigree files have rows 
            with the following format:

                ID (int > 0) Sire (int >= 0) Dam (int >= 0) Sex (M or F)

            Progeny with Sires/Dams labeled as 0 are founders.
        draw : boolean (optional) 
            If true, pedigree graph will be drawn (require networkx).
        """
        if file:
            self.build_pedigree(file)

            if draw:
                self.graph_pedigree(self.pedigree)


    def print_pedigree(self, pedigree=None):
        """
        Method that prints the pedigree in the following format:
            parent -> child_1 + child_2 ... + child_n

        Parameters
        ----------
        pedigree : {literal -> [literals]} (optional)
            Pedigree dict.  If none is supplied, self.pedigree will be used.
        """
        if not pedigree:
            pedigree = self.pedigree

        for p,cs in pedigree.items():
            cs = cs['cs']
            frm = cs[0] if len(cs) > 0 else ''
            for c in cs[1:]:
                frm += ' + ' + str(c)
            print(p, '->', frm)


    def print_ped_graph(self, ped_graph=None):
        if not ped_graph:
            ped_graph = self.ped_graph

        for node in ped_graph.nodes():
            out = node + ' -> ' + ' + '.join(ped_graph.successors(node))
            print(out)

    
    def build_pedigree(self, file):
        """
        Parameters
        ----------
        file : string
            String file name that contains pedigree data. It is expected that the
            file conains pedigree data, one individual per line. Refer to 
            constructor for format details.

        Returns
        -------
        pedigree : {literal: {'cs': [literal] 'sex': 'm'/'f'}}
            Dictionary rep (adj. list) of pedigree graph.
        """
        pedigree = {}
        inverted_ped = defaultdict(dict)
        ped_graph = nx.DiGraph()

        for r in self._file_iter(file):
            # build 'standard pedigree' (parent -> children)
            # also build ped graph for convenience
            pedigree[r[0]] = {'cs': [], 'sex': r[3].lower(), 'gen': int(r[4])}

            if int(r[1]) != 0:
                pedigree[r[1]]['cs'].append(r[0])
                ped_graph.add_edge(r[1], r[0])

            if int(r[2]) != 0:
                pedigree[r[2]]['cs'].append(r[0])
                ped_graph.add_edge(r[2], r[0])

            # also build 'inverted pedigree' (child -> parents)
            inverted_ped[r[0]] = {'s': r[1],'d': r[2]}

        self.pedigree = pedigree
        self.inverted_ped = inverted_ped
        self.ped_graph = ped_graph

        return pedigree


    def extract_sub_ped(self, num_founders=10, gens=10, max_offspring=2, start_gen=0):
        """
        """
        nodes = self.ped_graph.nodes()
        current_gen = [n for n in nodes if self.pedigree[n]['gen'] == start_gen]
        current_gen = set(sample(current_gen, num_founders))

        sub_ped = {}
        for i in range(gens):
            next_gen = set()
            for p in current_gen:
                sub_ped[p] = self.pedigree[p]
                if max_offspring > len(sub_ped[p]['cs']):
                    children = set(sub_ped[p]['cs'])
                else:
                    children = set(sample(sub_ped[p]['cs'], max_offspring))
                next_gen  |= children
                for child in children:
                    parents = self.ped_graph.predecessors(child) 
                    for parent in parents:
                        if parent not in sub_ped:
                            sub_ped[parent] = self.pedigree[parent]

            current_gen = next_gen

        self.pedigree = sub_ped
        return sub_ped


    def graph_pedigree(self, pedigree=None):
        """
        Method that graphs the given pedigree.

        Parameters
        ----------
        pedigree : {literal -> [literal]} (optional)
            Dictionary rep (adj. list) of pedigree graph. If none is provided,
            self.pedigree will be used.
        """
        if not pedigree:
            pedigree = self.pedigree

        G = nx.DiGraph()

        for p,cs in pedigree.items():
            G.add_node(p)

            for c in cs['cs']:
                G.add_edge(p, c)

        # change drawing method to change the shape of the graph
        pos = nx.draw_spectral(G)
        nx.draw(G, pos)
        plt.show()


    def _file_iter(self, file):
        """
        Custom iterator arbitrary files.  Method yields split lines of file,
        one at a time.

        Parameters
        ----------
        file : str
            Name of file that contains data (a pedigree, in this class).

        Yields 
        ------ 
        l : [str]
            Lines of file split on whitespace.
        """
        with open(file) as f:
            next(f)
            for line in f:
                l = line.split()
                yield l


if __name__ == '__main__':

    p = Pedigree('../data/test_data.txt')
    p.print_pedigree()
    #print()
    #print(p.pedigree)
    #p.print_ped_graph()
    p.extract_sub_ped(num_founders=3, gens=1, max_offspring=1)
    print()
    p.print_pedigree()
