import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys

from itertools import chain
from networkx import ancestors,DiGraph,Graph
from networkx import connected_component_subgraphs,topological_sort
from pedigree import Pedigree
from random import sample


class Query:
    
    def __init__(self, qtl):
        """
        Constructor for Query class. Query class extracts related group of variables
        from segregation networks to be used in MMAP queries.

        Parameters
        ----------
        qtl : QTL object
            The QTL object that contains the segregation network from which the 
            query variables will be extracted.
        """
        self.qtl = qtl


    def extract_allele_type(self, allele):
        """
        Method that extracts all variables of the same allele type:

        Parameters
        ----------
        allele : int
            The index in the QTL object of the allele to extract.

        Returns
        -------
        query_nodes : {literal : int}
            Dictionary that maps node names to UAI node indexes.
        """
        G = self.qtl.allele_graph

        query_nodes = {}
        for node in G.nodes():
            if str(self.qtl[node,2]) == str(allele):
                query_nodes[node] = self.qtl.index(node)

        return query_nodes


    def extract_within_range(self, person, rnge, alleles=[0], max_offspring=2):
        """
        Method that extracts families within rgne generations of person.

        Parameters
        ----------
        person : literal
            Progeny ID around which to extract nodes.
        rnge : tuple or int
            If rnge is a tuple, variables within rnge[0] generations before and 
            rnge[1] generations after person will be extracted. If rnge is a 
            tuple, rnge is converted to the tuple (rnge,rnge).
        alleles : [int]
            List that contains ints that correspond to allele indexes. Nodes with
            allele IDs in alleles will be extracted.
        max_offspring : int
            The maximum number of offspring to extract from each extracted parent.

        Returns
        -------
        query_nodes : {literal : int}
            Dictionary that maps node names to UAI node indexes.
        """
        if type(rnge) is tuple:
            suc,prd = rnge
        elif type(rnge) is int:
            suc,prd = rnge,rnge
        else:
            raise TypeError('Query.extract_within_range: rnge must be a tuple or int')

        if str(person) not in self.qtl.ped.ped_graph:
            raise ValueError('Query.extract_within_range: person must be in QTL')

        frontier    = self._get_starting({str(person)}, suc)
        extracted   = self._extract_from_ped(suc, prd, frontier, max_offspring)
        query_nodes = self._extract_from_seg_graph(extracted, alleles)

        return query_nodes


    def extract_random_person(self, rnge, alleles=[0], max_offspring=2):
        """
        Method that extracts random family within range of a randomly chose
        person. See extract_within_range for parameter descriptions.
        """
        person = None
        while not person or self.qtl[person,0] == 'S':
            person = sample(self.qtl.allele_graph.nodes(), 1)[0]
        return self.extract_within_range(self.qtl[person,0], rnge, alleles, max_offspring)


    def _get_starting(self, frontier, suc):
        PG = self.qtl.ped.ped_graph

        for i in range(suc):
            next = set(chain.from_iterable([PG.predecessors(p) for p in frontier]))
            if not next:
                break
            frontier = next

        return frontier


    def _extract_from_ped(self, suc, prd, frontier, max_offspring):
        PG        = self.qtl.ped.ped_graph
        extracted = set()

        for i in range(suc + prd):
            while frontier:
                n   = frontier.pop()
                new = set([x for x in PG.successors(n) if x not in frontier and x not in extracted])
                new = set(sample(new, max_offspring)) if len(new) > max_offspring else new

                frontier  |= new
                extracted |= set(chain.from_iterable([PG.predecessors(n) for n in new]))
                extracted.add(n)

        return extracted


    def _extract_from_seg_graph(self, extracted, alleles):
        G       = self.qtl.allele_graph
        q_nodes = {}

        for node in extracted:
            for a in alleles:
                mom = self.qtl.mk_nd(node, 0, a)
                dad = self.qtl.mk_nd(node, 1, a)

                msg = None
                if mom not in G and dad not in G:
                    msg = 'Query._extract_from_seg_graph: {}/{} not in seg graph'.format(mom, dad)
                elif mom not in G:
                    msg = 'Query._extract_from_seg_graph: {} not in seg graph'.format(mom)
                elif dad not in G:
                    msg = 'Query._extract_from_seg_graph: {} not in seg graph'.format(dad)

                if msg:
                    raise ValueError(msg)

                q_nodes[mom] = self.qtl.index(mom)
                q_nodes[dad] = self.qtl.index(dad)

        return q_nodes
            

if __name__ == '__main__':

    pass
