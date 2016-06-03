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
        self.qtl = qtl


    def extract_allele_type(self, allele):
        G = self.qtl.allele_graph

        query_nodes = {}
        for node in G.nodes():
            if str(self.qtl[node,2]) == str(allele):
                query_nodes[node] = self.qtl.index(node)

        return query_nodes


    def extract_within_range(self, person, rnge, alleles=[0], max_offspring=2):
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


    def extract_random_family(self, alleles=None, max_offspring=None):
        pass


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

                if mom not in G or dad not in G:
                    raise ValueError('Query._extract_from_seg_graph: node not in seg graph')

                q_nodes[mom] = self.qtl.index(mom)
                q_nodes[dad] = self.qtl.index(dad)

        return q_nodes
            

if __name__ == '__main__':

    pass
