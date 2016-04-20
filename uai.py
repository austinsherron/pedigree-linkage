import networkx as nx

from linkage import Linkage
from networkx import adjacency_matrix,topological_sort


class UAI:

	def __init__(self, linkage, out='out.uai'):
		self.linkage = linkage
		self.out = out

	
	def do_write(self, out=None):
		out = out if out else self.out

		print('MARKOV')


if __name__ == '__main__':

	G = nx.DiGraph()
	l = Linkage(mrk_file='../data/tiny_qtl.txt', ped_file='../data/tiny_data.txt')
	print(topological_sort(G))
