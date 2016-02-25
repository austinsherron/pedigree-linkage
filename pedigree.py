import matplotlib.pyplot as plt
import networkx as nx

from collections import defaultdict


class Pedigree:


	def __init__(self, file=None, draw=False):
		if file:
			self.build_pedigree(file)

			if draw:
				self.graph_pedigree(self.pedigree)


	def print_pedigree(self, pedigree=None):
		if not pedigree:
			pedigree = self.pedigree

		for p,cs in pedigree.items():
			frm = cs[0] if len(cs) > 0 else ''
			for c in cs[1:]:
				frm += ' + ' + str(c)
			print(p, '->', frm)

	
	def build_pedigree(self, file):
		"""
		Parameters
		----------
		file : string
			String file name that contains pedigree data. It is expected that the
			file conains pedigree data, one individual per line, with an id at index 0,
			and parents at indices 1 and 2.

		Returns
		-------
		pedigree : {literal -> [literal]}
			Dictionary rep (adj. list) of pedigree graph.
		"""
		pedigree = defaultdict(list)
		inverted_ped = defaultdict(dict)

		for r in self._file_iter(file):
			pedigree[r[0]] = []

			if int(r[1]) != 0:
				pedigree[r[1]].append(r[0])

			if int(r[2]) != 0:
				pedigree[r[2]].append(r[0])

			inverted_ped[r[0]] = {'s': r[1],'d': r[2]}

		self.pedigree = pedigree
		self.inverted_ped = inverted_ped

		return pedigree


	def graph_pedigree(self, pedigree=None):
		"""
		Parameters
		----------
		pedigree : {literal -> [literal]}
			Dictionary rep (adj. list) of pedigree graph.
		"""
		if not pedigree:
			pedigree = self.pedigree

		G = nx.DiGraph()

		for p,ans in pedigree.items():
			G.add_node(p)

			for an in ans:
				G.add_edge(an, p)

		pos = nx.draw_spectral(G)
		nx.draw(G, pos)
		plt.show()


	def _file_iter(self, file):
		"""
		Custom iterator for pedigree file.

		Parameters
		----------
		file : string
			String file name that contains pedigree data. It is expected that the
			file conains pedigree data, one individual per line, with an id at index 0,
			and parents at indices 1 and 2.
		"""
		with open(file) as f:
			next(f)
			for line in f:
				l = line.split()
				yield l


if __name__ == '__main__':

	p = Pedigree('../QMSim_Mac/r_ex01/p1_data_001.txt')
	print(p.print_pedigree())
