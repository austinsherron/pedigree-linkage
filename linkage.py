from pedigree import Pedigree


class Linkage(Pedigree):

	def __init__(self, mrk_file=None, ped_file=None, pedigree=None, *args):
		if pedigree:
			self.pedigree = pedigree
		elif ped_file:
			Pedigree.__init__(self, ped_file)
		else:
			raise ValueError('Linkage.__init__: Linkage object requires a pedigree')

		if len(args) <= 0:
			args = [1,2]

		if mrk_file:
			self.build_allele_assignments(mrk_file, *args)
			self.build_allele_network()


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
				'cs': cs
			}

		self.allele_network = allele_network
		return allele_network


	def build_allele_assignments(self, mrk_file, *args):
		allele_assigns = {}
		for r in self._file_iter(mrk_file):
			id = r[0]
			s_allele,d_allele = self.get_alleles(r, *args)
			allele_assigns[id] = {'s': s_allele, 'd': d_allele}

		self.allele_assigns = allele_assigns
		return allele_assigns
		

	def get_alleles(self, r, *args):
		return tuple(r[i] for i in args)


if __name__ == '__main__':

	p = Pedigree('../QMSim_Mac/r_ex01/p1_data_001.txt')
	l = Linkage(mrk_file='../QMSim_Mac/r_ex01/p1_qtl_001.txt', pedigree=p.pedigree)
	l.print_allele_network(sorting=sorted)
