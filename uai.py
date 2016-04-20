from linkage import Linkage


class UAI:

	def __init__(self, linkage, out='out.uai'):
		self.linkage = linkage
		self.out = out

	
	def do_write(self, out=None):
		out = out if out else self.out


if __name__ == '__main__':

	l = Linkage(mrk_file='../data/tiny_qtl.txt', ped_file='../data/tiny_data.txt')
