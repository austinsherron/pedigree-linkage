import configparser as cp
import sys

#from pedigree import Pedigree
#from qtl import QTL
#from query import Query


if len(sys.argv) < 2:
		print('Error: this script require one command line argument:')
		print('   ini_file.ini: the path to a valid ini file')
		sys.exit()

## helpers #####################################################################

def get_params(parser, section):
	try:
		return(dict(parser.items(section)))
	except cp.NoSectionError:
		return None

def get_param(section, param, default, type=str):
	if not section or param not in section:
		return default
	else:
		return type(section[param])

## parse ini file ##############################################################

parser = cp.ConfigParser()
with open(sys.argv[1]) as f:
	parser.readfp(f)

## get params ##################################################################

input   = get_params(parser, 'input')
output  = get_params(parser, 'output')
extract = get_params(parser, 'extract')
query   = get_params(parser, 'query')
options = get_params(parser, 'options')

## input params ################################################################

if not input:
	print('Error: [input] section required in ini file')
	sys.exit()

if 'pedigree' not in input:
	print('Error: "pedigree" param required in [input] section of ini file')
	sys.exit()

ped_file = input['pedigree']

if 'allele_freq' not in input:
	print('Error: "allele_freq" param required in [input] section of ini file')
	sys.exit()

frq_file = input['allele_freq']

if 'genotype' not in input:
	print('Error: "genotype" param required in [input] section of ini file')
	sys.exit()

qtl_file = input['genotype']

if 'qtl_pos' not in input:
	print('Error: "qtl_pos" param required in [input] section of ini file')
	sys.exit()

pos_file = input['qtl_pos']

## output params ###############################################################

network_output_file = get_param(output, 'model_output', 'out.uai')
evid_output_file    = get_param(output, 'evidence_output', 'out.uai.evid')
query_output_file   = get_param(output, 'query_output', 'out.uai.query')
variable_file       = get_param(output, 'variable_output', None)

## extract params ##############################################################

if extract:
	num_founders  = get_param(extract, 'num_founders', 5, int)
	num_gens      = get_param(extract, 'num_generations', 5, int)
	max_offspring = get_param(extract, 'max_offspring', 2, int)
	start_gen     = get_param(extract, 'start_generation', 0, int)

## query params ################################################################

if query:
	if 'query_range' in query:
		query_range = list(map(int, query['query_range'].split(' ')))
		if len(query_range) > 2:
			print('Error: "query_range" [query] param can have at most 2 values')
			sys.exit()
		query_range = query_range[0] if len(query_range) == 1 else tuple(query_range)
	else:
		query_range = 2

	if 'alleles' in query:
		alleles = list(map(int, query['alleles'].split(' ')))
	else:
		alleles = [0]

	if 'max_offspring' in query:
		max_offspring_query = query['max_offspring']
		max_offspring_query = float('inf') if max_offspring_query == 'inf' else int(max_offspring_query)
	else:
		max_offspring_query = 2

## option params ###############################################################

num_alleles = get_param(options, 'num_alleles', 1, int)
prob        = get_param(options, 'prob', 0.5, float)

## main ########################################################################

#ped = Pedigree(ped_file)
#
#if extract:
#	ped.extract_sub_ped(num_founders, num_gens, max_offspring, start_gen)
#
#qtl = QTL(ped, qtl_file, frq_file, pos_file, num_alleles)
#uai = UAI(qtl)
#
#with open(network_output_file, 'w') as f:
#	sys.stdout = f
#	uai.write()
#
#with open(evid_output_file, 'w') as f:
#	sys.stdout = f
#	uai.observe(prob)
#
#if variable_file:
#	with open(variable_file, 'w') as f:
#		sys.stdout = f
#		qtl.print_var_info()
#
#if query:
#	qur = Query(qtl)
#	with open(query_output_file, 'w') as f:
#		sys.stdout = f
#		vars = qur.extract_random_person(query_range, alleles, max_offspring_query).items()
#		print(len(vars), end=' ')
#		for v,vi in vars:
#			print(vi, end=' ')
#		print()
