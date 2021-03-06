# Using `main.py` 

## Command line Parameters

### Parameter: `ini_file`  

Required    : true
Type        : `str`
Flag        : `None`
Description : a path to a valid ini file used for generating UAI models, evidence,
		  	  and query files

## `ini` Files

Ini files should have the `.ini` extension.

### Sections

#### Input Section

Tag         : `[input]`
Required    : true
Description : this section contains paths to files that contain data necessary for
			  creating models

##### Input Parameters

###### Parameter: `pedigree`

Tag			: `pedigree`
Required    : true
Type        : `str`
Description : a path to a valid data (pedigree) file generated by QMSim, or any 
			  file in the same format

###### Parameter: `allele_freq`

Tag			: `allele_freq`
Required    : true
Type        : `str`
Description : a path to a valid allele frequency file generated by QMSim, or any 
			  file in the same format

###### Parameter: `genotype`

Tag			: `genotype`
Required    : true
Type        : `str`
Description : a path to a valid genotype file generated by QMSim, or any file in
		      the same format

###### Parameter: `qtl_pos`

Tag			: `qtl_pos`
Required    : true
Type        : `str`
Description : a path to a valid QTL position file generated by QMSim, or any file 
			  in the same format

#### Output Section

Tag         : `[output]`
Required    : false
Description : this section contains paths to files to which model data will be 
			  written

##### Output Parameters

###### Parameter: `model_output`

Tag			: `model_output`
Required    : false
Default     : `out.uai`
Type        : `str`
Description : the name/path of a file (doesn't need to exist yet) to which the
			  generated UAI model will be written

###### Parameter: `evidence_output`

Tag			: `evidence_output`
Required    : false
Default     : `out.uai.evid`
Type        : `str`
Description : the name/path of a file (doesn't need to exist yet) to which the
			  generated UAI evidence will be written

###### Parameter: `query_output`

Tag			: `query_output`
Required    : false
Default     : `out.uai.query` if `[query]` section is specified, `None` otherwise
Type        : `str`
Description : the name/path of a file (doesn't need to exist yet) to which query
			  file will be written

###### Parameter: `variable_output`

Tag			: `variable_output`
Required    : false
Default     : `None`
Type        : `str`
Description : the name/path of a file (doesn't need to exist yet) to which 
			  information about the variables in the generated model will
			  be written

#### Options Section

Tag         : `[options]`
Required    : false
Description : this section contains miscellaneous parameters used for building models

##### Options Parameters

###### Parameter: `num_alleles`

Tag			: `num_alleles`
Required    : false
Default     : 1
Type        : `int`
Description : the number of QTL to include in the model

###### Parameter: `prob`      

Tag			: `prob`
Required    : false
Default     : 0.5
Type        : `float`
Description : the probabily that a node will be observed at increasing depths in the
			  model's corresponding segregation network

#### Extract Section

Tag         : `[extract]`
Required    : false
Description : this section contains parameters regarding sub-pedigree extraction;
			  no sub-pedigree will be extracted if the `[extract]` section isn't
			  included in the `ini` file

##### Extract Parameters

###### `num_founders`

Tag			: `num_founders`
Required    : false
Default     : 5 if `[extract]` section is specified, `None` otherwise
Type        : `int`
Description : the number of founders to include in an extracted sub-pedigree

###### `num_generations`

Tag			: `num_generations`
Required    : false
Default     : 5 if `[extract]` section is specified, `None` otherwise
Type        : `int`
Description : the number of generations to include in an extracted sub-pedigree

###### `max_offspring`

Tag			: `max_offspring`
Required    : false
Default     : 2 if `[extract]` section is specified, `None` otherwise
Type        : `int`
Description : the number of offspring to include from each parent in an extracted 
			  sub-pedigree

###### `start_generation`

Tag			: `start_generation`
Required    : false
Default     : 0 if `[extract]` section is specified, `None` otherwise
Type        : `int`
Description : the first generation to include in extracted sub-pedigree

#### Query Section

Tag         : `[query]`
Required    : false
Description : this section contains parameters that specify how a query family
			  will be extracted from the generated model; no query will be extracted
			  if the `[query]` section isn't included in the `ini` file

##### Query Parameters

###### Parameter: `query_range`

Tag			: `query_range`
Required    : false
Default     : 2 if `[query]` section is specified, `None` otherwise
Type        : `int` or 2 `int`s separated by a space
Description : the number of generations that will be extracted "around" a query node;
			  for example if node 5 is chosen (randomly) as the query node, and 
			  `query_range` is 2, nodes from 2 generations before node 5 and 2
			  generations after node 5 will be extracted; if `query_range` is 3 1,
			  nodes from 3 generations before node 5, and 1 generation after node 5
			  will be extracted

###### Parameter: `alleles`      

Tag			: `alleles`
Required    : false
Default     : 0 if `[query]` section is specified, `None` otherwise
Type        : `int` or space separated `int`s
Description : the QTL indexes to extract from the model

###### Parameter: `max_offspring`

Tag			: `max_offspring`
Required    : false
Default     : 2 if `[query]` section is specified, `None` otherwise
Type        : `int` or `str`
Description : the max number of children from each parent to include in the
			  query family; can by `inf` (which acts as max)

# Using `uai.py` 

## Command Line Parameters

### Parameter: `qtl_file`  

Required    : true
Type        : `str`
Flag        : `-qtl` or `--qtl_file`
Description : a path to a valid genotype file generated by QMSim, or any file in
		      the same format

### Parameter: `allele_frequency_file`

Required    : true
Type        : `str`
Flag        : `-frq` or `--allele_frequency_file`
Description : a path to a valid allele frequency file generated by QMSim, or any 
			  file in the same format

### Parameter: `ped_file`

Required    : true
Type        : `str`
Flag        : `-ped` or `--ped_file`
Description : a path to a valid data (pedigree) file generated by QMSim, or any 
			  file in the same format

### Parameter: `loci_pos_file`

Required    : true
Type        : `str`
Flag        : `-pos` or `--loc_pos_file`
Description : a path to a valid QTL position file generated by QMSim, or any file 
			  in the same format

### Parameter: `network_output_file`

Required    : false
Default     : `out.uai`
Type        : `str`
Flag        : `-no` or `--network_output_file`
Description : the name/path of a file (doesn't need to exist yet) to which the
			  generated UAI model will be written

### Parameter: `evid_output_file`

Required    : false
Default     : `out.uai.evid`
Type        : `str`
Flag        : `-eo` or `--evid_output_file`
Description : the name/path of a file (doesn't need to exist yet) to which the
			  generated UAI evidence will be written

### Parameter: `num_alleles`

Required    : false
Default     : 1
Type        : `int`
Flag        : `-na` or `--num_alleles`
Description : the number of alleles to include in the model

### Parameter: `prob_to_observe`

Required    : false
Default     : 0.5
Type        : `float`
Flag        : `-po` or `--prob_to_observe`
Description : the probabily that a node will observed at increasing depths in the
			  model's corresponding segregation network

### Parameter: `var_assigns`

Required    : false
Default     : False
Type        : `boolean`
Flag        : `-va` or `--var_assigns`
Description : if true, a file named `var_assigns` will be created that contains
			  information about the variables in the generated model

### Parameter: `num_founders`

Required    : false
Default     : `None`
Type        : `int`
Flag        : `-nf` or `--num_founders`
Description : the number of founders to include in an extracted sub-pedigree
Note		: can only be used in conjunction with `num_gens`, `max_offspring`,
  			  and `start_gen`

### Parameter: `num_gens`

Required    : false
Default     : `None`
Type        : `int`
Flag        : `-ng` or `--num_gens`
Description : the number of generations to include in an extracted sub-pedigree
Note		: can only be used in conjunction with `num_founders`, `max_offspring`,
  			  and `start_gen`

### Parameter: `max_offspring`

Required    : false
Default     : `None`
Type        : `int`
Flag        : `-mo` or `--max_offspring`
Description : the number of offspring to include from each parent in an extracted 
			  sub-pedigree
Note		: can only be used in conjunction with `num_founders`, `num_gens`,
  			  and `start_gen`

### Parameter: `start_gen`

Required    : false
Default     : `None`
Type        : `int`
Flag        : `-sg` or `--start_gen`
Description : the first generation to include in extracted sub-pedigree
Note		: can only be used in conjunction with `num_founders`, `num_gens`,
  			  and `max_offspring`
