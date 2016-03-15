#!/usr/bin/env python
'''USAGE: python python interact_scc.py --sccfiles1 organism_1.scc.files --sccfiles2 organism_2.scc.files --name1 organism_1 --name2 organism_2 --confidence 0
Calculates interaction metrics and writes a report based on output from SCC networks (calc_seeds.py)
'''

# Initialize all modules, functions, and compound dictionary
import sys
import math
import os
import pickle
import argparse

# Find the directory where interaction.py is being run from
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Create a string with the name of the initial directory you start running the script from
starting_directory = str(os.getcwd())

# Read in pickled reaction to compound name dictionary
comp_pckl_path = script_path + '/support/compound.pkl'
compound_dict = pickle.load(open(comp_pckl_path, 'rb'))

#---------------------------------------------------------------------------------------#

# Function for calculating competitive index between 2 lists of set seeds
def competition(seed1, seed2):
	
	len_seed1 = len(seed1)
	len_seed2 = len(seed2)
	
	intersection_seeds = list(set(seed1) & set(seed2))
	seed_overlap = len(intersection_seeds)
	
	seed1_only = list(set(seed1) - set(seed2))
	seed2_only = list(set(seed2) - set(seed1))

	# Calculate competition metrics
	comp_index1 = float(seed_overlap) / float(len_seed1)
	comp_index2 = float(seed_overlap) / float(len_seed2)
	comp_indices = [comp_index1, comp_index2]
	
	output_List = [intersection_seeds, seed1_only, seed2_only, comp_indices]
	
	return(output_List)


# Function for calculating cooperative index between seeds and nonseeds/sinks
def cooperation(seed1, seed2, nonseed1, nonseed2, sink1, sink2, intermediate1, intermediate2): 
		
	seed1_len = len(seed1)
	seed2_len = len(seed2)

	# With respect to organism 1
	intersection_seed1_nonseed2 = list(set(seed1) & set(nonseed2))
	intersection_seed1_nonseed2_len = len(intersection_seed1_nonseed2)
	intersection_seed1_sink2 = list(set(seed1) & set(sink2))
	intersection_seed1_sink2_len = len(intersection_seed1_sink2)
	intersection_seed1_intermediate2 = list(set(seed1) & set(intermediate2))	
	intersection_seed1_intermediate2_len = len(intersection_seed1_intermediate2)
	
	# With respect to organism 2
	intersection_seed2_nonseed1 = list(set(seed2) & set(nonseed1))
	intersection_seed2_nonseed1_len = len(intersection_seed2_nonseed1)
	intersection_seed2_sink1 = list(set(seed2) & set(sink1))	
	intersection_seed2_sink1_len = len(intersection_seed2_sink1)
	intersection_seed2_intermediate1 = list(set(seed2) & set(intermediate1))	
	intersection_seed2_intermediate1_len = len(intersection_seed2_intermediate1)

	# Calculate cooperation metrics
	nonseed_coop_index1 = float(intersection_seed1_nonseed2_len) / float(seed1_len) # As defined in Borenstein paper
	nonseed_coop_index2 = float(intersection_seed2_nonseed1_len) / float(seed2_len) # As defined in Borenstein paper
	sink_coop_index1 = float(intersection_seed1_sink2_len) / float(seed1_len) # Uses sink compounds instead on nonseeds
	sink_coop_index2 = float(intersection_seed2_sink1_len) / float(seed2_len) # Uses sink compounds instead on nonseeds
	intermediate_coop_index1 = float(intersection_seed1_intermediate2_len) / float(seed1_len) # Uses intermediate compounds instead on nonseeds
	intermediate_coop_index2 = float(intersection_seed2_intermediate1_len) / float(seed2_len) # Uses intermediate compounds instead on nonseeds

	cooperation_metric_list = [nonseed_coop_index1, nonseed_coop_index2, sink_coop_index1, sink_coop_index2, intermediate_coop_index1, intermediate_coop_index2]
	
	output_list = [intersection_seed1_nonseed2, intersection_seed1_sink2, intersection_seed2_nonseed1, intersection_seed2_sink1, cooperation_metric_list, intersection_seed1_intermediate2, intersection_seed2_intermediate1]
	
	return(output_list)
	
#---------------------------------------------------------------------------------------#

# Set up input files
parser = argparse.ArgumentParser(description='Calculate metabolic interaction from two sets of seed information.')

parser.add_argument('--name1', default='organism1', help='Name of first (meta)organism')
parser.add_argument('--name2', default='organism2', help='Name of second (meta)organism')

parser.add_argument('--files1', default='none', help='Directory of SCC network output for first (meta)organism')
parser.add_argument('--files2', default='none', help='Directory of SCC network output for second (meta)organism')

parser.add_argument('--confidence', default=0, help='Minimum confidence values for seeds and sinks to be considered in calculations')

parser.add_argument('--quiet', default='y', help='Turns on verbose output mode')


args = parser.parse_args()
name1 = args.name1
name2 = args.name2
sccfiles1 = args.files1
sccfiles2 = args.files2
confidence = args.confidence
quiet = args.quiet

if sccfiles1 ==  'none' or sccfiles2 ==  'none': sys.exit('WARNING: Missing input file(s), quitting')
if os.path.exists(sccfiles1) == False or os.path.exists(sccfiles2) == False: sys.exit('WARNING: Missing input file(s), quitting')

#---------------------------------------------------------------------------------------#

# Retrieve and read in all necessary files

# cd to first and open/read in all files
sccfiles1_directory = sccfiles1 + '/seeds'
os.chdir(sccfiles1_directory)

for line in open('seeds.tsv','r'):
	if len(line.split()) == 1:
		print 'WARNING: No confidence values provided in seed file 1\n'
		seed_list1 = [x.split()[0] for x in open('seeds.tsv','r').readlines()]
		confidence_provided_seed1 = 0
	else:
		seed_list1 = []
		confidence_provided_seed1 = 1
		for object in open('seeds.tsv','r'):
			if object.split()[1] >= confidence:
				seed_list1.append(object.split()[0])
	break

for line in open('sinks.tsv','r'):
	if len(line.split()) == 1:
		print 'WARNING: No confidence values provided in sink file 1\n'
		sink_list1 = [x.split()[0] for x in open('sinks.tsv','r').readlines()]
		confidence_provided_sink1 = 0
	else:
		sink_list1 = []
		confidence_provided_sink1 = 1
		for object in open('sinks.tsv','r'):
			if object.split()[1] >= confidence:
				sink_list1.append(object.split()[0])
	break	

for line in open('seedParameters.txt','r'):
	if len(line.split()) < 2: 
		continue

	elif line.split()[1] == 'Giant':
		org1_giant = line.split()[3]
		continue

	elif line.split()[0] == 'Minimum':
		org1_minimum = int(line.split()[3])
		continue

	elif line.split()[1] == 'Threshold:':
		org1_threshold = float(line.split()[2])
		break
	
	
nonseed_list1 = [x.split()[0] for x in open('nonseeds.tsv','r').readlines()]

intermediate_list1 = [x.split()[0] for x in open('intermediates.tsv','r').readlines()]


# go back to starting directory
os.chdir(starting_directory)


# then cd to second and do the same
sccfiles2_directory = sccfiles2 + '/seeds'
os.chdir(sccfiles2_directory)

for line in open('seeds.tsv','r'):
	if len(line.split()) == 1:
		print 'WARNING: No confidence values provided in seed file 1\n'
		seed_list2 = [x.split()[0] for x in open('seeds.tsv','r').readlines()]
		confidence_provided_seed2 = 0
	else:
		seed_list2 = []
		confidence_provided_seed2 = 1
		for object in open('seeds.tsv','r'):
			if object.split()[1] >= confidence:
				seed_list2.append(object.split()[0])
	break

for line in open('sinks.tsv','r'):
	if len(line.split()) == 1:
		print 'WARNING: No confidence values provided in sink file 1\n'
		sink_list2 = [x.split()[0] for x in open('sinks.tsv','r').readlines()]
		confidence_provided_sink2 = 0
	else:
		sink_list2 = []
		confidence_provided_sink2 = 1
		for object in open('sinks.tsv','r'):
			if object.split()[1] >= confidence:
				sink_list2.append(object.split()[0])
	break	


for line in open('seedParameters.txt','r'):
	if len(line.split()) < 2: continue

	elif line.split()[1] == 'Giant':
		org2_giant = line.split()[3]
		continue
	
	elif line.split()[0] == 'Minimum':
		org2_minimum = int(line.split()[3])
		continue

	elif line.split()[1] == 'Threshold:':
		org2_threshold = float(line.split()[2])
		break

	
nonseed_list2 = [x.split()[0] for x in open('nonseeds.tsv','r').readlines()]

intermediate_list2 = [x.split()[0] for x in open('intermediates.tsv','r').readlines()]


# Go back to starting directory
os.chdir(starting_directory)


#---------------------------------------------------------------------------------------#

# Calculate putative metabolic interactions and parse the output
competition_information = list(competition(seed_list1, seed_list2))
intersection_seeds = competition_information[0]
seed1_only = competition_information[1]
seed2_only = competition_information[2]
comp_indices = competition_information[3]
comp_index1v2 = comp_indices[0]
comp_index2v1 = comp_indices[1]

cooperation_information = list(cooperation(seed_list1, seed_list2, nonseed_list1, nonseed_list2, sink_list1, sink_list2, intermediate_list1, intermediate_list2))
intersection_seed1_nonseed2 = cooperation_information[0]
intersection_seed1_sink2 = cooperation_information[1]
intersection_seed2_nonseed1 = cooperation_information[2]
intersection_seed2_sink1 = cooperation_information[3]
cooperation_metric_list = cooperation_information[4]
intersection_seed1_intermediate2 = cooperation_information[5]
intersection_seed2_intermediate1 = cooperation_information[6]

nonseed_coop_index1v2 = cooperation_metric_list[0]
nonseed_coop_index2v1 = cooperation_metric_list[1]
sink_coop_index1v2 = cooperation_metric_list[2]
sink_coop_index2v1 = cooperation_metric_list[3]
intermediate_coop_index1v2 = cooperation_metric_list[4]
intermediate_coop_index2v1 = cooperation_metric_list[5]

# Calculate total nodes per network
node_total1 = len(seed_list1) + len(nonseed_list1)
node_total2 = len(seed_list2) + len(nonseed_list2)

#---------------------------------------------------------------------------------------#

# Write all output to strings

confidence_str = '''
Seed confidence threshold for comparison: {confidence}
	
-----------------------------------------------
'''.format(confidence=confidence)

organism_name_str = '''
Individual Network Information

Organism 1: {organism1}
Only Giant Component: {giant1}
Minimum Component Size: {component1}
Seed Threshold: {threshold1}

Organism 2: {organism2}
Only Giant Component: {giant2}
Minimum Component Size: {component2}
Seed Threshold: {threshold2}

-----------------------------------------------
'''.format(organism1 = name1, organism2 = name2, giant1 = org1_giant, giant2 = org2_giant, component1 = org1_minimum, component2 = org2_minimum, threshold1 = org1_threshold, threshold2 = org2_threshold)

network_summary_str = '''
Network component summaries

{organism1}
Seed count: {se1}
Sink count: {si1}
Intermediate count: {i1}
Node total: {n1}

{organism2}
Seed count: {se2}
Sink count: {si2}
Intermediate count: {i2}
Node total: {n2}

-----------------------------------------------
'''.format(se1 = str(len(seed_list1)), si1 = str(len(sink_list1)), i1 = str(len(intermediate_list1)), n1 = str(node_total1), se2 = str(len(seed_list2)), si2 = str(len(sink_list2)), i2 = str(len(intermediate_list2)), n2 = str(node_total2), organism1 = name1, organism2 = name2)


# Create readable strings of compounds for output
intersection_seeds_str = ''
intersection_seeds_count = 0
for object in intersection_seeds:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	intersection_seeds_str = intersection_seeds_str + '\n' + entry
	intersection_seeds_count += 1
if not intersection_seeds: 
	intersection_seeds_str = 'none'
	
seed1_only_str = ''
seed1_only_count = 0
for object in seed1_only:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	seed1_only_str = seed1_only_str + '\n' + entry
	seed1_only_count += 1
if not seed1_only: 
	seed1_only_str = 'none'

seed2_only_str = ''
seed2_only_count = 0
for object in seed2_only:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	seed2_only_str = seed2_only_str + '\n' + entry
	seed2_only_count += 1
if not seed2_only: 
	seed2_only_str = 'none'

competition_str = '''
Putative competition

{organism1} vs {organism2} competitive index: {comp1}
{organism2} vs {organism1} competitive index: {comp2}

{count1} overlapping seed compounds:
{overlap}

{count2} non-overlapping seed compounds in {organism1} but not {organism2}:
{nonoverlap1}

{count3} non-overlapping seed compounds in {organism2} but not {organism1}:
{nonoverlap2}

-----------------------------------------------
'''.format(count1=intersection_seeds_count, count2=seed1_only_count, count3=seed2_only_count, comp1=comp_index1v2, comp2=comp_index2v1, overlap=intersection_seeds_str, nonoverlap1=seed1_only_str, nonoverlap2=seed2_only_str, organism1 = name1, organism2 = name2)


intersection_seed1_sink2_str = ''
intersection_seed1_sink2_count = 0
for object in intersection_seed1_sink2:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	intersection_seed1_sink2_str = intersection_seed1_sink2_str + '\n' + entry
	intersection_seed1_sink2_count += 1
if not intersection_seed1_sink2: 
	intersection_seed1_sink2_str = 'none'

intersection_seed2_sink1_str = ''
intersection_seed2_sink1_count = 0
for object in intersection_seed2_sink1:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	intersection_seed2_sink1_str = intersection_seed2_sink1_str + '\n' + entry
	intersection_seed2_sink1_count += 1
if not intersection_seed2_sink1: 
	intersection_seed2_sink1_str = 'none'

sink_cooperation_str = '''
Putative cooperation

Seeds U Sinks:
Cooperative index of {organism1} with {organism2}: {coop_index1}
Cooperative index of {organism2} with {organism1}: {coop_index2}

{count1} compounds in seed of {organism1} and sink of {organism2}:
{coop_list1u2}

{count2} compounds in seed of {organism2} and sink of {organism1}:
{coop_list2u1}
'''.format(count1=intersection_seed1_sink2_count, count2=intersection_seed2_sink1_count, coop_index1=sink_coop_index1v2, coop_index2=sink_coop_index2v1, coop_list1u2=intersection_seed1_sink2_str, coop_list2u1=intersection_seed2_sink1_str, organism1 = name1, organism2 = name2)


intersection_seed1_nonseed2_str = ''
intersection_seed1_nonseed2_count = 0
for object in intersection_seed1_nonseed2:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	intersection_seed1_nonseed2_str = intersection_seed1_nonseed2_str + '\n' + entry
	intersection_seed1_nonseed2_count += 1
if not intersection_seed1_nonseed2: 
	intersection_seed1_nonseed2_str = 'none'

intersection_seed2_nonseed1_str = ''
intersection_seed2_nonseed1_count = 0
for object in intersection_seed2_nonseed1:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	intersection_seed2_nonseed1_str = intersection_seed2_nonseed1_str + '\n' + entry
	intersection_seed2_nonseed1_count += 1
if not intersection_seed2_nonseed1: 
	intersection_seed2_nonseed1_str = 'none'

nonseed_cooperation_str = '''
Seeds U Nonseeds:
Cooperative index of {organism1} with {organism2}: {coop_index1}
Cooperative index of {organism2} with {organism1}: {coop_index2}

{count1} compounds in seed of {organism1} and nonseed of {organism2}:
{coop_list1u2}

{count2} compounds in seed of {organism2} and nonseed of {organism1}:
{coop_list2u1}
'''.format(count1=intersection_seed1_nonseed2_count, count2=intersection_seed2_nonseed1_count, coop_index1=nonseed_coop_index1v2, coop_index2=nonseed_coop_index2v1, coop_list1u2=intersection_seed1_nonseed2_str, coop_list2u1=intersection_seed2_nonseed1_str, organism1 = name1, organism2 = name2)


intersection_seed1_intermediate2_str = ''
intersection_seed1_intermediate2_count = 0
for object in intersection_seed1_intermediate2:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	intersection_seed1_intermediate2_str = intersection_seed1_intermediate2_str + '\n' + entry
	intersection_seed1_intermediate2_count += 1
if not intersection_seed1_intermediate2: 
	intersection_seed1_intermediate2_str = 'none'

intersection_seed2_intermediate1_str = ''
intersection_seed2_intermediate1_count = 0
for object in intersection_seed2_intermediate1:
	entry = object + '\t' + str(compound_dict[object][0]).replace('_', ' ')
	intersection_seed2_intermediate1_str = intersection_seed2_intermediate1_str + '\n' + entry
	intersection_seed2_intermediate1_count += 1
if not intersection_seed2_intermediate1: 
	intersection_seed2_intermediate1_str = 'none'

intermediate_cooperation_str = '''
Seeds U Intermediate:
Cooperative index of {organism1} with {organism2}: {coop_index1}
Cooperative index of {organism2} with {organism1}: {coop_index2}

{count1} compounds in seed of {organism1} and an intermediate of {organism2}:
{coop_list1u2}

{count2} compounds in seed of {organism2} and an intermediate of {organism1}:
{coop_list2u1}
'''.format(count1=intersection_seed1_intermediate2_count, count2=intersection_seed2_intermediate1_count, coop_index1=intermediate_coop_index1v2, coop_index2=intermediate_coop_index2v1, coop_list1u2=intersection_seed1_intermediate2_str, coop_list2u1=intersection_seed2_intermediate1_str, organism1 = name1, organism2 = name2)

#---------------------------------------------------------------------------------------#

outfile_name = name1 + '_U_' + name2 + '.report'

# Print out the results to the screen
if quiet == 'n':
	print(confidence_str)
	print(organism_name_str)
	print(network_summary_str)
	print(competition_str)
	print(sink_cooperation_str)
	print(nonseed_cooperation_str)
	print(intermediate_cooperation_str)
	print('-----------------------------------------------\n\nWriting results to ' + outfile_name + '\n')

with open(outfile_name, 'w') as outfile:
	outfile.write(confidence_str)
	outfile.write(organism_name_str)
	outfile.write(network_summary_str)
	outfile.write(competition_str)
	outfile.write(sink_cooperation_str)
	outfile.write(nonseed_cooperation_str)
	outfile.write(intermediate_cooperation_str)

print('Done.\n')