#!/usr/bin/env python
'''USAGE: python python interactions.py --interactions interaction.files --p_value n.s.
Calculates potential metabolic crosstalk within communities of bacteria based on output from bigSMALL
'''

# Initialize all modules, functions, and compound dictionary
import sys
import math
import os
import pickle
import argparse

#---------------------------------------------------------------------------------------#

# Find the directory where interaction.py is being run from
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Create a string with the name of the initial directory you start running the script from
starting_directory = str(os.getcwd())

#---------------------------------------------------------------------------------------#

# Define some functions

# Function to read in all species cominations from interaction file
def read_files(file_names):

	names = []
	omit = []

	for line in file_names:

		species_1 = str(line.split()[0])
		species_2 = str(line.split()[1])
		if species_1 in omit or species_2 in omit:
			continue

		if os.path.exists(species_1) == False:
			print('WARNING: ' + species_1 + ' does not exist. Omitting combinations.')
			omit.append(species_1)
			continue
		elif os.path.exists(species_2) == False:
			print('WARNING: ' + species_2 + ' does not exist. Omitting combinations.')
			omit.append(species_2)
			continue

		entry = [species_1, species_2]
		names.append(entry)
	
	print('\n')
	return names


# Function to get the organism name from the bigSMALL parameters file
def read_parameters(parameters):

	for line in parameters:

		line = line.split()

		if line[0] == 'KO expression file:':
			file_name = str(line[1])

		elif line[0] == 'Graph name:':
			graph_name = str(line[1])

	if graph_name != 'organism':
		name = graph_name
	else:
		name = file_name

	return name


# Reads importance files, converts to catagories, and generates a dictionary
def convert_scores(importance_scores, p_cutoff):

	score_list = []
	entry_list = []

	for line in compound_scores:
		
		line = line.split()
		if line[0] == 'Compound_code': continue
		
		score_list.append(float(line[2]))

		p_value = line[3]
		if line[3] == '<0.01':
			p_value = 0.01
		elif line[3] == '<0.05':
			p_value = 0.05
		else:
			p_value = 1

		entry_list.append([line[0],line[1],line[2],p_value])
		
	min_level = min(score_list) / 3
	output_hi = [min(score_list), (min(score_list) - min_level - 0.001)]
	output_med = [(min(score_list) - min_level), (min(score_list) - ( 2 * min_level) - 0.001)]
	output_lo = [(min(score_list) - ( 2 * min_level)),  -0.001]

	max_level = max(score_list) / 3
	input_hi = [max(score_list), (max(score_list) - min_level + 0.001)]
	input_med = [(max(score_list) - max_level), (max(score_list) - ( 2 * max_level) + 0.001)]
	input_lo = [(max(score_list) - ( 2 * max_level)),  0.001]

	score_dictionary = {}
	for index in entry_list:

		if index[3] > p_cutoff: continue

		if index[2] > 0.0:

			if input_hi[1] <= index[2] >= input_hi[0]:
				score_dictionary[index[0]] = [index[1], 3.0]
			elif input_hi[1] <= index[2] >= input_hi[0]:
				score_dictionary[index[0]] = [index[1], 2.0]
			elif input_hi[1] <= index[2] >= input_hi[0]:
				score_dictionary[index[0]] = [index[1], 1.0]

		elif index[2] < 0.0:

			if output_hi[1] >= index[2] <= output_hi[0]:
				score_dictionary[index[0]] = [index[1], -3.0]
			elif output_hi[1] >= index[2] <= output_hi[0]:
				score_dictionary[index[0]] = [index[1], -2.0]
			elif output_hi[1] >= index[2] <= output_hi[0]:
				score_dictionary[index[0]] = [index[1], -1.0]

		else:
			score_dictionary[index[0]] = [index[1], 0.0]


	return score_dictionary


# Function for calculating likely metabolic interaction
def interaction(score_dict1, score_dict2):
	
	all_compounds = list(set(score_dict1.keys() + score_dict2.keys()))
	
	interaction_list = []

	for index in all_compounds:
		
		name1 = 'place_holder'
		name2 = 'place_holder'

		try:
			name1 = score_dict1[index][0]
			score1 = float(score_dict1[index][1])
		except keyError:
			score1 = 0.0
			
		try:
			name2 = score_dict2[index][0]
			score2 = float(score_dict2[index][1])
		except keyError:
			score2 = 0.0

		if name1 == 'place_holder':
			name = name2
		else:
			name = name1

		interaction_score = str(score1 + score2)
	
		entry = [index, name, str(score1), str(score2), interaction_score]
		interaction_list.append(entry)

	return interaction_list


# Function to write data to output file
def write_output(header, out_data, p_cutoff, file_name):

	with open(file_name, 'w') as outfile: 
		
		p_value = str(p_cutoff) + '\n'
		outfile.write(header)
			
		for index in out_data:
			index = index.append(p_value)
			outfile.write('\t'.join(index))	
	
#---------------------------------------------------------------------------------------#

# Set up arguments
parser = argparse.ArgumentParser(description='Calculate metabolic pair-wise interactions of species from the output of bigSMALL.')
parser.add_argument('--interactions', default='none', help='2 column list of species interactors with bigSMALL output for each (directories)')
parser.add_argument('--p_value', default='n.s.', help='Minimum p-value for metabolites to be considered in calculations')

args = parser.parse_args()
interactions = args.interactions
p_value = float(args.p_value)

if interactions ==  'none': sys.exit('WARNING: Missing input file, quitting')
if p_value != 'n.s.' and p_value < 0.0: sys.exit('WARNING: p-value cutoff is less than 0, quitting')

#---------------------------------------------------------------------------------------#

# Retrieve and read in the necessary files
print('\nReading in interaction file.\n')
all_interations = read_files(interactions)

interaction_count = len(all_interations)
current = 1

for index in all_interations:

	print('Calculating interaction ' + str(current) + ' of ' + str(interaction_count) + '.')

	# Retrieve individual species information
	species_1 = index[0]
	species_2 = index[1]
	os.chdir(species_1)
	name_1 = read_parameters(open('parameters.txt','r'))
	scores_1 = convert_scores(open('importances.tsv','r'), p_value)
	os.chdir(starting_directory)
	os.chdir(species_2)
	name_2 = read_parameters(open('parameters.txt','r'))
	scores_2 = convert_scores(open('importances.tsv','r'), p_value)
	os.chdir(starting_directory)
	
	# Calculate putative metabolic interactions and parse the output
	print('Calculating putative interaction of ' + name_1 + ' and ' + name_2 + '.\n')
	crosstalk = interaction(scores_1, scores_2)
	current += 1
		
	# Write output tables and summary to files
	head = 'Compound_code\tCompound_name\tSpecies_1_score\tSpecies_2_score\tInteraction_score\tp_value\n'
	file = name_1 + 'AND' + name_2 + '.interaction.tsv'
	write_output(head, crosstalk, p_value, file)


print('Done.\n')


