#!/usr/bin/env python
'''USAGE: python crosstalk.py interaction.files --p n.s. --norm n
Multi-level inference of substrate competition and cooperation between transcriptome-informed genome-scale models
Calculates putative community-level and pair-wise metabolic interactions between species from aggregated bigSMALL analysis
'''

# Initialize all modules, functions, and compound dictionary
import sys
import numpy
import os
import argparse

#---------------------------------------------------------------------------------------#

# Find the directory where interaction.py is being run from
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Create a string with the name of the initial directory you start running the script from
starting_directory = str(os.getcwd())

#---------------------------------------------------------------------------------------#

# Set up arguments
parser = argparse.ArgumentParser(description='Calculate metabolic pair-wise and community-level interactions of species from the output of bigSMALL.')
parser.add_argument('input_file')
parser.add_argument('--p', default='n.s.', help='Minimum p-value for metabolites to be considered in calculations')
parser.add_argument('--norm', default='n', help='Normalize each metabolic model to total transcript recruited to each (y or n)')

args = parser.parse_args()
interactions = args.input_file
p_value = args.p
normalize = args.norm

if os.stat(interactions).st_size == 0 : sys.exit('WARNING: Input file empty, quitting')
if p_value != 'n.s.' and p_value < 0.0: sys.exit('WARNING: Invalid p-value cutoff, quitting')
if normalize != 'n' and normalize != 'y': sys.exit('WARNING: Invalid normalization response, quitting')

print('\n')

#---------------------------------------------------------------------------------------#

# Define functions

# Function to read in all species cominations from interaction file
def read_files(files):

	community = []
	for line in files:
		species = line.strip()
		if os.path.exists(species) == False:
			print('WARNING: ' + species + ' does not exist. Omitting combination.')
			continue
		community.append(species)

	interactions = []
	for index1 in community:
		for index2 in community:
			if index1 == index2:
				continue
			else:
				interactions.append([index1, index2])

	return interactions


# Reads importance files, applying p-value filter, normalizes score to reads, and generates a dictionary for compound names and compound scores
def read_scores(importance_scores, p_cutoff, norm):

	score_dictionary = {}
	name_dictionary = {}

	for line in importance_scores:
		
		line = line.split()
		if line[0] == 'Compound_code': continue
		
		compound_code = str(line[0])
		compound_name = str(line[1])
		
		score = float(line[2])
		if score == 0.0: continue
		
		p_value = str(line[3]).lstrip('<')
		if p_value == 'n.s.': p_value = 1
		p_value = float(p_value)
		if p_value > p_cutoff: continue

		score_dictionary[compound_code] = score
		name_dictionary[compound_code] = compound_name	

	if norm == 'y':
		final_score_dictionary = {}
		score_sum = sum([abs(x) for x in score_dictionary.values()])
		for index in score_dictionary.keys():
			final_score = score_dictionary[index] / score_sum
			final_score_dictionary[index] = [name_dictionary[index], final_score]
	else:
		final_score_dictionary = {}
		for index in score_dictionary.keys():
			final_score_dictionary[index] = [name_dictionary[index], score_dictionary[index]]
	
	return final_score_dictionary


# Function for calculating edges of metabolic competition
def single_interaction(score_dict_1, score_dict_2):
	
	all_compounds = list(set(score_dict_1.keys() + score_dict_2.keys()))
	
	interaction_dictionary = {}
	for index in all_compounds:

		# Consider only those compounds in the network of both organisms
		try:
			score_1 = float(score_dict_1[index][1])
			name = str(score_dict_1[index][0])
		except KeyError:
			continue
		try:
			score_2 = float(score_dict_2[index][1])
		except KeyError:
			continue

		# Determine type and strength of interaction
		if score_1 < 0 or score_2 < 0:
			temp_score_1 = 2**abs(score_1)
			temp_score_2 = 2**abs(score_2)
			ratio = min([(temp_score_1 / temp_score_2), (temp_score_2 / temp_score_1)])
			magnitude = (2**abs(score_1)) + (2**abs(score_2))
			interaction = -numpy.log2(ratio * magnitude)
		else:
			temp_score_1 = score_1
			temp_score_2 = score_2
			ratio = min([(temp_score_1 / temp_score_2), (temp_score_2 / temp_score_1)])
			magnitude = (2**score_1) + (2**score_2)
			interaction = numpy.log2(ratio * magnitude)

		interaction_dictionary[index] = [name, score_1, score_2, ratio, magnitude, interaction]

	return interaction_dictionary


# Calculate cumulative importance of each compound across the groups of models tested
def community_demand(community_dict, member_dict):

	for compound in member_dict.keys():
		name = member_dict[compound][0]
		score = member_dict[compound][1]
		
		if score < 0:
			score = -(2**abs(score))
		else:
			score = 2**abs(score)
	
		try:
			cumulative_score = community_dict[compound][1] + score
			if score > 0:
				consumption = community_dict[compound][2] + score
				production = community_dict[compound][3]
			elif score < 0:
				production = community_dict[compound][3] + score
				consumption = community_dict[compound][2]
			community_dict[compound] = [name, cumulative_score, consumption, production]
		except KeyError:
			if score > 0:
				consumption = score
				production = 0
			elif score < 0:
				consumption = 0
				production = score
			community_dict[compound] = [name, score, consumption, production]

	return community_dict


# Calculates the percentile for the given type of score
def calc_percentile(output_dictionary, type_index):

	current_list = []
	for index in output_dictionary.keys():
		current_list.append(output_dictionary[index][type_index])

	per_90 = [numpy.percentile(current_list, 10), numpy.percentile(current_list, 90)]
	per_80 = [numpy.percentile(current_list, 20), numpy.percentile(current_list, 80)]
	per_70 = [numpy.percentile(current_list, 30), numpy.percentile(current_list, 70)]
	per_60 = [numpy.percentile(current_list, 40), numpy.percentile(current_list, 60)]

	for index in output_dictionary.keys():
		if output_dictionary[index][type_index] < per_90[0] or output_dictionary[index][type_index] > per_90[1]:
			output_dictionary[index] = output_dictionary[index] + ['90']
			continue
		elif output_dictionary[index][type_index] < per_80[0] or output_dictionary[index][type_index] > per_80[1]:
			output_dictionary[index] = output_dictionary[index] + ['80']
			continue
		elif output_dictionary[index][type_index] < per_70[0] or output_dictionary[index][type_index] > per_70[1]:
			output_dictionary[index] = output_dictionary[index] + ['70']
			continue
		elif output_dictionary[index][type_index] < per_60[0] or output_dictionary[index][type_index] > per_60[1]:
			output_dictionary[index] = output_dictionary[index] + ['60']
			continue
		else:
			output_dictionary[index] = output_dictionary[index] + ['50']

	return output_dictionary


# Function to write data to output file
def write_output(header, output_dictionary, file_name, type_output):

	with open(file_name, 'w') as outfile:
		outfile.write(header)
		
		if type_output == 'single':

			for index in output_dictionary.keys():
				
				name = output_dictionary[index][0]
				score_1 = output_dictionary[index][1]
				score_2 = output_dictionary[index][2]
				ratio = output_dictionary[index][3]
				magnitude = output_dictionary[index][4]
				interaction = output_dictionary[index][5]
				percentile = output_dictionary[index][6]

				entry = '\t'.join([str(index), str(name), str(score_1), str(score_2), str(round(float(ratio), 3)), str(round(float(magnitude), 3)), str(round(float(percentile), 3))]) + '\n'
				outfile.write(entry)

		elif type_output == 'community':

			for index in output_dictionary.keys():
			
				name = output_dictionary[index][0]
				score = output_dictionary[index][1]
				consumption = output_dictionary[index][2]
				production = output_dictionary[index][3]
				percentile = output_dictionary[index][4]

				# Transform scores back to log2
				if score == 0.0:
					score = 0.0
				elif score < 0.0:
					score = numpy.log2(abs(score)) * -1
				else:
					score = numpy.log2((score))

				if consumption == 0.0:
					consumption = 0.0
				else:
					consumption = numpy.log2(consumption)

				if production == 0.0:
					production = 0.0
				else:
					production = numpy.log2(abs(production)) * -1

				entry = '\t'.join([str(index), str(name), str(score), str(round(float(consumption), 3)), str(round(float(production), 3)), str(round(float(percentile), 3))]) + '\n'
				outfile.write(entry)


#---------------------------------------------------------------------------------------#

# Worflow

# Retrieve and read in the necessary files
interactions = open(interactions, 'r')
interactions_list = read_files(interactions)
if not os.path.exists('community.files'):	
	os.makedirs('community.files')
community_dictionary = {}
community = []
current = 0
for index in interactions_list:

	os.chdir(index[0])
	scores_1 = read_scores(open('importances.tsv','r'), p_value, normalize)
	os.chdir(starting_directory)
	os.chdir(index[1])
	scores_2 = read_scores(open('importances.tsv','r'), p_value, normalize)
	os.chdir(starting_directory)

	current += 1
	print('Calculating interaction ' + str(current) + ' of ' + str(len(interactions_list)) + '.')
	interaction = single_interaction(scores_1, scores_2)
	interaction = calc_percentile(interaction, 5)

	if not str(index[0]) in community:
		community_dictionary = community_demand(community_dictionary, scores_1)
		community.append(str(index[0]))
	if not str(index[1]) in community:
		community_dictionary = community_demand(community_dictionary, scores_2)
		community.append(str(index[1]))

	header = 'compound_code\tcompound_name\tscore_1\tscore_2\tratio\tmagnitude\tinteraction_score\tpercentile\n'
	file_name = str('community.files/' + index[0]) + '.and.' + str(index[1]) + '.interaction.txt'
	write_output(header, interaction, file_name, 'single')

# Write cumulative scores to a file
community_dictionary = calc_percentile(community_dictionary, 1)
header = 'compound_code\tcompound_name\tcumulative_metabolite_score\tconsumption_score\tproduction_score\tpercentile\n'
file_name = 'community.files/community_importance.tsv'
write_output(header, community_dictionary, file_name, 'community')
print('Done\n')