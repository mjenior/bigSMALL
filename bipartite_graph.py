#!/usr/bin/env python
'''USAGE: bipartite_graph.py KO_expressionfile --name organism_name --min 0 --degree 0 --iters 1
The function of this script is to convert lists of genes to unweighted, 
directed graphs and compute importance of each compound to metabolism 
based on the expression of surrounding enzyme nodes.
'''

# Written by Matthew Jenior, University of Michigan, Schloss Laboratory, 2016

# Our equation for metabolite score is a variation of Eigenvector centrality and incorporates both 
# expression of local enzymes and degree of connectedness of each compound node in the calculation.  
# Relative importance of each surrounding compound node is then calculated by dividing the sum of 
# surrounding transcripts (the eigenvalue) by the number of edges connected to the node of interest.  
# This is repeated respective to incoming, outgoing, and combined edges.  A Negative Binomial
# distribution of simulated transcript abundance is then created and repeatedly subsampled from  
# to generate confidence interval for to compare the measured values against.

# Dependencies:  
# The script itself needs to be run from from a directory containing the /support/ sub-directory
# The only argument is a 2 column matrix text file containing a column of KO codes with corresponding expression
# Example:
# K00045		0
# K03454		4492
# K10021		183
# ...

# Generate files:  A new directory in ./ ending in ".bipartite.files" that contains all output including:
	# A 2 column directed, bipartite network file of compounds and enzymes
	# A text file containing reference errors thrown during the translation of KOs to chemical equations
	# A text file containing user defined parameters
	# List of unique compound nodes
	# List of unique enzymes nodes
	# 4 files with information derived from the network:
	#	1.  Input metabolite score (including confidence interval when applicable)
	#	2.  Output metabolite score (including confidence interval when applicable)
	#	3.  Description of network topology for each node (indegree and outdegree)

#---------------------------------------------------------------------------------------#		

# Import python modules
import sys
import os
import pickle
import math
import argparse
import random
import numpy
import time

#---------------------------------------------------------------------------------------#		

# Start timer
start = time.time()

#---------------------------------------------------------------------------------------#		

# Define arguments
parser = argparse.ArgumentParser(description='Generate bipartite metabolic models and calculates importance of substrate nodes based on gene expression.')
parser.add_argument('input_file')
parser.add_argument('--name', default='organism', help='Organism or other name for KO+expression file (default is organism)')
parser.add_argument('--min', default=0, help='minimum importance value to report')
parser.add_argument('--degree', default=0, help='minimum degree value to report')
parser.add_argument('--iters', default=1, help='iterations for random distribution subsampling')
args = parser.parse_args()

# Assign variables
KO_input_file = str(args.input_file)
file_name = str(args.name)
min_importance = int(args.min)
min_degree = int(args.degree)
iterations = int(args.iters)

#---------------------------------------------------------------------------------------#			

# Check for user error
if KO_input_file == 'input_file':
	print 'No KO+expression file provided. Aborting.'
	sys.exit()
elif os.stat(KO_input_file).st_size == 0:
	print('Empty input file provided. Aborting.')
	sys.exit()
elif min_importance < 0:
	print 'Invalid importance minimum. Aborting.'
	sys.exit()
elif min_degree < 0:
	print 'Invalid degree minimum. Aborting.'
	sys.exit()
elif iterations < 0:
	print 'Invalid iteration value. Aborting.'
	sys.exit()
	
#---------------------------------------------------------------------------------------#			

# Define functions

# Create a dictionary for transcript value associated with its KO
def transcription_dictionary(KO_file):
	
	seq_total = 0  # Total number of reads
	seq_max = 0  # Highest single number of reads
	transcript_dict = {}  # Dictionary for transcription
	
	KO_list = []
	for index in KO_file:
		index_split = index.split()
		seq_total += float(index_split[1])
		KO_list.append(str(index_split[0]))
		
		if not str(index_split[0]) in transcript_dict.keys():
			transcript_dict[str(index_split[0])] = float(index_split[1])
			if float(index_split[1]) > seq_max: seq_max = float(index_split[1])
		else:
			transcript_dict[str(index_split[0])] = transcript_dict[str(index_split[0])] + float(index_split[1])
			if transcript_dict[str(index_split[0])] > seq_max: seq_max = transcript_dict[str(index_split[0])]
			continue
	
	KO_list = list(set(KO_list))
	
	return transcript_dict, KO_list, seq_total, seq_max


# Translates a list of KOs to the bipartite graph
def translateKO(KOs, ko_dict, reaction_dict):

	# Set some starting points
	triedCountKO = 0
	excludedCountKO = 0
	triedCountReact = 0
	excludedCountReact = 0
	totalIncludedReact = 0
	
	compound_list = []
	enzyme_list = []
	network_list = []
	
	# Open a file to write the bipartite graph to
	graph = open('bipartite.graph', 'w')

	# Nested loops to convert the KO list to a directed graph of input and output compounds
	# Outside loop finds the biochemical reactions corresponding the the given KO	
	print('Translating KEGG orthologs to bipartite enzyme-to-compound graph...\n')
	
	errorfile = open('key_error.log', 'w')
	
	with open('bipartite.graph', 'w') as graph:
	
		for line in KOs:
	
			current_ko = line.strip('ko:')
			triedCountKO += 1
			
			try:
				reaction_number = ko_dict[current_ko]
			except KeyError:
				errorString = 'WARNING: ' + str(current_ko) + ' not found in KO-to-Reaction dictionary. Omitting.\n'
				errorfile.write(errorString)
				excludedCountKO += 1
				continue 
	
			# Inner loop translates the reaction codes to collections of input and output compounds
			for index in reaction_number:
				triedCountReact += 1
				try:
					reaction_collection = reaction_dict[index]
				except KeyError:
					errorString = 'WARNING: ' + str(index) + ' not found in Reaction-to-Compound dictionary. Omitting.\n'
					errorfile.write(errorString)
					excludedCountReact += 1
					continue
		
				# The innermost loop creates two columns of input and output compounds, incorporating reversibility information
				for x in reaction_collection:
				
					totalIncludedReact += 1
					
					# Split reaction input and output as well as the list of compounds with each
					reaction_info = x.split(':')
					input_compounds = reaction_info[0].split('|')
					output_compounds = reaction_info[2].split('|')
					rev = reaction_info[1].split('|')
						
					for input_index in input_compounds:
						network_list.append([str(input_index), str(current_ko)])
						graph.write(''.join([str(input_index), '\t', str(current_ko), '\n']))
						
						if rev == 'R':
							network_list.append([str(current_ko), str(input_index)])
							graph.write(''.join([str(current_ko), '\t', str(input_index), '\n']))
							
						compound_list.append(str(input_index))
						enzyme_list.append(str(current_ko))
			
			
					for output_index in output_compounds:
						network_list.append([str(current_ko), str(output_index)])
						graph.write(''.join([str(current_ko), '\t', str(output_index), '\n']))
						
						if rev == 'R':
							network_list.append([str(output_index), str(current_ko)])
							graph.write(''.join([str(output_index), '\t', str(current_ko), '\n']))
							
						compound_list.append(str(output_index))
						enzyme_list.append(str(current_ko))
	
	
	error_string = '''KOs successfully translated to Reactions: {KO_success}
KOs unsuccessfully translated to Reactions: {KO_failed}

Reactions successfully translated to Compounds: {Reaction_success}
Reactions unsuccessfully translated to Compounds: {Reaction_failed}
'''.format(KO_success = str(triedCountKO - excludedCountKO), KO_failed = str(excludedCountKO), Reaction_success = str(triedCountReact - excludedCountReact), Reaction_failed = str(excludedCountReact))

	errorfile.write(error_string)
	errorfile.close()
	
	network_list = [list(x) for x in set(tuple(x) for x in network_list)]  # List of unique edges (KOs and compounds)
	compound_list = list(set(compound_list))  # List of unique compounds
	enzyme_list = list(set(enzyme_list))  # List of unique enzymes
	
	print('Done.\n')
	
	return network_list, compound_list, enzyme_list


# Compile surrounding input and output node transcripts into a dictionary, same for degree information
def network_dictionaries(network, transcript_dictionary):

	# Open blank dictionaries to populate with compounds and their corresponding transcription
	input_dictionary = {}
	output_dictionary = {}
	indegree_dictionary = {}
	outdegree_dictionary = {}

	for edge_info in network:
	
		# Output, 'C' is useful because it is at the beginning of every compound code
		if edge_info[1][0] == 'C':
		
			# Fill output score dictionary
			if not edge_info[1] in output_dictionary.keys():
				try:
					current_transcript_count = transcript_dictionary[edge_info[0]]
				except KeyError:
					current_transcript_count = 0
		
				output_dictionary[edge_info[1]] = [current_transcript_count]
			else:
				output_dictionary[edge_info[1]].append(transcript_dictionary[edge_info[0]])
			
			# Fill indegree dictionary	
			if not edge_info[1] in indegree_dictionary.keys():
				indegree_dictionary[edge_info[1]] = 1
			else:	
				indegree_dictionary[edge_info[1]] = indegree_dictionary[edge_info[1]] + 1
			
			continue 
				
				
		# Input
		if edge_info[0][0] == 'C':
		
			if not edge_info[0] in input_dictionary.keys():
				try:
					current_transcript_count = transcript_dictionary[edge_info[1]]
				except KeyError:
					current_transcript_count = 0

				input_dictionary[edge_info[0]] = [current_transcript_count]
			
			else:
				input_dictionary[edge_info[0]].append(transcript_dictionary[edge_info[1]])
			
			# Fill outdegree dictionary
			if not edge_info[0] in outdegree_dictionary.keys():
				outdegree_dictionary[edge_info[0]] = 1
			else:
				outdegree_dictionary[edge_info[0]] = outdegree_dictionary[edge_info[0]] + 1
	
	all_degree_dictionary = {}
	all_transcript_dictionary = {}
	for index in transcript_dictionary.keys():
		
		try:
			input_score = sum(input_dictionary[index])
		except KeyError:
			input_score = 0
		try:	
			output_score = sum(output_dictionary[index])
		except KeyError:
			output_score = 0	
		all_transcript_dictionary[index] = [index, input_score, output_score]
		
		try:
			outdegree = sum(outdegree_dictionary[index])
		except KeyError:
			outdegree = 0
		try:	
			indegree = sum(indegree_dictionary[index])
		except KeyError:
			indegree = 0
		all_degree_dictionary[index] = [index, indegree, outdegree]
		
	return all_transcript_dictionary, all_degree_dictionary


# Calculate input and output scores and well as degree of each compound node
def calculate_importance(transcript_dict, degree_dict, compound_dict, min_score, min_deg):
	
	# Open blank lists for all calculated values
	inputscore_list = []
	outputscore_list = []
	degree_list = []
		
	# Calculate cumulative scores for all compounds as inputs or outputs
	for index in compound_dict.keys():
		
		compound = compound_dict[index]
		
		try:
			indegree = int(degree_dict[index][1])
		except KeyError:
			indegree = 0		
		try:
			outdegree = int(degree_dict[index][2])
		except KeyError:
			outdegree = 0	
			
		
		try:
			input_scores = [int(x) for x in transcript_dict[index][1]]
		except KeyError:
			input_scores = [0]
		if outdegree == 0:
			final_score = 0
		else:
			final_score = float(sum(input_scores) / outdegree)
			if indegree == 0:
				final_score = final_score * (1 + outdegree) * outdegree
		if final_score >= min_score:
			inputscore_list.append([compound, index, str(final_score)])
		
		
		try:
			output_scores = [int(x) for x in transcript_dict[index][2]]
		except KeyError:
			output_scores = [0]
		if indegree == 0:
			final_score = 0
		else:
			final_score = float(sum(input_scores) / indegree)
			if outdegree == 0:
				final_score = final_score * (1 + indegree) * indegree
		if final_score >= min_score:
			outputscore_list.append([compound, index, str(final_score)])

					
		if indegree >= min_deg and outdegree >= min_deg:
			degree_list.append([compound, index, str(indegree), str(outdegree)])

	return inputscore_list, outputscore_list, degree_list
	

# Perform iterative Monte Carlo simulation to create confidence interval for compound importance values
def monte_carlo_sim(network, kos, iterations, compound_dict, min_importance, min_degree, seq_total, seq_max):
	
	gene_count = len(kos)
	probability = 1.0 / gene_count
	
	distribution = list(numpy.random.negative_binomial(1, probability, seq_total))  # Negative Binomial distribution
	distribution = [i for i in distribution if i < seq_max]  # screen for transcript mapping greater than largest value actually sequenced
	
	input_distribution_dict = {}
	output_distribution_dict = {}
	
	increment = 100.0 / float(iterations)
	progress = 0.0
	sys.stdout.write('\rProgress: ' + str(progress) + '%')
	sys.stdout.flush()
	
	for index in compound_dict.keys():
		input_distribution_dict[index] = []
		output_distribution_dict[index] = []
	
	for current in range(0, iterations):
			
		sim_transcriptome = random.sample(distribution, gene_count)

		sim_transcript_dict = {}
		for index in range(0, gene_count):
			sim_transcript_dict[kos[index]] = sim_transcriptome[index]
		
		substrate_dict, degree_dict = network_dictionaries(network, sim_transcript_dict)
		inputscore_list, outputscore_list, degree_list = calculate_importance(substrate_dict, degree_dict, compound_dict, min_importance, min_degree)
		
		# Make dictionaries of scores for each compound for each direction
		for index in inputscore_list:
			input_distribution_dict[index[1]].append(float(index[2]))
		for index in outputscore_list:
			output_distribution_dict[index[1]].append(float(index[2]))
		
		progress += increment
		sys.stdout.write('\rProgress: ' + str(progress) + '%')
		sys.stdout.flush()
	
	print '\n'
		
	# Compile the scores for each compound and take the mean and standard deviation
	input_interval_list = []
	output_interval_list = []
	for index in compound_dict.keys():

		input_current_mean = float("%.3f" % (numpy.mean(input_distribution_dict[index])))
		input_current_std = float("%.3f" % (numpy.std(input_distribution_dict[index])))
		input_interval_list.append([index, input_current_mean, input_current_std])

		output_current_mean = float("%.3f" % (numpy.mean(output_distribution_dict[index])))
		output_current_std = float("%.3f" % (numpy.std(output_distribution_dict[index])))
		output_interval_list.append([index, output_current_mean, output_current_std])

	return input_interval_list, output_interval_list


# Assesses measured values against confidence interval from Monte Carlo simulation 
def confidence_interval(importance, interval):

	labeled_confidence = []

	for index in range(0, len(importance)):
		
		current_importance = importance[index]
		
		if float(importance[index][2]) > float(interval[index][1]):
		
			if float(importance[index][2]) > (float(interval[index][1]) + float(interval[index][2])):
			
				if float(importance[index][2]) > (float(interval[index][1]) + (float(interval[index][2]) * 2)):
				
					if float(importance[index][2]) > (float(interval[index][1]) + (float(interval[index][2]) * 3)):
						current_importance.extend((interval[index][1], interval[index][2], '+', '***'))
						labeled_confidence.append(current_importance)		
					else:
						current_importance.extend((interval[index][1], interval[index][2], '+', '**'))
						labeled_confidence.append(current_importance)
				else:
					current_importance.extend((interval[index][1], interval[index][2], '+', '*'))
					labeled_confidence.append(current_importance)
			else:
				current_importance.extend((interval[index][1], interval[index][2], '+', 'n.s.'))
				labeled_confidence.append(current_importance)
				
		elif float(importance[index][2]) < float(interval[index][1]):

			if float(importance[index][2]) < (float(interval[index][1]) - float(interval[index][2])):
			
				if float(importance[index][2]) < (float(interval[index][1]) - (float(interval[index][2]) * 2)):
					
					if float(importance[index][2]) < (float(interval[index][1]) - (float(interval[index][2]) * 3)):
						current_importance.extend((interval[index][1], interval[index][2], '-', '***'))
						labeled_confidence.append(current_importance)
					else:
						current_importance.extend((interval[index][1], interval[index][2], '-', '**'))
						labeled_confidence.append(current_importance)
				else:
					current_importance.extend((interval[index][1], interval[index][2], '-', '*'))
					labeled_confidence.append(current_importance)
			else:
				current_importance.extend((interval[index][1], interval[index][2], '-', 'n.s.'))
				labeled_confidence.append(current_importance)
	
	return labeled_confidence


# Function to write data to output files	
def write_output(header, out_data, file_name):

	with open(file_name, 'w') as outfile: 
		
		outfile.write(header)
			
		for index in out_data:
			index = [str(x) for x in index]
			index[-1] = str(index[-1]) + '\n'
			outfile.write('\t'.join(index))


##########################################################################################		
#																						 #
# 									  DO THE WORK										 #
#																						 #
##########################################################################################		


if file_name != 'organism':
	print '\nImputing metabolism for ' + file_name + '\n'

# Read in and create dictionary for transcription
with open(KO_input_file, 'r') as KO_file:
	transcript_dict, KO_list, total, max = transcription_dictionary(KO_file)

#---------------------------------------------------------------------------------------#		

# Determine starting directory
starting_directory = str(os.getcwd())
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Create and navigate to new output directory
directory = str(os.getcwd()) + '/' + file_name + '.bipartite.files'
if not os.path.exists(directory):	
	os.makedirs(directory)
os.chdir(directory)

#---------------------------------------------------------------------------------------#		

# Write parameters to a file

with open('parameters.txt', 'w') as parameter_file:
	outputString = '''User Defined Parameters
KO expression file: {ko}
Graph name: {name}
Minimum compound importance: {imp}
Minimum edges per node: {deg}
Monte Carlo simulation iterations: {iter}
'''.format(ko=str(KO_input_file), name=str(file_name), imp=str(min_importance), deg=str(min_degree), iter=str(iterations))
	parameter_file.write(outputString)

#---------------------------------------------------------------------------------------#		

# Create a dictionary of KO expression scores and load KEGG dictionaries

print('\nReading in KEGG dictionaries...\n')
# Read in pickled KO to reaction dictionary
ko_reactionpkl_path = script_path + '/support/ko_reaction.pkl'
ko_dictionary = pickle.load(open(ko_reactionpkl_path, 'rb'))

# Read in pickled reaction to reaction_mapformula dictionary
#reaction_mapformulapkl_path = script_path + '/support/reaction_mapformula.pkl'
reaction_mapformulapkl_path = script_path + '/support/reaction_mapformula_nonrev.pkl'
reaction_dictionary = pickle.load(open(reaction_mapformulapkl_path, 'rb'))

# Read in pickled compound name dictionary
compoundpkl_path = script_path + '/support/compound.pkl'
compound_dictionary = pickle.load(open(compoundpkl_path, 'rb'))
print('Done.\n')

#---------------------------------------------------------------------------------------#	

# Call translate function and separate output lists
reaction_graph, compound_list, enzyme_list = translateKO(KO_list, ko_dictionary, reaction_dictionary)

# Write compounds and enzymes to files
with open('compound.lst', 'w') as compound_file:
	for index in compound_list:
		compound_file.write(''.join([index, '\n']))

with open('enzyme.lst', 'w') as enzyme_file:
	for index in enzyme_list:
		enzyme_file.write(''.join([index, '\n']))
			
# Calculate actual importance scores for each compound in the network
print 'Calculating compound node connectedness and metabolite scores...\n'
transcript_dict, degree_dict = network_dictionaries(reaction_graph, transcript_dict)
inputscore_list, outputscore_list, degree_list = calculate_importance(transcript_dict, degree_dict, compound_dictionary, min_importance, min_degree)
print 'Done.\n'

#---------------------------------------------------------------------------------------#		

# Calculate simulated importance values if specified
if iterations > 1:

	print 'Comparing to simulated transcript distribution...\n'
	input_interval_list, output_interval_list = monte_carlo_sim(reaction_graph, enzyme_list, iterations, compound_dictionary, min_importance, min_degree, total, max)
	print '\nDone.\n'
	
	# Write all the calculated data to files
	print 'Writing score data with Monte Carlo simulation to files...\n'
	outname = file_name + '.input_score.monte_carlo.txt'
	final_output = confidence_interval(inputscore_list, input_interval_list)
	write_output('Compound_name	Compound_code	Input_metabolite_score	Simulated_Mean	Simulated_Std_Dev	Relationship_to_Mean	Significance\n', final_output, outname)
	
	outname = file_name + '.output_score.monte_carlo.txt'
	final_output = confidence_interval(outputscore_list, output_interval_list)
	write_output('Compound_name	Compound_code	Output_metabolite_score	Simulated_Mean	Simulated_Std_Dev	Relationship_to_Mean	Significance\n', final_output, outname)
	print 'Done.\n'
	
else:
	print 'Writing score data to files...\n' 
	outname = file_name + '.input_score.txt'
	write_output('Compound_name	Compound_code	Input_metabolite_score\n', inputscore_list, outname)

	outname = file_name + '.output_score.txt'
	write_output('Compound_name	Compound_code	Output_metabolite_score\n', outputscore_list, outname)
	print 'Done.\n'

#---------------------------------------------------------------------------------------#		

# Write network topology info to files
print 'Writing degree information to files...\n' 
outname = file_name + '.topology.txt'
write_output('Compound_name	Compound_code	Indegree	Outdegree\n', degree_list, outname)
print 'Done.\n'

# Return to the directory the script was called to
os.chdir(starting_directory)	

# Report time if iterations are performed
end = time.time()
if iterations > 10:
	duration = str(int(end - start))
	print '\nCompleted in ' + duration + ' seconds.\n'
else :
	print '\n'
	
print 'Output files located in: ' + directory + '\n\n'
