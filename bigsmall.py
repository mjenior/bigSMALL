#!/usr/bin/env python
'''USAGE: bigsmall.py KO_expressionfile --name organism_name --iters 1000
The function of this script is to convert lists of genes to unweighted, 
directed graphs and compute importance of each compound to metabolism 
based on the expression of surrounding enzyme nodes.
'''

# Written by Matthew Jenior, University of Michigan, Schloss Laboratory, 2016

# Dependencies:  
# The script itself needs to be run from from a directory containing the /support/ sub-directory
# The only argument is a 2 column matrix text file containing a column of KO codes with corresponding expression
# Example:
# K00045		0
# K03454		4492
# K10021		183
# ...
# Knnnnn 		n

# Generate files:  A new directory in ./ ending in ".bipartite.files" that contains all output including:
	# A 2 column directed, bipartite network file of compounds and enzymes
	# A text file containing reference errors thrown during the translation of KOs to chemical equations
	# A text file containing user defined parameters
	# List of unique compound nodes
	# List of unique enzymes nodes
	# Table containing network topology
	# Table containing importance values and significance (when applicable)

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
import datetime

#---------------------------------------------------------------------------------------#		

# Start timer
start = time.time()

#---------------------------------------------------------------------------------------#		

# User defined arguments
parser = argparse.ArgumentParser(description='Generate bipartite metabolic models and calculates importance of substrate nodes based on gene expression.')
parser.add_argument('input_file')
parser.add_argument('--name', default='default', help='Organism or other name for KO+expression file (default is organism)')
parser.add_argument('--iters', default='1000', help='Number of iterations of probability distribution for score comparison')
args = parser.parse_args()

# Assign variables
KO_input_file = str(args.input_file)
file_name = str(args.name)
iterations = int(args.iters)

#---------------------------------------------------------------------------------------#			

# Check if the user fucked it up
if KO_input_file == 'input_file':
	print('No KO+expression file provided. Aborting.')
	sys.exit()
elif os.stat(KO_input_file).st_size == 0:
	print('Empty input file provided. Aborting.')
	sys.exit()
elif file_name == '':
	print('Invalid names argument provided. Aborting.')
	sys.exit()
elif iterations < 0:
	print('Invalid iterations value. Aborting.')
	sys.exit()

# Make sure no spaces are in the name argument
file_name = file_name.replace(' ', '_')
	
#---------------------------------------------------------------------------------------#			

# Define all the functions!

# Create a dictionary for transcript value associated with its KO
def transcription_dictionary(KO_file):
	
	seq_total = 0  # Total number of reads
	seq_max = 0  # Highest single number of reads
	transcript_dict = {}  # Dictionary for transcription
	
	for line in KO_file:
		entry = line.split()
		
		ko = str(entry[0]).strip('ko:')
		expression = float(entry[1])
		
		seq_total += expression
		
		if not ko in transcript_dict.keys():
			transcript_dict[ko] = expression
		else:
			transcript_dict[ko] = transcript_dict[ko] + expression
		
		if transcript_dict[ko] > seq_max: seq_max = transcript_dict[ko]
	
	return transcript_dict, seq_total, seq_max


# Translates a list of KOs to the bipartite graph
def network_dictionaries(KOs, ko_dict, reaction_dict):

	# Set some starting points
	triedCountKO = 0
	excludedCountKO = 0
	triedCountReact = 0
	excludedCountReact = 0
	totalIncludedReact = 0
	
	network_list = []
	compound_lst = []
	KO_lst = []
	
	ko_input_dict = {}
	ko_output_dict = {}

	# Nested loops to convert the KO list to a directed graph of input and output compounds
	# Outside loop finds the biochemical reactions corresponding the the given KO	
	print('Translating KEGG orthologs to bipartite enzyme-to-compound graph...\n')
	
	with open('key_error.log', 'w') as errorfile:

		for current_ko in KOs:
	
			triedCountKO += 1
			
			if not current_ko in ko_input_dict:
				ko_input_dict[current_ko] = []
				ko_output_dict[current_ko] = []
			
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
				KO_lst.append(current_ko)
				for x in reaction_collection:
				
					totalIncludedReact += 1
					
					# Split reaction input and output as well as the list of compounds with each
					reaction_info = x.split(':')
					input_compounds = reaction_info[0].split('|')
					output_compounds = reaction_info[2].split('|')
					rev = reaction_info[1].split('|')
						
					for input_index in input_compounds:
						network_list.append([str(input_index), str(current_ko)])
						ko_input_dict[current_ko].append(str(input_index))
						
						if rev == 'R':
							network_list.append([str(current_ko), str(input_index)])
							ko_output_dict[current_ko].append(str(input_index))	
							
						compound_lst.append(str(input_index))		
			
					for output_index in output_compounds:
						network_list.append([str(current_ko), str(output_index)])
						ko_output_dict[current_ko].append(str(output_index))
						
						if rev == 'R':
							network_list.append([str(output_index), str(current_ko)])
							ko_input_dict[current_ko].append(str(output_index))
							
						compound_lst.append(str(output_index))
								
		error_string = '''KOs successfully translated to Reactions: {KO_success}
KOs unsuccessfully translated to Reactions: {KO_failed}

Reactions successfully translated to Compounds: {Reaction_success}
Reactions unsuccessfully translated to Compounds: {Reaction_failed}
'''.format(KO_success = str(triedCountKO - excludedCountKO), KO_failed = str(excludedCountKO), Reaction_success = str(triedCountReact - excludedCountReact), Reaction_failed = str(excludedCountReact))
		errorfile.write(error_string)
	
	network_list = [list(x) for x in set(tuple(x) for x in network_list)]  # List of unique edges (KOs and compounds)
	compound_lst = list(set(compound_lst))
	KO_lst = list(set(KO_lst))
	
	errorfile.close()
	print('Done.\n')
	
	return network_list, ko_input_dict, ko_output_dict, compound_lst, KO_lst


# Compile surrounding input and output node transcripts into a dictionary, same for degree information
def compile_transcripts(transcript_dictionary, ko_input_dict, ko_output_dict, compound_lst, KO_lst):

	compound_transcript_dict = {}
	compound_degree_dict = {}
	for compound in compound_lst:
		compound_transcript_dict[compound] = [0, 0] # [input, output]
		compound_degree_dict[compound] = [0, 0] # [indegree, outdegree]
		
	for ko in KO_lst:
	
		transcription = transcript_dictionary[ko]
		
		input_compounds = ko_input_dict[ko]
		output_compounds = ko_output_dict[ko]
		
		# Add transcription
		for compound in input_compounds:
			compound_transcript_dict[compound][0] = compound_transcript_dict[compound][0] + transcription
			compound_degree_dict[compound][1] = compound_degree_dict[compound][1] + 1
		
		for compound in output_compounds:
			compound_transcript_dict[compound][1] = compound_transcript_dict[compound][1] + transcription
			compound_degree_dict[compound][0] = compound_degree_dict[compound][0] + 1
	
	return compound_transcript_dict, compound_degree_dict


# Calculate input and output scores and well as degree of each compound node
def calculate_score(compound_transcript_dict, compound_degree_dict, compound_name_dict, compound_lst):
	
	score_dict = {}
	degree_dict = {}
		
	# Calculate metabolite scores integrating input and output reactions weightings
	for compound in compound_lst:
	
		score_dict[compound] = []
		degree_dict[compound] = []
		
		compound_name = compound_name_dict[compound]
		indegree = compound_degree_dict[compound][0]
		outdegree = compound_degree_dict[compound][1]
		input_transcription = compound_transcript_dict[compound][0]
		output_transcription = compound_transcript_dict[compound][1]	
		
		if outdegree == 0.0:
			input_score = 0.0
		else:
			input_score = input_transcription / outdegree

		if indegree == 0.0:
			output_score = 0.0
		else:
			output_score = output_transcription / indegree
		
		score_difference = input_score - output_score

		if score_difference < 1 and score_difference > -1:
			final_score = 0.0
		elif score_difference <= -1:
			final_score = math.log(abs(score_difference), 2) * -1
		else:
			final_score = math.log(score_difference, 2)

		final_score = float("%.3f" % final_score)

		score_dict[compound].extend((compound_name, final_score))
		degree_dict[compound].extend((compound_name, indegree, outdegree))	
					
	return score_dict, degree_dict

	
# Perform iterative simulation to create confidence interval for compound importance values
def probability_distribution(ko_input_dict, ko_output_dict, degree_dict, kos, compound_name_dict, seq_total, seq_max, compound_lst, transcription_dict, iterations):
	
	# Screen transcript distribution for those KOs included in the metabolic network
	transcript_distribution = []
	for index in kos:
		transcript_distribution.append(int(transcription_dict[index]))

	print 'Permuting transcript distributions...\n'
	increment = 100.0 / float(iterations) 
	progress = 0.0
	sys.stdout.write('\rProgress: ' + str(progress) + '%')
	sys.stdout.flush() 
	all_distributions = set()
	for index in range(iterations):

		# Generate bootstrapped transcript distributions
		transcript_distribution = random.sample(transcript_distribution, len(kos))

		if not tuple(transcript_distribution) in all_distributions:
			all_distributions.add(tuple(transcript_distribution))
			progress += increment
			progress = float("%.3f" % progress)
			sys.stdout.write('\rProgress: ' + str(progress) + '%')
			sys.stdout.flush() 

	all_distributions = list(all_distributions)
	sys.stdout.write('\rDone.                       \n\n')
	
	distribution_dict = {}
	for compound in compound_lst:
		distribution_dict[compound] = []

	# Memory intensive
	print 'Calculating importance scores for ' + str(iterations) + ' probability distributions...\n'
	progress = 0.0
	sys.stdout.write('\rProgress: ' + str(progress) + '%')
	sys.stdout.flush() 
	for index in all_distributions:

		current_distribution = list(index)

		sim_transcript_dict = {}
		for index in range(len(kos)):
			sim_transcript_dict[kos[index]] = current_distribution[index]

		substrate_dict, degree_dict = compile_transcripts(sim_transcript_dict, ko_input_dict, ko_output_dict, compound_lst, kos)
		score_dict, degree_dict = calculate_score(substrate_dict, degree_dict, compound_name_dict, compound_lst)
		
		# Make dictionaries of scores for each compound for each direction
		for compound in compound_lst:
			distribution_dict[compound].append(score_dict[compound][1])

		progress += increment
		progress = float("%.3f" % progress)
		sys.stdout.write('\rProgress: ' + str(progress) + '%')
		sys.stdout.flush() 

	sys.stdout.write('\rDone.                       \n\n')


	# Code for distribution testing purposes (uncomment when necessary)
	#with open('test_distribution.txt', 'w') as sample_dist_file:
	#    dist1 = list(distribution_dict['C00369']) # Starch - change as needed
	#    dist2 = list(distribution_dict['C00118']) # D-Glyceraldehyde_3-phosphate - change as needed
	#    for index in range(0, len(dist1)):
	#		entry1 = str(dist1[index])
	#		entry2 = str(dist2[index])
	#		sample_entry = entry1 + '\t' + entry2 + '\n'
	#		sample_dist_file.write(sample_entry)


	print 'Calculating summary statistics of each importance score distribution...\n'
	# Compile the scores for each compound and find the median and standard deviation
	interval_lst = []

	for compound in compound_lst:

		# Get the distribution
		current_dist = list(set(distribution_dict[compound]))

		# Calculate median
		current_median = numpy.median(current_dist)

		# McGill et al. (1978). Variations of Box Plots. The American Statistician, 32:1, 12-16.
		lower_iqr, upper_iqr = numpy.percentile(current_dist, [25, 75])
		numerator = 1.25 * abs(upper_iqr - lower_iqr)
		denominator = 1.35 * math.sqrt(len(current_dist))
		range_factor = numerator / denominator
		range_95 = 1.6 * range_factor
		range_99 = 1.9 * range_factor
		lower_95 = float("%.3f" % (current_median - range_95))
		upper_95 = float("%.3f" % (current_median + range_95))
		lower_99 = float("%.3f" % (current_median - range_99))
		upper_99 = float("%.3f" % (current_median + range_99))

		interval_lst.append([compound, lower_99, lower_95, current_median, upper_95, upper_99])

	print 'Done.\n'
	return interval_lst


# Compare randomized confidence intervals and format final data structures
def confidence_interval(score_dict, interval_lst, degree_dict):

	labeled_confidence = []
	sig_count = 0

	for index in interval_lst:
		
		current_compound = index[0]
		current_name = score_dict[current_compound][0]
		current_indegree = degree_dict[current_compound][1]
		current_outdegree = degree_dict[current_compound][2]
		current_score = float(score_dict[current_compound][1])
		
		current_median = float(index[3])
		current_simlower_95conf = float(index[2])
		current_simupper_95conf = float(index[4])
		current_simlower_99conf = float(index[1])
		current_simupper_99conf = float(index[5])
		
		if current_score > current_median:
			current_relationship = 'above'
		elif current_score < current_median:
			current_relationship = 'below'
		else:
			current_relationship = 'none'

		if current_score > current_simupper_95conf:
			if current_score > current_simupper_99conf:
				current_sig = '<0.01'
				sig_count += 1
			else:
				current_sig = '<0.05'
				sig_count += 1
		elif current_score < current_simlower_95conf:
			if current_score < current_simlower_99conf:
				current_sig = '<0.01'
				sig_count += 1
			else:
				current_sig = '<0.05'
				sig_count += 1
		else:
			current_sig = 'n.s.'

		labeled_confidence.append([current_compound, current_name, current_score, current_sig])	

	print('Detected ' + str(sig_count) + ' significant of ' + str(len(interval_lst)) + ' total metabolites.\n')

	return labeled_confidence


# Function to write lists to files	
def write_list(header, out_lst, file_name):

	with open(file_name, 'w') as out_file: 
		
		if not header == 'none': out_file.write(header)
			
		for index in out_lst:
			index = [str(x) for x in index]
			index[-1] = str(index[-1]) + '\n'
			out_file.write('\t'.join(index))

	out_file.close()

# Specialized version of previous function
def write_list_short(header, out_lst, file_name):

	with open(file_name, 'w') as out_file: 
		
		if not header == 'none': out_file.write(header)
			
		for index in out_lst:
			index = [str(x) for x in index]
			index[-1] = str(index[-1]) + '\n'
			out_file.write(''.join(index))

	out_file.close()
			

# Function to write dictionaries to files (next 2 functions are similar)	
def write_dictionary(header, out_dict, file_name):

	all_keys = out_dict.keys()
	
	with open(file_name, 'w') as out_file: 
		
		if not header == 'none': out_file.write(header)
			
		for index in all_keys:
			elements = out_dict[index]
			elements.insert(0, index)
			elements = [str(x) for x in elements]
			elements[-1] = elements[-1] + '\n'
			out_file.write('\t'.join(elements))

	out_file.close()

def write_dictionary_short(header, out_dict, file_name):

	all_keys = out_dict.keys()
	
	with open(file_name, 'w') as out_file: 
		
		if not header == 'none': out_file.write(header)
			
		for index in all_keys:
			entry = index + '\t' + str(out_dict[index]) + '\n'
			out_file.write(entry)

	out_file.close()

def write_dictionary_list(header, out_dict, file_name):

	all_keys = out_dict.keys()
	
	with open(file_name, 'w') as out_file: 
		
		if not header == 'none': out_file.write(header)
			
		for index in all_keys:
			entry = index + '\t' + ','.join(out_dict[index]) + '\n'
			out_file.write(entry)

	out_file.close()


##########################################################################################		
#											 											 #
#									Do The Analysis!								 	 #
#											 											 #
##########################################################################################		


# Citation text
print '''\nbigSMALL v1.0
Released: 12/01/2016

by
Matthew L. Jenior

Department of Microbiology & Immunology
University of Michigan
mljenior@umich.edu

When using, please cite:
Jenior, M.L., Leslie, J.L., Young, V.B., Schloss, P.D. (2016). Clostridium difficile colonizes alternative nutrient niches during infection across distinct murine gut environments. Biorxiv preprint: http://biorxiv.org/content/early/2016/12/07/092304

Distributed under the GNU General Public License\n\n'''

#---------------------------------------------------------------------------------------#		

# Print organism name to screen to track progress in case of loop
if file_name != 'default':
	print '\nImputing metabolism for ' + file_name + '\n'
else:
	current_time = datetime.datetime.now().time()
	current_time = current_time.strftime('%s/%d/%m/%Y')
	current_time = current_time.replace('/','_')
	current_time = current_time.replace('-','')
	file_name = current_time

# Read in and create dictionary for expression
with open(KO_input_file, 'r') as KO_file:
	transcript_dict, total, seq_max = transcription_dictionary(KO_file)
all_KO_lst = transcript_dict.keys()

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
compound_name_dictionary = pickle.load(open(compoundpkl_path, 'rb'))
print('Done.\n')

#---------------------------------------------------------------------------------------#	

# Call translate function and separate output lists
reaction_graph, ko_input_dict, ko_output_dict, compound_lst, KO_lst = network_dictionaries(all_KO_lst, ko_dictionary, reaction_dictionary)

#---------------------------------------------------------------------------------------#	

# Write compounds and enzymes to files
write_list_short('none', compound_lst, 'metabolite.lst')
write_list_short('none', KO_lst, 'enzyme.lst')

# Write network to a two column matrix for use in Neo4j or R
write_list('none', reaction_graph, 'graph.tsv')

#---------------------------------------------------------------------------------------#	

# Calculate actual importance scores for each compound in the network
print 'Calculating compound node connectedness and metabolite scores...\n'
compound_transcript_dict, compound_degree_dict = compile_transcripts(transcript_dict, ko_input_dict, ko_output_dict, compound_lst, KO_lst)
score_dict, degree_dict = calculate_score(compound_transcript_dict, compound_degree_dict, compound_name_dictionary, compound_lst)
print 'Done.\n'

#---------------------------------------------------------------------------------------#		

# Calculate simulated importance values if specified
if iterations >= 1:
	interval_lst = probability_distribution(ko_input_dict, ko_output_dict, degree_dict, KO_lst, compound_name_dictionary, total, seq_max, compound_lst, transcript_dict, iterations)
	final_data = confidence_interval(score_dict, interval_lst, degree_dict)

	# Write all the calculated data to files
	print 'Writing score data with probability distributions to a file...\n'
	outname = 'importances.tsv'
	write_list('Compound_code\tCompound_name\tMetabolite_score\tp_value\n', final_data, outname)
	outname = 'confidence_intervals.tsv'
	write_list('Compound_code\tLower_99_Interval\tLower_95_Interval\tSim_Median\tUpper_95_Interval\tUpper_99_Interval\n', interval_lst, outname)
	print 'Done.\n'

# If simulation not performed, write only scores calculated from measured expression to files	
else:
	print 'Writing score data to a file...\n' 
	outname = 'importances.tsv'
	write_dictionary_short('Compound_code\tCompound_name\tMetabolite_score\n', score_dict, outname)
	print 'Done.\n'

print 'Writing network topology and transcipt counts to files...\n'
outname = 'topology.tsv'
write_dictionary('Compound_code\tCompound_name\tIndegree\tOutdegree\n', degree_dict, outname)
outname = 'KO_mapping.tsv'
write_dictionary_short('KO_code\tTranscripts\n', transcript_dict, outname)
outname = 'input_metabolites.tsv'
write_dictionary_list('KO_code\tCompound_codes\n', ko_input_dict, outname)
outname = 'output_metabolites.tsv'
write_dictionary_list('KO_code\tCompound_codes\n', ko_output_dict, outname)
print 'Done.\n'

#---------------------------------------------------------------------------------------#		

# Wrap everything up

# Report time if iterations are performed
end = time.time()
if end > 10:
	duration = str(int(end - start))
	print '\nCompleted in ' + duration + ' seconds.\n'
else :
	print '\n'
	
print 'Output files located in: ' + directory + '\n\n'		

# Define calculation selection with a string
if iterations > 1:
	iter_str = 'yes'
else:
	iter_str = 'no'

time_unit = 'seconds'
if int(duration) >= 120:
	duration = int(duration) / 60
	time_unit = 'minutes'
if int(duration) >= 120:
	duration = int(duration) / 60
	time_unit = 'hours'

# Write parameters to a file
with open('parameters.txt', 'w') as parameter_file:
	outputString = '''User Defined Parameters
KO expression file: {ko}
Graph name: {name}
KEGG ortholog nodes: {kos}
Substrate nodes: {substrate}
Probability distribution generated: {iter}
Permutations: {perms}
Duration: {time} {tunit}
'''.format(ko=str(KO_input_file), name=str(file_name), iter=iter_str, kos=str(len(KO_lst)), substrate=str(len(compound_lst)), perms=str(iterations), time=str(duration), tunit=time_unit)
	parameter_file.write(outputString)

# Return to the directory the script was called to
os.chdir(starting_directory)	

