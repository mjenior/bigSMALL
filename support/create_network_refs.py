#!/usr/bin/env python3
'''USAGE: python create_network_refs.py kegg_directory
This script creates pickles for all KEGG datasets needed to create genome-scale metabolic network files'''

# On Axiom, KEGG files are located in /mnt/EXT/Schloss-data/kegg/kegg

import sys
import pickle

#---------------------------------------------------------------------------------------#		

# Define where KEGG files are located
kegg = str(sys.argv[1]).rstrip('/')
ko_react = kegg + '/genes/ko/ko_reaction.list'
mapformula = kegg + '/ligand/reaction/reaction_mapformula.lst'
substrate = kegg + '/ligand/compound/compound'

#---------------------------------------------------------------------------------------#		

# Read in and make dictionaries for each reference file needed in the network generation pipeline.


# Create dictionary for converting KO #'s to biochemical reactions #'s

# Reference file excerpt: 
# ko:K00001	rn:R00623
# ko:K00001	rn:R00754
# ko:K00001	rn:R02124
# ko:K00001	rn:R04805
# ko:K00001	rn:R04880
# ko:K00001	rn:R05233
# ko:K00001	rn:R05234
# ko:K00001	rn:R06917
# ko:K00001	rn:R06927
# ko:K00001	rn:R07105
# ko:K00001	rn:R08281

# Output dictionary:
# K00001 = ['R00623', 'R00754', 'R02124', ...]
	
print('\nCreating KO to Reaction dictionary...')

with open(ko_react,'r') as ko:
	
	ko_dict = {}
	
	for line in ko:
		temp = line.split()
		if not temp[0].strip('ko:') in ko_dict.keys():
			ko_dict[temp[0].strip('ko:')] = [temp[1].strip('rn:')]
		else:
			ko_dict[temp[0].strip('ko:')].append(temp[1].strip('rn:'))

print('Complete.\n')	

	
# Create dictionary for converting reaction #'s to input and output compound codes
# Composed of the input and output of each posible compound pairing in a reaction separated by an indicator for if the reaction in reversible.

# Reference file excerpt: 
# R00005: 00330: C01010 => C00011
# R00005: 00791: C01010 => C00011
# R00005: 01100: C01010 <=> C00011
# R00005: 01120: C01010 <=> C00011
# R00006: 00770: C00022 => C00900
# R00008: 00362: C06033 => C00022
# R00008: 00660: C00022 => C06033
# R00008: 01120: C06033 <=> C00022

# Output dictionary:
# R00008 = ['C06033|C00022:N:C06041|C06035', 'C06033|C00022:R:C06041']

print('Creating Reaction to Formula dictionary...')

with open(mapformula,'r') as reaction_map:
	
	reaction_dict = {}

	for line in reaction_map:
	
		temp = line.split(':')
		reactionID = str(temp[0])
		reactionFormula = str(temp[2]).strip()
		
		# Convert directional symbol to a single letter code
		reactionFormula = reactionFormula.replace(' <=> ', ':R:')
		reactionFormula = reactionFormula.replace(' => ', ':N:')
		reactionFormula = reactionFormula.replace(' <= ', ':NR:')
	
		reversibility_reactionFormula = reactionFormula.split(':')[1]
	
		# Check if KEGG reaction is listed right to left
		if reversibility_reactionFormula == 'NR':
			input_reactionFormula = reactionFormula.split(':')[2]
			output_reactionFormula = reactionFormula.split(':')[0]
			reversibility_reactionFormula = 'N'
		
		else:
			input_reactionFormula = reactionFormula.split(':')[0]
			output_reactionFormula = reactionFormula.split(':')[2]
	
		inputString_reactionFormula = input_reactionFormula.replace(' + ', '|')
		outputString_reactionFormula = output_reactionFormula.replace(' + ', '|')
	
	
		# Put all the formula information together in one string
		completeString_reactionFormula = inputString_reactionFormula + ':' + reversibility_reactionFormula + ':' + outputString_reactionFormula
		
		# Add entries to the dictionary, and append to those which already exist
		if not reactionID in reaction_dict.keys():
			reaction_dict[reactionID] = [completeString_reactionFormula]
		else:
			reaction_dict[reactionID].append(completeString_reactionFormula)


# Create essentially the same dictionary as above, but without reversibility information
with open(mapformula,'r') as reaction_map_nr:
	
	reaction_dict_nonrev = {}

	for line in reaction_map_nr:
	
		temp = line.split(':')
		reactionID = str(temp[0])
		reactionFormula = str(temp[2]).strip()
		
		# Convert directional symbol to a single letter code
		reactionFormula = reactionFormula.replace(' <=> ', ':N:')
		reactionFormula = reactionFormula.replace(' => ', ':N:')
		reactionFormula = reactionFormula.replace(' <= ', ':NR:')
	
		reversibility_reactionFormula = reactionFormula.split(':')[1]
	
		# Check if KEGG reaction is listed right to left
		if reversibility_reactionFormula == 'NR':
			input_reactionFormula = reactionFormula.split(':')[2]
			output_reactionFormula = reactionFormula.split(':')[0]
			reversibility_reactionFormula = 'N'
		
		else:
			input_reactionFormula = reactionFormula.split(':')[0]
			output_reactionFormula = reactionFormula.split(':')[2]
	
		inputString_reactionFormula = input_reactionFormula.replace(' + ', '|')
		outputString_reactionFormula = output_reactionFormula.replace(' + ', '|')
	
		# Put all the formula information together in one string
		completeString_reactionFormula = inputString_reactionFormula + ':' + reversibility_reactionFormula + ':' + outputString_reactionFormula
		
		# Add entries to the dictionary, and append to those which already exist
		if not reactionID in reaction_dict_nonrev.keys():
			reaction_dict_nonrev[reactionID] = [completeString_reactionFormula]
		else:
			reaction_dict_nonrev[reactionID].append(completeString_reactionFormula)

print('Complete.\n')


print('Creating Compound Code to Chemical Name dictionary...')

with open(substrate,'r') as compounds:
	
	compound_dict = {}
	
	for line in compounds:
		
		line = line.split()
		
		if line[0] == 'ENTRY':
			entry = line[1]
			continue
			
		elif line[0] == 'NAME':	
			name = '_'.join(line[1:]).rstrip(';')
			if not entry in compound_dict.keys():
				compound_dict[entry] = name


		# Need to add substrate classifications to names, separated by ;
		
print('Complete.\n')

#---------------------------------------------------------------------------------------#		

# Write a pickle for each of the dictionaries created above.
print('Writing dictionaries to pkl files...')

# Creates a pickle of the dictionary for translating KOs into lists of biochemical reactions.
with open('ko_reaction.pkl','wb') as outfile:
	pickle.dump(ko_dict, outfile)
	
# Creates a pickle of the dictionary for translating biochemical reactions into lists of input and output compounds.
with open('reaction_mapformula.pkl','wb') as outfile:
	pickle.dump(reaction_dict, outfile)

with open('reaction_mapformula_nonrev.pkl','wb') as outfile:
	pickle.dump(reaction_dict_nonrev, outfile)	

# Creates a pickle of the dictionary for KEGG compound codes to their common names.
with open('compound.pkl','wb') as outfile:
	pickle.dump(compound_dict, outfile)
	
print('Complete.\n')

