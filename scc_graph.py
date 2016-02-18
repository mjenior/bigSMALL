#!/usr/bin/env python
'''USAGE: python scc_graph.py KOfile --giant False --component 5 --confidence 0.2
Calculates SCC genome-scale metabolic networks based on KOs and derives lists of 
seed and sink compounds (along with several other files listed below).
'''

# Written by Matthew Jenior, University of Michigan, Schloss Laboratory, 2014-2015
# Compatible with python 2.7

# Graph traversal and seed calculation algorithms adopted from NetSeed python by Rogan Carr, 2014
	# Using python version of SCC graph travesal and seed calculation code

# The function of this code is to convert lists of genes to unweighted, directed graphs,
	# and use this data structure to impute information about the overall metabolism of the organism(s).

# Dependencies:  first argument is a text file containing a single column of KO codes 
	# that are to be included in the metabolic network.  Additional arguments to 
	# specify properties of the output files are as follows:
		# --giant, default=False, Only focus on the largest component or not
		# --component, default=5, Minimum threshold for SCC size to be included
		# --confidence, default=0.2, Minimum confidence value allowed for seeds

# Generate files:  A new directory in ./ ending in ".scc.files" that contains all other output including:
	# A 2 column network file of input and output compounds forming the given metabolic network after KO translation
	# A text file containing reference errors thrown during the translation of KOs to chemical equations
	# A "seeds" directory contained files with information about strongly connected components from the graph
		# nodes.txt - a list of all nodes included in the final network
		# nonseeds_names.tsv - Chemical names of all non-seed compounds
		# nonseeds.tsv - KEGG code for non-seed compounds
		# nonsinks_names.tsv - Chemical names of all non-sink compounds
		# nonsinks.tsv - KEGG code for non-sink compounds
		# prunedNetwork.tsv - New network file containing only those nodes in the SCC network
		# pruned.txt - List of nodes from original network that were not strongly connected
		# seedgroups.tsv - Lists of seeds with those that are included in cyclic components together
		# seedParameters.txt - User defined parameters for stringency of seed calculation
		# seeds_names.tsv - Chemical names of all seed compounds
		# seeds.tsv - KEGG code for seed compounds
		# sinkgroups.tsv - Lists of sinks with those that are included in cyclic components together
		# sinks_names.tsv - Chemical names of all sink compounds
		# sinks.tsv - KEGG code for sink compounds
		
#---------------------------------------------------------------------------------------#		

# Import all the necessary python modules
from __future__ import print_function
import sys
import os
import pickle
import argparse
import imp
import itertools

starting_directory = str(os.getcwd())
# Get the correct path for the reference files
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Function for reading a tab-separated file of edges
def readGraph(graphFile):

	graph={}
	readGraph = open(graphFile, 'r')
	
	for line in readGraph:
		vals = line.strip().split("\t");
		
		if len(vals[0]) == 0 or vals[0][0] == '#':
			continue
	
		elif len(vals) != 2:
			print("Bad line: " + line)
			continue
			
		graph.setdefault(vals[0],[]).append(vals[1])
        
	return graph
	
# Function for creating seed output files from dictionaries
def dictPrint(seedData, fileName):
	
	seedFile = open(fileName,'w')
	outputList = []
	
	for key in seedData.keys():
		seedFile.write(str(key) + "\t" + str(seedData[key]))
		outputList.append(str(key))
		seedFile.write('\n')

	seedFile.close()
	return outputList

# Function to write group files
def dictPrintGroup(seedData, fileName):
	
	seedFile = open(fileName,'w')
	
	for key in seedData.keys():
		newLine = '\t'.join(seedData[key])
		seedFile.write(newLine)
		seedFile.write('\n')

	seedFile.close()

# Function for creating seed output files from sets/lists
def textPrint(seedData, fileName):
	
	seedFile = open(fileName, 'w')
	
	for x in list(seedData):
		seedFile.write(x)
		seedFile.write('\n')

	seedFile.close()
	

# Function to created pruned network file
def prunedNetFile(pruned, original, newName, seeds, sinks):
		
	prunedFile = open(newName, 'w')
	inputList = []
	outputList = []
	
	for index in original:
	
		temp = index.split()
		
		if not temp[0] in pruned and not temp[1] in pruned:
			prunedFile.write(index)
	
	original.close()	
	prunedFile.close()
				

#---------------------------------------------------------------------------------------#		

# Functions for graph traversal and SCC identification. Adapted NetSeed python by Rogan Carr PhD., 2014
# calculate_seeds()
# depth_first_search()
# dfs_get_neighbors()

# Calculate the seed set
def calculate_seeds(graph,onlyGiant=None,minComponentSize=None,seedThreshold=None):

    #  Construct the reverse graph
    reverse_graph={}

    # Calculate node order
    nodes=set()
    nodeOrder=[];
    
    for source_node in graph:

        # Add it to the list of nodes
        if not source_node in nodes:
            nodes.add(source_node)
            nodeOrder.append(source_node)

        for target_node in graph[source_node]:
            if not target_node in reverse_graph.keys():
                reverse_graph[target_node]=[source_node]
            else:
                reverse_graph[target_node].append(source_node)

            # Add it to the list of nodes
            if not target_node in nodes:
                nodes.add(target_node)
                nodeOrder.append(target_node)
    
    
    # Return all the nodes that were not pruned    
    Pruned=set()
    # Trim (all or some) Small components
    if onlyGiant or minComponentSize:
        
        # Use DFS to determine the connected components
        nComponents,nodeToComponent = depth_first_search(graph,nodeOrder,reverse_graph)[1:3]
        
        # Calculate the size of each component
        componentSize = [0]*nComponents        

        for node in nodeOrder:
            componentSize[nodeToComponent[node]]+=1

        # If we only keep the giant component
        #  Then set the minimum size to the giant c's size
        if onlyGiant:
            minComponentSize=0
            for c in xrange(nComponents):
                if componentSize[c] > minComponentSize:
                    minComponentSize=componentSize[c]
        
        # Delete all nodes not in the minimal component

        #  Identify the components to trim
        componentsToTrim = set()
        for c in xrange(nComponents):
            if componentSize[c] < minComponentSize:
                componentsToTrim.add(c)
        
        #  Trim the nodes
        #    Note the .keys() function is necessary because we are deleting keys!
        #    Also, we add to Pruned, remove frome nodeOrder
        for node in graph.keys():      		
            if nodeToComponent[node] in componentsToTrim:
                del graph[node]

                if node in reverse_graph:
                    del reverse_graph[node]
                Pruned.add(node)

                if node in nodeOrder:
                    nodeOrder.remove(node)
        
        #  Capture the only-pointed-to nodes
        #    Note the .keys() function is necessary because we are deleting keys!
        for node in reverse_graph.keys():
            if nodeToComponent[node] in componentsToTrim:
                del reverse_graph[node]
                Pruned.add(node)
                if node in nodeOrder:
                    nodeOrder.remove(node)

    # Identify the Strongly-Connected Components (SCC)
    finishOrder,nComponents = depth_first_search(graph,nodeOrder)[0:2]
    nSCC,nodeToSCC = depth_first_search(reverse_graph,finishOrder)[1:3]
    
    #  Calculate the size of each scc
    sccSize = [0]*nSCC
    for node in nodeOrder:
        sccSize[nodeToSCC[node]]+=1
    
    
    # Determine which SCCs are sources
    #   Initialize to all SCCs
    SeedGroups={}
    SinkGroups={}
    for i in xrange(nSCC):
        SeedGroups[i]=[]
        SinkGroups[i]=[]
    #   Remove SCCs that are targets of others
    for source_node in graph:
        for target_node in graph[source_node]:
            if nodeToSCC[target_node] != nodeToSCC[source_node]:
                if nodeToSCC[target_node] in SeedGroups:
                    del SeedGroups[nodeToSCC[target_node]]
                
                if nodeToSCC[source_node] in SinkGroups:
                	del SinkGroups[nodeToSCC[source_node]]
					
	
    # Assemble the seeds
    Seeds={}
    nonSeeds={}
    for node in nodeOrder:
        seed_rank=0 
        
        if nodeToSCC[node] in SeedGroups:
            seed_rank = 1 / float(sccSize[nodeToSCC[node]])
            
        if seed_rank > seedThreshold:
            Seeds[node]=seed_rank
            SeedGroups[nodeToSCC[node]].append(node)
        else:
            # This can be nonzero
            nonSeeds[node]=seed_rank # need to create a statment to omit sinks from nonseed

    # Remove seed groups that are defunct due to their size
    #    Note the .keys() function is necessary because we are deleting keys!
    for group in SeedGroups.keys():   													
        if len(SeedGroups[group]) < sccSize[group]:
            if len(SeedGroups[group]) == 0 and 1/float(sccSize[group]) <= seedThreshold:
                print("Removing Seed Group", group, "of size", sccSize[group])
                del SeedGroups[group]
            else:
                print("Serious problem in SeedGroup", group, ": Expected", sccSize[group], "and found", len(SeedGroups[group]),".")

	
	# Assemble the sinks
    Sinks={}
    nonSinks={}
    for node in nodeOrder:
        sink_rank=0
        if nodeToSCC[node] in SinkGroups:
            sink_rank = 1 / float(sccSize[nodeToSCC[node]])
            
        if sink_rank > seedThreshold:
            Sinks[node]=sink_rank
            SinkGroups[nodeToSCC[node]].append(node)
        else:
            # This can be nonzero
            nonSinks[node]=sink_rank
            
    # Remove sink groups that are defunct due to their size
    for group in SinkGroups.keys():
        if len(SinkGroups[group]) < sccSize[group]:
            if len(SinkGroups[group]) == 0 and 1/float(sccSize[group]) <= seedThreshold:
                print("Removing Sink Group", group, "of size", sccSize[group])
                del SinkGroups[group]
            else:
                print("Serious problem in SinkGroup", group, ": Expected", sccSize[group], "and found", len(SinkGroups[group]),".")
     
            
    return [Seeds, SeedGroups, nonSeeds, Pruned, nodes, Sinks, SinkGroups, nonSinks]
        
# Function for depth-first graph traversal
def depth_first_search(graph,nodeOrder,reverse_graph=[]):
    
    # To return
    component=-1
    nodeToComponent={}
    finishOrder=[]
    
    # Our search space
    stack=[]

    # Working
    visited={}
    
    # Search every node
    for node in reversed(nodeOrder):
        if node in visited:
            continue
        
        # We must be on a new component if we haven't been here
        component+=1

        # Put it on our search space
        stack.append(node)

        # We visit each node twice
        while stack:
            vertex = stack.pop()

            # If it's visited, roll it up
            if vertex in visited:
                if visited[vertex]<2:
                    finishOrder.append(vertex)
                    visited[vertex]=2
            else:
                # Assign the component
                nodeToComponent[vertex] = component
                # Put it on our search space
                stack.append(vertex)
                
                # Mark it as visited
                visited[vertex]=1

                dfs_get_neighbors(graph,visited,vertex,stack)
                # Get the neighobrs for the revers graph if it's defined
                if reverse_graph:
                    dfs_get_neighbors(reverse_graph,visited,vertex,stack)

    # convert to a 1-based list
    component+=1
    return [finishOrder, component, nodeToComponent]

# Function to identify neighbors in SCC network
def dfs_get_neighbors(graph,visited,node,stack):
    if node in graph:
        for next_node in graph[node]:
            if next_node not in visited:
                stack.append(next_node)
                
#---------------------------------------------------------------------------------------#		
	
# Set the parameters you find important

parser = argparse.ArgumentParser(description='Calculate set seeds from KO lists.')
parser.add_argument('newfile')
parser.add_argument('--giant', default=False, help='Only focus on the largest component or not')
parser.add_argument('--component', default=5, help='Minimum threshold for SCC size to be included')
parser.add_argument('--confidence', default=0.2, help='Minimum confidence value allowed for seeds')
args = parser.parse_args()

#---------------------------------------------------------------------------------------#		

# Load in dictionaries for converting KOs to compound graph

# Read in pickled glycan list - DATABASE HAS IMPROVED, NO LONGER INCLUDED
#glycanpkl_path = script_path + '/support/glycan.pkl'
#glycanpkl_file = open(glycanpkl_path, 'rb')
#glycan_list = pickle.load(glycanpkl_file)
#glycanpkl_file.close()

# Read in pickled KO to reaction dictionary
ko_reactionpkl_path = script_path + '/support/ko_reaction.pkl'
ko_reactionpkl_file = open(ko_reactionpkl_path, 'rb')
ko_dict = pickle.load(ko_reactionpkl_file)
ko_reactionpkl_file.close()

# Read in pickled reaction to reaction_mapformula dictionary
reaction_mapformulapkl_path = script_path + '/support/reaction_mapformula.pkl'
#reaction_mapformulapkl_path = script_path + '/support/reaction_mapformula_nonrev.pkl'
reaction_mapformulapkl_file = open(reaction_mapformulapkl_path, 'rb')
reaction_dict = pickle.load(reaction_mapformulapkl_file)
reaction_mapformulapkl_file.close()

# Read in pickled reaction to compound dictionary
compoundpkl_path = script_path + '/support/compound.pkl'
compoundpkl_file = open(compoundpkl_path, 'rb')
compound_dict = pickle.load(compoundpkl_file)
compoundpkl_file.close()

#---------------------------------------------------------------------------------------#

# Set up files and directories

# Take in input KO file
infile = open(args.newfile,'r')
outfile_name = str(args.newfile).split('/')[-1]
    
# Navigate to KO file directory in case your not already there
new_directory = os.path.split(os.path.abspath(outfile_name))[0] + '/'
os.chdir(new_directory)

# Create a new directory for output files and navigate to it
directory = str(os.getcwd()) + '/' + outfile_name + '.scc.files'
if not os.path.exists(directory):	
	os.makedirs(directory)
os.chdir(directory)
print('\nOutput files located in: ' + directory)

# Create file for compound graph
outfile_name = 'network.' + outfile_name
outfile = open(outfile_name,'w')

#---------------------------------------------------------------------------------------#		

# Translate organism's KO list into an input and output compound graph

# Controls if complex glycan containing reactions are omitted from the graph (0 or 1)
# Since the KEGG glycan database annotations have improved since Borenstein et al. 2008, not omitted by default
glycan_switch = 0 

# Create file for reporting key errors
errorfile = open('RefErrorLog.txt', 'w')

excludedCountKO = 0
triedCountKO = 0
excludedCountReact = 0
triedCountReact = 0
totalIncludedReact = 0
glycanExcluded = 0

# Nested loops to finally convert the KO list to a directed graph of input and output compounds	
for line in infile:
	current_ko = str(line.split()[-1]).strip('ko:')
	triedCountKO += 1
	
	try:
		reaction_number = ko_dict[current_ko]
	except KeyError:
		errorString = 'WARNING: ' + str(current_ko) + ' not found in KO dictionary.'
		errorfile.write(errorString)
		errorfile.write('\n')
		excludedCountKO += 1
		continue
		
	for index in reaction_number:
		marked = 0
		triedCountReact += 1
		try:
			reaction_collection = reaction_dict[index]
		except KeyError:
			errorString = 'WARNING: ' + str(index) + ' not found in reaction dictionary.'
			errorfile.write(errorString)
			errorfile.write('\n')
			excludedCountReact += 1
			continue
			
		# Filter out those reactions that contain compounds from Glycans KEGG database
		if glycan_switch != 0:
			for x in reaction_collection:
				reaction_info = x.split(':')	
			
				if str(reaction_info[0]) in glycan_list or str(reaction_info[2]) in glycan_list:
					marked = 1
					glycanExcluded += 1
				else: continue
		
		if marked == 1: continue
		
		for x in reaction_collection:
			totalIncludedReact += 1
			# Spit reaction input and output as well as the list of compounds with each
			reaction_info = x.split(':')
			input_compounds = reaction_info[0].split('|')
			output_compounds = reaction_info[2].split('|')
			
			# Create both the the forward and reverse combinations of compounds
			F_compound_list = [input_compounds, output_compounds]
			R_compound_list = [output_compounds, input_compounds]
			
			# Calculate all combinations of forward and reverse reactions using itertools
			F_combination_list = list(itertools.product(*F_compound_list))
			R_combination_list = list(itertools.product(*R_compound_list))
			
			if str(reaction_info[1]) == 'R':
			
				for index in F_combination_list:
					outfile.write('\t'.join([str(index[0]), str(index[1])]))
					outfile.write('\n')
					
				for index in R_combination_list:
					outfile.write('\t'.join([str(index[0]), str(index[1])]))
					outfile.write('\n')
					
			elif str(reaction_info[1]) == 'N':
				for index in F_combination_list:
					outfile.write('\t'.join([str(index[0]), str(index[1])]))
					outfile.write('\n')	

infile.close()
outfile.close()
errorfile.close()

# Calculate and print translation success
includedCountKO = triedCountKO - excludedCountKO
print('\nOmitting ' + str(excludedCountKO) + ' untranslatable KEGG orthologs')

excludedAllReact = excludedCountReact + glycanExcluded
includedCountReact =  triedCountReact - excludedAllReact
print('Omitting ' + str(excludedAllReact) + ' untranslatable Reaction groups\n')

#---------------------------------------------------------------------------------------#		

# Use the directed compound graphs to calculate seeds and non-seeds

# Read in the files and get the properly formatted graph
graph = readGraph(outfile_name)

# Calculate the seed information
onlyGiant = args.giant
minComponentSize = args.component
seedThreshold = args.confidence
Seeds, SeedGroups, nonSeeds, Pruned, Nodes, Sinks, SinkGroups, nonSinks = calculate_seeds(graph, onlyGiant, minComponentSize, seedThreshold)

#---------------------------------------------------------------------------------------#		

# Write ReadMe file with explanations of output files in the new scc.files directory
readmefile = open('ReadMe.txt', 'w')
readme_str = '''nodes.txt - a list of all nodes included in the final network
nonseeds_names.tsv - Chemical names of all non-seed compounds
nonseeds.tsv - KEGG code for non-seed compounds, with associated confidence values
nonsinks_names.tsv - Chemical names of all non-sink compounds, with associated confidence values
nonsinks.tsv - KEGG code for non-sink compounds, with associated confidence values
prunedNetwork.tsv - New network file containing only those nodes in the SCC network
pruned.txt - List of nodes from original network that were not strongly connected
seedgroups.tsv - Lists of seeds with those that are included in cyclic components together
seedParameters.txt - User defined parameters for stringency of seed calculation
seeds_names.tsv - Chemical names of all seed compounds, with associated confidence values
seeds.tsv - KEGG code for seed compounds, with associated confidence values
sinkgroups.tsv - Lists of sinks with those that are included in cyclic components together
sinks_names.tsv - Chemical names of all sink compounds, with associated confidence values
sinks.tsv - KEGG code for sink compounds, with associated confidence values
intermediates_names.tsv - Chemical names of all intermediate compounds, with associated confidence values
intermediates.tsv - KEGG code for intermediate compounds, with associated confidence values'''
readmefile.write(readme_str)
readmefile.close()

# Create and navigate to a new directory for seed information
directory = str(os.getcwd()) + '/seeds'
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

# Write seed parameters to a file
parameterFile = open('seedParameters.txt', 'w')

outputString = '''Seed Calculation Parameters

Input Graph Name: {name}
Only Giant Component: {giant}
Minimum Component Size: {component}
Seed Threshold: {threshold}

-------------------------------------

Translation Results

Reactions attempted during KO translation: {tryko}
Reactions excluded during KO translation: {exko}
Reactions included during KO translation: {inko}

Compound combinations attempted during reaction translation: {tryrxn}
Compound combinations excluded during reaction translation: {exrxn}
Compound combinations included during reaction translation: {inrxn}

-------------------------------------

Output Information

Total Nodes: {node}
Total Pruned Nodes: {prune}
Total Seeds: {seed}
Total Seed Groups: {seedgroups}
Total Non-seeds: {nonseed}
Total Sinks: {sink}
Total Sink Groups: {sinkgroups}
Total Non-sinks: {nonsink}
'''.format(name=str(outfile_name), giant=str(onlyGiant), component=str(minComponentSize), threshold=str(seedThreshold), tryko=str(triedCountKO), exko=str(excludedCountKO), inko=str(includedCountKO), tryrxn=str(triedCountReact), exrxn=str(excludedCountReact), inrxn=str(includedCountReact), node=str(len(Nodes)), prune=str(len(Pruned)), seed=str(len(Seeds)), seedgroups=str(len(SeedGroups)), nonseed=str(len(nonSeeds)), sink=str(len(Sinks)), sinkgroups=str(len(SinkGroups)), nonsink=str(len(nonSinks)))

parameterFile.write(outputString)
parameterFile.close()


# Print most important seed information
print('\nPruning ' + str(len(Pruned)) + ' nodes\n')
print('Total Seeds: ' + str(len(Seeds)))
print('Total Sinks: ' + str(len(Sinks)))
print('Total Intermediates: ' + str(len(nonSeeds)-len(Sinks)))
print('Total SCC Nodes: ' + str(len(nonSeeds)+len(Seeds)))


seedsOut = 'seeds.tsv'
seedgroupsOut = 'seedgroups.tsv'
nonseedsOut = 'nonseeds.tsv'
prunedOut = 'pruned.txt'
nodesOut = 'nodes.txt'
sinksOut = 'sinks.tsv'
sinkgroupsOut = 'sinkgroups.tsv'
nonsinksOut = 'nonsinks.tsv'
prunedNet = 'prunedNetwork.tsv'

seed_list = dictPrint(Seeds, seedsOut)
dictPrintGroup(SeedGroups, seedgroupsOut)
nonseed_list = dictPrint(nonSeeds, nonseedsOut)
textPrint(Pruned, prunedOut)
textPrint(Nodes, nodesOut)
sink_list = dictPrint(Sinks, sinksOut)
dictPrintGroup(SinkGroups, sinkgroupsOut)
nonsink_list = dictPrint(nonSinks, nonsinksOut)
intermediatesFile = open('intermediates.tsv','w')
names_intermediatesFile = open('intermediates_names.tsv','w')

final_outfile_name = '../' + outfile_name
full_network = open(final_outfile_name, 'r')
pruned_list = list(Pruned)

# Add seed and sink lists
prunedNetFile(pruned_list, full_network, prunedNet, seed_list, sink_list)

#---------------------------------------------------------------------------------------#		

# Translate seeds, non-seeds, sinks, and non-sinks to compounds names

seedFile = open(seedsOut,'r')
nonseedFile = open(nonseedsOut,'r')
names_seedsOut = seedsOut.replace('.tsv', '_names.tsv')
names_nonseedsOut = nonseedsOut.replace('.tsv', '_names.tsv')
names_seedFile = open(names_seedsOut,'w')
names_nonseedFile = open(names_nonseedsOut,'w')
sinkFile = open(sinksOut,'r')
nonsinkFile = open(nonsinksOut,'r')
names_sinksOut = sinksOut.replace('.tsv', '_names.tsv')
names_nonsinksOut = nonsinksOut.replace('.tsv', '_names.tsv')
names_sinkFile = open(names_sinksOut,'w')
names_nonsinkFile = open(names_nonsinksOut,'w')

for index in seedFile:
	entry = '\t'.join([compound_dict[str(index.split()[0])][0], str(index.split()[0]), str(index.split()[1]) + '\n'])
	names_seedFile.write(entry)

temp_nonseed_list = []
nonseed_confidence_dict = {}
for index in nonseedFile:
	entry = '\t'.join([compound_dict[str(index.split()[0])][0], str(index.split()[0]), str(index.split()[1]) + '\n'])
	temp_nonseed_list.append(str(index.split()[0]))
	names_nonseedFile.write(entry)
	nonseed_confidence_dict[str(index.split()[0])]=str(index.split()[1])

temp_sink_list = []	
for index in sinkFile:
	entry = '\t'.join([compound_dict[str(index.split()[0])][0], str(index.split()[0]), str(index.split()[1]) + '\n'])
	temp_sink_list.append(str(index.split()[0]))
	names_sinkFile.write(entry)

for index in nonsinkFile:
	entry = '\t'.join([compound_dict[str(index.split()[0])][0], str(index.split()[0]), str(index.split()[1]) + '\n'])
	names_nonsinkFile.write(entry)	

intermediates = list(set(temp_nonseed_list) - set(temp_sink_list))
for index in intermediates:
	entry = str(index) + '\t' + str(nonseed_confidence_dict[index]) + '\n'
	intermediatesFile.write(entry)
	entry = '\t'.join([str(compound_dict[index][0]), str(index), str(nonseed_confidence_dict[index]) + '\n'])
	names_intermediatesFile.write(entry)	
	
seedFile.close()
nonseedFile.close()
names_seedFile.close()
names_nonseedFile.close()
sinkFile.close()
nonsinkFile.close()
names_sinkFile.close()
names_nonsinkFile.close()
intermediatesFile.close()
names_intermediatesFile.close()

os.chdir(starting_directory)

print('\nDone.\n')
