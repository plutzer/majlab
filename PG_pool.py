#PG_pool script

import numpy as np
import pandas as pd
import networkx
import sys

def simplify_protein_IDs(dataframe):
    '''
        Takes a proteinGroups dataframe and returns a proteinGroups dataframe with simplified protein IDs.
        dataframe: pandas Dataframe with a 'Protein IDs' column
        Output: pandas Dataframe with the 'Protein IDs' column simplified
    '''
    dataframe['Protein IDs'] = [protein_group.split(';')[0] for protein_group in dataframe['Protein IDs']]
    return dataframe

def merge_pgs(dataframe1,dataframe2 = None):
    '''
        Merges two proteinGroups files on simplified Protein IDs, and rectifies other columns needed for downstream analysis.
            If no second dataframe is given, the function will perform filtering and column simplification steps on just the 
            first dataframe. Otherwise, the first dataframe is assumed to already be rectified and the filtering/simplification
            will only be done on the second dataframe.
        dataframe1: pandas dataframe for first proteinGroups file
        dataframe2: pandas dataframe for second proteinGroups file. Default is null if just the first dataframe is being filtered. 
        Output: pandas dataframe for the merged proteinGroups file, or filtered dataframe.
    '''
    # If no second dataframe is given, just filter and fix the Protein IDs column of the first dataframe.
    if dataframe2 is None:
        # Use the simplify function to fix the Protein IDs column.
        dataframe1 = simplify_protein_IDs(dataframe1)
        #Filter out any rows where the protein group was identified only by type
        dataframe1 = dataframe1[dataframe1['Score'] > 0]
        # Return the filtered and fixed dataframe
        return dataframe1
    # If a second dataframe is given, it is assumed that the first dataframe is already filtered and rectified.
    else:
        # Simplify the Protein IDs column.
        dataframe2 = simplify_protein_IDs(dataframe2)
        #Filter out any rows where the protein group was identified only by type
        dataframe2 = dataframe2[dataframe2['Score'] > 0]
        
        # Merge the columns unique to the second dataframe with the first dataframe using protein IDs
        combined_df = dataframe1.merge(dataframe2[list(np.setdiff1d(dataframe2.columns,dataframe1.columns))+['Protein IDs']],how='outer',on='Protein IDs')
        
        ##### Next fix the columns of the new dataframe. I'm only fixing the columns that matter for downstream SAINT analysis.
        
        # Fix Score (max of scores for a given protein group)
        max_scores = []
        for group in combined_df['Protein IDs']:
            max_scores.append(np.max(list(dataframe1[dataframe1['Protein IDs'] == group]['Score'])+list(dataframe2[dataframe2['Protein IDs'] == group]['Score'])))
        combined_df['Score'] = max_scores
        
        # Fix Reverse (+ if either dataframe has a +)
        reverse = []
        for group in combined_df['Protein IDs']:
            if ('+' in list(dataframe1[dataframe1['Protein IDs'] == group]['Reverse'])) or ('+' in list(dataframe2[dataframe2['Protein IDs'] == group]['Reverse'])):
                reverse.append('+')
            else:
                reverse.append('')
        combined_df['Reverse'] = reverse
        
        # Fix Contaminant (+ if either dataframe has a +)
        contaminant = []
        for group in combined_df['Protein IDs']:
            if ('+' in list(dataframe1[dataframe1['Protein IDs'] == group]['Potential contaminant'])) or ('+' in list(dataframe2[dataframe2['Protein IDs'] == group]['Potential contaminant'])):
                contaminant.append('+')
            else:
                contaminant.append('')
                
        # Fix Sequence length (should be the same - still use max just in case)
        max_seqlength = []
        for group in combined_df['Protein IDs']:
            max_seqlength.append(np.max(list(dataframe1[dataframe1['Protein IDs'] == group]['Sequence length'])+list(dataframe2[dataframe2['Protein IDs'] == group]['Sequence length'])))
        combined_df['Sequence length'] = max_seqlength
        # Fix Identified by site (+ if either dataframe has a + ... just in case: shouldn't matter because I'm filtering these out to start)
        combined_df['Only identified by site'] = ['' for item in range(len(combined_df))]
        
        # Return the combined dataframe
        return combined_df

def generate_linked_groups(groupedIDs):
    '''
        Takes a list of protein groups and generates linked protein groups based on shared proteins.
        Generates an adjacency matrix for protein groups where:
            0 = no shared proteins between groups
            1 = shared proteins between groups
        And then finds connected subgraphs for that network. Each subgraph is a set of linked groups.

        groupedIDs: list of strings - protein IDs column from protein groups file. Proteins within a group should be separated by ';'
        returns: list of lists of ints - each element is a list of indexes protein groups in a linked group
    '''
    adjMat = np.zeros((len(groupedIDs),len(groupedIDs))) #Initializing the adjacency matrix
    
    #Iterate over each protein group, and find other protein groups that contain the same proteins.
    for group_index in range(len(groupedIDs)):
        for protein in list(groupedIDs)[group_index].split(';'):
            
            # also in is a list of indexes of protein groups that share a protein with this group.
            also_in  = list(np.where([protein in group.split(';') for group in groupedIDs])[0])
            
            # Iterate through also_in and update the adjacency matrix
            for index in also_in:
                if group_index != index:
                    adjMat[group_index,index] = 1
    
    # Generate a list of indexes for component subgraphs
    graph = networkx.from_numpy_matrix(adjMat)
    return [list(graph.subgraph(c).copy().nodes()) for c in networkx.connected_components(graph)]

def add_FDR(mergedPG,colname,func=np.max):
	'''
		This function calculates FDR using the decoy approach for linkages of protein groups and
		adds it the the mergedPG dataframe. The FDR is calculated for linkages of protein groups
		that share some peptide with other groups in a given column.

		mergedPG: a proteinGroups pandas dataframe
		colname: the name of the column in mergedPG to use for generating linkage groups
				For example: 'Protein IDs', 'Majority protein IDs'
		func: function to use for consolidating scores for linked protein groups. Max is
				default because high evidence for any protein group in a linkage is poteintially
				grounds to maintain the entire linkage.

	'''
	# Generate linked groups on the Protein IDs column (this shouldn't generate any linkages)
	linked_groups = generate_linked_groups(mergedPG[colname])

	# Get the linked scores (again, these should just be the same scores as in the original dataframe)
	linked_scores = [func(list(mergedPG.iloc[indexes]['Score'])) for indexes in linked_groups]

	# Denote whether each linked group is a decoy. This is a boolean list, same indexes as the other 2 lists. Use Reverse column...
	linked_decoys = []
	for indexes in linked_groups:
	    boolval = False
	    indexes_subset = mergedPG.iloc[indexes]
	    if '+' in list(indexes_subset['Reverse']): #sum(indexes_subset['Reverse']) != 0:
	        boolval = True
	    linked_decoys.append(boolval)
	
	# Make a dataframe out of the lists
	linkages_info = pd.DataFrame(zip(linked_groups,linked_scores,linked_decoys),columns = ['groups','scores','decoys'])

	# Initialize the FDRs list
	FDRs = list(np.zeros([len(mergedPG)]))

	# Compute FDR using the decoy approach
	linkage_FDRs = []
	for i in range(len(linked_groups)): #Iterating over the indexes of the linked groups list itself
	    # subset linkages_info to only elements that have a score less than than the current linked_group
	    subset = linkages_info[linkages_info['scores'] >= linked_scores[i]]
	    if len(subset) < 1:
	    	FDR = 0
	    else:
	    	FDR = np.sum((subset['decoys']))/(len(subset)-np.sum((subset['decoys'])))
	    linkage_FDRs.append(FDR)
	    for j in linked_groups[i]:
	    	FDRs[j] = FDR

	linkages_info['FDR'] = linkage_FDRs
	#linkages_info.to_csv('linkages_info.csv') # in case you need to debug

	#Add the FDRs column to the datarame.
	mergedPG['FDR'] = FDRs

	return mergedPG


####################################################
################ MAIN SCRIPT #######################
####################################################

#Read in the experimental design file from the command line
expt_des_filename = sys.argv[1]
expt_design = pd.read_csv(expt_des_filename)

#Get list of paths - column needs to be named 'PG paths'
pathslist = list(expt_design['PG paths'].unique())

# Use the merge_pgs function to combine together all of the paths
merged_PG = merge_pgs(pd.read_csv(pathslist[0],delimiter='\t'))
for i in range(1,len(pathslist)):
	merged_PG = merge_pgs(merged_PG,pd.read_csv(pathslist[i],delimiter='\t'))

# If there's no Gene names column in the protein groups, add it using the Protein IDs
# TO DO: what if one file has this column but others don't?
names = False
#if not ('Gene names' in merged_PG.columns):
#	names = False
if not names:
	merged_PG['Gene names'] = [item for item in merged_PG['Protein IDs']]

# Fill in NaNs in the MS/MS count columns
for col in merged_PG.columns:
	if 'MS/MS count' in col:
		merged_PG[col] = merged_PG[col].fillna(0)

# Compute the FDR and add it to the dataframe
merged_PG = add_FDR(merged_PG,'Protein IDs')

# Apply a 1% cutoff on the computed FDR
merged_PG = merged_PG[merged_PG['FDR'] <= 0.01]

# Write out the file to txt format for input into SAINT
merged_PG.to_csv('proteinGroups.txt',sep='\t')




####################################################

# Things to test: 
# does the simplify function work
# 	is the merging step going as expected
# 	is the score column correct for merged protein groups
# 	is the reverse column correct
# 	is the contaminant column correct
# 	is the sequence length correct
# 	is the id by site correct
# 	is the FDR being generated corectly
