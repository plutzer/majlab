import numpy as np
import pandas as pd
import networkx

def linked_FDRs(pathslist,func = np.max):
	'''
		Takes in a list of paths to seperately searched maxquant protein groups output files,
		and aggregates them into a single combined protein groups file with an extra column for
		the FDR.

		pathslist: list of strings - paths to individually searched protein groups files.
		func: function that can work on a list of doubles - will be the function used to combine
			protein group scores from sets of linked groups.
		returns: pandas dataframe of combined protein groups that can be writen out.
	'''

	# Make a combined protein groups file.
	combined_df = pd.concat([pd.read_csv(filename,delimiter = '\t') for filename in pathslist])


	# Generate a list of linked protein groups for this file using the protein IDs column.
	linked_groups = generate_linked_groups(combined_df['Protein IDs'])

	# Get a score for each linked group. This is a list with indexes corresponing to linked groups.
	linked_scores = [func(list(combined_df.iloc[indexes]['Score'])) for indexes in linked_groups]

	# Denote whether each linked group is a decoy. This is a boolean list, same indexes.
	linked_decoys = ['REV_' in ''.join(list(combined_df.iloc[indexes]['Protein IDs'])) for indexes in linked_groups]

	# Make a dataframe out of the linked columns
	linkages_info = pd.DataFrame(zip(linked_groups,linked_scores,linked_decoys),columns = ['groups','scores','decoys'])
	
	#Initialize a list for the FDRs
	FDRs = list(np.zeros([len(combined_df)]))
	

	for i in range(len(linked_groups)): #Iterating over the indexes of the linked groups list itself
	    # subset linkages_info to only elements that have a score less than than the current linked_group
	    subset = linkages_info[linkages_info['scores'] < linked_scores[i]]
	    FDR = np.sum((subset['decoys']))/(len(subset)-np.sum((subset['decoys'])))
	    for j in linked_groups[i]:
	    	FDRs[j] = FDR

	# Add the FDRs column to the dataframe.
	combined_df['FDR'] = FDRs

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
            also_in  = list(np.where([protein in group for group in groupedIDs])[0])
            
            # Iterate through also_in and update the adjacency matrix
            for index in also_in:
                if group_index != index:
                    adjMat[group_index,index] = 1
    
    # Generate a list of indexes for component subgraphs
    graph = networkx.from_numpy_matrix(adjMat)
    return [list(graph.subgraph(c).copy().nodes()) for c in networkx.connected_components(graph)]