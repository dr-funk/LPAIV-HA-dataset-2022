# a script that checks how many substitutions are necessary to reach an A stretch over a given threshold for all sequences of a dataset
# 05/2023 Mathis Funk

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import sys

# pointing the script towards where the data is. if a path is provided on the command line, take that, else take the hard-coded one
if len(sys.argv)>1:
    data_path = sys.argv[1]
else:
    data_path='' # input path here to hardcode it

# specifying which variable to split by
if len(sys.argv)>2:
    var2_name = sys.argv[2]
else:
    var2_name = 'none' # hardcode which variable to split the data by. set to "none" by default so that if no variable is provided, no splitting is done

if var2_name not in ('region', 'species', 'region-species','none'):
    sys.exit('second variable has to be one of: region, species, region-species, none')

# setting up parameters
results = {}
maxresult = 0
# setting up the title of the figure

data_files = ['america-poultry','eurasia-poultry','america-ans_cha','eurasia-ans_cha'] # these are the names of the subfolders in data/ and of the start of the sequence files they contain

## matching up long form names with the shorthand used in the script. the longer version is used to annotate axes and legends in the graphs
# first the values of the second variable
abr_dict={'ans_cha':'aquatic wild birds', 
        'poultry':'terrestrial poultry', 
        'eurasia':'AEO',    
        'america':'Americas',
        'america-poultry':'terrestrial poultry from the Americas',
        'america-ans_cha':'aquatic wild birds from the Americas',
        'eurasia-poultry':'terrestrial poultry from the AEO',
        'eurasia-ans_cha':'aquatic wild birds from the AEO'}
# second the names of the second variable
var2_longform={'none':'all sequences',
                'region':'split by region of isolation',
                'species':'split by species of isolation',
                'region-species':'split by region and species of isolation'}

# establishing the values taken by the second variable
var2_lists = {'region-species':data_files, 
              'region':['america','eurasia'],
              'species':['poultry','ans_cha'],
              'none':['all']}

palette_dict = {'region-species':['#f06e73', '#74abe3', '#DC050C', '#1965b0'],
                'region':['#DC050C','#1965b0'],
                'species':['#74abe3','#1965b0'],
                'none':'Blues_r'}


## this makes sure each value is present in the results dictionnary. this is important to have all analyses have the same x values, even if that particular value was not encountered in the analysis
def pad_results(results, maxresult):
    for var2 in data_files:
        for t in range(1, 17):
            currentsubtype = f'H{t}'
            for i in range (0, maxresult+1): # for the stretch analyses x values are handled differently since both 0 and 1 are grouped into <2
                if i not in results[currentsubtype][var2].keys():
                    results[currentsubtype][var2][i]=0
    return results

## this turns the dictionary into a dataframe for plotting
# here we iterate through the different levels of the result dictionary and build a list of lists containing the values of each level. this can then easily be turned into a dataframe
def dict_to_df(results, abr_dict):
    resultdf = pd.DataFrame([[k1, abr_dict[k2], abr_dict[k2.split('-')[0]], abr_dict[k2.split('-')[1]], k3, results[k1][k2][k3], results[k1][k2][k3]]
                        for k1 in results.keys()
                            for k2 in results[k1].keys()
                                for k3 in results[k1][k2].keys()], 
                        columns = ['subtype','region-species','region','species', var1_name,'count','percentage'])
    return resultdf

## this normalizes the dataframe based on the 2nd variable
# this automatically normalizes by region, species, region-species or all, as apropriate
def normalize_df(resultdf, var2_name):
    plotdf_norm = resultdf.copy() # important, deep copy the original dataframe
    for t in range(1, 17):
        currentsubtype = f'H{t}'
        if var2_name != 'none':
            for var2 in var2_list:
                plotdf_norm.loc[(plotdf_norm['subtype'] == currentsubtype) & (plotdf_norm[var2_name]==abr_dict[var2]), 'percentage'] = resultdf.loc[(resultdf['subtype'] == currentsubtype) & (resultdf[var2_name]==abr_dict[var2]), 'percentage']/resultdf.loc[(resultdf['subtype'] == currentsubtype) & (resultdf[var2_name]==abr_dict[var2]), "percentage"].sum()*100
        else:
            plotdf_norm.loc[plotdf_norm['subtype'] == currentsubtype, 'percentage'] = resultdf.loc[resultdf['subtype'] == currentsubtype, 'percentage']/resultdf.loc[resultdf['subtype'] == currentsubtype, "percentage"].sum()*100    
    plotdf_norm.sort_values(by=[var1_name], inplace=True) # sorting is crucial for stretch analyses, they will be plotted as categorical data
    return plotdf_norm

## get the bar labels from the dataframe and mask big values
def get_barlabels(plotdf_norm, t):
    labels = []
    if var2_name != 'none':
        # collapsing the counts by var2 to sum up both species in region mode for example
        collapsed_df = plotdf_norm[plotdf_norm['subtype']==f'H{t}'].groupby([var1_name, var2_name], as_index=False)['count'].sum()
        # pulling raw counts from the collapsed dataframe
        for var2_value in var2_list:
            counts = collapsed_df.loc[(collapsed_df[var2_name]==abr_dict[var2_value]), ['count']].values.tolist()
            # cleaning up the labels by ignoring zeros and pulling the values out of the individual lists
            label_text = [count[0] if count != [0] else '' for count in counts]
            labels.append(label_text)
    else:
        # collapsing the counts 
        collapsed_df = plotdf_norm[plotdf_norm['subtype']==f'H{t}'].groupby([var1_name], as_index=False)['count'].sum()
        # pulling raw counts from the collapsed dataframe
        counts = collapsed_df['count'].values.tolist()
        # cleaning up the labels by ignoring zeros and pulling the values out of the individual lists
        label_text = [count if count != 0 else '' for count in counts]
        labels = label_text
    return labels

## running the analysis
# make dictionnaries of allowed codons and count sequences
# this makes a lot of dictionaries, one for each way of splitting the data. this is important to restrict codon usage by species/region/both later on
allowed_codons ={}
for t in range (1, 17): # iterating through subtypes
    subtype = f'H{t}'
    if subtype not in allowed_codons.keys():
        allowed_codons[subtype]={}
        for i in ['all','america','eurasia','poultry','ans_cha']: 
            allowed_codons[subtype][i] = [[],[],[],[]] # setting up list of lists with one empty list for each position, these will contained all the allowed codons
    for var2 in data_files:
        allowed_codons[subtype][var2] = [[],[],[],[]]
        for record in SeqIO.parse(f'{data_path}/by_region-species/{var2}/{var2}_h{t}_unique_CSregion.fasta', 'fasta'):
            for c in range(0,4): # iterate through each codon and append it to every dictionary it needs to be in
                # add to region-species dict, no splicing of var2 needed, dictionary key is whole var2
                if str(record.seq).upper()[31+3*c:31+3*c+3] not in allowed_codons[subtype][var2][c]:
                    allowed_codons[subtype][var2][c].append(str(record.seq).upper()[31+3*c:31+3*c+3])
                # add to all dict, dictionary key is the same no matter the origin of the sequence
                if str(record.seq).upper()[31+3*c:31+3*c+3] not in allowed_codons[subtype]['all'][c]:
                    allowed_codons[subtype]['all'][c].append(str(record.seq).upper()[31+3*c:31+3*c+3])
                # add to region dict, dictionary key is the first part of the the region_species variable var2
                if str(record.seq).upper()[31+3*c:31+3*c+3] not in allowed_codons[subtype][var2.split('-')[0]][c]:
                    allowed_codons[subtype][var2.split('-')[0]][c].append(str(record.seq).upper()[31+3*c:31+3*c+3])
                # add to species dict, dictionary key is the second part of the the region_species variable var2
                if str(record.seq).upper()[31+3*c:31+3*c+3] not in allowed_codons[subtype][var2.split('-')[1]][c]:
                    allowed_codons[subtype][var2.split('-')[1]][c].append(str(record.seq).upper()[31+3*c:31+3*c+3])

# now we run the actual analysis
## intializing variables
var2_list = var2_lists[var2_name]

## getting number of subs
# here the sequence of allowed codons is crucial as it determines the length of the A stretch. therefore we compare all 12 nts of our sequence to all 12 codons of the target.
for threshold in range(3, 11): # varying the threshold from 3 to 10
    var1_name = f'substitutions to have at least {threshold} consecutive As'
    print(f'running {var1_name} on {var2_name if var2_name != "none" else "all sequences"}')
    for var2 in var2_list:
        for t in range(1,17):
            currentsubtype = f'H{t}'
            # build all valid combinations of codons for each subtype
            possible_loops = []
            for c1 in allowed_codons[currentsubtype][var2][0]:
                for c2 in allowed_codons[currentsubtype][var2][1]:
                    for c3 in allowed_codons[currentsubtype][var2][2]:
                        for c4 in allowed_codons[currentsubtype][var2][3]:
                            if c4 in ['AGA','AGG','CGT','CGA','CGC','CGG']: # postion P1 is fixed to R, don't even bother with any other codons
                                if f'{c1}{c2}{c3}{c4}'.find('A'*threshold) != -1: # checking that there's enough As in the loop
                                    possible_loops.append(f'{c1}{c2}{c3}{c4}') # keeping the valid combination
            if currentsubtype not in results.keys():
                results[currentsubtype]={}
            for file in ('america-poultry','america-ans_cha','eurasia-poultry','eurasia-ans_cha'): # we iterate through all 4 possible files of sequences and check if we need to look at the sequences it contains
                if var2_name == "none" or var2 in file: # the second part makes sure var2 "eurasia" will look at files "eurasia-poultry" and "eurasia-ans_cha" but not the american files. first part is an escape in case we want to look at all sequences
                    results[currentsubtype][file]={0:0}
                    if possible_loops == []: # no sequences fulfil the requirement, skip counting mutations
                        results[currentsubtype][file]={}
                    else:
                        for record in SeqIO.parse(f'{data_path}/by_region-species/{file}/{file}_h{t}_unique_CSregion.fasta', 'fasta'): # iterate through sequences of currently considered file
                            seq = str(record.seq).upper()[31:43]
                            muts = [] # list which will keep track of the number of differences between sequence and each valid loop
                            if seq in possible_loops: # no mutations needed, the sequence is already a valid loop
                                muts = [0]
                            else:
                                for loop in possible_loops: # iterate through all valid loops
                                    to_target = 0 # initiate counter
                                    for c in range(0, len(seq)): # go through each nucleotide of the sequence
                                        if seq[c] != loop[c]: # and compare it to the nucleotide of the current valid loop
                                            to_target += 1
                                    muts.append(to_target) # append the number of differences between sequence loop and valid loop
                                    if to_target == 1: # can't need less than 1 mutation when in this loop, no need to scan the rest
                                        break
                            if min(muts) in results[currentsubtype][file].keys():
                                results[currentsubtype][file][min(muts)]+=1
                            else:
                                results[currentsubtype][file][min(muts)]=1
                            if min(muts) > maxresult:
                                maxresult = min(muts)

    ## making sure each value is present in each dataset
    results = pad_results(results, maxresult)

    ## reshaping the dictionary to a list and building a dataframe
    resultdf = dict_to_df(results, abr_dict)

    ## normalizing the data                        
    plotdf_norm = normalize_df(resultdf, var2_name)

    ## plotting
    ## initializing the figure
    plt.rcParams.update({'font.size': 16})
    sns.set_style("whitegrid")
    fig = plt.figure(figsize =(15,15))
    for t in range(1, 17): # iterating through subtypes
        label_text = get_barlabels(plotdf_norm, t) # getting bar labels
        ax = plt.subplot2grid((4,4), (math.floor((t-1)/4),(t-1) % 4)) # this positions the figure correctly in the 4x4 grid
        # next line squashes the data together, e.g., combine america-poultry and america-ans_cha if plotting by region or combining all subsets if plotting without splitting
        sns.barplot(data=plotdf_norm[plotdf_norm['subtype']==f'H{t}'].groupby([var1_name, var2_name], as_index=False)['percentage'].sum() if var2_name != 'none' else plotdf_norm[plotdf_norm['subtype']==f'H{t}'].groupby([var1_name], as_index=False)['percentage'].sum(), 
                    ax=ax, 
                    x=var1_name, 
                    y='percentage',
                    hue=var2_name if var2_name != 'none' else None, 
                    hue_order=[abr_dict[var2] for var2 in var2_list] if var2_name != 'none' else None, # making sure the bars are plotted in the right order to make sure they get the right color
                    palette=palette_dict[var2_name]) # grabbing the right color palette
        # bar labels
        bar = 0
        for c in ax.containers:
            ax.bar_label(c, labels=label_text[bar] if var2_name != 'none' else label_text, fontsize=12) # label_text is a list of lists if the data is split, but only a normal list if data isn't split
            bar += 1  
        ax.set(ylim = (0,100))
        ax.grid(False)
        ax.set(xlabel=None)
        ax.set(ylabel=None)
        if var2_name != 'none': # no legend gets plotted if there is a single dataseries, avoiding an error by not trying to remove a legend that doesn't exist
            ax.get_legend().remove()
        sns.despine()
    fig.suptitle(f'{var1_name}, {var2_longform[var2_name]}') # setting a big ol' title
    fig.tight_layout()
    fig.savefig(f"{data_path}/muts_to_{threshold}_A_stretch_by{var2_name}.pdf")

