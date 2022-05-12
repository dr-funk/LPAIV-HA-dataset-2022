# a script that checks how many substitutions are necessary to reach an MBCS for all sequences of a dataset
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

results = {}
maxresult = 0
tribasic = True
consecutive = False
tetrabasic = False
var1_name = f'substitutions to {"tri" if tribasic == True else "tetra"}basic {"non-" if consecutive == False else ""}consecutive site'

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
        # collapsing the counts by var2 to sum up both species in region mode for example
        collapsed_df = plotdf_norm[plotdf_norm['subtype']==f'H{t}'].groupby([var1_name], as_index=False)['count'].sum()
        # pulling raw counts from the collapsed dataframe
        counts = collapsed_df['count'].values.tolist()
        print(counts)
        # cleaning up the labels by ignoring zeros and pulling the values out of the individual lists
        label_text = [count if count != 0 else '' for count in counts]
        labels = label_text
    return labels

## intializing variables
results = {}
maxresult = 0
var2_list = var2_lists[var2_name]
target = [['AAA','AAG','AGA','AGG'],['AAA','AAG','AGA','AGG'],['AAA','AAG','AGA','AGG'],['AGA','AGG']] # list of lists of codons acceptable at each position
print(f'running {var1_name} on {var2_name}')

for var2 in data_files:
    for t in range(1, 17):
        currentsubtype = f'H{t}'
        if currentsubtype not in results.keys():
            results[currentsubtype]={}
        results[currentsubtype][var2]={0:0}
        infile = f"{data_path}/by_region-species/{var2}/{var2}_h{t}_unique_CSregion.fasta"
        for record in SeqIO.parse(infile, 'fasta'):
            sequence = str(record.seq).upper()
            mutsminpos = [0,0,0,0] # one counter per position, each starting on 0
            sequence = sequence[31:43] # slicing out sequence of interest
            for n in range (0,4): # iterating through positions
                codon = sequence[3*n:3*n+3] # getting the codon for the position
                mutstotarget = [0 for item in target[n]] # one counter per valid codon at this position, all starting at 0
                if codon not in target[n]: # if the codon is in the list of valid codons, counter will remain at 0, else:
                    for c in range (0, len(target[n])): # iterate through each valid codon
                        for i in range(0,3): # iterate through each position of the valid codon
                            if codon[i] != target[n][c][i] and target[n][c][i] != 'X': # comparing sequence codon position with valid codon position. second part allows wildcard positions IN THE TARGET CODON by including an X instead of one of the nucleotide, everything will match this. 
                                mutstotarget[c] += 1 # if different, 1 substitution is required
                mutsminpos[n]=min(mutstotarget) # the minimum number of substitutions at  position n is the minimum number of mutations required to reach any target codon valid at position n
                # minimum number of substitutions required for this sequence is the sum of the minimum number required at each position of the target motif
                # we force an R in P1 which means this position (index 3 in the list) has to ALLWAYS be checked. one other position is allowed to stay non-substituted since we are looking for a tribasic site. this will be the position requiring the highest amount of substitutions
                mutsmin = min(
                    sum(mutsminpos), # odds of this being the single lowest value are null to be honest, since the 3 subsequent terms will alway be lower or equal. still, keeping this and removing the 3 subsequent ones allows to enforce tetrabasic sites
                    sum(mutsminpos) - mutsminpos[0],
                    sum(mutsminpos) - mutsminpos[1],
                    sum(mutsminpos) - mutsminpos[2])
            if mutsmin in results[currentsubtype][var2].keys():
                results[currentsubtype][var2][mutsmin]+=1
            else:
                results[currentsubtype][var2][mutsmin]=1
            if mutsmin > maxresult:
                maxresult = mutsmin

## making sure each value is present in each dataset
results = pad_results(results, maxresult)

## reshaping the dictionary to a list and building a dataframe
resultdf = dict_to_df(results, abr_dict)

## normalizing the data                        
plotdf_norm = normalize_df(resultdf, var2_name)

## plotting
## initializing the figure
sns.color_palette("Set2",len(var2_list))
plt.rcParams.update({'font.size': 16})
sns.set_style("whitegrid")
fig = plt.figure(figsize =(15,15))
for t in range(1, 17):
    ## get bar labels
    label_text = get_barlabels(plotdf_norm, t)
    ax = plt.subplot2grid((4,4), (math.floor((t-1)/4),(t-1) % 4))
    sns.barplot(data=plotdf_norm[plotdf_norm['subtype']==f'H{t}'].groupby([var1_name, var2_name], as_index=False)['percentage'].sum() if var2_name != 'none' else plotdf_norm[plotdf_norm['subtype']==f'H{t}'].groupby([var1_name], as_index=False)['percentage'].sum(), 
                ax=ax, 
                x=var1_name, 
                y='percentage',
                hue=var2_name if var2_name != 'none' else None, 
                hue_order=[abr_dict[var2] for var2 in var2_list] if var2_name != 'none' else None, # making sure the bars are plotted in the right order to make sure they get the right color
                palette=palette_dict[var2_name])
    bar = 0
    for c in ax.containers:
        ax.bar_label(c, labels=label_text[bar] if var2_name != 'none' else label_text, fontsize=12)
        bar += 1  
    ax.set(ylim = (0,100))
    ax.grid(False)
    ax.set(xlabel=None)
    ax.set(ylabel=None)
    if var2_name != 'none':
        ax.get_legend().remove()
    sns.despine()
fig.suptitle(f'{var1_name}, {var2_longform[var2_name]}')
fig.tight_layout()
fig.savefig(f"{data_path}/muts_to_{'tri' if tribasic == True else 'tetra'}basic_by{var2_name}_{'non-cons' if consecutive == False else 'cons'}.pdf")

