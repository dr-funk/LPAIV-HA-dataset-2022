# a script that performs a number of analyses on sequence data and outputs a single file containing separate tabs presenting the graphs of the different analyses
# 05/2023 Mathis Funk
# this script is pretty modular, it sets up shared functions at the start, but runs the analyses one after the other completely individually
# analyses can be removed/added freely, in which case they will disappear/appear on the final file
#   in case you remove the first analysis, make sure to have another one include the plotlyjs library in the html file, or the graphs will not be displayed (see line 266)
#   in case you add an analysis, make sure to update the fig_nb dictionary to set the title to use in the tab of the final file

from Bio import SeqIO
import math
import re
import sys
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from jinja2 import Environment, FileSystemLoader

# pointing the script towards where the data is. if a path is provided on the command line, take that, else take the hard-coded one
if len(sys.argv)>1:
    data_path = sys.argv[1]
else:
    data_path='' # input path here to hardcode it

# set up the jinja environment
env = Environment(loader=FileSystemLoader(data_path))
template = env.get_template('jinja_template.html') # name of the jinja template file, needs to be in the data folder or a path specified in the line above instead of data_path
var1_names =[] # this list will get populated with all performed analyses
var2_names = ['all sequences','region','species','region-species'] # these are the different second variables to split the data by
data_files = ['america-poultry','eurasia-poultry','america-ans_cha','eurasia-ans_cha'] # these are the names of the subfolders in data/ and of the start of the sequence files they contain
graphs = {} # this is where the plotly graphs will be stored

# defining the colors of the different categories
colors ={"america-poultry": "#FDDBC7",
        "america-ans_cha": "#DC050C",
        "eurasia-poultry": "#74abe3",
        "eurasia-ans_cha": "#1965b0",}

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
var2_longform={'all sequences':'all poultry/wild bird sequences',
                'region':'split by region of isolation',
                'species':'split by species of isolation',
                'region-species':'split by region and species of isolation'}

# establishing the values taken by the second variable
var2_lists = {'region-species':data_files, 
              'region':['america','eurasia'],
              'species':['poultry','ans_cha'],
              'all sequences':['all']}

# attributing names and numbers to the figures. these will be used to annotate the tabs
fig_nb = {'number of adenines':'Figure S3a',
            'number of purines':'Figure S3c',
            'longest adenine stretch':'Figure S3b',
            'longest purine stretch':'Figure S3d',
            'substitutions to tribasic site':'Figure S3e',
            'substitutions to adenine stretch':'Figure S3f'}

# setting the mask threshold, any value above this will not get a bar label in the final figure
maskthreshold = 9

## setting up the subtype annotations on the figures
# coordinates have to be specified by hand, there are 4 rows and 4 lines
annotations_x =[0.19,0.46,0.74,1]
annotations_y =[1,0.73,0.45,0.16]
# making a list matching up the annotation text with the correct row and line
annotations_list = []
for t in range(16,0,-1):
    annotations_list.append(dict(text=f'H{t}', x=annotations_x[(t-1) % 4], y=annotations_y[math.floor((t-1)/4)], xref='paper', yref='paper', showarrow=False))

## this makes sure each value is present in the results dictionnary. this is important to have all analyses have the same x values, even if that particular value was not encountered in the analysis
def pad_results(results, maxresult, var1_name):
    for var2 in data_files:
        for t in range(1, 17):
            currentsubtype = f'H{t}'
            for i in range (2 if ' longest' in var1_name else 0, maxresult+1): # for the stretch analyses x values are handled differently since both 0 and 1 are grouped into <2
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
        if var2_name != 'all sequences':
            for var2 in var2_list:
                plotdf_norm.loc[(plotdf_norm['subtype'] == currentsubtype) & (plotdf_norm[var2_name]==abr_dict[var2]), 'percentage'] = resultdf.loc[(resultdf['subtype'] == currentsubtype) & (resultdf[var2_name]==abr_dict[var2]), 'percentage']/resultdf.loc[(resultdf['subtype'] == currentsubtype) & (resultdf[var2_name]==abr_dict[var2]), "percentage"].sum()*100
        else:
            plotdf_norm.loc[plotdf_norm['subtype'] == currentsubtype, 'percentage'] = resultdf.loc[resultdf['subtype'] == currentsubtype, 'percentage']/resultdf.loc[resultdf['subtype'] == currentsubtype, "percentage"].sum()*100    
    plotdf_norm.sort_values(by=[var1_name], inplace=True) # sorting is crucial for stretch analyses, they will be plotted as categorical data
    return plotdf_norm

## this sets up the offset, ie which bars will stack and which will not
def get_offset(var2_name):
    if var2_name != 'all sequences': # if the data is split by region and species there will be no offset (""), plotly will not stack anything as intended
        offset={'america-poultry':f'{"america" if var2_name != "species" else ""}{"poultry" if var2_name != "region" else ""}',
                'america-ans_cha':f'{"america" if var2_name != "species" else ""}{"ans_cha" if var2_name != "region" else ""}',
                'eurasia-poultry':f'{"eurasia" if var2_name != "species" else ""}{"poultry" if var2_name != "region" else ""}',
                'eurasia-ans_cha':f'{"eurasia" if var2_name != "species" else ""}{"ans_cha" if var2_name != "region" else ""}'}
    else: # in case the data is split, explicitly set the same arbitrary offset for all bars to ensure stacking
        offset={'america-poultry':'1',
                'america-ans_cha':'1',
                'eurasia-poultry':'1',
                'eurasia-ans_cha':'1'}
    return offset

## this finds the right values to shift the columns up when stacking so they don't overlap. top to bottom is america poultry, eurasia poultry, america ans_cha, eurasia ans_cha
def get_plotbase(var2_name, t, plotdf_norm):
    plotbase={'america-poultry':None, 
                'america-ans_cha':None if var2_name != 'region' else plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['america-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'],
                'eurasia-poultry':None if var2_name != 'species' else plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['america-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'],
                'eurasia-ans_cha':None}
    if var2_name == 'region':
        plotbase['eurasia-ans_cha'] = plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['eurasia-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage']
    if var2_name == 'species':
        plotbase['eurasia-ans_cha'] = plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['america-ans_cha']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage']
    if var2_name == 'all sequences':
        plotbase ={'america-poultry':None,
                    'eurasia-poultry':plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['america-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'],
                    'america-ans_cha':plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['america-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'].reset_index(drop=True).add(plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['eurasia-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'].reset_index(drop=True)),
                    'eurasia-ans_cha':plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['america-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'].reset_index(drop=True).add(plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['eurasia-poultry']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'].reset_index(drop=True)).add(plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict['america-ans_cha']) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'].reset_index(drop=True))}
    return plotbase

## get the bar labels from the dataframe and mask big values
def get_barlabels(plotdf_norm, maskthreshold, regionspecies, t):
    # pulling raw counts and percentages
    textdf =plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict[regionspecies]) & (plotdf_norm['subtype']==f'H{t}'), ['count', 'percentage']].copy()
    # masking numbers with low count and low percentage
    textdf.mask((textdf['count'] > maskthreshold) | (textdf['percentage']<50), inplace=True)
    # masking zeros
    textdf.mask(textdf == 0, inplace=True)
    text=textdf['count']
    return text

## actually makes the graph using plotly graph objects
def make_graph(plotdf_norm, var1_name, var2_name, abr_dict):
    # setting up a grid of 4x4 subplots
    fig =  make_subplots(rows=4, cols=4, x_title=var1_name, y_title='percent of sequences')

    # set up x-offsets
    offset = get_offset(var2_name)
    # iterating through subtypes
    for t in range(1, 17):
        # getting the correct y-offsets
        plotbase = get_plotbase(var2_name, t, plotdf_norm)    

        # drawing graphs. we always iterate through region-species since we want to draw a bar for all combinations, regardless of what they are split by
        for regionspecies in data_files: 

            # pull out the bar label data and mask anything that isn't below 10   
            textdf = get_barlabels(plotdf_norm, maskthreshold, regionspecies, t)

            # time to actually draw the graph now that we have all the pieces. the order of customdata fields is the same as the column order of plotdf
            fig.add_trace(go.Bar(name=abr_dict[regionspecies], 
                                x=[i for i in range(maxresult+1)] if ' longest' not in var1_name else [i if i != 1 else 0 for i in range(1, maxresult+1)], 
                                y=plotdf_norm.loc[(plotdf_norm['region-species']==abr_dict[regionspecies]) & (plotdf_norm['subtype']==f'H{t}'), 'percentage'],
                                customdata=plotdf_norm[(plotdf_norm['region-species']==abr_dict[regionspecies]) & (plotdf_norm['subtype']==f'H{t}')], 
                                hovertemplate='<b>subtype: %{customdata[0]}</b><br>region: %{customdata[2]} <br>species: %{customdata[3]} <br>count: %{customdata[5]}<br>percentage: %{customdata[6]:.2f}%<extra></extra>',
                                offsetgroup=offset[regionspecies],
                                base= plotbase[regionspecies],
                                marker_color=[colors[regionspecies],]*16,
                                text=textdf,
                                textposition='auto',
                                showlegend = True if t == 1 else False,
                                legendgroup=regionspecies
                                ), 
                        row=math.floor((t-1)/4)+1,
                        col=(t-1) % 4+1
                        )

    # some formatting and adding a title
    fig.update_layout(template='simple_white',
                    hoverlabel={'bgcolor':'white',
                                'font_size':16},
                                autosize=True,
                    title={'text': f"{var1_name} in terrestrial poultry and wild aquatic bird LPAIV sequences {f'{var2_longform[var2_name]}' if var2_name !='all sequences' else ''}" ,
                            'y':0.99,
                            'x':0.5,
                            'xanchor': 'center',
                            'yanchor': 'top'},
                    margin={'b':60,'l':60,'r':25,'t':40})

    # adding the subtype annotations
    for annotation in annotations_list:
        fig.add_annotation(annotation)

    # formatting axes. stretch analysis again need to be handled separately because of the <2 category spanning both 0 and 1 values
    if ' longest' not in var1_name:
        fig.update_xaxes(showticklabels=True,
                        tickmode = 'linear',
                        tick0 = 0,
                        dtick= 1)
    else:
        fig.update_xaxes(showticklabels=True,
                    type="category",
                    tickmode = 'array',
                    tickvals = [i if i != 1 else 0 for i in range(1, maxresult+1)],
                    ticktext = [str(i) if i != 1 else '<2' for i in range(1, maxresult+1)])
    fig.update_yaxes(range=[0, 100])
    if var1_name == 'substitutions to adenine stretch':
        fig.update_xaxes(range=[-0.5,7.5])
    
    return fig

## A count
var1_name = 'number of adenines'
var1_names.append(var1_name)
graphs[var1_name]=[]
for var2_name in var2_names:
    print(f'running {var1_name} on {var2_name}')
    ## intializing variables
    results = {}
    maxresult = 0 # keeping track of highest encountered value
    var2_list = var2_lists[var2_name]

    ## getting A count
    for var2 in data_files:
        for t in range(1, 17):
            currentsubtype = f"H{t}"
            if currentsubtype not in results.keys():
                results[currentsubtype]={}
            results[currentsubtype][var2]={0:0}
            infile = f"{data_path}/by_region-species/{var2}/{var2}_h{t}_unique_CSregion.fasta" # this is the path of the file containing all sequences for this region-species combo of this subtype
            for record in SeqIO.parse(infile, "fasta"): # iterating through each sequence
                seqslice = str(record.seq)[31:43].upper() # slicing up the sequence to get only P4-P1, the region of interest
                anumber = seqslice.upper().count("A") # simply counting how many As we have
                if anumber > maxresult:
                    maxresult=anumber
                if anumber in results[currentsubtype][var2].keys():
                    results[currentsubtype][var2][anumber]+=1 
                else:
                    results[currentsubtype][var2][anumber]=1

    ## making sure each value is present in each dataset
    results = pad_results(results, maxresult, var1_name)

    ## reshaping the dictionary to a list and building a dataframe
    resultdf = dict_to_df(results, abr_dict)

    ## normalizing the data                        
    plotdf_norm = normalize_df(resultdf, var2_name)

    ## plotting the figure
    fig =  make_graph(plotdf_norm, var1_name, var2_name, abr_dict)

    ## adding the html string to the graphs. The plotly library scripts are huge and only need to be present once. we only include them with the first graph to avoid massive file sizes
    graphs[var1_name].append([var2_name,fig.to_html(full_html=False, include_plotlyjs=True if var2_name == 'all sequences' and var1_name == 'number of adenines' else False).strip('<div>').rstrip('</div>')])

## adenine stretch
var1_name = 'longest adenine stretch'
var1_names.append(var1_name)
graphs[var1_name]=[]
for var2_name in var2_names:
    print(f'running {var1_name} on {var2_name}')
    ## intializing variables
    results = {}
    maxresult = 0 # keeping track of highest encountered value
    var2_list = var2_lists[var2_name]

    ## getting longest A-stretches
    for var2 in data_files:
        for t in range(1, 17):
            currentsubtype = f"H{t}"
            if currentsubtype not in results.keys():
                results[currentsubtype]={}
            results[currentsubtype][var2]={0:0}
            infile = f"{data_path}/by_region-species/{var2}/{var2}_h{t}_unique_CSregion.fasta"
            for record in SeqIO.parse(infile, "fasta"):
                seqslice = str(record.seq)[31:43].upper()
                astretch = len(max(re.compile("(A+A)*").findall(seqslice.upper()))) # using regex to find all stretches of 2 or more As, then taking the longest one via max(), then getting it's length via len()
                if astretch in results[currentsubtype][var2].keys():
                    results[currentsubtype][var2][astretch]+=1
                else:
                    results[currentsubtype][var2][astretch]=1
                if astretch > maxresult:
                    maxresult = astretch
  
    ## making sure each value is present in each dataset
    results = pad_results(results, maxresult, var1_name)

    ## reshaping the dictionary to a list and building a dataframe
    resultdf = dict_to_df(results, abr_dict)

    ## normalizing the data                        
    plotdf_norm = normalize_df(resultdf, var2_name)

    ## setting up some stuff for plotting
    offset = get_offset(var2_name)

    ## plotting the figure
    fig =  make_graph(plotdf_norm, var1_name, var2_name, abr_dict)

    # fig.show()
    graphs[var1_name].append([var2_name,fig.to_html(full_html=False, include_plotlyjs=False).strip('<div>').rstrip('</div>')])

## purine count
var1_name = 'number of purines'
var1_names.append(var1_name)
graphs[var1_name]=[]
for var2_name in var2_names:
    print(f'running {var1_name} on {var2_name}')
    ## intializing variables
    results = {}
    maxresult = 0 # keeping track of biggest variable
    var2_list = var2_lists[var2_name]

    ## getting purine count
    for var2 in data_files:
        for t in range(1, 17):
            currentsubtype = f"H{t}"
            if currentsubtype not in results.keys():
                results[currentsubtype]={}
            results[currentsubtype][var2]={0:0}
            infile = f"{data_path}/by_region-species/{var2}/{var2}_h{t}_unique_CSregion.fasta"
            for record in SeqIO.parse(infile, "fasta"):
                seqslice = str(record.seq)[31:43].upper() 
                anumber = seqslice.upper().count("A")+seqslice.upper().count('G') # simply take the number of As and the number of Gs this time
                if anumber > maxresult:
                    maxresult=anumber
                if anumber in results[currentsubtype][var2].keys():
                    results[currentsubtype][var2][anumber]+=1 
                else:
                    results[currentsubtype][var2][anumber]=1

    ## making sure each value is present in each dataset
    results = pad_results(results, maxresult, var1_name)

    ## reshaping the dictionary to a list and building a dataframe
    resultdf = dict_to_df(results, abr_dict)

    ## normalizing the data                        
    plotdf_norm = normalize_df(resultdf, var2_name)

    ## plotting the figure
    fig =  make_graph(plotdf_norm, var1_name, var2_name, abr_dict)

    graphs[var1_name].append([var2_name,fig.to_html(full_html=False, include_plotlyjs=False).strip('<div>').rstrip('</div>')])


## Purine stretch
var1_name = 'longest purine stretch'
var1_names.append(var1_name)
graphs[var1_name]=[]
for var2_name in var2_names:
    print(f'running {var1_name} on {var2_name}')

    ## intializing variables
    results = {}
    maxresult = 0
    var2_list = var2_lists[var2_name]

    ## getting longest purine stretches
    for var2 in data_files:
        for t in range(1, 17):
            currentsubtype = f"H{t}"
            if currentsubtype not in results.keys():
                results[currentsubtype]={}
            results[currentsubtype][var2]={0:0}
            infile = f"{data_path}/by_region-species/{var2}/{var2}_h{t}_unique_CSregion.fasta"
            for record in SeqIO.parse(infile, "fasta"):
                seqslice = str(record.seq)[31:43].upper() 
                astretch = len(max(re.compile("([AG]+[AG])*").findall(seqslice.upper()))) # same approach as above, but both A and G are valid characters to start/continue the stretch
                if astretch in results[currentsubtype][var2].keys():
                    results[currentsubtype][var2][astretch]+=1
                else:
                    results[currentsubtype][var2][astretch]=1
                if astretch > maxresult:
                    maxresult = astretch

    ## making sure each value is present in each dataset
    results = pad_results(results, maxresult, var1_name)

    ## reshaping the dictionary to a list and building a dataframe
    resultdf = dict_to_df(results, abr_dict)

    ## normalizing the data                        
    plotdf_norm = normalize_df(resultdf, var2_name)

    ## setting up some stuff for plotting
    offset = get_offset(var2_name)

    ## plotting the figure
    fig =  make_graph(plotdf_norm, var1_name, var2_name, abr_dict)

    # fig.show()
    graphs[var1_name].append([var2_name,fig.to_html(full_html=False, include_plotlyjs=False).strip('<div>').rstrip('</div>')])

## muts to basic
var1_name = 'substitutions to tribasic site'
var1_names.append(var1_name)
graphs[var1_name]=[]
for var2_name in var2_names:

    ## intializing variables
    results = {}
    maxresult = 0
    var2_list = var2_lists[var2_name]
    target = [['AAA','AAG','AGA','AGG'],['AAA','AAG','AGA','AGG'],['AAA','AAG','AGA','AGG'],['AGA','AGG']] # list of lists of codons acceptable at each position
    print(f'running {var1_name} on {var2_name}')

    ## getting A count
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
    results = pad_results(results, maxresult, var1_name)

    ## reshaping the dictionary to a list and building a dataframe
    resultdf = dict_to_df(results, abr_dict)

    ## normalizing the data                        
    plotdf_norm = normalize_df(resultdf, var2_name)

    ## plotting the figure
    fig =  make_graph(plotdf_norm, var1_name, var2_name, abr_dict)

    graphs[var1_name].append([var2_name,fig.to_html(full_html=False, include_plotlyjs=False).strip('<div>').rstrip('</div>')])

## muts to adenine stretch
var1_name = 'substitutions to adenine stretch'
var1_names.append(var1_name)
graphs[var1_name]={}
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
for var2_name in var2_names:
    print(f'running {var1_name} on {var2_name}')
    ## intializing variables
    results = {}
    maxresult = 0
    var2_list = var2_lists[var2_name]
    graphs[var1_name][var2_name]=[]

    ## getting number of subs
    # here the sequence of allowed codons is crucial as it determines the length of the A stretch. therefore we compare all 12 nts of our sequence to all 12 codons of the target.
    for threshold in range(3, 11): # varying the threshold from 3 to 10
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
                    if var2_name == "all sequences" or var2 in file: # the second part makes sure var2 "eurasia" will look at files "eurasia-poultry" and "eurasia-ans_cha" but not the american files. first part is an escape in case we want to look at all sequences
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
        results = pad_results(results, maxresult, var1_name)

        ## reshaping the dictionary to a list and building a dataframe
        resultdf = dict_to_df(results, abr_dict)

        ## normalizing the data                        
        plotdf_norm = normalize_df(resultdf, var2_name)

        ## plotting the figure
        fig =  make_graph(plotdf_norm, var1_name, var2_name, abr_dict)

        graphs[var1_name][var2_name].append([threshold,fig.to_html(full_html=False, include_plotlyjs=False).strip('<div>').rstrip('</div>')])

print('making figure')
with open(f'{data_path}/multifig_jinja.html', 'w') as fj:
    #calling jinja and passing along all our variables
    fj.write(template.render(var2_names=var2_names, graphs=graphs, var1_names =var1_names, fig_nb=fig_nb, var2_longform=var2_longform))
