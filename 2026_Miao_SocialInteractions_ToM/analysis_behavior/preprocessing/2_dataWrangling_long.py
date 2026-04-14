# Behavior data wrangling
# Author: Zizhuang Miao

# This script extracts behavioral rating data, both continuous and intermittent, obtained from the online study, and put them into single long-format data sheets. 
# One data file will be generated for one type of questions collected (e.g., one for social interactions, another for theory of mind).
# Meanwhile, because some behavioral metrics are related to quality control (proportions of valid data, for example),
# this script also generates a new quality control file at the last code chunk.

# %%
# get and define metadata
import pandas as pd
import os
import numpy as np
import glob

rootDir = 'C:\\'
dataDir = os.path.join(rootDir, 'data')
outputDir = os.path.join(rootDir, 'derivedData')

# %%
# ______________________________________________
# Social interactions and ToM
# ______________________________________________

# %% 
# generate the backbone dataframe
parInfo = pd.read_csv(os.path.join(rootDir, 'participants_cleaned_int.csv'))
parList = parInfo['parID']
numericID = parInfo['numericID']

outputDf_social = pd.DataFrame()
outputDf_tom = pd.DataFrame()

frameRate = []
sampling_interval = []

for j, par in enumerate(parList):
    dataFile = os.path.join(dataDir, par+'_group3.1-3.2-social-tom.csv')
    raw = pd.read_csv(dataFile)
    frameRate.append(raw.loc[0, 'frameRate'])
    rawData = raw.loc[~np.isnan(raw['fixation.started'])]    # discard other rows
    rawData = rawData.reset_index()
    rawData = rawData.drop(index=[0, 19, 38, 57])    # discard practice trials
    rawData = rawData.reset_index()
    parDf = pd.DataFrame()
    for i in range(len(par)-1, -1, -1):
        if par[i].isnumeric():
            parID = int(par[i])
            break
    if parID % 2 == 0:
        narrativeList = [7]*9+[8]*9+[5]*9+[6]*9+[3]*9+[4]*9+[1]*9+[2]*9
    else:
        narrativeList = [8]*9+[7]*9+[6]*9+[5]*9+[4]*9+[3]*9+[2]*9+[1]*9
    situationList = list(range(1,10))*8
    runList = [1]*18+[2]*18+[3]*18+[4]*18

    # get sampling interval
    # because the first few sampling intervals are not stable, I will use the later ones
    sampling_interval_obtained = False
    timeDiff = 83    # default value in case it cannot be found in the continuous rating data
    for i in range(len(rawData)):
        if not np.isnan(rawData.loc[i, 'fixation.started']) and not sampling_interval_obtained:
            con = 0
            exec('con = ' + rawData.loc[i, 'continuousResponses_audio'])
            try:
                timeDiff = con[5][1] - con[4][1]
                sampling_interval.append(timeDiff)
                sampling_interval_obtained = True
            except:
                sampling_interval_obtained = False

    for i in range(len(rawData)):
        #parDf.loc[i, 'narrative'] = ''
        #parDf.loc[i, 'situation'] = ''    
        #parDf.loc[i, 'run'] = 
        # I need to modify the experiment file to save those information
        parDf.loc[i, 'narrative'] = narrativeList[i]
        parDf.loc[i, 'situation'] = situationList[i]
        parDf.loc[i, 'run'] = runList[i]
        parDf.loc[i, 'trial'] = i%18+1
        parDf.loc[i, 'con_rating'] = rawData.loc[i, 'continuousResponses_audio'] if i <= 35 else rawData.loc[i, 'continuousResponses_text']
        parDf.loc[i, 'int_rating'] = rawData.loc[i, 'intQuestionSlider.response']

        # calculate the proportions of valid responses in continuous rating
        # note that the data was stored as strings, although they look like a python list
        numberCon = (rawData.loc[i, 'audioPre.stopped'] - rawData.loc[i, 'audioPre.started'])*1000 // timeDiff if i <= 35 else (rawData.loc[i, 'textPre.stopped'] - rawData.loc[i, 'textPre.started'])*1000 // timeDiff
        numberValid = 0
        con = 0    # avoid error warnings
        exec('con=' + parDf.loc[i, 'con_rating'])
        for k in range(len(con)):
            if con[k][0] != 999:
                numberValid += 1
        parDf.loc[i, 'con_rating_valid'] = numberValid/numberCon


    parDf['ID'] = numericID[j]
    parDf = parDf[['ID', 'run', 'trial', 'narrative', 'situation', 'con_rating', 'int_rating', 'con_rating_valid']]
    if raw.loc[0, 'exp_question'] == 'social':
        outputDf_social = pd.concat([outputDf_social, parDf])
    elif raw.loc[0, 'exp_question'] == 'tom':
        outputDf_tom = pd.concat([outputDf_tom, parDf])

# %%
# save long-formatted data
outputDf_social.to_csv(os.path.join(outputDir, 'social_long.csv'), index=False)
outputDf_tom.to_csv(os.path.join(outputDir, 'tom_long.csv'), index=False)

# %%
# calculate some summary statistics for quality control
quality_catch = pd.read_csv(os.path.join(rootDir, 'preprocess', 'outputs', 'qualityControl_catchQuestion_social-tom.csv'))
output_all = pd.concat([outputDf_social, outputDf_tom])

quality_output = pd.DataFrame()
for i, par in enumerate(np.unique(output_all['ID'])):
    parName = parInfo.loc[parInfo['numericID']==par, 'parID']
    quality_output.loc[i, 'parID'] = parName.iloc[0]
    quality_output.loc[i, 'numericID'] = par
    for column in ['time_taken', 'exp_question', 'catch_correct', 'catch_accuracy', 'reaction_time_mean', 'reaction_time_median', 'reaction_time_minimum', 'rt_mean_z', 'rt_median_z', 'rt_min_z']:
        a = quality_catch.loc[quality_catch['participant']==parName.iloc[0], column]
        quality_output.loc[i, column] = a.iloc[0]
    dataDf = output_all.loc[output_all['ID']==par]
    quality_output.loc[i, 'int_range'] = max(dataDf['int_rating']) - min(dataDf['int_rating'])
    quality_output.loc[i, 'int_std'] = dataDf['int_rating'].std()
    quality_output.loc[i, 'con_valid_min'] = min(dataDf['con_rating_valid'])
    quality_output.loc[i, 'con_valid_max'] = max(dataDf['con_rating_valid'])
    quality_output.loc[i, 'con_valid_mean'] = dataDf['con_rating_valid'].mean()
    
    idx = list(parList).index(parName.iloc[0])
    quality_output.loc[i, 'frame_rate'] = frameRate[idx]
    quality_output.loc[i, 'sampling_interval'] = sampling_interval[idx]

# save results
quality_output.to_csv(os.path.join(rootDir, 'preprocess', 'outputs', 'quality_all_social-tom_cleaned.csv'), index=False)

#%%
print(parList[-2:-1])