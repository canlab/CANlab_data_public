## Reorganize continuous ratings
# Author: Zizhuang Miao

# This script is used to reorganize continuous rating data separately for each trial, which could then be used to obtain group summary statistics which will be applied to fMRI data.

# The script will be dedicated to processing social interaction ratings. However, the same methods (and codes) could be easily used for theory of mind ratings.

# We generate long-type data files for each trial.

# get long type data
import pandas as pd
import numpy as np
import os
import codecs

rootDir = 'C:\\'
rawData = pd.read_csv(os.path.join(rootDir, 'derivedData', 'social_long.csv'))

narrativeList = list(range(1,9))
situationList = list(range(1,10))

stimuliDuration = pd.read_csv(os.path.join(rootDir, 'derivedData', 'durationInfo.csv'))

for n in narrativeList:
    for s in situationList:
        dat = rawData.loc[(rawData['narrative']==n) & (rawData['situation']==s)]
        dat = dat.reset_index()
        
        maxColumn = stimuliDuration.loc[(stimuliDuration.narrative==n) & (stimuliDuration.situation==s), 'maximum_column'].iloc[0]
        outputDic = {'ID': []}   # a dictionary storing outputs

        for t in range(0, int(maxColumn), 230):    
            outputDic[str(t)] = []    # create keys

        for i in range(len(dat)):
            if dat.loc[i, 'con_rating_valid'] < 0.4:
                continue
            outputDic['ID'].append(dat.loc[i, 'ID'])
            ratings = []
            exec('ratings = ' + dat.loc[i, 'con_rating'])

            # -- for loop below is an implementation of linear interpolation method --
            idx = 0    # which ratings to record now
            for t in range(0, int(maxColumn), 230):
                outputDic[str(t)].append(np.nan)
                if idx == 0:
                    if t >= ratings[0][1] - 115:   # the first viable rating is not interpolated
                        outputDic[str(t)][-1] = (ratings[idx][0]) if ratings[idx][0]!=999 else np.nan
                        idx += 1
                elif idx < len(ratings):    # the intermediate ratings are linearly interpolated 
                    if t >= ratings[idx-1][1]:
                        while not (ratings[idx-1][1] <= t and ratings[idx][1] > t):
                            idx += 1
                            if idx >= len(ratings):
                                break
                        if idx < len(ratings):
                            outputDic[str(t)][-1] = (ratings[idx-1][0] + (ratings[idx][0] - ratings[idx-1][0])*(t-ratings[idx-1][1])/(ratings[idx][1] - ratings[idx-1][1])) if ratings[idx-1][0]!=999 else ratings[idx][0]
                        else:
                            outputDic[str(t)][-1] = (ratings[idx-1][0])    # the last rating is not interpolated
                    else:
                        continue
                elif idx == len(ratings):
                    continue
        
        # generate a DataFrame for this trial
        outputDf = pd.DataFrame(outputDic)
        outputDf = outputDf.replace(999, np.nan)
        outputDf = outputDf.melt(id_vars='ID', var_name='time', value_name='rating')
        outputDf['time'] = round(outputDf['time'].astype(int)/1000, 2)
        # store the dataFrame
        outputDf.to_csv(os.path.join(rootDir, 'derivedData', 'continuous_social_1a2a3a4a_long', f"narrative{n}situation{s}.csv"), index=False)