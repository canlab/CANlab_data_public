# Generate time course of manual annotations
# Zizhuang Miao

# This script can resample the time course of experimenters' annotations for social interactions and theory of mind demands into the same resolution as the online participants' ratings -- 230 ms
# It allows for direct comparisons between them

import pandas as pd
import numpy as np
import os

rootDir = 'C:\\'
dataDir = 'C:\\'
outputDir = 'C:\\'

stimuliDuration = pd.read_csv(os.path.join(rootDir, 'derivedData', 'durationInfo.csv'))

for n in range(1,9):
    for s in range(1,10):
        annotation = pd.read_csv(os.path.join(dataDir, f"Narrative{n}Situation{s}_transcript_social_GPT.csv"))
        # add time points into the annotation to match the stimuli used in the online studies
        if n >= 5:    # audio runs
            # add 1.6 seconds at the end
            annotation = pd.concat([annotation, pd.DataFrame({'onset': [annotation.iloc[-1, 0]+1.6], 'socialInteraction': [annotation.iloc[-1, 3]]})])
            # and 0.3 seconds at the start
            annotation['onset'] += 0.3
            annotation = pd.concat([pd.DataFrame({'onset': [0], 'socialInteraction': [annotation.iloc[0, 3]]}), annotation])
            annotation = annotation.reset_index()
        else:    # text runs
            # add 1.5 seconds at the end
            annotation = pd.concat([annotation, pd.DataFrame({'onset': [annotation.iloc[-1, 0]+1.5], 'socialInteraction': [annotation.iloc[-1, 3]]})])
            annotation = annotation.reset_index()
        maxColumn = stimuliDuration.loc[(stimuliDuration.narrative==n) & (stimuliDuration.situation==s), 'maximum_column'].iloc[0]
        
        outputDic = {'time': [], 'rating': []}
        idx = 0
        for t in range(0, int(maxColumn), 230):
            outputDic['time'].append(round(t/1000, 2))
            if idx == 0:
                if t >= annotation.loc[0, 'onset']*1000 - 115:   # the first viable rating is not interpolated
                    outputDic['rating'].append(annotation.loc[0, 'socialInteraction']*100)
                    idx += 1
            elif idx < len(annotation):    # the intermediate ratings are linearly interpolated 
                if t >= annotation.loc[idx-1, 'onset']*1000 :
                    while not (annotation.loc[idx-1, 'onset']*1000 <= t and annotation.loc[idx, 'onset']*1000 > t):
                        idx += 1
                        if idx >= len(annotation):
                            break
                    if idx < len(annotation):
                        outputDic['rating'].append((annotation.loc[idx-1, 'socialInteraction'] + (annotation.loc[idx, 'socialInteraction'] - annotation.loc[idx-1, 'socialInteraction'])*(t/1000-annotation.loc[idx-1, 'onset'])/(annotation.loc[idx, 'onset'] - annotation.loc[idx-1, 'onset']))*100)
                    else:
                        outputDic['rating'].append(annotation.loc[idx-1, 'socialInteraction']*100)    # the last rating is not interpolated
                else:
                    continue
            elif idx == len(annotation):
                continue
        if len(outputDic['time']) < len(outputDic['rating']):
            while not len(outputDic['time']) == len(outputDic['rating']):
                outputDic['time'].append(outputDic['time'][-1] + 0.23)
        if len(outputDic['time']) > len(outputDic['rating']):
            while not len(outputDic['time']) == len(outputDic['rating']):
                outputDic['rating'].append(outputDic['rating'][-1])

        outputDf = pd.DataFrame(outputDic)
        outputDf.to_csv(os.path.join(outputDir, f'narrative{n}situation{s}.csv'), index=False)
