# Quality check on the online data collected
# Author: Zizhuang Miao

# This script extracts quality-related data/metrics from the raw behavior data of the online study (ratings on written narratives)
# Data of interest includes: number of correct answers on catch questions, variabilities in intermittent ratings, sampling rate/frame rate in continuous ratings, and valid data entries in continuous ratings

# %%
# read data
import pandas as pd
import os
import numpy as np
import glob
from scipy.stats import zscore

rootDir = 'C:\\'
dataDir = os.path.join(rootDir, 'data')
parInfo = pd.read_csv(os.path.join(rootDir, 'participants_dataAvail.csv'), encoding='utf-8')
parList = parInfo['parID']
timeTaken = parInfo['time_taken']


# %%
# extract catch trial accuracies
## define columns that are to be stored in the outputs
catchCorrect = []
catchAccuracy = []
threshold = 7.5    # the threshold over which to judge each catch question answer as correct or not
exp_question = []
rt_mean = []    # the average reaction time in intermittent ratings
rt_median = []  # the median ...
rt_min = []     # the minimum ...
sampling_interval = []    
frameRate = []  # different participants have different frame rates, as well as sampling rates in continuous ratings

for par in parList:
    correct = []
    dataFile = os.path.join(dataDir, par+'_group3.1-3.2-social-tom.csv')
    rawData = pd.read_csv(dataFile)
    exp_question.append(rawData.loc[0, 'exp_question'])
    sampling_interval_obtained = False
    frameRate.append(rawData.loc[0, 'frameRate'])
    rt = []    # all reaction time
    for i in range(len(rawData)):
        if not np.isnan(rawData.loc[i, 'catchQuestionSlider.response']):
            response = rawData.loc[i, 'catchQuestionSlider.response']
            trueAnswer = rawData.loc[i, 'catch_value']
            if trueAnswer == 0:
                if response <= threshold:
                    correct.append(1)
                else:
                    correct.append(0)
            elif trueAnswer == 50:
                if response >= 50-threshold and response <= 50+threshold:
                    correct.append(1)
                else:
                    correct.append(0)
            elif trueAnswer == 100:
                if response >= 100-threshold:
                    correct.append(1)
                else:
                    correct.append(0)
        if not np.isnan(rawData.loc[i, 'intQuestionSlider.rt']):
            rt.append(rawData.loc[i, 'intQuestionSlider.rt'])
        if not np.isnan(rawData.loc[i, 'fixation.started']) and not sampling_interval_obtained:
            con = 0
            exec('con = ' + rawData.loc[i, 'continuousResponses_audio'])
            try:
                timeDiff = con[2][1] - con[1][1]
                sampling_interval.append(timeDiff)
                sampling_interval_obtained = True
            except:
                sampling_interval_obtained = False
    catchCorrect.append(correct)
    catchAccuracy.append(sum(correct)/4)
    rt_mean.append(np.mean(rt))
    rt_median.append(np.median(rt))
    rt_min.append(min(rt))

rt_mean_z = zscore(rt_mean)
rt_median_z = zscore(rt_median)
rt_min_z = zscore(rt_min)
outputDf = pd.DataFrame({'participant': parList, 'time_taken': timeTaken, 'exp_question': exp_question, 'catch_correct': catchCorrect, 'catch_accuracy': catchAccuracy, 'reaction_time_mean': rt_mean, 'reaction_time_median': rt_median, 'reaction_time_minimum':rt_min, 'rt_mean_z':rt_mean_z, 'rt_median_z': rt_median_z, 'rt_min_z': rt_min_z, 'sampling_interval': sampling_interval, 'frame_rate':frameRate})
outputDf.to_csv(os.path.join(rootDir, 'preprocess', 'outputs', 'qualityControl_catchQuestion_social-tom.csv'), index=False)

#%%
situationList = list(range(1,10))*8
print(situationList)