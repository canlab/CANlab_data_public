## first-level GLM analysis for social interaction (SInt) within participants' ratings 
# Author: Zizhuang Miao

# Main regressors:
# rating: the rating period in every trial (the last 8 seconds)
# SInt: medians of online participants' ratings of how much social interactions at each moment (-1 is not at all and 1 is definitely yes), mean-centered in each modalities
# all the beta maps from the main regressors were stored

# Nuisance regressors:
# 24 Motion parameters
# csf mean activity
# linear and quadartic trend
# discrete cosine basis functions (as a high pass filter at 1/128 Hz)
# indicator function of each spike (defined as larger than 3 SD from the global mean or from the previous volume)
# indicator function of each time point that should not be modeled (the first few points in each trial)

# BOLD data were first spatial smoothed using a Gaussian kernel (full-width-at-half-maximum = 6mm) before being regressed

import pandas as pd
import numpy as np
import os
import glob
import shutil
import matplotlib.pyplot as plt
import nibabel as nib
from nltools.data import Brain_Data, Design_Matrix
from nltools.stats import regress, zscore
from nltools.file_reader import onsets_to_dm
from nltools.external import glover_hrf
import time      # monitor running time
import psutil    # monitor RAM usage
import sys       # read the command line argument (job ID)

# ---------------------------------------------------------------------
#                            function
# ---------------------------------------------------------------------
def make_motion_covariates(mc, tr):
    '''Create motion covariates regressors from realignment parameters
    Args:
        mc: (pd.DataFrame) realignment parameters
        tr: (float) repetition time
    Returns:
        mcReg: (Design_Matrix) instance that contains all 24 motion covariates
    '''
    z_mc = zscore(mc)
    all_mc = pd.concat([z_mc, z_mc**2, z_mc.diff(), z_mc.diff()**2], axis=1)
    all_mc.fillna(value=0, inplace=True)
    mcReg = Design_Matrix(all_mc, sampling_freq=1/tr)
    mcReg.columns = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'trans_x2',
       'trans_y2', 'trans_z2', 'rot_x2', 'rot_y2', 'rot_z2', 'trans_x3', 'trans_y3',
       'trans_z3', 'rot_x3', 'rot_y3', 'rot_z3', 'trans_x4', 'trans_y4', 'trans_z4',
       'rot_x4', 'rot_y4', 'rot_z4']
    return mcReg

# ---------------------------------------------------------------------
#                 parameters (change per time running)
# ---------------------------------------------------------------------
dataDir = ''    # smoothed preprocessed BOLD data directory
covDir = ''    # fmriprep output directory that contain confound files
eventsDir = 'data/all_design_matrices/designs_social'    # directory that contains design matrix files
outputDir = ''

process = psutil.Process()


##### MORDIFY THE FOLLOWING LINE WITH BATCH JOB ARRANGEMENT AND RUNNING GOALS
jobID = int(sys.argv[1]) - 1
subList = pd.read_csv('cleanedSubjects.csv')['subjects'].tolist()
subList = subList[jobID*8:(jobID+1)*8] if jobID < 11 else subList[jobID*8:]    # 12 jobs, 8 or less subjects per job
#####
runList = ['1', '2', '3', '4']

tr = 0.46    # TR = 0.46s
nDummy = 6   # number of dummy volumes at the start of each run
start_time = time.time()
i = 1

##### BAD RUNS
badruns = pd.read_csv(os.path.join(eventsDir, 'problematic_runs_narratives.csv'))

# ---------------------------------------------------------------------
#                            main loop
# ---------------------------------------------------------------------
for sub in ['sub-0131']:
    print('***********************************', flush=True)
    print(f'Subject #{i} in this job', flush=True)
    allBold = Brain_Data()    # a Brain_Data instance containging the BOLD data of all runs of a single subject 
    allDm = Design_Matrix(sampling_freq = 1/tr)    # a multi-run design matrix
    stats = []
    data = []
    nodataCount = 0    # if all events files were not available, don't do regressions
    
    for run in runList:
        # skip bad runs
        if sub in list(badruns['sub']):
            runInorNot = ('run-'+run == badruns.loc[badruns['sub']==sub, 'run'])
            if runInorNot.any():
                nodataCount += 1
                continue

        # find run length (unit: TR)
        boldFile = os.path.join(dataDir, f'{sub}_ses-02_task-narratives_acq-mb8_run-{run}_space-MNI152NLin2009cAsym_desc-preproc_bold_smoothed-fwhm6mm.nii.gz')
        if not os.path.isfile(boldFile):
            print(f'There is no bold file for {sub} run-{run}!')
            nodataCount += 1
            continue
        numberTr = nib.load(boldFile).shape[-1]    # run length (dummy volumes already excluded before saving this smoothed file)

        # load events file as design matrix
        onsetsFile = os.path.join(eventsDir, f'{sub}_run0{run}.csv')
        if not os.path.isfile(onsetsFile):
            print(f'There is no events file for {sub} run-0{run}!')
            nodataCount += 1
            continue
        
        designs = pd.read_csv(onsetsFile)
        dm = Design_Matrix(designs, sampling_freq=1/tr)
        dmCon = dm.convolve(columns=['SInt_audio', 'rating']) if run in ['1', '2'] else dm.convolve(columns=['SInt_text', 'rating'])
        
        # nuisance regressors
        # motion covariates
        covFile = os.path.join(covDir, sub, 'ses-02', 'func', f'{sub}_ses-02_task-narratives_acq-mb8_run-{run}_desc-confounds_timeseries.tsv')
        cov = pd.read_csv(covFile, sep='\t')
        cov = cov.loc[nDummy:]
        cov = cov.reset_index(drop=True)
        mc = cov[['trans_x','trans_y','trans_z','rot_x', 'rot_y', 'rot_z']]
        mcReg = make_motion_covariates(mc, tr=tr)
        # csf signals
        csf = cov[['csf']]
        csf = Design_Matrix(csf, sampling_freq=1/tr)
        covs = Design_Matrix(pd.concat([mcReg, csf], axis=1), sampling_freq=1/tr)
        # intercept
        covs = covs.add_poly(order=0)
        
        data = Brain_Data(boldFile)
        print(f'{sub} run-0{run}: BOLD data loaded', data.shape())
        spikes = data.find_spikes(global_spike_cutoff=3, diff_spike_cutoff=3)
        spikes = Design_Matrix(spikes.iloc[:,1:], sampling_freq=1/tr)
        # rename the spikes so that they won't be in the same regressors in the all-run dm
        for col in spikes.columns:
            spikes = spikes.rename(columns={col:f'{run}_{col}'})

        # runwise grand mean scaling
        grand_mean = np.mean(data.mean().data)
        data.data = 100 * data.data / grand_mean
        
        allBold = allBold.append(data)
        print(f"Current data size: {allBold.shape()}", flush=True)

        # high-pass filter
        dmCon = dmCon.add_dct_basis(duration=128)
        for col in dmCon.columns:
            if col[:3] == 'cos':
                 dmCon = dmCon.rename(columns={col:f'{run}_{col}'})

        # design matrix containing all event regressors and nuisance regressors
        dmCon = Design_Matrix(pd.concat([dmCon, covs, spikes], axis=1), sampling_freq=1/tr)
        # append the design matrix to all run matrix, create separate columns for the covariates
        allDm = allDm.append(dmCon, axis=0, unique_cols=covs.columns)
    
    print(f'{sub}: All design matrices concatenated')

    # regression
    if nodataCount == 4:
        print(f'No events or BOLD data for {sub}!')
        continue
    allBold.X = allDm
    stats = allBold.regress()

    # save beta maps and standard error maps
    indSocInt_audio = np.where(allDm.columns=='SInt_audio_c0')
    indSocInt_text = np.where(allDm.columns=='SInt_text_c0')
    indRating = np.where(allDm.columns=='rating_c0')
    stats['beta'][indSocInt_audio].write(os.path.join(outputDir, f'{sub}_SInt_audio_beta.nii.gz'))
    stats['beta'][indSocInt_text].write(os.path.join(outputDir, f'{sub}_SInt_text_beta.nii.gz'))
    stats['beta'][indRating].write(os.path.join(outputDir, f'{sub}_rating_beta.nii.gz'))
    
    print(f'{sub}: Writing done!')
    print(f"Current RAM usage: {process.memory_info().rss / 1024 / 1024:.2f} MB")
   
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds, or {elapsed_time/60:.2f} minutes")
    
    i += 1