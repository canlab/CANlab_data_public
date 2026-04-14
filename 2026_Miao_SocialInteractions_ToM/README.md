# Scripts and behavioral data for *Common and distinct neural correlates of social interaction processing and theory of mind in narratives*

Zizhuang Miao

[![DOI](https://zenodo.org/badge/902574045.svg)](https://doi.org/10.5281/zenodo.18315908)

This repo hosts the codes and behavioral data accompanying the article titled "Common and distinct neural correlates of social interaction processing and theory of mind in narratives". Readers should be able to use and adpat these codes to replicate experiments, analyses, and results.

## Software dependency
The analyses were performed in Python 3.9 and MATLAB R2022b. Besides many commonly used packages (e.g., `pandas` and `numpy` in Python), I used the following open-source packages:
+ [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
+ Several packages in [CANLab toolbox](https://github.com/canlab)
+ [NLTools](https://nltools.org/)
+ [Nilearn](https://nilearn.github.io/stable/index.html)

Additional packages needed will be mentioned where necessary.

## Script navigation
Below is a high-level overview of the naming and functions of the scripts in this repo:

+ `stimuli`: The audios and texts of narratives used in the experiments; each .wav and .txt file was used in one trial. Additionally, `all_narratives` contains the text documents for each of the eight narratives. Note that there are two more audio files than text files, which were used in the practice trials of the online experiment.
+ `exp_online_*`: The scripts, condition files, additional stimuli, and dependencies (PsychJS) needed to run the online experiment on Pavlovia.org. Note that both .psyexp files in PsychoPy and .js files are provided, but there are slight differences in how the same functions are implemented because of the innate differences between PsychoPy and JavaScript; the .js files were used for data collection.
+ `exp_fmri`: The scripts and condition file used to run the experiment during fMRI scanning. Note that Psychtoolbox (MATLAB) is needed to run the scripts.
+ `analysis_behavior`: Scripts used to analyze online participants' continuous rating data. 
  + `preprocessing`: Scripts for organizing and cleaning continuous rating data. They should be run in the sequence of the numbers in file names.
  + `analyzing`: Scripts for calculating metrics for online particpiants' continous ratings and generating reports. Figure 2 in the paper can be reproduced using those scripts.
+ `analysis_fmri`: Scripts used to analyze and visualize fMRI data
  + `preprocessing`: Preprocessing of fMRI data was done using a Docker image of fmriprep on a computing cluster. Specific parameters for fmriprep can be read from the shell script.
  + `first_level_GLM`: The scripts used to build design matrices and run first-level GLM. If you would like to reproduce first-level GLM analysis, all you need is the .py files, subject list (`cleanedSubjects.csv`), and the design matrices (detailed in the next section). We provide scripts for two GLM models (only social interactions, and both social interactions and theory of mind) as examples; other models were run by simple variants of these scripts.
  + `second_level_GLM`: Scripts in this folder take the beta maps from first-level GLMs as inputs, and generate all main results in the paper. Several types of analyses were performed:
    + Group level t-tests. One example t test analysis was provided in `social_tom.mlx`, based on SPM and CANLab toolbox. The core functions were `ttest()` on `fmri_data` objects. It can be easily adapted to other t tests in the paper. In the same script, we also provide codes that compare t maps to resting-state functional networks.
    + Bayes factor analysis. It was done in `bayes_factor.mlx`. Related to Figures 6 and S14.
    + Mask-based analysis: extracting beta values within brain masks and performing analysis on the group level. `mask_based_analysis.mlx` contains codes to replacate results in Figures 3b, 4b, 5c, 6b, S4, S5, and S14.
  + `second_level_GLM/brain_visualization.mlx` provides example codes that show volumetric data on the brain surface. Related to Figures 3-6, S2-S11, and S13.

## Data navigation
Two types of data are shared:
+ Ratings and annotations of the narratives (in `ratings_*` folders and `all_ratings_annotations`). In each of these files, ratings or annotations were resampled to the sample time points of the online studies for continuous ratings. The narrative word closest to each time point were included in the "word" column. For the three continuous ratings, column `_ratings_mad` indicates median absolute deviations, `_validPar_trial` indicates the number of participants included in that trial, `_validPar_time` indicates the number of participants with valid rating at that time point, and `_included` indicates whether that rating should be included in analyzing fMRI data (1 = included, 0 = not included).
+ Design matrices for performing first-level GLMs (in `all_design_matrices`).