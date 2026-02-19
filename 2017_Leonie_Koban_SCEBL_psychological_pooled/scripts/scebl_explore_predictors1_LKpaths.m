basedir = '/Users/Leonie/Dropbox/WORK/PROJECTS/A_2017_Leonie_Koban_SCEBL_psychological_pooled';
cd(basedir)

addpath(fullfile(basedir, 'scripts'));

datadir = fullfile(basedir, 'data');
load(fullfile(datadir, 'SCEBL_pooleddata_canlabdataset.mat'));

list_variables(scebl_pooleddata)

%%
xname = 'LOT_R' ;
yname = 'ExpectancyEffectonPain  (med temp)';

scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%
yname = 'SocialInfluenceEffectonPain (med temp)';
scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%
yname = 'CSEffectonPain  (med temp)' ;
scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%


xname = 'BIS' ;
yname = 'ExpectancyEffectonPain  (med temp)';

scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%
yname = 'SocialInfluenceEffectonPain (med temp)';
scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%
yname = 'CSEffectonPain  (med temp)' ;
scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%


xname = 'BAS_Drive' ;
yname = 'ExpectancyEffectonPain  (med temp)';

scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%
yname = 'SocialInfluenceEffectonPain (med temp)';
scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

%%
yname = 'CSEffectonPain  (med temp)' ;
scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

