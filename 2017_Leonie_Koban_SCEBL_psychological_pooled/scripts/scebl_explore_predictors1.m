basedir = '/Users/tor/Dropbox (Personal)/A_2017_Leonie_Koban_SCEBL_psychological_pooled';
basedir = '/Users/tor/Documents/Code_Repositories/CANlab_data_private/2017_Leonie_Koban_SCEBL_psychological_pooled';
datadir = fullfile(basedir, 'data');
resultsdir = fullfile(basedir, 'results');

load(fullfile(datadir, 'SCEBL_pooleddata_canlabdataset.mat'));

cd(basedir)

addpath(fullfile(basedir, 'scripts'));

list_variables(scebl_pooleddata)

print_summary(scebl_pooleddata)

%% list Ns by study version

wh_expt = get_var(scebl_pooleddata, 'InstructionVersion');
[indic, x] = condf2indic(wh_expt);

N = sum(indic)';


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

