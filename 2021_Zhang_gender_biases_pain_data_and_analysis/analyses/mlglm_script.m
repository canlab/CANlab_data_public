basedir = '/Users/yoni/Dropbox/empathy for pain/ml_glm_analysis';
cd(basedir);

PSPI = readtable('PSPI at stimuli level.xlsx');
PSPI.Stimulusbl = PSPI.VideoClips; % so the join works below

SelfReportedPain = readtable('target self-reports.xlsx');
SelfReportedPain.VideoClips = SelfReportedPain.BaselineSTIMULI; % to make the join work

% DROP VARS I DON'T NEED NOW, MAY NEED LATER THO
stimtabl = join(PSPI, SelfReportedPain);
stimtabl.BaselineSTIMULI = [];
stimtabl.Genderbl = [];
stimtabl.IntensityBl = [];
stimtabl.frame = [];
stimtabl.VideoClips = [];
stimtabl.intensity = [];

cd rating


%% 

% corr of self-reported pain w/ pain judgements
corr([stimtabl.OPR_0_5_ stimtabl.AFF_0_16_ stimtabl.VAS_0_10_ stimtabl.SEN_0_16_ stimtabl.PSPIaverage stimtabl.PSPIpeak])

% could do a mediation:  self report -> PSPI -> observed pain, w/ gender
% moderator on paths



%% set up cell arrays for ml glm

subjs = filenames('empathy*1.txt');

subjs([10 40 41]) = [];  % SOMETHING IS OFF W/ NUM 10 and 40 41

X={};Y={};Gender={};
for i=1:length(subjs)

    tabl = readtable(subjs{i});
    tabl = tabl(:, {'Subject', 'TrialNum', 'Stimulusbl', 'Intensitybl', 'Genderbl', 'VASratingbl'}); % just to get relevant cols, to simplify. unnecesary.

    % BEST TO IGNORE INTENSITY FOR NOW, WAS CHOSEN SEMI-ARBITRARILY
    
    % for now, just look at high intensity
    %tabl = tabl( char(tabl.Intensitybl)=='M' ,:);
  %{  
    % convert intensity to 1-2-3
    low_inds = categorical(tabl.Intensitybl)=='L';
    med_inds = categorical(tabl.Intensitybl)=='M';
    high_inds = categorical(tabl.Intensitybl)=='H';
    tabl.Intensity_num = sum([low_inds med_inds*2 high_inds*3], 2);
    %}
    tabl = join(tabl, stimtabl);

    myX = [];
    myX(:,1) = grp2idx(categorical(tabl.Genderbl)); % 1 is female, 2 is male 
    myX(:,2) = tabl.VAS_0_10_;
 %   myX(:,3) = tabl.PSPIpeak;
    X{i} = myX;
    
    Gender{i} = grp2idx(categorical(tabl.Genderbl)); % 1 is female, 2 is male 
    PSPIs{i} = tabl.PSPIpeak; 
    observed{i} = tabl.VASratingbl;
    VAS{i} = tabl.VAS_0_10_;
    
    males = categorical(tabl.Genderbl)=='Male';
    acc(i, 1) = corr(tabl.VASratingbl(~males), tabl.VAS_0_10_(~males));
    acc(i, 2) = corr(tabl.VASratingbl(males), tabl.VAS_0_10_(males));
    mae(i, 1) = mean(abs(tabl.VASratingbl(~males) - tabl.VAS_0_10_(~males)));
    mae(i, 2) = mean(abs(tabl.VASratingbl(males) - tabl.VAS_0_10_(males)));    
    % also add PSPI for this targ

    Y{i} = tabl.VASratingbl;
end

%% empathic accuracy across subjects, divided by gender

figure; violinplot(acc, 'xlabel', {'Female targs' 'Male targs'}); title('corr of self-reported and observed pain');
figure; violinplot(mae, 'xlabel', {'Female targs' 'Male targs'}); title('MAE of self-reported and observed pain');
%figure; plot(acc')

[~,p]=ttest(acc(:,1)-acc(:,2))

%% the ml glm

% is target gender predicting of pain ratings, controlling for
% self-reported pain (VAS)?
glmfit_multilevel(Y, X, [], 'names', {'Intrcpt' 'Gender' 'TargVAS'});
% 
% 2nd-level B01
% 	Intrcpt	Gender	Targ report	
% Coeff	15.01	-2.09	3.61	
% STE	1.79	0.82	0.20	
% t	8.38	-2.56	18.52	
% Z	6.39	-2.20	8.13	
% p	0.0000	0.0139	0.0000	

%% is target gender predicting of pain ratings, controlling for
% self-reported pain (VAS)?  and when we also control for PSPI peak?
glmfit_multilevel(Y, X, [], 'names', {'Intrcpt' 'Gender' 'TargVAS' 'PSPIpeak'});
% % 
% 2nd-level B01
% 	Intrcpt	Gender	TargVAS	PSPIpeak	
% Coeff	9.13	1.99	2.05	2.49	
% STE	1.78	0.82	0.15	0.15	
% t	5.13	2.43	13.68	16.73	
% Z	4.39	2.07	8.13	8.13	
% p	0.0000	0.0191	0.0000	0.0000	

%% natural next question:  is the effect of target gender on pain ratings mediated by PSPI, controlling for self-report?

%% multilevel mediation

[paths, toplevelstats, firstlevelstats] = mediation(Gender, Y, PSPIs, 'plots', 'names', {'Gender' 'Pain ratings' 'PSPI'}, 'bootstrapfirst', 'bootsamples', 10000, 'verbose');

%{
for HIGH intensity, using PSPI average

Multi-level model
	a	b	c'	c	ab	
Coeff	-0.83	2.09	-0.85	-2.72	-0.92	
STE	0.01	0.44	1.31	1.21	0.30	
t (~N)	-83.53	4.71	-0.64	-2.25	-3.03	
Z	-8.21	4.20	-0.64	-2.17	-2.86	
p	0.0000	0.0000	0.5225	0.0300	0.0043	



for M intensity, using PSPI peak -- tho code threw a lot of warnings!
________________________________________

Multi-level model
	a	b	c'	c	ab	
Coeff	-1.24	-0.37	2.27	2.71	0.29	
STE	0.05	0.35	1.06	0.98	0.36	
t (~N)	-26.78	-1.07	2.15	2.77	0.80	
Z	-8.21	-1.05	2.08	2.64	0.80	
p	0.0000	0.2920	0.0377	0.0083	0.4261	


%}

%% are women more expressive of the same pain?

create_figure('gend differences in expression of pain');
males = categorical(stimtabl.gender)=='Male';
scatter(stimtabl.VAS_0_10_(males), stimtabl.PSPIpeak(males), 50,'r'); h=lsline; 
scatter(stimtabl.VAS_0_10_(~males), stimtabl.PSPIpeak(~males), 50,'b'); h=lsline;
xlabel('targ self report'), ylabel('PSPI peak')
h(1).Color = 'r'; h(2).Color = 'b';

create_figure('PSPI by gender'); violinplot([stimtabl.PSPIpeak(~males) stimtabl.PSPIpeak(males)],'xlabel', {'Female targs' 'Male targs'}); 
title('PSPI peak')

%% are women more expressive of the same pain CONTROLLING FOR SELF-REPORTED PAIN?

[B,DEV,STATS] = glmfit(stimtabl.VAS_0_10_, stimtabl.PSPIpeak);  % regresss VAS out of PSPI
stimtabl.PSPIresidVAS_0_10_ = STATS.resid;

create_figure('PSPI by gender resid self-report'); violinplot([stimtabl.PSPIresidVAS_0_10_(~males) stimtabl.PSPIresidVAS_0_10_(males)],'xlabel', {'Female targs' 'Male targs'}); 
title('PSPI peak')

[~,p]=ttest2(stimtabl.PSPIresidVAS_0_10_(~males), stimtabl.PSPIresidVAS_0_10_(males))