basedir = '/Users/tor/Documents/Code_Repositories/CANlab_data_private/2017_Leonie_Koban_SCEBL_psychological_pooled';
datadir = fullfile(basedir, 'data');
resultsdir = fullfile(basedir, 'results');

load(fullfile(datadir, 'SCEBL_pooleddata_canlabdataset.mat'));

%% Fix a couple of things

n = size(scebl_pooleddata.Subj_Level.data, 1);
scebl_pooleddata.Subj_Level.id = {};
for i = 1:n
    scebl_pooleddata.Subj_Level.id{i, 1} = num2str(i);
end

scebl_pooleddata.Subj_Level.type = repmat({'numeric'}, 1, size(scebl_pooleddata.Subj_Level.data, 2));

% print summary

print_summary(scebl_pooleddata)

% save text output
write_text(scebl_pooleddata)

% save mat file with changes
save SCEBL_pooleddata_canlabdataset scebl_pooleddata
