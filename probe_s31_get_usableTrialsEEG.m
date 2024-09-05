clc;clear
close all;

%% parameters
sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'29';'30';'31';'32';'33';'34';};

% for pp = 1
for pp = 1:length(sub_list)
% for pp = [6 15]
    % may need to exclude pp6 & pp15 for further data analysis
    

param = getSubjParam(pp);
 
%% load data
load([param.pathnewlen, '/processed_data/' 'epoched_data_eeg' '__', param.subjName], 'data');

%% keep only channels of interest
cfg = [];
cfg.channel = {'EEG'};
data = ft_preprocessing(cfg, data);

%% remove bad ICA components
load([param.pathnewlen, '/saved_data/' 'ICAcomponents', '__' param.subjName], 'comp2rem','ica');
cfg           = [];
cfg.component = comp2rem;
data          = ft_rejectcomponent(cfg, ica, data);

%% find bad trials
data.trialinfo(:,end+1) = 1:length(data.trial);
trials_old = data.trialinfo(:,end);

% % browse the data first
% cfg = [];
% cfg.viewmode = 'vertical';
% databrow = ft_databrowser(cfg, data);

cfg = [];
cfg.method = 'summary';
cfg.channel = {'EEG'};
data = ft_rejectvisual(cfg, data);
databrow = ft_databrowser(cfg, data);

cfg.preproc.bpfilter = 'yes';
cfg.preproc.bpfreq = [4 7];
data = ft_rejectvisual(cfg, data);

cfg.preproc.bpfilter = 'yes';
cfg.preproc.bpfreq = [8 30];
data = ft_rejectvisual(cfg, data);
trials_new = data.trialinfo(:,end);

% metrics
trl2keep = ismember(trials_old, trials_new);

propkeep(pp) = mean(trl2keep)

%% save
save([param.pathnewlen, '/processed_data/' 'usableTrials', '__' param.subjName], 'trl2keep');

end

