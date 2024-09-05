clc;clear 
close all

%% parameters
   
for laplacian    = 1;
for window       = 1; % 1 = 300 ms, 2 = 500 ms
for chanclusters = 1; % 0 = PO7/8,  1 = PO7/8 + neighbouring channels.


for pp           = [1:23];    

plotResults      = 0;
plotResults_topo = 0;
plotResults_cvsi = 0;
colormap2use     = fliplr(brewermap(100, 'RdBu'));

%% load data
param = getSubjParam(pp);

load([param.pathnewlen, '/processed_data/' 'epoched_data_eeg' '__', param.subjName], 'data');

%% keep only channels of interest
cfg         = [];
cfg.channel = {'EEG'};
data        = ft_preprocessing(cfg, data);

%% remove bad ICA components
load([param.pathnewlen, '/saved_data/' 'ICAcomponents', '__' param.subjName], 'comp2rem','ica');
cfg             = [];
cfg.component   = comp2rem;
data            = ft_rejectcomponent(cfg, ica, data);

%% remove bad trials
load([param.pathnewlen, '/processed_data/' 'usableTrials', '__' param.subjName], 'trl2keep');

cfg         = [];
cfg.trials  = trl2keep;
data        = ft_selectdata(cfg, data);

%% surface laplacian?
if laplacian
    cfg         = [];
    cfg.elec    = ft_read_sens('standard_1020.elc');
    data        = ft_scalpcurrentdensity(cfg, data);
end

%% get time-frequency response
if window == 1      windowsize = 0.3; % 300 ms sliding time window
elseif window == 2  windowsize = 0.5; % 300 ms sliding time window
end

cfg             = [];
cfg.method      = 'mtmconvol';
cfg.keeptrials  = 'yes';
cfg.taper       = 'hanning';
cfg.foi         = [2:1:40];
cfg.pad         = 10; 
cfg.toi         = [data.time{1}(1) + (windowsize / 2) : .01 : data.time{1}(end) - (windowsize / 2)]; % steps of 10 ms. 
cfg.t_ftimwin   = ones(1,length(cfg.foi))*windowsize;
tfr             = ft_freqanalysis(cfg, data);
clear data


%% split into three conditions
imper100 = ismember(tfr.trialinfo(:,2), 4);
high100  = ismember(tfr.trialinfo(:,2), 1);
med80    = ismember(tfr.trialinfo(:,2), 2);
low60    = ismember(tfr.trialinfo(:,2), 3);

valid = ismember(tfr.trialinfo(:,3), 1);
invalid  = ismember(tfr.trialinfo(:,3), 2);

% reverse cueside for invalid trials
% left     = ismember(tfr.trialinfo(:,1), [41]);
% right    = ismember(tfr.trialinfo(:,1), [42]);
trialinfo = tfr.trialinfo;
for t = 1:length(trialinfo)
    if trialinfo(t,3) == 1 % valid
        trialinfo(t,4) = trialinfo(t,1)-40; % valid
    elseif trialinfo(t,3) == 2 % invalid
        trialinfo(t,4) = 3-(trialinfo(t,1)-40); % valid
    end
end
tfr.trialinfo = trialinfo; clear trialinfo

left     = ismember(tfr.trialinfo(:,4), [1]);
right    = ismember(tfr.trialinfo(:,4), [2]);



%% time-freq for condition comparisons + put into a structure we can later plot
timefreq = [];
timefreq.time = tfr.time*1000;
timefreq.freq = tfr.freq;
timefreq.label = tfr.label;
timefreq.dimord = 'chan_freq_time';

timefreq.imper100 = squeeze(mean(tfr.powspctrm(imper100,:,:,:)));
timefreq.high100 = squeeze(mean(tfr.powspctrm(high100,:,:,:)));
timefreq.med80 = squeeze(mean(tfr.powspctrm(med80,:,:,:)));
timefreq.low60 = squeeze(mean(tfr.powspctrm(low60,:,:,:)));

timefreq.imper_vs_high      = ((timefreq.imper100 - timefreq.high100) ./ (timefreq.imper100 + timefreq.high100)) * 100;
timefreq.high_vs_low        = ((timefreq.high100  - timefreq.low60) ./ (timefreq.high100  + timefreq.low60))   * 100;


%% get contra-ipsi--for all channels--averaged over trials
chL = []; chR = [];
chL = match_str(tfr.label, {'Fp1','AF3','AF7','F1','F3','F5','F7','FC1','FC3','FC5','FT7','C1','C3','C5','T7','CP1','CP3','CP5','TP7','P1','P3','P5','P7','P9','PO3','PO7','O1'});
chR = match_str(tfr.label, {'Fp2','AF4','AF8','F2','F4','F6','F8','FC2','FC4','FC6','FT8','C2','C4','C6','T8','CP2','CP4','CP6','TP8','P2','P4','P6','P8','P10','PO4','PO8','O2'});

data_chL = tfr.powspctrm(:,chL,:,:);
data_chR = tfr.powspctrm(:,chR,:,:);

ci_allchan = [];
ci_allchan.time = timefreq.time;
ci_allchan.freq = timefreq.freq;
ci_allchan.dimord = 'time_freq_chan';
ci_allchan.label = {'imper100-cvsi','high100-cvsi','med80-cvsi','low60-cvsi'};

sel = imper100;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(1,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = high100;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(2,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = med80;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(3,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = low60;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(4,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;


ci_allchan.data(5,:,:,:) = ci_allchan.data(1,:,:,:) - ci_allchan.data(2,:,:,:); %im100-infor100
ci_allchan.data(6,:,:,:) = ci_allchan.data(2,:,:,:) - ci_allchan.data(4,:,:,:); %infor100-infor60

clear sel contra ipsi


%% get contra - ipsi--averaged over channel & trials
chL = []; chR = [];
chL = match_str(tfr.label, 'PO7'); 
% chL = [chL, chL];
chR = match_str(tfr.label, 'PO8');
% chR = [chR, chR];

% if chanclusters
% chL = match_str(tfr.label, {'PO7','PO3','P7','P5','P3','O1'});
% chR = match_str(tfr.label, {'PO8','PO4','P8','P6','P4','O2'});
% end

if chanclusters
chL = []; chR = [];
chL = match_str(tfr.label, {'O1','PO7','PO3','P9','P7','P5','P3','P1'});
chR = match_str(tfr.label, {'O2','PO8','PO4','P10','P8','P6','P4','P2'});
end

data_chL = squeeze(mean(tfr.powspctrm(:,chL,:,:),2));
data_chR = squeeze(mean(tfr.powspctrm(:,chR,:,:),2));

ci = [];
ci.time = timefreq.time;
ci.freq = timefreq.freq;
ci.dimord = 'time_freq_chan';
ci.label = {'imper100-cvsi','high100-cvsi','med80-cvsi','low60-cvsi'};

sel = imper100;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(1,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(1,:,:) = contra; cii.data(1,:,:) = ipsi; clear contra ipsi

sel = high100;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(2,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(2,:,:) = contra; cii.data(2,:,:) = ipsi; clear contra ipsi

sel = med80;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(3,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(3,:,:) = contra; cii.data(3,:,:) = ipsi; clear contra ipsi

sel = low60;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(4,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(4,:,:) = contra; cii.data(4,:,:) = ipsi; clear contra ipsi

ci.data(5,:,:) = ci.data(1,:,:) - ci.data(2,:,:); %im100-infor100
ci.data(6,:,:) = ci.data(2,:,:) - ci.data(4,:,:); %infor100-infor60



%% save
if laplacian toadd1 = '_laplacian'; else toadd1 = ''; end
if window == 1 toadd2 = ''; elseif window == 2 toadd2 = '_500mswindow'; end
if chanclusters toadd3 = '_chanclusters'; else toadd3 = ''; end

save([param.pathnewlen, 'saved_data/cue_cor_timefreq', toadd1,toadd2,toadd3, '__' param.subjName], 'timefreq','ci','cic','cii','ci_allchan');


end
end
end
end
