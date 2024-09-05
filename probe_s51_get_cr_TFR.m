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

left     = ismember(tfr.trialinfo(:,1), [41]);
right    = ismember(tfr.trialinfo(:,1), [42]);

valid    = ismember(tfr.trialinfo(:,3), 1);
invalid  = ismember(tfr.trialinfo(:,3), 2);


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

timefreq.valid = squeeze(mean(tfr.powspctrm(valid,:,:,:)));
timefreq.invalid = squeeze(mean(tfr.powspctrm(invalid,:,:,:)));

timefreq.validim100 = squeeze(mean(tfr.powspctrm(valid&imper100,:,:,:)));
timefreq.validinfor100 = squeeze(mean(tfr.powspctrm(valid&high100,:,:,:)));
timefreq.validinfor80 = squeeze(mean(tfr.powspctrm(valid&med80,:,:,:)));
timefreq.validinfor60 = squeeze(mean(tfr.powspctrm(valid&low60,:,:,:)));
timefreq.validinfor40 = squeeze(mean(tfr.powspctrm(invalid&low60,:,:,:)));
timefreq.validinfor20 = squeeze(mean(tfr.powspctrm(invalid&med80,:,:,:)));

timefreq.imper_vs_high      = ((timefreq.imper100 - timefreq.high100) ./ (timefreq.imper100 + timefreq.high100)) * 100;
timefreq.high_vs_low        = ((timefreq.high100  - timefreq.low60) ./ (timefreq.high100  + timefreq.low60))   * 100;
timefreq.valid_vs_invalid   = ((timefreq.valid    - timefreq.invalid) ./ (timefreq.valid    + timefreq.invalid)) * 100;


%% multiplot
if plotResults
cfg = [];
cfg.figure = 'gcf';
cfg.layout = 'easycapM1.mat';
cfg.style = 'straight';
cfg.marker = 'off';
cfg.zlim = [-20 20];

% figure; cfg.parameter = 'imper_vs_high'; ft_multiplotTFR(cfg, timefreq); title(cfg.parameter);
% figure; cfg.parameter = 'high_vs_low';   ft_multiplotTFR(cfg, timefreq); title(cfg.parameter);

drawnow;
end

%% singleplot

if plotResults
cfg = [];
cfg.figure = 'gcf';
cfg.zlim = [-30 30];

figure;
cfg.channel = 'AFz';
subplot(4,2,1); cfg.parameter = 'imper_vs_high'; ft_singleplotTFR(cfg, timefreq);  title(['imperative-vs-high', ' - AFz']);
subplot(4,2,2); cfg.parameter = 'high_vs_low';   ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - AFz']);
cfg.channel = 'Pz';
subplot(4,2,3); cfg.parameter = 'imper_vs_high'; ft_singleplotTFR(cfg, timefreq);  title(['imperative-vs-high', ' - Pz']);
subplot(4,2,4); cfg.parameter = 'high_vs_low';   ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - Pz']);
cfg.channel = 'C3';
subplot(4,2,5); cfg.parameter = 'imper_vs_high'; ft_singleplotTFR(cfg, timefreq);  title(['imperative-vs-high', ' - C3']);
subplot(4,2,6); cfg.parameter = 'high_vs_low';   ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - C3']);
cfg.channel = 'Oz';
subplot(4,2,7); cfg.parameter = 'imper_vs_high'; ft_singleplotTFR(cfg, timefreq);  title(['imperative-vs-high', ' - Oz']);
subplot(4,2,8); cfg.parameter = 'high_vs_low';   ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - Oz']);
colormap(colormap2use);
end

%% running topos
if plotResults_topo
cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.comment = 'no';
cfg.style = 'straight';
cfg.marker = 'off';
cfg.zlim = [-30 30];

stepsize = 250; times1 = [0:stepsize:2000]; n = length(times1);

cfg.parameter = 'imper_vs_high';

figure; 
for t = 1:n-1
    cfg.xlim = times1(t:t+1);
    subplot(3,n-1,t);            cfg.ylim = [4 7];   ft_topoplotER(cfg, timefreq); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
    subplot(3,n-1,t+n-1);        cfg.ylim = [8 12];  ft_topoplotER(cfg, timefreq);
    subplot(3,n-1,t+(n-1)*2);    cfg.ylim = [13 30]; ft_topoplotER(cfg, timefreq);
end
colormap(colormap2use);

cfg.parameter = 'high_vs_low'; 

figure; 
for t = 1:n-1
    cfg.xlim = times1(t:t+1);
    subplot(3,n-1,t);            cfg.ylim = [4 7];   ft_topoplotER(cfg, timefreq); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
    subplot(3,n-1,t+n-1);        cfg.ylim = [8 12];  ft_topoplotER(cfg, timefreq);
    subplot(3,n-1,t+(n-1)*2);    cfg.ylim = [13 30]; ft_topoplotER(cfg, timefreq);
end
colormap(colormap2use);
end


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
ci_allchan.label = {'imper100-cvsi','high100-cvsi','med80-cvsi','low60-cvsi',...
            'valid-cvsi','invalid-cvsi',...
            'validim100-cvsi','validinfor100-cvsi','validinfor80-cvsi','validinfor60-cvsi','validinfor40-cvsi','validinfor20-cvsi',...
             'imper-vs-high','high-vs-low','valid-vs-invalid'};

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

sel = valid;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(5,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = invalid;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(6,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = valid&imper100;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(7,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = valid&high100;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(8,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = valid&med80;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(9,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = valid&low60;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(10,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = invalid&low60;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(11,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;

sel = invalid&med80;
contra = squeeze(   mean(data_chR(sel&left,:,:,:)) + mean(data_chL(sel&right,:,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:,:)) + mean(data_chR(sel&right,:,:,:))  ) ./ 2;
ci_allchan.data(12,:,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;


ci_allchan.data(13,:,:,:) = ci_allchan.data(1,:,:,:) - ci_allchan.data(2,:,:,:); %im100-infor100
ci_allchan.data(14,:,:,:) = ci_allchan.data(2,:,:,:) - ci_allchan.data(4,:,:,:); %infor100-infor60
ci_allchan.data(15,:,:,:) = ci_allchan.data(5,:,:,:) - ci_allchan.data(6,:,:,:); %valid-invalid

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
ci.label = {'imper100-cvsi','high100-cvsi','med80-cvsi','low60-cvsi',...
            'valid-cvsi','invalid-cvsi',...
            'validim100-cvsi','validinfor100-cvsi','validinfor80-cvsi','validinfor60-cvsi','validinfor40-cvsi','validinfor20-cvsi',...
             'imper-vs-high','high-vs-low','valid-vs-invalid'};

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

sel = valid;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(5,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(5,:,:) = contra; cii.data(5,:,:) = ipsi; clear contra ipsi

sel = invalid;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(6,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(6,:,:) = contra; cii.data(6,:,:) = ipsi; clear contra ipsi

sel = valid&imper100;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(7,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(7,:,:) = contra; cii.data(7,:,:) = ipsi; clear contra ipsi

sel = valid&high100;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(8,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(8,:,:) = contra; cii.data(8,:,:) = ipsi; clear contra ipsi

sel = valid&med80;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(9,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(9,:,:) = contra; cii.data(9,:,:) = ipsi; clear contra ipsi

sel = valid&low60;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(10,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(10,:,:) = contra; cii.data(10,:,:) = ipsi; clear contra ipsi

sel = invalid&low60;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(11,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(11,:,:) = contra; cii.data(11,:,:) = ipsi; clear contra ipsi

sel = invalid&med80;
contra = squeeze(   mean(data_chR(sel&left,:,:)) + mean(data_chL(sel&right,:,:))  ) ./ 2;
ipsi = squeeze(   mean(data_chL(sel&left,:,:)) + mean(data_chR(sel&right,:,:))  ) ./ 2;
ci.data(12,:,:) = ((contra-ipsi) ./ (contra+ipsi)) * 100;
cic.data(12,:,:) = contra; cii.data(12,:,:) = ipsi; clear contra ipsi


ci.data(13,:,:) = ci.data(1,:,:) - ci.data(2,:,:); %im100-infor100
ci.data(14,:,:) = ci.data(2,:,:) - ci.data(4,:,:); %infor100-infor60
ci.data(15,:,:) = ci.data(5,:,:) - ci.data(6,:,:); %valid-invalid


%% plot cvsi

if plotResults_cvsi
cfg = [];
cfg.colorbar = 'no';
cfg.zlim = [-20 20];
cfg.parameter = 'data';
cfg.figure = 'gcf';

figure;
for ch2pl = 1:6
    cfg.channel = ci.label(ch2pl);
    subplot(2,4,ch2pl); ft_singleplotTFR(cfg, ci);
end
colormap(colormap2use);

fsel = ci.freq >= 8 & ci.freq <= 12;

figure; hold on; title('alpha - contra vs ipsi');
plot(ci.time, squeeze(mean(ci.data(1,fsel,:),2)), 'color', [0, 0, 1]);
plot(ci.time, squeeze(mean(ci.data(2,fsel,:),2)), 'color', [1, 0, 0]);
plot(ci.time, squeeze(mean(ci.data(3,fsel,:),2)), 'color', [1, 0, 1]);
plot(ci.time, squeeze(mean(ci.data(4,fsel,:),2)), 'color', [0, 1, 1]);
plot(xlim, [0,0], ':k');
xlim([ci.time(1), ci.time(end)]);
legend(ci.label([1,2,3,4]));

figure; hold on; title('alpha - contra vs ipsi');
plot(ci.time, squeeze(mean(ci.data(5,fsel,:),2)), 'color', [0, 0, 1]);
plot(ci.time, squeeze(mean(ci.data(6,fsel,:),2)), 'color', [1, 0, 0]);
plot(xlim, [0,0], ':k');
xlim([ci.time(1), ci.time(end)]);
legend(ci.label([5,6]));
end

%% save
if laplacian toadd1 = '_laplacian'; else toadd1 = ''; end
if window == 1 toadd2 = ''; elseif window == 2 toadd2 = '_500mswindow'; end
if chanclusters toadd3 = '_chanclusters'; else toadd3 = ''; end

save([param.pathnewlen, 'saved_data/timefreq', toadd1,toadd2,toadd3, '__' param.subjName], 'timefreq','ci','cic','cii','ci_allchan');

end
end
end
end
