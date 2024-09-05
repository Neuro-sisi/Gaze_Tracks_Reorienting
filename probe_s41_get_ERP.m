clc;clear ;close all
ft_defaults

%% parameters
sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'29';'30';'31';'32';'33';'34';};

for laplacian    = 0; % defalt--0
for chancluster  = 0; % defalt--0; 1-averaged over cluster channels; 0-single channel PO7/8
% for timebl     = [1 2 3]; % defalt--1; 1-precue baseline: -0.2-0; 2-prefixation baseline: -1.2--1; 3-preprobe baseline: 1.3--1.5
for timebl     = 3; % defalt--1; 1-precue baseline: -0.2-0; 2-prefixation baseline: -1.2--1; 3-preprobe baseline: 1.3--1.5

for pp           = 1:length(sub_list)
% for pp = 1

plotResults_topo = 0;
plotResults_sc   = 0;
plotResults_ci   = 0;

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

%% baseline correct
% figure; hold on; plot(data.time{1}, data.trial{1}(1,:), 'k');
cfg                 = [];
cfg.demean          = 'yes';
if timebl==1
cfg.baselinewindow  = [-.2 0]; % 200 ms pre-cue baseline
elseif timebl==2
cfg.baselinewindow  = [-1.3 -1]; % 200 ms pre-fixation baseline
elseif timebl==3
cfg.baselinewindow  = [1.3 1.5]; % 200 ms pre-probe baseline
end    
data                = ft_preprocessing(cfg, data);

% plot(data.time{1}, data.trial{1}(1,:), 'm');

% lpfilter
cfg          = [];
cfg.lpfilter = 'yes';
cfg.lpfreq   = 30;
data         = ft_preprocessing(cfg, data);
% plot(data.time{1}, data.trial{1}(1,:), 'b');

%% get ERP
cfg             = [];
cfg.keeptrials  = 'yes';
tl              = ft_timelockanalysis(cfg, data); % data in tl.trial matrix

tl.time         = tl.time*1000; % s to ms

%% split into three conditions
imper100 = ismember(tl.trialinfo(:,2), 4);
high100  = ismember(tl.trialinfo(:,2), 1);
med80    = ismember(tl.trialinfo(:,2), 2);
low60    = ismember(tl.trialinfo(:,2), 3);

left     = ismember(tl.trialinfo(:,1), [41]);
right    = ismember(tl.trialinfo(:,1), [42]);

valid = ismember(tl.trialinfo(:,3), 1);
invalid  = ismember(tl.trialinfo(:,3), 2);

%% ERP per condition, put into a structure we can later plot
erp = [];
erp.time = tl.time;
erp.label = tl.label;
erp.dimord = 'chan_time';

erp.imper100 = squeeze(mean(tl.trial(imper100,:,:)));
erp.high100 = squeeze(mean(tl.trial(high100,:,:)));
erp.med80 = squeeze(mean(tl.trial(med80,:,:)));
erp.low60 = squeeze(mean(tl.trial(low60,:,:)));

erp.valid = squeeze(mean(tl.trial(valid,:,:)));
erp.invalid = squeeze(mean(tl.trial(invalid,:,:)));

erp.validim100 = squeeze(mean(tl.trial(valid&imper100,:,:)));
erp.validinfor100 = squeeze(mean(tl.trial(valid&high100,:,:)));
erp.validinfor80 = squeeze(mean(tl.trial(valid&med80,:,:)));
erp.validinfor60 = squeeze(mean(tl.trial(valid&low60,:,:)));
erp.validinfor40 = squeeze(mean(tl.trial(invalid&low60,:,:)));
erp.validinfor20 = squeeze(mean(tl.trial(invalid&med80,:,:)));

erp.imper_vs_high = erp.imper100 - erp.high100;
erp.high_vs_low   = erp.high100  - erp.low60;
erp.valid_vs_invalid   = erp.valid  - erp.invalid;


%% running topos
if plotResults_topo

cfg = [];
cfg.layout = 'easycapM1.mat';
cfg.comment = 'no';
cfg.style = 'straight';
cfg.marker = 'off';
cfg.zlim = [-5 5]; if laplacian cfg.zlim = [-1 1]*1e-3; end

    stepsize = 250; times1 = [-250:stepsize:2500]; n = length(times1);
    figure; 
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(4,n-1,t);           cfg.parameter = 'imper100'; ft_topoplotER(cfg, erp); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
        subplot(4,n-1,t+n-1);       cfg.parameter = 'high100';  ft_topoplotER(cfg, erp);
        subplot(4,n-1,t+(n-1)*2);   cfg.parameter = 'med80';    ft_topoplotER(cfg, erp);
        subplot(4,n-1,t+(n-1)*3);   cfg.parameter = 'low60';    ft_topoplotER(cfg, erp);
    end
    colormap(colormap2use);

    figure;
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(2,n-1,t);           cfg.parameter = 'imper_vs_high'; ft_topoplotER(cfg, erp); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
        subplot(2,n-1,t+n-1);       cfg.parameter = 'high_vs_low';   ft_topoplotER(cfg, erp);
    end
    colormap(colormap2use);
    drawnow;
end

%% ERP timecourses in selected channels

if plotResults_sc
figure;
for ch = 1:3
    if ch == 1         choi = ismember(erp.label, 'Oz'); title2use = 'Oz';
    elseif ch == 2     choi = ismember(erp.label, 'Pz'); title2use = 'Pz';
    elseif ch == 3     choi = ismember(erp.label, 'Fz'); title2use = 'Fz';
    end
    
    subplot(2,3,ch); hold on; title(title2use);
    plot(erp.time, erp.imper100(choi,:), 'b');
    plot(erp.time, erp.high100(choi,:), 'r');
    plot(erp.time, erp.med80(choi,:), 'm');
    plot(erp.time, erp.low60(choi,:), 'c');
    legend('imper100','high100','med80','med60');
    xlim([erp.time(1), erp.time(end)]);
        
    subplot(2,3,ch+3); hold on; title(title2use);
    plot(erp.time, erp.imper_vs_high(choi,:), 'g');
    plot(erp.time, erp.high_vs_low(choi,:), 'y');
    plot(xlim, [0,0], '--k');
    legend('imper-vs-high','high-vs-low');
    xlim([erp.time(1), erp.time(end)]);    
end
drawnow;
end


%% get contra & ipsi & cvsi erps
chL = []; chR = [];
chL = ismember(erp.label, {'Fp1','AF3','AF7','F1','F3','F5','F7','FC1','FC3','FC5','FT7','C1','C3','C5','T7','CP1','CP3','CP5','TP7','P1','P3','P5','P7','P9','PO3','PO7','O1'});
chR = ismember(erp.label, {'Fp2','AF4','AF8','F2','F4','F6','F8','FC2','FC4','FC6','FT8','C2','C4','C6','T8','CP2','CP4','CP6','TP8','P2','P4','P6','P8','P10','PO4','PO8','O2'});

ci_allchan = [];
ci_allchan.time = erp.time;
ci_allchan.label = {'imper100-contra','imper100-ipsi','imper100-cvsi',...
            'high100-contra','high100-ipsi','high100-cvsi',...
            'med80-contra','med80-ipsi','med80-cvsi',...
            'low60-contra','low60-ipsi','low60-cvsi',...
            'valid-contra','valid-ipsi','valid-cvsi',...
            'invalid-contra','invalid-ipsi','invalid-cvsi',...         
            'validim100-contra','validim100-ipsi','validim100-cvsi',...
            'validinfor100-contra','validinfor100-ipsi','validinfor100-cvsi',...
            'validinfor80-contra','validinfor80-ipsi','validinfor80-cvsi',...
            'validinfor60-contra','validinfor60-ipsi','validinfor60-cvsi',...
            'validinfor40-contra','validinfor40-ipsi','validinfor40-cvsi',...
            'validinfor20-contra','validinfor20-ipsi','validinfor20-cvsi',...
            };

ci_allchan.data(1,:,:) =  (squeeze(mean(tl.trial(imper100&left,chR,:))) + squeeze(mean(tl.trial(imper100&right,chL,:)))    ) / 2;
ci_allchan.data(2,:,:) =  (squeeze(mean(tl.trial(imper100&left,chL,:))) + squeeze(mean(tl.trial(imper100&right,chR,:)))    ) / 2;
ci_allchan.data(3,:,:) =  ci_allchan.data(1,:,:) - ci_allchan.data(2,:,:); 
ci_allchan.data(4,:,:) =  (squeeze(mean(tl.trial(high100&left,chR,:))) + squeeze(mean(tl.trial(high100&right,chL,:)))    ) / 2;
ci_allchan.data(5,:,:) =  (squeeze(mean(tl.trial(high100&left,chL,:))) + squeeze(mean(tl.trial(high100&right,chR,:)))    ) / 2;
ci_allchan.data(6,:,:) =  ci_allchan.data(4,:,:) - ci_allchan.data(5,:,:); 
ci_allchan.data(7,:,:) =  (squeeze(mean(tl.trial(med80&left,chR,:))) + squeeze(mean(tl.trial(med80&right,chL,:)))    ) / 2;
ci_allchan.data(8,:,:) =  (squeeze(mean(tl.trial(med80&left,chL,:))) + squeeze(mean(tl.trial(med80&right,chR,:)))    ) / 2;
ci_allchan.data(9,:,:) =  ci_allchan.data(7,:,:) - ci_allchan.data(8,:,:); 
ci_allchan.data(10,:,:) =  (squeeze(mean(tl.trial(low60&left,chR,:))) + squeeze(mean(tl.trial(low60&right,chL,:)))    ) / 2;
ci_allchan.data(11,:,:) =  (squeeze(mean(tl.trial(low60&left,chL,:))) + squeeze(mean(tl.trial(low60&right,chR,:)))    ) / 2;
ci_allchan.data(12,:,:) =  ci_allchan.data(10,:,:) - ci_allchan.data(11,:,:); 
ci_allchan.data(13,:,:) =  (squeeze(mean(tl.trial(valid&left,chR,:))) + squeeze(mean(tl.trial(valid&right,chL,:)))    ) / 2;
ci_allchan.data(14,:,:) =  (squeeze(mean(tl.trial(valid&left,chL,:))) + squeeze(mean(tl.trial(valid&right,chR,:)))    ) / 2;
ci_allchan.data(15,:,:) =  ci_allchan.data(13,:,:) - ci_allchan.data(14,:,:); 
ci_allchan.data(16,:,:) =  (squeeze(mean(tl.trial(invalid&left,chR,:))) + squeeze(mean(tl.trial(invalid&right,chL,:)))    ) / 2;
ci_allchan.data(17,:,:) =  (squeeze(mean(tl.trial(invalid&left,chL,:))) + squeeze(mean(tl.trial(invalid&right,chR,:)))    ) / 2;
ci_allchan.data(18,:,:) =  ci_allchan.data(16,:,:) - ci_allchan.data(17,:,:); 
ci_allchan.data(19,:,:) =  (squeeze(mean(tl.trial(valid&imper100&left,chR,:))) + squeeze(mean(tl.trial(valid&imper100&right,chL,:)))    ) / 2;
ci_allchan.data(20,:,:) =  (squeeze(mean(tl.trial(valid&imper100&left,chL,:))) + squeeze(mean(tl.trial(valid&imper100&right,chR,:)))    ) / 2;
ci_allchan.data(21,:,:) =  ci_allchan.data(19,:,:) - ci_allchan.data(20,:,:); 
ci_allchan.data(22,:,:) =  (squeeze(mean(tl.trial(valid&high100&left,chR,:))) + squeeze(mean(tl.trial(valid&high100&right,chL,:)))    ) / 2;
ci_allchan.data(23,:,:) =  (squeeze(mean(tl.trial(valid&high100&left,chL,:))) + squeeze(mean(tl.trial(valid&high100&right,chR,:)))    ) / 2;
ci_allchan.data(24,:,:) =  ci_allchan.data(22,:,:) - ci_allchan.data(23,:,:); 
ci_allchan.data(25,:,:) =  (squeeze(mean(tl.trial(valid&med80&left,chR,:))) + squeeze(mean(tl.trial(valid&med80&right,chL,:)))    ) / 2;
ci_allchan.data(26,:,:) =  (squeeze(mean(tl.trial(valid&med80&left,chL,:))) + squeeze(mean(tl.trial(valid&med80&right,chR,:)))    ) / 2;
ci_allchan.data(27,:,:) =  ci_allchan.data(25,:,:) - ci_allchan.data(26,:,:); 
ci_allchan.data(28,:,:) =  (squeeze(mean(tl.trial(valid&low60&left,chR,:))) + squeeze(mean(tl.trial(valid&low60&right,chL,:)))    ) / 2;
ci_allchan.data(29,:,:) =  (squeeze(mean(tl.trial(valid&low60&left,chL,:))) + squeeze(mean(tl.trial(valid&low60&right,chR,:)))    ) / 2;
ci_allchan.data(30,:,:) =  ci_allchan.data(28,:,:) - ci_allchan.data(29,:,:); 
ci_allchan.data(31,:,:) =  (squeeze(mean(tl.trial(invalid&low60&left,chR,:))) + squeeze(mean(tl.trial(invalid&low60&right,chL,:)))    ) / 2;
ci_allchan.data(32,:,:) =  (squeeze(mean(tl.trial(invalid&low60&left,chL,:))) + squeeze(mean(tl.trial(invalid&low60&right,chR,:)))    ) / 2;
ci_allchan.data(33,:,:) =  ci_allchan.data(31,:,:) - ci_allchan.data(32,:,:); 
ci_allchan.data(34,:,:) =  (squeeze(mean(tl.trial(invalid&med80&left,chR,:))) + squeeze(mean(tl.trial(invalid&med80&right,chL,:)))    ) / 2;
ci_allchan.data(35,:,:) =  (squeeze(mean(tl.trial(invalid&med80&left,chL,:))) + squeeze(mean(tl.trial(invalid&med80&right,chR,:)))    ) / 2;
ci_allchan.data(36,:,:) =  ci_allchan.data(34,:,:) - ci_allchan.data(35,:,:); 


%% get averaged contra and ipsi--channel cluster/single channel
if chancluster
chL = []; chR = [];
chL = ismember(erp.label, {'P1','P3','P5','P7','P9','PO3','PO7','O1'});
chR = ismember(erp.label, {'P2','P4','P6','P8','P10','PO4','PO8','O2'});

ci = [];
ci.time = erp.time;
ci.label = {'imper100-contra','imper100-ipsi','imper100-cvsi',...
            'high100-contra','high100-ipsi','high100-cvsi',...
            'med80-contra','med80-ipsi','med80-cvsi',...
            'low60-contra','low60-ipsi','low60-cvsi',...
            'valid-contra','valid-ipsi','valid-cvsi',...
            'invalid-contra','invalid-ipsi','invalid-cvsi',...         
            'validim100-contra','validim100-ipsi','validim100-cvsi',...
            'validinfor100-contra','validinfor100-ipsi','validinfor100-cvsi',...
            'validinfor80-contra','validinfor80-ipsi','validinfor80-cvsi',...
            'validinfor60-contra','validinfor60-ipsi','validinfor60-cvsi',...
            'validinfor40-contra','validinfor40-ipsi','validinfor40-cvsi',...
            'validinfor20-contra','validinfor20-ipsi','validinfor20-cvsi',...
            };

ci.data(1,:) =  squeeze(    mean(squeeze(mean(tl.trial(imper100&left,chR,:)))) + mean(squeeze(mean(tl.trial(imper100&right,chL,:))))    ) / 2;
ci.data(2,:) =  squeeze(    mean(squeeze(mean(tl.trial(imper100&left,chL,:)))) + mean(squeeze(mean(tl.trial(imper100&right,chR,:))))    ) / 2;
ci.data(3,:) =  ci.data(1,:) - ci.data(2,:); 
ci.data(4,:) =  squeeze(    mean(squeeze(mean(tl.trial(high100&left,chR,:)))) + mean(squeeze(mean(tl.trial(high100&right,chL,:))))    ) / 2;
ci.data(5,:) =  squeeze(    mean(squeeze(mean(tl.trial(high100&left,chL,:)))) + mean(squeeze(mean(tl.trial(high100&right,chR,:))))    ) / 2;
ci.data(6,:) =  ci.data(4,:) - ci.data(5,:); 
ci.data(7,:) =  squeeze(    mean(squeeze(mean(tl.trial(med80&left,chR,:)))) + mean(squeeze(mean(tl.trial(med80&right,chL,:))))    ) / 2;
ci.data(8,:) =  squeeze(    mean(squeeze(mean(tl.trial(med80&left,chL,:)))) + mean(squeeze(mean(tl.trial(med80&right,chR,:))))    ) / 2;
ci.data(9,:) =  ci.data(7,:) - ci.data(8,:); 
ci.data(10,:) =  squeeze(    mean(squeeze(mean(tl.trial(low60&left,chR,:)))) + mean(squeeze(mean(tl.trial(low60&right,chL,:))))    ) / 2;
ci.data(11,:) =  squeeze(    mean(squeeze(mean(tl.trial(low60&left,chL,:)))) + mean(squeeze(mean(tl.trial(low60&right,chR,:))))    ) / 2;
ci.data(12,:) =  ci.data(10,:) - ci.data(11,:); 
ci.data(13,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&left,chR,:)))) + mean(squeeze(mean(tl.trial(valid&right,chL,:))))    ) / 2;
ci.data(14,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&left,chL,:)))) + mean(squeeze(mean(tl.trial(valid&right,chR,:))))    ) / 2;
ci.data(15,:) =  ci.data(13,:) - ci.data(14,:); 
ci.data(16,:) =  squeeze(    mean(squeeze(mean(tl.trial(invalid&left,chR,:)))) + mean(squeeze(mean(tl.trial(invalid&right,chL,:))))    ) / 2;
ci.data(17,:) =  squeeze(    mean(squeeze(mean(tl.trial(invalid&left,chL,:)))) + mean(squeeze(mean(tl.trial(invalid&right,chR,:))))    ) / 2;
ci.data(18,:) =  ci.data(16,:) - ci.data(17,:); 
ci.data(19,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&imper100&left,chR,:)))) + mean(squeeze(mean(tl.trial(valid&imper100&right,chL,:))))    ) / 2;
ci.data(20,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&imper100&left,chL,:)))) + mean(squeeze(mean(tl.trial(valid&imper100&right,chR,:))))    ) / 2;
ci.data(21,:) =  ci.data(19,:) - ci.data(20,:); 
ci.data(22,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&high100&left,chR,:)))) + mean(squeeze(mean(tl.trial(valid&high100&right,chL,:))))    ) / 2;
ci.data(23,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&high100&left,chL,:)))) + mean(squeeze(mean(tl.trial(valid&high100&right,chR,:))))    ) / 2;
ci.data(24,:) =  ci.data(22,:) - ci.data(23,:); 
ci.data(25,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&med80&left,chR,:)))) + mean(squeeze(mean(tl.trial(valid&med80&right,chL,:))))    ) / 2;
ci.data(26,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&med80&left,chL,:)))) + mean(squeeze(mean(tl.trial(valid&med80&right,chR,:))))    ) / 2;
ci.data(27,:) =  ci.data(25,:) - ci.data(26,:); 
ci.data(28,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&low60&left,chR,:)))) + mean(squeeze(mean(tl.trial(valid&low60&right,chL,:))))    ) / 2;
ci.data(29,:) =  squeeze(    mean(squeeze(mean(tl.trial(valid&low60&left,chL,:)))) + mean(squeeze(mean(tl.trial(valid&low60&right,chR,:))))    ) / 2;
ci.data(30,:) =  ci.data(28,:) - ci.data(29,:); 
ci.data(31,:) =  squeeze(    mean(squeeze(mean(tl.trial(invalid&low60&left,chR,:)))) + mean(squeeze(mean(tl.trial(invalid&low60&right,chL,:))))    ) / 2;
ci.data(32,:) =  squeeze(    mean(squeeze(mean(tl.trial(invalid&low60&left,chL,:)))) + mean(squeeze(mean(tl.trial(invalid&low60&right,chR,:))))    ) / 2;
ci.data(33,:) =  ci.data(31,:) - ci.data(32,:); 
ci.data(34,:) =  squeeze(    mean(squeeze(mean(tl.trial(invalid&med80&left,chR,:)))) + mean(squeeze(mean(tl.trial(invalid&med80&right,chL,:))))    ) / 2;
ci.data(35,:) =  squeeze(    mean(squeeze(mean(tl.trial(invalid&med80&left,chL,:)))) + mean(squeeze(mean(tl.trial(invalid&med80&right,chR,:))))    ) / 2;
ci.data(36,:) =  ci.data(34,:) - ci.data(35,:); 


if plotResults_ci
% figure;
% subplot(2,4,1); plot(ci.time, ci.data(1:2,:));   legend(ci.label(1:2)); title('imper100'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,2); plot(ci.time, ci.data(4:5,:));   legend(ci.label(4:5)); title('high100'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,3); plot(ci.time, ci.data(7:8,:));   legend(ci.label(7:8)); title('med80'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,4); plot(ci.time, ci.data(10:11,:)); legend(ci.label(10:11)); title('low60'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,5); plot(ci.time, ci.data(3,:));     legend(ci.label(3)); hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,6); plot(ci.time, ci.data(6,:));     legend(ci.label(6)); hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,7); plot(ci.time, ci.data(9,:));     legend(ci.label(9));  hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,8); plot(ci.time, ci.data(12,:));    legend(ci.label(9));  hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);

figure;
plot(ci.time, ci.data([3:3:36],:)); legend(ci.label([3:3:36])); title('cvsi');
hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
end

drawnow;

else

% % get contra and ipsi--single channel
chL = []; chR = [];
chL = ismember(erp.label, 'PO7');
chR = ismember(erp.label, 'PO8');

ci = [];
ci.time = erp.time;
ci.label = {'imper100-contra','imper100-ipsi','imper100-cvsi',...
            'high100-contra','high100-ipsi','high100-cvsi',...
            'med80-contra','med80-ipsi','med80-cvsi',...
            'low60-contra','low60-ipsi','low60-cvsi',...
            'valid-contra','valid-ipsi','valid-cvsi',...
            'invalid-contra','invalid-ipsi','invalid-cvsi',...         
            'validim100-contra','validim100-ipsi','validim100-cvsi',...
            'validinfor100-contra','validinfor100-ipsi','validinfor100-cvsi',...
            'validinfor80-contra','validinfor80-ipsi','validinfor80-cvsi',...
            'validinfor60-contra','validinfor60-ipsi','validinfor60-cvsi',...
            'validinfor40-contra','validinfor40-ipsi','validinfor40-cvsi',...
            'validinfor20-contra','validinfor20-ipsi','validinfor20-cvsi',...
            };

ci.data(1,:) =  squeeze(    mean(tl.trial(imper100&left,chR,:)) + mean(tl.trial(imper100&right,chL,:))    ) / 2;
ci.data(2,:) =  squeeze(    mean(tl.trial(imper100&left,chL,:)) + mean(tl.trial(imper100&right,chR,:))    ) / 2;
ci.data(3,:) =  ci.data(1,:) - ci.data(2,:); 
ci.data(4,:) =  squeeze(    mean(tl.trial(high100&left,chR,:)) + mean(tl.trial(high100&right,chL,:))    ) / 2;
ci.data(5,:) =  squeeze(    mean(tl.trial(high100&left,chL,:)) + mean(tl.trial(high100&right,chR,:))    ) / 2;
ci.data(6,:) =  ci.data(4,:) - ci.data(5,:); 
ci.data(7,:) =  squeeze(    mean(tl.trial(med80&left,chR,:)) + mean(tl.trial(med80&right,chL,:))    ) / 2;
ci.data(8,:) =  squeeze(    mean(tl.trial(med80&left,chL,:)) + mean(tl.trial(med80&right,chR,:))    ) / 2;
ci.data(9,:) =  ci.data(7,:) - ci.data(8,:); 
ci.data(10,:) =  squeeze(    mean(tl.trial(low60&left,chR,:)) + mean(tl.trial(low60&right,chL,:))    ) / 2;
ci.data(11,:) =  squeeze(    mean(tl.trial(low60&left,chL,:)) + mean(tl.trial(low60&right,chR,:))    ) / 2;
ci.data(12,:) =  ci.data(10,:) - ci.data(11,:); 
ci.data(13,:) =  squeeze(    mean(tl.trial(valid&left,chR,:)) + mean(tl.trial(valid&right,chL,:))    ) / 2;
ci.data(14,:) =  squeeze(    mean(tl.trial(valid&left,chL,:)) + mean(tl.trial(valid&right,chR,:))    ) / 2;
ci.data(15,:) =  ci.data(13,:) - ci.data(14,:); 
ci.data(16,:) =  squeeze(    mean(tl.trial(invalid&left,chR,:)) + mean(tl.trial(invalid&right,chL,:))    ) / 2;
ci.data(17,:) =  squeeze(    mean(tl.trial(invalid&left,chL,:)) + mean(tl.trial(invalid&right,chR,:))    ) / 2;
ci.data(18,:) =  ci.data(16,:) - ci.data(17,:); 
ci.data(19,:) =  squeeze(    mean(tl.trial(valid&imper100&left,chR,:)) + mean(tl.trial(valid&imper100&right,chL,:))    ) / 2;
ci.data(20,:) =  squeeze(    mean(tl.trial(valid&imper100&left,chL,:)) + mean(tl.trial(valid&imper100&right,chR,:))    ) / 2;
ci.data(21,:) =  ci.data(19,:) - ci.data(20,:); 
ci.data(22,:) =  squeeze(    mean(tl.trial(valid&high100&left,chR,:)) + mean(tl.trial(valid&high100&right,chL,:))    ) / 2;
ci.data(23,:) =  squeeze(    mean(tl.trial(valid&high100&left,chL,:)) + mean(tl.trial(valid&high100&right,chR,:))    ) / 2;
ci.data(24,:) =  ci.data(22,:) - ci.data(23,:); 
ci.data(25,:) =  squeeze(    mean(tl.trial(valid&med80&left,chR,:)) + mean(tl.trial(valid&med80&right,chL,:))    ) / 2;
ci.data(26,:) =  squeeze(    mean(tl.trial(valid&med80&left,chL,:)) + mean(tl.trial(valid&med80&right,chR,:))    ) / 2;
ci.data(27,:) =  ci.data(25,:) - ci.data(26,:); 
ci.data(28,:) =  squeeze(    mean(tl.trial(valid&low60&left,chR,:)) + mean(tl.trial(valid&low60&right,chL,:))    ) / 2;
ci.data(29,:) =  squeeze(    mean(tl.trial(valid&low60&left,chL,:)) + mean(tl.trial(valid&low60&right,chR,:))    ) / 2;
ci.data(30,:) =  ci.data(28,:) - ci.data(29,:); 
ci.data(31,:) =  squeeze(    mean(tl.trial(invalid&low60&left,chR,:)) + mean(tl.trial(invalid&low60&right,chL,:))    ) / 2;
ci.data(32,:) =  squeeze(    mean(tl.trial(invalid&low60&left,chL,:)) + mean(tl.trial(invalid&low60&right,chR,:))    ) / 2;
ci.data(33,:) =  ci.data(31,:) - ci.data(32,:); 
ci.data(34,:) =  squeeze(    mean(tl.trial(invalid&med80&left,chR,:)) + mean(tl.trial(invalid&med80&right,chL,:))    ) / 2;
ci.data(35,:) =  squeeze(    mean(tl.trial(invalid&med80&left,chL,:)) + mean(tl.trial(invalid&med80&right,chR,:))    ) / 2;
ci.data(36,:) =  ci.data(34,:) - ci.data(35,:); 

if plotResults_ci
% figure;
% subplot(2,4,1); plot(ci.time, ci.data(1:2,:));   legend(ci.label(1:2)); title('imper100'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,2); plot(ci.time, ci.data(4:5,:));   legend(ci.label(4:5)); title('high100'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,3); plot(ci.time, ci.data(7:8,:));   legend(ci.label(7:8)); title('med80'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,4); plot(ci.time, ci.data(10:11,:)); legend(ci.label(10:11)); title('low60'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,5); plot(ci.time, ci.data(3,:));     legend(ci.label(3)); hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,6); plot(ci.time, ci.data(6,:));     legend(ci.label(6)); hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,7); plot(ci.time, ci.data(9,:));     legend(ci.label(9));  hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
% subplot(2,4,8); plot(ci.time, ci.data(12,:));    legend(ci.label(9));  hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);

figure;
plot(ci.time, ci.data([3:3:36],:)); legend(ci.label([3:3:36])); title('cvsi');
hold on; plot(xlim, [0,0], '--k'); xlim([ci.time(1), ci.time(end)]);
end

drawnow;

end


%% save
if laplacian toadd1 = '_laplacian'; else toadd1 = ''; end
if chancluster toadd2 = ''; else toadd2 = '_sinch'; end
if timebl==1 toadd3 = '_precuebl'; elseif timebl==2 toadd3 = '_prefixbl'; elseif timebl==3 toadd3 = '_preprobebl'; end
save([param.pathnewlen, 'saved_data/erp', toadd1,toadd2,toadd3,'__' param.subjName], 'erp','ci','ci_allchan');

clear erp ci ci_allchan
close all
end
end
end
end