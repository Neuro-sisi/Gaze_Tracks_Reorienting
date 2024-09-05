clc;clear
close all

% % Run fastICA together
% % parameters
sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'29';'30';'31';'32';'33';'34';};

% for pp = 1
for pp = 1:length(sub_list)
    
param = getSubjParam(pp);


% % load data
load([param.pathnewlen, '/processed_data/' 'epoched_data_eeg' '__', param.subjName], 'data');

% % quick clean before proceding?
% cfg = [];
% cfg.method = 'summary';
% cfg.channel = {'EEG'};
% cfg.keepchannels = 'yes';
% data = ft_rejectvisual(cfg, data);


% % cut out EEG data only
cfg = [];
cfg.keeptrials = 'yes';
cfg.channel = {'EEG'};
d_eeg = ft_timelockanalysis(cfg, data);


% % run ica
cfg = [];
cfg.method = 'fastica';
ica = ft_componentanalysis(cfg, d_eeg);


% % look at components
cfg           = [];
cfg.component = [1:length(ica.label)];      
cfg.layout    = 'biosemi64.lay'; %'easycapM1.lay';
cfg.comment   = 'no';
figure(pp); ft_topoplotIC(cfg, ica)
colormap('jet')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, [['cueva_eeg_topoicacompo_sub' num2str(pp)], '.jpg'])


% % save appended data before badchannel replacement
save(fullfile([param.pathnewlen, '/processed_data/' 'icad_data_eeg', '__', param.subjName]));

end



%% Manually inspect bad components
clear;clc
close all;

% % Do correlation between ICA components & EOG
sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'29';'30';'31';'32';'33';'34';};

for pp = 1:length(sub_list)
    
param = getSubjParam(pp);


% % load data
load([param.pathnewlen, '/processed_data/' 'icad_data_eeg' '__', param.subjName]);


% % correlate ica timecourses with measured eog

cfg = [];
cfg.keeptrials = 'yes';
d_ica = ft_timelockanalysis(cfg, ica);
cfg.channel = {'HEOG','VEOG'};
d_eog = ft_timelockanalysis(cfg, data);

% % check eog in all trials
% figure; 
% subplot(1,2,1); plot(d_eog.time, squeeze(d_eog.trial))
% 
% alltrl = 1:size(d_eog.trial,1);
% badtrl = find(sum(squeeze(d_eog.trial)>2000,2));
% oktrl = ~ismember(alltrl, badtrl);
% 
% % remove bad ones from eog and ica, before checking their correlation below
% cfg = [];
% cfg.trials = oktrl;
% d_ica = ft_selectdata(cfg, d_ica);
% cfg.trials = oktrl;
% d_eog = ft_selectdata(cfg, d_eog);
% 
% subplot(1,2,2); plot(d_eog.time, squeeze(d_eog.trial))

% now do the correlation
x1 = [];
x2 = [];
y = [];
correlations1 = [];
correlations2 = [];
x1 = d_eog.trial(:,1,:); % h-eog
x2 = d_eog.trial(:,2,:); % v-eog
for c = 1:size(d_ica.trial,2)
y = d_ica.trial(:,c,:); % components
correlations1(c) = corr(y(:), x1(:));
correlations2(c) = corr(y(:), x2(:));
end

figure; 
subplot(2,1,1); bar(1:c, abs(correlations1),'r'); title('HEOG correlations with component timecourses');   
subplot(2,1,2); bar(1:c, abs(correlations2),'r'); title('VEOG correlations with component timecourses');   
xlabel('comp #');

drawnow;

% Open the ICA topomap for each pp to check the distribution
% % check which components are bad:
comp2rem = input('bad components are: ');


% % check

% plot before
ch2pl = ismember(data.label, 'AFz');
figure;
perm = randperm(length(data.trial));
for x = 1:16
    subplot(4,4,x); hold on;plot(data.time{perm(x)}, data.trial{perm(x)}(ch2pl,:), 'k'); xlim([data.time{1}(1), data.time{1}(end)]);
end

% remove
cfg           = [];
cfg.component = comp2rem;
data = ft_rejectcomponent(cfg, ica,data);

% plot after
for x = 1:16
    subplot(4,4,x); plot(data.time{perm(x)}, data.trial{perm(x)}(ch2pl,:), 'm');  xlim([data.time{1}(1), data.time{1}(end)]);
end
title(['pp' num2str(pp)]);

checknot = input('are you sure[1-yes, 0-no]: ');


% % save
% remove unnecessary field to save disk space
ica = rmfield(ica, 'trial');

% now actually saves
save([param.pathnewlen, '/saved_data/' 'ICAcomponents', '__' param.subjName], 'comp2rem','ica');


% close plots
close all

end