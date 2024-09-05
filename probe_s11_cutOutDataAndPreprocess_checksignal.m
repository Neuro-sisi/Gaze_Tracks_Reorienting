clc;clear
close all

%% Check bad channels

%% Probe-locked analysis

%% add fieldtrip & path
% addpath '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/Data_ana/EEGdata_Ana_FvE/toolbox/fieldtrip-20201023'
ft_defaults

addpath '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/Data_ana/EEGdata_Ana_FvE/';

%% S1: append data & check bad channelss

% % set loops
sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'29';'30';'31';'32';'33';'34';};
% sub'27' not included due to lacking data


% % parameters
for pp = 1:length(sub_list)
% for pp = 1

% cut out relative to retrocue
values2use  = [41,42]; % retrocue
prestim     = -1.5;  % from X before retrocue
poststim    = 2.5; % until X sec after 1s after probe


% % correct for trigger mess-up in pp1
if pp==1 values2use = [41,43]; end % in pp 1 all triggers +1 -- will correct for it after epoching, by again subtracting 1 from trialinfo.


% % get relevant file directories and subject info
param = getSubjParam(pp);

eegfile1 = [param.eeg1];
eegfile2 = [param.eeg2];
eegfile3 = [param.eeg3];
eegfile4 = [param.eeg4];


% % check datafile
checkch = 0;

if checkch
% % check header
hdr1 = ft_read_header(eegfile1);
hdr2 = ft_read_header(eegfile2);
hdr3 = ft_read_header(eegfile3);
hdr4 = ft_read_header(eegfile4);

% check events
event1 = ft_read_event(eegfile1);
sel1 = find(strcmp({event1.type}, 'STATUS'));
event1 = event1(sel1);

event2 = ft_read_event(eegfile2);
sel2 = find(strcmp({event2.type}, 'STATUS'));
event2 = event2(sel2);

event3 = ft_read_event(eegfile3);
sel3 = find(strcmp({event3.type}, 'STATUS'));
event3 = event3(sel3);

event4 = ft_read_event(eegfile4);
sel4 = find(strcmp({event4.type}, 'STATUS'));
event4 = event4(sel4);

figure(98); 
subplot(1,4,1);
plot([event1.sample], [event1.value], '.');
xlabel('sample in recording'); ylabel('trigger value');

subplot(1,4,2);
plot([event2.sample], [event2.value], '.');
xlabel('sample in recording'); ylabel('trigger value');

subplot(1,4,3);
plot([event3.sample], [event3.value], '.');
xlabel('sample in recording'); ylabel('trigger value');

subplot(1,4,4);
plot([event4.sample], [event4.value], '.');
xlabel('sample in recording'); ylabel('trigger value');
end

% % epoch

% ds 1
cfg = [];
cfg.dataset = eegfile1;
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = values2use;
cfg.trialdef.prestim    = -prestim;
cfg.trialdef.poststim   = poststim;
cfg = ft_definetrial(cfg);
cfg.reref = 'yes';
cfg.refchannel = {'EXG1','EXG2'}; % reref to avg of both mastoids
cfg.demean = 'yes';
data1 = ft_preprocessing(cfg);

% ds 2
cfg = [];
cfg.dataset = eegfile2;
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = values2use;
cfg.trialdef.prestim    = -prestim;
cfg.trialdef.poststim   = poststim;
cfg = ft_definetrial(cfg);
cfg.reref = 'yes';
cfg.refchannel = {'EXG1','EXG2'}; % reref to avg of both mastoids
cfg.demean = 'yes';
data2 = ft_preprocessing(cfg);

% ds 3
cfg = [];
cfg.dataset = eegfile3;
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = values2use;
cfg.trialdef.prestim    = -prestim;
cfg.trialdef.poststim   = poststim;
cfg = ft_definetrial(cfg);
cfg.reref = 'yes';
cfg.refchannel = {'EXG1','EXG2'}; % reref to avg of both mastoids
cfg.demean = 'yes';
data3 = ft_preprocessing(cfg);

% ds 4
cfg = [];
cfg.dataset = eegfile4;
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = values2use;
cfg.trialdef.prestim    = -prestim;
cfg.trialdef.poststim   = poststim;
cfg = ft_definetrial(cfg);
cfg.reref = 'yes';
cfg.refchannel = {'EXG1','EXG2'}; % reref to avg of both mastoids
cfg.demean = 'yes';
data4 = ft_preprocessing(cfg);


% % append
cfg = [];
data = ft_appenddata(cfg, data1, data2, data3, data4);
clear data1 data2 data3 data4


% % correction for trigger mess-up in pp1
if pp==1 data.trialinfo(data.trialinfo(:,1)==43,1) = 42; end % now correct trigger 43 back to 42.


% % create bipolar EMGs and EOG
eog1 = ismember(data.label,  'EXG3');    eog2 = ismember(data.label,  'EXG4'); % bipolar EOG chan 
eog3 = ismember(data.label,  'EXG5');    eog4 = ismember(data.label,  'EXG6'); % bipolar EOG chan 2  

data.label(end+1:end+2) = {'HEOG','VEOG'};
for trl = 1:size(data.trial,2)
    data.trial{trl}(end+1,:) = data.trial{trl}(eog1,:) - data.trial{trl}(eog2,:); %HEOG
    data.trial{trl}(end+1,:) = data.trial{trl}(eog3,:) - data.trial{trl}(eog4,:); %VEOG
end


% % convert to 10-10 system labels
cfg = [];
cfg.layout = 'biosemi64.lay'; % assumes this is standard
layout = ft_prepare_layout(cfg);
data.label(1:64) = layout.label(1:64);


% % keep only channels of interest
cfg = [];
cfg.channel = {'EEG','HEOG','VEOG'};
data = ft_selectdata(cfg, data);


% % plot some data to check

if checkch
figure(99); 
subplot(2,2,1); plot(data.time{1}, data.trial{1}(1:64,:)); title('trial 1');
subplot(2,2,2);  plot(data.time{100}, data.trial{100}(1:64,:)); title('trial 100'); 
legend(data.label(1:64));
subplot(2,2,3); plot(data.time{1}, data.trial{1}(65:end,:));
subplot(2,2,4);  plot(data.time{100}, data.trial{100}(65:end,:)); 
legend(data.label(65:end));

figure(pp); 
for ch = 1:64
subplot(8,8,ch); hold on;
plot(data.time{1}, data.trial{1}(ch,:), 'b');
plot(data.time{100}, data.trial{100}(ch,:), 'r');
title(data.label(ch));
xlim([-1 1]); 
ylim([-100 100]);
end
drawnow;
saveas(gcf, [['cueva_eeg_examtrial_sub' num2str(pp)], '.jpg'])
end

% % save appended data before badchannel replacement
% save([param.pathnewlen, '/processed_data/' 'appended_data_eeg', '__', param.subjName], 'data');
% clear data
% close all
% end pp loop
% end



%% replace bad channel & resample
% sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'29';'30';'31';'32';'33';'34';};

% % parameters
% for pp = 1:length(sub_list)
%     
% param = getSubjParam(pp);
% 
% 
% % % load data
% load([param.pathnewlen, '/processed_data/' 'appended_data_eeg' '__', param.subjName], 'data');


% % interpolate noisy channels

% noisy channel P2
if ismember(pp, [3 5 8 9 11 12 13 16 17 18 20]) % noisy channel B25/P2
    badchan = ismember(data.label, 'P2');
    replacechans =  ismember(data.label, {'CP2', 'P4'});
    for trl = 1:size(data.trial,2)
        data.trial{trl}(badchan,:) = mean(data.trial{trl}(replacechans,:));
    end      
end

% noisy channel PO4
if ismember(pp, [2,5,6,8,16,17,18,19,20,23]) % noisy channel PO4
    badchan = ismember(data.label, 'PO4');
    if pp==17
    replacechans =  ismember(data.label, {'POz', 'O2'});
    elseif pp==18
    replacechans =  ismember(data.label, {'PO8', 'P4'});
    else
    replacechans =  ismember(data.label, {'PO8', 'POz'});
    end
    for trl = 1:size(data.trial,2)
        data.trial{trl}(badchan,:) = mean(data.trial{trl}(replacechans,:));
    end      
end

% noisy channel O2
if ismember(pp, [1,8,18,20]) 
    badchan = ismember(data.label, 'O2');
    replacechans =  ismember(data.label, {'PO8', 'Oz'});
    for trl = 1:size(data.trial,2)
        data.trial{trl}(badchan,:) = mean(data.trial{trl}(replacechans,:));
    end      
end

% noisy channel Pz
if ismember(pp, [16,22]) % noisy channel Pz
    badchan = ismember(data.label, 'Pz');
    replacechans =  ismember(data.label, {'P1', 'CPz'});
    for trl = 1:size(data.trial,2)
        data.trial{trl}(badchan,:) = mean(data.trial{trl}(replacechans,:));
    end      
end

% noisy channel POz
if ismember(pp, [18,22]) % noisy channel Pz
    badchan = ismember(data.label, 'POz');
    replacechans =  ismember(data.label, {'PO3', 'Oz'});
    for trl = 1:size(data.trial,2)
        data.trial{trl}(badchan,:) = mean(data.trial{trl}(replacechans,:));
    end      
end

% noisy channel TP8
if ismember(pp, [2]) 
    badchan = ismember(data.label, 'TP8');
    replacechans =  ismember(data.label, {'CP6', 'T8'});
    for trl = 1:size(data.trial,2)
        data.trial{trl}(badchan,:) = mean(data.trial{trl}(replacechans,:));
    end      
end

% noisy channel FC5
if ismember(pp, [19]) 
    badchan = ismember(data.label, 'FC5');
    replacechans =  ismember(data.label, {'FT7', 'FC3'});
    for trl = 1:size(data.trial,2)
        data.trial{trl}(badchan,:) = mean(data.trial{trl}(replacechans,:));
    end      
end


% % check again
if ismember(pp, [5,6])
    figure;
    for ch = 1:64
        subplot(8,8,ch); hold on;
        plot(data.time{1}, data.trial{1}(ch,:), 'b');
        plot(data.time{100}, data.trial{100}(ch,:), 'r');
        title(data.label(ch));
        xlim([-1 1]); ylim([-100 100]);
    end
end


% % resample manually
if data.fsample == 1024 resamplefactor = 4; elseif data.fsample == 512 resamplefactor = 2; end; % 1024, to 256.
data.fsample = data.fsample ./ resamplefactor;
for trl = 1:length(data.trial)
data.time{trl} = data.time{trl}(1:resamplefactor:end);
data.trial{trl} = data.trial{trl}(:,1:resamplefactor:end);
end

% cfg = [];
% cfg.resamplefs = 200;
% cfg.detrend = 'yes';
% data = ft_resampledata(cfg, data);


% % add condition order to trial info
behdata = load(param.beh);
param.blockorder = behdata.deg_matrix_allp(:,19); % condi: 1-none100%, 2-none80%, 3-none60%, 4-im100%.
param.cuevalidity = behdata.deg_matrix_allp(:,7); % 7-cuetype_list1Valid 2Invalid
data.trialinfo(:,2) = param.blockorder;
data.trialinfo(:,3) = param.cuevalidity;


% % save epoched data
save([param.pathnewlen, '/processed_data/' 'epoched_data_eeg', '__', param.subjName], 'data');

clear data
close all

end
