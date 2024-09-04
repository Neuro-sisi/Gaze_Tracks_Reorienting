%% Eye-tracking data analysis--Cue Validity Experiment

% 24 July 2023 Sisi Wang

% 1 Apr 2024, correction
% the cue side is locked to test item, reverse cueside after cue for invalid cue trials

%% Step3-- gaze-shift calculation

% time-locked to retrocue onset: 0
% baseline: -1500--1000 ms
% artifact rejection: 0-1000 ms 

%% start clean
clear; clc; close all;

% % add fieldtrip & path
addpath '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/Eyedata_Analysis/toolbox/fieldtrip-20201023'
ft_defaults

addpath '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/recue_valid/Data_ana/Eyedata_Ana/';

%% set loops
% sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';};
% sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'30';'31';'32';'33';'34';};
% sub27 & 29 not included due to eyetracker failer for one part

% for pp = 1:length(sub_list)
for pp = [1:16 18 20:24]

    oneOrTwoD       = 1; oneOrTwoD_options = {'_1D','_2D'};
    plotResults     = 0; % plot-1

    
    %% load epoched data of this participant data and concattenate the three parts
    param = getSubjParam(pp);
    x1 = load([param.resupath, 'results/epoched_data/CueValidity_v1_', param.subjName, '_condi1'], 'eyedata');
    x2 = load([param.resupath, 'results/epoched_data/CueValidity_v1_', param.subjName, '_condi2'], 'eyedata');
    x3 = load([param.resupath, 'results/epoched_data/CueValidity_v1_', param.subjName, '_condi3'], 'eyedata');
    x4 = load([param.resupath, 'results/epoched_data/CueValidity_v1_', param.subjName, '_condi4'], 'eyedata');

    % append 
    cfg = [];
    eyedata = ft_appenddata(cfg, x1.eyedata, x2.eyedata, x3.eyedata, x4.eyedata);
    clear x*
    
    %% add relevant behavioural file data
    behdata = load(param.log);
    eyedata.trialinfo(:,2) = behdata.deg_matrix_allc(:,8); % cueside-relative to test: 1-left, 2-right.
    eyedata.trialinfo(:,3) = behdata.deg_matrix_allc(:,19); % condi: 1-none100%, 2-none80%, 3-none60%, 4-im100%.
    eyedata.trialinfo(:,4) = behdata.deg_matrix_allc(:,7); % '7-cuetype_list1Valid2Invalid'

    % reverse cueside for invalid trials
    trialinfo = eyedata.trialinfo;
    for t = 1:size(trialinfo,1)
        if trialinfo(t,4) == 1
            trialinfo(t,5) = trialinfo(t,2); % cueside--relative to cue
        elseif trialinfo(t,4) == 2
            trialinfo(t,5) = 3-trialinfo(t,2); % reverse invalid trials
        end
    end
    eyedata.trialinfo = trialinfo; clear trialinfo
    

    %% only keep channels of interest
    cfg = [];
    cfg.channel = {'eyeX','eyeY'}; % only keep x & y axis
    eyedata = ft_selectdata(cfg, eyedata); % select x & y channels

    %% reformat such that all data in single matrix of trial x channel x time
    cfg = [];
    cfg.keeptrials = 'yes';
    tl = ft_timelockanalysis(cfg, eyedata); % realign the data: from trial*time cells into trial*channel*time?

    %% pixel to degree
    [dva_x, dva_y] = frevede_pixel2dva(squeeze(tl.trial(:,1,:)), squeeze(tl.trial(:,2,:)));
    tl.trial(:,1,:) = dva_x;
    tl.trial(:,2,:) = dva_y;

    %% selection vectors for conditions

    % cued item location
    cueL = ismember(tl.trialinfo(:,5), [1]);
    cueR = ismember(tl.trialinfo(:,5), [2]);

    % cue validity: 1-4
    NonE100 = ismember(tl.trialinfo(:,3), [1]);
    NonE80 = ismember(tl.trialinfo(:,3), [2]);
    NonE60 = ismember(tl.trialinfo(:,3), [3]);
    Im100 = ismember(tl.trialinfo(:,3), [4]);


    %% channels
    chX = ismember(tl.label, 'eyeX');
    chY = ismember(tl.label, 'eyeY');

    %% get gaze shifts using our custom function
    cfg = [];
    data_input = squeeze(tl.trial);
    time_input = tl.time*1000;

    if oneOrTwoD == 1         [shiftsX, velocity, times]             = PBlab_gazepos2shift_1D(cfg, data_input(:,chX,:), time_input);
    elseif oneOrTwoD == 2     [shiftsX,shiftsY, peakvelocity, times] = PBlab_gazepos2shift_2D(cfg, data_input(:,chX,:), data_input(:,chY,:), time_input);
    end

    %% get saccade shift data for trial-cateogrizing
    data_shift.shift = shiftsX;
    data_shift.time = times;
    data_shift.trialinfo = eyedata.trialinfo;

    % save
    save([param.resupath, 'results/saved_data/cor_saccadeshift', oneOrTwoD_options{oneOrTwoD} '_', param.subjName], 'data_shift');
    clear data_shift

    %% select usable gaze shifts
    minDisplacement = 0;
    maxDisplacement = 1000;

    if oneOrTwoD == 1     saccadesize = abs(shiftsX);
    elseif oneOrTwoD == 2 saccadesize = abs(shiftsX+shiftsY*1i);
    end
    shiftsL = shiftsX<0 & (saccadesize>minDisplacement & saccadesize<maxDisplacement);
    shiftsR = shiftsX>0 & (saccadesize>minDisplacement & saccadesize<maxDisplacement);

    %% get relevant contrasts out
    saccade = [];
    saccade.time = times;
    saccade.label = {'all','NonE100','NonE80','NonE60','Im100'};

    for selection = 1:5
        if     selection == 1  sel = ones(size(cueL));
        elseif selection == 2  sel = NonE100;
        elseif selection == 3  sel = NonE80;
        elseif selection == 4  sel = NonE60;
        elseif selection == 5  sel = Im100;
        end
        saccade.toward(selection,:) =  (mean(shiftsL(sel&cueL,:)) + mean(shiftsR(sel&cueR,:))) ./ 2;
        saccade.away(selection,:)  =   (mean(shiftsL(sel&cueR,:)) + mean(shiftsR(sel&cueL,:))) ./ 2;
    end

    % add towardness field
    saccade.effect = (saccade.toward - saccade.away);


    %% smooth and turn to Hz
    integrationwindow = 100; % window over which to integrate saccade counts
    saccade.toward = smoothdata(saccade.toward,2,'movmean',integrationwindow)*1000; % *1000 to get to Hz, given 1000 samples per second.
    saccade.away   = smoothdata(saccade.away,2,  'movmean',integrationwindow)*1000;
    saccade.effect = smoothdata(saccade.effect,2,'movmean',integrationwindow)*1000;

    %% plot
    if plotResults
        figure;    for sp = 1:5 subplot(2,3,sp); hold on; plot(saccade.time, saccade.toward(sp,:), 'r'); plot(saccade.time, saccade.away(sp,:), 'b'); title(saccade.label(sp)); legend({'toward','away'}); end
        figure;    for sp = 1:5 subplot(2,3,sp); hold on; plot(saccade.time, saccade.effect(sp,:), 'k'); plot(xlim, [0,0], '--k');                    title(saccade.label(sp)); legend({'effect'}); end
        figure;                                   hold on; plot(saccade.time, saccade.effect([1:5],:)); plot(xlim, [0,0], '--k'); legend(saccade.label([1:5]));
        drawnow;
    end

    %% also get as function of saccade size - identical as above, except with extra loop over saccade size.
    binsize = 0.5;
    halfbin = binsize/2;

    saccadesize = [];
    saccadesize.dimord = 'chan_freq_time';
    saccadesize.freq = halfbin:0.1:6-halfbin; % shift sizes, as if "frequency axis" for time-frequency plot
    saccadesize.time = times;
    saccadesize.label = {'all','NonE100','NonE80','NonE60','Im100'};

    cnt = 0;
    for sz = saccadesize.freq;
        cnt = cnt+1;
        shiftsL = [];
        shiftsR = [];
        shiftsL = shiftsX<-sz+halfbin & shiftsX > -sz-halfbin; % left shifts within this range
        shiftsR = shiftsX>sz-halfbin  & shiftsX < sz+halfbin; % right shifts within this range

        for selection = 1:5
            if     selection == 1  sel = ones(size(cueL));
            elseif selection == 2  sel = NonE100;
            elseif selection == 3  sel = NonE80;
            elseif selection == 4  sel = NonE60;
            elseif selection == 5  sel = Im100;
            end
            saccadesize.toward(selection,cnt,:) = (mean(shiftsL(cueL&sel,:)) + mean(shiftsR(cueR&sel,:))) ./ 2;
            saccadesize.away(selection,cnt,:) =   (mean(shiftsL(cueR&sel,:)) + mean(shiftsR(cueL&sel,:))) ./ 2;
        end
    end
    
    % add towardness field
    saccadesize.effect = (saccadesize.toward - saccadesize.away);


    %% smooth and turn to Hz
    integrationwindow = 100; % window over which to integrate saccade counts
    saccadesize.toward = smoothdata(saccadesize.toward,3,'movmean',integrationwindow)*1000; % *1000 to get to Hz, given 1000 samples per second.
    saccadesize.away   = smoothdata(saccadesize.away,3,  'movmean',integrationwindow)*1000;
    saccadesize.effect = smoothdata(saccadesize.effect,3,'movmean',integrationwindow)*1000;

    if plotResults
        cfg = [];
        cfg.parameter = 'effect';
        cfg.figure = 'gcf';
        cfg.zlim = [-.1 .1];
        figure;
        for chan = 1:5
            cfg.channel = chan;
            subplot(2,3,chan); ft_singleplotTFR(cfg, saccadesize);
        end
        colormap('jet');
        drawnow;
    end

    %% save
    save([param.resupath, 'results/saved_data/cor_saccadeEffects', oneOrTwoD_options{oneOrTwoD} '_', param.subjName], 'saccade','saccadesize');

    %% close loops
end % end pp loop


