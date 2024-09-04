%% Eye-tracking data analysis--Cue Validity Experiment

% 24 July 2023 Sisi Wang

% 1 Apr 2024, correction
% the cue side is locked to test item, reverse cueside after cue for invalid cue trials

%% Step3b--grand average plots of gaze-shift (saccade) results

% time-locked to retrocue onset: 0
% baseline: -1500--1000 ms
% artifact rejection: 0-1000 ms

%% start clean
clear; clc; close all;

% % add fieldtrip & path
% addpath '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/Eyedata_Analysis/toolbox/fieldtrip-20201023'
% ft_defaults

% addpath '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/Data_ana/Eyedata_Ana/';
% cd '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/Data_ana/Eyedata_Ana/';
addpath '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/recue_valid/Data_ana/Eyedata_Ana/';
cd '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/recue_valid/Data_ana/Eyedata_Ana/';


% % set loops
% sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';};
% sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'30';'31';'32';'33';'34';};
% sub27 & 29 not included due to eyetracker failer for one part

% select subs to average
% pp2do = 1:length(sub_list);
pp2do =  [1:16 18 20:24];

oneOrTwoD       = 1;        oneOrTwoD_options = {'_1D','_2D'};
nsmooth         = 25;
% nsmooth         = 50; % smooth data 
plotSinglePps   = 1;
xlimtoplot      = [-200 1500];


%% load and aggregate the data from all pp
s = 0;
for pp = pp2do
    s = s+1;

    % get participant data
    param = getSubjParam(pp);

    % load
    disp(['getting data from participant ', param.subjName]);
    load([param.resupath, 'results/saved_data/cor_saccadeEffects', oneOrTwoD_options{oneOrTwoD} '_', param.subjName], 'saccade','saccadesize');

    % smooth?
    if nsmooth > 0
        for x1 = 1:size(saccade.toward,1)
            saccade.toward(x1,:)  = gsmooth(squeeze(saccade.toward(x1,:)), nsmooth);
            saccade.away(x1,:)    = gsmooth(squeeze(saccade.away(x1,:)), nsmooth);
            saccade.effect(x1,:)  = gsmooth(squeeze(saccade.effect(x1,:)), nsmooth);
        end
        % also smooth saccadesize data over time.
        for x1 = 1:size(saccadesize.toward,1)
            for x2 = 1:size(saccadesize.toward,2)
                saccadesize.toward(x1,x2,:) = gsmooth(squeeze(saccadesize.toward(x1,x2,:)), nsmooth);
                saccadesize.away(x1,x2,:)   = gsmooth(squeeze(saccadesize.away(x1,x2,:)), nsmooth);
                saccadesize.effect(x1,x2,:) = gsmooth(squeeze(saccadesize.effect(x1,x2,:)), nsmooth);
            end
        end
    end

    % put into matrix, with pp as first dimension
    d1(s,:,:) = saccade.toward;
    d2(s,:,:) = saccade.away;
    d3(s,:,:) = saccade.effect;

    d4(s,:,:,:) = saccadesize.toward;
    d5(s,:,:,:) = saccadesize.away;
    d6(s,:,:,:) = saccadesize.effect;
end


% % make GA for the saccadesize fieldtrip structure data, to later plot as "time-frequency map" with fieldtrip. For timecourse data, we directly plot from d structures above.
saccadesize.toward = squeeze(mean(d4));
saccadesize.away   = squeeze(mean(d5));
saccadesize.effect = squeeze(mean(d6));


%% saccade cluster test -- [0-1500ms]
% timerange = xlimtoplot;
timerange = [0 1500];
timediff1 = abs(saccade.time-timerange(1));
pnt_start = find(timediff1==min(timediff1));
timediff2 = abs(saccade.time-timerange(2));
pnt_end = find(timediff2==min(timediff2));
clear timediff1 timediff2

data_cond1 = squeeze(d3(:,5,pnt_start:pnt_end));
data_cond2 = squeeze(d3(:,2,pnt_start:pnt_end));
data_cond3 = squeeze(d3(:,3,pnt_start:pnt_end));
data_cond4 = squeeze(d3(:,4,pnt_start:pnt_end));
data_zero = zeros(size(data_cond1)); % compare with 0

statcfg.xax = saccade.time(pnt_start:pnt_end);
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.025; % two-sided (0.05 if one-sided)
statcfg.nsub = s;
statcfg.statMethod = 'montecarlo'; % or 'analytic'; 

% cluster test
stat1 = frevede_ftclusterstat1D(statcfg, data_cond1, data_zero);
stat2 = frevede_ftclusterstat1D(statcfg, data_cond2, data_zero);
stat3 = frevede_ftclusterstat1D(statcfg, data_cond3, data_zero);
stat4 = frevede_ftclusterstat1D(statcfg, data_cond4, data_zero);

signegclurange1 = statcfg.xax(find(stat1.negclusterslabelmat==1));
signegclurange2 = statcfg.xax(find(stat2.negclusterslabelmat==1));
signegclurange3 = statcfg.xax(find(stat3.negclusterslabelmat==1));
signegclurange4 = statcfg.xax(find(stat4.negclusterslabelmat==1));
sigposclurange1 = statcfg.xax(find(stat1.posclusterslabelmat==1));
sigposclurange2 = statcfg.xax(find(stat2.posclusterslabelmat==1));
sigposclurange3 = statcfg.xax(find(stat3.posclusterslabelmat==1));
sigposclurange4 = statcfg.xax(find(stat4.posclusterslabelmat==1));

save(fullfile('CueVa_v1_cor_saccade_effect_permutation_results_1500.mat'))


%% saccade cluster test -- [0-1000ms]
% timerange = xlimtoplot;
timerange = [0 1000];
timediff1 = abs(saccade.time-timerange(1));
pnt_start = find(timediff1==min(timediff1));
timediff2 = abs(saccade.time-timerange(2));
pnt_end = find(timediff2==min(timediff2));
clear timediff1 timediff2

data_cond1 = []; data_cond2 = [];  data_cond3 = []; data_cond4 = []; 
data_cond1 = squeeze(d3(:,5,pnt_start:pnt_end));
data_cond2 = squeeze(d3(:,2,pnt_start:pnt_end));
data_cond3 = squeeze(d3(:,3,pnt_start:pnt_end));
data_cond4 = squeeze(d3(:,4,pnt_start:pnt_end));
data_zero = zeros(size(data_cond1)); % compare with 0

statcfg.xax = saccade.time(pnt_start:pnt_end);
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.025; % two-sided (0.05 if one-sided)
statcfg.nsub = s;
statcfg.statMethod = 'montecarlo'; % or 'analytic'; 

% cluster test
stat1 = frevede_ftclusterstat1D(statcfg, data_cond1, data_zero);
stat2 = frevede_ftclusterstat1D(statcfg, data_cond2, data_zero);
stat3 = frevede_ftclusterstat1D(statcfg, data_cond3, data_zero);
stat4 = frevede_ftclusterstat1D(statcfg, data_cond4, data_zero);

signegclurange1 = statcfg.xax(find(stat1.negclusterslabelmat==1));
signegclurange2 = statcfg.xax(find(stat2.negclusterslabelmat==1));
signegclurange3 = statcfg.xax(find(stat3.negclusterslabelmat==1));
signegclurange4 = statcfg.xax(find(stat4.negclusterslabelmat==1));
sigposclurange1 = statcfg.xax(find(stat1.posclusterslabelmat==1));
sigposclurange2 = statcfg.xax(find(stat2.posclusterslabelmat==1));
sigposclurange3 = statcfg.xax(find(stat3.posclusterslabelmat==1));
sigposclurange4 = statcfg.xax(find(stat4.posclusterslabelmat==1));

save(fullfile('CueVa_v1_cor_saccade_effect_permutation_results_1000.mat'))


%% plot grand average

%% Paper plot--V2

%% plot grand average data patterns of interest, with error bars
% title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
title2use = {'Im100%','100%','80%','60%'};
% col2use = [0,0,1; 1,0,0; 1,0,1; 0,1,1;];
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];

% saccade right vs. left cues--color
figure('Name','paperplot_CueVa_v1_cor_saccade_LeftvsRight_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    if sp==1 con=5; elseif sp==2 con=2; elseif sp==3 con=3; elseif sp==4 con=4; end
    subplot(1,4,sp); hold on; title(title2use{sp});
    p1 = frevede_errorbarplot(saccade.time, squeeze(d1(:,con,:)), col2use(sp,:), 'se');
    p2 = frevede_errorbarplot(saccade.time, squeeze(d2(:,con,:)), [0.5,0.5,0.5], 'se');
    plot(xlim, [0,0], '--k');xlim(xlimtoplot); ylim([0 1.5]);
    xline(0,'--k','Cue Onset');
    ylabel('Saccade Rate'); xlabel('Time (ms)'); 
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
end
legend([p1, p2], {'toward','away'},'FontSize',12,'AutoUpdate','off','Box','off'); 

saveas(gcf, 'paperplot_CueVa_v1_cor_saccade_LeftvsRight_sepcon_s22_v2', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueVa_v1_cor_saccade_LeftvsRight_sepcon_s22_v2', 'epsc')


%saccade effect (towardness): toward-away seperate
figure('Name','paperplot_CueVa_v1_cor_saccade_effect_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    if sp==1 con=5; elseif sp==2 con=2; elseif sp==3 con=3; elseif sp==4 con=4; end
    subplot(1,4,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(saccade.time, squeeze(d3(:,con,:)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlimtoplot); ylim([-0.5 1]);
    xline(0,'--k','Cue Onset');
    ylabel('Saccade Effect'); xlabel('Time (ms)'); 
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
end
legend({'toward minus away'},'FontSize',12,'AutoUpdate','off','Box','off'); 
saveas(gcf, 'paperplot_CueVa_v1_cor_saccade_effect_sepcon_color_s22_v2', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueVa_v1_cor_saccade_effect_sepcon_color_s22_v2', 'epsc')


%saccade effect (towardness): overlay
figure('Name','paperplot_CueVa_v1_cor_saccade_effect_overlay_4con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])
hold on;
p2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,5,:)), [col2use(1,:)], 'se');
p3 = frevede_errorbarplot(saccade.time, squeeze(d3(:,2,:)), [col2use(2,:)], 'se');
p4 = frevede_errorbarplot(saccade.time, squeeze(d3(:,3,:)), [col2use(3,:)], 'se');
p5 = frevede_errorbarplot(saccade.time, squeeze(d3(:,4,:)), [col2use(4,:)], 'se');
legend([p2,p3,p4,p5], title2use{1:4},'FontSize',12,'AutoUpdate','off','Box','off'); 
plot(xlim, [0,0], '--k');xlim(xlimtoplot);ylim([-0.5 1]);
xline(0,'--k','Cue Onset');
ylabel('Saccade Effect'); xlabel('Time (ms)');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
% title('Saccade Effect'); 

saveas(gcf, 'paperplot_CueVa_v1_cor_saccade_effect_overlay_4con_s22_v2', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueVa_v1_cor_saccade_effect_overlay_4con_s22_v2', 'epsc')



%% bar plot--averaged saccade rate
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];

timeproon = 0;
timerange = [200+timeproon 600+timeproon];
timediff1 = abs(saccade.time-timerange(1));
pnt_start = find(timediff1==min(timediff1));
timediff2 = abs(saccade.time-timerange(2));
pnt_end = find(timediff2==min(timediff2));
clear timediff1 timediff2

% average ERP across time
aveerp = [];
aveerp(:,1) = mean(squeeze(d3(:,5,pnt_start:pnt_end)),2);
aveerp(:,2) = mean(squeeze(d3(:,2,pnt_start:pnt_end)),2);
aveerp(:,3) = mean(squeeze(d3(:,3,pnt_start:pnt_end)),2);
aveerp(:,4) = mean(squeeze(d3(:,4,pnt_start:pnt_end)),2);

[~,P_12,~,STATS_12] = ttest(aveerp(:,1),aveerp(:,2));
[~,P_23,~,STATS_23] = ttest(aveerp(:,2),aveerp(:,3));
[~,P_34,~,STATS_34] = ttest(aveerp(:,3),aveerp(:,4));

% bar plot
figure('Name','paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 400 300])
set(gcf,'Position',[0 0 200 450])
hold on,
b = bar(1:4, mean(aveerp(:,:)));
b.LineWidth = 2;
b.FaceColor = 'flat';
for c = 1:4
    b.CData(c,:) = col2use(c,:);
end
hold on,
errorbar(1:4, mean(aveerp(:,:)), std(aveerp(:,:))./sqrt((size(aveerp,1))), '.k', 'LineWidth',2);
% xticks(1:4); xticklabels({'Imper100%','Infor100%','Infor80%','Infor60%'});
xticks(1:4); xticklabels({'Im100%','100%','80%','60%'});
ylim([0 0.6]);
ylabel('Saccade Rate'); xlabel('Block Type');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3);

% saveas(gcf, 'paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600_s22', 'jpg')
% set(gcf,'renderer','painters')
% saveas(gcf, 'paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600_s22', 'epsc')
% clear aveerp b

saveas(gcf, 'paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600_new_s22', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600_new_s22', 'epsc')
clear aveerp b


%% bar plot--averaged saccade rate--normalize mean
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];

timeproon = 0;
timerange = [200+timeproon 600+timeproon];
timediff1 = abs(saccade.time-timerange(1));
pnt_start = find(timediff1==min(timediff1));
timediff2 = abs(saccade.time-timerange(2));
pnt_end = find(timediff2==min(timediff2));
clear timediff1 timediff2

% average ERP across time
aveerp = [];
aveerp(:,1) = mean(squeeze(d3(:,5,pnt_start:pnt_end)),2);
aveerp(:,2) = mean(squeeze(d3(:,2,pnt_start:pnt_end)),2);
aveerp(:,3) = mean(squeeze(d3(:,3,pnt_start:pnt_end)),2);
aveerp(:,4) = mean(squeeze(d3(:,4,pnt_start:pnt_end)),2);

% normalize averaged gaze with mean across 4 conditions
aveerp_norm = aveerp/mean(mean(aveerp));

[~,P_12,~,STATS_12] = ttest(aveerp_norm(:,1),aveerp_norm(:,2));
[~,P_23,~,STATS_23] = ttest(aveerp_norm(:,2),aveerp_norm(:,3));
[~,P_34,~,STATS_34] = ttest(aveerp_norm(:,3),aveerp_norm(:,4));

% bar plot
figure('Name','paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600_norm','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])
hold on,
b = bar(1:4, mean(aveerp_norm(:,:)));
b.LineWidth = 2;
b.FaceColor = 'flat';
for c = 1:4
    b.CData(c,:) = col2use(c,:);
end
hold on,
errorbar(1:4, mean(aveerp_norm(:,:)), std(aveerp_norm(:,:))./sqrt((size(aveerp_norm,1))), '.k', 'LineWidth',2);
xticks(1:4); xticklabels({'Im100%','100%','80%','60%'});
ylim([0 3]); 
% yticklabels([0:1:3]);
ylabel('Normalized Saccade Rate'); xlabel('Cue Type');
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);

saveas(gcf, 'paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600_s22', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueVa_v1_cor_cuelocked_saccade_effect_4con_bar200600_norm_s22', 'epsc')
clear aveerp b


%% saccade size X time
colormap2use        = fliplr(brewermap(100, 'RdBu'));
% title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
title2use = {'Im100%','100%','80%','60%'};

cfg = [];
cfg.parameter = 'effect';
cfg.figure = 'gcf';
cfg.zlim = [-.1 .1];
cfg.xlim = xlimtoplot;
cfg.colorbar = 'off';

figure('Name','CueVa_v1_cor_saccadesize_sepcon','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 1600 250])
set(gcf,'Position',[0 0 1400 250])
for c = 1:4
    if c== 1 chan = 5; elseif c== 2 chan = 2; elseif c== 3 chan = 3; elseif c== 4 chan = 4; end
    cfg.channel = chan;
    subplot(1,4,c); ft_singleplotTFR(cfg, saccadesize);
    title(title2use{c}); 
    xline(0,'--k','Cue Onset');
    ylabel('Pixels'); xlabel('Time (ms)'); 
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
% colormap('jet');
colormap(colormap2use);

saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_nobar', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_nobar', 'epsc')

% saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_bar', 'jpg')
% set(gcf,'renderer','painters')
% saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_bar', 'epsc')


%% saccade size X time -- smaller plots
colormap2use        = fliplr(brewermap(100, 'RdBu'));
title2use = {'Imp100%','100%','80%','60%'};

cfg = [];
cfg.parameter = 'effect';
cfg.figure = 'gcf';
cfg.zlim = [-.1 .1];
cfg.xlim = xlimtoplot;
cfg.colorbar = 'off';

figure('Name','CueVa_v1_cor_saccadesize_sepcon','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 1600 250])
set(gcf,'Position',[0 0 1100 230])
for c = 1:4
    if c== 1 chan = 5; elseif c== 2 chan = 2; elseif c== 3 chan = 3; elseif c== 4 chan = 4; end
    cfg.channel = chan;
    subplot(1,4,c); ft_singleplotTFR(cfg, saccadesize);
    title(title2use{c}); 
    xline(0,'--k','Cue Onset');
    ylabel('Pixels'); xlabel('Time (ms)'); 
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
end
% colormap('jet');
colormap(colormap2use);

saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_nobar_new', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_nobar_new', 'epsc')

% saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_bar_new', 'jpg')
% set(gcf,'renderer','painters')
% saveas(gcf, 'paperplot_CueVa_v1_cor_saccadesize_sepcon_s22_bar_new', 'epsc')

