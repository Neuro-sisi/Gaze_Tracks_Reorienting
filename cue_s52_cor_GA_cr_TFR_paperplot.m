%% Plot for papers

clc;clear;close all


% % parameters
laplacian           = 1;
chanclusters        = 1; % 0 = PO7/PO8, 1 = PO7+5, PO8+5
window              = 1; % 1 = 300 ms, 2 = 500 ms


pp2do               = [1:23];


plotRunningTopos    = 0;
plotSingleSubs      = 0;
plotsingleChan      = 0;
colormap2use        = fliplr(brewermap(100, 'RdBu'));

xlim2plot           = [-200 1500];
xlim2plot_probe     = [1300 2400];


% % load
if laplacian    toadd1 = '_laplacian';     else toadd1 = ''; end
if window == 1  toadd2 = ''; elseif window == 2 toadd2 = '_500mswindow'; end;
if chanclusters toadd3 = '_chanclusters';  else toadd3 = ''; end;

s = 0;
for pp = pp2do;
    s = s+1;
    param = getSubjParam(pp);

    disp(['getting data from ' param.subjName]);

    load([param.pathnewlen, 'saved_data/cue_cor_timefreq', toadd1,toadd2,toadd3, '__' param.subjName], 'timefreq','ci','cic','cii','ci_allchan');


    d1(s,:,:,:) = timefreq.imper_vs_high;
    d2(s,:,:,:) = timefreq.high_vs_low;
    d3(s,:,:,:) = ci.data;
    d4(s,:,:,:) = timefreq.imper100;
    d5(s,:,:,:) = timefreq.high100;
    d6(s,:,:,:) = timefreq.med80;
    d7(s,:,:,:) = timefreq.low60;

    d8(s,:,:,:) = cic.data;
    d9(s,:,:,:) = cii.data;

end


%% put back into structure
ci.label{1,5} = 'im-high';
ci.label{1,6} = 'high-low';
timefreq.imper_vs_high = squeeze(mean(d1));
timefreq.high_vs_low   = squeeze(mean(d2));
ci.data                = squeeze(mean(d3));
timefreq.imper100      = squeeze(mean(d4));
timefreq.high100       = squeeze(mean(d5));
timefreq.med80         = squeeze(mean(d6));
timefreq.low60         = squeeze(mean(d7));
cic.data                = squeeze(mean(d8));
cii.data                = squeeze(mean(d9));

save(fullfile('CueValid_v1_cor_cuelocked_alpha_results_timefreq.mat'))


%% Permutation test

%% later alpha cluster test--cuelocked 4con
fsel = ci.freq >= 8 & ci.freq <= 12;

data_cond1 = squeeze(mean(d3(:,1,fsel,:),3));
data_cond2 = squeeze(mean(d3(:,2,fsel,:),3));
data_cond3 = squeeze(mean(d3(:,3,fsel,:),3));
data_cond4 = squeeze(mean(d3(:,4,fsel,:),3));
data_zero = zeros(size(data_cond1)); % compare with 0

statcfg.xax = ci.time;
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.025; % two-sided
statcfg.nsub = s;
statcfg.statMethod = 'montecarlo'; % or 'analytic'; 

% cluster test
stat1 = frevede_ftclusterstat1D(statcfg, data_cond1, data_zero);
stat2 = frevede_ftclusterstat1D(statcfg, data_cond2, data_zero);
stat3 = frevede_ftclusterstat1D(statcfg, data_cond3, data_zero);
stat4 = frevede_ftclusterstat1D(statcfg, data_cond4, data_zero);

sigclurange1 = statcfg.xax(find(stat1.negclusterslabelmat==1));
sigclurange2 = statcfg.xax(find(stat2.negclusterslabelmat==1));
sigclurange3 = statcfg.xax(find(stat3.negclusterslabelmat==1));
sigclurange4 = statcfg.xax(find(stat4.negclusterslabelmat==1));

save(fullfile('CueValid_v1_cor_cuelocked_alpha_results_timefreq_permtest_alpha.mat'))


%% later Theta & Beta cluster test--cuelocked 4con
fsel1 = ci.freq >= 3 & ci.freq <= 7;
fsel2 = ci.freq >= 13 & ci.freq <= 30;

% define time range
tpnt_s = 136; % 0-1500 ms
tpnt_e = 286;

% define permtest params
statcfg.xax = ci.time(tpnt_s:tpnt_e);
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.025; % two-sided
statcfg.nsub = s;
statcfg.statMethod = 'montecarlo'; % or 'analytic'; 

% prep data--theta
data_cond11 = squeeze(mean(d3(:,1,fsel1,tpnt_s:tpnt_e),3));
data_cond12 = squeeze(mean(d3(:,2,fsel1,tpnt_s:tpnt_e),3));
data_cond13 = squeeze(mean(d3(:,3,fsel1,tpnt_s:tpnt_e),3));
data_cond14 = squeeze(mean(d3(:,4,fsel1,tpnt_s:tpnt_e),3));
data_zero = zeros(size(data_cond11)); % compare with 0

% cluster test
stat11 = frevede_ftclusterstat1D(statcfg, data_cond11, data_zero);
stat12 = frevede_ftclusterstat1D(statcfg, data_cond12, data_zero);
stat13 = frevede_ftclusterstat1D(statcfg, data_cond13, data_zero);
stat14 = frevede_ftclusterstat1D(statcfg, data_cond14, data_zero);

sigclurange11 = statcfg.xax(find(stat11.negclusterslabelmat==1));
sigclurange12 = statcfg.xax(find(stat12.negclusterslabelmat==1));
sigclurange13 = statcfg.xax(find(stat13.negclusterslabelmat==1));
sigclurange14 = statcfg.xax(find(stat14.negclusterslabelmat==1));


% prep data--beta
data_cond21 = squeeze(mean(d3(:,1,fsel2,tpnt_s:tpnt_e),3));
data_cond22 = squeeze(mean(d3(:,2,fsel2,tpnt_s:tpnt_e),3));
data_cond23 = squeeze(mean(d3(:,3,fsel2,tpnt_s:tpnt_e),3));
data_cond24 = squeeze(mean(d3(:,4,fsel2,tpnt_s:tpnt_e),3));
data_zero = zeros(size(data_cond21)); % compare with 0

% cluster test
stat21 = frevede_ftclusterstat1D(statcfg, data_cond21, data_zero);
stat22 = frevede_ftclusterstat1D(statcfg, data_cond22, data_zero);
stat23 = frevede_ftclusterstat1D(statcfg, data_cond23, data_zero);
stat24 = frevede_ftclusterstat1D(statcfg, data_cond24, data_zero);

sigclurange21 = statcfg.xax(find(stat21.negclusterslabelmat==1));
sigclurange22 = statcfg.xax(find(stat22.negclusterslabelmat==1));
sigclurange23 = statcfg.xax(find(stat23.negclusterslabelmat==1));
sigclurange24 = statcfg.xax(find(stat24.negclusterslabelmat==1));

save(fullfile('CueValid_v1_cor_cuelocked_alpha_results_timefreq_permtest_thetabeta.mat'))



%% plot ------------------------------------------------------------------------------------

%% Plot Cue-locked

%% Plot lateralized alpha time-series--Cuelocked
fsel = ci.freq >= 8 & ci.freq <= 12;
% title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
title2use = {'Im100%','100%','80%','60%'};

% contra vs. ipsi alpha
figure('Name','paperplot_CueValid_v1_cor_cuelocked_alpha_conipsi_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    subplot(1,4,sp); hold on; title(title2use{sp});
    [p1] = frevede_errorbarplot(ci.time, squeeze(mean(d8(:,sp,fsel,:),3)), [1,0.6,0], 'se');
    [p2] = frevede_errorbarplot(ci.time, squeeze(mean(d9(:,sp,fsel,:),3)), [0,1,1], 'se');
    xlim(xlim2plot); xline(0,'--k','Cue Onset');
    %   plot(xlim, [0,0], '--k');  ylim([-10 10]);
    set(gca, 'YDir','reverse'); ylabel('Power'), xlabel('Time (ms)');
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
end
legend([p1, p2], {'Contralateral','Ipsilateral'},'FontSize',12,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_conipsi_sepcon_s22', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_conipsi_sepcon_s22', 'epsc')


%lateralized alpha: seperate conditions
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];
figure('Name','paperplot_CueValid_v1_cor_cuelocked_alpha_later_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    subplot(1,4,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp,fsel,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-10 5]);
    xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
end
% legend({'Lateralized Alpha'},'FontSize',8,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_sepcon_s22', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_sepcon_s22', 'epsc')


%later alpha: overlay
% col2use = [0,0,1; 1,0,0; 1,0,1; 0,1,1;];
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];
figure('Name','paperplot_CueValid_v1_cor_cuelocked_alpha_later_overlay_4con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])
hold on;
% for sp = 1:4
p1 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,1,fsel,:),3)), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,2,fsel,:),3)), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,3,fsel,:),3)), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,4,fsel,:),3)), [col2use(4,:)], 'se');
% end
legend([p1,p2,p3,p4],title2use{1:4},'FontSize',12,'AutoUpdate','off','Box','off');
plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-10 5]);
xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
title('Lateralized Alpha');

saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_overlay_4con_s22', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_overlay_4con_s22', 'epsc')


%% Bar plot--averaged lateralized alpha
fsel = ci.freq >= 8 & ci.freq <= 12;
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];

timeproon = 0;
timerange = [400+timeproon 800+timeproon];
%     timerange = [300+timeproon 800+timeproon];

timediff1 = abs(ci.time-timerange(1));
pnt_start = find(timediff1==min(timediff1));
timediff2 = abs(ci.time-timerange(2));
pnt_end = find(timediff2==min(timediff2));
clear timediff1 timediff2

% average TFA across time
aveerp = [];
aveerp(:,1) = mean(squeeze(mean(d3(:,1,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,2) = mean(squeeze(mean(d3(:,2,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,3) = mean(squeeze(mean(d3(:,3,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,4) = mean(squeeze(mean(d3(:,4,fsel,pnt_start:pnt_end),4)),2);

% bar plot
figure('Name','CueValid_v1_cor_cuelocked_alpha_later_ave_4con_bar400800','NumberTitle','off', 'Color','white'),
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
ylim([-10 0]);set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Block Type');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3);

% saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_ave_4con_bar400800_s23', 'jpg');
% set(gcf,'renderer','painters')
% saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_ave_4con_bar400800_s23', 'epsc')
% clear aveerp b

saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_ave_4con_bar400800_new_s23', 'jpg');
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_alpha_later_ave_4con_bar400800_new_s23', 'epsc')
clear aveerp b


%% Plot Lateralized Theta & Alpha & Beta time-series--Cuelocked
% title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
title2use = {'Im100%','100%','80%','60%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];

fsel1 = ci.freq >= 3 & ci.freq <= 7;
fsel = ci.freq >= 8 & ci.freq <= 12;
fsel2 = ci.freq >= 13 & ci.freq <= 30;

% Lateralized Theta: 
% seperate conditions
figure('Name','paperplot_CueValid_v1_cor_cuelocked_thalbe_later_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    subplot(1,4,sp); hold on; title(title2use{sp});
    p1 = frevede_errorbarplot_dashline(ci.time, squeeze(mean(d3(:,sp,fsel1,:),3)), [col2use(sp,:)], 'se');
    p2 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp,fsel,:),3)), [col2use(sp,:)], 'se');
    p3 = frevede_errorbarplot_dotline(ci.time, squeeze(mean(d3(:,sp,fsel2,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-10 5]);
    xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
end
legend([p1,p2,p3],{'Lateralized Theta','Lateralized Alpha','Lateralized Beta'},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_thalbe_later_sepcon_s22', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValid_v1_cor_cuelocked_thalbe_later_sepcon_s22', 'epsc')


%% plot spectrum -- c vs i -- cue-locked
cfg             = [];
cfg.colorbar    = 'no';
cfg.zlim        = [-8 8];
cfg.xlim        = xlim2plot;
cfg.parameter   = 'data';
cfg.figure      = 'gcf';

ci.dimord = 'chan_freq_time';
% title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
title2use = {'Im100%','100%','80%','60%'};

% compare: 4 cue validity conditions
figure('Name','paperplot_CueValid_v1_cor_tfa_cvsi_spectrum_cuelocked_4con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for channel = 1:4
    cfg.channel     = ci.label(channel);
    subplot(1,4,channel);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
    xline(0,'--k','Cue Onset')
    title(title2use{channel});
    ylabel('Frequency (Hz)'); xlabel('Time (ms)');
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
    colorbar;
end

% saveas(gcf, 'paperplot_CueValid_v1_cor_tfa_cvsi_spectrum_cuelocked_4con_s23', 'jpg');
% set(gcf,'renderer','painters')
% saveas(gcf, 'paperplot_CueValid_v1_cor_tfa_cvsi_spectrum_cuelocked_4con_s23', 'epsc')

saveas(gcf, 'paperplot_CueValid_v1_cor_tfa_cvsi_spectrum_cuelocked_4con_s23_bar', 'jpg');
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValid_v1_cor_tfa_cvsi_spectrum_cuelocked_4con_s23_bar', 'epsc')



