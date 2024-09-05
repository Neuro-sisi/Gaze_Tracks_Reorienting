clc;clear 
close all

%% parameters
% loops
for laplacian           = 1;
for chanclusters        = 1; % 0 = PO7/PO8, 1 = PO7+5, PO8+5
for window              = 1; % 1 = 300 ms, 2 = 500 ms


pp2do               = [1:23];

maskbysignT         = 0;

plotRunningTopos    = 0;
plotSingleSubs      = 0;
plotsingleChan      = 0;
colormap2use        = fliplr(brewermap(100, 'RdBu'));
xlim2plot           = [-500 2500];

%% load
if laplacian    toadd1 = '_laplacian';     else toadd1 = ''; end
if window == 1  toadd2 = ''; elseif window == 2 toadd2 = '_500mswindow'; end;
if chanclusters toadd3 = '_chanclusters';  else toadd3 = ''; end;

s = 0;
for pp = pp2do;
    s = s+1;
    param = getSubjParam(pp);
    
    disp(['getting data from ' param.subjName]);
    
    load([param.pathnewlen, 'saved_data/timefreq', toadd1,toadd2,toadd3, '__' param.subjName], 'timefreq','ci');
    
    
    d1(s,:,:,:) = timefreq.imper_vs_high;
    d2(s,:,:,:) = timefreq.high_vs_low;    
    d3(s,:,:,:) = ci.data;    
    d4(s,:,:,:) = timefreq.valid;
    d5(s,:,:,:) = timefreq.invalid;    
    d6(s,:,:,:) = timefreq.validim100;
    d7(s,:,:,:) = timefreq.validinfor100;
    d8(s,:,:,:) = timefreq.validinfor80;
    d9(s,:,:,:) = timefreq.validinfor60;
    d10(s,:,:,:) = timefreq.validinfor40;
    d11(s,:,:,:) = timefreq.validinfor20;
    d12(s,:,:,:) = timefreq.imper100;
    d13(s,:,:,:) = timefreq.high100;
    d14(s,:,:,:) = timefreq.med80;
    d15(s,:,:,:) = timefreq.low60;
    d16(s,:,:,:) = timefreq.valid_vs_invalid;
    
end


%% put back into structure
timefreq.imper_vs_high = squeeze(mean(d1));
timefreq.high_vs_low   = squeeze(mean(d2));
ci.data                = squeeze(mean(d3));
timefreq.valid         = squeeze(mean(d4));
timefreq.invalid       = squeeze(mean(d5));
timefreq.validim100    = squeeze(mean(d6));
timefreq.validinfor100 = squeeze(mean(d7));
timefreq.validinfor80  = squeeze(mean(d8));
timefreq.validinfor60  = squeeze(mean(d9));
timefreq.validinfor40  = squeeze(mean(d10));
timefreq.validinfor20  = squeeze(mean(d11));
timefreq.imper100      = squeeze(mean(d12));
timefreq.high100       = squeeze(mean(d13));
timefreq.med80         = squeeze(mean(d14));
timefreq.low60         = squeeze(mean(d15));
timefreq.valid_vs_invalid = squeeze(mean(d16));


if maskbysignT
   h = ttest(d1);  timefreq.imper_vs_high = timefreq.imper_vs_high.*squeeze(h);
   h = ttest(d2);  timefreq.high_vs_low   = timefreq.high_vs_low.*squeeze(h);
   h = ttest(d3);  ci.data                = ci.data.*squeeze(h);    
end

%% plot ------------------------------------------------------------------------------------

%% singleplot--plot single channel
if plotsingleChan
cfg             = [];
cfg.zlim        = [-15 15];
cfg.xlim        = xlim2plot;
cfg.figure      = 'gcf';

cfg.parameter   = 'imper_vs_high'; 
figure;
subplot(3,3,1); cfg.channel = 'AFz'; ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - AFz']);
subplot(3,3,2); cfg.channel = 'Fz'; ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - Fz']);
subplot(3,3,3); cfg.channel = 'FCz'; ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - FCz']);
subplot(3,3,4); cfg.channel = 'Cz'; ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - Cz']);
subplot(3,3,5); cfg.channel = 'CPz'; ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - CPz']);
subplot(3,3,6); cfg.channel = 'Pz';  ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - Pz']);
subplot(3,3,7); cfg.channel = 'POz';  ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - POz']);
subplot(3,3,8); cfg.channel = 'Oz';  ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - Oz']);
subplot(3,3,9); cfg.channel = 'C3';  ft_singleplotTFR(cfg, timefreq);  title(['imper-vs-high', ' - C3']);
colormap(colormap2use);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_compare_impervshigh_s23', toadd1,toadd2,toadd3], '.jpg']);

cfg.parameter   = 'high_vs_low'; 
figure;
subplot(3,3,1); cfg.channel = 'AFz'; ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - AFz']);
subplot(3,3,2); cfg.channel = 'Fz'; ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - Fz']);
subplot(3,3,3); cfg.channel = 'FCz'; ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - FCz']);
subplot(3,3,4); cfg.channel = 'Cz'; ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - Cz']);
subplot(3,3,5); cfg.channel = 'CPz'; ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - CPz']);
subplot(3,3,6); cfg.channel = 'Pz';  ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - Pz']);
subplot(3,3,7); cfg.channel = 'POz';  ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - POz']);
subplot(3,3,8); cfg.channel = 'Oz';  ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - Oz']);
subplot(3,3,9); cfg.channel = 'C3';  ft_singleplotTFR(cfg, timefreq);  title(['high-vs-low', ' - C3']);
colormap(colormap2use);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_compare_highvslow_s23', toadd1,toadd2,toadd3], '.jpg']);

cfg.parameter   = 'valid_vs_invalid'; 
figure;
subplot(3,3,1); cfg.channel = 'AFz'; ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - AFz']);
subplot(3,3,2); cfg.channel = 'Fz'; ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - Fz']);
subplot(3,3,3); cfg.channel = 'FCz'; ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - FCz']);
subplot(3,3,4); cfg.channel = 'Cz'; ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - Cz']);
subplot(3,3,5); cfg.channel = 'CPz'; ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - CPz']);
subplot(3,3,6); cfg.channel = 'Pz';  ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - Pz']);
subplot(3,3,7); cfg.channel = 'POz';  ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - POz']);
subplot(3,3,8); cfg.channel = 'Oz';  ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - Oz']);
subplot(3,3,9); cfg.channel = 'C3';  ft_singleplotTFR(cfg, timefreq);  title(['valid-vs-invalid', ' - C3']);
colormap(colormap2use);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_compare_validvsinvalid_s23', toadd1,toadd2,toadd3], '.jpg']);

end


%% singleplot--plot single channel
if plotsingleChan

figure;
set(gcf,'Position',[0 0 1800 300])

cfg             = [];
cfg.zlim        = [-10 10];
% cfg.xlim        = xlim2plot;
cfg.xlim        = [-400 2500];
cfg.figure      = 'gcf';
cfg.baseline    = [-0.4 0];
cfg.baselinetype= 'db';

subplot(1,6,1); cfg.parameter = 'validim100'; cfg.channel = {'AFz';'Fz';'F1';'F2';'F3';'F4'}; ft_singleplotTFR(cfg, timefreq);  title(['validim100', ' - Frontal']);
subplot(1,6,2); cfg.parameter = 'validinfor100'; cfg.channel = {'AFz';'Fz';'F1';'F2';'F3';'F4'}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor100', ' - Frontal']);
subplot(1,6,3); cfg.parameter = 'validinfor80'; cfg.channel = {'AFz';'Fz';'F1';'F2';'F3';'F4'}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor80', ' - Frontal']);
subplot(1,6,4); cfg.parameter = 'validinfor60'; cfg.channel = {'AFz';'Fz';'F1';'F2';'F3';'F4'}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor60', ' - Frontal']);
subplot(1,6,5); cfg.parameter = 'validinfor40'; cfg.channel = {'AFz';'Fz';'F1';'F2';'F3';'F4'}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor40', ' - Frontal']);
subplot(1,6,6); cfg.parameter = 'validinfor20'; cfg.channel = {'AFz';'Fz';'F1';'F2';'F3';'F4'}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor20', ' - Frontal']);
colormap(colormap2use);


figure;
set(gcf,'Position',[0 0 1800 300])

cfg             = [];
cfg.zlim        = [-10 10];
% cfg.xlim        = xlim2plot;
cfg.xlim        = [-400 2500];
cfg.figure      = 'gcf';
cfg.baseline    = [-0.4 0];
cfg.baselinetype= 'db';

subplot(1,6,1); cfg.parameter = 'validim100'; cfg.channel = {'Cz';'C1';'C2';}; ft_singleplotTFR(cfg, timefreq);  title(['validim100', ' - Central']);
subplot(1,6,2); cfg.parameter = 'validinfor100'; cfg.channel = {'Cz';'C1';'C2';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor100', ' - Central']);
subplot(1,6,3); cfg.parameter = 'validinfor80'; cfg.channel = {'Cz';'C1';'C2';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor80', ' - Central']);
subplot(1,6,4); cfg.parameter = 'validinfor60'; cfg.channel = {'Cz';'C1';'C2';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor60', ' - Central']);
subplot(1,6,5); cfg.parameter = 'validinfor40'; cfg.channel = {'Cz';'C1';'C2';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor40', ' - Central']);
subplot(1,6,6); cfg.parameter = 'validinfor20'; cfg.channel = {'Cz';'C1';'C2';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor20', ' - Central']);
colormap(colormap2use);


figure;
set(gcf,'Position',[0 0 1800 300])

cfg             = [];
cfg.zlim        = [-10 10];
% cfg.xlim        = xlim2plot;
cfg.xlim        = [-400 2500];
cfg.figure      = 'gcf';
cfg.baseline    = [-0.4 0];
cfg.baselinetype= 'db';

subplot(1,6,1); cfg.parameter = 'validim100'; cfg.channel = {'Pz';'P1';'P2';'P3';'P4';'POz';}; ft_singleplotTFR(cfg, timefreq);  title(['validim100', ' - Parietal']);
subplot(1,6,2); cfg.parameter = 'validinfor100'; cfg.channel = {'Pz';'P1';'P2';'P3';'P4';'POz';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor100', ' - Parietal']);
subplot(1,6,3); cfg.parameter = 'validinfor80'; cfg.channel = {'Pz';'P1';'P2';'P3';'P4';'POz';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor80', ' - Parietal']);
subplot(1,6,4); cfg.parameter = 'validinfor60'; cfg.channel = {'Pz';'P1';'P2';'P3';'P4';'POz';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor60', ' - Parietal']);
subplot(1,6,5); cfg.parameter = 'validinfor40'; cfg.channel = {'Pz';'P1';'P2';'P3';'P4';'POz';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor40', ' - Parietal']);
subplot(1,6,6); cfg.parameter = 'validinfor20'; cfg.channel = {'Pz';'P1';'P2';'P3';'P4';'POz';}; ft_singleplotTFR(cfg, timefreq);  title(['validinfor20', ' - Parietal']);
colormap(colormap2use);

end


%% running topos
if plotRunningTopos
cfg             = [];
cfg.layout      = 'easycapM1.mat';
cfg.comment     = 'no';
cfg.style       = 'straight';
cfg.marker      = 'off';
cfg.zlim        = [-20 20];
cfg.xlim        = xlim2plot;
cfg.figure      = 'gcf';

stepsize = 250; times1 = [-750:stepsize:2000]; n = length(times1);

cfg.parameter   = 'imper_vs_high'; 
figure; 
for t = 1:n-1
    cfg.xlim = times1(t:t+1);
    subplot(3,n-1,t);          cfg.ylim = [4 7];   ft_topoplotER(cfg, timefreq); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
    subplot(3,n-1,t+n-1);      cfg.ylim = [8 12];  ft_topoplotER(cfg, timefreq);
    subplot(3,n-1,t+(n-1)*2);  cfg.ylim = [13 30]; ft_topoplotER(cfg, timefreq);
end
colormap(colormap2use);

cfg.parameter   = 'high_vs_low'; 
figure; 
for t = 1:n-1
    cfg.xlim = times1(t:t+1);
    subplot(3,n-1,t);          cfg.ylim = [4 7];   ft_topoplotER(cfg, timefreq); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
    subplot(3,n-1,t+n-1);      cfg.ylim = [8 12];  ft_topoplotER(cfg, timefreq);
    subplot(3,n-1,t+(n-1)*2);  cfg.ylim = [13 30]; ft_topoplotER(cfg, timefreq);
end
colormap(colormap2use);
end

%% cvsi--for all conditions & all freqs
% time-locked to cue onset
cfg             = [];
cfg.colorbar    = 'no';
cfg.zlim        = [-7 7];
cfg.xlim        = xlim2plot;
cfg.parameter   = 'data';
cfg.figure      = 'gcf';

% compare: 4 cue validity conditions
figure('Name','CueValidity_v1_tfa_spectrum_probe_cuevalid','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 300])
for channel = 1:4
    cfg.channel     = ci.label(channel);
    subplot(1,4,channel);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
end
ylabel('Frequency'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_cuevalid_s23', toadd1,toadd2,toadd3], '.jpg']);


% time-locked to probe onset
cfg             = [];
cfg.colorbar    = 'no';
cfg.zlim        = [-7 7];
cfg.xlim        = [1300,2500];
cfg.parameter   = 'data';
cfg.figure      = 'gcf';

% compare: valid vs invalid
figure('Name','CueValidity_v1_tfa_spectrum_probe_vainva','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 300])
for channel = 5:6
    cfg.channel     = ci.label(channel);
    subplot(1,2,channel-4);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
end
ylabel('Frequency'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_vainva_s23', toadd1,toadd2,toadd3], '.jpg']);

% compare: 6 cue validity conditions
figure('Name','CueValidity_v1_tfa_spectrum_probe_6con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 300])
for channel = 7:12
    cfg.channel     = ci.label(channel);
    subplot(1,6,channel-6);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
end
ylabel('Frequency'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_6con_s23', toadd1,toadd2,toadd3], '.jpg']);


%% timecourses--lateralized alpha (difference)--Time-locked to cue onset
fsel = ci.freq >= 8 & ci.freq <= 12;

figure('Name','CueValidity_v1_later_alpha_probe_cue_valid4con_overlay','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 1800 900])
hold on; title('alpha - contra vs ipsi');
[m_plot1] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,1,fsel,:),3)), [0, 0, 1], 'se');
[m_plot2] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,2,fsel,:),3)), [1, 0, 0], 'se');
[m_plot3] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,3,fsel,:),3)), [1, 0, 1], 'se');
[m_plot4] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,4,fsel,:),3)), [0, 1, 1], 'se');
plot(xlim, [0,0], ':k');
legend([m_plot1,m_plot2,m_plot3,m_plot4] , ci.label([1:4]));
xlim(xlim2plot);
set(gca, 'YDir','reverse') 
ylabel('Percentage Change'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
% ylim([1 3]);
% saveas(gcf, 'CueValidity_v1_alpha_overlay_s23', 'jpg')

saveas(gcf, [['CueValidity_v1_later_alpha_probe_cue_valid4con_overlay_s23', toadd1,toadd2,toadd3], '.jpg']);



%% timecourses--lateralized alpha (difference)--valid 6cons
xlim2plot_probe = [1000 2500];
fsel = ci.freq >= 8 & ci.freq <= 12;

figure('Name','CueValidity_v1_later_alpha_probe_probe_valid6con_overlay','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 1800 900])
hold on; title('alpha - contra vs ipsi');
[m_plot1] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,7,fsel,:),3)), [0, 0, 1], 'se');
[m_plot2] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,8,fsel,:),3)), [1, 0, 0], 'se');
[m_plot3] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,9,fsel,:),3)), [0, 1, 0], 'se');
[m_plot4] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,10,fsel,:),3)), [1, 1, 0], 'se');
[m_plot5] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,11,fsel,:),3)), [0, 1, 1], 'se');
[m_plot6] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,12,fsel,:),3)), [1, 0, 1], 'se');
plot(xlim, [0,0], ':k');
legend([m_plot1,m_plot2,m_plot3,m_plot4,m_plot5,m_plot6] , ci.label([7:12]));
legend("boxoff");
xlim(xlim2plot_probe);
set(gca, 'YDir','reverse') 
ylabel('Percentage Change'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
% ylim([1 3]);
% saveas(gcf, 'CueValidity_v1_alpha_overlay_s23', 'jpg')

saveas(gcf, [['CueValidity_v1_later_alpha_probe_probe_valid6con_overlay_s23', toadd1,toadd2,toadd3], '.jpg']);


figure('Name','CueValidity_v1_later_alpha_probe_probe_valid6con_singletrial','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 1800 900])
hold on; title('alpha - contra vs ipsi');
[m_plot1] = plot(ci.time, squeeze(mean(d3(:,7,fsel,:),3)), 'b');
[m_plot2] = plot(ci.time, squeeze(mean(d3(:,8,fsel,:),3)), 'r');
[m_plot3] = plot(ci.time, squeeze(mean(d3(:,9,fsel,:),3)), 'g');
[m_plot4] = plot(ci.time, squeeze(mean(d3(:,10,fsel,:),3)), 'y');
[m_plot5] = plot(ci.time, squeeze(mean(d3(:,11,fsel,:),3)), 'c');
[m_plot6] = plot(ci.time, squeeze(mean(d3(:,12,fsel,:),3)), 'm');
xlim(xlim2plot_probe);
plot(xlim, [0,0], ':k');
set(gca, 'YDir','reverse') 
saveas(gcf, [['CueValidity_v1_later_alpha_probe_probe_valid6con_singletrial_s23', toadd1,toadd2,toadd3], '.jpg']);


%% timecourses--lateralized alpha (difference)--valid 6cons
xlim2plot_probe = xlim2plot;
fsel = ci.freq >= 8 & ci.freq <= 12;

figure('Name','CueValidity_v1_later_alpha_probe_cue_valid6con_overlay','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 1800 900])
hold on; title('alpha - contra vs ipsi');
[m_plot1] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,7,fsel,:),3)), [0, 0, 1], 'se');
[m_plot2] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,8,fsel,:),3)), [1, 0, 0], 'se');
[m_plot3] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,9,fsel,:),3)), [0, 1, 0], 'se');
[m_plot4] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,10,fsel,:),3)), [1, 1, 0], 'se');
[m_plot5] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,11,fsel,:),3)), [0, 1, 1], 'se');
[m_plot6] = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,12,fsel,:),3)), [1, 0, 1], 'se');
plot(xlim, [0,0], ':k');
legend([m_plot1,m_plot2,m_plot3,m_plot4,m_plot5,m_plot6] , ci.label([7:12]));
xlim(xlim2plot_probe);
set(gca, 'YDir','reverse') 
ylabel('Percentage Change'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
% ylim([1 3]);
% saveas(gcf, 'CueValidity_v1_alpha_overlay_s23', 'jpg')

saveas(gcf, [['CueValidity_v1_later_alpha_probe_cue_valid6con_overlay_s23', toadd1,toadd2,toadd3], '.jpg']);


figure('Name','CueValidity_v1_later_alpha_probe_cue_valid6con_singletrial','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 1800 900])
hold on; title('alpha - contra vs ipsi');
[m_plot1] = plot(ci.time, squeeze(mean(d3(:,7,fsel,:),3)), 'b');
[m_plot2] = plot(ci.time, squeeze(mean(d3(:,8,fsel,:),3)), 'r');
[m_plot3] = plot(ci.time, squeeze(mean(d3(:,9,fsel,:),3)), 'g');
[m_plot4] = plot(ci.time, squeeze(mean(d3(:,10,fsel,:),3)), 'y');
[m_plot5] = plot(ci.time, squeeze(mean(d3(:,11,fsel,:),3)), 'c');
[m_plot6] = plot(ci.time, squeeze(mean(d3(:,12,fsel,:),3)), 'm');
xlim(xlim2plot_probe);
plot(xlim, [0,0], ':k');
set(gca, 'YDir','reverse') 
saveas(gcf, [['CueValidity_v1_later_alpha_probe_cue_valid6con_singletrial_s23', toadd1,toadd2,toadd3], '.jpg']);


%% Non-lateralized tfa plot
plotsingleChan_nonlater = 0;
if plotsingleChan_nonlater

for selch_area = [1 2 3];
% for selch_area = [1];

if selch_area == 1
selch_lb = {'AFz','Fz','F1','F2','F3','F4'};
toadd4 = '_front';
elseif selch_area == 2
selch_lb = {'Cz','C1','C2'};
toadd4 = '_cent';
elseif selch_area == 3
selch_lb = {'Pz','P1','P2','POz'};
toadd4 = '_post';
end
selch = match_str(timefreq.label,selch_lb);

ave = [];
ave.time = timefreq.time;
ave.freq = timefreq.freq;
ave.dimord = 'chan_freq_time';
ave.label = {'imper100-ave','high100-ave','med80-ave','low60-ave',...
            'valid-ave','invalid-ave',...
            'validim100-ave','validinfor100-ave','validinfor80-ave','validinfor60-ave','validinfor40-ave','validinfor20-ave',...
             'imper-vs-high','high-vs-low','valid-vs-invalid'};

ave.data(1,:,:) = squeeze(mean(timefreq.imper100(selch,:,:)));
ave.data(2,:,:) = squeeze(mean(timefreq.high100(selch,:,:)));
ave.data(3,:,:) = squeeze(mean(timefreq.med80(selch,:,:)));
ave.data(4,:,:) = squeeze(mean(timefreq.low60(selch,:,:)));
ave.data(5,:,:) = squeeze(mean(timefreq.valid(selch,:,:)));
ave.data(6,:,:) = squeeze(mean(timefreq.invalid(selch,:,:)));
ave.data(7,:,:) = squeeze(mean(timefreq.validim100(selch,:,:)));
ave.data(8,:,:) = squeeze(mean(timefreq.validinfor100(selch,:,:)));
ave.data(9,:,:) = squeeze(mean(timefreq.validinfor80(selch,:,:)));
ave.data(10,:,:) = squeeze(mean(timefreq.validinfor60(selch,:,:)));
ave.data(11,:,:) = squeeze(mean(timefreq.validinfor40(selch,:,:)));
ave.data(12,:,:) = squeeze(mean(timefreq.validinfor20(selch,:,:)));
ave.data(13,:,:) = squeeze(mean(timefreq.imper_vs_high(selch,:,:)));
ave.data(14,:,:) = squeeze(mean(timefreq.high_vs_low(selch,:,:)));
ave.data(15,:,:) = squeeze(mean(timefreq.valid_vs_invalid(selch,:,:)));


% % plot averaged tfa--for all conditions & all freqs
% time-locked to cue onset
cfg             = [];
cfg.colorbar    = 'no';
% cfg.zlim        = [-0.000001 0.000001];
cfg.xlim        = xlim2plot;
cfg.parameter   = 'data';
cfg.figure      = 'gcf';
cfg.baseline    = [-0.4 0];
% cfg.baseline    = [1.1 1.5];
cfg.baselinetype= 'db';

% compare: 4 cue validity conditions
figure('Name','CueValidity_v1_tfa_spectrum_probe_ave_cuevalid','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 300])
for channel = 1:4
    cfg.channel     = ave.label(channel);
    subplot(1,4,channel);    ft_singleplotTFR(cfg, ave);
    colormap(colormap2use);
end
ylabel('Frequency'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_ave_cuevalid_s23', toadd1,toadd2,toadd3,toadd4], '.jpg']);


% time-locked to probe onset
cfg             = [];
cfg.colorbar    = 'no';
% cfg.zlim        = [-0.000001 0.000001];
cfg.xlim        = [1300,2500];
cfg.parameter   = 'data';
cfg.figure      = 'gcf';
cfg.baseline    = [-0.4 0];
% cfg.baseline    = [1.1 1.5];
cfg.baselinetype= 'db';

% compare: valid vs invalid
figure('Name','CueValidity_v1_tfa_spectrum_probe_ave_vainva','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 300])
for channel = 5:6
    cfg.channel     = ave.label(channel);
    subplot(1,2,channel-4);    ft_singleplotTFR(cfg, ave);
    colormap(colormap2use);
end
ylabel('Frequency'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_ave_vainva_s23', toadd1,toadd2,toadd3,toadd4], '.jpg']);

% compare: 6 cue validity conditions
figure('Name','CueValidity_v1_tfa_spectrum_probe_6con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1800 300])
for channel = 7:12
    cfg.channel     = ave.label(channel);
    subplot(1,6,channel-6);    ft_singleplotTFR(cfg, ave);
    colormap(colormap2use);
end
ylabel('Frequency'), xlabel('Time')
set(gca, 'FontSize', 10)
set(gca,'LineWidth',2);
saveas(gcf, [['CueValidity_v1_tfa_spectrum_probe_ave_6con_s23', toadd1,toadd2,toadd3,toadd4], '.jpg']);

end
end


end
end
end

% close all
