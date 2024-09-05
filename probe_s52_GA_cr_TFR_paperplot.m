%% Plot for papers

clc;clear
close all


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

    load([param.pathnewlen, 'saved_data/timefreq', toadd1,toadd2,toadd3, '__' param.subjName], 'timefreq','ci','cic','cii','ci_allchan');


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

    d17(s,:,:,:) = cic.data;
    d18(s,:,:,:) = cii.data;

end


%% put back into structure
% calculate valid-invalid difference for 80 & 60, respectively
for ss = 1:s
    dvainva_80(s,:,:,:)   = ((d8(s,:,:,:) - d11(s,:,:,:)) ./ (d8(s,:,:,:) + d11(s,:,:,:))) * 100;
    dvainva_60(s,:,:,:)   = ((d9(s,:,:,:) - d10(s,:,:,:)) ./ (d9(s,:,:,:) + d10(s,:,:,:))) * 100;
end

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
cic.data                = squeeze(mean(d17));
cii.data                = squeeze(mean(d18));

timefreq.vainva_80  = squeeze(mean(dvainva_80));
timefreq.vainva_60  = squeeze(mean(dvainva_60));

save(fullfile('recue_valid_results_timefreq.mat'))


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
statcfg.clusterStatEvalaluationAlpha = 0.01; % two-sided
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

save(fullfile('recue_valid_results_timefreq_permtest_cuelocked_alpha.mat'))


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

save(fullfile('recue_valid_results_timefreq_permtest_cuelocked_thetabeta.mat'))


%% later alpha cluster test--probelocked 6con
fsel = ci.freq >= 8 & ci.freq <= 12;
% pnt_start = 266; % 1300ms
pnt_start = 286; % 1500ms
pnt_end = 370; % 2400ms

data_cond111 = []; data_cond112 = []; 
data_cond111 = squeeze(mean(d3(:,7,fsel,pnt_start:pnt_end),3)); % im100%
data_cond112 = squeeze(mean(d3(:,8,fsel,pnt_start:pnt_end),3)); % in100%
data_cond1=[]; data_cond2=[]; data_cond3=[]; data_cond4=[]; 
data_cond1 = squeeze(mean(d3(:,9,fsel,pnt_start:pnt_end),3)); % 80-valid
data_cond2 = squeeze(mean(d3(:,12,fsel,pnt_start:pnt_end),3)); % 80-invalid
data_cond3 = squeeze(mean(d3(:,10,fsel,pnt_start:pnt_end),3)); % 60-valid
data_cond4 = squeeze(mean(d3(:,11,fsel,pnt_start:pnt_end),3)); % 60-invalid
data_zero = zeros(size(data_cond1)); % compare with 0

statcfg.xax = ci.time(pnt_start:pnt_end);
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.025; % two-sided
statcfg.nsub = s;
statcfg.statMethod = 'montecarlo'; % or 'analytic'; 

% cluster test
% compare with zeros
stat111=[]; stat112=[]; stat1=[]; stat2=[]; stat3=[]; stat4=[];
stat111 = frevede_ftclusterstat1D(statcfg, data_cond111, data_zero);
stat112 = frevede_ftclusterstat1D(statcfg, data_cond112, data_zero);
stat1 = frevede_ftclusterstat1D(statcfg, data_cond1, data_zero);
stat2 = frevede_ftclusterstat1D(statcfg, data_cond2, data_zero);
stat3 = frevede_ftclusterstat1D(statcfg, data_cond3, data_zero);
stat4 = frevede_ftclusterstat1D(statcfg, data_cond4, data_zero);

% compare conditions
stat12=[];stat34=[];
stat12 = frevede_ftclusterstat1D(statcfg, data_cond1, data_cond2);
stat34 = frevede_ftclusterstat1D(statcfg, data_cond3, data_cond4);

sigposclurange12=[]; 
sigposclurange12 = statcfg.xax(find(stat12.posclusterslabelmat==1));
signegclurange1=[]; signegclurange2=[]; 
signegclurange1 = statcfg.xax(find(stat1.negclusterslabelmat==1));
signegclurange2 = statcfg.xax(find(stat2.negclusterslabelmat==1));

save(fullfile('recue_valid_results_timefreq_permtest_probelocked_alpha.mat'))


%% plot ------------------------------------------------------------------------------------

%% Plot Cue-locked

%% Plot lateralized alpha time-series--Cuelocked
fsel = ci.freq >= 8 & ci.freq <= 12;
title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};

% contra vs. ipsi alpha
figure('Name','paperplot_CueValidity_v1_cuelocked_alpha_conipsi_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    subplot(1,4,sp); hold on; title(title2use{sp});
    [p1] = frevede_errorbarplot(ci.time, squeeze(mean(d17(:,sp,fsel,:),3)), [1,0.6,0], 'se');
    [p2] = frevede_errorbarplot(ci.time, squeeze(mean(d18(:,sp,fsel,:),3)), [0,1,1], 'se');
    xlim(xlim2plot); xline(0,'--k','Cue Onset');
    %   plot(xlim, [0,0], '--k');  ylim([-10 10]);
    set(gca, 'YDir','reverse'); ylabel('Power'), xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
legend([p1, p2], {'Contralateral','Ipsilateral'},'FontSize',8,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValidity_v1_cuelocked_alpha_conipsi_sepcon_s22', 'jpg')


%lateralized alpha: seperate conditions
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];
figure('Name','paperplot_CueValidity_v1_cuelocked_alpha_later_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    subplot(1,4,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp,fsel,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-12 5]);
    xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
% legend({'Lateralized Alpha'},'FontSize',8,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValidity_v1_cuelocked_alpha_later_sepcon_s22', 'jpg')


%later alpha: overlay
% col2use = [0,0,1; 1,0,0; 1,0,1; 0,1,1;];
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];
figure('Name','paperplot_CueValidity_v1_cuelocked_alpha_later_overlay_4con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])
hold on;
% for sp = 1:4
p1 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,1,fsel,:),3)), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,2,fsel,:),3)), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,3,fsel,:),3)), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,4,fsel,:),3)), [col2use(4,:)], 'se');
% end
legend([p1,p2,p3,p4],title2use{1:4},'FontSize',8,'AutoUpdate','off','Box','off');
plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-12 5]);
xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
title('Lateralized Alpha');

saveas(gcf, 'paperplot_CueValidity_v1_cuelocked_alpha_later_overlay_4con_s22', 'jpg')


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
figure('Name','CueValidity_v1_cuelocked_alpha_later_ave_4con_bar400800','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])
hold on,
b = bar(1:4, mean(aveerp(:,:)));
b.LineWidth = 2;
b.FaceColor = 'flat';
for c = 1:4
    b.CData(c,:) = col2use(c,:);
end
hold on,
errorbar(1:4, mean(aveerp(:,:)), std(aveerp(:,:))./sqrt((size(aveerp,1))), '.k', 'LineWidth',2);
%     xticks(1:4); xticklabels({'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'});
xticks(1:4); xticklabels({'Imper100%','Infor100%','Infor80%','Infor60%'});
ylim([-10 2]);set(gca, 'YDir','reverse');
ylabel('Mean Percentage Change'); xlabel('Cue Type');
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);

saveas(gcf, [['paperplot_CueValidity_v1_cuelocked_alpha_later_ave_4con_bar400800_s23', toadd1,toadd2,toadd3], '.jpg']);
clear aveerp b


%% Plot Lateralized Theta & Alpha & Beta time-series--Cuelocked

title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];

fsel1 = ci.freq >= 3 & ci.freq <= 7;
fsel = ci.freq >= 8 & ci.freq <= 12;
fsel2 = ci.freq >= 13 & ci.freq <= 30;

% Lateralized Theta: 
% seperate conditions
figure('Name','paperplot_CueValidity_v1_cuelocked_thalbe_later_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    subplot(1,4,sp); hold on; title(title2use{sp});
    p1 = frevede_errorbarplot_dashline(ci.time, squeeze(mean(d3(:,sp,fsel1,:),3)), [col2use(sp,:)], 'se');
    p2 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp,fsel,:),3)), [col2use(sp,:)], 'se');
    p3 = frevede_errorbarplot_dotline(ci.time, squeeze(mean(d3(:,sp,fsel2,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-12 5]);
    xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
legend([p1,p2,p3],{'Lateralized Theta','Lateralized Alpha','Lateralized Beta'},'FontSize',8,'AutoUpdate','off','Box','off');
% legend({'Lateralized Oscillations'},'FontSize',8,'AutoUpdate','off','Box','off');
% saveas(gcf, 'paperplot_CueValidity_v1_cuelocked_thalbe_later_sepcon_s22', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValidity_v1_cuelocked_thalbe_later_sepcon_s22', 'epsc')



%% Plot Averaged Center Theta time-series--Cuelocked

title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];

fsel1 = ci.freq >= 3 & ci.freq <= 7;
% choose electrode
e2plot = match_str(timefreq.label, 'Fz'); 

% % Averaged Theta--no BC: 
% % seperate conditions
% figure('Name','paperplot_CueValidity_v1_cuelocked_centertheta_overlay','NumberTitle','off', 'Color','white'),
% set(gcf,'Position',[0 0 400 300])
% 
% hold on; title('Center Theta')
% p1 = frevede_errorbarplot(ci.time, squeeze(mean(d12(:,e2plot,fsel1,:),3)), [col2use(1,:)], 'se');
% p2 = frevede_errorbarplot(ci.time, squeeze(mean(d13(:,e2plot,fsel1,:),3)), [col2use(2,:)], 'se');
% p3 = frevede_errorbarplot(ci.time, squeeze(mean(d14(:,e2plot,fsel1,:),3)), [col2use(3,:)], 'se');
% p4 = frevede_errorbarplot(ci.time, squeeze(mean(d15(:,e2plot,fsel1,:),3)), [col2use(4,:)], 'se');
% plot(xlim, [0,0], '--k'); xlim(xlim2plot); 
% xline(0,'--k','Cue Onset');xlabel('Time (ms)');set(gca, 'YDir','reverse');
% % ylim([-12 5]);ylabel('Percentage Change'); 
% set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
% 
% legend([p1,p2,p3,p4],title2use{1:4},'FontSize',8,'AutoUpdate','off','Box','off');
% saveas(gcf, 'paperplot_CueValidity_v1_cuelocked_centertheta_overlay_s22', 'jpg')


tpnt_s = 116; % -200
tpnt_e = 286; % 1500
tpnt_0 = 136; % 0

% baseline correction:
for ss = 1:size(d12,1)
     theta = squeeze(mean(d12(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d12(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,1,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
     theta = squeeze(mean(d13(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d13(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,2,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
     theta = squeeze(mean(d14(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d14(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,3,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
     theta = squeeze(mean(d15(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d15(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,4,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
end


% Averaged Theta--BC: 
% seperate conditions
figure('Name','paperplot_CueValidity_v1_cuelocked_centertheta_overlay','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])

hold on; title('Center Theta')
p1 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,1,:)), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,2,:)), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,3,:)), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,4,:)), [col2use(4,:)], 'se');
plot(xlim, [0,0], '--k'); xlim([-200 1000]); 
xline(0,'--k','Cue Onset');xlabel('Time (ms)');set(gca, 'YDir','reverse');
% ylim([-12 5]);ylabel('Percentage Change'); 
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);

legend([p1,p2,p3,p4],title2use{1:4},'FontSize',8,'AutoUpdate','off','Box','off');


saveas(gcf, 'paperplot_CueValidity_v1_cuelocked_centertheta_overlay_s22', 'jpg')



%% plot topos--cuelocked
plotRunningTopos = 1;

if plotRunningTopos
    cfg             = [];
    cfg.layout      = 'easycapM1.mat';
    cfg.comment     = 'no';
    cfg.style       = 'straight';
    cfg.marker      = 'off';
%     cfg.zlim        = [-10 10];
    cfg.xlim        = xlim2plot;
    cfg.figure      = 'gcf';

    stepsize = 250; times1 = [0:stepsize:2500]; n = length(times1);

    cfg.parameter   = 'low60';
    figure(1);
    set(gcf,'Position',[0 0 1600 600])
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(3,n-1,t);          cfg.ylim = [4 7];   ft_topoplotER(cfg, timefreq); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
        subplot(3,n-1,t+n-1);      cfg.ylim = [8 12];  ft_topoplotER(cfg, timefreq);
        subplot(3,n-1,t+(n-1)*2);  cfg.ylim = [13 30]; ft_topoplotER(cfg, timefreq);
    end
    colormap(colormap2use);
    colorbar

    cfg.parameter   = 'vainva_80';
    figure(2);
    set(gcf,'Position',[0 0 1600 600])
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(3,n-1,t);          cfg.ylim = [4 7];   ft_topoplotER(cfg, timefreq); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
        subplot(3,n-1,t+n-1);      cfg.ylim = [8 12];  ft_topoplotER(cfg, timefreq);
        subplot(3,n-1,t+(n-1)*2);  cfg.ylim = [13 30]; ft_topoplotER(cfg, timefreq);
    end
    colormap(colormap2use);
    colorbar

    cfg.parameter   = 'vainva_60';
    figure(3);
    set(gcf,'Position',[0 0 1600 600])
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(3,n-1,t);          cfg.ylim = [4 7];   ft_topoplotER(cfg, timefreq); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
        subplot(3,n-1,t+n-1);      cfg.ylim = [8 12];  ft_topoplotER(cfg, timefreq);
        subplot(3,n-1,t+(n-1)*2);  cfg.ylim = [13 30]; ft_topoplotER(cfg, timefreq);
    end
    colormap(colormap2use);
    colorbar

end


%% Probelocked

%% Plot alpha--Probelocked
% contra vs. ipsi alpha
fsel = ci.freq >= 8 & ci.freq <= 12;
title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%','Informative 40%','Informative 20%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];

figure('Name','paperplot_CueValidity_v1_probelocked_alpha_conipsi_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 250])
for sp = 1:6
    subplot(1,6,sp); hold on; title(title2use{sp});
    [p1] = frevede_errorbarplot(ci.time, squeeze(mean(d17(:,sp+6,fsel,:),3)), [1,0.6,0], 'se');
    [p2] = frevede_errorbarplot(ci.time, squeeze(mean(d18(:,sp+6,fsel,:),3)), [0,1,1], 'se');
    xlim(xlim2plot_probe); xline(1500,'--k','Probe Onset');
    %   plot(xlim, [0,0], '--k');  ylim([-10 10]);
    set(gca, 'YDir','reverse'); ylabel('Power'), xlabel('Time (ms)');
    set(gca, 'FontSize', 10); set(gca,'LineWidth',2);
end
legend([p1, p2], {'Contralateral','Ipsilateral'},'FontSize',8,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_conipsi_sepcon_s22', 'jpg')


%lateralized alpha: seperate conditions
figure('Name','paperplot_CueValidity_v1_probelocked_alpha_later_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 250])
for sp = 1:6
    subplot(1,6,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp+6,fsel,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
    xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 10); set(gca,'LineWidth',2);
end
legend({'Lateralized Alpha'},'FontSize',8,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_sepcon_s22', 'jpg')


% overlay
% col2use = [0,0,1; 1,0,0; 1,0,1; 0,1,1;];
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];
figure('Name','paperplot_CueValidity_v1_probelocked_alpha_later_overlay_6con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 400])
hold on;
% for sp = 1:6
p1 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,7,fsel,:),3)), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,8,fsel,:),3)), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,9,fsel,:),3)), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,10,fsel,:),3)), [col2use(4,:)], 'se');
p5 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,11,fsel,:),3)), [col2use(5,:)], 'se');
p6 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,12,fsel,:),3)), [col2use(6,:)], 'se');
% end
legend([p1,p2,p3,p4,p5,p6],title2use{1:6},'FontSize',6,'AutoUpdate','off','Box','off');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
set(gca, 'FontSize', 10); set(gca,'LineWidth',2);
title('Lateralized Alpha');

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_overlay_6con_s22', 'jpg')



%% Plot alpha--Probelocked--new sequence
% contra vs. ipsi alpha
fsel = ci.freq >= 8 & ci.freq <= 12;
title2use = {'Imperative 100%','Informative 100%','Informative 80%-Valid','Informative 80%-Invalid','Informative 60%-Valid','Informative 60%-Invalid'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1;0.72 0.27 1; 0.9,0.5,0.1;0.56 0.2 0.02;];


% alpha--contra & ipsilateral
figure('Name','paperplot_CueValidity_v1_probelocked_alpha_conipsi_sepcon_newseq','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 200])
for sp = 1:6
    if sp==1 con=7; elseif sp==2 con=8; elseif sp==3 con=9; elseif sp==4 con=12;elseif sp==5 con=10;elseif sp==6 con=11; end
    subplot(1,6,sp); hold on; title(title2use{sp});
    [p1] = frevede_errorbarplot(ci.time, squeeze(mean(d17(:,con,fsel,:),3)), [0,0,0], 'se');
    [p2] = frevede_errorbarplot(ci.time, squeeze(mean(d18(:,con,fsel,:),3)), [0.7,0.7,0.7], 'se');
    xlim(xlim2plot_probe); xline(1500,'--k','Test Onset');
    %   plot(xlim, [0,0], '--k');  ylim([-10 10]);
    set(gca, 'YDir','reverse'); ylabel('Power'), xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
legend([p1, p2], {'Contralateral','Ipsilateral'},'FontSize',8,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_conipsi_sepcon_newseq_s22', 'jpg')


%lateralized alpha: seperate conditions
figure('Name','paperplot_CueValidity_v1_probelocked_alpha_later_sepcon_newseq','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 200])
for sp = 1:6
    if sp==1 con=7; elseif sp==2 con=8; elseif sp==3 con=9; elseif sp==4 con=12;elseif sp==5 con=10;elseif sp==6 con=11; end
    subplot(1,6,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, squeeze(mean(d3(:,con,fsel,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
    xline(1500,'--k','Test Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
legend({'Lateralized Alpha'},'FontSize',8,'AutoUpdate','off','Box','off');
saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_sepcon_newseq_s22', 'jpg')


% overlay
% col2use = [0,0,1; 1,0,0; 1,0,1; 0,1,1;];
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];
figure('Name','paperplot_CueValidity_v1_probelocked_alpha_later_overlay_6con_newseq','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])
hold on;
% for sp = 1:6
p1 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,7,fsel,:),3)), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,8,fsel,:),3)), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,9,fsel,:),3)), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,12,fsel,:),3)), [col2use(4,:)], 'se');
p5 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,10,fsel,:),3)), [col2use(5,:)], 'se');
p6 = frevede_errorbarplot(ci.time, squeeze(mean(d3(:,11,fsel,:),3)), [col2use(6,:)], 'se');
% end
legend([p1,p2,p3,p4,p5,p6],title2use{1:6},'FontSize',6,'AutoUpdate','off','Box','off');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Test Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
title('Lateralized Alpha');

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_overlay_6con_newseq_s22', 'jpg')



%% Plot alpha--Probelocked--new layout--4plots

% contra vs. ipsi alpha
fsel = ci.freq >= 8 & ci.freq <= 12;
title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1;0.72 0.27 1; 0.9,0.5,0.1;0.56 0.2 0.02;];

%lateralized alpha: seperate conditions
figure('Name','paperplot_CueValidity_v1_probelocked_alpha_later_sepcon_4plots','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:2
    subplot(1,4,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp+6,fsel,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
    xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    legend({'Valid'},'FontSize',8,'AutoUpdate','off','Box','off');
end
hold on;

subplot(1,4,3); hold on; title(title2use{3});
p1=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,9,fsel,:),3)), [col2use(3,:)], 'se');
p2=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,12,fsel,:),3)), [col2use(4,:)], 'se');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

subplot(1,4,4); hold on; title(title2use{4});
p1=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,10,fsel,:),3)), [col2use(5,:)], 'se');
p2=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,11,fsel,:),3)), [col2use(6,:)], 'se');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_sepcon_4plots_s22', 'jpg')


%% Probelocked alpha--bar plot--averaged lateralized alpha
fsel = ci.freq >= 8 & ci.freq <= 12;
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1;0.72 0.27 1; 0.9,0.5,0.1;0.56 0.2 0.02;];

timeproon = 1500;
timerange = [400+timeproon 800+timeproon];
%     timerange = [300+timeproon 800+timeproon];

timediff1 = abs(ci.time-timerange(1));
pnt_start = find(timediff1==min(timediff1));
timediff2 = abs(ci.time-timerange(2));
pnt_end = find(timediff2==min(timediff2));
clear timediff1 timediff2

% average TFA across time
aveerp = [];
aveerp(:,1) = mean(squeeze(mean(d3(:,7,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,2) = mean(squeeze(mean(d3(:,8,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,3) = mean(squeeze(mean(d3(:,9,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,4) = mean(squeeze(mean(d3(:,12,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,5) = mean(squeeze(mean(d3(:,10,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,6) = mean(squeeze(mean(d3(:,11,fsel,pnt_start:pnt_end),4)),2);

% bar plot
figure('Name','CueValidity_v1_probelocked_alpha_later_ave_6con_bar400800','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 280 400])
hold on,
b = bar(1:6, mean(aveerp(:,:)));
b.LineWidth = 2;
b.FaceColor = 'flat';
for c = 1:6
    b.CData(c,:) = col2use(c,:);
end
hold on,
errorbar(1:6, mean(aveerp(:,:)), std(aveerp(:,:))./sqrt((size(aveerp,1))), '.k', 'LineWidth',2);
%     xticks(1:6); xticklabels({'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'});
xticks(1:6); xticklabels({'Imp100%','100%','80%-valid','80%-invalid','60%-valid','60%-invalid'});
ylim([-6 2]);set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Cue Type');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3);

% saveas(gcf, [['paperplot_CueValidity_v1_probelocked_alpha_later_ave_6con_bar400800_s23', toadd1,toadd2,toadd3], '.jpg']);
% clear aveerp b

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_ave_6con_bar400800_s23', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_ave_6con_bar400800_s23', 'epsc')


%% Probelocked alpha--bar plot--averaged lateralized alpha--decending sequence
fsel = ci.freq >= 8 & ci.freq <= 12;
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1;0.72 0.27 1; 0.9,0.5,0.1;0.56 0.2 0.02;];

timeproon = 1500;
timerange = [400+timeproon 800+timeproon];
%     timerange = [300+timeproon 800+timeproon];

timediff1 = abs(ci.time-timerange(1));
pnt_start = find(timediff1==min(timediff1));
timediff2 = abs(ci.time-timerange(2));
pnt_end = find(timediff2==min(timediff2));
clear timediff1 timediff2

% average TFA across time
aveerp = [];
aveerp(:,1) = mean(squeeze(mean(d3(:,7,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,2) = mean(squeeze(mean(d3(:,8,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,3) = mean(squeeze(mean(d3(:,9,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,4) = mean(squeeze(mean(d3(:,10,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,5) = mean(squeeze(mean(d3(:,11,fsel,pnt_start:pnt_end),4)),2);
aveerp(:,6) = mean(squeeze(mean(d3(:,12,fsel,pnt_start:pnt_end),4)),2);

% bar plot
figure('Name','CueValidity_v1_probelocked_alpha_later_ave_6con_bar400800','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 280 400])
hold on,
b = bar(1:6, mean(aveerp(:,:)));
b.LineWidth = 2;
b.FaceColor = 'flat';
for c = 1:6
    b.CData(c,:) = col2use(c,:);
end
hold on,
errorbar(1:6, mean(aveerp(:,:)), std(aveerp(:,:))./sqrt((size(aveerp,1))), '.k', 'LineWidth',2);
xticks(1:6); xticklabels({'Imp100%','100%','80%','60%','40%','20%'});
ylim([-6 2]);set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Cue Reliability');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3);

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_ave_6con_bar400800_newseq_s23', 'jpg')
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValidity_v1_probelocked_alpha_later_ave_6con_bar400800_newseq_s23', 'epsc')



%% Plot theta & beta--Probelocked--new layout--4plots
title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1;0.72 0.27 1; 0.9,0.5,0.1;0.56 0.2 0.02;];

fsel1 = ci.freq >= 4 & ci.freq <= 8;
fsel2 = ci.freq >= 13 & ci.freq <= 30;

%lateralized Theta
figure('Name','paperplot_CueValidity_v1_probelocked_theta_later_sepcon_4plots','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:2
    subplot(1,4,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp+6,fsel1,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
    xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    legend({'Valid'},'FontSize',8,'AutoUpdate','off','Box','off');
end
hold on;

subplot(1,4,3); hold on; title(title2use{3});
p1=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,9,fsel1,:),3)), [col2use(3,:)], 'se');
p2=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,12,fsel1,:),3)), [col2use(4,:)], 'se');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

subplot(1,4,4); hold on; title(title2use{4});
p1=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,10,fsel1,:),3)), [col2use(5,:)], 'se');
p2=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,11,fsel1,:),3)), [col2use(6,:)], 'se');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_theta_later_sepcon_4plots_s22', 'jpg')



%lateralized Beta
figure('Name','paperplot_CueValidity_v1_probelocked_beta_later_sepcon_4plots','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:2
    subplot(1,4,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, squeeze(mean(d3(:,sp+6,fsel2,:),3)), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
    xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
    ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    legend({'Valid'},'FontSize',8,'AutoUpdate','off','Box','off');
end
hold on;

subplot(1,4,3); hold on; title(title2use{3});
p1=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,9,fsel2,:),3)), [col2use(3,:)], 'se');
p2=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,12,fsel2,:),3)), [col2use(4,:)], 'se');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

subplot(1,4,4); hold on; title(title2use{4});
p1=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,10,fsel2,:),3)), [col2use(5,:)], 'se');
p2=frevede_errorbarplot(ci.time, squeeze(mean(d3(:,11,fsel2,:),3)), [col2use(6,:)], 'se');
plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-10 5]);
xline(1500,'--k','Probe Onset');set(gca, 'YDir','reverse');
ylabel('Percentage Change'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_beta_later_sepcon_4plots_s22', 'jpg')



%% Plot Averaged Center Theta & Posterior Alpha time-series--Probelocked

title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;];

tpnt_s = 266; % 1300
tpnt_e = 370; % 2300
tpnt_0 = 286; % 1500

% Center Theta
fsel1 = ci.freq >= 3 & ci.freq <= 7;
% choose electrode
e2plot = match_str(timefreq.label, 'Fz'); 

% baseline correction:
theta = []; theta_bc = [];
for ss = 1:size(d12,1)
     theta = squeeze(mean(d12(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d12(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,1,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
     theta = squeeze(mean(d13(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d13(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,2,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
     theta = squeeze(mean(d14(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d14(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,3,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
     theta = squeeze(mean(d15(ss,e2plot,fsel1,tpnt_s:tpnt_e),3));
     bl = mean(squeeze(mean(d15(ss,e2plot,fsel1,tpnt_s:tpnt_0),3)));
     theta_bc(ss,4,:) = ((theta-bl) ./(theta+bl))*100; 
     clear theta bl;
end

% Averaged Theta--BC: 
% seperate conditions
figure('Name','paperplot_CueValidity_v1_probelocked_centertheta_overlay','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])

hold on; title('Center Theta')
p1 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,1,:)), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,2,:)), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,3,:)), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(theta_bc(:,4,:)), [col2use(4,:)], 'se');
plot(xlim, [0,0], '--k'); xlim([1300 2400]); 
xline(0,'--k','Cue Onset');xlabel('Time (ms)');set(gca, 'YDir','reverse');
% ylim([-12 5]);ylabel('Percentage Change'); 
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2,p3,p4],title2use{1:4},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_centertheta_overlay_s22', 'jpg')


% Posterior Alpha
fsel = ci.freq >= 8 & ci.freq <= 12;
% choose electrode
e2plot = match_str(timefreq.label, {'O1','PO7','PO3','P9','P7','P5','P3','P1','O2','PO8','PO4','P10','P8','P6','P4','P2','Pz','POz','Oz'}); 

% baseline correction:
alpha = []; alpha_bc = [];
for ss = 1:size(d12,1)
     alpha = mean(squeeze(mean(d12(ss,e2plot,fsel,tpnt_s:tpnt_e),3)),1);
     bl = mean(mean(squeeze(mean(d12(ss,e2plot,fsel,tpnt_s:tpnt_0),3)),1));
     alpha_bc(ss,1,:) = ((alpha-bl) ./(alpha+bl))*100; 
     clear alpha bl;
     alpha = mean(squeeze(mean(d13(ss,e2plot,fsel,tpnt_s:tpnt_e),3)),1);
     bl = mean(mean(squeeze(mean(d13(ss,e2plot,fsel,tpnt_s:tpnt_0),3)),1));
     alpha_bc(ss,2,:) = ((alpha-bl) ./(alpha+bl))*100; 
     clear alpha bl;
     alpha = mean(squeeze(mean(d14(ss,e2plot,fsel,tpnt_s:tpnt_e),3)),1);
     bl = mean(mean(squeeze(mean(d14(ss,e2plot,fsel,tpnt_s:tpnt_0),3)),1));
     alpha_bc(ss,3,:) = ((alpha-bl) ./(alpha+bl))*100; 
     clear alpha bl;
     alpha = mean(squeeze(mean(d15(ss,e2plot,fsel,tpnt_s:tpnt_e),3)),1);
     bl = mean(mean(squeeze(mean(d15(ss,e2plot,fsel,tpnt_s:tpnt_0),3)),1));
     alpha_bc(ss,4,:) = ((alpha-bl) ./(alpha+bl))*100; 
     clear alpha bl;
end

% Averaged alpha--BC: 
% seperate conditions
figure('Name','paperplot_CueValidity_v1_probelocked_postalpha_overlay','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])

hold on; title('Posterior Alpha')
p1 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(alpha_bc(:,1,:)), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(alpha_bc(:,2,:)), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(alpha_bc(:,3,:)), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time(tpnt_s:tpnt_e), squeeze(alpha_bc(:,4,:)), [col2use(4,:)], 'se');
plot(xlim, [0,0], '--k'); xlim([1300 2400]); 
xline(0,'--k','Cue Onset');xlabel('Time (ms)');set(gca, 'YDir','reverse');
% ylim([-12 5]);ylabel('Percentage Change'); 
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
legend([p1,p2,p3,p4],title2use{1:4},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, 'paperplot_CueValidity_v1_probelocked_postalpha_overlay_s22', 'jpg')




%% plot spectrum

%% plot spectrum -- c vs i -- cue-locked
xlim2plot           = [-200 1500];
xlim2plot_probe     = [1300 2500];

cfg             = [];
cfg.colorbar    = 'no';
cfg.zlim        = [-7 7];
cfg.xlim        = xlim2plot;
cfg.parameter   = 'data';
cfg.figure      = 'gcf';

title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};

% compare: 4 cue validity conditions
figure('Name','paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 250])
for channel = 1:4
    cfg.channel     = ci.label(channel);
    subplot(1,4,channel);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
    xline(0,'--k','Cue Onset')
    title(title2use{channel});
    ylabel('Frequency (Hz)'); xlabel('Time (ms)');
    set(gca, 'FontSize', 10); set(gca,'LineWidth',2);
    %     colorbar;
end
% saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_s23', toadd1,toadd2,toadd3, '_colorbar'], '.jpg']);
saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_s23', toadd1,toadd2,toadd3], '.jpg']);


%% plot spectrum -- c vs i -- cue-locked--new
cfg             = [];
cfg.colorbar    = 'no';
cfg.zlim        = [-5 5];
cfg.xlim        = xlim2plot;
cfg.parameter   = 'data';
cfg.figure      = 'gcf';

title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};

% compare: 4 cue validity conditions
figure('Name','paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 400])
for channel = 1:4
    cfg.channel     = ci.label(channel);
    subplot(1,4,channel);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
    xline(0,'--k','Cue Onset')
    title(title2use{channel});
    ylabel('Frequency (Hz)'); xlabel('Time (ms)');
    set(gca, 'FontSize', 10); set(gca,'LineWidth',2);
    %     colorbar;
end
% saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_new_s23', toadd1,toadd2,toadd3, '_colorbar'], '.jpg']);
saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_new_s23', toadd1,toadd2,toadd3], '.jpg']);



%% plot spectrum -- c vs i -- cue-locked--v2_reorient
cfg             = [];
cfg.colorbar    = 'no';
cfg.zlim        = [-5 5];
cfg.xlim        = xlim2plot;
cfg.parameter   = 'data';
cfg.figure      = 'gcf';

title2use = {'Im100%','100%','80%','60%'};

% compare: 4 cue validity conditions
figure('Name','paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_v2_reorient','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for channel = 1:4
    cfg.channel     = ci.label(channel);
    subplot(1,4,channel);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
    xline(0,'--k','Cue Onset','LineWidth',2); 
%     yline(8,'--k'); yline(12,'--k');
    title(title2use{channel});
    ylabel('Frequency (Hz)'); xlabel('Time (ms)');
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);
end
colorbar;

saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_v2_reorient_s23', toadd1,toadd2,toadd3, '_colorbar'], '.jpg']);
set(gcf,'renderer','painters')
saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_v2_reorient_s23', toadd1,toadd2,toadd3, '_colorbar'], '.epsc']);

% saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_v2_reorient_s23', toadd1,toadd2,toadd3], '.jpg']);
% set(gcf,'renderer','painters')
% saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_cuelocked_4con_v2_reorient_s23', toadd1,toadd2,toadd3], '.epsc']);



%% plot spectrum -- c vs i -- probe-locked
cfg             = [];
cfg.colorbar    = 'no';
cfg.zlim        = [-7 7];
cfg.xlim        = xlim2plot_probe;
cfg.parameter   = 'data';
cfg.figure      = 'gcf';

title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%','Informative 40%','Informative 20%'};

% compare: 6 cue validity conditions
figure('Name','paperplot_CueValidity_v1_tfa_cvsi_spectrum_probelocked_6con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 250])
for channel = 7:12
    cfg.channel     = ci.label(channel);
    subplot(1,6,channel-6);    ft_singleplotTFR(cfg, ci);
    colormap(colormap2use);
    xline(1500,'--k','Probe Onset')
    title(title2use{channel-6});
    ylabel('Frequency (Hz)'); xlabel('Time (ms)');
    set(gca, 'FontSize', 10); set(gca,'LineWidth',2);
    %     colorbar;
end
% saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_probelocked_6con_s23', toadd1,toadd2,toadd3, '_colorbar'], '.jpg']);
saveas(gcf, [['paperplot_CueValidity_v1_tfa_cvsi_spectrum_probelocked_6con_s23', toadd1,toadd2,toadd3], '.jpg']);


%% 80% valid vs. invalid
%% Latency-analysis--Jackknife--80% valid vs. invalid
fsel = ci.freq >= 8 & ci.freq <= 12;

% average TFA across time over frequency level: change subxconxfreqxtime to subxconxtime
avealpha = [];
avealpha(:,1,:) = squeeze(mean(d3(:,9,fsel,:),3)); % 80valid
avealpha(:,2,:) = squeeze(mean(d3(:,12,fsel,:),3)); % 80invalid

%  % S1: calculate mean Saccade leave-one-out method
sub_num = size(avealpha,1);

saccade_jk = [];
for s = 1:sub_num
    for con = 1:2
        if s == 1
            saccade_jk(:,s,con) = mean(squeeze(avealpha([s+1:sub_num],con,:)),1);
        elseif s == sub_num
            saccade_jk(:,s,con) = mean(squeeze(avealpha([1:sub_num-1],con,:)),1);
        elseif s < sub_num && s > 1
            saccade_jk(:,s,con) = mean(squeeze(avealpha([1:s-1 s+1:sub_num],con,:)),1);
        end
    end
end


% % Jackknife--S2: find peak latency during a specific time range--50%
% define onset & offset latency 50% peak

% srate = 1; % sampling rate
srate = 9.9705; % sampling rate
% bl = abs(saccade.time(1)); % baseline length: ms
bl = abs(ci.time(1)); % baseline length: ms
method = 1; % 1-min; 2-max; min(for negative component)/max(for positive component)
% relativeperc = 0.3; % relative percentage of peak amplitude to detect the component onset 
relativeperc = 0.5; % relative percentage of peak amplitude to detect the component onset 
timerg = [1500 2340]; % interested time window of ERP component
% find min(for negative component) / max(for positive component) during the time window for each condition
for s = 1:size(saccade_jk,2) % sub_num
    for c = 1:size(saccade_jk,3) % condition_num
        grandave = saccade_jk(round((bl+timerg(1))/srate):round((bl+timerg(2))/srate),s,c);
        if method == 1
            peak_val = min(grandave);
        elseif method == 2
            peak_val = max(grandave);
        end
        grandave_diff = abs(grandave-relativeperc*peak_val); % calculate absulute difference of ERP & 50% of peak value
        
        peak_val_loc = find(grandave==peak_val); % find location of peak value
        
        % component onset & offset
        peak_val_loc_50onset = find(grandave_diff==min(grandave_diff(1:peak_val_loc))); % find the nearest location of the 50% of peak value before peak latency
        peak_val_loc_50offset = find(grandave_diff==min(grandave_diff(peak_val_loc:end))); % find the nearest location of the 50% of peak value before peak latency
        area_diff = peak_val_loc_50offset-peak_val_loc_50onset;

        saccade_jk_peak_val(s,c) = peak_val;

        saccade_jk_peak_lat(s,c) = (peak_val_loc-1)*srate+timerg(1); % get peak latency
        saccade_jk_peak_onset_lat(s,c) = (peak_val_loc_50onset-1)*srate+timerg(1); % get peak latency
        saccade_jk_peak_offset_lat(s,c) = (peak_val_loc_50offset-1)*srate+timerg(1); % get peak latency
        saccade_jk_area50_lat(s,c) = area_diff; % get 50% area latency
        clear grandave peak_val peak_val_loc peak_val_loc_50onset peak_val_loc_50offset grandave_diff area_diff
    end
end


% S3: Retrieved estamation (Smulders, F. T. Y. (2010))
% calculated retrieved latency for each sub
for c = 1:2
    for s = 1:sub_num
        saccade_jk_peak_lat_re(s,c) = sub_num*mean(saccade_jk_peak_lat(:,c),1)-(sub_num-1)*saccade_jk_peak_lat(s,c); % n*aver(lat)-(n-1)*lat_persub
        saccade_jk_peak_onset_lat_re(s,c) = sub_num*mean(saccade_jk_peak_onset_lat(:,c),1)-(sub_num-1)*saccade_jk_peak_onset_lat(s,c); % n*aver(lat)-(n-1)*lat_persub
        saccade_jk_peak_offset_lat_re(s,c) = sub_num*mean(saccade_jk_peak_offset_lat(:,c),1)-(sub_num-1)*saccade_jk_peak_offset_lat(s,c); % n*aver(lat)-(n-1)*lat_persub
        saccade_jk_area50_lat_re(s,c) = sub_num*mean(saccade_jk_area50_lat(:,c),1)-(sub_num-1)*saccade_jk_area50_lat(s,c); % n*aver(lat)-(n-1)*lat_persub
    end
end
 
% do t-test
[~,p_peak_lat_re_12,~,stats_peak_lat_re_12] = ttest(saccade_jk_peak_lat_re(:,1),saccade_jk_peak_lat_re(:,2));

[~,p_peak_onset_lat_re_12,~,stats_peak_onset_lat_re_12] = ttest(saccade_jk_peak_onset_lat_re(:,1),saccade_jk_peak_onset_lat_re(:,2));

[~,p_peak_offset_lat_re_12,~,stats_peak_offset_lat_re_12] = ttest(saccade_jk_peak_offset_lat_re(:,1),saccade_jk_peak_offset_lat_re(:,2));

[~,p_area50_lat_re_12,~,stats_area50_lat_re_12] = ttest(saccade_jk_area50_lat_re(:,1),saccade_jk_area50_lat_re(:,2));

save(fullfile('CueVa_v1_probe_alpha_later_latency_50%_results_Smulders_s23_80%vavsinva.mat'))


