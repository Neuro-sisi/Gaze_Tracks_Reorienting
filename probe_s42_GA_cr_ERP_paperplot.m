%% Plot ERPs for paper

clc;clear
close all

% % parameters
laplacian           = 0;
chancluster         = 0;
% for timebl        = [1 2 3]; % 1-precuebl: -0.2-0; 2-prefixbl: -1.2--1; 3-preprobebl: 1.3--1.5
timebl              = [3]; % 1-precuebl: -0.2-0; 2-prefixbl: -1.2--1; 3-preprobebl: 1.3--1.5

pp2do               = [1:23];

nsmooth             = 50; % 0-no smooth, 20-smooth erp

plotRunningTopos    = 1;
plot_sech           = 1;
plot_sech_probe     = 1;

colormap2use        = fliplr(brewermap(100, 'RdBu'));


xlim2plot           = [-200 1500];
xlim2plot_probe     = [1300 2400];


%% load
if laplacian toadd1 = '_laplacian'; else toadd1 = ''; end
if chancluster toadd2 = ''; else toadd2 = '_sinch'; end
if timebl==1 toadd3 = '_precuebl'; elseif timebl==2 toadd3 = '_prefixbl'; elseif timebl==3 toadd3 = '_preprobebl'; end

s = 0;
for pp = pp2do;
    s = s+1;
    param = getSubjParam(pp);

    disp(['getting data from ' param.subjName]);

    load([param.pathnewlen, 'saved_data/erp', toadd1,toadd2,toadd3, '__' param.subjName], 'erp','ci','ci_allchan');

    if nsmooth > 0
        for x1 = 1:size(erp.imper100,1);
            erp.imper100(x1,:) = gsmooth(squeeze(erp.imper100(x1,:)), round(nsmooth/4));
            erp.high100(x1,:) = gsmooth(squeeze(erp.high100(x1,:)), round(nsmooth/4));
            erp.med80(x1,:) = gsmooth(squeeze(erp.med80(x1,:)), round(nsmooth/4));
            erp.low60(x1,:) = gsmooth(squeeze(erp.low60(x1,:)), round(nsmooth/4));

            erp.valid(x1,:) = gsmooth(squeeze(erp.valid(x1,:)), round(nsmooth/4));
            erp.invalid(x1,:) = gsmooth(squeeze(erp.invalid(x1,:)), round(nsmooth/4));

            erp.validim100(x1,:) = gsmooth(squeeze(erp.validim100(x1,:)), round(nsmooth/4));
            erp.validinfor100(x1,:) = gsmooth(squeeze(erp.validinfor100(x1,:)), round(nsmooth/4));
            erp.validinfor80(x1,:) = gsmooth(squeeze(erp.validinfor80(x1,:)), round(nsmooth/4));
            erp.validinfor60(x1,:) = gsmooth(squeeze(erp.validinfor60(x1,:)), round(nsmooth/4));
            erp.validinfor40(x1,:) = gsmooth(squeeze(erp.validinfor40(x1,:)), round(nsmooth/4));
            erp.validinfor20(x1,:) = gsmooth(squeeze(erp.validinfor20(x1,:)), round(nsmooth/4));

            erp.imper_vs_high(x1,:) = gsmooth(squeeze(erp.imper_vs_high(x1,:)), round(nsmooth/4));
            erp.high_vs_low(x1,:) = gsmooth(squeeze(erp.high_vs_low(x1,:)), round(nsmooth/4));
            erp.valid_vs_invalid(x1,:) = gsmooth(squeeze(erp.valid_vs_invalid(x1,:)), round(nsmooth/4));
        end
        for x1 = 1:size(ci.data,1);
            ci.data(x1,:) = gsmooth(squeeze(ci.data(x1,:)), round(nsmooth/4));
        end
        for x1 = 1:size(ci_allchan.data,1); % conditions
            for c1 = 1:size(ci_allchan.data,2); % channels
                ci_allchan.data(x1,c1,:) = gsmooth(squeeze(ci_allchan.data(x1,c1,:)), round(nsmooth/4));
            end
        end
    end

    d1(s,:,:)  = erp.imper100;
    d2(s,:,:)  = erp.high100;
    d3(s,:,:)  = erp.med80;
    d4(s,:,:)  = erp.low60;
    d5(s,:,:)  = erp.valid;
    d6(s,:,:)  = erp.invalid;
    d7(s,:,:)  = erp.validim100;
    d8(s,:,:)  = erp.validinfor100;
    d9(s,:,:)  = erp.validinfor80;
    d10(s,:,:) = erp.validinfor60;
    d11(s,:,:) = erp.validinfor40;
    d12(s,:,:) = erp.validinfor20;
    d13(s,:,:) = erp.imper_vs_high;
    d14(s,:,:) = erp.high_vs_low;
    d15(s,:,:) = erp.valid_vs_invalid;
    d16(s,:,:) = ci.data;
    d17(s,:,:,:) = ci_allchan.data;

end


%% put back into structure
% calculate valid-invalid difference for 80 & 60, respectively
for ss = 1:s
    dvainva_80(s,:,:)   = d9(s,:,:) - d12(s,:,:);
    dvainva_60(s,:,:)   = d10(s,:,:) - d11(s,:,:);
end

erp.imper100        = squeeze(mean(d1));
erp.high100         = squeeze(mean(d2));
erp.med80           = squeeze(mean(d3));
erp.low60           = squeeze(mean(d4));
erp.valid           = squeeze(mean(d5));
erp.invalid         = squeeze(mean(d6));
erp.validim100      = squeeze(mean(d7));
erp.validinfor100   = squeeze(mean(d8));
erp.validinfor80    = squeeze(mean(d9));
erp.validinfor60    = squeeze(mean(d10));
erp.validinfor40    = squeeze(mean(d11));
erp.validinfor20    = squeeze(mean(d12));
erp.imper_vs_high   = squeeze(mean(d13));
erp.high_vs_low     = squeeze(mean(d14));
erp.valid_vs_invalid= squeeze(mean(d15));
ci.data             = squeeze(mean(d16));
ci_allchan.data     = squeeze(mean(d17));

erp.vainva_80  = squeeze(mean(dvainva_80));
erp.vainva_60  = squeeze(mean(dvainva_60));


%% plot ---------------------------------------------------------------------------------------


%% Cue-locked ERPs

%% Lateralized ERP
title2use = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];

% % 1. lateralized ERP--contra & ipsi
figure('Name','CueValidity_v1_ERP_cuelocked_cvsi_4con','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])

for sp = 1:4
    if sp == 1, con = 0; elseif sp == 2, con = 3; elseif sp == 3, con = 6; elseif sp == 4, con = 9; end
    subplot(1,4,sp); hold on; title(title2use{sp});
    p1 = frevede_errorbarplot(ci.time, d16(:,1+con,:), [0, 0, 0], 'se');
    p2 = frevede_errorbarplot(ci.time, d16(:,2+con,:), [0.7,0.7,0.7], 'se');
    plot(xlim, [0,0], '--k');xlim(xlim2plot); ylim([-2 5]);
    xline(0,'--k','Cue Onset');
    ylabel('Amplitude'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
legend([p1, p2], {'Contralateral','Ipsilateral'},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_cvsi_4con_s23', toadd1,toadd2,toadd3], '.jpg']);


% % 2. lateralized ERP--contra-ipsi seperate plots
figure('Name','CueValidity_v1_ERP_cuelocked_cmi_4con_sepcon','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 1600 300])
for sp = 1:4
    subplot(1,4,sp); hold on; title(title2use{sp});
    frevede_errorbarplot(ci.time, d16(:,3*sp,:), [col2use(sp,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-2 1]);
    xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
    ylabel('Amplitude'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
end
legend({'Lateralized ERP'},'FontSize',8,'AutoUpdate','off','Box','off');

saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_cmi_4con_sepcon_s23', toadd1,toadd2,toadd3], '.jpg']);


% % 3.lateralized ERP--contra-ipsi overlay
% overlay
% col2use = [0,0,1; 1,0,0; 1,0,1; 0,1,1;];
col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];
figure('Name','CueValidity_v1_ERP_cuelocked_cmi_4con_overlay','NumberTitle','off', 'Color','white'),
set(gcf,'Position',[0 0 400 300])
hold on;
p1 = frevede_errorbarplot(ci.time, d16(:,3,:), [col2use(1,:)], 'se');
p2 = frevede_errorbarplot(ci.time, d16(:,6,:), [col2use(2,:)], 'se');
p3 = frevede_errorbarplot(ci.time, d16(:,9,:), [col2use(3,:)], 'se');
p4 = frevede_errorbarplot(ci.time, d16(:,12,:), [col2use(4,:)], 'se');
legend([p1,p2,p3,p4],title2use{1:4},'FontSize',6,'AutoUpdate','off','Box','off');
plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-2 1]);
xline(0,'--k','Cue Onset');set(gca, 'YDir','reverse');
ylabel('Amplitude'); xlabel('Time (ms)');
set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
title('Lateralized ERP');

saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_cmi_4con_overlay_s23', toadd1,toadd2,toadd3], '.jpg']);


%% Cue-locked ERPs

%% Averaged ERP--center electrodes
if plot_sech

    title2uselabel = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};

    % plot timeseries
    figure;
    set(gcf,'Position',[0 0 1600 600])
    for ch = 1:8
        if     ch == 1     choi = ismember(erp.label, 'Oz'); title2use = 'Oz';
        elseif ch == 2     choi = ismember(erp.label, 'POz'); title2use = 'POz';
        elseif ch == 3     choi = ismember(erp.label, 'Pz'); title2use = 'Pz';
        elseif ch == 4     choi = ismember(erp.label, 'CPz'); title2use = 'CPz';
        elseif ch == 5     choi = ismember(erp.label, 'Cz'); title2use = 'Cz';
        elseif ch == 6     choi = ismember(erp.label, 'FCz'); title2use = 'FCz';
        elseif ch == 7     choi = ismember(erp.label, 'Fz'); title2use = 'Fz';
        elseif ch == 8     choi = ismember(erp.label, 'AFz'); title2use = 'AFz';
        end

        subplot(2,4,ch); hold on; title(title2use);
        p1 = frevede_errorbarplot(ci.time, d1(:,choi,:), [col2use(1,:)], 'se');
        p2 = frevede_errorbarplot(ci.time, d2(:,choi,:), [col2use(2,:)], 'se');
        p3 = frevede_errorbarplot(ci.time, d3(:,choi,:), [col2use(3,:)], 'se');
        p4 = frevede_errorbarplot(ci.time, d4(:,choi,:), [col2use(4,:)], 'se');
        plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-10 4]);
        xline(0,'--k','Cue Onset');
        ylabel('Amplitude'); xlabel('Time (ms)');
        set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    end
    legend([p1,p2,p3,p4],title2uselabel{1:4},'FontSize',6,'AutoUpdate','off','Box','off');

    saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_centerelec_s23', toadd1,toadd2,toadd3], '.jpg']);

end



%% Cue-locked ERPs
%% Averaged ERP topo--all timerange

if plotRunningTopos
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.marker = 'off';
    cfg.zlim = [-10 10]; if laplacian cfg.zlim = [-.1 .1]*1e-3; end

    %     stepsize = 200; times1 = [0:stepsize:1000]; n = length(times1);
    stepsize = 400; times1 = [0:stepsize:800]; n = length(times1);

    figure;
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(4,n-1,t);                 cfg.parameter = 'imper100'; ft_topoplotER(cfg, erp); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
        subplot(4,n-1,t+n-1);         cfg.parameter = 'high100';  ft_topoplotER(cfg, erp);
        subplot(4,n-1,t+(n-1)*2);   cfg.parameter = 'med80';    ft_topoplotER(cfg, erp);
        subplot(4,n-1,t+(n-1)*3);   cfg.parameter = 'low60';    ft_topoplotER(cfg, erp);
    end
    %     colorbar
    colormap(colormap2use);

    %     saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_topo_01000_s23', toadd1,toadd2,toadd3], '.jpg']);
    saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_topo_0800_s23', toadd1,toadd2,toadd3], '.jpg']);

end

%% Cue-locked topo--300-800ms
if plotRunningTopos
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.marker = 'off';
    cfg.zlim = [-4 4];
    if laplacian cfg.zlim = [-.1 .1]*1e-3; end

    %     stepsize = 200; times1 = [0:stepsize:1000]; n = length(times1);
    %     stepsize = 400; times1 = [0:stepsize:800]; n = length(times1);

    %     figure('Name','CueValidity_v1_ERP_cuelocked_aveERP_topo_400800','NumberTitle','off', 'Color','white'),
    figure('Name','CueValidity_v1_ERP_cuelocked_aveERP_topo_300800','NumberTitle','off', 'Color','white'),
    set(gcf,'Position',[0 0 1600 300])

    %     cfg.xlim = [400 800];
    cfg.xlim = [300 800];
    subplot(1,4,1); cfg.parameter = 'imper100'; ft_topoplotER(cfg, erp); title(title2uselabel{1});
    subplot(1,4,2); cfg.parameter = 'high100';   ft_topoplotER(cfg, erp); title(title2uselabel{2});
    subplot(1,4,3); cfg.parameter = 'med80';      ft_topoplotER(cfg, erp); title(title2uselabel{3});
    subplot(1,4,4); cfg.parameter = 'low60';       ft_topoplotER(cfg, erp); title(title2uselabel{4});
    %         colorbar
    colormap(colormap2use);

    %     saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_topo_400800_cobar_s23', toadd1,toadd2,toadd3], '.jpg']);
    %     saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_topo_400800_s23', toadd1,toadd2,toadd3], '.jpg']);
    %         saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_topo_300800_cobar_s23', toadd1,toadd2,toadd3], '.jpg']);
    saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_topo_300800_s23', toadd1,toadd2,toadd3], '.jpg']);

end


%% Cue-locked ERP
%% Averaged ERP--POz (based on topo distribution)
if plot_sech

    title2uselabel = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
    col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];


    % % seperate plots
    figure('Name','CueValidity_v1_ERP_cuelocked_aveERP_POz_sepcon','NumberTitle','off', 'Color','white'),
    set(gcf,'Position',[0 0 1600 300])

    choi = ismember(erp.label, 'POz');

    for sp = 1:4
        if sp == 1, dd = d1; elseif sp == 2, dd = d2; elseif sp == 3, dd = d3; elseif sp == 4, dd = d4; end
        subplot(1,4,sp); hold on; title(title2uselabel{sp});
        frevede_errorbarplot(ci.time, dd(:,choi,:), [col2use(sp,:)], 'se');
        plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-8 4]);
        xline(0,'--k','Cue Onset');
        ylabel('Amplitude'); xlabel('Time (ms)');
        set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    end
    legend({'ERP-POz'},'FontSize',6,'AutoUpdate','off','Box','off');

    saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_POz_sepcon_s23', toadd1,toadd2,toadd3], '.jpg']);



    % plot overlay
    figure('Name','CueValidity_v1_ERP_cuelocked_aveERP_POz_overlay','NumberTitle','off', 'Color','white'),
    set(gcf,'Position',[0 0 400 300])

    choi = ismember(erp.label, 'POz');

    hold on; title('ERP-POz');
    p1 = frevede_errorbarplot(ci.time, d1(:,choi,:), [col2use(1,:)], 'se');
    p2 = frevede_errorbarplot(ci.time, d2(:,choi,:), [col2use(2,:)], 'se');
    p3 = frevede_errorbarplot(ci.time, d3(:,choi,:), [col2use(3,:)], 'se');
    p4 = frevede_errorbarplot(ci.time, d4(:,choi,:), [col2use(4,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot); ylim([-8 4]);
    xline(0,'--k','Cue Onset');
    ylabel('Amplitude'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);

    legend([p1,p2,p3,p4],title2uselabel{1:4},'FontSize',6,'AutoUpdate','off','Box','off');

    saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_POz_overlay_s23', toadd1,toadd2,toadd3], '.jpg']);

end


%% Cue-locked ERP
%% plot bar of averaged ERPs--P3 component
if plot_sech

    col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1; 0.9,0.5,0.1;0.72 0.27 1;0.56 0.2 0.02;];

    choi = ismember(erp.label, 'POz');

    timeproon = 0;
    %     timerange = [400+timeproon 800+timeproon];
    timerange = [300+timeproon 800+timeproon];

    timediff1 = abs(ci.time-timerange(1));
    pnt_start = find(timediff1==min(timediff1));
    timediff2 = abs(ci.time-timerange(2));
    pnt_end = find(timediff2==min(timediff2));
    clear timediff1 timediff2

    % average ERP across time
    aveerp = [];
    aveerp(:,1) = mean(squeeze(d1(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,2) = mean(squeeze(d2(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,3) = mean(squeeze(d3(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,4) = mean(squeeze(d4(:,choi,pnt_start:pnt_end)),2);

    % bar plot
    %     figure('Name','CueValidity_v1_ERP_cuelocked_aveERP_POz_bar400800','NumberTitle','off', 'Color','white'),
    figure('Name','CueValidity_v1_ERP_cuelocked_aveERP_POz_bar300800','NumberTitle','off', 'Color','white'),
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
    ylim([-4 2]);
    ylabel('Mean Amplitude'); xlabel('Cue Type');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);

    %     saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_POz_bar400800_s23', toadd1,toadd2,toadd3], '.jpg']);
    saveas(gcf, [['paperplot_CueValidity_v1_ERP_cuelocked_aveERP_POz_bar300800_s23', toadd1,toadd2,toadd3], '.jpg']);
    clear aveerp b

end




%% Probe-locked--ERP
%% Averaged ERP--center electrodes
% !! redo baseline correction (pre-test baseline)

if plot_sech_probe

    title2uselabel = {'Imperative 100%','Informative 100%','Informative 80%-Valid','Informative 80%-Invalid','Informative 60%-Valid','Informative 60%-Invalid'};
    col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1;0.72 0.27 1; 0.9,0.5,0.1;0.56 0.2 0.02;];

    figure;
    set(gcf,'Position',[0 0 1600 600])
    for ch = 1:8
        if     ch == 1     choi = ismember(erp.label, 'Oz'); title2use = 'Oz';
        elseif ch == 2     choi = ismember(erp.label, 'POz'); title2use = 'POz';
        elseif ch == 3     choi = ismember(erp.label, 'Pz'); title2use = 'Pz';
        elseif ch == 4     choi = ismember(erp.label, 'CPz'); title2use = 'CPz';
        elseif ch == 5     choi = ismember(erp.label, 'Cz'); title2use = 'Cz';
        elseif ch == 6     choi = ismember(erp.label, 'FCz'); title2use = 'FCz';
        elseif ch == 7     choi = ismember(erp.label, 'Fz'); title2use = 'Fz';
        elseif ch == 8     choi = ismember(erp.label, 'AFz'); title2use = 'AFz';
        end

        subplot(2,4,ch); hold on; title(title2use);
%         p1 = frevede_errorbarplot(ci.time, d7(:,choi,:), [col2use(1,:)], 'se');
%         p2 = frevede_errorbarplot(ci.time, d8(:,choi,:), [col2use(2,:)], 'se');
%         p3 = frevede_errorbarplot(ci.time, d9(:,choi,:), [col2use(3,:)], 'se');
%         p4 = frevede_errorbarplot(ci.time, d12(:,choi,:), [col2use(4,:)], 'se');
        p5 = frevede_errorbarplot(ci.time, d10(:,choi,:), [col2use(5,:)], 'se');
        p6 = frevede_errorbarplot(ci.time, d11(:,choi,:), [col2use(6,:)], 'se');
        plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-4 16]);
        xline(1500,'--k','Test Onset');
        ylabel('Amplitude'); xlabel('Time (ms)');
        set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    end
    legend([p1,p2,p3,p4,p5,p6],title2uselabel{1:6},'FontSize',6,'AutoUpdate','off','Box','off');

    saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_centerelec_s23', toadd1,toadd2,toadd3], '.jpg']);

end



%% Probe-locked--ERP
%% Averaged ERP topo--all timerange

if plotRunningTopos
    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.marker = 'off';
    cfg.zlim = [-5 5];
    if laplacian cfg.zlim = [-.1 .1]*1e-3; end

    %         stepsize = 200; times1 = [0:stepsize:1000]; n = length(times1);
    stepsize = 200; times1 = [1500:stepsize:2300]; n = length(times1);

    figure(1);
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(6,n-1,t);                 cfg.parameter = 'validim100';     ft_topoplotER(cfg, erp); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
        subplot(6,n-1,t+n-1);         cfg.parameter = 'validinfor100';  ft_topoplotER(cfg, erp);
        subplot(6,n-1,t+(n-1)*2);   cfg.parameter = 'validinfor80';    ft_topoplotER(cfg, erp);
        subplot(6,n-1,t+(n-1)*3);   cfg.parameter = 'validinfor20';    ft_topoplotER(cfg, erp);
        subplot(6,n-1,t+(n-1)*4);   cfg.parameter = 'validinfor60';    ft_topoplotER(cfg, erp);
        subplot(6,n-1,t+(n-1)*5);   cfg.parameter = 'validinfor40';    ft_topoplotER(cfg, erp);
    end
    colorbar
    colormap(colormap2use);

      cfg.zlim = [-0.5 0.5];
  figure(2);
    for t = 1:n-1
        cfg.xlim = times1(t:t+1);
        subplot(2,n-1,t);                 cfg.parameter = 'vainva_80';     ft_topoplotER(cfg, erp); title([num2str(times1(t)), ' to ', num2str(times1(t+1))]);
         subplot(2,n-1,t+n-1);         cfg.parameter = 'vainva_60';  ft_topoplotER(cfg, erp);
   end
    colorbar
    colormap(colormap2use);

%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_topo_01000_s23', toadd1,toadd2,toadd3], '.jpg']);
%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_topo_0800_s23', toadd1,toadd2,toadd3], '.jpg']);

end


%% Probe-locked--ERP
% %% topo--400-800ms
%% topo--300-600ms
if plotRunningTopos

    title2uselabel = {'Imperative 100%','Informative 100%','Informative 80%-Valid','Informative 80%-Invalid','Informative 60%-Valid','Informative 60%-Invalid'};

    cfg = [];
    cfg.layout = 'easycapM1.mat';
    cfg.comment = 'no';
    cfg.style = 'straight';
    cfg.marker = 'off';
    cfg.zlim = [-7 7];
    if laplacian cfg.zlim = [-.1 .1]*1e-3; end

    figure('Name','CueValidity_v1_ERP_probelocked_aveERP_topo_400800','NumberTitle','off', 'Color','white'),
    set(gcf,'Position',[0 0 1600 200])

    bl = 1500;
    cfg.xlim = [300+bl 600+bl];
    subplot(1,6,1); cfg.parameter = 'validim100'; ft_topoplotER(cfg, erp); title(title2uselabel{1});
    subplot(1,6,2); cfg.parameter = 'validinfor100';   ft_topoplotER(cfg, erp); title(title2uselabel{2});
    subplot(1,6,3); cfg.parameter = 'validinfor80';      ft_topoplotER(cfg, erp); title(title2uselabel{3});
    subplot(1,6,4); cfg.parameter = 'validinfor20';       ft_topoplotER(cfg, erp); title(title2uselabel{4});
    subplot(1,6,5); cfg.parameter = 'validinfor60';       ft_topoplotER(cfg, erp); title(title2uselabel{5});
    subplot(1,6,6); cfg.parameter = 'validinfor40';       ft_topoplotER(cfg, erp); title(title2uselabel{6});
    colorbar
    colormap(colormap2use);

    %     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_topo_400800_cobar_s23', toadd1,toadd2,toadd3], '.jpg']);
%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_topo_400800_s23', toadd1,toadd2,toadd3], '.jpg']);
    saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_topo_300600_colorbar_s23', toadd1,toadd2,toadd3], '.jpg']);
%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_topo_300600_s23', toadd1,toadd2,toadd3], '.jpg']);

end


%% Probe-locked--ERP
% %% Averaged ERP--Pz (based on topo distribution)
%% Averaged ERP--CPz (based on topo distribution)
if plot_sech_probe

    title2uselabel1 = {'Imperative 100%','Informative 100%','Informative 80%','Informative 60%'};
    title2uselabel = {'Imperative 100%','Informative 100%','Informative 80%-Valid','Informative 80%-Invalid','Informative 60%-Valid','Informative 60%-Invalid'};
    col2use = [0,0.3,0.9; 0,0.6,0.6; 0.9,0.6,1;0.72 0.27 1; 0.9,0.5,0.1;0.56 0.2 0.02;];

%     choi = ismember(erp.label, 'Pz');
    choi = ismember(erp.label, 'CPz');

%     % % seperate plots
% %     figure('Name','CueValidity_v1_ERP_probelocked_aveERP_Pz_sepcon','NumberTitle','off', 'Color','white'),
%     figure('Name','CueValidity_v1_ERP_probelocked_aveERP_Pz_sepcon','NumberTitle','off', 'Color','white'),
%     set(gcf,'Position',[0 0 1600 200])
%     for sp = 1:6
%         if sp == 1, dd = d7; elseif sp == 2, dd = d8; elseif sp == 3, dd = d9; elseif sp == 4, dd = d12; elseif sp == 5, dd = d10; elseif sp == 6, dd = d11; end
%         subplot(1,6,sp); hold on; title(title2uselabel{sp});
%         frevede_errorbarplot(ci.time, dd(:,choi,:), [col2use(sp,:)], 'se');
%         plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-4 16]);
%         xline(1500,'--k','Test Onset');
%         ylabel('Amplitude'); xlabel('Time (ms)');
%         set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
%     end
% %     legend({'ERP-Pz'},'FontSize',6,'AutoUpdate','off','Box','off');
% %     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_Pz_sepcon_s23', toadd1,toadd2,toadd3], '.jpg']);
%     legend({'ERP-CPz'},'FontSize',6,'AutoUpdate','off','Box','off');
%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_CPz_sepcon_s23', toadd1,toadd2,toadd3], '.jpg']);


    % % seperate plots--4plots
%     figure('Name','CueValidity_v1_ERP_probelocked_aveERP_Pz_sepcon_4plots','NumberTitle','off', 'Color','white'),
    figure('Name','CueValidity_v1_ERP_probelocked_aveERP_CPz_sepcon_4plots','NumberTitle','off', 'Color','white'),
    set(gcf,'Position',[0 0 1600 300])

    for sp = 1:2
        if sp == 1, dd = d7; elseif sp == 2, dd = d8; elseif sp == 3, dd = d9; elseif sp == 4, dd = d12; elseif sp == 5, dd = d10; elseif sp == 6, dd = d11; end
        subplot(1,4,sp); hold on; title(title2uselabel1{sp});
        frevede_errorbarplot(ci.time, dd(:,choi,:), [col2use(sp,:)], 'se');
        plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-4 12]);
        xline(1500,'--k','Test Onset');
        ylabel('Amplitude'); xlabel('Time (ms)');
        set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
        legend({'Valid'},'FontSize',8,'AutoUpdate','off','Box','off');
    end
    hold on;
    subplot(1,4,3); hold on; title(title2uselabel1{3});
    p1=frevede_errorbarplot(ci.time, d9(:,choi,:), [col2use(3,:)], 'se');
    p2=frevede_errorbarplot(ci.time, d12(:,choi,:), [col2use(4,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-4 12]);
    xline(1500,'--k','Test Onset');
    ylabel('Amplitude'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

    subplot(1,4,4); hold on; title(title2uselabel1{4});
    p1=frevede_errorbarplot(ci.time, d10(:,choi,:), [col2use(5,:)], 'se');
    p2=frevede_errorbarplot(ci.time, d11(:,choi,:), [col2use(6,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-4 12]);
    xline(1500,'--k','Test Onset');
    ylabel('Amplitude'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);
    legend([p1,p2],{'Valid','Invalid'},'FontSize',8,'AutoUpdate','off','Box','off');

%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_Pz_sepcon_4plots_s23', toadd1,toadd2,toadd3], '.jpg']);
saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_CPz_sepcon_4plots_s23', toadd1,toadd2,toadd3], '.jpg']);
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValidity_v1_ERP_probelocked_aveERP_CPz_sepcon_4plots_s23', 'epsc')



    % plot overlay
%     figure('Name','CueValidity_v1_ERP_probelocked_aveERP_Pz_overlay','NumberTitle','off', 'Color','white'),
    figure('Name','CueValidity_v1_ERP_probelocked_aveERP_CPz_overlay','NumberTitle','off', 'Color','white'),
    set(gcf,'Position',[0 0 400 300])
    hold on; title('ERP-Reporting Stage');
    p1 = frevede_errorbarplot(ci.time, d7(:,choi,:), [col2use(1,:)], 'se');
    p2 = frevede_errorbarplot(ci.time, d8(:,choi,:), [col2use(2,:)], 'se');
    p3 = frevede_errorbarplot(ci.time, d9(:,choi,:), [col2use(3,:)], 'se');
    p4 = frevede_errorbarplot(ci.time, d12(:,choi,:), [col2use(4,:)], 'se');
    p5 = frevede_errorbarplot(ci.time, d10(:,choi,:), [col2use(5,:)], 'se');
    p6 = frevede_errorbarplot(ci.time, d11(:,choi,:), [col2use(6,:)], 'se');
    plot(xlim, [0,0], '--k'); xlim(xlim2plot_probe); ylim([-4 16]);
    xline(1500,'--k','Test Onset');
    ylabel('Amplitude'); xlabel('Time (ms)');
    set(gca, 'FontSize', 12); set(gca,'LineWidth',3);

    legend([p1,p2,p3,p4,p5,p6],title2uselabel{1:6},'FontSize',6,'AutoUpdate','off','Box','off');

%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_Pz_overlay_s23', toadd1,toadd2,toadd3], '.jpg']);
saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_CPz_overlay_s23', toadd1,toadd2,toadd3], '.jpg']);
set(gcf,'renderer','painters')
saveas(gcf, 'paperplot_CueValidity_v1_ERP_probelocked_aveERP_CPz_overlay_s23', 'epsc')

end


%% Probe-locked--ERP
%% Pz & CPz-erp permutation cluster test
% choi = ismember(erp.label, 'Pz');
choi = ismember(erp.label, 'CPz');
pnt_start = 718;
pnt_end = 1024; % 1300-2500ms

data_cond1 = squeeze(d9(:,choi,pnt_start:pnt_end)); % 80-valid
data_cond2 = squeeze(d12(:,choi,pnt_start:pnt_end)); % 80-invalid
data_cond3 = squeeze(d10(:,choi,pnt_start:pnt_end)); % 60-valid
data_cond4 = squeeze(d11(:,choi,pnt_start:pnt_end)); % 60-invalid
data_zero = zeros(size(data_cond1)); % compare with 0

statcfg.xax = ci.time(pnt_start:pnt_end);
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.025; % two-sided
statcfg.nsub = s;
statcfg.statMethod = 'montecarlo'; % or 'analytic'; 

% cluster test
stat12=[];stat34=[];
stat12 = frevede_ftclusterstat1D(statcfg, data_cond1, data_cond2);
stat34 = frevede_ftclusterstat1D(statcfg, data_cond3, data_cond4);

signegclurange12=[]; signegclurange34=[]; sigposclurange12=[]; sigposclurange34=[];
signegclurange12 = statcfg.xax(find(stat12.negclusterslabelmat==1));
sigposclurange12 = statcfg.xax(find(stat12.posclusterslabelmat==1));
signegclurange34 = statcfg.xax(find(stat34.negclusterslabelmat==1));
sigposclurange34 = statcfg.xax(find(stat34.posclusterslabelmat==1));


%% Probe-locked--ERP
%% Posterior cluster-erp permutation cluster test
% choi = ismember(erp.label, 'Pz');
choi = ismember(erp.label, {'CPz'});
pnt_start = 718;
pnt_end = 1024; % 1300-2500ms

% data_cond1 = squeeze(d9(:,choi,pnt_start:pnt_end)); % 80-valid
% data_cond2 = squeeze(d12(:,choi,pnt_start:pnt_end)); % 80-invalid
% data_cond3 = squeeze(d10(:,choi,pnt_start:pnt_end)); % 60-valid
% data_cond4 = squeeze(d11(:,choi,pnt_start:pnt_end)); % 60-invalid
data_cond1 = squeeze(mean(d9(:,choi,pnt_start:pnt_end),2)); % 80-valid
data_cond2 = squeeze(mean(d12(:,choi,pnt_start:pnt_end),2)); % 80-invalid
data_cond3 = squeeze(mean(d10(:,choi,pnt_start:pnt_end),2)); % 60-valid
data_cond4 = squeeze(mean(d11(:,choi,pnt_start:pnt_end),2)); % 60-invalid
data_zero = zeros(size(data_cond1)); % compare with 0

statcfg.xax = ci.time(pnt_start:pnt_end);
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.025; % two-sided
statcfg.nsub = s;
statcfg.statMethod = 'montecarlo'; % or 'analytic'; 

% cluster test
stat12=[];stat34=[];
stat12 = frevede_ftclusterstat1D(statcfg, data_cond1, data_cond2);
stat34 = frevede_ftclusterstat1D(statcfg, data_cond3, data_cond4);

signegclurange12=[]; signegclurange34=[]; sigposclurange12=[]; sigposclurange34=[];
signegclurange12 = statcfg.xax(find(stat12.negclusterslabelmat==1));
sigposclurange12 = statcfg.xax(find(stat12.posclusterslabelmat==1));
signegclurange34 = statcfg.xax(find(stat34.negclusterslabelmat==1));
sigposclurange34 = statcfg.xax(find(stat34.posclusterslabelmat==1));



%% Probe-locked--ERP
%% plot bar of averaged ERPs--P3 component
if plot_sech_probe

%     choi = ismember(erp.label, 'Pz');
    choi = ismember(erp.label, 'CPz');

    timeproon = 1500;
    timerange = [300+timeproon 600+timeproon];
    timediff1 = abs(ci.time-timerange(1));
    pnt_start = find(timediff1==min(timediff1));
    timediff2 = abs(ci.time-timerange(2));
    pnt_end = find(timediff2==min(timediff2));
    clear timediff1 timediff2

    % average ERP across time
    aveerp = [];
    aveerp(:,1) = mean(squeeze(d7(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,2) = mean(squeeze(d8(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,3) = mean(squeeze(d9(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,4) = mean(squeeze(d12(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,5) = mean(squeeze(d10(:,choi,pnt_start:pnt_end)),2);
    aveerp(:,6) = mean(squeeze(d11(:,choi,pnt_start:pnt_end)),2);

    [~,P_34,~,STATS_34] = ttest(aveerp(:,3),aveerp(:,4));
    [~,P_56,~,STATS_56] = ttest(aveerp(:,5),aveerp(:,6));

    % bar plot
%     figure('Name','CueValidity_v1_ERP_probelocked_aveERP_Pz_bar300600','NumberTitle','off', 'Color','white'),
    figure('Name','CueValidity_v1_ERP_probelocked_aveERP_CPz_bar300600','NumberTitle','off', 'Color','white'),
    set(gcf,'Position',[0 0 400 400])
    hold on,
    b = bar(1:6, mean(aveerp(:,:)));
    b.LineWidth = 2;
    b.FaceColor = 'flat';
    for c = 1:6
        b.CData(c,:) = col2use(c,:);
    end
    hold on,
    errorbar(1:6, mean(aveerp(:,:)), std(aveerp(:,:))./sqrt((size(aveerp,1))), '.k', 'LineWidth',2);
    %     xticks(1:6); xticklabels({'Imperative 100%','Informative 100%','Informative 80%-Valid','Informative 80%-Invalid','Informative 60%-Valid','Informative 60%-Invalid'});
    xticks(1:6); xticklabels({'Im100%','100%','80%-Valid','80%-Invalid','60%-Valid','60%-Invalid'});
    ylim([0 15]);
    ylabel('Mean Amplitude'); xlabel('Cue Type');
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3);

%     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_Pz_bar300600_s23', toadd1,toadd2,toadd3], '.jpg']);
    saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_CPz_bar300600_s23', toadd1,toadd2,toadd3], '.jpg']);
     saveas(gcf, [['paperplot_CueValidity_v1_ERP_probelocked_aveERP_CPz_bar300600_s23', toadd1,toadd2,toadd3], '.epsc']);
   clear aveerp b

end



