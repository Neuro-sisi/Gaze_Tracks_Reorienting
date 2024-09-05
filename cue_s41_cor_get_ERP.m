%% Cue-locked ERP calculation

% correction 1 Apr 2024, Sisi Wang
% reverse cueside for invalid cue trials after cue

%%
clc;clear ;close all
ft_defaults

%% parameters
sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'28';'29';'30';'31';'32';'33';'34';};

for laplacian    = 0; % defalt--0
    for chancluster  = 0; % defalt--0; 1-averaged over cluster channels; 0-single channel PO7/8
        % for timebl     = [1 2 3]; % defalt--1; 1-precue baseline: -0.2-0; 2-prefixation baseline: -1.2--1; 3-preprobe baseline: 1.3--1.5
        for timebl     = 1; % defalt--1; 1-precue baseline: -0.2-0; 2-prefixation baseline: -1.2--1; 3-preprobe baseline: 1.3--1.5

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

                valid = ismember(tl.trialinfo(:,3), 1);
                invalid  = ismember(tl.trialinfo(:,3), 2);

                % reverse cueside for invalid trials
                % left     = ismember(tl.trialinfo(:,1), [41]);
                % right    = ismember(tl.trialinfo(:,1), [42]);
                trialinfo = tl.trialinfo;
                for t = 1:length(trialinfo)
                    if trialinfo(t,3) == 1 % valid
                        trialinfo(t,4) = trialinfo(t,1)-40; % valid
                    elseif trialinfo(t,3) == 2 % invalid
                        trialinfo(t,4) = 3-(trialinfo(t,1)-40); % valid
                    end
                end
                tl.trialinfo = trialinfo; clear trialinfo

                left     = ismember(tl.trialinfo(:,4), [1]);
                right    = ismember(tl.trialinfo(:,4), [2]);



                %% ERP per condition, put into a structure we can later plot
                erp = [];
                erp.time = tl.time;
                erp.label = tl.label;
                erp.dimord = 'chan_time';

                erp.imper100 = squeeze(mean(tl.trial(imper100,:,:)));
                erp.high100 = squeeze(mean(tl.trial(high100,:,:)));
                erp.med80 = squeeze(mean(tl.trial(med80,:,:)));
                erp.low60 = squeeze(mean(tl.trial(low60,:,:)));

                erp.imper_vs_high = erp.imper100 - erp.high100;
                erp.high_vs_low   = erp.high100  - erp.low60;



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


                end


                %% save
                if laplacian toadd1 = '_laplacian'; else toadd1 = ''; end
                if chancluster toadd2 = ''; else toadd2 = '_sinch'; end
                if timebl==1 toadd3 = '_precuebl'; elseif timebl==2 toadd3 = '_prefixbl'; elseif timebl==3 toadd3 = '_preprobebl'; end
                save([param.pathnewlen, 'saved_data/cue_cor_erp', toadd1,toadd2,toadd3,'__' param.subjName], 'erp','ci','ci_allchan');

                clear erp ci ci_allchan
                close all
            end
        end
    end
end