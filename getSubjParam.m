function param = getSubjParam(pp)

% 4-12-2023 Freek van Ede
% 4-13-2023 Sisi Wang

%% participant-specific notes

%% set path and pp-specific file locations

param.behpath = '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/recue_valid/Data_ana/Behdata_Ana/behdata_new/';
param.eyepath = '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/recue_valid/recuevalid_v1_EEG/Eyedata_raw/rename/';
param.resupath = '/Volumes/sisBNU4Tnew/VUAm_2023/Research_VUAm/recue_valid/Data_ana/Eyedata_Ana/';

% param.behpath = '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/Data_ana/Behdata_Ana/behdata_new/';
% param.eyepath = '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/recuevalid_v1_EEG/Eyedata_raw/rename/';
% param.resupath = '/home/swa100/Documents/sisi_exp/Research_VUAm/recue_valid/Data_ana/Eyedata_Ana/';

sub_list = {'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';};

% pp = 1:length(sub_list);

if pp < 10 % sub_No < 10, add 0 to sub_No
    param.subjName = ['pp0' num2str(pp)];
    param.log =  [param.behpath, ['RecueValid_v1_beh_sub' sub_list{pp} '_allcondi.mat']];
    param.eds1 = [param.eyepath, ['recuevalid_v1_condi_1_sub' sub_list{pp} '.asc']];
    param.eds2 = [param.eyepath, ['recuevalid_v1_condi_2_sub' sub_list{pp} '.asc']];
    param.eds3 = [param.eyepath, ['recuevalid_v1_condi_3_sub' sub_list{pp} '.asc']];
    param.eds4 = [param.eyepath, ['recuevalid_v1_condi_4_sub' sub_list{pp} '.asc']];

else % sub_No >= 10
    param.subjName = ['pp' num2str(pp)];
    param.log =  [param.behpath, ['RecueValid_v1_beh_sub' sub_list{pp} '_allcondi.mat']];
    param.eds1 = [param.eyepath, ['recuevalid_v1_condi_1_sub' sub_list{pp} '.asc']];
    param.eds2 = [param.eyepath, ['recuevalid_v1_condi_2_sub' sub_list{pp} '.asc']];
    param.eds3 = [param.eyepath, ['recuevalid_v1_condi_3_sub' sub_list{pp} '.asc']];
    param.eds4 = [param.eyepath, ['recuevalid_v1_condi_4_sub' sub_list{pp} '.asc']];

end

end


