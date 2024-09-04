function param = getSubjParam_sisi_distraction(pp)

%% participant-specific notes

%% set path and pp-specific file locations
param.path = 'E:\VUAm_2023\Research_VUAm\EyeAna_code\EyeData_ana_FreekNew\';

% experiment 2, with distractors before, during, and after retrocue.
if pp == 1     
    param.subjName = 'pp01';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub46_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub46_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub46_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub46_3.asc']; 
end
if pp == 2    
    param.subjName = 'pp02';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub47_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub47_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub47_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub47_3.asc']; 
end
if pp == 3   
    param.subjName = 'pp03';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub48_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub48_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub48_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub48_3.asc']; 
end
if pp == 4      
    param.subjName = 'pp04';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub49_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub49_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub49_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub49_3.asc']; 
end
if pp == 5      
    param.subjName = 'pp05';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub50_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub50_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub50_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub50_3.asc']; 
end
if pp == 6      
    param.subjName = 'pp06';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub83_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub83_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub83_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub83_3.asc']; 
end
if pp == 7      
    param.subjName = 'pp07';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub84_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub84_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub84_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub84_3.asc']; 
end
if pp == 8     
    param.subjName = 'pp08';   
    param.log =  [param.path, 'rawdata_beh\Orient_retrocue_pv3_varCDI_sd_sub85_block1_part3.mat']; 
    param.eds1 = [param.path, 'rawdata_eye\dis_v3_part1_sub85_1.asc'];  
    param.eds2 = [param.path, 'rawdata_eye\dis_v3_part2_sub85_2.asc']; 
    param.eds3 = [param.path, 'rawdata_eye\dis_v3_part3_sub85_3.asc']; 
end


