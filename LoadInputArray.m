function [IA_raw, IA_num, IA_VarNames] = LoadInputArray(Config)
    %% Load InputArray
    InputArrayFileName = [Config.InputArray '.xlsx'];

    if Config.ReadNewInputArray
        [IA_num,IA_txt,IA_raw] = xlsread(InputArrayFileName); %#ok<ASGLU>
        save('InputArray.mat','IA_num','IA_txt','IA_raw')
    elseif ~Config.ReadNewInputArray
        load('InputArray.mat','IA_raw','IA_num')
    end

   IA_VarNames = IA_raw(2:end,1);
end