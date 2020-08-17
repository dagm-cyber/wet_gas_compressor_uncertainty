function [InputNames,OutputNames] = InputOutputNames1(IA_raw,FP, Compressor, totCompNames)

totCompNames = strrep(totCompNames,'-','_');

%% input names
InputNames = {IA_raw{2:end,1}}';

    
%% OutPut
CompPerfNames = strcat(repmat({'Compressor_PVTsim_'},length(fields(Compressor.PVTsim)),1),fields(Compressor.PVTsim)); %Names of compressor performance data


VarNames = fields(FP.Sep.PVTsim.Result); %Find names of all variables
VarNames = VarNames(~ismember(VarNames,'outputUnits'));% remove outputUnits
VarNames = VarNames(~ismember(VarNames,'GERG'));% remove GERG
VarNames = VarNames(~ismember(VarNames,'GERG_Mod'));% remove GERG_Mod
FP_Sep_Names = {};%Create empty cell
counter = 1;
for i = 1:length(VarNames) %insert component names 
    if ismember(VarNames(i),{'mixComposition' 'gasComposition' 'oilComposition' 'waterComposition'})
        CompNames = strcat('FP_Sep_PVTsim_',repmat(VarNames(i),length(totCompNames),1),totCompNames');
        FP_Sep_Names = [FP_Sep_Names; CompNames]; %insert comp names in bottom of variable cell array
        counter = counter+length(totCompNames);
    else
        FP_Sep_Names(counter,:) = strcat('FP_Sep_PVTsim_',VarNames(i)); %insert variable name in bottom of array
        counter = counter +1;
    end
end

VarNames = fields(FP.Orifice.PVTsim.Result); %Find names of all variables
VarNames = VarNames(~ismember(VarNames,'outputUnits'));% remove outputUnits
VarNames = VarNames(~ismember(VarNames,'GERG'));% remove GERG
VarNames = VarNames(~ismember(VarNames,'GERG_Mod'));% remove GERG_Mod
FP_Orifice_Names = {};%Create empty cell
counter = 1;
for i = 1:length(VarNames) %insert component names 
    if ismember(VarNames(i),{'mixComposition' 'gasComposition' 'oilComposition' 'waterComposition'})
        CompNames = strcat('FP_Orifice_PVTsim_',repmat(VarNames(i),length(totCompNames),1),totCompNames');
        FP_Orifice_Names = [FP_Orifice_Names; CompNames]; %insert comp names in bottom of variable cell array
        counter = counter+length(totCompNames);
    else
        FP_Orifice_Names(counter,:) = strcat('FP_Orifice_PVTsim_',VarNames(i)); %insert variable name in bottom of array
        counter = counter+1;
    end
    
end

VarNames = fields(FP.In.PVTsim.Result); %Find names of all variables
VarNames = VarNames(~ismember(VarNames,'outputUnits'));% remove outputUnits
VarNames = VarNames(~ismember(VarNames,'GERG'));% remove GERG
VarNames = VarNames(~ismember(VarNames,'GERG_Mod'));% remove GERG_Mod
FP_In_Names = {};%Create empty cell
counter = 1;
for i = 1:length(VarNames) %insert component names 
    if ismember(VarNames(i),{'mixComposition' 'gasComposition' 'oilComposition' 'waterComposition'})
        CompNames = strcat('FP_In_PVTsim_',repmat(VarNames(i),length(totCompNames),1),totCompNames');
        FP_In_Names = [FP_In_Names; CompNames]; %insert comp names in bottom of variable cell array
        counter = counter+length(totCompNames);
    else
        FP_In_Names(counter,:) = strcat('FP_In_PVTsim_',VarNames(i)); %insert variable name in bottom of array
        counter = counter + 1;
    end
end

VarNames = fields(FP.Out.PVTsim.Result); %Find names of all variables
VarNames = VarNames(~ismember(VarNames,'outputUnits'));% remove outputUnits
VarNames = VarNames(~ismember(VarNames,'GERG'));% remove GERG
VarNames = VarNames(~ismember(VarNames,'GERG_Mod'));% remove GERG_Mod
FP_Out_Names = {};%Create empty cell
counter = 1;

for i = 1:length(VarNames) %insert component names 
    if ismember(VarNames(i),{'mixComposition' 'gasComposition' 'oilComposition' 'waterComposition'})
        CompNames = strcat('FP_Out_PVTsim_',repmat(VarNames(i),length(totCompNames),1),totCompNames');
        FP_Out_Names = [FP_Out_Names; CompNames]; %insert comp names in bottom of variable cell array
        counter = counter+length(totCompNames);
    else
        FP_Out_Names(counter,:) = strcat('FP_Out_PVTsim_',VarNames(i)); %insert variable name in bottom of array
        counter = counter + 1;
    end
end

%% Concatinate to one output name array
OutputNames =...
[CompPerfNames;
FP_Sep_Names;
FP_Orifice_Names;
FP_In_Names;
FP_Out_Names];

