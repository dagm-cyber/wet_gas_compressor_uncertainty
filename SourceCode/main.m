
%clear all
%close all
%clfall

%% Config
tic
Config.EOS = 'PR_HV'; %'SRK_HV' or 'SRK_CPA' or 'PR_HV' or 'PR_CPA' or 'PR_Classic' or SRK_Classic %select EOS and polar components model  
Config.InputArray = 'InputArrayWet'; %select 'InputArrayWet' for wet or 'InputArrayDry' for dry case
Config.UseConstLoadComposition = false; % if true --> simulate 0 uncertainty in composition
Config.SensitivityOnOff = true; %false --> turn off sensitivity analysis 
Config.ReadNewInputArray = true; %false can be used to speed up calculation if no change has been done to input array and it is already loaded into memory
Config.UseMeanInputAsFirstIteration = true; %true to use expectation values as first input to monte carlo simulations. OK for large number of iterations.
Config.UseErrorInOrificeDischargeCoeff = true; %true to include unsertainty in orifice discharge coefficient, except for input mean calculation
Config.UseAfterNormComp = true; %Use composition after normalization in sensitivity analysis if true 
Config.GasCompRange = 1:14; %location of gas composition in input array
Config.OilCompRange = 25:38; %location of oil composition in input array
Config.LastFluidIndex = 48; %location of last compositional data in input array

rng('default'); %default seed, the same random numbers are produced as if you restarted MATLAB
n = 10000; %number of perturbations
h =  findobj('type','figure');

%% Load InputArray
[IA_raw, IA_num, IA_VarNames] = LoadInputArray(Config);
IA = IA_num(:,2); 

%% Perturbations
[IM] = Perturbations(n, IA_num,Config);

%% Monte Carlo PVTsim
[FP, LoopCalcData,IS,Constants, totCompNames] = PVTsimLoopCalc_IsolateF(IM,Config,'MonteCarlo Progresss',false);

%% Compressor Performace Equation Node 
[Compressor] = CompressorPerformanceEquationNode(FP,IS,Constants);
[Compressor] = CompressorPerformanceStatistics(Compressor);

%% MonteCarlo Plots
MonteCarloPlots(Compressor,IM,FP,LoopCalcData,Config,totCompNames)

%% Monte Carlo uncertainty table
[MCM_UncTable, MCM_UncTableRound] = MonteCarloUncertaintyTable(Compressor, Config); %Create uncertainty table
writetable(MCM_UncTable,['MCM_UncTable_' Config.EOS '_' Config.InputArray '_n=' num2str(n) '.csv']) %Write uncertaity table to file
writetable(MCM_UncTableRound,['MCM_UncTableRound_' Config.EOS '_' Config.InputArray '_n=' num2str(n) '.csv'])%Write rounded uncertaity table to file
MCM_UncTableRoundExtract = MCM_UncTableRound(:,{'Performance_Parameters','meanValuesRound','STDRound','STD_PercRound','CoverageIntervalRound'}); %Extract main values to uncertainty table
writetable(MCM_UncTableRoundExtract,['MCM_UncTableRoundExtract_' Config.EOS '_' Config.InputArray '_n=' num2str(n) '.csv'])%Write exstract uncertaity table to file


%% sensitivity PVTsim
if Config.SensitivityOnOff %no sensitivity analysis if false
%% Step up and step down all inputs given in the input array, and generate step up and step down matrix
    [IM_Up] = StepUpDown(IA, 'Up', 0.5);
    [IM_Down] = StepUpDown(IA, 'Down', 0.5);

%% Use composition after normalization in sensitivity analysis if true 
if Config.UseAfterNormComp 
    [IM_Up] = StepUpDownNormalizeComp(IM_Up, Config);
    [IM_Down] = StepUpDownNormalizeComp(IM_Down, Config);
end

%% Calculate sensitivity
[SensitivityMatrix, SensitivityRelativeMatrix, SigmaNormSensitivityMatrix, totCompNames] = ...
    sensitivityPVTsim(IA,IM_Up,IM_Down,IA_raw,Config,Compressor,FP,IM);

%% Input Output variable names
[InputNames,OutputNames] = InputOutputNames1(IA_raw,FP,Compressor,totCompNames); %collect input and output variables names into cell arrays
%% sensitivity table
SensitivityRelative_Table = array2table(SensitivityRelativeMatrix,'RowNames',InputNames','VariableNames',OutputNames');%create sensitivity table
save('AppSensitivity','SensitivityRelative_Table','InputNames','OutputNames') %save sensitivity table, inputnames and outputnames as .mat file

SigmaNormSensitivityMatrix_Table = array2table(SigmaNormSensitivityMatrix,'RowNames',InputNames','VariableNames',OutputNames');%create Sigma normalized sensitivity matrix table
save('SigmaNormSensitivityMatrix_Table','InputNames','OutputNames') %save Sigma normalized sensitivity matrix table, inputnames and outputnames as .mat file

writetable(SigmaNormSensitivityMatrix_Table,['SigmaNormSensitivityMatrix_Table_' Config.EOS '_' Config.InputArray '_n=' num2str(n) '.csv'],'WriteRowNames',true)%Create .CSV of the SigmaNormSensitivityMatrix 

%% Sensitivity Plots
[correlation_CompressorVsInput] = SensitivityPlots(SensitivityRelative_Table,SigmaNormSensitivityMatrix_Table,IM,FP,Compressor,InputNames,OutputNames); %create sensitivity plots

%% Correlation Tables
correlation_CompressorVsInput_Table = array2table(correlation_CompressorVsInput,'RowNames',InputNames','VariableNames',...
    fieldnames(Compressor.PVTsim)'); %create correlation table of compressor performance parameters compared with input variables
writetable(correlation_CompressorVsInput_Table,['correlation_CompressorVsInput_Table_' Config.EOS '_'...
    Config.InputArray '_n=' num2str(n) '.csv',],'WriteRowNames',true); %create .csv of correlations

%%Top 10 sensitivity
[Top10Sens, Top10SensTable] = Top10_Sensitivity(SigmaNormSensitivityMatrix_Table,InputNames); % Collect top 10 sensitivities
writetable(Top10SensTable,['Top10SigmaNormalizedSensitivityTable' Config.EOS '_'...
    Config.InputArray '_n=' num2str(n) '.csv']); %create .csv of top 10 sensitivities

end
toc


