function [SensitivityMatrix, SensitivityRelativeMatrix, SigmaNormSensitivityMatrix, totCompNames] = sensitivityPVTsim(IA,IM_Up,IM_Down,IA_raw,Config,Compressor,FP,IM)
%% Calculate sensitivity by step up and step down of all input parameters and calculate derivatives of all output parameters

%% no change
[FP_IA, LoopCalcData_IA,IS_IA,Constants,totCompNames] = PVTsimLoopCalc_IsolateF(IA,Config,'Sensitivity No Change progress',true);
[Compressor_noChange] = CompressorPerformanceEquationNode(FP_IA,IS_IA,Constants); %calculate compressor performance parameters
[OutMatrix_NoChange] = OutMatrix(FP_IA,Compressor_noChange); %Create OutMatrix for the no change case
%% Step up 
[FP_Up, LoopCalcData_UP,IS_Up] = PVTsimLoopCalc_IsolateF(IM_Up,Config,'Sensitivity Step Up progress',true); %loop through the step up input matrix 
[Compressor_StepUp] = CompressorPerformanceEquationNode(FP_Up,IS_Up,Constants); %calculate compressor performance parameters
[OutMatrix_Up] = OutMatrix(FP_Up,Compressor_StepUp); %Create OutMatrix for the step up case

%% step down
[FP_Down, LoopCalcData_Down,IS_Down] = PVTsimLoopCalc_IsolateF(IM_Down,Config,'Sensitivity Step Down progress',true); %loop through the step down input matrix 
[Compressor_StepDown] = CompressorPerformanceEquationNode(FP_Down,IS_Down,Constants); %calculate compressor performance parameters
[OutMatrix_Down] = OutMatrix(FP_Down, Compressor_StepDown);%Create OutMatrix for the step down case

%% Delta matixes
DeltaIn = diag(IM_Up - IM_Down);
DeltaOut = OutMatrix_Up - OutMatrix_Down;

%% Sensitivity matrix and relative sensitivity matrix
SensitivityMatrix = (DeltaOut./DeltaIn')';
SensitivityMatrix(isnan(SensitivityMatrix)) = 0; %set sensitivity to 0 for nan that has arrisen from DeltIn = 0, or from nonexisting output parameters (exsample wet Re_M for dry case)
SensitivityRelativeMatrix = (IA'./OutMatrix_NoChange)'.*SensitivityMatrix;

%% calculate sigma normalized sensitivity
[OutMatrixAll] = OutMatrix(FP, Compressor);
sigmaY = std(OutMatrixAll,0,2);
sigmaX = std(IM,0,2);
sX_div_SY = sigmaX*(sigmaY.^(-1))';
SigmaNormSensitivityMatrix = sX_div_SY.*SensitivityMatrix;



