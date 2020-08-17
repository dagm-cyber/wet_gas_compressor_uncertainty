function [MCM_UncTable, MCM_UncTableRound] = MonteCarloUncertaintyTable(Compressor,Config)

%% Create uncertainty tables(non round and rounded to N-significant numbers) for the compressor performance parameters

Arrays = struct2cell(Compressor.PVTsim(:)');
STD = struct2cell(Compressor.PVTsimSTD(:)');
STD_Rel = struct2cell(Compressor.PVTsimSTD_Rel(:)');
STD_Perc = struct2cell(Compressor.PVTsimSTD_Perc(:)');
CoverageLower = struct2cell(Compressor.PVTsimLower(:)');
CoverageUpper = struct2cell(Compressor.PVTsimUpper(:)');
Performance_Parameters = fieldnames(Compressor.PVTsim);
CoverageInterval = strcat(repmat('[',length(CoverageLower),1), string(CoverageLower), repmat('-',length(CoverageLower),1),...
    string(CoverageUpper), repmat(']',length(CoverageLower),1));
for i=1:length(Performance_Parameters)
    meanValues(i,1) = mean(Arrays{i});
end

sign = 4; %significant numbers
signSTD = 2;

meanValuesRound = round(meanValues,sign,'significant');
STDRound = round(cell2mat(STD),signSTD,'significant');
STD_RelRound = round(cell2mat(STD_Rel),signSTD,'significant');
STD_PercRound = round(cell2mat(STD_Perc),signSTD,'significant');
CoverageLowerRound = round(cell2mat(CoverageLower),sign,'significant');
CoverageUpperRound = round(cell2mat(CoverageUpper),sign,'significant');
CoverageIntervalRound = strcat(repmat('[',length(CoverageLower),1), string(round(cell2mat(CoverageLower),sign,'significant')),...
    repmat('-',length(CoverageLower),1),string(round(cell2mat(CoverageUpper),sign,'significant')), repmat(']',length(CoverageLower),1));

if Config.UseMeanInputAsFirstIteration
    for i = 1:length(Performance_Parameters)
        NoPertOutput(i) = Arrays{i}(1);
    end
    NoPertOutput = NoPertOutput';
    NoPertOutputRound = round(NoPertOutput,sign,'significant');
    MCM_UncTable = table(Performance_Parameters,NoPertOutput, meanValues,STD, STD_Rel, STD_Perc, CoverageLower, CoverageUpper,CoverageInterval);
    MCM_UncTableRound = table(Performance_Parameters,NoPertOutputRound, meanValuesRound,STDRound, STD_RelRound, STD_PercRound, CoverageLowerRound, CoverageUpperRound,CoverageIntervalRound);
else
    MCM_UncTable = table(Performance_Parameters, meanValues,STD, STD_Rel, STD_Perc, CoverageLower, CoverageUpper, CoverageInterval);
    MCM_UncTableRound = table(Performance_Parameters, meanValuesRound,STDRound, STD_RelRound, STD_PercRound, CoverageLowerRound, CoverageUpperRound,CoverageIntervalRound);
end

