function [] = MonteCarloPlots(Compressor,IM,FP,LoopCalcData,Config,totCompNames)


%% Compressor performance spread histogram

h = findobj('type','figure');
figure(length(h)+1)
tg = uitabgroup; % tabgroup
plotRow = 3;
plotCol = 3;
subplotnumber = plotRow*plotCol;
var = struct2cell(Compressor.PVTsim(:)');
varNames = fieldnames(Compressor.PVTsim);
Counter = 0;
for k = 1:ceil(length(var)/subplotnumber)
    thistab = uitab(tg,'Title',varNames{Counter+1}(1: min(cellfun('length', varNames(Counter+1)),19)));
    axes('Parent',thistab); % somewhere to plot
    for i = 1:subplotnumber
        Counter = Counter + 1;
        subplot(plotRow,plotCol,i)
        try
            histogram(var{Counter})
            title(varNames{Counter},'Interpreter', 'none')
            ylabel('Count');
%             set(gca,'XTick',1:length(InputNames))
%             set(gca,'xticklabel',InputNames)%,'fontsize',8,'FontWeight','bold')
%             set(groot, 'DefaultAxesTickLabelInterpreter', 'none') % 'latex' or 'none'
%             xtickangle(45)
            %grid on
            catch
        end
    end
end
    

%% convergence plot

var = struct2cell(Compressor.PVTsim(:)');
for i = 1:length(fields(Compressor.PVTsim))
    CumMean(i,:) = CumulativeMean(var{i});
end    
h = findobj('type','figure');
figure(length(h)+1)
tg = uitabgroup; % tabgroup
plotRow = 3;
plotCol = 3;
subplotnumber = plotRow*plotCol;
varNames = fieldnames(Compressor.PVTsim);
Counter = 0;

for k = 1:ceil(length(CumMean(:,1))/subplotnumber)
    thistab = uitab(tg,'Title',varNames{Counter+1}(1: min(cellfun('length', varNames(Counter+1)),19)));
    axes('Parent',thistab); % somewhere to plot
    for i = 1:subplotnumber
        Counter = Counter + 1;
        subplot(plotRow,plotCol,i)
        try
            plot(CumMean(Counter,:))
            title([varNames{Counter} ' convergence'],'Interpreter', 'none')
            xlabel('Number of Calculations');
            catch
        end
    end
end
 

%% STD of comperssor performance parameters

h = findobj('type','figure');
figure(length(h)+1)
var = struct2array(Compressor.PVTsimSTD_Perc(:)');
varNames = fieldnames(Compressor.PVTsimSTD_Perc);
bar(var)
ylabel('StandardDeviation in Percent');
somenames = varNames;
set(gca,'XTick',1:length(somenames))
set(gca,'xticklabel',somenames)
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
xtickangle(45)

%% Code Performance to se if code slows down due to memory usage or other

h = findobj('type','figure');
figure(length(h)+1)
scatter(1:1:length(LoopCalcData.AccululativeIterationTime),LoopCalcData.AccululativeIterationTime)
title('Code Performance')
ylabel('AccululativeIterationTime [s]')
xlabel('Iteration number [-]')

%% DEVIATION EXPECTED VALUE VS MONTE CARLO MEAN

if Config.UseMeanInputAsFirstIteration
h = findobj('type','figure');
figure(length(h)+1)

var = struct2cell(Compressor.PVTsim(:)');
var = cell2mat(var);
varNames = fieldnames(Compressor.PVTsim);
varMean = mean(var,2);
varExp =  var(:,1);
somenames = varNames;
deviationsFromExpectedValues = 100*(varMean./varExp-1)';
bar(deviationsFromExpectedValues)
    ylabel('Deviation in Percent');
    title('DEVIATION FROM MEAN EXPECTED INPUT VALUES')
    set(gca,'XTick',1:length(somenames))    
    set(gca,'xticklabel',somenames)
    set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
    xtickangle(45)
    
end

%% Total Composition uncertainty
h = findobj('type','figure');
figure(length(h)+1)

LoadComp = cell2mat(FP.Sep.PVTsim.Result.mixComposition(:)');
LoadCompNames = totCompNames;

LoadCompAvg = mean(LoadComp')';
LoadCompStd = std(LoadComp')';
STDPercent = 100*LoadCompStd./LoadCompAvg;
bar(STDPercent)
title('Mix Composition')
ylabel('Standard deviation in Percent');
set(gca,'XTick',1:length(LoadCompAvg))
set(gca,'xticklabel',LoadCompNames)
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
xtickangle(45)


%% hp vs flow coeff
h = findobj('type','figure');
figure(length(h)+1)
s = scatter(Compressor.PVTsim.FlowCoeffDry, Compressor.PVTsim.PolyHead);
 xlim([0.025 0.075])
 ylim([10 40])
title('Polytropic Head')
ylabel(['Polytropic Head' ' [kJ/kg]']);
xlabel(['Flow Coff' ' [-]']);
s.Marker = 'o';
s.SizeData = 10;
transp = 0.1;
markColor = 'b';
s.MarkerFaceColor = markColor;
s.MarkerFaceAlpha = transp;
s.MarkerEdgeColor = markColor;
s.MarkerEdgeAlpha = transp;


 
%% hp vs flow
h = findobj('type','figure');
figure(length(h)+1)
s = scatter(3600*FP.In.PVTsim.Result.totalMassFlow./FP.In.PVTsim.Result.mixDensity, Compressor.PVTsim.PolyHead);
 xlim([0 4000])
 ylim([10 40])
title('Polytropic Head')
ylabel(['Polytropic Head' ' [kJ/kg]']);
xlabel(['Flow' ' [$\frac{Am^3}{h}$]'], 'Interpreter','latex');%'$\hat{\psi}$','Interpreter','latex'
s.Marker = 'o';
s.SizeData = 10;
transp = 0.1;
markColor = 'b';
s.MarkerFaceColor = markColor;
s.MarkerFaceAlpha = transp;
s.MarkerEdgeColor = markColor;
s.MarkerEdgeAlpha = transp;

%% hp vs flow
h = findobj('type','figure');
figure(length(h)+1)
s = scatter(3600*FP.In.PVTsim.Result.totalMassFlow./FP.In.PVTsim.Result.mixDensity, Compressor.PVTsim.PolyEff_thermo);
 xlim([0 4000])
 ylim([0 1])
title('Polytropic Efficiency Thermodynamic')
ylabel(['Polytropic Efficiency Thermodynamic' ' [-]']);
xlabel(['Flow' ' [$\frac{Am^3}{h}$]'], 'Interpreter','latex');%'$\hat{\psi}$','Interpreter','latex'
s.Marker = 'o';
s.SizeData = 10;
transp = 0.1;
markColor = 'b';
s.MarkerFaceColor = markColor;
s.MarkerFaceAlpha = transp;
s.MarkerEdgeColor = markColor;
s.MarkerEdgeAlpha = transp;


