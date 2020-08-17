function [correlation] = SensitivityPlots(T,SNT,IM,FP,Compressor,InputNames,OutputNames)

%% Choose what to plot
CorrelationsScatterOnOff = true;
PieOnOff = false;
SigmaNormalizedDerivatives =true;
barOutputOnOff = true;
barInputOnOff = false;

%% 
f = waitbar(0,'SensitivityPlots Compressor Performance');


%% Scatter plots incl correlations
if CorrelationsScatterOnOff
waitbar(0.2,f,'Scatter Plots including correlations');
h = findobj('type','figure');
figure(length(h)+1)
tg = uitabgroup; % tabgroup
plotRow = 3;
plotCol = 3;
subplotnumber = plotRow*plotCol;
Y = struct2cell(Compressor.PVTsim(:)');
Y_Names = fieldnames(Compressor.PVTsim);
X = IM;
X_Names = InputNames;
NumberOfPlots = length(X_Names)*length(Y_Names); 
Counter = 1;
X_index = repmat([1:length(X_Names)]',[NumberOfPlots/length(X_Names),1]);

OnesVec = ones(1,length(X_Names))';
Y_index = OnesVec;
for p = 1:length(Y_Names)-1 %Create Y index vector
    nTimesOnes = (p+1)*OnesVec;
    Y_index = [Y_index; nTimesOnes];
end
for k = 1:ceil((length(Y)*length(X(:,1)))/subplotnumber) %loop over tabs
    thistab = uitab(tg,'Title',Y_Names{Y_index(Counter)}(1: min(cellfun('length', Y_Names(Y_index(Counter))),19)));
    axes('Parent',thistab); % somewhere to plot
    for i = 1:subplotnumber %loop over subplots
        subplot(plotRow,plotCol,i)
        try
            correlation(X_index(Counter),Y_index(Counter)) = corr(X(X_index(Counter),:)',Y{Y_index(Counter)}');
            scatter(X(X_index(Counter),:),Y{Y_index(Counter)})
            title([Y_Names{Y_index(Counter)} ' vs ' X_Names{X_index(Counter)} ' correlation= ' mat2str(correlation(X_index(Counter),Y_index(Counter)),4)],'Interpreter', 'none')
            ylabel(Y_Names{Y_index(Counter)},'Interpreter', 'none');
            xlabel(X_Names{X_index(Counter)},'Interpreter', 'none');

            Counter = Counter+1;
            waitbar(Counter/subplotnumber,f,'Scatter Plots including correlations');
        catch
        end
    end
end
end

%% pie plot sensitivity
if PieOnOff
waitbar(0.2,f,'Sensitivity Pie Plot');
h = findobj('type','figure');
figure(length(h)+1)
tg = uitabgroup; % tabgroup
subplotnumber = 4;
varIdx = 0;
vartot =  length(OutputNames);
for k = 1:ceil(length(OutputNames)/subplotnumber)
    thistab = uitab(tg,'Title',OutputNames{varIdx+1}(1:10)); % build iith tab
    axes('Parent',thistab); % somewhere to plot
    for i = 1:subplotnumber
        waitbar(varIdx/vartot,f,'Sensitivity Pie Plot');
        varIdx = varIdx + 1;
        subplot(2,2,i)
        try
            plotData = abs(table2array(T(:,varIdx)));
            plotData = plotData(:,1)./sum(plotData(~isnan(plotData)));
            plotData(isnan(plotData))=0; %replace nan with 0
            newStr = strrep(InputNames,'_',' ');
            pie(plotData,newStr);
            title(['Sensitivity ' OutputNames{varIdx}],'Interpreter', 'none')
            catch
        end
    end
end
end

%% SigmaNormalizedDerivatives bar plot
if SigmaNormalizedDerivatives

waitbar(0.2,f,'SigmaNormalizedDerivatives');
h = findobj('type','figure');
figure(length(h)+1)
tg = uitabgroup; % tabgroup
subplotnumber = 3;
varIdx = 0;
vartot =  length(OutputNames);
for k = 1:ceil(length(OutputNames)/subplotnumber)
    thistab = uitab(tg,'Title',OutputNames{varIdx+1}(1:10)); % build iith tab
    axes('Parent',thistab); % somewhere to plot
    for i = 1:subplotnumber
        waitbar(varIdx/vartot,f,'SigmaNormalizedDerivatives');
        varIdx = varIdx + 1;
        subplot(3,1,i)
        try
            bar(table2array(SNT(:,varIdx)))
            title(['Sigma Normalized Derivatives' OutputNames{varIdx}],'Interpreter', 'none')
            ylabel('Sigma Norm Deriv');
            set(gca,'XTick',1:length(InputNames))
            set(gca,'xticklabel',InputNames)%,'fontsize',8,'FontWeight','bold')
            set(groot, 'DefaultAxesTickLabelInterpreter', 'none') % 'latex' or 'none'
            xtickangle(45)
            grid on
            catch
        end
    end
end

end

%% Relative sensitivity vs all Output vaiables sensitivity bar Plot
    % access column
if barOutputOnOff 
waitbar(0.2,f,'SensitivityPlots One Plot per Output vaiables');
h = findobj('type','figure');
figure(length(h)+1)
tg = uitabgroup; % tabgroup
subplotnumber = 3;
varIdx = 0;
vartot =  length(OutputNames);
for k = 1:ceil(length(OutputNames)/subplotnumber)
    thistab = uitab(tg,'Title',OutputNames{varIdx+1}(1:10)); % build iith tab
    axes('Parent',thistab); % somewhere to plot
    for i = 1:subplotnumber
        waitbar(varIdx/vartot,f,'SensitivityPlots All Output vaiables');
        varIdx = varIdx + 1;
        subplot(3,1,i)
        try
            bar(table2array(T(:,varIdx)))
            title(['Relative sensitivity ' OutputNames{varIdx}],'Interpreter', 'none')
            ylabel('Relative sensitivity');
            set(gca,'XTick',1:length(InputNames))
            set(gca,'xticklabel',InputNames)%,'fontsize',8,'FontWeight','bold')
            set(groot, 'DefaultAxesTickLabelInterpreter', 'none') % 'latex' or 'none'
            xtickangle(45)
            grid on
            catch
        end
    end
end
end

%% Relative sensitivity vs all Input vaiables sensitivity bar Plot
    %access row
if barInputOnOff
    h = findobj('type','figure');
    figure(length(h)+1)
    tg = uitabgroup; % tabgroup
    subplotnumber = 4;
    varIdx = 0;
    vartot =  length(InputNames);
    for k = 1:ceil(length(InputNames)/subplotnumber)
        thistab = uitab(tg,'Title',InputNames{varIdx+1}(1: min(cellfun('length', InputNames(varIdx+1)),19))); % build iith tab
        axes('Parent',thistab); % somewhere to plot
        for i = 1:subplotnumber
            waitbar(varIdx/vartot,f,'SensitivityPlots one plot per Input vaiables');
            varIdx = varIdx + 1;
            subplot(2,2,i)
            try
            bar(table2array(T(varIdx,:)))
            title(['Relative sensitivity ' InputNames{varIdx}],'Interpreter', 'none')
            set(gca,'XTick',1:length(OutputNames))
            set(gca,'xticklabel',OutputNames)
            set(groot, 'DefaultAxesTickLabelInterpreter', 'none')
            xtickangle(45)
            catch
            end
        end
    end
end
    
      close(f)




