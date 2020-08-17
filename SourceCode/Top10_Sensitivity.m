function [Top10Sens, Top10SensTable] = Top10_Sensitivity(SigmaNormSensitivityMatrix_Table,InputNames);
%% Collect top 10 sensitivities

Comp = table2cell(SigmaNormSensitivityMatrix_Table(:,1:23));
Comp = cell2mat(Comp);

Fields = fields(SigmaNormSensitivityMatrix_Table);
CompressorFields = Fields(1:23);

for  i = 1:length(Comp(1,:))
    B = Comp(:,i);
    [B_sort,I] = sort(abs(B),'descend');
    NewNameList = InputNames(I);
    NewVar = B(I);
    Top10Sens(i).Top10Var = NewVar(1:10);
    Top10Sens(i).Top10VarNames = {NewNameList{1:10}}';
    Top10Sens(i).PerformanceParameter = CompressorFields(i);
end

Top10SensTableCell = {};

for i = 1:length(Comp(1,:))
    
    Top10SensTableCell = [Top10SensTableCell {Top10Sens(i).Top10VarNames{:}}' num2cell(Top10Sens(i).Top10Var(:))]; 
end
RepNames = {};
for i = 1:length(Comp(1,:))
    RepNames = [RepNames CompressorFields(i) CompressorFields(i)];
end

Top10SensTableCell = [RepNames; Top10SensTableCell];
Top10SensTable = cell2table(Top10SensTableCell);


