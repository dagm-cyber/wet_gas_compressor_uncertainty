function [OutMatrix] = OutMatrix(FP,Compressor)

%% Compressor Performance

Comp = struct2cell(Compressor.PVTsim);
Comp = cell2mat(Comp);

%% Loop calc flash points  


A(1).FP = struct2cell(FP.Sep.PVTsim.Result);
A(2).FP = struct2cell(FP.Orifice.PVTsim.Result);
A(3).FP = struct2cell(FP.In.PVTsim.Result);
A(4).FP = struct2cell(FP.Out.PVTsim.Result);


%create Matrix for Flash points 
Array = [];
for k = 1:length(A)
    for i = 1:length(A(k).FP)
        if istable(A(k).FP{i})
           Array = Array;  
        elseif isdouble(A(k).FP{i})
            Array = [Array; A(k).FP{i}];
        elseif iscell(A(k).FP{i})
             Array = [Array; cell2mat(A(k).FP{i})];
        %elseif isstruct(A(k).FP{i})
         %    Array = [Array; cell2mat(A(k).FP{i})];
        end
    end
end
Loop = Array;
OutMatrix = [Comp; Loop];

