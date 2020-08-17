function [StepMatrix] = StepUpDown(IA, StepType, StepPercent)

if ~(isequal(StepType,'Up') || isequal(StepType,'Down'))
    error('StepType must be eighter Up or Down')
end


if isequal(StepType,'Down')
    StepPercent = -StepPercent;
end


IA_Rep = repmat(IA,1,length(IA)); %Repeat input array (crate matrix)

NewDiag = IA.*(1+StepPercent/100); %Create new diagonal

DiagRemoved = IA_Rep - diag(diag(IA_Rep)); %Remover diagonal

StepMatrix = DiagRemoved + diag(NewDiag); %Insert new diagonal