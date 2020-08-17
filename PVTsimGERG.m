function [Output]=PVTsimGERG(Output,Comp,P,T,FlashGERG,FlashInputGERG,FlashOutputGERG,NcompGERG,GERG_idx,GERG_C10pussIdx,i)

%Collect composition to be used in Gerg
if NcompGERG == length(Comp)
    GERG_Comp = Comp;
else
    GERG_Comp = [Comp(GERG_idx );sum(Comp(GERG_C10pussIdx))];
    GERG_Comp = norm(GERG_Comp);
end

%% Performe Gerg flash
FlashOutputGERG = PVTsimFlashCalc(GERG_Comp,P,T,FlashGERG,FlashInputGERG,FlashOutputGERG);

%% Collect Gerg flash properties
[Output.GERG] =  PVTsimFlashCalcProp(Output.GERG,FlashOutputGERG,i);

%% Collect Gerg flash phase compositions
[Output.GERG] =  PVTsimFlashCalcComp(Output.GERG,FlashOutputGERG,NcompGERG,i);

end

%% Normalizatrion function
  function X = norm(X)
    X=X./sum(X);
    end
