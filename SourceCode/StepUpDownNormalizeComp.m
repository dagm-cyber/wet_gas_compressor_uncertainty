function [IM] = StepUpDownNormalizeComp(IM, Config)


%Find gas and oil comp
Gas_X = IM(Config.GasCompRange,:);
Oil_X = IM(Config.OilCompRange,:);
    
% Normalize
if any(Gas_X(:)) %No normalization if no gas composition i present
    Gas_X = Gas_X./sum(Gas_X);
end
if any(Oil_X(:)) %No normalization if no oil composition i present
    Oil_X = Oil_X./sum(Oil_X);
end
%Replace gas and oil comp with norm compositions
IM(Config.GasCompRange,:) = 100*Gas_X;
IM(Config.OilCompRange,:) = 100*Oil_X;