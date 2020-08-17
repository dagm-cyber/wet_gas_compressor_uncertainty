function [CumMean] = CumulativeMean(array)

wascolumn = false;

if iscolumn(array)
    wascolumn = true;
    array = array';
end


CumMean = zeros(1,length(array));
    for i = 1:length(array)
        CumMean(i) = mean(array(1:i));
    end

if wascolumn
    CumMean = CumMean';
end