function results = meanAndError(dataX,dataY)
timepoints = unique(dataX);
results = zeros(length(timepoints),3);
    for i = 1:length(timepoints)
       filter = dataX == timepoints(i);
       results(i,1) = timepoints(i);
       results(i,2) = nanmean(dataY(filter));
       results(i,3) = nanstd(dataY(filter));
    end
end

