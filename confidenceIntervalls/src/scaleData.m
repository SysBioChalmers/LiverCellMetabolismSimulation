function [dataYScaled, scalingfactor] = scaleData(dataX, dataY)
    %asumes growth is nr1
    valueNumbers = unique(dataX(:,2));
    dataYScaled = dataY;
    scalingfactor = zeros(length(valueNumbers),1);
    
    %Set growth as scale 1
    filter = dataX(:,2) == 1;
    growthScale = mean(dataY(filter));
    scalingfactor(1) = 1;
    
    for i = 2:length(valueNumbers)
        filter = dataX(:,2) == i;
        scalingfactor(i) = mean(dataY(filter))/growthScale;
        dataYScaled(filter) = (dataY(filter)/scalingfactor(i));
    end
end


