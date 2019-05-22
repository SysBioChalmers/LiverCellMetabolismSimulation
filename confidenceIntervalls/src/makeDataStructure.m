function [dataX, dataY, metlabels] = makeDataStructure(condition)

[dataX1, dataY1] =  extractGrowthData(condition);
[dataX2, dataY2, mets] =  extractMetData(condition);

metlabels = [{'growth'};mets];
dataX = [dataX1;dataX2];
dataY = [dataY1;dataY2];

end

function [dataX, dataY] =  extractGrowthData(condition)
%Growth data
    growthData = importdata('growthData.txt');
    dataX = growthData.data;
    dataX(dataX(:,1)==0,:) = [];
    dataX(dataX(:,2)>0,2) = 1;
    
    dataX(:,1) = dataX(:,1)-min(dataX(:,1)); %move to t=0

    if strcmp(condition, '0mM')
        dataX(dataX(:,2)==1,:) = [];
    else
        dataX(dataX(:,2)==0,:) = [];
    end
    dataY = dataX(:,3);
    dataY = dataY/10^6;
    dataX(:,3) = [];
    dataX(:,2) = 1;
    
end

function [dataX, dataY, mets] =  extractMetData(condition)
    metData = importdata('metData.txt');
    data = metData.data;
    data(:,1) = data(:,1)-min(data(:,1)); %move to t=0
    
    mets = metData.textdata(2:end,1);
    [mets, a, b] = unique(mets);
    data = [data b+1];
    wells = [4 7 10 13 16];
    if strcmp(condition, '6mM')
        wells = wells + 1;
    elseif strcmp(condition, '22mM')
        wells = wells + 2;
    end
    data(not(ismember(data(:,2),wells)),:) = [];
    
    dataY = data(:,3);
    dataX = data(:,[1 4]);
end

