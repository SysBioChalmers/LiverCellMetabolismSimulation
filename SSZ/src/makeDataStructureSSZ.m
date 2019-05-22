function [dataX, dataY, metlabels] = makeDataStructureSSZ(condition, includeMets)
timePoints = [24.75 45.25];
[dataX1, dataY1, mets1] =  extractGrowthData(condition);
[dataX2, dataY2, mets2] =  extractMetData(condition, timePoints, includeMets);
includeMets = setdiff(includeMets, mets2);
[dataX3, dataY3, mets3] =  extractGlcPyrLacData(condition, timePoints, includeMets);
dataX3(:,2) = dataX3(:,2) + length(mets2);
 
%move to time 0
dataX1(:,1) = dataX1(:,1) - min(dataX1(:,1));
dataX2(:,1) = dataX2(:,1) - min(dataX2(:,1));
dataX3(:,1) = dataX3(:,1) - min(dataX3(:,1));

metlabels = [mets1;mets2;mets3];
dataX = [dataX1;dataX2;dataX3];
dataY = [dataY1;dataY2;dataY3];

end

function [dataX, dataY, mets] =  extractGrowthData(condition)
    condition = str2num(condition);
    %Growth data
    growthData = importdata('growthSSZ.txt');
    dataX = growthData.data;
    
    %remove later timepoints
    dataX(dataX(:,1)>46,:) = [];
    
    %Remove time 0
    dataX(dataX(:,1)==0,:) = [];
    
    %remove NaN data
    dataX(isnan(dataX(:,3)),:) = [];
    
    %keep data from selected condition
    dataX = dataX(dataX(:,2) == condition,:);
    
    dataY = dataX(:,3);
    dataY = dataY/10^6; %unit million cells
    dataX(:,3) = [];
    dataX(:,2) = 1;
    mets = {'growth'};
end

function [dataX, dataY, mets] =  extractMetData(condition, timePoints, includeMets)
    metData = IO('SSZdata.tsv');
    colNames = metData(1,:);
    
    %Remove mets
    metNames = metData(:,ismember(colNames,'met'));
    metData(not(ismember(metNames, includeMets)),:) = [];
    
    %Extract condition
    con = metData(:,ismember(colNames,'condition'));
    
    %renameCondition at time 0
    filter = contains(metData(:,ismember(colNames,'time')), '0');
    con(filter) = {condition};
    metData(not(ismember(con, condition)),:) = [];
    
    
    %Extract timepoints
    time = cell2nummat(metData(:,ismember(colNames,'time')));
    metData(not(ismember(time,timePoints)),:) = [];
    time(not(ismember(time,timePoints)),:) = [];
    time = time-min(time); %move to t=0
    
    %Extract values
    value = cell2nummat(metData(:,ismember(colNames,'value')));
    value = value/1000;

    %Make met ids
    
    metNames = metData(:,ismember(colNames,'met'));
    [mets, a, b] = unique(metNames);
    metId = b+1;
    
    dataY = value;
    dataX = [time metId];
end

function [dataX, dataY, mets] = extractGlcPyrLacData(condition, timePoints, includeMets)
    metData = IO('carbs.txt');
    colNames = metData(1,:);

    %remove NaN data
    metData(contains(metData(:,5), 'NaN'), :) = []; 
    
    %Remove mets
    metNames = metData(:,ismember(colNames,'met'));
    metData(not(ismember(metNames, includeMets)),:) = [];
    
    %Extract condition
    con = metData(:,ismember(colNames,'condition'));
    metData(not(ismember(con, condition)),:) = [];
    
    %Extract timepoints
    time = cell2nummat(metData(:,ismember(colNames,'time')));
    metData(not(ismember(time,timePoints)),:) = [];
    time(not(ismember(time,timePoints)),:) = [];
    time = time-min(time); %move to t=0
    
    %Extract values
    value = cell2nummat(metData(:,ismember(colNames,'value')));
    value = value/1000;

    %Make met ids
    metNames = metData(:,ismember(colNames,'met'));
    [mets, a, b] = unique(metNames);
    metId = b+1;
    
    dataY = value;
    dataX = [time metId];
end