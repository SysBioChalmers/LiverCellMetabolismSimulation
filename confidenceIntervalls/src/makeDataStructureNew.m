function [dataX, dataY, metlabels] = makeDataStructureNew(condition, includeMets)

timePoints = [23.25, 30.25, 47.5];  

poolData = true;

[dataX1, dataY1, mets1] =  extractGrowthData(condition, poolData);
[dataX2, dataY2, mets2] =  extractMetData(condition, timePoints, includeMets);
includeMets = setdiff(includeMets, mets2);
[dataX3, dataY3, mets3] =  extractGlcPyrLacData(condition, timePoints, includeMets);
dataX3(:,2) = dataX3(:,2) + length(mets2);

metlabels = [mets1;mets2;mets3];
dataX = [dataX1;dataX2;dataX3];
dataY = [dataY1;dataY2;dataY3];

end

function [dataX, dataY, mets] =  extractGrowthData(condition, poolData)
%Growth data
    growthData = importdata('growthData.txt');
    dataX = growthData.data;
    
    %remove time 0
    dataX(dataX(:,1)==0,:) = [];
    
    if poolData == true
        %Pool later time points (for 0 mM condition)
        dataX(dataX(:,1)>48,1) = 48;
    else
        dataX(dataX(:,1)>48,:) = [];
    end
    
    %dataX(:,1) = dataX(:,1) - 23.75;
    dataX(:,1) = dataX(:,1)-min(dataX(:,1)); %move to t=0

    if strcmp(condition, '0')
        dataX(dataX(:,2)~=0,:) = [];
    elseif strcmp(condition, '6')
        if poolData == true
            dataX(dataX(:,2)==0,:) = [];
        else
            dataX(dataX(:,2)~=6,:) = [];
        end
    elseif strcmp(condition, '22')
        if poolData == true
            dataX(dataX(:,2)==0,:) = [];
        else
            dataX(dataX(:,2)~=22,:) = [];
        end
    end
    
    dataY = dataX(:,3);
    dataY = dataY/10^6; %unit million cells
    dataX(:,3) = [];
    dataX(:,2) = 1;
    mets = {'growth'};
end

function [dataX, dataY, mets] =  extractMetData(condition, timePoints, includeMets)
    metData = IO('growthDataNew.tsv');
    colNames = metData(1,:);
    
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

function [dataX, dataY, mets] =  extractGlcPyrLacData(condition, timePoints, includeMets)
    metData = IO('../data/hepg2-serialized.txt');
    colNames = metData(1,:);
    
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