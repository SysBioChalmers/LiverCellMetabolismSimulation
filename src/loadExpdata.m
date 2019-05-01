function [expData, removalPoints, growth] = loadExpdata(folder, organism, condition)
metaboliteMap = IO([folder '/metaboliteMap.tsv']);
metTranslation = containers.Map(metaboliteMap(:,1),metaboliteMap(:,2));

pathing = [folder '/' organism '-serialized.txt'];
expData = containers.Map;

%Met data
metData = IO(pathing);
colNames = metData(1,:);
conditionCol = ismember(colNames, 'condition');
metCol = ismember(colNames, 'met');
timeCol = ismember(colNames, 'time');
valueCol = ismember(colNames, 'value');
timeAndValue = or(timeCol, valueCol);

for i = 2:size(metData, 1)
    curCondtion = metData{i,conditionCol};
    if strcmp(condition, curCondtion)
        currentMet = metData{i, metCol};
        currentMet = metTranslation(currentMet);
        curVal = cell2nummat(metData(i,timeAndValue))';
        curVal = [curVal; 0]; %for a standard deviation slot
        if not(expData.isKey(currentMet))
            expData(currentMet) = curVal;
        else
            expData(currentMet) = [expData(currentMet) curVal];
        end
    end
end

addStdValues = true;
if addStdValues
    stdFile = importdata([folder '/standardDeviations.txt']);
    stdMets = stdFile.textdata;
    stdData = stdFile.data;
    for i = 1:length(stdMets)
        curData = expData(stdMets{i});
        curData(3,:) = stdData(i);
        expData(stdMets{i}) = curData;
    end
end

%Growth data
pathing = [folder '/' organism '-growthData.txt'];
growthData = IO(pathing);
conditionCol = ismember(colNames, 'condition');
growthData = growthData(2:end,:);

timeCol = ismember(colNames, 'time');
valueCol = ismember(colNames, 'value');
timeAndValue = or(timeCol, valueCol);

matchingCondition = ismember(growthData(:,conditionCol), condition);


growth = cell2nummat(growthData(matchingCondition, timeAndValue));

%Sampling points
removalPoints = importdata([folder '/' organism '-' 'volume.txt']); %sample removal points
end

