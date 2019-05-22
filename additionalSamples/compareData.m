close all
%fileName = 'SSZdata.tsv';
fileName = 'growthData.tsv';


rawData = IO(['output/' fileName]);

colData = rawData(1,:);
data = rawData(2:end,:);

mets = data(:,ismember(colData, 'met'));
metTypes = unique(mets);

conditions = data(:,ismember(colData, 'condition'));
conditionTypes = unique(conditions);

for i = 1:length(metTypes)
    subplot(4,6,i)
    hold all
    curMet = ismember(mets, metTypes{i});
    
    for j = 1:length(conditionTypes)
        curCondition = ismember(conditions, conditionTypes{j});
        curData = data(and(curMet, curCondition),:);
        curTime = cell2nummat(curData(:, ismember(colData, 'time')));
        curValue = cell2nummat(curData(:, ismember(colData, 'value')));
        scatter(curTime, curValue, 'fill')
    end
    title(metTypes{i})
    curYlim = ylim;
    ylim([0 curYlim(2)])
end


condtions = {'hepg2-0mm-mets.txt';
             'hepg2-6mm-mets.txt';
             'hepg2-22mm-mets.txt'};
conditionNames = {'0'; '6'; '22'};
         
for j = 1:length(condtions)
    rawData = IO(['../data/' condtions{j}]);
    mets = rawData(:,1);
    colData = rawData(1,:);
    for i = 1:length(metTypes)
        curMet = ismember(mets, metTypes{i});
        curData = rawData(curMet,:);
        curTime = cell2nummat(curData(:, ismember(colData, 'Time')));
        curValue = cell2nummat(curData(:, ismember(colData, 'Concentration')));       
        subplot(4,6,i)
        hold all
        scatter(curTime, curValue)
    end
end

subplot(4,6,1)
hold all
legend([conditionTypes; conditionNames])
legend boxoff