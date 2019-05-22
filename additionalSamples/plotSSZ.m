close all
data = 'SSZdata.tsv';

rawData = IO('output/SSZdata.tsv');

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
    if i == 1
        legend(conditionTypes)
        legend boxoff
    end
end
