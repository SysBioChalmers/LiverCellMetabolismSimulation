function overlayData(dict, conditionType, width)
hold all
rawData = IO('additionalSamples/output/growthData.tsv');

reverseMap = IO('data/metaboliteMap.tsv');
metTranslation = containers.Map(reverseMap(:,2),reverseMap(:,1));

height = ceil(length(dict)/width);

colData = rawData(1,:);
data = rawData(2:end,:);

mets = data(:,ismember(colData, 'met'));

conditions = data(:,ismember(colData, 'condition'));
curCondition = ismember(conditions, conditionType);

for i = 1:length(dict)
    subplot(height,width,i)
    hold all
    curMetName = metTranslation(dict{i});
    curMet = ismember(mets, curMetName);
    curData = data(and(curMet, curCondition),:);
    curTime = cell2nummat(curData(:, ismember(colData, 'time')));
    curValue = cell2nummat(curData(:, ismember(colData, 'value')));
    scatter(curTime, curValue, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',[0 0 0])
    curYlim = ylim;
    if curYlim(2) < max(curValue)
        ylim([0 max(curValue)])
    end
end


end