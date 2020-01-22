selectedCondition = [0 6 22]; %0 6 22
selectedTime = [23.25; 30.25; 47.5]; %23.25 30.25 47.5

addpath('../src')
serumData = IO('serumAA.txt');
mets = strrep(serumData(:,1), 'L-', '');
mets = strrep(mets, 'Isoleucine', 'Iso-leucine');
values = str2double(serumData(:,2:3));

cultureData = IO('../confidenceIntervalls/growthDataNew.tsv');
cultureMets = cultureData(2:end,1);
cultureTime = str2double(cultureData(2:end,3));
cultureCondition = str2double(cultureData(2:end,4));
cultureValues = str2double(cultureData(2:end,6));


conditionFilter = ismember(cultureCondition, selectedCondition);
timeFilter = ismember(cultureTime, selectedTime);
combined = and(conditionFilter, timeFilter);

compareValues = zeros(size(values,1), 2);
for i = 1:size(compareValues, 1)
    metFilter = ismember(cultureMets, mets{i});
    data =cultureValues(and(combined, metFilter));
    compareValues(i,1) = mean(data);
    compareValues(i,2) = std(data);
end
hold all
X = compareValues(:,1);
Y = values(:,1);
XN = -compareValues(:,2);
XP = compareValues(:,2);
YN = -values(:,2);
YP = values(:,2);

plot([0 1000], [0 1000], 'k-')
errorbar(X, Y, YN, YP, XN, XP, '.k')
scatter(X, Y, 'filled')
text(X+10, Y+15, mets)

xlabel('Culture [uM]')
ylabel('Blood [uM]')
axis equal
xlim([0 1000])
ylim([0 1000])
