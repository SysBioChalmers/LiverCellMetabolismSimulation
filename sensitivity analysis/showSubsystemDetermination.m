function showSubsystemDetermination(model, result, tresh)
minRxns = 4;
upShiftText = 0.06;
barColor = [141 206 255]/256;

hold all
allSubs = model.subSystems(result(:,1));
uniqueSubs = unique(allSubs);

for i = 1:length(uniqueSubs)
    affected = ismember(allSubs, uniqueSubs{i});
    if sum(affected)<minRxns
        allSubs(affected) = {'Other'};
    end
end
uniqueSubs = unique(allSubs);

diff = abs(result(:,4)-result(:,3));
determined = (diff./abs(result(:,2)))<tresh;

analysis = zeros(length(uniqueSubs),2);

for i = 1:length(uniqueSubs)
    affected = determined(ismember(allSubs, uniqueSubs{i}));
    analysis(i,1) = sum(affected);
    analysis(i,2) = length(affected);
end

determinedRate = analysis(:,1)./analysis(:,2);
[determinedRate, idx] = sort(determinedRate);
uniqueSubs = uniqueSubs(idx);
analysis = analysis(idx,:);

barh(determinedRate, 0.9, 'facecolor', barColor, 'edgecolor', 'none')
yticklabels([])

for i =1:length(analysis)
    text(0.02,i+upShiftText, uniqueSubs{i}) 
	text(1.02,i+upShiftText, num2str(analysis(i,2))) 
end

plot([1 1], [0 length(analysis)+0.5], 'k')

xlim([0 1.1])
ylim([0.5 length(analysis)+0.5])

end

