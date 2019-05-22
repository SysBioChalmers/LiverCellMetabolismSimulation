clc
load('model/genericHuman2')
addpath('src')

close all
celltype = {
    'hepG2'
    'hepG2' 
    %'huh7'
    'hepG2'    
    %'huh7'
    };

conditions = {
    '0'
    '6'
    %22
    '22'
    %0
    };

conditionName = {
    'hepG2 0mM'
    'hepG2 6mM'
    %'huh7 22mM'
    'hepG2 22mM'  
    %'huh7 0mM'
    };

fluxProfile = [
    2
    2
    %2
    2
    %2
    ];

allGrowth = {};
for i = 1:length(conditions)
    [expData, volumePoints, growthdat] = loadExpdata('data', celltype{i}, conditions{i});
    allGrowth{i} = growthdat;
end    

%%
plotColors = get(gca,'ColorOrder');
hold all


assumedDepletion = 48;
plot(assumedDepletion * [1 1], [0 4*10^6], 'k--', 'HandleVisibility','off')

for i = 1:length(allGrowth)
    curGrowth = allGrowth{i};
    curGrowth(curGrowth(:,1)>100,:) = [];
    scatter(curGrowth(:,1), curGrowth(:,2), 'fill', 'markerfacecolor', plotColors(i,:))
end

tvals1 = linspace(0, assumedDepletion);
tvals2 = linspace(assumedDepletion, 80);
results = zeros(length(conditions),2);
for i = 1:length(allGrowth)
    curGrowth = allGrowth{i};
    curGrowth(curGrowth(:,1)>100,:) = [];    
    filter1 = curGrowth(:,1) <= assumedDepletion;
    KM = polyfit(curGrowth(filter1,1), log(curGrowth(filter1,2)), 1);
    plot(tvals1, exp(KM(2))*exp(tvals1*KM(1)), 'color', plotColors(i,:), 'linewidth', 3);
    
    cellAtDepletion = KM(2) + assumedDepletion*KM(1);
    filter2 = curGrowth(:,1) > assumedDepletion;
    mdlL = fitlm(curGrowth(filter2,1)-assumedDepletion, log(curGrowth(filter2,2))-cellAtDepletion,'Intercept',false);
    slope = mdlL.Coefficients.Estimate(1);
    if slope<0
       slope = 0; 
    end
    plot(tvals2, exp(cellAtDepletion) * exp((tvals2-assumedDepletion)*slope), 'color', plotColors(i,:), 'linewidth', 3);
    results(i,1) = KM(1);
    results(i,2) = slope;
end


legend(conditionName, 'location', 'nw')
legend boxoff
xlabel('time [h]')
ylabel('cells')
xlim([0 80])

results