clc
close all
addpath('src')
global massPerCell
massPerCell = 1000 * 426.8 * 10^-12; %mg per cell

[expData, volumePoints, growthdat] = loadExpdata('../data', 'hepg2', '22mm');

growthdat((end):end,:) = [];

%Identify growth rates:
timeSpan = [10 43];
%growthdat(end,:) = [];
growthRates = expFit(growthdat, timeSpan);

timeSpan = [timeSpan 100];
growthRates = [growthRates; 0];

%%

%growthRates = [0.0124 0.0286 0.0173];
figure()
%Identify fluxes:

[metabolites, crap] = loadFluxes('../fluxvalues', 'hepg2-0mm-.txt');

initalX = growthdat(1,2);
result = zeros(length(metabolites), length(growthRates));

sizeOfPlot = ceil(sqrt(length(metabolites)));

for i = 1:length(metabolites)
    subplot(sizeOfPlot,sizeOfPlot,i);
    curMet = metabolites{i};
    if isKey(expData, curMet)
        data = expData(curMet)';
        data(:,2:3) = data(:,2:3)/1000;
        data(end,:) = [];
    else
       data = []; 
    end

    title(curMet)
    if size(data,1)>2
        fluxes = metFit(data, timeSpan, growthRates, initalX*10^-6, volumePoints);
    else
        fluxes = zeros(1,length(growthRates));
    end
    result(i,:) = fluxes;
    xlabel('h')
    ylabel('mM')
end

for i =1:length(metabolites)
    curMet = metabolites{i};
    fprintf('%s', curMet)
    for j = 1:size(result,2)
        fprintf('\t%f', 10^-3*result(i,j)/massPerCell);        
    end
    fprintf('\n')
end
