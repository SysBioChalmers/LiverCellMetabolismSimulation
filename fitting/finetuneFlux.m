clc
close all
addpath('src')
tolerance = 0.1;

initalX = y(1);
%if subscript dimention missmatch switch these lines
%result = zeros(length(fluxMets), length(mu));
result = zeros(length(fluxMets), length(mu)+1);


sizeOfPlot = ceil(sqrt(length(fluxMets)));

tLim = max(breakPoints);

for i = 1:length(fluxMets)
    fluxMets{i}
    subplot(sizeOfPlot,sizeOfPlot,i);
    curMet = fluxMets{i};
    if isKey(expData, curMet)
        data = expData(curMet)';
        data(data(:,1)>tLim,:)= []; %dont fit beond scope
    else
       data = []; 
    end
    
    title(curMet)
    if size(data,1)>2
        prior = fluxValues(i,:);
        fluxOut = metTune(data, breakPoints', [mu; 0], initalX, volumePoints, prior, tolerance);
    else
        fluxOut = fluxValues(i,:);
    end
    result(i,:) = fluxOut;
    xlabel('h')
    ylabel('mM')
end

for i =1:length(fluxMets)
    curMet = fluxMets{i};
    fprintf('%s', curMet)
    for j = 1:size(result,2)
        fprintf('\t%s', num2str(result(i,j)));        
    end
    fprintf('\n')
end
