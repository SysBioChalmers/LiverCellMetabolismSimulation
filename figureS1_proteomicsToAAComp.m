addpath('proteomicsToAAComp/src')


[sequenceProt, sequences, labels] = loadSequenceData('proteomicsToAAComp/data/AAComp.txt');

% proteomicsData = {
%     'data/platletUniprot.txt'
%     'data/plasmaUniprot.txt'
%     'data/brainUniprot.txt'
%     'data/heartUniprot.txt'
%     'data/liverUniprot.txt'
%     
%     };
proteomicsData = {
    'proteomicsToAAComp/data/Wisniewski/hepg2Uniprot.txt'
    'proteomicsToAAComp/data/Wisniewski/primaryUniprot.txt'
    'proteomicsToAAComp/data/paxDb/Uniprot_Human_liver_Chinese_2010.txt'
    'proteomicsToAAComp/data/paxDb/Uniprot_PA_liver_201308.txt'    
    'proteomicsToAAComp/data/paxDb/Uniprot_liver_Wilhelm_2014_Maxquant.txt'    
    };

%compareDistributions(proteomicsData)

proteomicsLabels = {
        'HEP-G2'
        'Liver'
        'Liver (paxDB 2010)'
        'Liver (paxDB 2013)'        
        'Liver (paxDB 2014)'
    };

proteomicsSymbols = {'s', 'o', '^',  'p', 'v'};

composition = zeros(length(proteomicsData), length(labels));
for i = 1:length(proteomicsData)
    composition(i,:) = getCompositionFromProteomics(proteomicsData{i}, sequenceProt, sequences);
end


samples = randomSample(1000, sequences, 'logNormal');


clf

hold all
meanValue = mean(samples);
[meanValue indx] = sort(meanValue);
sortedSamples = samples(:,indx);
sortedLabels = labels(indx);
sortedComposition = composition(:,indx);

makeBoxplot(sortedSamples, sortedLabels, true)

tickVals = 1:length(labels);
plotHandels = zeros(length(proteomicsLabels),1);
for i = 1:length(proteomicsLabels)
   plotHandels(i) = scatter(sortedComposition(i,:), tickVals, 50,proteomicsSymbols{i},'filled', 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
end

%result = gemAA(sortedLabels);
%hgem = scatter(result, tickVals, 50, 'o', 'filled', 'MarkerFaceColor', [0.4 1 0.4], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
result = humanLiverAA(sortedLabels);
hHuman = scatter(result, tickVals, 50, 'o', 'filled', 'MarkerFaceColor', [1 0.4 0.4], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
%legend([plotHandels; hgem; hHuman], [proteomicsLabels; 'Infant'; 'Human Liver'], 'location', 'se')
legend([plotHandels; hHuman], [proteomicsLabels; 'Human Liver'], 'location', 'se')


xlim([0 0.13])


molWeight = getMolWeight(sortedLabels);
meanAA = sum(meanValue.*molWeight') 

masses = sortedSamples .* repmat(molWeight', size(sortedSamples,1), 1);
sumOfMass = sum(masses,2);
figure()
histogram(sumOfMass)


