halfLife = importdata('proteinHalflife/proteinHalflife.txt');
abundance = importdata('proteinHalflife/abundanceData.txt');
lengthMatrix = importdata('proteinHalflife/lengthMatrix.txt');
lengthMatrix.id = lengthMatrix.textdata(:,1);
lengthMatrix.lngth = sum(lengthMatrix.data,2);

symbolMap = tdfread('proteinHalflife/symbol.txt');
annotation = tdfread('proteinHalflife/annotation.txt');

uniportConversion = importdata('proteinHalflife/ID_conversion_key.txt');
uniprotMap = containers.Map(uniportConversion.textdata(:,4), uniportConversion.textdata(:,1));


halfLifeId = halfLife.textdata(2:end,1);
abundanceId = abundance.textdata(2:end,1);
halfLifeVal = halfLife.data;
abundanceVal = abundance.data;

%remove zeroes
abundanceId(abundanceVal==0) = [];
abundanceVal(abundanceVal==0) = [];

%add and remove duplicates
for i = length(abundanceId):-1:1
    duplicateMatches = find(ismember(abundanceId, abundanceId{i}));
    duplicateMatches(end) = [];
    if not(isempty(duplicateMatches))
        abundanceVal(duplicateMatches(end)) = abundanceVal(duplicateMatches(end)) + abundanceVal(i);
        abundanceVal(i) = [];
        abundanceId(i) = [];
    end    
end

%normalize abundance
totalAbundance = sum(abundanceVal);
abundanceVal = abundanceVal/totalAbundance;

%Make maps
abundanceMap = containers.Map(abundanceId,abundanceVal);
lengthMap = containers.Map(lengthMatrix.id(2:end),lengthMatrix.lngth);
symbolMap = containers.Map(cellstr(symbolMap.Uniprot_ID), cellstr(symbolMap.Gene_Symbol));
annotationMap = containers.Map(cellstr(annotation.Uniprot_ID), cellstr(annotation.Annotation));

%Remove halflife with missing abundance
halfLifeVal = halfLifeVal(ismember(halfLifeId, abundanceId));
halfLifeId = halfLifeId(ismember(halfLifeId, abundanceId));

totalMass = 0;
for i = 1:length(abundanceId)
    if isKey(lengthMap,abundanceId{i})
        totalMass = totalMass + abundanceVal(i).*lengthMap(abundanceId{i});
    end
end

%%
AAConcentration = 5.1;
polymerizationCost = 4.3;

set(gca,'DefaultTextFontSize',20)
for i = 1:length(halfLifeId)
    halfLifeAbudnance(i) = abundanceMap(halfLifeId{i});
    halfLifeLengths(i) = lengthMap(halfLifeId{i});
end
sumOfAbundance = sum(halfLifeAbudnance.*halfLifeLengths);
massCoverage = sumOfAbundance/totalMass;


degradationRate = log(2)./halfLifeVal;

absoluteAbundances = AAConcentration * halfLifeAbudnance .* halfLifeLengths/totalMass; 
synthesisRate = degradationRate' .* absoluteAbundances;
proteinCost = polymerizationCost * synthesisRate;
totalSynthesis = sum(synthesisRate);
measuredCost = sum(proteinCost);
totalCost = measuredCost/massCoverage;

[proteinCost, indx] = sort(proteinCost, 2, 'descend');
matchId = halfLifeId(indx);
sortedDegradation = degradationRate(indx);

costFraction = proteinCost/measuredCost;
cumCost = cumsum(proteinCost);

xVals = (1:length(costFraction))/length(costFraction);
hold all
plot(xVals, cumCost, '-', 'linewidth', 3)

geneToDisp = 8;

set(findall(gcf,'-property','FontSize'),'FontSize',15)
ylim([0 inf])
xlabel('fraction of meassured proteome')
ylabel('expenditure mmol ATP/gdw/h')

plot(xVals(1:geneToDisp), cumCost(1:geneToDisp), 'o', 'linewidth', 2)


for i = 1:geneToDisp
    symbol = symbolMap(matchId{i});
    text(0.005*i + 0.02, cumCost(1) -0.1 + i*0.06, sprintf('%s (%2.1f%%)',symbol, 100*costFraction(i)), 'fontsize', 10)
end

fprintf('Rank\tUniprot\tCost %%\tCost mmol ATP/h\tTurnover rate\tAbundance\tLength\tSymbol\tAnnotation\n')

for i = 1:length(matchId)
    symbol = symbolMap(matchId{i});
    dscrpt = annotationMap(matchId{i});
    costPercent = 100*costFraction(i);
    A = abundanceMap(matchId{i});
    L = lengthMap(matchId{i});
    D = sortedDegradation(i);
    fprintf('%i\t%s\t%2.5f\t%2.5f\t%2.5f\t%2.5f\t%2.5f\t%s\t%s\n', i, matchId{i}, costPercent, proteinCost(i), D, A, L, symbol, dscrpt)
end



coverage = sum(halfLifeAbudnance);


%%
%Statistics
figure()
celltype = 'hepg2';
conditions = 22;
conditionName = 'hepG2 22mM';
fluxProfile = 2;

GAM = 48;
M = 1;
model = setupBiomass(model, GAM, M);

flux = zeros(length(model.rxns),length(conditions));

if and(strcmp(celltype, 'hepg2'), fluxProfile == 2)
    fileName = ['confidenceIntervalls\output\hepg2-' num2str(conditions) '.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    fluxIn = cell2nummat(raw(2:end,2))/1000;
else
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', celltype, conditions);
    fluxIn = fluxValues(:,fluxProfile)/1000;
end

model = bindFBA(model, fluxMets, fluxIn);

solution = solveLinMin(model,1);
rxnsWithFlux = solution.x>0;

allGenes = model.genes';
fluxGenes = getInvolvedGenes(model, rxnsWithFlux);

totalInModel = 0;
totalInFluxModel = 0;
totalEC = 0;
for i = 1:length(matchId)
    geneId = uniprotMap(matchId{i});
    if ismember(geneId, allGenes)
        totalInModel = totalInModel +  proteinCost(i);
        if ismember(geneId, fluxGenes)
            totalInFluxModel = totalInFluxModel +  proteinCost(i);
        end
    else
        dscrpt = annotationMap(matchId{i});
        containsEC = strfind(annotationMap(matchId{i}), '(EC ');
        if not(isempty(containsEC))
            totalEC = totalEC + proteinCost(i);
        end       
    end
end

notQuantified = totalCost-measuredCost;
notInModel = measuredCost-totalInModel;
inmodelButNoFlux = totalInModel-totalInFluxModel;

fprintf('Not quantified\t%2.2f\n', notQuantified)
fprintf('Not in model\t%2.2f\n', notInModel)
fprintf('In model without flux\t%2.2f\n', inmodelButNoFlux)
fprintf('With flux\t%2.2f\n', totalInFluxModel)
%fprintf('Has EC\t%2.2f\n', totalEC)

data = [notInModel, nan;
        inmodelButNoFlux, nan;
        totalInFluxModel, nan;
        notQuantified, nan];

labels = {'Not in model'
          'In model without flux'
          'With flux'
          'Not quantified'
          };
      
    
bar(data','stacked', 'linestyle', 'none')
legend(labels)
legend boxoff
ylabel('ATP expenditure [mmol ATP/gdw/h]')
ylim([0 4])
xlim([0.25 1.75])
