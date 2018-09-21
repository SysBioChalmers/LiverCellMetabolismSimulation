clc
load('genericHuman2')
addpath('src')
[id, exchangeRxns] =getExchangeRxns(model);
model = configureSMatrix(model, 90, 'HumanGrowth', 'human_growthMaintainance[c]');
model = configureSMatrix(model, 6, 'HumanGrowth', 'human_protein_pool[c]');
model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_RNAPool[c]');
model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_DNAPool[c]');
model = configureSMatrix(model, 1, 'HumanGrowth', 'glycogen[c]');
model = configureSMatrix(model, 0, 'lipidPool', 'fattyAcidPool[c]');
model = configureSMatrix(model, 0.3, 'HumanGrowth', 'lipidPool[c]');
%FA composition http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0074283


%https://www.google.com/patents/US20060014223 1.2 µg/mL culture
%Product should be 1000 * 10^-6 * 120 * 4 * 10^-3 * 1000 ~ 0.5 mg

%Set Non-GAM
model = setParam(model, 'lb', 'human_ATPMaintainance', 1);
[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');
nrDists = size(fluxValues,2)-1;

allSolutions = zeros(length(model.rxns), nrDists);
mu = zeros(1,nrDists);
for i = 2:nrDists
    model = bindFBA(model, fluxMets, fluxValues(:,i)/1000);
    objectiveFunction = 'HumanGrowth';
    model = setParam(model, 'obj', objectiveFunction, 1);
    model = setParam(model, 'ub', objectiveFunction, 0.05);
    solution = solveLinMin(model);
    allSolutions(:,i) = solution.x;
    mu(i) = solution.x(findIndex(model.rxns,objectiveFunction));
end
tresh = 10^-6;
results = allSolutions(exchangeRxns,:);

exchange=sum(abs(results),2)>tresh;
results = results(exchange,:);
[res, id] = sort(mean(results,2));
results = results(id,:);

eq = constructEquations(model, exchangeRxns(exchange));
eq = eq(id);
for i = 1:length(eq)
    
   fprintf('%s', eq{i}) 
   for j = 1:length(results(i,:))
       fprintf('\t%f', results(i,j)) 
   end
   fprintf('\n') 
end

%%
close all
fluxMets = {'O2[s]', 'H+[s]',  'CO2[s]', 'HCO3-[s]', 'NH3[s]', 'urea[s]', 'Pi[s]', 'sulfate[s]'};
fluxNr = getBounds(model, fluxMets);
n = ceil(sqrt(length(fluxMets)));
for i = 1:length(fluxMets)
    subplot(n,n,i)
    hold all
    flux = allSolutions(fluxNr(i),:);
    if sign(sum(flux)) ==-1
        plot(-flux,'ro-', 'linewidth', 2);
    else
        plot(flux,'bo-', 'linewidth', 2);
    end
    plot([0 nrDists], [0 0], 'k-')
    title(fluxMets{i})
end

%%
dataFiles = {'hepg2-0mm-.txt', 'hepg2-6mm-.txt', 'hepg2-22mm-.txt', 'huh7-0mm-.txt', 'huh7-22mm-.txt'};
dataFolders = {'fluxvalues', 'fluxvalues', 'fluxvalues', 'fluxvalues', 'fluxvalues', 'fluxvalues'};
startPoint = [-1 0 0 -1 0];
%dataFiles = {'hepg2-6mm-.txt', 'hepg2-22mm-.txt'};

results = [];
for i = 1:length(dataFiles)    
    [fluxMets, fluxValues] = loadFluxes(dataFolders{i}, dataFiles{i});
    nrDists = size(fluxValues,2)-1;
    allSolutions = zeros(length(model.rxns), nrDists);
    for j = 1:nrDists
        model = bindFBA(model, fluxMets, fluxValues(:,j)/1000);
        objectiveFunction = 'HumanGrowth';
        model = setParam(model, 'obj', objectiveFunction, 1);
        model = setParam(model, 'ub', objectiveFunction, 0.05);
        solution = solveLinMin(model);
        allSolutions(:,j) = solution.x;    
    end
    results{i} = allSolutions;
end
    
%%

serialized = [];

close all
fluxMets = {'O2[s]', 'NH3[s]', 'biomass[s]'};
fluxNr = getBounds(model, fluxMets);
n = ceil(sqrt(length(fluxMets)));
    for i = 1:length(results)
        allSolutions = results{i};
        for j = 1:length(fluxMets)
            subplot(n,n,j)
            hold all
            flux = allSolutions(fluxNr(j),:);
            x = startPoint(i):(length(flux) -1 + startPoint(i));
            
            if sign(sum(flux)) ==-1
                plot(x, -flux, 'linewidth', 2, 'DisplayName', dataFiles{i});
            else
                plot(x, flux, 'linewidth', 2, 'DisplayName', dataFiles{i});
            end
            %plot([0 nrDists], [0 0], 'k-')
            title(fluxMets{j})
            ylim([0 inf])
        end
    serialized = [serialized;allSolutions(fluxNr,:)'];
        
    end 
legend(gca,'show')
figure
plot(serialized(:,3), -serialized(:,1),'o')

