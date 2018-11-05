clc
load('model/genericHuman2')
addpath('src')
clf
[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');
fluxData = fluxValues(:,2);
model = setupBiomass(model, 150, 0.5);
model = bindFBA(model, fluxMets, fluxData/1000);
modifiedMetNames = modifyMetNames(model);
solution = solveLinMin(model,1)
mu = -solution.f;

fluxMets = strrep(fluxMets, '[s]', '');

fluxMets(findIndex(fluxMets, 'cystine')) = {'cysteine'};
fluxData(findIndex(fluxMets, 'cysteine')) = 2 * fluxData(findIndex(fluxMets, 'cysteine'));

fluxData(ismember(fluxMets,'tryptophan')) = [];
fluxMets(ismember(fluxMets,'tryptophan')) = [];



essential = {'histidine'
            'isoleucine'
            'leucine'
            'lysine'
            'methionine'
            'phenylalanine'
            'threonine'
            'valine'};
        
        

rxnIndx = findIndex(model.rxns, 'human_proteinPool');
metIndex = model.S(:,rxnIndx);
aminoAcids = model.metNames(metIndex<0);
aminoAcidStochiometry = model.S(findIndex(modifiedMetNames,'human_protein_pool[c]'),findIndex(model.rxns,'HumanGrowth'));
metProtein = aminoAcidStochiometry * metIndex(metIndex<0);

%Metabolite pool:
rxnIndx = findIndex(model.rxns, 'metabolitePool');
aminoPoolStochiometry = model.S(findIndex(modifiedMetNames,'metabolitePool[c]'),findIndex(model.rxns,'HumanGrowth'));
freeMetAmount = model.S(metIndex<0,rxnIndx) * aminoPoolStochiometry;

metAmount = freeMetAmount + metProtein;
metAmount = 1000* metAmount;

mapedFluxes = zeros(length(metAmount),1);

for i = 1:length(fluxMets)
    curIndx = findIndex(aminoAcids,fluxMets{i});
    if isempty(curIndx)
       disp(fluxMets{i});
    else
        mapedFluxes(curIndx) = fluxData(i);
    end
end

%mu = min(-mapedFluxes(ismember(aminoAcids, essential))./metAmount(ismember(aminoAcids, essential)));
%mu = 0.033;

ub = round(max(metAmount)*1.02,1);
hold all
plot([0 ub], [0 0], 'k-');
plot([0 ub], [0 -mu*ub], 'color', [0.8 0.8 0.8], 'linewidth', 2);

text(ub*0.8, -mu*ub, sprintf('µ = %2.3f',mu))


lineMets = {'glutamine', 'alanine', 'glutamate'};

for i = 1:length(lineMets)
    curMet = findIndex(aminoAcids, lineMets{i});
    x = metAmount(curMet);
    y1 = -mu*x;
    y2 = mapedFluxes(curMet);
    plot(x * [1 1], [y1 y2], 'color', [0.8 0.8 0.8]);
    text(x, y1 + 0.5 * (y2-y1), sprintf('%2.2f',(y2-y1)));
end


flippedLocation = {'cysteine', 'methionine', 'arginine', 'isoleucine', 'phenylalanine', 'valine', 'leucine'};

fprintf('\n')
for i = 1:length(metAmount)
    if ismember(aminoAcids{i}, essential) 
        scatter(metAmount(i), mapedFluxes(i), 50, 'filled', 'o', 'MarkerFaceColor', [216 85 39]/256, 'LineWidth', 0.001);
    else
        scatter(metAmount(i), mapedFluxes(i), 50, 'filled', 'o', 'MarkerFaceColor', [17 115 187]/256, 'LineWidth', 0.001);
    end
    
    if ismember(aminoAcids{i}, flippedLocation) 
        text(metAmount(i)-0.004, mapedFluxes(i)-0.002, aminoAcids{i}, 'fontsize', 15, 'Rotation',45, 'HorizontalAlignment', 'right', 'color', [0.5 0.5 0.5]);        
    else
        text(metAmount(i)+0.002, mapedFluxes(i)+0.002, aminoAcids{i}, 'fontsize', 15, 'Rotation',45, 'color', [0.5 0.5 0.5]);        
    end
    fprintf('%s\t%2.2f\n', aminoAcids{i}, mapedFluxes(i)+mu*metAmount(i))
end

%plot cystine
cysteineMet = ismember(aminoAcids, 'cysteine');
scatter(metAmount(cysteineMet), mapedFluxes(cysteineMet)/2, 50, 'filled', 'o', 'MarkerFaceColor', [115 115 115]/256, 'LineWidth', 0.001);
plot(metAmount(cysteineMet) * [1 1], mapedFluxes(cysteineMet) * [0.5 1], 'color', [115 115 115]/256)
% fluxMets(findIndex(fluxMets, 'cystine')) = {'cysteine'};
% fluxData(findIndex(fluxMets, 'cysteine')) = 2 * fluxData(findIndex(fluxMets, 'cysteine'));


xlabel('Biomass')
ylabel('Flux')
%ylim([-0.11 0.06])
xlim([0 ub])
ylim(1000*[-0.15 0.0505])
set(gca,'FontSize', 15, 'FontWeight', 'bold');

