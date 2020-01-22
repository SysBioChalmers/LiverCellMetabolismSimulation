clc
load('model/genericHuman2')
addpath('src')

color = get(gca,'colororder');

tresh = 10^-6;

cellType = 'hepg2';
condition = '22';
setErrorBounds = false;

if strcmp(cellType, 'hepg2')
    fileName = ['confidenceIntervalls\output\hepg2-' num2str(condition) '.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    fluxValues = cell2nummat(raw(2:end,2));
    fluxData = fluxValues/1000;
    fluxError = cell2nummat(raw(2:end,3:4));
else
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', cellType, condition);
    fluxData = fluxValues(:,2)/1000;
end

model = setupBiomass(model, 48, 1);
model = bindFBA(model, fluxMets, fluxData);

%Prevent cycle around proline
model.ub(findIndex(model.rxns, 'HMR_4963')) = 0;

solution = solveLinMin(model,1);
flux = solution.x;

constrained = abs(flux)<tresh;
model.ub(constrained) = 0;
model.lb(constrained) = 0;

constrained = flux>tresh;
model.ub(constrained) = flux(constrained)+tresh;
model.lb(constrained) = max([-flux(constrained)-tresh model.lb(constrained)]')';

constrained = -flux>tresh;
model.lb(constrained) = flux(constrained)-tresh;
model.ub(constrained) = min([-flux(constrained)-tresh model.ub(constrained)]')';


exMets = getBounds(model, fluxMets);
model.ub(exMets) = 1000;
model.lb(exMets) = -1000;


% model.lb(ismember(model.rxns, {'HMR_4931', 'FreeGlutamateTransport'})) = -1000;
% model.ub(ismember(model.rxns, {'HMR_4931', 'FreeGlutamateTransport'})) = 1000;

model.ub(findIndex(model.rxns, 'human_ATPMaintainance')) = 1000;

%%
%Glutamate model

gluModel = model;
gluRxn = getBounds(model, {'glutamate[s]'});
gluFlux = fluxData(ismember(fluxMets, 'glutamate[s]'));

concentrations = linspace(0.5, 1, 50)';
results = zeros(length(concentrations), 2);

%Without BCAA constraint
for i = 1:length(concentrations)
    gluModel.ub(gluRxn) = gluFlux*concentrations(i)+tresh;
    gluModel.lb(gluRxn) = gluFlux*concentrations(i)-tresh;
    
    solution = solveLin(gluModel,1);
    results(i, 1) = -solution.f;
end

%With BCAA constraint
scale = 0.999;
bcatModel = model;
bcatModel.lb(findIndex(model.rxns, 'HMR_3777')) = flux(findIndex(model.rxns, 'HMR_3777'))*scale;
bcatModel.lb(findIndex(model.rxns, 'HMR_6923')) = flux(findIndex(model.rxns, 'HMR_6923'))*scale;
bcatModel.lb(findIndex(model.rxns, 'HMR_3747')) = flux(findIndex(model.rxns, 'HMR_3747'))*scale;    
    
for i = 1:length(concentrations)
    bcatModel.ub(gluRxn) = gluFlux*concentrations(i)+tresh;
    bcatModel.lb(gluRxn) = gluFlux*concentrations(i)-tresh;    

    solution = solveLin(bcatModel,1);
    results(i, 2) = -solution.f;
end

hold all
plot(concentrations, results(:,1)/max(results(:)), '-', 'linewidth', 2)

included = results(:,2)>tresh;
plot(concentrations(included), results(included,2)/max(results(:)), '-', 'linewidth', 2)
plot([0 1], [0 1], 'k-', 'HandleVisibility', 'off')
text(0.92, 0.9, 'X=Y')

ylim([0 inf])
growthConcentrations = results>tresh;
%trend = polyfit(concentrations(growthConcentrations), results(growthConcentrations), 1);
%plot(concentrations, polyval(trend, concentrations), 'r--');
xlabel('Relative glutamate excretion')
ylabel('Relative growth rate')
%ylim([0 0.04])

xlim([0.5 1])
legend({'unconstrained BCAT', 'constrained BCAT'}, 'location', 'se')
legend boxoff

%%
figure()
%Nucleotide model
nucleotideModel = model;

nucleotideNames = {'human_DNAPool', 'human_RNAPool'};
nucleotideRxns = ismember(model.rxns, nucleotideNames);
nucleotideFlux = flux(ismember(model.rxns, nucleotideNames));
gluRxn = getBounds(model, {'glutamate[s]'});
gluFlux = flux(gluRxn);

objectiveFunction = 'HumanGrowth';
concentrations = linspace(0.5, 1, 50)';
results = zeros(length(concentrations), 6);

bcatModel = model;
bcatModel.lb(findIndex(model.rxns, 'HMR_3777')) = flux(findIndex(model.rxns, 'HMR_3777'))*scale;
bcatModel.lb(findIndex(model.rxns, 'HMR_6923')) = flux(findIndex(model.rxns, 'HMR_6923'))*scale;
bcatModel.lb(findIndex(model.rxns, 'HMR_3747')) = flux(findIndex(model.rxns, 'HMR_3747'))*scale;    

mitoRxn = findIndex(model.rxns, 'HMR_5067');
mitoflux = flux(mitoRxn);

for i = 1:length(concentrations)
    nucleotideModel.ub(nucleotideRxns) = nucleotideFlux*concentrations(i)+tresh;
    nucleotideModel.lb(nucleotideRxns) = nucleotideFlux*concentrations(i)-tresh;
    
    nucleotideModel = setParam(nucleotideModel, 'obj', objectiveFunction, 1);    
    solution = solveLin(nucleotideModel, 1);
    results(i, 1) = -solution.f;
    
    nucleotideModel = setParam(nucleotideModel, 'obj', gluRxn, 1);
    solution = solveLin(nucleotideModel, 1);   
    results(i, 2) = solution.x(gluRxn); 
    
    nucleotideModel = setParam(nucleotideModel, 'obj', gluRxn, -1);
    solution = solveLin(nucleotideModel, 1);
    results(i, 3) = solution.x(gluRxn);
    
    bcatModel.ub(nucleotideRxns) = nucleotideFlux*concentrations(i)+tresh;
    bcatModel.lb(nucleotideRxns) = nucleotideFlux*concentrations(i)-tresh;    
    bcatModel = setParam(bcatModel, 'obj', gluRxn, 1);
    solution = solveLin(bcatModel, 1);   
    
    if length(solution.x)>1
        results(i, 4) = solution.x(gluRxn); 

        bcatModel = setParam(bcatModel, 'obj', gluRxn, -1);
        solution = solveLin(bcatModel, 1);
        results(i, 5) = solution.x(gluRxn);    
        solution.x(mitoRxn)
%         tmpModel = bcatModel;
%         tmpModel.lb(mitoRxn) = mitoflux-tresh;
%         tmpModel.ub(mitoRxn) = mitoflux+tresh;
%         
%         tmpModel = setParam(tmpModel, 'obj', gluRxn, 1);
%         solution = solveLin(tmpModel, 1);
%         results(i, 6) = solution.x(gluRxn);           
    end

end

hold all
%plot(concentrations, results(:,1)/results(end,1), 'k-', 'linewidth', 2)
X = [concentrations; flip(concentrations)];
Y = [results(:,2); flip(results(:,3))]/gluFlux;
fill(X, Y, color(1,:), 'linestyle', 'none', 'facealpha', 0.4)


Y = [results(:,4); flip(results(:,5))]/gluFlux;
fill(X, Y, color(2,:), 'linestyle', 'none', 'facealpha', 0.4)

Y = [results(:,6); flip(results(:,5))]/gluFlux;
fill(X, Y, color(2,:), 'linestyle', 'none', 'facealpha', 0.4)


plot([0 1], [0 1], 'k-', 'HandleVisibility', 'off')
text(0.92, 0.9, 'X=Y')


xlabel('Relative nucleotide synthesis rate')
ylabel('Relative glutamate excretion rate')
%ylim([0 0.04])

xlim([0.5 1])
ylim([-0.5 1.5])
legend({'unconstrained BCAT', 'constrained BCAT'}, 'location', 'NE')
legend boxoff

%%
figure()
glnModel = model;
glnModel.lb(findIndex(model.rxns, 'HMR_3777')) = flux(findIndex(model.rxns, 'HMR_3777'))*scale;
glnModel.lb(findIndex(model.rxns, 'HMR_6923')) = flux(findIndex(model.rxns, 'HMR_6923'))*scale;
glnModel.lb(findIndex(model.rxns, 'HMR_3747')) = flux(findIndex(model.rxns, 'HMR_3747'))*scale;    

concentrations = linspace(0.5, 1, 50)';
glnRxn = getBounds(model, {'glutamine[s]'});
glnFlux = fluxData(ismember(fluxMets, 'glutamine[s]'));
results = zeros(length(concentrations), 3);

level = 0.8;
gluRxn = getBounds(model, {'glutamate[s]'});
gluFlux = flux(gluRxn);

for i = 1:length(concentrations)
    glnModel.ub(glnRxn) = glnFlux*concentrations(i)+tresh;
    glnModel.lb(glnRxn) = glnFlux*concentrations(i)-tresh;

    glnModel1 = glnModel;
    solution = solveLin(glnModel1,1);
    results(i, 1) = -solution.f;
    
    glnModel2 = glnModel;
    glnModel2.lb(gluRxn) = level*gluFlux - tresh;
    glnModel2.ub(gluRxn) = level*gluFlux + tresh;
    
    solution = solveLin(glnModel2,1);
    results(i, 2) = -solution.f;
    
    glnModel3 = glnModel;    
    glnModel3.lb(gluRxn) = gluFlux - tresh;
    glnModel3.ub(gluRxn) = gluFlux + tresh;
        
    solution = solveLin(glnModel3,1);
    results(i, 3) = -solution.f;    
end

% %With BCAA constraint
% scale = 0.999;
% GluGLNmodel = model;
% GluGLNmodel.ub(gluRxn) = 0;
%     
% for i = 1:length(concentrations)
%     GluGLNmodel.ub(glnRxn) = glnFlux*concentrations(i)+tresh;
%     GluGLNmodel.lb(glnRxn) = glnFlux*concentrations(i)-tresh; 
% 
%     solution = solveLin(GluGLNmodel,1);
%     results(i, 2) = -solution.f;
% end

figure()
hold all
plot(concentrations, results(:,1)/max(results(:)), '-', 'linewidth', 2)
plot(concentrations, results(:,2)/max(results(:)), '-', 'linewidth', 2)
plot(concentrations, results(:,3)/max(results(:)), '-', 'linewidth', 2)

plot([0 1], [0 1], 'k-', 'HandleVisibility', 'off')
text(0.92, 0.9, 'X=Y')

xlim([0.5 1])
ylim([0 inf])
growthConcentrations = results>tresh;
%trend = polyfit(concentrations(growthConcentrations), results(growthConcentrations), 1);
%plot(concentrations, polyval(trend, concentrations), 'r--');
xlabel('Relative glutamine uptake')
ylabel('Relative growth rate')
%ylim([0 0.04])
legend({'Glu free', 'Glu low', 'Glu high'}, 'location', 'SW')
legend box off
