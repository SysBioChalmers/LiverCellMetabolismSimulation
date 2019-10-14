
addpath('src')
load('model/genericHuman2')
celltype = {
    'hepg2' 
    'hepg2'
    'hepg2'
%    'hepg2'    
%     'hepg2'
%     'hepg2'
%     'huh7'
%     'huh7'
    };

conditions = [
    0
    6
    22
%    22
%     6
%     22
%     0
%     22
    ];

conditionName = {
    'hepG2 0mM'
    'hepG2 6mM'
    'hepG2 22mM'
%    'hepG2 22mM (rapid)'   
%     'hepG2 6mM (gln dep)'
%     'hepG2 22mM (gln dep)'
%     'huh7 0mM'
%     'huh7 6 & 22mM'
    };

fluxProfile = [
    2
    2
    2
%    1
%     3
%     3
%     2
%     2
    ];
GAM = 48;
M = 1;
model = setupBiomass(model, GAM, M);

flux = zeros(length(model.rxns),length(conditions));

for i = 1:length(conditions)
    
    if and(strcmp(celltype{i}, 'hepg2'), fluxProfile(i) == 2)
        fileName = ['confidenceIntervalls\output\hepg2-' num2str(conditions(i)) '.tsv'];
        raw = IO(fileName);
        fluxMets = raw(2:end,1);
        fluxIn = cell2nummat(raw(2:end,2))/1000;
    else
        [fluxMets, fluxValues] = loadFluxes('fluxvalues', celltype{i}, conditions(i));
        fluxIn = fluxValues(:,fluxProfile(i))/1000;
    end
        
    model = bindFBA(model, fluxMets, fluxIn);
    solution = solveLinMin(model,1);
    flux(:,i) = solution.x;
end    


%%
close all
clf
%subplot(10,1,1:5);
hold on 
set(gca,'DefaultTextFontSize',14)
set(gca,'DefaultTextFontWeight','bold')

results = zeros(3, length(conditions));
for i = 1:length(conditions)
    [outflux, eqnsOut] = calculateATPValues(model, flux(:,i), 'in');
    results(:,i) = outflux;    
end

plotBarWithPercent(results, eqnsOut, conditionName)    
%xticks([])
%set(findall(gcf,'-property','FontSize'),'FontSize',15)

% pos = get(gca, 'Position');
% pos(1) = pos(1)*2.5;
% pos(3) = pos(3)*0.8;
% set(gca, 'Position', pos)

estimatedMaintainance = M;
plot([0.5 length(conditions)+0.5], estimatedMaintainance * [1 1], 'k-' ,'HandleVisibility','off')

ylabel('ATP synthesis [mmol/gdw/h]')

xlim([0.4 length(conditions)+0.5])

% subplot(10,1,6:10);
figure()
hold on  
set(gca,'DefaultTextFontSize',15)
set(gca,'DefaultTextFontWeight','bold')

results = zeros(5, length(conditions));
for i = 1:length(conditions)
    [outflux, eqnsOut] = calculateATPValues(model, flux(:,i), 'out');
    results(:,i) = outflux;    
end

array2table(100*results(:,:)./sum(results),'RowNames',eqnsOut,'VariableNames',matlab.lang.makeValidName(conditionName))

plotBarWithPercent(results, eqnsOut, conditionName)    
xlim([0.4 length(conditions)+1.5])
% pos = get(gca, 'Position');
% pos(1) = pos(1)*2.5;
% pos(3) = pos(3)*0.8;
% set(gca, 'Position', pos)
% 
estimatedProteinTurnover = 2.57;
plot([0.5 length(conditions)+0.5], estimatedProteinTurnover * [1 1], 'k-','HandleVisibility','off')

ylabel('ATP synthesis [mmol/gdw/h]')
%set(findall(gcf,'-property','FontSize'),'FontSize',15)
xlim([0.4 length(conditions)+0.5])
