addpath('src')
load('model/genericHuman2')
conditions = {
    'hepg2-0mm-.txt' 
    'hepg2-6mm-.txt'
    'hepg2-22mm-.txt'     
    'hepg2-6mm-.txt'
    'hepg2-22mm-.txt'
    'huh7-0mm-.txt'
    'huh7-22mm-.txt'
    };

conditionName = {
    'hepG2 0mM'
    'hepG2 6mM'
    'hepG2 22mM'    
    'hepG2 6mM (gln dep)'
    'hepG2 22mM (gln dep)'
    'huh7 0mM'
    'huh7 6 & 22mM'
    };

fluxProfile = [
    2
    2
    2
    3
    3
    2
    2
    2
    ];

model = setupBiomass(model, 150, 0.5);
metNames = modifyMetNames(model);
glutamateMet = findIndex(metNames, 'glutamate[c]');
glutamateRxns = find(model.S(glutamateMet,:));
glutamateStochiometry = full(model.S(glutamateMet,glutamateRxns));


groupNames{1} = 'BCAA';
reactionGroups{1} = {'HMR_3747', 'HMR_3777', 'HMR_6923'};

groupNames{2} = 'Biosynthesis';
reactionGroups{2} = {'HMR_4300', 'HMR_3903', 'HMR_4194', 'HMR_4034', 'HMR_4808', 'HMR_4406', 'HMR_4046', 'HMR_4658'};

groupNames{3} = 'Deamination of other amino acids';
reactionGroups{3} = {'HMR_6725',  'HMR_6768', 'HMR_4596'};

groupNames{4} = 'Biomass';
reactionGroups{4} = {'human_proteinPool', 'metabolitePool', 'HMR_4324'};

groupNames{5} = 'Serine synthesis';
reactionGroups{5} = {'HMR_3841'};

groupNames{6} = 'Glutamine synthesis';
reactionGroups{6} = {'HMR_3890'};

groupNames{7} = 'Excretion';
reactionGroups{7} = {'FreeGlutamateTransport', 'HMR_4931'};

groupNames{8} = 'Exchange Mitochondria';
reactionGroups{8} = {'HMR_3829', 'HMR_5122',  'HMR_3825'};

flux = zeros(length(model.rxns),length(conditions));

for i = 1:length(conditions)
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', conditions{i});
    model = bindFBA(model, fluxMets, fluxValues(:,fluxProfile(i))/1000);
    solution = solveLinMin(model,1);
    flux(:,i) = solution.x;
end    

growth = flux(find(model.c),:);

flux = flux * 1000;

%%
reactionsWithFlux = glutamateRxns(sum(abs(flux(glutamateRxns,:)),2)>10^-6);
results = zeros(length(reactionGroups), length(conditionName));

for i = 1:length(conditionName)
    fprintf('\t%s', conditionName{i});
end

fprintf('\n');



for i = 1:length(reactionGroups)   
   sumOfFlux = zeros(1, length(conditions));
   curGroup = reactionGroups{i};
   for j = 1:length(curGroup)
      curRxn = findIndex(model.rxns, curGroup{j});
      stochiometry = full(model.S(glutamateMet,curRxn)); 
      sumOfFlux = sumOfFlux + stochiometry * flux(curRxn,:);
      
      %double check that all reactions are counted
      reactionsWithFlux(ismember(reactionsWithFlux, curRxn)) = [];
   end
   
   results(i,:) = sumOfFlux;
   
   fprintf('%s', groupNames{i});
   for j = 1:length(sumOfFlux)
       fprintf('\t%2.2f', sumOfFlux(j));
   end
   fprintf('\n');
end


sumOfFlux = zeros(1, length(conditions));
for i = 1:length(reactionsWithFlux)
  curRxn = reactionsWithFlux(i);
  stochiometry = full(model.S(glutamateMet,curRxn)); 
  sumOfFlux = sumOfFlux + stochiometry * flux(curRxn,:);
end


fprintf('other');
for j = 1:length(sumOfFlux)
   fprintf('\t%2.2f', sumOfFlux(j));
end
fprintf('\n');

subplot(10,1,1:5);
positiveResults = results;
positiveResults(positiveResults<0) = 0;
plotBarWithPercent(positiveResults([1 2 3 8],:), groupNames([1 2 3 8]), conditionName) 
xlim([0 80])
xticks([])
%set(gca,'TickLength',[0 0])
set(gca,'box','off'); %remove the box

pos = get(gca, 'Position');
pos(1) = pos(1)*3;
pos(3) = pos(3)*0.7;
set(gca, 'Position', pos)

subplot(10,1,6:10);
negativeResults = -results;
negativeResults(negativeResults<0) = 0;
plotBarWithPercent(negativeResults(4:8,:), groupNames(4:8), conditionName) 
xlim([0 80])
%set(gca,'TickLength',[0 0])
set(gca,'box','off'); %remove the box

pos = get(gca, 'Position');
pos(1) = pos(1)*3;
pos(3) = pos(3)*0.7;
set(gca, 'Position', pos)
%%
figure()
color2 = [215 86 40]/256;
color1 = [93 155 211]/256;

xvals = linspace(0, 0.04);

subplot(10,1,1:5);
hold all

scatter(growth, results(2,:), 30, 'MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
scatter(growth, -results(4,:), 30, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)

dlm = fitlm(growth,results(2,:),'Intercept',false);
km = dlm.Coefficients.Estimate;
plot(xvals, xvals*km, 'color', color1)
text(0.03, 0.03*km, sprintf('r2=%2.2f\np=%2.2e', dlm.Rsquared.Ordinary, coefTest(dlm))) 

dlm = fitlm(growth, -results(4,:),'Intercept',false);
km = dlm.Coefficients.Estimate;
plot(xvals, xvals*km, 'color', color2)
text(0.03, 0.03*km, sprintf('r2=%2.2f\np=%2.2e', dlm.Rsquared.Ordinary, coefTest(dlm))) 

xlim([0 0.04])
ylim([0 30])
xticks([])
legend(groupNames{2}, groupNames{4})
legend boxoff
subplot(10,1,6:10);
hold all

scatter(growth, -results(7,:), 30, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)

dlm = fitlm(growth, -results(7,:));
km = dlm.Coefficients.Estimate;
plot(xvals, km(1) + xvals*km(2), 'color', color2)
text(0.03, km(1) + 0.03*km(2), sprintf('r2=%2.2f\np=%2.2e', dlm.Rsquared.Ordinary, coefTest(dlm))) 

legend(groupNames{7})
legend boxoff

 xlim([0 0.04])
 xlabel('growth rate')
% ylim([0 30])

%%
figure()
hold all
aspMet = findIndex(metNames, 'aspartate[c]');
%aspRxns = find(model.S(aspMet,:));
%rxnsWithFlux = sum(abs(flux(aspRxns,:)),2)>10^-6;
%aspRxns = aspRxns(rxnsWithFlux);
aspRxns = [229	294	336	530	8129 8145 5445 6705];

aspStochiometry = full(model.S(aspMet,aspRxns));
aspSink = flux(aspRxns,:);
for i = 1:size(aspSink,1)
    aspSink(i,:) = aspSink(i,:) .* aspStochiometry(i);
end

aspSink = -sum(aspSink);

scatter(growth, -results(8,:), 30, 'MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
scatter(growth, aspSink, 30, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.6)
xlim([0 0.04])
ylim([-10 35])

dlm = fitlm(growth, aspSink, 'Intercept',false);
km = dlm.Coefficients.Estimate;
plot(xvals, xvals*km, 'color', color2)

