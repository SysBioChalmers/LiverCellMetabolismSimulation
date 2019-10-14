clc
load('model/genericHuman2')
addpath('src')

NGAM = 1;
GAM = 48;

close all
celltype = {
    'hepG2'
    'hepG2' 
    'hepG2'    
    };

conditions = [
    0
    6
    22
    ];

fluxProfile = [
    2
    2
    2
    ];

result = containers.Map();
growthRate = zeros(length(conditions),1);
ATPrate = zeros(length(conditions),5);

objectiveFunction = 'HumanGrowth';
ATPsynthesis = 'human_ATPMaintainance';
ATPsynthRxn = findIndex(model.rxns, 'human_ATPMaintainance');

model = setupBiomass(model, 0, 0);
epsilon = 10^-6;
for i = 1:length(conditions)
    if strcmp(celltype{i}, 'hepG2')
        fileName = ['confidenceIntervalls\output\hepg2-' num2str(conditions(i)) '.tsv'];
        raw = IO(fileName);
        fluxMets = raw(2:end,1);
        fluxValues = cell2nummat(raw(2:end,2))/1000;
        fluxError = cell2nummat(raw(2:end,3:4))/1000;
    end
    
    %Calculate growth rate
    model = bindFBA(model, fluxMets, fluxValues);
    solution = solveLin(model,1);
    growthRate(i) = -solution.f;
    
    %maximize ATP at this growth rate
    lactateRxn = findIndex(fluxMets, 'L-lactate[s]');
    glucoseRxn = findIndex(fluxMets, 'glucose[s]');
    for j = 1:5
        tmp = fluxValues;
        
        switch j
            case 1
                tmp = fluxValues;
            case 2
                tmp(lactateRxn) = fluxError(lactateRxn,1);
            case 3
                tmp(lactateRxn) = fluxError(lactateRxn,2);
            case 4
                tmp(glucoseRxn) = fluxError(glucoseRxn,1);
            case 5
                tmp(glucoseRxn) = fluxError(glucoseRxn,2);      
        end
                
        model = bindFBA(model, fluxMets, tmp);
        
        model = setParam(model, 'lb', objectiveFunction, growthRate(i) - epsilon);
        model = setParam(model, 'ub', objectiveFunction, growthRate(i) + epsilon);   

        model.ub(ATPsynthRxn) = 1000;
        model.c = zeros(length(model.rxns),1);
        model.c(ATPsynthRxn) = 1;
        solution = solveLin(model,1);
        ATPrate(i,j) = -solution.f;
    end
    
end    


%%
hold all
area(xvals, GAM*xvals + NGAM, 'HandleVisibility', 'off', 'LineStyle', 'none', 'FaceColor', [0.8 0.8 0.8])
plot(growthRate, ATPrate(:,1), 'ko-', 'MarkerFaceColor', 'black')

for i = 2:size(ATPrate,2)
    scatter(growthRate, ATPrate(:,i), 'fill')
end
xlim([0 0.04])
ylim([0 8])

xvals = [0 0.04];
xlabel('Specific growth rate [1/h]')
ylabel('ATP synthesis capacity [mmol ATP/gdw/h]')

legend({'MLE', 'lb lactate', 'ub lactate', 'lb glucose', 'ub glucose'}, 'location', 'nw')
legend boxoff    

text(growthRate, ATPrate(:,1), {'0 mM', '6 mM', '22 mM'})

