cellLine = {
    'hepg2'
    'hepg2'
    'hepg2' 
    'huh7'
    'huh7'    
    'huh7'
    };

condition = {
    '0mm'
    '6mm'
    '22mm' 
    '0mm'
    '6mm'
    '22mm' 
    };

metabolites = {'alanine[s]'
    'arginine[s]'
    'asparagine[s]'
    'aspartate[s]'
    'cystine[s]'
    'glutamate[s]'
    'glutamine[s]'
    'glycine[s]'
    'histidine[s]'
    'isoleucine[s]'
    'leucine[s]'
    'lysine[s]'
    'methionine[s]'
    'ornithine[s]'
    'phenylalanine[s]'
    'serine[s]'
    'taurine[s]'
    'threonine[s]'
    'tyrosine[s]'
    'valine[s]'};

flux = zeros(length(metabolites),length(condition));


for i = 1:length(condition)
    expData = loadExpdata('../data', cellLine{i}, condition{i});
    
    
    for j = 1:length(metabolites)
        if isKey(expData, metabolites{j})
            curData = expData(metabolites{j});
            expDataFirst = curData(2,1);
            flux(j,i) = expDataFirst;
        end
    end
end    

hold all
xvalues = 1:length(metabolites);
average = mean(flux');
standardDev = std(flux');

plot(flux, 'ko')
errorbar(xvalues, average, standardDev, 'ro');

set(gca,'xtick', xvalues)
set(gca,'xticklabel', metabolites)
set(gca, 'XTickLabelRotation', 45);
figure
plot(standardDev./average, 'ko')

set(gca,'xtick', xvalues)
set(gca,'xticklabel', metabolites)
set(gca, 'XTickLabelRotation', 45);

