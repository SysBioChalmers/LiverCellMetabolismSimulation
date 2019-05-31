function plotGenes(A, B, proteinCodingTranscripts, geneConvertionList, geneList)
exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;
     
allGeneNames = {};
maxValue = 0;
for i = 1:length(geneList)
    hold all
    k = 1 + 2*(i-1);
    [data1, data2 ] = getGeneData(A, B, proteinCodingTranscripts, geneList{i});
    allGeneNames{i} = getGeneName(geneConvertionList, geneList{i});
    
    tot1 = sum(data1,1);
    tot2 = sum(data2,1);

    highVal = max([tot1 tot2]);
    
    bar(k, mean(tot1), 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
    bar(k+1, mean(tot2), 'FaceColor', exMap(3,:), 'EdgeColor', 'none')
    errorbar(k, mean(tot1), std(tot1), 'k')
    errorbar(k+1, mean(tot2), std(tot2), 'k')

    [h, p, ci, stats] = ttest2(tot1,tot2);

    text(k, highVal, sprintf('p=%2.2e',p));
    
    maxValue = max(maxValue, highVal);
end
allGeneNames
set(gca,'xtick',0.5+1:2:(2*length(allGeneNames)))
set(gca,'xticklabel',allGeneNames)
%xtickangle(45)
ylabel('TPM')

set(gca,'fontsize',12)
ylim([0 ceil(maxValue)*1.05])
xlim([0.5 k+1.5])
legend('Liver', 'HepG2')
legend('boxoff')
end

