addpath('RNAseq')
exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;

A = importdata('RNAseq/tisue.txt');
B = importdata('RNAseq/seq.txt');

C = IO('RNAseq/ID_conversion_key.txt');
proteinCodingTranscripts = C(:,2);
geneConvertionList = C(:,[1 5]);

%%
clf
geneList{1} = {'ENSG00000129596'}; %CDO1
geneList{2} = {'ENSG00000116761'}; %CTH
geneList{3} = {'ENSG00000244005'}; %NFS1
geneList{4} = {'ENSG00000105755'}; %ETHE1
geneList{5} = {'ENSG00000060982', 'ENSG00000105552'}; %BCAT 1 and 2
geneList{6} = {'ENSG00000167701', 'ENSG00000166123'}; %ALT 1 and 2
geneList{7} = {'ENSG00000118520', 'ENSG00000081181'}; %ARG 1 and 2



nRows = ceil(length(geneList)/2);

for i = 1:length(geneList)
    subplot(nRows,2,i)
    plotGenes(A, B, proteinCodingTranscripts, geneConvertionList, geneList{i})
end


%%
figure()     
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

    [h p] = ttest2(tot1,tot2);
    [p,h] = ranksum(tot1,tot2)
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



