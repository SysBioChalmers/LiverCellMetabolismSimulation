A = importdata('tisue.txt');
B = importdata('seq.txt');
exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;
C = importdata('ID_conversion_key.txt');
proteinCodingTranscripts = C.textdata(:,2);

%%
clf

hold all
%[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000183044');

[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000115419'); %GLS1

[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000117054'); %propionyl-CoA carboxylase

[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000244005'); %Cysteine desulfurase, sugests cytoplasmic isoform 

[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000105755'); %Cysteine desulfurase, sugests cytoplasmic isoform 

[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000021826'); %Carbamoyl-phosphate synthase 1

[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000130707'); %Argininosuccinate synthase
[data1, data2, transc1, transc2] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000126522'); %Argininosuccinate lyase 



	

transc1
maxValue = 0;
for i = 1:size(data1,1)
    k = 1 + 2*(i-1);

    tot1 = data1(i,:);
    tot2 = data2(i,:);

    highVal = max([tot1 tot2]);
    
    bar(k, mean(tot1), 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
    bar(k+1, mean(tot2), 'FaceColor', exMap(3,:), 'EdgeColor', 'none')
    errorbar(k, mean(tot1), std(tot1), 'k')
    errorbar(k+1, mean(tot2), std(tot2), 'k')

    [h p] = ttest2(tot1,tot2);
    text(k, highVal, sprintf('p=%2.2e',p));
    
    maxValue = max(maxValue, highVal);
end
set(gca,'xtick',1+1:2:(2*size(data1,1)))
set(gca,'xticklabel',transc1)
xtickangle(45)
ylabel('TPM')

set(gca,'fontsize',12)
ylim([0 ceil(maxValue)*1.05])
xlim([0.5 k+1.5])
legend('Liver', 'HepG2')
legend('boxoff')
