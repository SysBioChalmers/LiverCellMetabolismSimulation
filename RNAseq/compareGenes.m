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
hold all
[data1, data2 ] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000060982');

tot1 = sum(data1);
tot2 = sum(data2);

bar(1, mean(tot1), 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
bar(2, mean(tot2), 'FaceColor', exMap(3,:), 'EdgeColor', 'none')
errorbar(1, mean(tot1), std(tot1), 'k')
errorbar(2, mean(tot2), std(tot2), 'k')

[h p] = ttest2(tot1,tot2)
text(1.5, 25, num2str(p))

[data1, data2 ] = getGeneData(A, B, proteinCodingTranscripts, 'ENSG00000105552');
tot1 = sum(data1);
tot2 = sum(data2);

bar(3, mean(tot1), 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
bar(4, mean(tot2), 'FaceColor', exMap(3,:), 'EdgeColor', 'none')

errorbar(3, mean(tot1), std(tot1), 'k')
errorbar(4, mean(tot2), std(tot2), 'k')

[h p] = ttest2(tot1,tot2)
text(3.5, 25, num2str(p))


ylim([0 30])
xlim([0.5 4.5])
legend('Liver', 'HepG2')

