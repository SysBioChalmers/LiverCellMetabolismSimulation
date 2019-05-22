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
geneList = { 'ENSG00000167701', 'ENSG00000166123'}; %cytosolic vs mitochondrial GPT1 vs GPT2

geneList = {'ENSG00000120053','ENSG00000169154', 'ENSG00000125166'}; %Cytosolic, Cytosolic, Mitochondrial Aspartate
%geneList = {'ENSG00000162174','ENSG00000166183', 'ENSG00000070669'}; %Asparagine synthesis

%geneList = {'ENSG00000151012'}; %Glut-cys2 antiporter
geneList = {'ENSG00000001084','ENSG00000023909','ENSG00000100983'}; %GSH forming

%
geneList = {'ENSG00000137944', 'ENSG00000171097'}; %Kynurenine pathway

geneList = {'ENSG00000132330'}; %Selenocysteine lyase



geneList =  {'ENSG00000006756', 'ENSG00000157399', 'ENSG00000205667', 'ENSG00000062096', 'ENSG00000135677', 'ENSG00000101846', 'ENSG00000180801', 'ENSG00000100299', 'ENSG00000196562', 'ENSG00000164291', 'ENSG00000010404', 'ENSG00000137573', 'ENSG00000141012', 'ENSG00000113273', 'ENSG00000141337'}; %sulfatases
geneList =  {'ENSG00000196562', 'ENSG00000164291', 'ENSG00000141012', 'ENSG00000113273', 'ENSG00000141337'}; %sulfatases

geneList = {'ENSG00000115665', 'ENSG00000070214', 'ENSG00000129353'}; %High af choline transporter

geneList = {'ENSG00000167720', 'ENSG00000135094'}; %Pyr -> Serine

geneList = {'ENSG00000119640', 'ENSG00000170634'}; %Acetyl-P to Acetate predicts death in liver cancer

geneList = {'ENSG00000100116', 'ENSG00000131480', 'ENSG00000131471', 'ENSG00000149476'}; %pyruvate->glycine pathway

geneList = {'ENSG00000179761', 'ENSG00000123453', 'ENSG00000124713'}; %Sarcosine


geneList = {'ENSG00000136881'}; %Serine -> glycine THF

geneList = {'ENSG00000198380', 'ENSG00000113552', 'ENSG00000163281'}; %hexosamine pathway

geneList = {'ENSG00000139344'};

geneList = {'ENSG00000172890'}; %Leukotriene metabolism

geneList = {'ENSG00000198380', 'ENSG00000113552', 'ENSG00000163281', 'ENSG00000100522', 'ENSG00000013375', 'ENSG00000197355', 'ENSG00000095380'}; %hexosamine pathway

geneList = {'ENSG00000244005'}; %Cysteine desulfurase, mitochondrial

geneList = {'ENSG00000171314', 'ENSG00000092621','ENSG00000135069', 'ENSG00000146733'}; %Serine synthesis

geneList = {'ENSG00000198380', 'ENSG00000113552', 'ENSG00000163281', 'ENSG00000100522', 'ENSG00000013375', 'ENSG00000197355', 'ENSG00000095380'}; %hexosamine pathway

geneList = {'ENSG00000120053','ENSG00000116761', 'ENSG00000125166', 'ENSG00000160200', 'ENSG00000274276', 'ENSG00000115919'}; %Beatine



geneList = {'ENSG00000102996', 'ENSG00000151694', 'ENSG00000137845', 'ENSG00000143537', 'ENSG00000168615', 'ENSG00000107959', 'ENSG00000084073'}; %Matrix metallopeptidase 1

geneList = {'ENSG00000138078', 'ENSG00000109861',  'ENSG00000254986', 'ENSG00000085377', 'ENSG00000176393'}; %Exopeptidase

geneList = {'ENSG00000124767', 'ENSG00000063854'}; %Lactoylglutathione lyase

geneList = {'ENSG00000002587', 'ENSG00000125430',  'ENSG00000153976', 'ENSG00000162040', 'ENSG00000182601', 'ENSG00000249853'}; %HMR_7222, heparan sulfate proteoglycan

geneList = {'ENSG00000130348', 'ENSG00000257218',  'ENSG00000059691'}; 

geneList = {'ENSG00000103356', 'ENSG00000136628'}; 

%transaminases
geneList = {'ENSG00000092295', 'ENSG00000198959', 'ENSG00000124491', 'ENSG00000125780', 'ENSG00000163810', 'ENSG00000281886', 'ENSG00000104055', 'ENSG00000166948', 'ENSG00000159495'};

geneList = {'ENSG00000139631', 'ENSG00000144644'}; %Cystine oxidation CSAD GADL1
	

geneList = {'ENSG00000110887', 'ENSG00000172482'}; %Cysteine desulfurase, mitochondrial

geneList = {'ENSG00000198431', 'ENSG00000184470', 'ENSG00000197763'}; %HMR_3995 cystine metabolism



geneList = {'ENSG00000116761'}; %Cystathionine gamma-lyase 

geneList = {'ENSG00000120053', 'ENSG00000125166', 'ENSG00000100344', 'ENSG00000183044', 'ENSG00000109576'}; %GOT1 GOT2

geneList = {'ENSG00000120053', 'ENSG00000125166'}; %GOT1 GOT2



geneList = {'ENSG00000176974', 'ENSG00000182199'}; %Serine -> glycine THF



geneList = {'ENSG00000128311', 'ENSG00000128309'}; %MPST TST

geneList = {'ENSG00000081181', 'ENSG00000118520'}; %HMR_3816 Arginine to Ornithine

geneList = {'ENSG00000083123', 'ENSG00000248098'}; %BCKDHB for BCAA

geneList = {'ENSG00000111058', 'ENSG00000131069','ENSG00000154930'}; %BCKDHB for BCAA

geneList = {'ENSG00000114054', 'ENSG00000175198'}; %propanoyl-CoA[m] -> methylmalonyl-CoA[m]

geneList = {'ENSG00000119711', 'ENSG00000113492', 'ENSG00000183044'}; %Side reactions for valine metabolism

geneList = {'ENSG00000176715', 'ENSG00000146085'}; %Main reactions valine metabolism


geneList = {'ENSG00000172340', 'ENSG00000163541', 'ENSG00000136143'}; %Suc CoA

%geneList = {'ENSG00000124370'}; %Main reactions isoleucine metabolism


geneList = {'ENSG00000146085'}; %Main reactions valine metabolism

geneList = {'ENSG00000175198'}; %Main reactions valine metabolism

%Whole valine degradation pathway
geneList = {'ENSG00000060982', 'ENSG00000248098', 'ENSG00000137992', 'ENSG00000117054', 'ENSG00000122971', 'ENSG00000196177', 'ENSG00000127884'}; %Main reactions valine metabolism


geneList = {'ENSG00000100116', 'ENSG00000189221', 'ENSG00000131480'}; %Cysteine desulfurase, mitochondrial


%The whole sulfur metabolism pathway 'ENSG00000139531' is also used by the
%model

geneList = {'ENSG00000112972', 'ENSG00000134240'}; %HMG-CoA synthesis


geneList = {'ENSG00000117115', 'ENSG00000142619', 'ENSG00000142623', 'ENSG00000159339', 'ENSG00000256049'}; %Protein-arginine deiminase


geneList = {'ENSG00000167701', 'ENSG00000060982', 'ENSG00000118520', 'ENSG00000129596', 'ENSG00000104951'}; %Removed reactions due to diff expression
geneList = { 'ENSG00000135423', 'ENSG00000115419'}; %Kidney vs Liver glutaminase


geneList = {'ENSG00000122729'}; %Cytosolic aconitase ACO1

%
geneList ={'ENSG00000138413'}; %IDH1 cytosolic
%geneList = {'ENSG00000112972'}; %HMGCS1 hmg-coa
geneList = {'ENSG00000244005', 'ENSG00000105755'}; %Cysteine desulfurase, mitochondrial

geneList = {'ENSG00000116761'}; %CTH
geneList = {'ENSG00000129596'}; %CDO1
geneList = {'ENSG00000060982', 'ENSG00000105552'}; %BCAT 1 and 2


hold all

maxValue = 0;
for i = 1:length(geneList)
    k = 1 + 2*(i-1);
    [data1, data2 ] = getGeneData(A, B, proteinCodingTranscripts, geneList{i});

    tot1 = sum(data1,1);
    tot2 = sum(data2,1);

    highVal = max([tot1 tot2]);
    
    bar(k, mean(tot1), 'FaceColor', exMap(1,:), 'EdgeColor', 'none')
    bar(k+1, mean(tot2), 'FaceColor', exMap(3,:), 'EdgeColor', 'none')
    errorbar(k, mean(tot1), std(tot1), 'k')
    errorbar(k+1, mean(tot2), std(tot2), 'k')

    [h p] = ttest2(tot1,tot2);
    text(k, highVal, sprintf('p=%2.2e',p));
    
    maxValue = max(maxValue, highVal);
end

set(gca,'xtick',1+1:2:(2*length(geneList)))
set(gca,'xticklabel',geneList)
xtickangle(45)
ylabel('TPM')

set(gca,'fontsize',12)
ylim([0 ceil(maxValue)*1.05])
xlim([0.5 k+1.5])
legend('Liver', 'HepG2')
legend('boxoff')
