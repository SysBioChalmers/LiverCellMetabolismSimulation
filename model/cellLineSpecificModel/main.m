addpath('C:\Program Files\SBML\libSBML-5.10.2-libxml2-x64\bindings\matlab\matlab')
modelHep=importModel('Hep-G2.xml');

load('../genericHuman2.mat')

%ignore knocks that have no gene association
hepRxns = modelHep.rxns;
genRxns = model.rxns;
hepKnock = setdiff(genRxns, hepRxns);



%remove aquaporin mitochondrial anotation:
model.grRules{findIndex(model.rxns,'HMR_4888')} = '';
model.grRules{findIndex(model.rxns,'HMR_4951')} = '';

noRules = model.rxns(ismember(model.grRules,''));
hepKnock = setdiff(hepKnock, noRules);

length(intersect(hepKnock, noRules))




%Make exceptions for important reactions:
hepKnock(findIndex(hepKnock, 'HMR_4446')) = []; %Add back methylenetetrahydrofolate reductase or else methionine depletion
hepKnock(findIndex(hepKnock, 'HMR_6725')) = []; %Add back GOT for phenylpyruvate synthesis
hepKnock(findIndex(hepKnock, 'HMR_7757')) = []; %Add back DBI for succinyl-CoA transport (dead end metabolite)
hepKnock(findIndex(hepKnock, 'HMR_4199')) = []; %Add back Aminoacylase for N-acetyl-L-alanine[c] synthesis (dead end metabolite)
hepKnock(findIndex(hepKnock, 'HMR_8416')) = []; %Add back NAGS for N-acetylglutamate synthesis (dead end metabolite)
hepKnock(findIndex(hepKnock, 'HMR_5113')) = []; %Add back alanine transport, removed due to non existant gene TP250
hepKnock(findIndex(hepKnock, 'HMR_4790')) = []; %Add back alanine transport, removed due to non existant gene TP250
hepKnock(findIndex(hepKnock, 'HMR_4790')) = []; %Add back alanine transport, removed due to non existant gene TP250
hepKnock(findIndex(hepKnock, 'HMR_8684')) = []; %GOT2 excists in HepG2, removed due to orphan metabolites?


% for i = 1:length(hepKnock)
%    rxn =  findIndex(model.rxns,hepKnock{i});
%    genes = model.grRules(rxn);
%    if isempty(strfind(genes, 'ENSG'))
%        disp(genes)
%    end
% end

fid = fopen('hepKnock.txt','w');
fprintf(fid,'%s\n', hepKnock{:});
fclose(fid);

rxns = constructEquations(model, hepKnock);

fid = fopen('hepKnockRxns.txt','w');
fprintf(fid,'%s\n', rxns{:});
fclose(fid);
