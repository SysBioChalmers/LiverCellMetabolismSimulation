%This loads the model, makes some changes (adds rxns, makes the rxn
%directionality uniform etc.) and saves the result as a matlab struct for
%future use.

%load Raven Model
addpath('../src')
addpath('cellLineSpecificModel')
model = importExcelModel('HMRdatabase2_00.xls');

%Standardize the direction of exchange fluxes, negative flux = uptake
model = createConsistentReactionDirection(model);

%Make sure that NADPH is used in anabolism and NADH produced in catabolism
model = fixNADHDirectionality(model);

%Fix P/O ratio
model = configureSMatrix(model, 3, 'HMR_6916', 'H+[c]');
model = configureSMatrix(model, -3, 'HMR_6916', 'H+[m]');

%Add fat reactions
model = addFatExchange(model, 's', 'AddedFatExchange');

%Update Fatty acid reactions
model = configureSMatrix(model, 0, 'ApproxPhosphatidate', 'palmitate[c]');
molRatio=2*[0.050037794
0.19337822
0.501334366
0.095723088
0.041258871
0.118267661];
    
FA = {'myristic acid[c]'
'pentadecylic acid[c]'
'palmitate[c]'
'margaric acid[c]'
'oleate[c]'
'stearate[c]'};

for i = 1:length(molRatio)
    model = configureSMatrix(model, molRatio(i), 'ApproxPhosphatidate', FA{i});
end

%heparan sulfate synthesis reaction
lactRxn = createRXNStuct(model, 'heparanSulfateExport', 'heparan sulfate[g] => heparan sulfate[c]', 0, 1000, 'Putative Reactions');
model=addRxns(model,lactRxn,3,'c',true);
lactRxn = createRXNStuct(model, 'heparanSulfateSynthesis', 'heparan sulfate proteoglycan[g] => heparan sulfate[g] + 3-beta-D-glucuronosyl-3-beta-D-galactosyl-4-beta-D-galactosyl-O-beta-D-xylosylprotein[g]', 0, 1000, 'Putative Reactions');
model=addRxns(model,lactRxn,3,'c',false);

%Remove proton gradient requirement for peroxisome transport of pyruvate
model = configureSMatrix(model, 0, 'HMR_4930', 'H+[c]');
model = configureSMatrix(model, 0, 'HMR_4930', 'H+[p]');

%Add 5-oxoproline as an extracellular metabolite
lactRxn = createRXNStuct(model, 'PyroglutamicTransport', '5-oxoproline[s] <=> 5-oxoproline[c]', -1000, 1000, 'Transport, extracellular');
model=addRxns(model,lactRxn,3,'c',true);
lactRxn = createRXNStuct(model, 'PyroglutamicUptake', '5-oxoproline[s] <=> ', -1000, 1000, 'Transport, extracellular');
model=addRxns(model,lactRxn,3,'s',false);

%Remove old albumin
model = setParam(model, 'ub', 'HMR_5151', 0);
model = setParam(model, 'ub', 'HMR_5151', 0);

save('genericHuman', 'model')

%The compartment of this reaction is golgi apparatus and vasicles
%(proteinatlas.org) Blocking this reaction prevents a deamination cycle in
%cytoplasm that shares reaction steps with glycolysis and thereby becomes
%parsimonious. 
model.ub(findIndex(model.rxns, 'HMR_4299')) = 0;

%Aldo-keto reductase family 1 is downregulated 400 fold
%Also unclear how this gene has anything to do with sulfur metabolism
model = setParam(model, 'ub', 'HMR_4840', 0);

%Aminoadipate transporter only proceds in reverse direction (reactome.org)
model.ub(findIndex(model.rxns, 'HMR_8022')) = 0;

%Glutaminase is not cytosolic 
model = setParam(model, 'lb', 'HMR_9802', 0);
model = setParam(model, 'ub', 'HMR_9802', 0);

%GRHPR are not mitochondrial
model = setParam(model, 'ub', 'HMR_8779', 0);

% The function of KYAT1 is to break down tryptophane, the following is a
% moonlight activity of limited quantitative importance
model = setParam(model, 'lb', 'HMR_4196', 0); %oxoglutamarate
model = setParam(model, 'ub', 'HMR_4196', 0); %

%No support for this pantetheine reaction in human
model.ub(findIndex(model.rxns,'HMR_4716')) =0;

%these pantetheine relate reactions are irreversible 
model.ub(findIndex(model.rxns,'HMR_4730')) =0;
model.ub(findIndex(model.rxns,'HMR_4731')) =0;
model.lb(findIndex(model.rxns,'HMR_4679')) =0;
model.lb(findIndex(model.rxns,'HMR_4717')) =0;

%ATP is the cofactor of PPCS in human
model.eccodes{findIndex(model.rxns,'HMR_4723')} = 'EC:6.3.2.51';
model = configureSMatrix(model, 0, 'HMR_4723', 'CTP[c]');
model = configureSMatrix(model, 0, 'HMR_4723', 'CMP[c]');
model = configureSMatrix(model, 1, 'HMR_4723', 'ATP[c]');
model = configureSMatrix(model, -1, 'HMR_4723', 'AMP[c]');

%CTP is thought to operate in oposite direction eg pmid: 18406340  
model.lb(findIndex(model.rxns,'HMR_4964')) = -1000;
%Prevents free proton gradient:
model = configureSMatrix(model, 0, 'HMR_4964', 'H+[m]');
model = configureSMatrix(model, 0, 'HMR_4964', 'H+[c]');

%Acts on hydroxyproline (uniprot Q9UF12)
%model.ub(findIndex(model.rxns,'HMR_3838')) =0;

%No support for this reaction CIT, ISO antiporter,
%in the reference (PMID 14598172)
model.lb(findIndex(model.rxns,'HMR_4972')) =0;
model.ub(findIndex(model.rxns,'HMR_4972')) =0;

%Provides TCA with substrate protein atlas (ENSG00000183048)
%model.ub(findIndex(model.rxns,'HMR_4865')) =0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mitochondrial lactate dehydrogenase bypasses the pyruvate carrier and need
%for malate-aspartate shutle
model = setParam(model, 'ub', 'HMR_4280', 0);
model = setParam(model, 'lb', 'HMR_4280', 0);

%Remove pyr to dATP, dGTP
model = setParam(model, 'lb', 'HMR_4421', 0);
model = setParam(model, 'ub', 'HMR_4421', 0);
model = setParam(model, 'lb', 'HMR_4573', 0);
model = setParam(model, 'ub', 'HMR_4573', 0);
model = setParam(model, 'lb', 'HMR_4487', 0);
model = setParam(model, 'ub', 'HMR_4487', 0);
model = setParam(model, 'lb', 'HMR_4210', 0);
model = setParam(model, 'ub', 'HMR_4210', 0);
model = setParam(model, 'ub', 'HMR_4193', 0);
model = setParam(model, 'ub', 'HMR_4171', 0);

%The gene codes for transporter of malate not AKG 
model = setParam(model, 'lb', 'HMR_6330', 0);
model = setParam(model, 'ub', 'HMR_6330', 0);
model = setParam(model, 'lb', 'HMR_4851', 0);
model = setParam(model, 'ub', 'HMR_4851', 0);

%Gene codes for transport of malate not GSH
model = setParam(model, 'lb', 'HMR_6391', 0);
model = setParam(model, 'ub', 'HMR_6391', 0);

%Remove undocumented acetoacetate tranporter
model = setParam(model, 'lb', 'HMR_5426', 0);
model = setParam(model, 'ub', 'HMR_5426', 0);

%Remove undocumented HMG-CoA tranporter
model.lb(findIndex(model.rxns, 'HMR_1572')) = 0;
model.ub(findIndex(model.rxns, 'HMR_1572')) = 0;

%Remove undocumented threonine conversion
model.ub(findIndex(model.rxns, 'HMR_4284')) = 0;

%remove undocumented transport of AKG
model.lb(findIndex(model.rxns, 'HMR_6293')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6293')) = 0;

model.lb(findIndex(model.rxns, 'HMR_6289')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6289')) = 0;

model.lb(findIndex(model.rxns, 'HMR_6286')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6286')) = 0;

%Allow reversed IDH flux in the mitochondria
%Formally correct but creates an annoying loop
%Mainly relevant for C13 studies
%model.lb(findIndex(model.rxns, 'HMR_3958')) = -1000;

%propanoate directionality, else cells can produce atp from 
%BCAA->propanoate
model = setParam(model, 'ub', 'HMR_0153', 0);
model = setParam(model, 'ub', 'HMR_3797', 0);
model = setParam(model, 'lb', 'HMR_4459', 0);

%Prevent incorrect proline synthesis reversibility
model = setParam(model, 'ub', 'HMR_3806', 0);

%Prevent Glutamine->Alanine from unspecific Transaminase
model = setParam(model, 'ub', 'HMR_4197', 0);

%Free glutamate Transport, model infeasible if no glutamate transport available 
%extracellularly
lactRxn = createRXNStuct(model, 'FreeGlutamateTransport', 'glutamate[s] <=> glutamate[c]', -1000, 1000, 'Transport, extracellular');
model=addRxns(model,lactRxn,3,'c',false);

%Free lactate Transport
lactRxn = createRXNStuct(model, 'FreeLactateTransport', 'L-lactate[s] <=> L-lactate[c]', -1000, 1000, 'Transport, extracellular');
model=addRxns(model,lactRxn,3,'c',false);

%Allow export of alanine and serine from mitochondria
model = setParam(model, 'lb', 'HMR_5113', -1000);
model = setParam(model, 'lb', 'HMR_5114', -1000);

%Allow aspartate export from mitochondria
% lactRxn = createRXNStuct(model, 'aspartateExporter', 'aspartate[m] => aspartate[c]', 0, 1000, 'Putative Reactions');
% model=addRxns(model,lactRxn,3,'c',false);

%reconstruction of cysteine de sulfurase pathway:
%Significantly upregulated (ENSG00000244005)
lactRxn = createRXNStuct(model, 'cysteineDesulfurase', 'cysteine[c] + GSH[c] => S-Sulfanylglutathione[c] + alanine[c]', 0, 1000, 'Sulfur metabolism');
model=addRxns(model,lactRxn,3,'m',true);

%Significantly upregulated  (ENSG00000105755)
lactRxn = createRXNStuct(model, 'persulfideDioxygenase', 'S-Sulfanylglutathione[c] + O2[c] + H2O[c] => GSH[c] + sulfite[c]', 0, 1000, 'Sulfur metabolism');
model=addRxns(model,lactRxn,3,'c',false);

%Make the cell specific knock outs
model = hepG2Constrain(model);

%Manual curation of HEPG2 specific expression
%Probably no expression of Interleukin 4 induced 1, also unclear if
%methionine is its primary substrate
model = setParam(model, 'lb', 'HMR_5390', 0); %Interleukin 4 induced 1
model = setParam(model, 'ub', 'HMR_5390', 0); %

%CDO1 is not expressed in HEPG2
model = setParam(model, 'ub', 'HMR_3908', 0);

%The mitochondrial arginase is expressed 100 fold higher
model = setParam(model, 'ub', 'HMR_3816', 0);

%Probably no mitochondrial BCAT in HEPG2:
model = setParam(model, 'lb', 'HMR_3778', 0);
model = setParam(model, 'ub', 'HMR_3778', 0);
model = setParam(model, 'lb', 'HMR_3744', 0);
model = setParam(model, 'ub', 'HMR_3744', 0);
model = setParam(model, 'lb', 'HMR_3765', 0);
model = setParam(model, 'ub', 'HMR_3765', 0);

%Proably no expression of cytsosolic GPT
model = setParam(model, 'lb', 'HMR_3899', 0); %cytosolic Ala - pyr
model = setParam(model, 'ub', 'HMR_3899', 0); %cytosolic Ala - pyr

%Reactions required for heparan sulfate synthesis reaction
model.ub(findIndex(model.rxns, 'HMR_7220'))=1000;
model.ub(findIndex(model.rxns, 'HMR_7221'))=1000;
model.ub(findIndex(model.rxns, 'HMR_7222'))=1000;

model = calculateAllMissingMetformulas(model);

%Overwrite raven model
save('genericHuman2', 'model')
