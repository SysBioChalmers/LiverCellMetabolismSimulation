exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;

load('genericHuman2.mat')
%A = importdata('tisue.txt');
A = importdata('seq.txt');

C = importdata('ID_conversion_key.txt');
proteinCodingTranscripts = C.textdata(:,2);
geneConvertionList = C.textdata(:,[1 5]);

%%%%%%%%%%%%%%%%%%%%%%%
cytoplasmMitMap = mapCytoToMito(model);
cytGenes = model.grRules(cytoplasmMitMap(:,1));
mitGenes = model.grRules(cytoplasmMitMap(:,2));
fprintf('Affected Rxns %1.0f\n', length(cytoplasmMitMap));


%%%%%%%%%%%%%%%%%%%%%%%
withGenes = and(not(cellfun(@isempty,cytGenes)), not(cellfun(@isempty,mitGenes)));
cytoplasmMitMap= cytoplasmMitMap(withGenes,:);
cytGenes = model.grRules(cytoplasmMitMap(:,1));
mitGenes = model.grRules(cytoplasmMitMap(:,2));
fprintf('Affected Rxns with gene rules %1.0f\n', sum(withGenes));


%%%%%%%%%%%%%%%%%%%%%%%
areDifferent = zeros(length(cytoplasmMitMap),1);
for i = 1:length(cytoplasmMitMap)
    areDifferent(i) = not(strcmp(cytGenes{i}, mitGenes{i}));
end
cytoplasmMitMap= cytoplasmMitMap(areDifferent==true,:);
cytGenes = model.grRules(cytoplasmMitMap(:,1));
mitGenes = model.grRules(cytoplasmMitMap(:,2));
fprintf('Have different genes in mit and cyt %1.0f\n', sum(withGenes));

%%%%%%%%%%%%%%%%%%%%%%%%
areUnique = strcat(cytGenes,mitGenes);
[areUnique,identifyerIndex,geneGroupId] = unique(areUnique);
fprintf('Unique gene combos %1.0f\n', length(areUnique));

% for i = 1:length(mitRxns)
%     fprintf('%s\t%s\n', model.rxns{cytRxns(i)}, model.rxns{mitRxns(i)});
% end

%%
%Plot genes
geneExpression = zeros(length(identifyerIndex),3);
geneSTD = zeros(length(identifyerIndex),3);

for i = 1:length(identifyerIndex)
    curCyt = parseGenes(cytGenes{identifyerIndex(i)});
    curMit = parseGenes(mitGenes{identifyerIndex(i)});
    both = intersect(curCyt,curMit);
    curCyt = setdiff(curCyt,both);
    curMit = setdiff(curMit,both);
    
    tot = [0 0];
    for j = 1:length(curCyt)
        [data, transcripts] = extractGeneData(A, proteinCodingTranscripts, curCyt{j});
        tot = tot + sum(data,1);   
    end
    geneExpression(i,1) = mean(tot);
    geneSTD(i,1) = std(tot);
    
    
    tot = [0 0];
    for j = 1:length(curMit)
        [data, transcripts] = extractGeneData(A, proteinCodingTranscripts, curMit{j});
        tot = tot + sum(data,1);
    end  
    geneExpression(i,2) = mean(tot);
    geneSTD(i,2) = std(tot);    
    
    tot = [0 0];
    for j = 1:length(both)
        [data, transcripts] = extractGeneData(A, proteinCodingTranscripts, both{j});
        tot = tot + sum(data,1);
    end
    geneExpression(i,3) = mean(tot);
    geneSTD(i,3) = std(tot);       
end


%%
% for i = 1:length(identifyerIndex)
%     subplot(6,7,i)
%     hold all
%     for j = 1:3
%         bar(j, geneExpression(i,j), 'FaceColor', exMap(j,:), 'EdgeColor', 'none')
%         errorbar(j, geneExpression(i,j), geneSTD(i,j), 'k')
%     end
% end
stdFactor = 2;
diffFactor = 3;
TPMcutOf = 1;

affiliation = zeros(length(identifyerIndex),1);
for i = 1:length(identifyerIndex)
    %cytosolic
    A = geneExpression(i,1) - stdFactor*geneSTD(i,1);
    B = geneExpression(i,2) + stdFactor*geneSTD(i,2);
    C = geneExpression(i,3) + stdFactor*geneSTD(i,3);
    
    if min(geneExpression(i,1:2))>TPMcutOf
        if A>(B + C)*diffFactor
            affiliation(i) = 1;
        end

        A = geneExpression(i,1) + stdFactor*geneSTD(i,1);
        B = geneExpression(i,2) - stdFactor*geneSTD(i,2);
        C = geneExpression(i,3) + stdFactor*geneSTD(i,3);    
        if B>(A + C)*diffFactor
            affiliation(i) = 2;
        end   
    end
end

distinctCompartmentalization = find(affiliation);

fprintf('Distinct compartmentalization %1.0f\n', length(distinctCompartmentalization));
fprintf('-Cytosolic %1.0f\n', sum(affiliation==1));
fprintf('-Mitochondrial %1.0f\n', sum(affiliation==2));


%%
% 
rxnEquations = constructEquations(model, cytoplasmMitMap(:,1));

for i = 1:length(distinctCompartmentalization)
    subplot(3,3,i)
    hold all
    curId = distinctCompartmentalization(i);
    for j = 1:3
        bar(j, geneExpression(curId,j), 'FaceColor', exMap(j,:), 'EdgeColor', 'none')
        errorbar(j, geneExpression(curId,j), geneSTD(curId,j), 'k')
        ylabel('TPM')
        titles = rxnEquations(find(geneGroupId==curId));
        xlabel(titles{1})
    end
end

%%

geneGroups = find(affiliation==1);
for i = 1:length(geneGroups)
    affectedReactions = find(geneGroupId == geneGroups(i));
    for j = 1:length(affectedReactions)
        curRxn = affectedReactions(j);
        fprintf('%s\t%s\n', rxnEquations{curRxn}, cytGenes{curRxn}) 
    end
    fprintf('-\n')
end

fprintf('---------------\n')
rxnEquations = constructEquations(model, cytoplasmMitMap(:,2));
geneGroups = find(affiliation==2);
for i = 1:length(geneGroups)
    affectedReactions = find(geneGroupId == geneGroups(i));
    for j = 1:length(affectedReactions)
        curRxn = affectedReactions(j);
        fprintf('%s\t%s\n', rxnEquations{curRxn}, mitGenes{curRxn}) 
    end
    fprintf('-\n')
end
    
curGenes = mitGenes(identifyerIndex(affiliation==2));
