function [C, labels, rxnStart] = generateBiPartiteFromMets(model, smallSolution, mets, fluxTresh)

    S = full(model.S);

    labels = modifyMetNames(model);
    metMatch = ismember(labels, mets);
    
    includedReactions = sum(abs(S(metMatch,:)),1) > 0;
    
    includedReactions(smallSolution<fluxTresh) = false;
    
    S(:,not(includedReactions)) = 0;
    S(not(metMatch),:) = 0;
        
    SPos = sign(S)>0;
    SNeg = sign(S)<0;
    
    %Nonconected Mets
    nonConected = sum(abs(S),2) == 0;
    %printRemoval(labels, nonConected)
    
    S(nonConected,:)=[];
    labels(nonConected) = [];
        
    rxnLabels = makeRxnLabels(model, true);
    
    %Nonconected reactions
    nonConected = sum(abs(S),1) == 0;
    %printRemoval(model.rxns, nonConected)
    %printRemoval(constructEquations(model), nonConected)
    S(:,nonConected) = [];
    smallSolution(nonConected) = [];
    rxnLabels(nonConected) = [];   
    
    rxnStart = length(labels);
    labels = [labels; rxnLabels];
    
    
    weightedS = S .* repmat(smallSolution', size(S,1),1);
    
    
    SPos = weightedS; 
    SNeg = weightedS;        
    
    Q1 = zeros(size(S,1));
    SPos(SPos<0) = 0;
    SNeg(SNeg>0) = 0;
    SNeg = -SNeg;
    
    Q4 = zeros(size(S,2));
    C = [Q1 SNeg; SPos' Q4];   
    
end

function rxnLabels = makeRxnLabels(model, printRxn)
    rxnLabels = model.rxns;
    if printRxn == true
        eqs = constructEquations(model);

        for i = 1:length(eqs)
            %Pre process equations
            equationName = rxnLabels{i};
            equationName = '';
            curEq = strrep(eqs{i},'<=>','=>');
            reaPro = strsplit(curEq, '=>');
            reactants = lineBrakeString(reaPro{1});
            products = lineBrakeString(reaPro{2});
            
            curEq = sprintf('%s\n%s=>\n%s', equationName, reactants, products);
            rxnLabels{i} = curEq;
        end
    end
end

function result = lineBrakeString(data)
    elementsPerRow = 3;
    result = [];
    data = strsplit(data, ' + ');
    intervals = 1:elementsPerRow:length(data);
    
    for i = 1:length(intervals)
        maxInter = intervals(i)+elementsPerRow-1;
        maxInter = min(maxInter, length(data));
        curInterval = intervals(i):maxInter;
        curLine = strjoin(data(curInterval), ' + ');
        result = [result curLine '\n'];
    end
end
