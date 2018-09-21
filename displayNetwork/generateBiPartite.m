function [C, labels, rxnStart, exchangeMets, biomassMets] = generateBiPartite(model, smallSolution, metCutOf, biomass, banMets, whiteMets)
    fluxTresh = 10^-5;
    
    [model, smallSolution, exchangeMets, biomassMets] = preProcess(model, smallSolution, biomass, fluxTresh);    
    S = full(model.S);
    Z = sign(abs(S));  
    labels = model.metNames;
    
    S(ismember(labels, banMets),:) = 0;
    exchangeMets(ismember(labels, banMets)) = false;
    biomassMets(ismember(labels, banMets)) = false;
    
    %Calculate networks scores
    metConnect = sum(Z,2);    
    
    [metConnect, indx] = sort(metConnect);
    S = S(indx,:);
    labels = labels(indx);
    exchangeMets = exchangeMets(indx);
    biomassMets = biomassMets(indx);
    
    firstElement = find(metConnect<metCutOf);
    firstElement = firstElement(end)+1;
    
    SPos = sign(S)>0;
    SNeg = sign(S)<0;
    
    currentDepth = metConnect(firstElement);
    currentI = 1;
    
    for i = firstElement:length(labels)
        if currentDepth<metConnect(i)
            currentDepth = metConnect(i);
            
        end
        currentI = i-1;        
        
        
        if not(ismember(labels{i},whiteMets))
            rxns = find(S(i,:));
            rxnsP = rxns(S(i,rxns)>0);
            rxnsN = rxns(S(i,rxns)<0);

            mapP = zeros(1,size(S,2));
            mapN = zeros(1,size(S,2));
%             for j = 1:length(rxnsP)
%                 if sum(SPos(1:currentI,rxnsP(j))) == 0
%                     mapP(rxnsP(j)) = 1;
%                 end
%             end         
% 
%             for j = 1:length(rxnsN)            
%                 if sum(SNeg(1:currentI,rxnsN(j))) == 0
%                     mapN(rxnsN(j)) = 1;
%                 end                              
%             end     
            
%             if sum(mapP) == 0
%                 [crap, idx] = min(smallSolution(rxnsP));
%                 mapP(rxnsP(idx)) = 1;
%             end
%             if sum(mapN) == 0   
%                 [crap, idx] = min(smallSolution(rxnsN));
%                 mapN(rxnsN(idx)) = 1;
%             end  

            S(i,:) = S(i,:) .* (mapP + mapN);
        end
    end    
    
    
    
    %Nonconected Mets
    nonConected = sum(abs(S),2) == 0;
    printRemoval(labels, nonConected)
    
    exchangeMets(nonConected) = [];
    S(nonConected,:)=[];
    labels(nonConected) = [];
    biomassMets(nonConected) = [];
        
    rxnLabels = makeRxnLabels(model, true);
    
    %Nonconected reactions
    nonConected = sum(abs(S),1) == 0;
    printRemoval(model.rxns, nonConected)
    printRemoval(constructEquations(model), nonConected)
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

function [model, smallSolution, exchangeMets, biomassMets] = preProcess(model, smallSolution, biomass, fluxTresh)  
    [id, exchangeRxns] = getExchangeRxns(model);
    
    exchangeMets = sum(model.S(:,exchangeRxns),2) > 0;
    %labels(exchangeMets>0)
    
    model.S(:,exchangeRxns) = [];
    model.rxns(exchangeRxns) = [];
    smallSolution(exchangeRxns) = [];

    %Remove small fluxes
    smalRxns = smallSolution<fluxTresh;
    model.S(:,smalRxns) = [];
    model.rxns(smalRxns) = [];
    smallSolution(smalRxns) = [];    
    
    biomassIndx = ismember(model.rxns, biomass);

    biomassMets = sum(abs(model.S(:,biomassIndx)),2) > 0;    
    model.S(:, biomassIndx) = [];
    model.rxns(biomassIndx) = [];
    smallSolution(biomassIndx) = [];   
end

function printRemoval(list, filter)
    fprintf('Removing unconected objects:\n')
    list = list(filter);
    for i = 1:length(list)
        fprintf('%s\n', list{i})
    end
    fprintf('\n')
end

function rxnLabels = makeRxnLabels(model, printRxn)
    rxnLabels = model.rxns;
    if printRxn == true
        eqs = constructEquations(model);

        for i = 1:length(eqs)
            %Pre process equations
            curEq = strrep(eqs{i},'<=>','=>');
            curEq = strrep(curEq,'[c]','');
            reaPro = strsplit(curEq, '=>');
            reactants = lineBrakeString(reaPro{1});
            products = lineBrakeString(reaPro{2});
            curEq = sprintf('%s\n%s=>\n%s', rxnLabels{i}, reactants, products);
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
