function [C, labels, rxnStart, sourceMets, sinkMets] = generateShortestDistance(model, smallSolution, sources, sinks, banMets, banRxns)
    fluxTresh = 10^-5;   
    S = full(model.S);
    removeExchange = true;
    
    %get 
    banRxns =  getNrs(model.rxns, banRxns);
    S(:,banRxns) = 0;    
    
    banMets =  getNrs(model.metNames, banMets);
    S(banMets,:) = 0;
    
    if removeExchange
       [crap, id] = getExchangeRxns(model);
       S(:,id) = 0;
    end
    
    rxnStart = length(model.metNames);
    sourceMets = getNrs(model.metNames, sources);
    sinkMets =  getNrs(model.metNames, sinks);
    sourceAndSink = [sourceMets; sinkMets];

    
    %Generate Weighted bipartiet graph
    C = makebipartite(S, smallSolution);
    absC = sign(C);
    rxnLabels = makeRxnLabels(model);
    labels = [model.metNames;rxnLabels];

    %Calculate shortest distance to source
    NodesToKeepA = getNodesToKeep(absC', sinkMets, sourceAndSink, labels);

    %Calculate shortest distance from source to a node 
    NodesToKeepB = getNodesToKeep(absC, sourceMets, sourceAndSink);
    NodesToKeep = unique([NodesToKeepA NodesToKeepB]);
    
    filterNodes = true(length(labels),1);
    filterNodes(NodesToKeep) = false;
    filterNodes(sourceMets) = false;
    filterNodes(sinkMets) = false;    

    
    C(filterNodes,:) = 0;
    C(:,filterNodes) = 0;    
 
    A = zeros(length(labels),1);
    B = zeros(length(labels),1);
    A(sourceMets) = 1;
    B(sinkMets) = 1;
    sourceMets = A;
    sinkMets = B;
    
    %Trim Graph
    [C, labels, rxnStart, sourceMets, sinkMets] = trimGraph(C, labels, rxnStart, sourceMets, sinkMets);
    
end

function [C, labels, rxnStart, sourceMets, sinkMets] = trimGraph(C, labels, rxnStart, sourceMets, sinkMets)
    for i = length(labels):-1:1
       if and(sum(C(i,:)) == 0, sum(C(:,i)) == 0)
           C(i,:) = [];
           C(:,i) = [];
           sourceMets(i) = [];
           sinkMets(i) = [];
           labels(i) = [];
           if i<=rxnStart
              rxnStart = rxnStart-1; 
           end
       end
        
    end
end

function NodesToKeep = getNodesToKeep(C, sources, sinks, labels)
    NodesToKeep = [];
    
    for i = 1:length(sources)
        curLen = inf;
        curPath = [];
        for j=1:length(sinks)
            if sources(i) ~= sinks(j)
                %[crap, sp] = dijkstra(C, sources(i), sinks(j));
                [cost, sp] = breadthSearch(C, sources(i), sinks(j));
                
                if cost<curLen
                    curPath = [[] sp(:)'];
                    curLen = cost;
                elseif cost == curLen
                    curPath = [curPath sp(:)']; %add all nodes of equal distance
                end
            end
        end
        NodesToKeep = [NodesToKeep curPath];                
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



function result = getNrs(list, sources)
result = zeros(length(sources),1);
    for i = 1:length(sources)
        result(i) = findIndex(list, sources{i});
    end
end

function C = makebipartite(S, smallSolution)
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

function rxnLabels = makeRxnLabels(model)
    rxnLabels = constructEquations(model);
        for i = 1:length(rxnLabels)
            %Pre process equations
            curEq = strrep(rxnLabels{i},'<=>','=>');
            curEq = strrep(curEq,'=>','&rarr;');
            curEq = strrep(curEq,'[c]','');
            rxnLabels{i} = curEq;
        end
end
