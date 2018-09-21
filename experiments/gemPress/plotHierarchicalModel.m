function [ output_args ] = plotHierarchicalModel(model, solution, mets, graphDirection)
if nargin <4
    graphDirection = 'right';
end
    
Z = abs(sign(model.S));
[crap, biomassRxn] = max(sum(Z,1));
metNames = modifyMetNames(model);

metNr = zeros(length(mets),1);
for i = 1:length(mets)
    metNr(i) = find(ismember(metNames, mets{i}));
end

plotSubgraph(model, solution, biomassRxn, metNr, graphDirection);
    
end


function plotSubgraph(model, solution, biomassRxn, mets, graphDirection)
    edgeLabelOn = false;
    layeredGraph = true;
    nodeList = [];
    s = [];
    t = [];
    w = [];
    
    reactionMap = [];

    for i = 1:length(mets)
        affectedRxns = find(model.S(mets(i),:));
        metName = model.metNames{mets(i)};
        nodeList = [nodeList; metName];
        reactionMap = [reactionMap; -1000];
        metNode = size(nodeList,1);

        for j = 1:length(affectedRxns)
            if affectedRxns(j) == biomassRxn
                rxnLabel = cellstr(sprintf('Biomass(%s)', metName));
                reactionMap = [reactionMap; -1000];
                nodeList = [nodeList; rxnLabel];
                rxnNode = length(nodeList);
            elseif not(ismember(affectedRxns(j), reactionMap))
                rxnLabel = constructEquations(model, affectedRxns(j));
                reactionMap = [reactionMap; affectedRxns(j)];
                nodeList = [nodeList; rxnLabel]; 
                rxnNode = length(nodeList);
            else
                rxnNode = find(reactionMap == affectedRxns(j));
            end
                reactionWeight = full(model.S(mets(i), affectedRxns(j)));
                direction = sign(reactionWeight);
                reactionWeight = reactionWeight * solution(affectedRxns(j));
                reactionWeight = reactionWeight * 1000;
                     
                if direction == 1
                    s = [s, rxnNode]; 
                    t = [t, metNode];                
                else
                    s = [s, metNode]; 
                    t = [t, rxnNode];
                end
                w = [w, abs(reactionWeight)];
                
        end
    end

    
     p=digraph(s,t, w);
     
     if layeredGraph
         H = plot(p,'Layout','layered', 'Direction', graphDirection, 'NodeLabel', '', 'ArrowSize', 15);
     else
         H = plot(p, 'NodeLabel', '', 'ArrowSize', 15);
     end
         
     if edgeLabelOn
         labeledge(H,s,t,w)
     end
     
     wNorm = p.Edges.Weight;
     wNorm = 5*wNorm/max(wNorm);
     H.LineWidth = wNorm;
     
    for j = 1:length(H.XData)
        rxnLabel = strrep(nodeList{j},'<=>', '=>');
        rxnLabel = strrep(rxnLabel,'=>', '\n=>\n');
        text(H.XData(j)+0.05, H.YData(j), sprintf(rxnLabel), 'fontSize', 6);
    end
    
    axis tight
    axis off
    %set(findobj(gcf, 'type','axes'), 'Visible','off');
    set(gcf,'color','w');
    set(findall(gcf,'type','text'),'FontSize',6,'fontWeight','bold')
    nodeList
end
