function [ output_args ] = plotHierarchicalModel(model, solution, mets, graphDirection)
if nargin <4
    graphDirection = 'right';
end
    
biomassRxn = find(contains(model.rxns,'HumanGrowth'));
metNames = modifyMetNames(model);

metNr = zeros(length(mets),1);
for i = 1:length(mets)
    metNr(i) = find(ismember(metNames, mets{i}));
end

plotSubgraph(model, solution, biomassRxn, metNr, graphDirection);
    
end


function plotSubgraph(model, solution, biomassRxn, mets, graphDirection)
    metNames = modifyMetNames(model);
    edgeLabelOn = false;
    layeredGraph = true;
    nodeList = [];
    fluxList = [];
    s = [];
    t = [];
    w = [];
    
    reactionMap = [];

    for i = 1:length(mets)
        affectedRxns = find(model.S(mets(i),:));
        metName = metNames{mets(i)};
        nodeList = [nodeList; metName];
        fluxList = [fluxList; 0];
        reactionMap = [reactionMap; -1000];
        metNode = size(nodeList,1);

        for j = 1:length(affectedRxns)
            reactionWeight = full(model.S(mets(i), affectedRxns(j)));
            reactionWeight = reactionWeight * solution(affectedRxns(j));
            reactionWeight = reactionWeight * 1000;            
            
            if affectedRxns(j) == biomassRxn
                rxnLabel = cellstr(sprintf('%s => Biomass', metName));
                reactionMap = [reactionMap; -1000];
                nodeList = [nodeList; rxnLabel];
                rxnNode = length(nodeList);
                fluxList = [fluxList; abs(reactionWeight)];
            elseif not(ismember(affectedRxns(j), reactionMap))
                rxnLabel = constructEquations(model, affectedRxns(j));
                reactionMap = [reactionMap; affectedRxns(j)];
                nodeList = [nodeList; rxnLabel]; 
                rxnNode = length(nodeList);
                fluxList = [fluxList; abs(1000*solution(affectedRxns(j)))];
            else
                rxnNode = find(reactionMap == affectedRxns(j));
            end

                direction = sign(reactionWeight);

                
                     
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
     
     wVal = p.Edges.Weight;
     wNorm = 5*wVal/max(wVal);
     H.LineWidth = wNorm;
     
    for j = 1:length(H.XData)
        rxnLabel = strrep(nodeList{j},'<=>', '=>');
        labeledArrow = sprintf('\n(%s) =>\n', num2str(abs(fluxList(j))));
        rxnLabel = strrep(rxnLabel,'=>', labeledArrow);
        text(H.XData(j)+0.05, H.YData(j), rxnLabel, 'fontSize', 6);
    end
    
    axis tight
    axis off
    %set(findobj(gcf, 'type','axes'), 'Visible','off');
    set(gcf,'color','w');
    set(findall(gcf,'type','text'),'FontSize',6,'fontWeight','bold')
    for i = 1:length(nodeList)
       fprintf('%s\n', nodeList{i}); 
    end
end

function [exchangeMets, production] = getExchangeMets(model)
    [id, exchangeRxns] = getExchangeRxns(model);
    exchangeMets = zeros(length(id),1);
    production = zeros(length(id),1);
    for i = 1:length(id)
        exchangeMets(i) =  find(model.S(:,exchangeRxns(i)));
        production(i) = -sign(model.S(exchangeMets(i), exchangeRxns(i)));
    end
end
