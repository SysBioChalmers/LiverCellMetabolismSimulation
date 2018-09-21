function [ output_args ] = plotHierarchicalModel(model, solution, mets)
Z = abs(sign(model.S));
[crap, biomassRxn] = max(sum(Z,1));
metNames = modifyMetNames(model);

for i = 1:length(mets)
    subplot(length(mets),1,i)
    metNr = find(ismember(metNames, mets{i}));
    plotSubgraph(model, solution, biomassRxn, metNr);
end

end


function plotSubgraph(model, solution, biomassRxn, startMet)

    affectedRxns = find(model.S(startMet,:));
    metList = model.metNames(startMet);
    s = [];
    t = [];
    w = [];
    for i = 1:length(affectedRxns)
        if affectedRxns(i) == biomassRxn
            rxnLabel = 'Biomass';
        else
            rxnLabel = constructEquations(model, affectedRxns(i));
            rxnLabel = strrep(rxnLabel,'<=>', '=>');
            rxnLabel = strrep(rxnLabel,'=>', '\n=>\n');
        end
            reactionWeight = full(model.S(startMet, affectedRxns(i)));
            direction = sign(reactionWeight);
            reactionWeight = reactionWeight * solution(affectedRxns(i));
            reactionWeight = reactionWeight * 1000;
            metList = [metList; rxnLabel];      
            if direction == 1
                s = [s, length(metList)]; 
                t = [t, 1];                
            else
                s = [s, 1]; 
                t = [t, length(metList)];
            end
            w = [w, abs(reactionWeight)];
    end

    
     p=digraph(s,t);
     
     H = plot(p,'Layout','layered', 'Direction','right','NodeLabel', '', 'ArrowSize', 15);
     labeledge(H,s,t,w)
     wNorm = cellfun(@str2num,H.EdgeLabel);
     wNorm = 5*wNorm/max(wNorm);
     H.LineWidth = wNorm;
     
    for i = 1:length(H.XData)
       text(H.XData(i)+0.05, H.YData(i), sprintf(metList{i}), 'fontSize', 6);
    end
    
    axis tight
    axis off
    %set(findobj(gcf, 'type','axes'), 'Visible','off');
    set(gcf,'color','w');
    xlim([0.5 10])   
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
