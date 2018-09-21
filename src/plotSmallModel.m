function plotSmallModel(model, metThreshold, rxnThreshold)

[cMatrix, labels] = generateCmatrix(model, metThreshold, rxnThreshold, false);

test = biograph(cMatrix, labels);
h = view(test);

for i = 1:length(h.Nodes)
    set(h.Nodes(i),'FontSize',10)
    set(h.Nodes(i),'Color',[0.2 0.7 1])    
    set(h.Nodes(i),'LineColor',[0.5 0.5 0.5]) 
end
set(h, 'ShowArrows', 'off')
set(h, 'Scale', 0.3)
end



function [C, labels] = generateCmatrix(model, metCutOf, rxnCutOf, r2r)
    S = sign(abs(model.S));   %We convert the S-matrix to binary form for topological analysis
    rxnConnect = sum(S,1);
    metConnect = sum(S,2);
    S(metConnect>metCutOf, :) = 0;
    S(:, rxnConnect>rxnCutOf) = 0;
    if r2r == true
        labels = model.rxns;
        C = S' * S;
    else
        metNames = modifyMetNames(model);
        labels = metNames;
        C = S * S';
    end
    
    C = triu(C);   
end
