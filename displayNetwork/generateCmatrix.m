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
        labels = model.metNames;
        C = S * S';
    end
    
    C = triu(C);   
end

