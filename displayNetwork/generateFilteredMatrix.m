function [C, labels] = generateFilteredMatrix(model, metCutOf, rxnCutOf)
    Z = full(sign(abs(model.S)));   %We convert the S-matrix to binary form for topological analysis    
    labels = [model.metNames];
    
    %Calculate networks scores
    rxnConnect = sum(Z,1);
    metConnect = sum(Z,2);    
    histogram(rxnConnect)
    
    %Remove biomass reaction (can later replace this with fake rxns)
    Z(:, rxnConnect>rxnCutOf) = 0;    
    
    
    [id, exchangeRxns] =getExchangeRxns(model);
    exchangeMets = sum(Z(:,exchangeRxns),2);
    %labels(exchangeMets>0)
    
 
    [metConnect, indx] = sort(metConnect);
    Z = Z(indx,:);
    labels = labels(indx)
    
    prevRxnUse = zeros(1,size(Z,2));
    for i = 1:length(labels)
        curRxns = Z(i,:);
        %curRxns(prevRxnUse
        
        prevRxnUse = prevRxnUse + curRxns;
        Z(i,:) = curRxns;
    end   

    C = Z * Z';
    C(eye(size(C))) = 0; 
    C = triu(C);
end

