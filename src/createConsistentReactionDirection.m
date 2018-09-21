function model = createConsistentReactionDirection(model)
% createConsistentReactionDirection
% Changes the S matrix so that all exchange reactions are negative for 
% uptake and positive for release. (by setting them to - in the S matrix)
%
%   model             a model structure
%   model             an updated model structure
%
%   Avlant Nilsson, 2016-05-16
%
    [exchangeRxns,exchangeRxnsIndexes] = getExchangeRxns(model,'both');
    for i = 1:length(exchangeRxnsIndexes)
        %Exchange reactions are reactions with only one entry in the column
        %in the S matrix
        if abs(sum(model.S(:, exchangeRxnsIndexes(i)))) ~= sum(abs(model.S(:, exchangeRxnsIndexes(i))))
           disp('Warning, antiport reactions are not handled by this script')
        else
            model.S(:, exchangeRxnsIndexes(i)) = -abs(model.S(:, exchangeRxnsIndexes(i)));
        end
    end

end

