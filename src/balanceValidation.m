function result = balanceValidation(model, solution)
    validateBalances = zeros(length(model.mets),1);
    for i = 1:length(model.mets)
        a=solution' .* model.S(i,:);
        validateBalances(i) = sum(a);
    end

    result = sum(abs(validateBalances));
end

