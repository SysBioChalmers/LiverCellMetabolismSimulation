function samples = randomSample(N, sequences, method)
    nrAA = size(sequences, 2);
    nrProteins = size(sequences,1);

    samples = zeros(N, nrAA);
    
    if strcmp(method, 'uni')
        randVector = rand(nrProteins,N)*100;
    elseif strcmp(method, 'logNormal')
        randVector = 10.^(randn(nrProteins,N)*1.0256 + 0.8561);
    end

    
    for i = 1:N
        currentMatrix = repmat(randVector(:,i), 1, nrAA);
        currentMatrix = currentMatrix.*sequences;
        currentDistribution = sum(currentMatrix);
        samples(i,:) = currentDistribution/sum(currentDistribution);
    end
end

