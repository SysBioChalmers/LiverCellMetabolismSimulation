function y = fitFlux(p,x,t,volume)
    tNr = x(:,1);
    valueNr = x(:,2);
    y = zeros(size(x,1),1);
    
    predictions = zeros(size(p,1), length(t));
    predictions(1,:) = p(1,1) * exp(p(1,2).*t);
    
    %iterative fit curve between each sample removal points
    predictions(2:end,1) = p(2:end,1);
    for i = 2:length(t)
        dt = t(i)-t(i-1);
        growthPart = predictions(1,i-1)./p(1,2) * (exp(p(1,2).*dt)-1);
        predictions(2:end,i) = predictions(2:end,i-1) + p(2:end,2)/volume(i) * growthPart;
    end
    for i = 1:size(x,1)
        y(i) = predictions(valueNr(i), tNr(i)==t);
    end
end

function predictMetaboliteConcentrations(p, t, volume)
    predictions = zeros(size(p,1)-1, length(t));

end