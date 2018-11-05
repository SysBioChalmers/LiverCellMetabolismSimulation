function solution = filterLowFlux(model, solution, threshold)
%filterLowFlux
% for each metabolite keep all reaction above the treshold
% model- a model structure
% solution - a solution vector
% threshold - the treshold in percent for each metabolite
solution = solution';

keep = zeros(length(solution),1);
S = full(model.S);

for i = 1:size(S,1)
    fluxes = S(i,:) .* solution;
    fluxSum = sum(fluxes(fluxes>0));
    metTresh = fluxSum*threshold;
    keepRxns = find(abs(fluxes)>metTresh);
    keep(keepRxns) = 1;
end
removeRatio = 1-sum(keep)/length(keep);
fprintf('filtered out %2.2f percent of the reactions\n', removeRatio*100);
solution = solution'.*keep;
end
