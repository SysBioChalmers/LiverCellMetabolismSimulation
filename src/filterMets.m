function solution = filterMets(model, solution, mets, threshold)
%filterLowFlux
% for each metabolite keep all reaction above the treshold
% model- a model structure
% solution - a solution vector
% the 
% threshold - the treshold in percent for each metabolite
solution = solution';

mustkeep = zeros(length(solution),1);
canRemove = zeros(length(solution),1);
S = full(model.S);

metIds = find(contains(model.metNames,mets));

for i = 1:length(metIds)
    fluxes = S(metIds(i),:) .* solution;
    fluxSum = sum(fluxes(fluxes>0));
    metTresh = fluxSum*threshold;
    mustkeep(abs(fluxes)>metTresh) = 1;
    canRemove(and(abs(fluxes)<metTresh, abs(fluxes)>0)) = 1;
end

canRemove(mustkeep==1) = 0;

%constructEquations(model, canRemove == 1)
removeRatio = sum(canRemove)/length(canRemove);
fprintf('filtered out %2.2f percent of the reactions\n', removeRatio*100);
solution = solution';
solution(canRemove==1) = 0;
end
